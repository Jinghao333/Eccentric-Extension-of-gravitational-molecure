"""One-folder-per-run, open-channel-only scientific storage for Xu Appendix A.

Canonical successful run:
    config.json  effective simulation inputs (written BEFORE computation)
    physics.h5   unique sampled Fourier coefficients; total open-channel sources;
                 physical complex outgoing amplitudes; partial Green overlaps.

No radial array, numerical error grids, orbital/atomic profiles, fluxes,
Wronskians or duplicated per-channel forcing multipoles are stored.
The radius grid is reconstructed deterministically from config; a checksum
warns if a different code version would silently change that grid.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import h5py
import numpy as np

SCHEMA = 'xu-fig8-open-v4'
GRID_RECIPE = 'geomspace_plus_inclusive_peri_apo_v1'


def radial_grid(config: dict[str, Any]) -> np.ndarray:
    """Exact numerical grid recipe, in units of the central mass coordinate."""
    minimum, maximum = float(config['r_min']), float(config['r_max'])
    base = np.geomspace(minimum, maximum, int(config['log_points']))
    e, p, M = float(config['e']), float(config['p']), float(config['M'])
    n_orb = int(config['orbit_interval_points'])
    if n_orb:
        r_p, r_a = M*p/(1+e), M*p/(1-e)
        left, right = max(minimum, r_p), min(maximum, r_a)
        if left < right:
            base = np.concatenate((base, np.linspace(left, right, n_orb)))
    values = np.unique(base)
    if len(values) < 4 or not np.all(np.isfinite(values)) or np.any(np.diff(values) <= 0):
        raise ValueError('invalid reconstructed radial grid')
    return values.astype('<f8', copy=False)


def grid_digest(values: np.ndarray) -> str:
    """Check the exact sampling grid without redundantly storing the grid."""
    array = np.ascontiguousarray(values, dtype='<f8')
    return hashlib.sha256(array.tobytes()).hexdigest()


def write_config_first(directory: Path, config: dict[str, Any]) -> None:
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=False)
    path = directory/'config.json'
    with path.open('x', encoding='utf-8') as f:
        json.dump(config, f, indent=2, sort_keys=True, allow_nan=False)
        f.write('\n')


def _matrix(h5: h5py.File, name: str, data: np.ndarray) -> None:
    data = np.asarray(data)
    if data.size and data.ndim == 2:
        h5.create_dataset(name, data=data, compression='gzip', compression_opts=3,
                          shuffle=True, chunks=(min(len(data), 16), data.shape[1]))
    else:
        h5.create_dataset(name, data=data)


def write_physics(path: Path, *, config: dict, forcing_cache: dict,
                  source_cache: dict[tuple[int,int,int], np.ndarray],
                  amplitude_cache: dict[tuple[int,int,int], complex],
                  sampled_zero_keys: list[tuple[int,int,int]],
                  partial_records: list[dict], provenance: dict | None = None,
                  quality: dict | None = None) -> None:
    """Atomic HDF5 writer with no duplicate (l,m,g) or open-channel arrays."""
    path = Path(path)
    if path.exists():
        raise FileExistsError(path)
    grid = radial_grid(config)
    n_r = grid.size
    ckeys = sorted(source_cache)
    pvals = sorted(partial_records, key=lambda x: (x['ell_f'], x['m_f'], x['g'], x['forcing_l']))
    # A batched forcing call may compute g values belonging to other sectors.
    # Persist ONLY coefficients actually needed by computed open channels.
    needed_forcing = {(int(p['forcing_l']), int(p['m_f'])-int(config['m_i']), int(p['g']))
                      for p in pvals}
    absent = needed_forcing - set(forcing_cache)
    if absent:
        raise ValueError(f'Missing required Fourier keys: {sorted(absent)[:8]}')
    fkeys = sorted(needed_forcing)
    zkeys = sorted(set(tuple(int(x) for x in key) for key in sampled_zero_keys))
    if ckeys != sorted(amplitude_cache):
        raise ValueError('Source and amplitude keys disagree')
    if set(zkeys) & set(ckeys):
        raise ValueError('An open channel cannot be both zero and computed')
    if len(fkeys)!=len(set(fkeys)) or len(ckeys)!=len(set(ckeys)):
        raise ValueError('Duplicate keys')
    for k in fkeys:
        obj = forcing_cache[k]
        values = np.asarray(obj.coefficient, dtype=np.complex128)
        if values.shape != (n_r,) or not np.all(np.isfinite(values)):
            raise ValueError(f'Invalid Fourier samples {k}')
        if not np.all(np.asarray(obj.quadrature_converged,dtype=bool)):
            raise ValueError(f'Unaccepted Fourier integral {k}')
    sources = np.empty((len(ckeys),n_r),dtype=np.complex128)
    amps = np.empty((len(ckeys),),dtype=np.complex128)
    for i,key in enumerate(ckeys):
        sources[i] = np.asarray(source_cache[key], dtype=np.complex128)
        amps[i] = complex(amplitude_cache[key])
    if not np.all(np.isfinite(sources)) or not np.all(np.isfinite(amps)):
        raise ValueError('Nonfinite source or outgoing amplitude')
    pkeys = np.asarray([(p['ell_f'],p['m_f'],p['g'],p['forcing_l']) for p in pvals],dtype='<i4').reshape(-1,4)
    Ip = np.asarray([complex(p['I_l_re'],p['I_l_im']) for p in pvals],dtype=np.complex128)
    if len({tuple(row) for row in pkeys}) != len(pkeys):
        raise ValueError('Duplicate partial-overlap key')
    if not np.all(np.isfinite(Ip)):
        raise ValueError('Nonfinite partial overlap')
    if any(tuple(row[:3]) not in source_cache for row in pkeys):
        raise ValueError('Partial overlap for unknown output channel')
    tmp = path.with_suffix(path.suffix+'.incomplete')
    if tmp.exists():
        tmp.unlink()
    try:
        with h5py.File(tmp,'w') as h5:
            h5.attrs.update({
                'schema': SCHEMA,
                'complete': False,
                'grid_recipe': GRID_RECIPE,
                'grid_sha256': grid_digest(grid),
                'n_radius': n_r,
                'fourier_basis': 'coordinate_time',
                'Fourier_sign': 'exp(+i*g*mean_anomaly)',
                'spherical_harmonics': 'Condon-Shortley',
                'Z_infinity_definition': 'physical A_out=exp(pi*eta/2)*(sum_l I_l)/W',
                'Z_infinity_kinematics': 'Xu KG flux momentum; NR Coulomb response',
            })
            if provenance:
                h5.attrs['provenance_json'] = json.dumps(provenance,sort_keys=True)
            if quality:
                h5.attrs['quality_json'] = json.dumps(quality,sort_keys=True)
            _matrix(h5,'forcing/keys_lmg',np.asarray(fkeys,dtype='<i4').reshape(-1,3))
            _matrix(h5,'forcing/u_lmg',np.asarray([
                np.asarray(forcing_cache[k].coefficient,dtype=np.complex128) for k in fkeys
            ],dtype=np.complex128).reshape(-1,n_r))
            _matrix(h5,'channels/keys_ell_m_g',np.asarray(ckeys,dtype='<i4').reshape(-1,3))
            _matrix(h5,'channels/S_total',sources)
            _matrix(h5,'channels/Z_infinity',amps)
            _matrix(h5,'channels/sampled_zero_keys',np.asarray(zkeys,dtype='<i4').reshape(-1,3))
            _matrix(h5,'overlap/keys_ell_m_g_l',pkeys)
            _matrix(h5,'overlap/I_l',Ip)
            h5.attrs['complete'] = True
        tmp.replace(path)
    except BaseException:
        if tmp.exists():
            tmp.unlink()
        raise


def load_primitives(directory: Path) -> dict:
    """Read and validate the v4 run; reconstruct the radial grid from config."""
    directory=Path(directory)
    config=json.loads((directory/'config.json').read_text(encoding='utf-8'))
    grid=radial_grid(config)
    with h5py.File(directory/'physics.h5','r') as h5:
        if h5.attrs.get('schema')!=SCHEMA or not bool(h5.attrs.get('complete')):
            raise ValueError('Not a complete v4 run')
        if (h5.attrs.get('grid_recipe')!=GRID_RECIPE or
            h5.attrs.get('grid_sha256')!=grid_digest(grid) or
            int(h5.attrs.get('n_radius'))!=len(grid)):
            raise ValueError('Recreated radial grid differs from the executed grid')
        result={
            'config':config, 'radius':grid,
            'forcing_keys':np.asarray(h5['forcing/keys_lmg']),
            'forcing':np.asarray(h5['forcing/u_lmg']),
            'channel_keys':np.asarray(h5['channels/keys_ell_m_g']),
            'source':np.asarray(h5['channels/S_total']),
            'Z_infinity':np.asarray(h5['channels/Z_infinity']),
            'zero_keys':np.asarray(h5['channels/sampled_zero_keys']),
            'partial_keys':np.asarray(h5['overlap/keys_ell_m_g_l']),
            'partial_I':np.asarray(h5['overlap/I_l']),
            'quality':json.loads(h5.attrs.get('quality_json','{}')),
            'provenance':json.loads(h5.attrs.get('provenance_json','{}')),
        }
    fk,ck,pk,zk=(result[x] for x in ('forcing_keys','channel_keys','partial_keys','zero_keys'))
    n_r=len(grid)
    if fk.shape[1:]!=(3,) or ck.shape[1:]!=(3,) or pk.shape[1:]!=(4,) or zk.shape[1:]!=(3,):
        raise ValueError('Invalid index shapes')
    if (result['forcing'].shape != (len(fk),n_r) or
        result['source'].shape != (len(ck),n_r) or
        result['Z_infinity'].shape != (len(ck),) or
        result['partial_I'].shape != (len(pk),)):
        raise ValueError('Inconsistent array shapes')
    for indices in (fk,ck,pk,zk):
        if len({tuple(map(int,row)) for row in indices})!=len(indices):
            raise ValueError('Duplicate saved keys')
    if set(map(tuple,ck))&set(map(tuple,zk)):
        raise ValueError('Open and zero sets overlap')
    if not all(np.all(np.isfinite(result[key])) for key in ('forcing','source','Z_infinity','partial_I')):
        raise ValueError('Nonfinite scientific data')
    if any(tuple(row[:3]) not in {tuple(v) for v in ck} for row in pk):
        raise ValueError('Partial Green overlap for absent channel')
    return result


def get_forcing(result: dict, l: int, m: int, g: int) -> tuple[np.ndarray,np.ndarray]:
    key=(l,m,g)
    for i,k in enumerate(result['forcing_keys']):
        if tuple(map(int,k))==key:
            return result['radius'],result['forcing'][i]
    raise KeyError(f'Fourier coefficient {key} was not computed for an open response')


def get_channel(result: dict, ell_f: int, m_f: int, g: int) -> dict:
    key=(ell_f,m_f,g)
    for i,k in enumerate(result['channel_keys']):
        if tuple(map(int,k))==key:
            partials={int(k4[3]):complex(value)
                      for k4,value in zip(result['partial_keys'],result['partial_I'])
                      if tuple(map(int,k4[:3]))==key}
            return {'r':result['radius'], 'S':result['source'][i],
                    'Z_infinity':complex(result['Z_infinity'][i]),'I_l':partials}
    raise KeyError(f'Channel {key} is not a computed open response')
