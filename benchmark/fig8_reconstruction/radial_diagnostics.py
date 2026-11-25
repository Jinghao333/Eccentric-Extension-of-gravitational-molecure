#!/usr/bin/env python3
"""On-demand radial response reconstruction from NONREDUNDANT primitives.

Physics provenance: r and u_lmg(r) come from physics.h5, while R_i(r),
Gaunt(l,m;li,mi;lf,mf) and the regular Coulomb Z_-(r) are recomputed with the
original hashed modules. NO Fourier quadrature or production Green overlap is
repeated. Cumulative grid integrals are qualitative diagnostics, NOT the
original adaptive overlaps.

By default no new file is created. --export-file saves an explicitly optional
analysis-only HDF5 view for chosen top channels.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import h5py
from scipy.integrate import cumulative_trapezoid
from one_point_fast import load_modules, file_sha256, MODULE_NAMES
from derive_quantities import derive


def profiles(run_directory:Path,module_dir:Path,*,top:int=5):
    if top<1:raise ValueError('top must be positive')
    result=derive(run_directory)
    cfg=result['run']['config']; raw=result['raw']
    # Refuse silent source/convention drift when recreating spatial sources.
    old_hashes=result['run'].get('source_sha256',{})
    for name in MODULE_NAMES:
        path=Path(module_dir)/(name+'.py')
        previous=old_hashes.get(name+'.py')
        if previous and (not path.is_file() or file_sha256(path)!=previous):
            raise ValueError(f'Physics module changed since the run: {name}.py')
    modules=load_modules(Path(module_dir))
    atom=modules['01_atom']; source=modules['06_radial_source']; response=modules['07_green_response']
    angular=modules['05_angular_coupling']
    state=atom.HydrogenicState(
        atom.AtomicParameters(cfg['M'],cfg['alpha']/cfg['M'],cfg['M']*cfg['M_b_over_M']),
        atom.QuantumNumbers.from_principal(cfg['N_i'],cfg['ell_i'],cfg['m_i']))
    r=np.asarray(raw['radius'],dtype=float)
    mu=cfg['alpha']/cfg['M']; omega=result['summary']['omega_orb']; Ei=result['summary']['E_i']
    common=2*mu**2*r*np.asarray(state.radial(r),dtype=float)
    forcing=result['forcing']
    backend=response.MpmathWhittakerBackend(cfg['mpmath_dps'])
    zcache={}
    ranked=sorted((row for row in result['channels'] if row['status']=='open'),
                  key=lambda x:abs(x['normalized_flux']),reverse=True)
    for row in ranked[:top]:
        ell,mf,g=row['ell_f'],row['m_f'],row['g'];m=mf-cfg['m_i']
        members=[p for p in result['partials'] if
                 (p['ell_f'],p['m_f'],p['g'])==(ell,mf,g)]
        ls=[p['forcing_l'] for p in members]
        partial=np.stack([common*forcing[(l,m,g)]*
            angular.gaunt_coupling(l,m,cfg['ell_i'],cfg['m_i'],ell,mf)
            for l in ls])
        source_total=np.sum(partial,axis=0)
        saved_source=result['source'][(ell,mf,g)]
        if not np.allclose(saved_source,source_total,rtol=1e-11,atol=1e-15):
            raise ValueError(f'Saved source differs from reconstructed Fourier/Gaunt source: {(ell,mf,g)}')
        label=(ell,g)
        if label not in zcache:
            descriptor=source.RadialOperatorDescriptor(cfg['M'],mu,ell,Ei,g,omega)
            params=response.WhittakerParameters.from_operator(descriptor)
            basis=response.CoulombHomogeneousBasis(params,backend)
            zcache[label]=np.asarray(basis.inner(r),dtype=np.complex128)
        z=zcache[label]
        integrands=partial*z[None,:]
        cumulative=cumulative_trapezoid(integrands,r,axis=1,initial=0.0)
        prodI=complex(row['I_re'],row['I_im'])
        yield dict(label=(ell,mf,g),r=r,forcing_l=np.asarray(ls),
                   source_partials=partial,source_total=source_total,z_minus=z,
                   overlap_integrands=integrands,cumulative_I_l_grid=cumulative,
                   cumulative_I_total_grid=cumulative.sum(axis=0),production_I=prodI)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('run_directory',type=Path)
    p.add_argument('--module-dir',type=Path,default=Path(__file__).resolve().parent/'physics')
    p.add_argument('--top',type=int,default=5)
    p.add_argument('--export-file',type=Path,help='OPTIONAL redundant view; default is no file')
    a=p.parse_args()
    data=list(profiles(a.run_directory,a.module_dir,top=a.top))
    for d in data:
        print(d['label'],'production I=',d['production_I'],
              'grid-trapezoid I=',d['cumulative_I_total_grid'][-1])
    if a.export_file:
        a.export_file.parent.mkdir(parents=True,exist_ok=True)
        with h5py.File(a.export_file,'w') as h5:
            h5.attrs['description']='Optional derived diagnostic view; not canonical numerical output'
            if data:
                h5.create_dataset('radius',data=data[0]['r'])
            for d in data:
                ell,mf,g=d['label']
                grp=h5.create_group(f'ell{ell}_m{mf}_g{g}')
                for k in ('forcing_l','source_partials','z_minus'):
                    grp.create_dataset(k,data=d[k],compression='gzip')
        print('Optional exported radial profile:',a.export_file)

if __name__=='__main__':main()
