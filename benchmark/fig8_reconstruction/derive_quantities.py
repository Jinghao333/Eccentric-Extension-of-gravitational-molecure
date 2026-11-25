#!/usr/bin/env python3
"""Reconstruct scientific observables from one config and an open-only HDF5.

Reads config.json + physics.h5. Does NOT redo orbit Fourier integration or radial
Green overlaps. Reconstructs Wronskians analytically and uses the stored complex
I_l to compute every Fig. 8 channel flux, all flux sums, interference, phases,
spectral properties, and orbital energy/charge algebraic balance.

No analysis products are written unless explicitly requested with --export-dir.
Derived exports are temporary views, NOT additional canonical run data.
"""
from __future__ import annotations

import argparse
import csv
from collections import defaultdict
import json
import math
from pathlib import Path
from typing import Any

import mpmath as mp
import numpy as np
from sympy.physics.wigner import wigner_3j

from primitive_store import load_primitives, SCHEMA


def _complex_pair(z: complex) -> dict[str, float]:
    return {'re': float(z.real), 'im': float(z.imag)}


def _sum(items) -> float:
    return math.fsum(float(v) for v in items)


def _ratio(n, d):
    return float(n/d) if d else None


def _eval_wronskian(mu: float, M: float, ell: int, delta: float, precision: int) -> complex:
    """Same canonical physical-r Whittaker Wronskian as 07_green_response.py.

    lambda=-2i sqrt(2mu Delta), kappa=2mu^2 M/lambda, and
    W=-lambda Gamma(2ell+2)/Gamma(ell+1-kappa).
    """
    if delta <= 0:
        raise ValueError('Outgoing Wronskian requires Delta>0')
    k_nr=math.sqrt(2*mu*delta)
    lam=complex(0.0,-2.0*k_nr)
    kappa=2.0*mu**2*M/lam
    with mp.workdps(precision):
        return complex(-mp.mpc(lam)*mp.gamma(2*ell+2)*mp.rgamma(ell+1-mp.mpc(kappa)))


def _coupled_l(ell_i: int, m_i: int, ell_f: int, m_f: int) -> tuple[int, ...]:
    """Exact allowed l for an equatorial Newtonian point source.

    Includes triangle, Gaunt, parity, and equatorial forcing support.
    """
    m=m_f-m_i
    result=[]
    for l in range(abs(ell_f-ell_i),ell_f+ell_i+1):
        if abs(m)>l or (l+ell_i+ell_f)%2 or (l+m)%2:
            continue
        if wigner_3j(l,ell_i,ell_f,0,0,0)==0:
            continue
        if wigner_3j(l,ell_i,ell_f,m,m_i,-m_f)==0:
            continue
        result.append(l)
    return tuple(result)


def derive(directory: Path) -> dict[str, Any]:
    """Return all algebraically reconstructible observables and diagnostics."""
    directory=Path(directory)
    data=load_primitives(directory)
    c=data['config']
    run={'config': c, 'output_schema':SCHEMA, 'completed':True,
         'source_sha256':data['provenance'].get('code_sha256',{}),
         'provenance':data['provenance'], 'quality':data['quality']}
    M=float(c['M']); mu=float(c['alpha'])/M
    mb=float(c['M_b_over_M']); q=float(c['q'])
    li=int(c['ell_i']); mi=int(c['m_i']); Ni=int(c['N_i'])
    orbit_e=float(c['e']); p=float(c['p'])
    Ei=-mu*float(c['alpha'])**2/(2*Ni**2)
    Omega=((1-orbit_e**2)/p)**1.5/M
    omega_i=mu+Ei
    precision=int(c['mpmath_dps'])
    divisor=q*q*mb
    g_values=[int(g) for g in c['g_values']]
    fkeys=[tuple(int(v) for v in row) for row in data['forcing_keys']]
    forcing={key: data['forcing'][j] for j,key in enumerate(fkeys)}
    channel_saved={tuple(map(int,k)):i for i,k in enumerate(data['channel_keys'])}
    zero_saved={tuple(map(int,k)) for k in data['zero_keys']}
    open_saved=set(channel_saved)
    if open_saved & zero_saved:
        raise ValueError('Computed open and sampled zero channels overlap')
    parts=defaultdict(list)
    rows=[]
    for key,value in zip(data['partial_keys'],data['partial_I']):
        ef,mf,g,l=(int(v) for v in key)
        z=complex(value)
        parts[(ef,mf,g)].append((l,z))
    part_rows=[]; pair_rows=[]; balance_rows=[]
    channel_rows=[]
    w_cache={}
    eligible=0
    expected_open=set()
    for ell in range(int(c['ell_f_max'])+1):
        for mf in range(-ell,ell+1):
            lm=_coupled_l(li,mi,ell,mf)
            m=mf-mi
            for g in g_values:
                delta=Ei+g*Omega
                structural=bool(lm) and (orbit_e!=0 or g==m)
                open_sector=delta>0
                threshold=delta==0
                is_eligible=bool(structural and open_sector)
                if is_eligible:
                    eligible+=1
                key=(ell,mf,g)
                if is_eligible:
                    expected_open.add(key)
                pvals=sorted(parts.get(key,[]))
                status=('open' if key in channel_saved else
                        'sampled_zero' if key in zero_saved else
                        'threshold' if structural and threshold else
                        'closed' if structural and delta<0 else 'excluded')
                if (pvals or key in channel_saved or key in zero_saved) and not is_eligible:
                    raise ValueError(f'Nonzero data exists for forbidden or closed channel {key}')
                if key in channel_saved:
                    missing_forcing=[(l,m,g) for l in lm if (l,m,g) not in forcing]
                    if missing_forcing:
                        raise ValueError(f'Missing cached forcing samples for {key}: {missing_forcing}')
                    if key in channel_saved and sorted(l for l,_ in pvals)!=list(lm):
                        raise ValueError(f'Partial-overlap catalogue incomplete or duplicated: {key}; expected {lm}')
                    if key in zero_saved and pvals:
                        raise ValueError(f'Zero channel has nonzero overlap metadata: {key}')
                I=sum((v for _,v in pvals),0j)
                k_nr=math.sqrt(2*mu*delta) if delta>0 else None
                omega_g=mu+delta
                k_kg=math.sqrt(omega_g**2-mu**2) if delta>0 else None
                eta=M*mu**2/k_nr if k_nr is not None else None
                W=None; A=0j
                if key in channel_saved:
                    wkey=(ell,g)
                    if wkey not in w_cache:
                        w_cache[wkey]=_eval_wronskian(mu,M,ell,delta,precision)
                    W=w_cache[wkey]
                    if W==0j:
                        raise ArithmeticError(f'Zero canonical Wronskian for {key}')
                    reconstructed=math.exp(math.pi*eta/2)*I/W
                    A=complex(data['Z_infinity'][channel_saved[key]])
                    if abs(A-reconstructed)>1e-10*max(abs(A),abs(reconstructed),1e-200):
                        raise ValueError(f'Saved Z_infinity is inconsistent with summed partial overlaps: {key}')
                power=abs(A)**2
                field=2*omega_g*k_kg*power if key in channel_saved else 0.0
                charge=-2*k_kg*power if key in channel_saved else 0.0
                orbital=2*g*Omega*k_kg*power if key in channel_saved else 0.0
                sum_abs=_sum(abs(z) for _,z in pvals)
                incoh=_sum(abs(z)**2 for _,z in pvals)
                coherence=_ratio(abs(I),sum_abs)
                enhancement=_ratio(abs(I)**2,incoh)
                for l,z in pvals:
                    part_rows.append(dict(ell_f=ell,m_f=mf,g=g,forcing_l=l,
                                          I_l_re=z.real,I_l_im=z.imag,
                                          abs_I_l=abs(z),phase_I_l=math.atan2(z.imag,z.real) if z else None))
                for j,(l,z) in enumerate(pvals):
                    for l2,z2 in pvals[j+1:]:
                        pair_rows.append(dict(ell_f=ell,m_f=mf,g=g,l1=l,l2=l2,
                                              interference=2*(z*z2.conjugate()).real))
                if key in channel_saved:
                    balance_rows.append(dict(ell_f=ell,m_f=mf,g=g,
                        orbital_minus_field_omega_charge=orbital-(field+omega_i*charge),
                        overlap_interference_identity=(abs(I)**2-incoh)-_sum(2*(z*z2.conjugate()).real for j,(_,z) in enumerate(pvals) for _,z2 in pvals[j+1:])))
                channel_rows.append(dict(ell_f=ell,m_f=mf,g=g,forcing_m=m,
                    status=status,delta_g=delta,omega_g=omega_g,
                    k_nr=k_nr,k_kg=k_kg,eta=eta,
                    allowed_forcing_l=lm,
                    I_re=I.real,I_im=I.imag,
                    W_re=W.real if W is not None else None,
                    W_im=W.imag if W is not None else None,
                    A_out_re=A.real,A_out_im=A.imag,
                    orbital_flux=orbital,normalized_flux=orbital/divisor,
                    field_energy_flux=field,charge_flux=charge,
                    relative_flux_fraction=None,
                    multipole_coherence=coherence,
                    multipole_enhancement=enhancement))
    if expected_open != (open_saved | zero_saved):
        missing=sorted(expected_open-(open_saved|zero_saved))
        unexpected=sorted((open_saved|zero_saved)-expected_open)
        raise ValueError(f'Incomplete or unexpected open channel coverage: missing={missing[:8]}, extra={unexpected[:8]}')
    total=_sum(r['normalized_flux'] for r in channel_rows)
    for row in channel_rows:
        row['relative_flux_fraction']=_ratio(row['normalized_flux'],total)
    by_dimension={}
    for name,dim in [('harmonics',('g',)),('angular',('ell_f',)),('azimuthal',('m_f',)),
                     ('ell_g',('ell_f','g')),('m_g',('m_f','g'))]:
        sums=defaultdict(list)
        for row in channel_rows:
            sums[tuple(row[k] for k in dim)].append(row['normalized_flux'])
        by_dimension[name]=[dict(**dict(zip(dim,key)),normalized_flux=_sum(vals),
                                 fraction=_ratio(_sum(vals),total))
                            for key,vals in sorted(sums.items())]
    cumul=0.0
    for row in by_dimension['harmonics']:
        cumul+=row['normalized_flux']
        row['cumulative_normalized_flux']=cumul
        row['cumulative_fraction']=_ratio(cumul,total)
    summary=dict(point=c.get('point_name','unnamed'),
                 alpha=c['alpha'],e=c['e'],p=c['p'],mu=mu,E_i=Ei,
                 omega_orb=Omega,r_peri=M*p/(1+orbit_e),
                 r_apo=M*p/(1-orbit_e),a=M*p/(1-orbit_e**2),
                 n_candidates=len(channel_rows),n_eligible=eligible,
                 n_open=sum(r['status']=='open' for r in channel_rows),
                 n_sampled_zero=sum(r['status']=='sampled_zero' for r in channel_rows),
                 normalized_flux=total,raw_orbital_flux=total*divisor,
                 n_unique_forcing_coefficients=len(fkeys),n_partial_overlaps=len(part_rows),
                 n_integral_groups=len({(int(v[0]),int(v[2])) for v in data['channel_keys']}),
                 warning='Derived from saved finite numerical primitives, not convergence-validated')
    return dict(run=run,raw=data,summary=summary,channels=channel_rows,
                partials=part_rows,interference_pairs=pair_rows,
                balance=balance_rows,**by_dimension,
                forcing=forcing,source={key:data['source'][i] for key,i in channel_saved.items()},
                Z_infinity={key:complex(data['Z_infinity'][i]) for key,i in channel_saved.items()})


def _csv(path: Path, rows: list[dict], columns: list[str] | None=None):
    if columns is None: columns=list(rows[0]) if rows else []
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=columns,extrasaction='ignore')
        w.writeheader()
        for row in rows:
            w.writerow({k:json.dumps(v) if isinstance(v,(tuple,list,dict)) else v
                        for k,v in row.items()})


def export_views(result: dict, folder: Path):
    """Optional, intentionally redundant analysis views, NOT canonical data."""
    folder=Path(folder); folder.mkdir(parents=True,exist_ok=True)
    for table in ('channels','partials','interference_pairs','balance',
                  'harmonics','angular','azimuthal','ell_g','m_g'):
        _csv(folder/(table+'.csv'),result[table])
    (folder/'summary.json').write_text(json.dumps(result['summary'],indent=2)+'\n',encoding='utf-8')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('run_directory',type=Path)
    p.add_argument('--export-dir',type=Path,help='OPTIONAL derived tables; not canonical run data')
    p.add_argument('--top',type=int,default=10)
    a=p.parse_args()
    result=derive(a.run_directory)
    print(json.dumps(result['summary'],indent=2))
    print('Top contributing channels (normalized flux):')
    top=sorted(result['channels'],key=lambda r:abs(r['normalized_flux']),reverse=True)
    for row in top[:a.top]:
        print(f"(ell_f,m_f,g)=({row['ell_f']},{row['m_f']},{row['g']})  "
              f"F={row['normalized_flux']:.12e}  coherence={row['multipole_coherence']}")
    if a.export_dir:
        export_views(result,a.export_dir)
        print('Optional analysis views:',a.export_dir)


if __name__=='__main__':
    main()
