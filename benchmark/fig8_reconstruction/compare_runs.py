#!/usr/bin/env python3
"""One-point optimized/reference equivalence check from minimal records.

Requires matching physical inputs, truncations and grid; only engine choices
may differ. Algebraically reconstructs all fluxes without solving again.
This is NOT a radial/grid/harmonic convergence test.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
from derive_quantities import derive


def compare(optimized: Path, reference: Path, *, rtol: float=1e-4,
            atol: float=0.0) -> dict:
    a=derive(optimized); b=derive(reference)
    ac=dict(a['run']['config']);bc=dict(b['run']['config'])
    for c in (ac,bc):
        c.pop('forcing_engine',None);c.pop('radial_engine',None)
    if any(a['run'].get('source_sha256',{}).get(n+'.py') != b['run'].get('source_sha256',{}).get(n+'.py') for n in ('01_atom','02_kepler_orbit','03_tidal_field','04_fourier_forcing','05_angular_coupling','06_radial_source','07_green_response','08_flux','09_response_channel_catalogue')):
        raise ValueError('Original physics module hashes differ; not a pure scheduling regression')
    if ac!=bc or len(a['raw']['radius'])!=len(b['raw']['radius']) or not all(a['raw']['radius']==b['raw']['radius']):
        raise ValueError('Physics, numerical cutoffs and radial grids must match')
    aa={(r['ell_f'],r['m_f'],r['g']):r for r in a['channels']}
    bb={(r['ell_f'],r['m_f'],r['g']):r for r in b['channels']}
    if aa.keys()!=bb.keys():
        raise ValueError('Channel catalogues disagree')
    failures=[]
    reference_total=b['summary']['normalized_flux']
    for key in sorted(aa):
        left,right=aa[key],bb[key]
        if left['status']!=right['status']:
            failures.append(dict(channel=key,issue='different status'))
        elif abs(left['normalized_flux']-right['normalized_flux']) > atol+rtol*max(
            abs(right['normalized_flux']),abs(reference_total)*1e-10):
            failures.append(dict(channel=key,issue='flux discrepancy',
                 optimized=left['normalized_flux'],reference=right['normalized_flux']))
    diff=abs(a['summary']['normalized_flux']-reference_total)
    return dict(number_channels=len(aa),channel_failures=failures,
        reference_flux=reference_total,optimized_flux=a['summary']['normalized_flux'],
        normalized_absolute_difference=diff,
        total_pass=diff<=atol+rtol*abs(reference_total),
        note='Method regression only; convergence deferred')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('optimized',type=Path);p.add_argument('reference',type=Path)
    p.add_argument('--rtol',type=float,default=1e-4)
    p.add_argument('--atol',type=float,default=0.0)
    p.add_argument('--output',type=Path,help='Optional analysis report (not canonical run data)')
    a=p.parse_args()
    result=compare(a.optimized,a.reference,rtol=a.rtol,atol=a.atol)
    print(json.dumps(result,indent=2))
    if a.output:a.output.write_text(json.dumps(result,indent=2)+'\n')
    if result['channel_failures'] or not result['total_pass']:raise SystemExit(1)


if __name__=='__main__':main()
