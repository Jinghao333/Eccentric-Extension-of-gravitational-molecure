"""Manufactured, analytic closure of the minimal two-file scientific record."""
from __future__ import annotations
import json
import math
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import pytest
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from primitive_store import write_physics, write_config_first, load_primitives, radial_grid, SCHEMA
from derive_quantities import derive, _eval_wronskian


def _fixture(tmp_path):
    c=json.loads((Path(__file__).resolve().parents[1]/'configs'/'fig8_point.json').read_text())
    c.update(ell_f_max=1,g_values=[1],log_points=4,orbit_interval_points=0)
    q=c['q'];r=radial_grid(c)
    s=lambda ell:SimpleNamespace(coefficient=q*np.ones(4,dtype=complex),
        estimated_absolute_error=np.zeros(4),evaluations=np.full(4,8),
        quadrature_converged=np.ones(4,dtype=bool))
    # (ell_f,m_f,g)=(1,1,1) receives l=0 and l=2 with m=0;
    # (ell_f,m_f,g)=(0,0,1) receives l=1 with m=-1.
    forcing={(0,0,1):s(0),(2,0,1):s(2),(1,-1,1):s(1),
             (1,1,1):s(1),(3,-1,1):s(3),(1,-2,1):s(1),(3,-2,1):s(3)}
    partial=[dict(ell_f=1,m_f=1,g=1,forcing_l=0,I_l_re=1.0,I_l_im=0.0),
             dict(ell_f=1,m_f=1,g=1,forcing_l=2,I_l_re=-.4,I_l_im=.0),
             dict(ell_f=0,m_f=0,g=1,forcing_l=1,I_l_re=.2,I_l_im=.1)]
    # Additional forcing keys required for all eligible modes of ell_f<=1:
    from derive_quantities import _coupled_l
    for ell in range(2):
        for mf in range(-ell,ell+1):
            for l in _coupled_l(1,1,ell,mf):
                forcing.setdefault((l,mf-1,1),s(l))
    groups=[dict(ell_f=0,g=1,n_active=1,radial_neval=16,
                 radial_error_norm=1e-10,elapsed_seconds=.01),
            dict(ell_f=1,g=1,n_active=1,radial_neval=16,
                 radial_error_norm=1e-10,elapsed_seconds=.01)]
    Ei=-.001; Omega=((1-.4**2)/50.)**1.5; delta=Ei+Omega; mu=.2
    def outgoing(ell,I):
        W=_eval_wronskian(mu,c['M'],ell,delta,c['mpmath_dps'])
        k=math.sqrt(2*mu*delta)
        return math.exp(math.pi*c['M']*mu**2/(2*k))*I/W
    write_config_first(tmp_path/'run',c)
    write_physics(tmp_path/'run'/'physics.h5',config=c,forcing_cache=forcing,
                  source_cache={(1,1,1):np.ones(len(r),dtype=complex),
                                (0,0,1):np.ones(len(r),dtype=complex)*.2},
                  amplitude_cache={(1,1,1):outgoing(1,.6),
                                   (0,0,1):outgoing(0,.2+.1j)},
                  sampled_zero_keys=[(1,-1,1)],partial_records=partial)
    return c


def test_primitive_closure_and_destructive_interference(tmp_path):
    c=_fixture(tmp_path)
    raw=load_primitives(tmp_path/'run')
    assert raw['radius'].shape==(4,)
    assert len(raw['partial_I'])==3
    assert 'groups' not in raw  # group work statistics are not saved
    result=derive(tmp_path/'run')
    target=next(row for row in result['channels'] if (row['ell_f'],row['m_f'],row['g'])==(1,1,1))
    assert math.isclose(target['I_re'],.6,rel_tol=1e-13)
    assert math.isclose(target['multipole_coherence'],.6/1.4,rel_tol=1e-13)
    assert math.isclose(target['multipole_enhancement'],(.6)**2/(1+.4**2),rel_tol=1e-13)
    cross=[r for r in result['interference_pairs'] if (r['ell_f'],r['m_f'],r['g'])==(1,1,1)]
    assert len(cross)==1
    assert math.isclose(cross[0]['interference'],-.8,rel_tol=1e-13)
    Ei=result['summary']['E_i']; Omega=result['summary']['omega_orb'];mu=result['summary']['mu']
    delta=Ei+Omega;k=math.sqrt(2*mu*delta);kg=math.sqrt((mu+delta)**2-mu**2)
    W=_eval_wronskian(mu,c['M'],1,delta,c['mpmath_dps'])
    amplitude=math.exp(math.pi*c['M']*mu**2/(2*k))*.6/W
    expected=2*Omega*kg*abs(amplitude)**2/(c['q']**2*c['M_b_over_M'])
    assert math.isclose(target['normalized_flux'],expected,rel_tol=1e-12)
    assert all(abs(x['orbital_minus_field_omega_charge'])<1e-10 for x in result['balance'])


def test_missing_forcing_is_detected(tmp_path):
    _fixture(tmp_path)
    import h5py
    with h5py.File(tmp_path/'run'/'physics.h5','a') as h:
        h['forcing/keys_lmg'][0]=[11,11,11]
    with pytest.raises(ValueError,match='Missing cached forcing samples'):
        derive(tmp_path/'run')
