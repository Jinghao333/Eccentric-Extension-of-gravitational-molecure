"""Real bundled-physics checks, not mocks. No original repository needed."""
from __future__ import annotations
import importlib
import json
import math
from pathlib import Path
import subprocess
import sys
import numpy as np
from scipy.integrate import quad

ROOT=Path(__file__).resolve().parents[1]
PHYSICS=ROOT/'physics'
from one_point_fast import load_modules

def physics():
    return load_modules(PHYSICS)


def test_all_nine_bundled_modules_import_without_external_repository():
    for i,name in enumerate(('01_atom','02_kepler_orbit','03_tidal_field','04_fourier_forcing',
                              '05_angular_coupling','06_radial_source','07_green_response',
                              '08_flux','09_response_channel_catalogue')):
        module=physics()[name]
        assert Path(module.__file__).resolve().parent==PHYSICS


def test_atom_normalization_bound_energy_and_orbit_scales():
    atom=physics()['01_atom'];kepler=physics()['02_kepler_orbit']
    s=atom.HydrogenicState(atom.AtomicParameters(1.,.2,1e-4),atom.QuantumNumbers.from_principal(2,1,1))
    assert math.isclose(s.energy,-.001,rel_tol=1e-14)
    # Closed-form integral of the production field normalization.
    integral=quad(lambda r:2*.2**2*r*r*s.radial(r)**2,0,900,epsabs=1e-14)[0]
    assert math.isclose(integral,1e-4,rel_tol=1e-8)
    orbit=kepler.KeplerOrbit(1,50,.4)
    assert math.isclose(orbit.periapsis_radius,50/1.4)
    assert math.isclose(orbit.apoapsis_radius,50/.6)
    assert orbit.omega_orb>0
    E=np.array([.25,2.0,3.9,5.6]);st=orbit.state_from_eccentric_anomaly(E)
    assert np.allclose(st.radius,orbit.semi_major_axis*(1-.4*np.cos(E)))
    assert np.allclose(st.x**2+st.y**2,st.radius**2)
    # Second half is NOT folded into first half.
    assert st.true_anomaly[-1]>math.pi


def test_exact_angular_and_tidal_symmetries():
    ang=physics()['05_angular_coupling']
    tidal=physics()['03_tidal_field']
    orbit=physics()['02_kepler_orbit'].KeplerOrbit(1,50,.4)
    # Real Condon-Shortley Gaunt selection and exact orthogonality.
    assert ang.gaunt_coupling(1,1,1,1,0,0)==0
    assert math.isclose(ang.gaunt_coupling(0,0,1,1,1,1),1/math.sqrt(4*math.pi),rel_tol=1e-12)
    state=orbit.state_from_eccentric_anomaly(.9)
    assert tidal.tidal_multipole(state,float(state.radius)*.75,1,1,.001)==0j
    assert tidal.tidal_multipole(state,float(state.radius)*1.2,1,1,.001)!=0j


def test_circular_fourier_lock_and_eccentric_second_half():
    forcing=physics()['04_fourier_forcing']
    kepler=physics()['02_kepler_orbit']
    circ=kepler.KeplerOrbit(1,50,0)
    u1=forcing.forcing_coefficient(2,2,2,75,circ,.001)
    u2=forcing.forcing_coefficient(2,2,3,75,circ,.001)
    assert u1.value!=0j and u2.value==0j
    orb=kepler.KeplerOrbit(1,50,.4)
    u=forcing.forcing_coefficient(2,1,1,55,orb,.001)
    v=forcing.forcing_coefficient(2,-1,-1,55,orb,.001)
    # Complex spherical harmonics satisfy u_(l,-m,-g)=(-1)^m conj(u_lmg).
    assert abs(v.value+u.value.conjugate())<1e-12
    assert u.converged and v.converged


def test_analytic_whittaker_wronskian_matches_numerical_derivative():
    response=physics()['07_green_response']
    op=physics()['06_radial_source'].RadialOperatorDescriptor(1.,.2,0,-.001,1,.002177528874665)
    p=response.WhittakerParameters.from_operator(op)
    basis=response.CoulombHomogeneousBasis(p,response.MpmathWhittakerBackend(28))
    r=62.;h=1e-3
    m0=basis.inner(r); w0=basis.outer(r)
    dm=(basis.inner(r+h)-basis.inner(r-h))/(2*h)
    dw=(basis.outer(r+h)-basis.outer(r-h))/(2*h)
    # Wronskian convention W=M*W' - M'*W (physical r).
    measured=m0*dw-dm*w0
    assert abs(measured-basis.analytic_wronskian())/abs(measured)<2e-7


def test_real_end_to_end_optimized_vs_reference(tmp_path):
    cfg=json.loads((ROOT/'configs/smoke_real.json').read_text())
    cfg.update(log_points=15,orbit_interval_points=6,r_max=220,mpmath_dps=22)
    base=tmp_path/'base.json';base.write_text(json.dumps(cfg))
    runs=[]
    for fengine,rengine in [('gauss_batched_fallback','grouped_quad_vec'),('adaptive_cached','reference_per_channel')]:
        folder=tmp_path/(fengine+'_'+rengine)
        proc=subprocess.run([sys.executable,str(ROOT/'one_point_fast.py'),
              '--config',str(base),'--forcing',fengine,'--radial',rengine,
              '--output',str(folder)],cwd=ROOT,text=True,capture_output=True,timeout=90)
        assert proc.returncode==0,proc.stderr+'\n'+proc.stdout
        assert sorted(p.name for p in folder.iterdir())==['config.json','physics.h5']
        runs.append(folder)
    sys.path.insert(0,str(ROOT))
    from derive_quantities import derive
    from analysis_dashboard.core import load_run,run_overview,channel_atlas,channel_detail
    results=[derive(path) for path in runs]
    a,b=(x['summary']['normalized_flux'] for x in results)
    assert a>0 and math.isclose(a,b,rel_tol=5e-5)
    assert results[0]['summary']['n_open']==1
    dashboard_run=load_run(runs[0])
    assert run_overview(dashboard_run)['totals']
    assert channel_atlas(dashboard_run)['top']
    detail=channel_detail(dashboard_run,0,0,1)
    assert detail['response']
