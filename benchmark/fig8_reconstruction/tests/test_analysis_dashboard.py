"""Fast tests of analysis algebra and read-only dashboard endpoints.

The fixture has manufactured forcing/amplitude, NOT a Xu Fig.8 physics result.
"""
from __future__ import annotations
import json
import math
from pathlib import Path
import sys
from types import SimpleNamespace
from threading import Thread
from http.server import ThreadingHTTPServer
from urllib.request import urlopen

import numpy as np
import pytest
from scipy.integrate import simpson

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from primitive_store import write_physics, write_config_first,radial_grid
from derive_quantities import _eval_wronskian
from analysis_dashboard.core import (radial_field,orbital_parameters,
    shell_and_scales,channel_atlas,channel_detail,load_run,gaunt)
from analysis_dashboard.server import create_handler

CONFIG=json.loads((Path(__file__).resolve().parents[1]/'configs/fig8_point.json').read_text())


def test_production_field_norm_and_cloud_scales():
    cfg=CONFIG
    r=np.linspace(0,1000,9000)
    mu=cfg['alpha']/cfg['M']; Mb=cfg['M_b_over_M']*cfg['M']
    mass=simpson(2*mu*mu*r*r*radial_field(r,cfg)**2,x=r)
    assert mass==pytest.approx(Mb,rel=3e-6)
    scales=shell_and_scales(cfg)
    assert scales['density_peak']==pytest.approx(50)
    assert scales['shell_peak']==pytest.approx(100)
    assert scales['mean_radius']==pytest.approx(125)
    assert scales['r50']==pytest.approx(116.772722,abs=.0001)
    assert scales['r90']==pytest.approx(199.83974,abs=.0001)


def test_orbital_period_and_closed_threshold_labels():
    p=orbital_parameters(CONFIG)
    assert p['period']==pytest.approx(2885.465896816684)
    assert p['periapsis']==pytest.approx(50/1.4)
    assert p['apoapsis']==pytest.approx(50/.6)
    assert p['E_i']==pytest.approx(-.001)


@pytest.fixture
def small_manufactured_run(tmp_path):
    cfg=dict(CONFIG)
    cfg.update(point_name='SYNTHETIC_TEST_NOT_PHYSICAL',
               ell_f_max=0,g_values=[1],log_points=20,orbit_interval_points=0,
               r_min=.1,r_max=350.0)
    r=radial_grid(cfg)
    coef=(3.2e-4+1.2e-4j)*np.exp(-r/180.)
    C=gaunt(1,-1,1,1,0,0)
    mu=cfg['alpha']/cfg['M']; M=cfg['M'];Mb=cfg['M_b_over_M']*M
    source=2*mu**2*r*radial_field(r,cfg)*C*coef
    # Deliberately arbitrary I_l; its existence does not imply numerical
    # equivalence to the adaptive physical integral. This tests storage algebra.
    I=complex(1.2e-7,-.43e-7)
    E=-mu*cfg['alpha']**2/(2*cfg['N_i']**2)
    omega=((1-cfg['e']**2)/cfg['p'])**1.5/M
    delta=E+omega
    W=_eval_wronskian(mu,M,0,delta,cfg['mpmath_dps'])
    eta=M*mu**2/math.sqrt(2*mu*delta)
    amp=math.exp(math.pi*eta/2)*I/W
    forcing=SimpleNamespace(coefficient=coef,quadrature_converged=np.ones(len(r),bool))
    folder=tmp_path/'manufactured'
    write_config_first(folder,cfg)
    write_physics(folder/'physics.h5',config=cfg,
                  forcing_cache={(1,-1,1):forcing},
                  source_cache={(0,0,1):source},amplitude_cache={(0,0,1):amp},
                  sampled_zero_keys=[],partial_records=[dict(ell_f=0,m_f=0,g=1,
                        forcing_l=1,I_l_re=I.real,I_l_im=I.imag)])
    return folder


def test_channel_algebra_and_profiles(small_manufactured_run):
    result=load_run(small_manufactured_run)
    atlas=channel_atlas(result)
    assert len(atlas['top'])==1
    assert atlas['top'][0]['fraction']==pytest.approx(1)
    assert atlas['top'][0]['ell']==0
    assert atlas['top'][0]['m']==0
    assert len(atlas['angular_modes'])==1
    assert atlas['angular_modes'][0]['ell']==0
    assert atlas['angular_modes'][0]['m']==0
    assert atlas['angular_modes'][0]['power']==pytest.approx(atlas['total_normalized_flux'])
    assert atlas['angular_modes'][0]['share']==pytest.approx(1)
    assert sum(row['share'] for row in atlas['angular_modes'])==pytest.approx(1)
    deep=channel_detail(result,0,0,1)
    assert deep['response']['source_relative_difference']<1e-12
    assert deep['response']['multipole_coherence']==pytest.approx(1)
    assert len(deep['radial']['overlap_total_abs'])==20
    assert deep['response']['A_out']['abs']>0
    assert deep['response']['grid_minus_production_I'] is not None
    # All details read-only; no diagnostic files written to run.
    assert sorted(p.name for p in small_manufactured_run.iterdir())==['config.json','physics.h5']


def test_http_api_offline(small_manufactured_run):
    http=ThreadingHTTPServer(('127.0.0.1',0),create_handler(small_manufactured_run.parent))
    thread=Thread(target=http.serve_forever,daemon=True);thread.start()
    root=f'http://127.0.0.1:{http.server_port}'
    try:
        def fetch(path):
            with urlopen(root+path,timeout=18) as f:return json.load(f)
        assert 'manufactured' in fetch('/api/runs')['runs']
        o=fetch('/api/overview?run=manufactured')
        assert o['scales']['density_peak']==pytest.approx(50)
        assert o['totals']['orbital_power']>0
        a=fetch('/api/atlas?run=manufactured')
        assert len(a['top'])==1
        d=fetch('/api/channel?run=manufactured&ell=0&m=0&g=1')
        assert d['response']['multipole_coherence']==pytest.approx(1)
        geom=fetch('/api/geometry?run=manufactured&n=90')
        assert len(geom['rho'])==90
        assert len(geom['rho'][0])==90
    finally:
        http.shutdown();thread.join(timeout=3);http.server_close()
