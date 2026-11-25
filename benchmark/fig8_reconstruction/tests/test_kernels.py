"""Fast, independent tests of the optimization kernels (no production modules)."""
from __future__ import annotations
import math
from pathlib import Path
import sys
from dataclasses import dataclass
from types import SimpleNamespace

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from batched_quadrature import fourier_at_radius, radial_overlap_batch
from one_point_fast import build_grid, make_cached_provider, read_config


class SyntheticOrbit:
    eccentricity = 0.4
    omega_orb = 0.002178
    periapsis_time = 0.0
    periapsis_radius = 35.714285714285715
    apoapsis_radius = 83.33333333333333
    def state_from_eccentric_anomaly(self,E):
        return SimpleNamespace(E=np.asarray(E))


def splits(_orbit, _radius):
    return SimpleNamespace(eccentric_anomaly_crossings=(1.0,5.2))


def test_batched_fourier_harmonic_uses_full_2pi_orbit():
    orb = SyntheticOrbit()
    # A manufactured signal u_m(E)=exp[-i*m*(E-e sin E)]:
    # the exact Fourier coefficient is delta_{g,m}.
    def potential(state,r,l,m,q):
        E = state.E
        mean=E-orb.eccentricity*np.sin(E)
        return q*np.exp(-1j*m*mean)
    batch = fourier_at_radius(orb,potential,splits,2,3,tuple(range(0,7)),50.0,0.001,
                              initial_order=16,max_order=128,epsabs=1e-12,epsrel=1e-11)
    target = np.zeros(7,dtype=complex)
    target[3] = 0.001
    assert batch.accepted
    assert np.max(abs(batch.values-target)) < 1e-11
    assert batch.evaluations > 0


def test_composite_radial_overlap_integrates_partial_modes_coherently():
    r=np.linspace(0.2,3.2,21)
    # CubicSpline interpolates these cubics exactly; inner=1.
    parts=np.stack([r**2+1j*r, 2*r**3, r*0+1],axis=0)
    result=radial_overlap_batch(r,parts,lambda x:1+0j,epsabs=1e-11,epsrel=1e-11,
                                breakpoints=(1.3,2.5))
    expected=np.array([(r[-1]**3-r[0]**3)/3 +1j*(r[-1]**2-r[0]**2)/2,
                      (r[-1]**4-r[0]**4)/2, r[-1]-r[0]])
    assert result.success
    assert np.allclose(result.partial_integrals,expected,rtol=1e-10,atol=1e-11)
    assert abs(np.sum(result.partial_integrals)-np.sum(expected))<1e-10


def test_grid_reproduces_original_reference_count():
    data=read_config(Path(__file__).resolve().parents[1]/'configs'/'fig8_point.json')
    r=build_grid(data,SyntheticOrbit())
    assert len(r)==240
    assert r[0]==data['r_min']
    assert r[-1]==data['r_max']


def test_adaptive_adapter_caches_identical_input():
    from dataclasses import dataclass
    calls=[]
    @dataclass(frozen=True,slots=True)
    class FakeAdapter:
        orbit: object
        mass_ratio: float
        output_basis: str
        settings: object = None
        available_degrees: frozenset | None = None
        spectrum_convergence: object = None
        def sample(self,l,m,g,r):
            calls.append((l,m,g))
            return SimpleNamespace(coefficient=np.asarray(r)*0+g,
                                   quadrature_converged=np.array([True]*len(r)))
    C=make_cached_provider(FakeAdapter,SimpleNamespace(),SimpleNamespace())
    p=C(orbit=SyntheticOrbit(),mass_ratio=.001,output_basis='coordinate_time',
        g_values=(1,2),engine='adaptive_cached',available_degrees=frozenset({1,2}))
    grid=np.linspace(1,3,9)
    a=p.sample(1,0,1,grid)
    assert a is p.sample(1,0,1,grid)
    assert p.sample(1,0,2,grid) is not a
    assert calls==[(1,0,1),(1,0,2)]
    assert p._stats['cache_hits']==1


def test_periapsis_tangency_is_not_duplicate_integration_interval():
    orb=SyntheticOrbit()
    def tangent(_o,_r):return SimpleNamespace(eccentric_anomaly_crossings=(0.0,))
    def potential(state,r,l,m,q):return np.ones_like(state.E,dtype=complex)*q
    batch=fourier_at_radius(orb,potential,tangent,0,0,(0,1),35.0,.001,
                            initial_order=16,max_order=128)
    assert batch.accepted
    assert np.isclose(batch.values[0],.001,rtol=1e-10)
    assert abs(batch.values[1])<1e-10
