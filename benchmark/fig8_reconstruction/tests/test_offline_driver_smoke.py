"""Offline wiring test for a complete output, using disposable fake physics.

This proves that the run orchestration, IO, grouping and cache contracts work.
It does NOT test physical agreement with the real Xu modules; use the provided
compare_runs.py on your original checkout for that regression.
"""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace
from enum import StrEnum
import math
import sys
import json

import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from one_point_fast import run, read_config, MODULE_NAMES
from primitive_store import load_primitives
from derive_quantities import derive, export_views, _eval_wronskian


@dataclass(frozen=True)
class Params:
    M:float; mu:float; M_b:float

@dataclass(frozen=True)
class Qnums:
    N:int;ell:int;m:int
    @classmethod
    def from_principal(cls,N,ell,m):return cls(N,ell,m)

class Atom:
    def __init__(self,parameters,quantum_numbers):
        self.parameters=parameters; self.quantum_numbers=quantum_numbers
        self.energy=-.001
    def unit_radial(self,r):return np.exp(-np.asarray(r)/100.0)
    def radial(self,r):return self.unit_radial(r)*math.sqrt(self.parameters.M_b)

class Orbit:
    def __init__(self,M,p,e,tp,varpi):
        self.central_mass=M; self.p=p;self.eccentricity=e
        self.periapsis_time=tp;self.periapsis_azimuth=varpi
        self.omega_orb=((1-e*e)/p)**1.5/M
        self.periapsis_radius=35.0
        self.apoapsis_radius=80.0
        self.semi_major_axis=58.0
    def state_from_eccentric_anomaly(self,E):
        E=np.asarray(E)
        v=E-self.eccentricity*np.sin(E)
        rad=58*(1-self.eccentricity*np.cos(E))
        return SimpleNamespace(E=E,mean_anomaly=v,true_anomaly=v,
                               radius=rad,azimuth=v,radial_velocity=np.zeros_like(E),
                               angular_velocity=np.ones_like(E)*self.omega_orb,
                               x=rad*np.cos(v),y=rad*np.sin(v))

class Basis(StrEnum):
    COORDINATE_TIME='coordinate_time'; MEAN_ANOMALY='mean_anomaly'

@dataclass
class ForceSample:
    l:int;m:int;g:int;radial_points:object;mass_ratio:float
    coefficient:object;estimated_absolute_error:object;evaluations:object
    quadrature_variable:str;basis:object;native_basis:object
    conversion_phase:complex;branch_structures:tuple;quadrature_converged:object

@dataclass(frozen=True)
class QuadSettings:
    epsabs:float;epsrel:float;limit:int

@dataclass(frozen=True)
class Adapter:
    orbit:object;mass_ratio:float;output_basis:object
    settings:object=None;available_degrees:object=None;spectrum_convergence:object=None
    def sample(self,l,m,g,r):
        r=np.asarray(r)
        return ForceSample(l,m,g,r,self.mass_ratio,
            self.mass_ratio*np.ones_like(r,dtype=complex),
            np.zeros_like(r),np.ones_like(r,dtype=int),'eccentric_anomaly',
            self.output_basis,Basis.MEAN_ANOMALY,1+0j,tuple([None]*len(r)),
            np.ones_like(r,dtype=bool))

@dataclass
class Label:
    ell_f:int;m_f:int;g:int

@dataclass
class Meta:
    label:Label
    def source_kwargs(self):return dict(ell_f=self.label.ell_f,m_f=self.label.m_f,g=self.label.g)
    @property
    def forcing_m(self):return self.label.m_f-1
    @property
    def delta_g(self):return .001
    @property
    def spectral_sector(self):return 'OPEN'
    @property
    def equatorial_parity_allowed(self):return True
    @property
    def allowed_l_values(self):return (1,)
    @property
    def equatorial_l_values(self):return (1,)
    @property
    def gaunt_values(self):return ((1,.5),)
    @property
    def eligible_for_flux_computation(self):return True
    @property
    def rejection_reasons(self):return ()

class Catalogue:
    def __init__(self,g_values):self.records=tuple(Meta(Label(0,0,g)) for g in g_values)
    def flux_records(self):return self.records
    def summary(self):return dict(candidate=len(self.records),flux_inputs=len(self.records))

class FakeStageIV:
    def __init__(self,radii,state,provider,ell_f,m_f,g):
        s=provider.sample(1,-1,g,radii)
        self.ell_f=ell_f;self.m_f=m_f;self.g=g;self.required_m=-1
        self.total_source=np.asarray(s.coefficient)*np.ones_like(radii)
        self.available_source_sum=self.total_source
        self.partial_sources=np.array([self.total_source])
        self.forcing_coefficients=np.array([s.coefficient])
        self.l_values=np.array([1]);self.missing_l_values=np.array([],dtype=int)
        self.gaunt_coefficients=np.array([.5]);self.available_absolute_error_bound=np.zeros_like(radii)
        self.inherited_forcing_metadata=(s,)
        self.completeness_status='complete';self.driven_energy=.001
        self.radial_operator=SimpleNamespace(ell_f=ell_f,g=g)

@dataclass
class RSettings:
    r_min:float;r_max:float;epsabs:float;epsrel:float;limit:int
    breakpoints:tuple=()

class MBackend:
    def __init__(self,dps):self.dps=dps

class Whittaker:
    @classmethod
    def from_operator(cls,operator):
        delta=-.001+operator.g*((1-.4**2)/50.)**1.5
        return SimpleNamespace(wave_number=math.sqrt(2*.2*delta),
                               ell_f=operator.ell_f,delta=delta)

class Cbasis:
    def __init__(self,p,b):self.parameters=p
    def require_invertible(self):return
    def analytic_wronskian(self):return _eval_wronskian(.2,1.,self.parameters.ell_f,self.parameters.delta,35)
    def inner(self,r):return 1+0j

class FMode(StrEnum):XU_KG_KINEMATIC='kg'
class FStatus(StrEnum):OPEN='open'

def fake_flux(amplitude,**kwargs):
    mu=.2;omega=kwargs['omega_orb'];delta=-.001+kwargs['g']*omega
    knr=math.sqrt(2*mu*delta);kkg=math.sqrt((mu+delta)**2-mu**2)
    power=abs(amplitude)**2
    return SimpleNamespace(status=FStatus.OPEN,
        orbital_energy_loss=2*kwargs['g']*omega*kkg*power,
        field_energy_flux=2*(mu+delta)*kkg*power,
        charge_flux=-2*kkg*power,k_nr=knr,k_kg=kkg,
        kinematic_factor_orbital=2*kwargs['g']*omega*kkg)


def test_driver_saves_only_nonredundant_primitives(monkeypatch,tmp_path):
    modules={
        '01_atom':SimpleNamespace(AtomicParameters=Params,QuantumNumbers=Qnums,HydrogenicState=Atom),
        '02_kepler_orbit':SimpleNamespace(KeplerOrbit=Orbit),
        '03_tidal_field':SimpleNamespace(),
        '04_fourier_forcing':SimpleNamespace(StageIIIForcingAdapter=Adapter,FourierBasis=Basis,
             QuadratureSettings=QuadSettings),
        '05_angular_coupling':SimpleNamespace(),
        '06_radial_source':SimpleNamespace(scalar_source_channel=lambda r,initial_state,forcing_provider,forcing_basis,ell_f,m_f,g:
            FakeStageIV(r,initial_state,forcing_provider,ell_f,m_f,g)),
        '07_green_response':SimpleNamespace(MpmathWhittakerBackend=MBackend,
             ResponseIntegrationSettings=RSettings,WhittakerParameters=Whittaker,
             CoulombHomogeneousBasis=Cbasis),
        '08_flux':SimpleNamespace(coulomb_eta=lambda M,mu,k:M*mu**2/k,
             FluxKinematicsMode=FMode,FluxChannelStatus=FStatus,
             channel_flux_from_amplitude=fake_flux),
        '09_response_channel_catalogue':SimpleNamespace(
             build_response_channel_catalogue=lambda **kwargs:Catalogue(kwargs['g_values'])),
    }
    mod_dir=tmp_path/'modules';mod_dir.mkdir()
    for name in MODULE_NAMES:
        (mod_dir/(name+'.py')).write_text('mock for orchestration test\n')
        monkeypatch.setitem(sys.modules,name,modules[name])
    c=read_config(Path(__file__).resolve().parents[1]/'configs'/'fig8_point.json')
    c.update(g_values=[1,2],ell_f_max=0,log_points=12,orbit_interval_points=0,
             forcing_engine='adaptive_cached',radial_engine='grouped_quad_vec')
    out=tmp_path/'result'
    result=run(c,mod_dir,out)
    assert result['number_computed_open']==2
    assert result['normalized_fig8_flux']>0
    # Two authoritative files only: everything else is reconstructed on demand.
    assert sorted(f.name for f in out.iterdir()) == ['config.json','physics.h5']
    assert json.loads((out/'config.json').read_text())['g_values']==[1,2]
    raw=load_primitives(out)
    assert raw['radius'].shape==(12,)
    assert raw['forcing_keys'].shape==(2,3)
    assert raw['partial_keys'].shape==(2,4)
    derived=derive(out)
    assert derived['summary']['n_open']==2
    assert math.isclose(derived['summary']['normalized_flux'],
                        result['normalized_fig8_flux'],rel_tol=1e-12)
    assert math.isclose(sum(row['normalized_flux'] for row in derived['harmonics']),
                        result['normalized_fig8_flux'],rel_tol=1e-12)
    assert all(abs(row['orbital_minus_field_omega_charge']) < 1e-12
               for row in derived['balance'])
    # Optional materialized views go elsewhere, NOT into the canonical run.
    optional=tmp_path/'diagnostics'
    export_views(derived,optional)
    assert (optional/'channels.csv').is_file()
    assert sorted(f.name for f in out.iterdir()) == ['config.json','physics.h5']
