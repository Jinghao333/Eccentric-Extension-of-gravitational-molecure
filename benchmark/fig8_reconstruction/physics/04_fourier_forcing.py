"""Branch-aware Fourier coefficients of the Xu tidal field.

u^[M]_{lmg}=(2pi)^-1 ∫u_lm(E)(1-e cos E)exp(+ig(E-e sin E))dE.
u^[t]_{lmg}=exp(+igΩ t_p)u^[M]_{lmg}.
"""
from __future__ import annotations
from dataclasses import dataclass
from enum import StrEnum
import importlib
import numpy as np
from scipy.integrate import quad_vec

tidal=importlib.import_module('03_tidal_field')

class FourierBasis(StrEnum):
    MEAN_ANOMALY='mean_anomaly'
    COORDINATE_TIME='coordinate_time'

@dataclass(frozen=True)
class QuadratureSettings:
    epsabs:float=1e-10
    epsrel:float=1e-8
    limit:int=240

@dataclass(frozen=True)
class ForcingCoefficient:
    ell:int
    m:int
    g:int
    field_radius:float
    mass_ratio:float
    value:complex
    estimated_absolute_error:float
    evaluations:int
    quadrature_variable:str
    basis:str
    branch_structure:object
    converged:bool

@dataclass(frozen=True)
class ForcingSample:
    l:int
    m:int
    g:int
    radial_points:object
    mass_ratio:float
    coefficient:object
    estimated_absolute_error:object
    evaluations:object
    quadrature_variable:str
    basis:FourierBasis
    native_basis:FourierBasis
    conversion_phase:complex
    branch_structures:tuple
    quadrature_converged:object

def forcing_coefficient(ell,m,g,field_radius,orbit,mass_ratio,settings=None,*,quadrature_variable='eccentric_anomaly'):
    settings=settings or QuadratureSettings()
    structure=tidal.radial_branch_structure(orbit,float(field_radius))
    if orbit.eccentricity==0:
        u=tidal.tidal_multipole(orbit.state_from_eccentric_anomaly(0.0),field_radius,ell,m,mass_ratio)
        return ForcingCoefficient(ell,m,g,float(field_radius),mass_ratio,complex(u) if g==m else 0j,0.0,0,
                                  quadrature_variable,'mean_anomaly',structure,True)
    e=orbit.eccentricity
    if quadrature_variable not in ('eccentric_anomaly','mean_anomaly'):
        raise ValueError('unknown Fourier quadrature variable')
    if quadrature_variable=='eccentric_anomaly':
        split=sorted({float(x) for x in structure.eccentric_anomaly_crossings if 0<float(x)<2*np.pi})
        def integrand(E):
            state=orbit.state_from_eccentric_anomaly(E)
            Ma=E-e*np.sin(E)
            return tidal.tidal_multipole(state,field_radius,ell,m,mass_ratio)*np.exp(1j*g*Ma)*(1-e*np.cos(E))
    else:
        split=sorted({float(x) for x in structure.mean_anomaly_crossings if 0<float(x)<2*np.pi})
        def integrand(Ma):
            return tidal.tidal_multipole(orbit.state_from_mean_anomaly(Ma),field_radius,ell,m,mass_ratio)*np.exp(1j*g*Ma)
    total=0j;error=0.;evals=0;success=True
    for a,b in zip((0.,*split),(*split,2*np.pi)):
        if b<=a:continue
        value,err,info=quad_vec(integrand,a,b,epsabs=settings.epsabs,epsrel=settings.epsrel,
                                limit=settings.limit,full_output=True)
        total+=complex(value);error+=float(err);evals+=int(info.neval);success=success and bool(info.success)
    return ForcingCoefficient(ell,m,g,float(field_radius),mass_ratio,total/(2*np.pi),error/(2*np.pi),
                              evals,quadrature_variable,'mean_anomaly',structure,success)

@dataclass(frozen=True)
class StageIIIForcingAdapter:
    orbit:object
    mass_ratio:float
    output_basis:FourierBasis
    settings:QuadratureSettings|None=None
    available_degrees:frozenset[int]|None=None
    spectrum_convergence:object=None
    @property
    def omega_orb(self):return self.orbit.omega_orb
    @property
    def central_mass(self):return self.orbit.central_mass
    def sample(self,l,m,g,radial_points):
        radii=np.asarray(radial_points,dtype=float)
        if self.available_degrees is not None and l not in self.available_degrees:
            raise KeyError(f'missing forcing degree {l}')
        rs=[forcing_coefficient(l,m,g,float(r),self.orbit,self.mass_ratio,self.settings) for r in radii]
        phase=np.exp(1j*g*self.orbit.omega_orb*self.orbit.periapsis_time) if self.output_basis==FourierBasis.COORDINATE_TIME else 1.+0j
        return ForcingSample(l,m,g,radii,self.mass_ratio,
                             np.asarray([x.value for x in rs],complex)*phase,
                             np.asarray([x.estimated_absolute_error for x in rs],float),
                             np.asarray([x.evaluations for x in rs],int),
                             'eccentric_anomaly' if self.orbit.eccentricity else 'mean_anomaly',self.output_basis,
                             FourierBasis.MEAN_ANOMALY,complex(phase),tuple(x.branch_structure for x in rs),
                             np.asarray([x.converged for x in rs],bool))
