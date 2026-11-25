"""Whittaker Coulomb Green response with analytic Wronskian and retarded branch."""
from __future__ import annotations
from dataclasses import dataclass
from enum import StrEnum
from functools import lru_cache
import math
import mpmath as mp
import numpy as np
from scipy.integrate import quad_vec
from scipy.interpolate import CubicSpline
source_mod=__import__('importlib').import_module('06_radial_source')

class StageVResponseStatus(StrEnum):
    COMPLETE='complete'
    INCOMPLETE_SOURCE='incomplete_source'
    THRESHOLD_SPECIAL_TREATMENT='threshold_special_treatment'
    EXACT_BOUND_POLE_NONINVERTIBLE='exact_bound_pole_noninvertible'

class ThresholdSpecialTreatmentError(ValueError):pass
class ExactBoundPoleError(ArithmeticError):pass

@dataclass(frozen=True)
class WhittakerParameters:
    central_mass:float
    scalar_mass:float
    ell_f:int
    energy:float
    sector:object
    lambda_g:complex|None
    kappa_g:complex|None
    whittaker_nu:float
    wave_number:float|None
    @classmethod
    def from_values(cls,*,central_mass,scalar_mass,ell_f,energy):
        sector=source_mod.classify_driven_energy(energy)
        if energy==0:return cls(central_mass,scalar_mass,ell_f,energy,sector,None,None,ell_f+.5,None)
        if energy<0:
            lam=complex(2*math.sqrt(-2*scalar_mass*energy));k=None
        else:
            k=math.sqrt(2*scalar_mass*energy)
            lam=complex(0,-2*k)
        return cls(central_mass,scalar_mass,ell_f,energy,sector,lam,
                   complex(2*scalar_mass**2*central_mass/lam),ell_f+.5,k)
    @classmethod
    def from_operator(cls,op):
        return cls.from_values(central_mass=op.central_mass,scalar_mass=op.scalar_mass,
                               ell_f=op.ell_f,energy=op.energy)
    def require_available(self):
        if self.energy==0:raise ThresholdSpecialTreatmentError('Zero-energy continuum requires a dedicated solver')
    def dimensionless_coordinate(self,r):
        self.require_available();return self.lambda_g*np.asarray(r)
    def exact_pole_principal_number(self):
        if self.energy>=0:return None
        N=int(round(abs(self.kappa_g)))
        if N<self.ell_f+1:return None
        expected=-self.scalar_mass*(self.central_mass*self.scalar_mass)**2/(2*N**2)
        return N if self.energy==expected else None

@dataclass(frozen=True)
class MpmathWhittakerBackend:
    precision_dps:int=35
    def __post_init__(self):
        if self.precision_dps<20:raise ValueError('precision >=20 needed')
    @staticmethod
    @lru_cache(maxsize=50000)
    def _scalar(name,kappa,nu,z,prec):
        with mp.workdps(prec):
            kp=mp.mpc(kappa); nu=mp.mpf(nu);z=mp.mpc(z)
            a=nu-kp+mp.mpf('.5'); b=1+2*nu
            if name=='M':
                value=mp.exp(-z/2)*mp.power(z,nu+mp.mpf('.5'))*mp.hyp1f1(a,b,z)
            else:
                value=mp.exp(-z/2)*mp.power(z,nu+mp.mpf('.5'))*mp.hyperu(a,b,z)
            return complex(value)
    def _eval(self,name,kappa,nu,z):
        a=np.asarray(z,dtype=np.complex128)
        result=np.empty(a.shape,dtype=np.complex128)
        for i in np.ndindex(a.shape):
            result[i]=self._scalar(name,complex(kappa),float(nu),complex(a[i]),self.precision_dps)
        return complex(result) if result.ndim==0 else result
    def M(self,kappa,nu,z):return self._eval('M',kappa,nu,z)
    def W(self,kappa,nu,z):return self._eval('W',kappa,nu,z)

@dataclass(frozen=True)
class CoulombHomogeneousBasis:
    parameters:WhittakerParameters
    backend:MpmathWhittakerBackend
    inner_normalization:complex=1+0j
    outer_normalization:complex=1+0j
    @property
    def canonical(self):return self.inner_normalization==1+0j and self.outer_normalization==1+0j
    def inner(self,r):
        p=self.parameters
        return self.inner_normalization*self.backend.M(p.kappa_g,p.whittaker_nu,p.dimensionless_coordinate(r))
    def outer(self,r):
        p=self.parameters
        return self.outer_normalization*self.backend.W(p.kappa_g,p.whittaker_nu,p.dimensionless_coordinate(r))
    def analytic_wronskian(self):
        p=self.parameters;p.require_available()
        if p.exact_pole_principal_number() is not None:return 0j
        with mp.workdps(self.backend.precision_dps):
            W=-mp.mpc(self.inner_normalization)*mp.mpc(self.outer_normalization)*mp.mpc(p.lambda_g)*mp.gamma(2*p.ell_f+2)*mp.rgamma(p.ell_f+1-mp.mpc(p.kappa_g))
            return complex(W)
    def require_invertible(self):
        self.parameters.require_available()
        if self.analytic_wronskian()==0j:
            raise ExactBoundPoleError('Coulomb resolvent at an exact bound pole')

@dataclass(frozen=True)
class ResponseIntegrationSettings:
    r_min:float
    r_max:float
    epsabs:float
    epsrel:float
    limit:int
    breakpoints:tuple=()

@dataclass(frozen=True)
class SampledRadialSource:
    radii:np.ndarray
    partial_values:np.ndarray
    l_values:np.ndarray
    inherited_channel:object
    @classmethod
    def from_stage_iv(cls,channel):
        return cls(np.asarray(channel.radial_points),np.asarray(channel.partial_sources),
                   np.asarray(channel.l_values),channel)

@dataclass(frozen=True)
class AsymptoticOverlapResult:
    partial_overlaps:np.ndarray
    total_overlap:complex
    quadrature_error_bound:float

@dataclass(frozen=True)
class AsymptoticResponseResult:
    status:StageVResponseStatus
    basis:object
    overlap:object
    wronskian:complex|None
    basis_asymptotic_coefficient:complex|None
    operator:object
    source:object
    @property
    def canonical_asymptotic_coefficient(self):
        return self.basis_asymptotic_coefficient if self.basis.canonical else None


def solve_asymptotic_response(*,operator,source,settings,backend,inner_normalization=1+0j,outer_normalization=1+0j):
    if settings.r_min<source.radii[0] or settings.r_max>source.radii[-1]:raise ValueError('source outside sample domain')
    basis=CoulombHomogeneousBasis(WhittakerParameters.from_operator(operator),backend,inner_normalization,outer_normalization)
    basis.require_invertible()
    W=basis.analytic_wronskian()
    overlaps=[];errtot=0.0
    for partial in source.partial_values:
        spline=CubicSpline(source.radii,partial)
        value,err,info=quad_vec(lambda r:basis.inner(r)*spline(r),settings.r_min,settings.r_max,
                                 epsabs=settings.epsabs,epsrel=settings.epsrel,limit=settings.limit,
                                 points=settings.breakpoints or None,full_output=True)
        if not info.success:raise ArithmeticError(f'Green overlap quadrature failed: {info.message}')
        overlaps.append(complex(value));errtot+=float(err)
    Ip=np.asarray(overlaps,complex);I=complex(Ip.sum())
    return AsymptoticResponseResult(StageVResponseStatus.COMPLETE,basis,AsymptoticOverlapResult(Ip,I,errtot),
                                    W,I/W,operator,source)
