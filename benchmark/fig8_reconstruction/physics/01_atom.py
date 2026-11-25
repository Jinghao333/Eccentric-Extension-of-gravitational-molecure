"""Xu Appendix A gravitational atom; field normalization matches source paper.

The complex scalar radial field R has integral r²|R|² dr = M_b/(2 mu²).
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from scipy.special import eval_genlaguerre, gammaln
@dataclass(frozen=True)
class AtomicParameters:
    M: float
    mu: float
    M_b: float
    def __post_init__(self):
        if not all(np.isfinite(x) and x>0 for x in (self.M,self.mu,self.M_b)):
            raise ValueError('positive finite masses required')
    @property
    def alpha(self): return self.M*self.mu

@dataclass(frozen=True)
class QuantumNumbers:
    n_r: int
    ell: int
    m: int
    def __post_init__(self):
        if self.n_r<0 or self.ell<0 or abs(self.m)>self.ell: raise ValueError('invalid angular/radial quantum numbers')
    @property
    def N(self): return self.n_r+self.ell+1
    @classmethod
    def from_principal(cls,N,ell,m):
        if N<ell+1:raise ValueError('N must be >= ell+1')
        return cls(int(N-ell-1),int(ell),int(m))

def gravitational_coupling(M,mu): return M*mu

def newtonian_bound_energy(mu,alpha,N): return -mu*alpha**2/(2*N**2)

def unit_radial(r,*,central_mass,scalar_mass,quantum_numbers):
    """Unit cloud mass R(r), including the 1/mu field normalization."""
    q=quantum_numbers; nr,ell,N=q.n_r,q.ell,q.N
    mu=float(scalar_mass); alpha=central_mass*mu
    r=np.asarray(r,dtype=float)
    if np.any(r<0): raise ValueError('negative radius')
    beta=2*mu*alpha/N
    rho=beta*r
    logpref=(-.5*np.log(2)-np.log(mu)+.5*(3*np.log(beta)+gammaln(nr+1)-np.log(2*N)-gammaln(N+ell+1)))
    result=np.exp(logpref)*rho**ell*eval_genlaguerre(nr,2*ell+1,rho)*np.exp(-rho/2)
    return float(result) if result.ndim==0 else result

@dataclass(frozen=True)
class HydrogenicState:
    parameters: AtomicParameters
    quantum_numbers: QuantumNumbers
    @property
    def N(self):return self.quantum_numbers.N
    @property
    def alpha(self):return self.parameters.alpha
    @property
    def energy(self):return newtonian_bound_energy(self.parameters.mu,self.alpha,self.N)
    def unit_radial(self,r):return unit_radial(r,central_mass=self.parameters.M,scalar_mass=self.parameters.mu,quantum_numbers=self.quantum_numbers)
    def radial(self,r):return np.asarray(self.unit_radial(r))*np.sqrt(self.parameters.M_b)
