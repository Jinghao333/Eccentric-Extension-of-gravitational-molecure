"""Frozen Gaunt integral ∫Y_lm Y_ell_i,m_i conjugate(Y_ell_f,m_f)dΩ."""
from __future__ import annotations
import math
from functools import lru_cache
from dataclasses import dataclass
from sympy.physics.wigner import wigner_3j

@dataclass(frozen=True)
class SelectionReport:
    l:int
    m:int
    ell_i:int
    m_i:int
    ell_f:int
    m_f:int
    gaunt_value:float
    status:str

def candidate_forcing_degrees(ell_i,m_i,ell_f,m_f):
    m=m_f-m_i
    return tuple(l for l in range(abs(ell_f-ell_i),ell_f+ell_i+1)
                 if abs(m)<=l and (l+ell_i+ell_f)%2==0)

@lru_cache(maxsize=None)
def gaunt_coupling(l,m,ell_i,m_i,ell_f,m_f):
    if m+m_i!=m_f or l<0 or abs(m)>l or abs(m_i)>ell_i or abs(m_f)>ell_f or abs(l-ell_i)>ell_f or ell_f>l+ell_i or (l+ell_i+ell_f)%2:
        return 0.0
    pref=(-1)**m_f*math.sqrt((2*l+1)*(2*ell_i+1)*(2*ell_f+1)/(4*math.pi))
    return float(pref*wigner_3j(l,ell_i,ell_f,0,0,0)*wigner_3j(l,ell_i,ell_f,m,m_i,-m_f))

def selection_status(l,m,ell_i,m_i,ell_f,m_f):
    value=gaunt_coupling(l,m,ell_i,m_i,ell_f,m_f)
    return SelectionReport(l,m,ell_i,m_i,ell_f,m_f,value,'NONZERO_COUPLING' if value!=0.0 else 'ZERO_COUPLING')

def allowed_forcing_degrees(*,ell_i,m_i,ell_f,m_f):
    m=m_f-m_i
    return tuple(l for l in candidate_forcing_degrees(ell_i,m_i,ell_f,m_f)
                 if gaunt_coupling(l,m,ell_i,m_i,ell_f,m_f)!=0.0)
