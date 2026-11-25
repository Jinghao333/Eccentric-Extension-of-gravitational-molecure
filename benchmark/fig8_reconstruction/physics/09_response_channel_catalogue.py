"""Pure structural channel catalogue, before any expensive Fourier calculations."""
from __future__ import annotations
from dataclasses import dataclass
import importlib
import math
angular=importlib.import_module('05_angular_coupling')

@dataclass(frozen=True)
class ChannelLabel:
    ell_f:int
    m_f:int
    g:int

@dataclass(frozen=True)
class ResponseChannelMetadata:
    label:ChannelLabel
    forcing_m:int
    delta_g:float
    spectral_sector:str
    equatorial_parity_allowed:bool
    allowed_l_values:tuple
    equatorial_l_values:tuple
    gaunt_values:tuple
    eligible_for_flux_computation:bool
    rejection_reasons:tuple
    def source_kwargs(self):return dict(ell_f=self.label.ell_f,m_f=self.label.m_f,g=self.label.g)

@dataclass(frozen=True)
class ChannelCatalogue:
    records:tuple
    def flux_records(self):return tuple(r for r in self.records if r.eligible_for_flux_computation)
    def summary(self):
        return dict(candidate=len(self.records),flux_inputs=len(self.flux_records()),
                    closed=sum(r.spectral_sector=='CLOSED' for r in self.records),
                    threshold=sum(r.spectral_sector=='THRESHOLD' for r in self.records),
                    equatorial_allowed=sum(bool(r.equatorial_l_values) for r in self.records))

def build_response_channel_catalogue(*,point_name,initial_state,orbit,ell_f_max,g_values,forcing_l_max=None,equatorial=True):
    li=initial_state.quantum_numbers.ell;mi=initial_state.quantum_numbers.m
    Ei=initial_state.energy
    records=[]
    for lf in range(ell_f_max+1):
        for mf in range(-lf,lf+1):
            m=mf-mi
            allowed=angular.allowed_forcing_degrees(ell_i=li,m_i=mi,ell_f=lf,m_f=mf)
            eq_l=tuple(l for l in allowed if not equatorial or (l+m)%2==0)
            # Y_lm(theta=pi/2)=0 when l+m odd.
            parity=bool(eq_l)
            gaunts=tuple((l,angular.gaunt_coupling(l,m,li,mi,lf,mf)) for l in eq_l)
            missing=forcing_l_max is not None and any(l>forcing_l_max for l in eq_l)
            for g in g_values:
                delta=Ei+g*orbit.omega_orb
                sector='OPEN' if delta>0 else 'CLOSED' if delta<0 else 'THRESHOLD'
                temporal=(orbit.eccentricity!=0 or g==m)
                eligible=parity and temporal and not missing and sector=='OPEN'
                why=[]
                if not parity:why.append('ANGULAR_OR_EQUATORIAL_ZERO')
                if not temporal:why.append('CIRCULAR_LOCKING')
                if missing:why.append('MISSING_FORCING_DEGREE')
                if sector!='OPEN':why.append(sector)
                records.append(ResponseChannelMetadata(ChannelLabel(lf,mf,int(g)),m,delta,sector,
                             parity,allowed,eq_l,gaunts,eligible,tuple(why)))
    return ChannelCatalogue(tuple(records))
