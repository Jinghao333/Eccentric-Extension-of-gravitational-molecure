"""Xu outgoing scalar ionization flux, NR Coulomb response + KG momentum."""
from __future__ import annotations
from enum import StrEnum
from dataclasses import dataclass
import math

class FluxKinematicsMode(StrEnum):
    XU_KG_KINEMATIC='xu_kg_kinematic'
    NR='nr'
class FluxChannelStatus(StrEnum):
    OPEN='open'
    CLOSED='closed'
    THRESHOLD_UNRESOLVED='threshold_unresolved'
    INVALID_FREQUENCY_SECTOR='invalid_frequency_sector'

@dataclass(frozen=True)
class CarrierFrequencies:
    scalar_mass:float
    initial_energy:float
    g:int
    omega_orb:float
    omega:float
    delta_g:float
    omega_g:float

def carrier_frequencies(scalar_mass,initial_energy,g,omega_orb):
    return CarrierFrequencies(scalar_mass,initial_energy,g,omega_orb,
                              scalar_mass+initial_energy,initial_energy+g*omega_orb,
                              scalar_mass+initial_energy+g*omega_orb)

def nr_momentum(mu,delta_g):return math.sqrt(2*mu*delta_g)
def kg_momentum(mu,omega_g):return math.sqrt(omega_g**2-mu**2)
def coulomb_eta(M,mu,k_nr):return M*mu**2/k_nr

@dataclass(frozen=True)
class FluxChannelResult:
    status:FluxChannelStatus
    orbital_energy_loss:float|None
    field_energy_flux:float|None
    charge_flux:float|None
    k_nr:float|None
    k_kg:float|None
    kinematic_factor_orbital:float|None
    amplitude:complex|None

def channel_flux_from_amplitude(amplitude,*,scalar_mass,initial_energy,g,omega_orb,ell_f,m_f,
                                kinematics_mode,ell_i=None,m_i=None,**kwargs):
    f=carrier_frequencies(scalar_mass,initial_energy,g,omega_orb)
    if f.delta_g<0:return FluxChannelResult(FluxChannelStatus.CLOSED,0.,0.,0.,None,None,None,None)
    if f.delta_g==0:return FluxChannelResult(FluxChannelStatus.THRESHOLD_UNRESOLVED,None,None,None,None,None,None,None)
    if initial_energy>=0 or g<=0:
        return FluxChannelResult(FluxChannelStatus.INVALID_FREQUENCY_SECTOR,None,None,None,None,None,None,None)
    kn=nr_momentum(scalar_mass,f.delta_g);kg=kg_momentum(scalar_mass,f.omega_g)
    selected=kg if kinematics_mode==FluxKinematicsMode.XU_KG_KINEMATIC else kn
    power=abs(amplitude)**2
    field=2*f.omega_g*selected*power
    charge=-2*selected*power
    orbital=field+f.omega*charge
    return FluxChannelResult(FluxChannelStatus.OPEN,orbital,field,charge,kn,kg,2*g*omega_orb*selected,complex(amplitude))

def flux_channel(response,*,kinematics_mode):
    c=response.source.inherited_channel
    p=response.basis.parameters
    if response.status!='complete':raise ArithmeticError('incomplete response')
    eta=coulomb_eta(p.central_mass,p.scalar_mass,p.wave_number)
    amp=math.exp(math.pi*eta/2)*response.basis_asymptotic_coefficient*response.basis.outer_normalization
    return channel_flux_from_amplitude(amp,scalar_mass=p.scalar_mass,initial_energy=response.operator.initial_energy,
                                       g=response.operator.g,omega_orb=response.operator.omega_orb,
                                       ell_f=response.operator.ell_f,m_f=c.m_f,kinematics_mode=kinematics_mode,
                                       ell_i=c.ell_i,m_i=c.m_i)
