"""Scalar radial source; coherent l sum, radial Coulomb-operator descriptor."""
from __future__ import annotations
import importlib
from dataclasses import dataclass
from enum import StrEnum
import numpy as np
angular=importlib.import_module('05_angular_coupling')

class SpectralSector(StrEnum):
    BOUND='bound'
    THRESHOLD='threshold'
    CONTINUUM='continuum'

class SourceCompletenessStatus(StrEnum):
    COMPLETE='complete'
    INCOMPLETE_MISSING_FORCING='incomplete_missing_forcing'

def driven_energy(initial_energy,g,omega_orb):return initial_energy+g*omega_orb

def classify_driven_energy(E):
    if E<0:return SpectralSector.BOUND
    if E==0:return SpectralSector.THRESHOLD
    return SpectralSector.CONTINUUM

@dataclass(frozen=True)
class RadialOperatorDescriptor:
    central_mass:float
    scalar_mass:float
    ell_f:int
    initial_energy:float
    g:int
    omega_orb:float
    @property
    def energy(self):return driven_energy(self.initial_energy,self.g,self.omega_orb)
    @property
    def spectral_sector(self):return classify_driven_energy(self.energy)

@dataclass(frozen=True)
class PartialSourceContribution:
    l:int
    m:int
    g:int
    gaunt_coefficient:float
    selection_report:object
    forcing:object
    source:object
    absolute_error_bound:object

@dataclass(frozen=True)
class ScalarSourceChannel:
    ell_i:int
    m_i:int
    ell_f:int
    m_f:int
    g:int
    required_m:int
    driven_energy:float
    spectral_sector:SpectralSector
    radial_points:object
    required_l_values:object
    l_values:object
    missing_l_values:object
    gaunt_coefficients:object
    partial_contributions:tuple
    partial_sources:object
    forcing_coefficients:object
    available_source_sum:object
    total_source:object
    inherited_forcing_metadata:tuple
    spectrum_convergence:object
    available_absolute_error_bound:object
    completeness_status:SourceCompletenessStatus
    radial_operator:RadialOperatorDescriptor


def scalar_source_channel(radial_points,*,initial_state,forcing_provider,forcing_basis,ell_f,m_f,g):
    radii=np.asarray(radial_points,dtype=float)
    if forcing_provider.output_basis!=forcing_basis:raise ValueError('Fourier basis mismatch')
    if forcing_provider.central_mass!=initial_state.parameters.M:raise ValueError('central masses disagree')
    ell_i=initial_state.quantum_numbers.ell;m_i=initial_state.quantum_numbers.m
    m=m_f-m_i
    required=angular.allowed_forcing_degrees(ell_i=ell_i,m_i=m_i,ell_f=ell_f,m_f=m_f)
    available=forcing_provider.available_degrees
    missing=tuple(l for l in required if available is not None and l not in available)
    supplied=tuple(l for l in required if l not in missing)
    mu=initial_state.parameters.mu
    common=2*mu**2*radii*initial_state.radial(radii)
    partials=[];sources=[];fsamples=[];errors=[];gaunts=[]
    for l in supplied:
        report=angular.selection_status(l,m,ell_i,m_i,ell_f,m_f)
        ff=forcing_provider.sample(l,m,g,radii)
        gaunt=report.gaunt_value
        source=common*gaunt*np.asarray(ff.coefficient)
        err=np.abs(common*gaunt)*np.asarray(ff.estimated_absolute_error)
        partials.append(PartialSourceContribution(l,m,g,gaunt,report,ff,source,err))
        sources.append(source);fsamples.append(ff);errors.append(err);gaunts.append(gaunt)
    partial_arrays=np.stack(sources) if sources else np.empty((0,len(radii)),complex)
    total=np.sum(partial_arrays,axis=0)
    error=np.sum(np.stack(errors),axis=0) if errors else np.zeros(len(radii))
    operator=RadialOperatorDescriptor(initial_state.parameters.M,mu,ell_f,initial_state.energy,g,forcing_provider.omega_orb)
    return ScalarSourceChannel(ell_i,m_i,ell_f,m_f,g,m,operator.energy,operator.spectral_sector,
        radii,np.asarray(required,int),np.asarray(supplied,int),np.asarray(missing,int),np.asarray(gaunts,float),
        tuple(partials),partial_arrays,
        np.stack([x.coefficient for x in fsamples]) if fsamples else np.empty((0,len(radii)),complex),
        total,total if not missing else None,tuple(fsamples),forcing_provider.spectrum_convergence,error,
        SourceCompletenessStatus.COMPLETE if not missing else SourceCompletenessStatus.INCOMPLETE_MISSING_FORCING,operator)
