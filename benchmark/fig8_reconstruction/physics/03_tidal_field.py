"""Direct + indirect Newtonian multipoles in the accelerated primary frame."""
from __future__ import annotations
import importlib
from dataclasses import dataclass
import numpy as np
from scipy.special import sph_harm_y
KeplerOrbit=importlib.import_module('02_kepler_orbit').KeplerOrbit

@dataclass(frozen=True)
class RadialBranchStructure:
    kind: str
    eccentric_anomaly_crossings: tuple[float,...]
    mean_anomaly_crossings: tuple[float,...]

def radial_branch_structure(orbit: KeplerOrbit,field_radius:float):
    r=float(field_radius);rp=orbit.periapsis_radius;ra=orbit.apoapsis_radius
    if orbit.eccentricity==0:return RadialBranchStructure('circular',(),())
    if r<rp:return RadialBranchStructure('always_interior',(),())
    if r>ra:return RadialBranchStructure('always_exterior',(),())
    cosE=(1-r/orbit.semi_major_axis)/orbit.eccentricity
    first=np.arccos(np.clip(cosE,-1,1))
    second=2*np.pi-first
    return RadialBranchStructure('crossing',(float(first),float(second)),
                                 (float(first-orbit.eccentricity*np.sin(first)),float(second-orbit.eccentricity*np.sin(second))))

def _source_harmonic(state,ell,m):
    return np.conj(sph_harm_y(ell,m,np.pi/2,np.asarray(state.azimuth_wrapped)))

def direct_tidal_multipole(orbit_state,field_radius,ell,m,mass_ratio):
    field=np.asarray(field_radius,dtype=float);rstar=np.asarray(orbit_state.radius,dtype=float)
    mp=mass_ratio*orbit_state.central_mass
    radial=np.where(field<=rstar, field**ell/rstar**(ell+1),rstar**ell/field**(ell+1))
    return -4*np.pi*mp/(2*ell+1)*_source_harmonic(orbit_state,ell,m)*radial

def indirect_dipole_multipole(orbit_state,field_radius,m,mass_ratio):
    rstar=np.asarray(orbit_state.radius,dtype=float)
    return 4*np.pi*mass_ratio*orbit_state.central_mass/3*_source_harmonic(orbit_state,1,m)*np.asarray(field_radius)/rstar**2

def tidal_multipole(orbit_state,field_radius,ell,m,mass_ratio):
    if ell!=1:return direct_tidal_multipole(orbit_state,field_radius,ell,m,mass_ratio)
    field=np.asarray(field_radius,dtype=float);rstar=np.asarray(orbit_state.radius,dtype=float)
    # Exact cancellation in region field < rstar, including equality.
    mp=mass_ratio*orbit_state.central_mass
    # Evaluate only the physically permitted exterior formula to avoid singular 0*inf.
    exterior=-4*np.pi*mp/3*_source_harmonic(orbit_state,1,m)*(rstar/field**2-field/rstar**2)
    out=np.where(field>rstar,exterior,0.0j)
    return complex(out) if np.ndim(out)==0 else out
