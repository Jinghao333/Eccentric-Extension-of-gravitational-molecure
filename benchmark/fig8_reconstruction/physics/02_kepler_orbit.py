"""Prescribed primary-centred Newtonian Kepler orbit, signed anomalies."""
from __future__ import annotations
from dataclasses import dataclass
import math
import numpy as np
from scipy.optimize import newton

def _out(x):
    a=np.asarray(x);return float(a) if a.ndim==0 else a

def wrap_angle(x):return _out(np.mod(np.asarray(x),2*np.pi))

@dataclass(frozen=True)
class OrbitState:
    central_mass: float
    p: float
    eccentricity: float
    mean_anomaly: object
    eccentric_anomaly: object
    true_anomaly: object
    azimuth: object
    azimuth_wrapped: object
    radius: object
    radial_velocity: object
    angular_velocity: object
    x: object
    y: object

@dataclass(frozen=True)
class KeplerOrbit:
    central_mass: float
    p: float
    eccentricity: float
    periapsis_time: float = 0.0
    periapsis_azimuth: float = 0.0
    def __post_init__(self):
        if self.central_mass<=0 or self.p<=0 or not(0<=self.eccentricity<1):
            raise ValueError('bound Kepler orbit required')
    @property
    def semi_major_axis(self):return self.p*self.central_mass/(1-self.eccentricity**2)
    @property
    def omega_orb(self):return ((1-self.eccentricity**2)/self.p)**1.5/self.central_mass
    @property
    def periapsis_radius(self):return self.p*self.central_mass/(1+self.eccentricity)
    @property
    def apoapsis_radius(self):return self.p*self.central_mass/(1-self.eccentricity)
    @property
    def specific_angular_momentum(self):return self.central_mass*math.sqrt(self.p)
    def mean_anomaly_from_time(self,t):return _out(self.omega_orb*(np.asarray(t)-self.periapsis_time))
    def state_from_mean_anomaly(self,mean_anomaly,*,tolerance=1e-13,max_iterations=100):
        mean=np.asarray(mean_anomaly,dtype=float)
        sol=newton(lambda E:E-self.eccentricity*np.sin(E)-mean,mean,
                   fprime=lambda E:1-self.eccentricity*np.cos(E),tol=tolerance,maxiter=max_iterations)
        return self.state_from_eccentric_anomaly(sol)
    def state_from_time(self,t):return self.state_from_mean_anomaly(self.mean_anomaly_from_time(t))
    def state_from_eccentric_anomaly(self,eccentric_anomaly):
        E=np.asarray(eccentric_anomaly,dtype=float);e=self.eccentricity
        Ma=E-e*np.sin(E)
        # Signed atan2 retains the second half of orbit; unwrap across 2pi for returned true anomaly.
        Ewrap=np.mod(E,2*np.pi)
        nu=np.mod(np.arctan2(np.sqrt(1-e*e)*np.sin(Ewrap),np.cos(Ewrap)-e),2*np.pi)
        phi=self.periapsis_azimuth+nu
        phiwrap=np.mod(phi,2*np.pi)
        r=self.semi_major_axis*(1-e*np.cos(E))
        vr=np.sqrt(self.central_mass/self.semi_major_axis)*e*np.sin(E)/(1-e*np.cos(E))
        vphi=self.specific_angular_momentum/r**2
        return OrbitState(self.central_mass,self.p,e,_out(Ma),_out(E),_out(nu),
                          _out(phi),_out(phiwrap),_out(r),_out(vr),_out(vphi),
                          _out(r*np.cos(phiwrap)),_out(r*np.sin(phiwrap)))
