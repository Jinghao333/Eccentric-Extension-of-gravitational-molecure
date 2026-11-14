# Boson Cloud Binary Dynamics

This context names the physical and numerical concepts used to evolve a boson cloud in the time-dependent gravitational field of a binary black hole.

## Language

**Binary separation**:
The instantaneous coordinate distance between the primary and secondary black holes.
_Avoid_: Orbital radius

**COM radius**:
The distance of one black hole from the binary center of mass; the primary and secondary COM radii sum to the binary separation.
_Avoid_: Separation of a body

**Orbital angle**:
The signed azimuthal coordinate of the primary black hole in the center-of-mass frame; the secondary is displaced by pi.
_Avoid_: Orbital phase

**Accumulated orbital phase**:
A strictly increasing measure of completed orbital motion used to select averaging windows, including for counter-rotating motion. For prescribed eccentric motion this is the unwrapped mean anomaly, not the signed orbital angle.
_Avoid_: Orbital angle

**Fixed circular orbit**:
A moving circular trajectory whose separation and black-hole masses do not respond to field power, gravitational-wave losses, or accretion.
_Avoid_: Static binary

**Fixed eccentric orbit**:
A moving prescribed Kepler ellipse whose masses, eccentricity and semimajor axis remain fixed while instantaneous separation, position and velocities vary around the orbit.
_Avoid_: Static eccentric binary, eccentric inspiral

**Mean motion**:
The constant Kepler clock rate `n = 2 pi / P` of a fixed ellipse. It differs from the instantaneous angular speed on an eccentric orbit.
_Avoid_: Instantaneous orbital frequency

**Evolving circular orbit**:
A circular trajectory whose separation and masses respond to gravitational-wave losses, phase-averaged cloud power, and absorption-derived accretion.
_Avoid_: Free orbit

**Startup interval**:
The configured initial time interval over which orbital frequency, and optionally companion-potential strength, is smoothly introduced.
_Avoid_: Absorption ramp, backreaction ramp

**Horizon absorption**:
Removal of the cloud's local s-wave component by a compact sink centered on each black hole.
_Avoid_: Outer absorption, CAP

**Orbital feedback**:
The post-startup update of binary separation and masses from gravitational-wave losses, cloud power, and accretion.
_Avoid_: Field evolution

**Binary power**:
The inner-region expectation value of the orbital-motion derivative `R_dot partial_R V + theta_dot partial_theta V`; artificial companion-potential startup switching work is reported separately.
_Avoid_: Energy flux

**Complex absorbing potential**:
A negative-imaginary outer-domain potential included in the Hamiltonian to remove outgoing radiation.
_Avoid_: Damping mask, horizon absorption

**Damping mask**:
A multiplicative outer-domain attenuation applied to the field after a timestep, outside the Hamiltonian.
_Avoid_: Complex absorbing potential
