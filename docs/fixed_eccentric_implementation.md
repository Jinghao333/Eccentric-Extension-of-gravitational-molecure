# Fixed eccentric orbit implementation

## Scope and terminology

Stage I prescribes a bound Kepler ellipse while evolving the scalar field in
its time-dependent potential. “Fixed orbit” means that the component masses,
eccentricity, semimajor axis and mean motion do not respond to the cloud or to
gravitational-wave losses. It does **not** mean fixed source positions:
`R`, `theta`, `R_dot` and `theta_dot` change around the ellipse. Existing
circular inspiral remains available only at `e = 0`; eccentric inspiral is a
future orbit producer and is not implemented here.

All quantities use `G = c = M1_initial = 1`. The ordinary-binary COM
convention is

\[
\mathbf{x}_1=+\frac{M_2}{M_1+M_2}R(\cos\theta,\sin\theta,0),\qquad
\mathbf{x}_2=-\frac{M_1}{M_1+M_2}R(\cos\theta,\sin\theta,0).
\]

## Implementation ledger

| Category | Source file/function | Previous behavior/formula | New behavior/formula | Test | Status |
| --- | --- | --- | --- | --- | --- |
| A | `include/OrbitState.hpp`; `CrankNicolsonModes` ownership | Independently mutable `R`, `R1`, `R2`, `M1`, `M2`, `q`, `theta`, `phase`, `Omega`, `P` | One authoritative `OrbitState {M1,M2,R,theta,R_dot,theta_dot,phase}`; `q`, total mass, `R1`, `R2` and `|theta_dot|` are getters | `orbit_state`, circular regression suite | Complete |
| A | `F1`, `F2`, `orbital_phase`, potential builders | Read solver orbital scalars implicitly | Consume `const OrbitState&`; source radii are instantaneous `R1()` and `R2()` and angular factor remains `exp(i m_* theta)` with the existing secondary parity | `parallel_correctness`, `physics`, `eccentric_orbit` | Complete |
| A | absorption projector/profile builders | Cached positions derived from circular scalar state | Consume `const OrbitState&` and rebuild kernels/projectors at the instantaneous source positions; normalization and positive sink are unchanged | `fixed_orbit_absorption`, `spherical_absorption`, `eccentric_orbit` | Complete |
| A | isolated black hole | Displacement shared the old `R` field | Explicit exceptional source radius `orbit.R`; ordinary COM getters are not used | `isolated_bh` | Complete |
| B | `FixedKeplerOrbit::configure/evaluate` | Circular `Omega=sqrt(M/R^3)`, constant `R` | `a=R_peri/(1-e)`, `n=sqrt(M_T/a^3)`, `E-e sin E=M`, `R=a(1-e cos E)`, true anomaly from `atan2`; exact eccentric `R_dot` and `theta_dot` | `kepler_orbit` | Complete |
| B | `build_static_perturbation` | `(dV/dt)_orbit = theta_dot partial_theta V` | `(dV/dt)_orbit = R_dot partial_R V + theta_dot partial_theta V`; radial term enabled for prescribed eccentric motion while the legacy circular-inspiral diagnostic path is retained | `eccentric_orbit` finite difference and circular limit | Complete |
| B | companion switching derivative | Switching was implicit in `V=V1+lambda V2` and absent from power output | `dVdt_switch=lambda_dot V2`; `p_ramp_inst=<psi|dVdt_switch|psi>` remains separate from mechanical orbital power | `eccentric_orbit` switching finite difference | Complete |
| B | period users | `P=2 pi/Omega` | Prescribed eccentric period is `P0=2 pi/n`; instantaneous `|theta_dot|` is diagnostic only | `kepler_orbit`, `command_line`, dashboard Python test | Complete |
| C | `advance_orbit`, `update_orbit_and_ramp` | Circular angle/phase advanced by the previous ramped `Omega`, then endpoint `Omega` evaluated | Positive mean anomaly advances by `M_{k+1}=M_k+n g(t_k) dt`; endpoint geometry and velocities use `g(t_{k+1})`. The endpoint Hamiltonian convention is unchanged | `orbital_ramp`, `eccentric_orbit` | Complete |
| C | `step` | Hamiltonian assembly also controlled circular scheduling state | Orbit producer advances first; potential and absorption consume the endpoint `OrbitState`; CN propagation, mask, diagnostics, and post-field circular feedback retain their ordering | circular suite, `eccentric_orbit` | Complete |
| C | circular feedback | Circular GW/cloud/accretion update after propagation | Unchanged for `e=0`; `e>0 && evolve_orbit` is rejected | `physics`, `orbital_ramp`, configuration tests | Complete |
| D | `Config`, `main` | No eccentricity key; run length used circular period | `[binary].e` defaults to zero and requires `0<=e<1`; `[binary].R` is pericenter for `e>0`; run length uses `2 pi/n` | `config_and_io`, `command_line` | Complete |
| D | diagnostics | No explicit phase, radial velocity, or switching-work columns | Adds `phase_orb.txt`, `R_dot_orb.txt`, `p_ramp_inst.txt`; `Omega_orb.txt=|theta_dot|` remains compatible | `eccentric_orbit`, `command_line`, `restart` | Complete |
| D | checkpoint/restart | `GMCKPT01` lacked radial velocity and signed rate | `GMCKPT02` explicitly stores complete orbit fields; loader accepts v1; old circular folders backfill only missing new columns | `restart`, `eccentric_orbit` | Complete |
| D | initial-state scripts | Molecular generator implied a stationary rotating-frame state for any config | Atomic translation remains valid at the initial pericenter; molecular generator rejects `e!=0` with its circular limitation stated | Python initial-state test | Complete |
| D | dashboard | Used `2 pi/|Omega|` and circular theory for every run | Uses `a=R_peri/(1-e)` and `P0=2 pi sqrt(a^3/M_T)` for eccentric runs; hides circular-only theory/reference lines | Python dashboard test | Complete |
| D | documentation/example | Circular configuration only | Documents eccentric conventions, startup, power split, restrictions, and `examples/fixed_eccentric/config.toml` | documentation review | Complete |

## Kepler producer

`FixedKeplerOrbit` stores only fixed trajectory parameters: `e`, `a`, `n`, and
the direction sign `s`. The instantaneous state is written only to the solver's
authoritative `orbit`. For initial pericenter `R0`,

\[
a=\frac{R_0}{1-e},\qquad n=\sqrt{\frac{M_1+M_2}{a^3}},\qquad P_0=\frac{2\pi}{n}.
\]

The positive unwrapped `orbit.phase` is mean anomaly `M`. The bracketed Newton
solver preserves full revolutions while solving

\[
E-e\sin E=\mathcal M.
\]

At each endpoint it evaluates

\[
R=a(1-e\cos E),\quad
f=\operatorname{atan2}(\sqrt{1-e^2}\sin E,\cos E-e),\quad
\theta=s f,
\]

with `f` unwrapped consistently, and

\[
\dot R=g\frac{a n e\sin E}{1-e\cos E},\qquad
\dot\theta=g\,s\frac{n\sqrt{1-e^2}}{(1-e\cos E)^2}.
\]

Thus counter-rotation changes the signs of `theta` and `theta_dot`, while mean
anomaly remains strictly increasing for rolling averages.

## Startup and timestep ordering

For ramp duration `T=ramp_periods*P0`,

\[
g(t)=3x^2-2x^3,\quad x=t/T,
\]

saturated to zero before and one after the interval. If `T=0`, `g=1` from the
start. On timestep `k`, the code deliberately performs

\[
\mathcal M_{k+1}=\mathcal M_k+n g(t_k)\Delta t,
\]

then evaluates `R`, `theta`, `R_dot`, `theta_dot`, the optional potential
weight `lambda`, absorption state, and Hamiltonian at `t_{k+1}`. This is the
legacy discrete phase/endpoint-Hamiltonian convention, not an analytic ramp
integral and not a midpoint Hamiltonian.

If companion ramping is enabled, `lambda=g` and

\[
\dot\lambda=\frac{6x-6x^2}{T}
\]

inside the startup interval. Otherwise `lambda=1` and `dot(lambda)=0`.
Horizon absorption remains inactive until the end of startup, and circular
backreaction begins afterward with its existing activation/order. Neither
activation is folded into the Kepler trajectory.

## Potential derivative and power meanings

For multipole radial profile

\[
F_\ell(r,d)=\frac{\min(r,d)^\ell}{\max(r,d)^{\ell+1}},\qquad d_h=c_hR,
\]

the implementation uses

\[
\partial_R F_\ell=c_h\begin{cases}
-(\ell+1)r^\ell/d_h^{\ell+2},&r<d_h,\\
\ell d_h^{\ell-1}/r^{\ell+1},&r>d_h.
\end{cases}
\]

Mass weights, Gaunt/Wigner coefficients, mode ordering, secondary parity, and
sparse patterns are unchanged. For `m_*=m_col-m_row`,

\[
(\partial_tV)_{AB}^{orbit}=\dot R[\partial_RV_1+\lambda\partial_RV_2]_{AB}
+i m_*\dot\theta[V_1+\lambda V_2]_{AB}.
\]

Diagnostics mean:

- `p_bin_inst`: inner-region expectation of the orbital expression above.
- `p_bin`: unchanged rolling phase average of `p_bin_inst`.
- `p_ramp_inst`: inner-region expectation of `dot(lambda) V2`, kept separate
  because it is artificial startup work.
- `dVdt`: the in-memory orbital derivative matrix only.
- `dV_bin_dt`: legacy `Re <dpsi/dt|V|psi>` diagnostic; it is not the complete
  explicit potential derivative.

## Ownership and dependency boundary

`CrankNicolsonModes::orbit` is the single instantaneous authority. Only the
trajectory/startup update, legacy circular feedback, isolated-geometry setup,
and restart reconstruction write it. Hamiltonian, absorption, diagnostics and
checkpoint output read it. `FixedKeplerOrbit` owns fixed elements but never a
shadow instantaneous state. The scalar solver continues to own time, `dt`,
ramp configuration, field vectors, sparse operators and rolling averages.

This boundary deliberately permits a future eccentric-inspiral producer to
update the same `OrbitState` and any required evolving elements without changing
the gravitational multipole or scalar solver interfaces. No eccentric
backreaction equations are included in Stage I.

## Restart and compatibility

`GMCKPT02` writes named fields in a fixed order rather than raw-serializing the
C++ struct: time, angle, unwrapped phase, separation, radial velocity, signed
angular velocity, both masses, timestep and rolling binary power, followed by
the field and rolling-average state. `GMCKPT01` remains loadable; its unsigned
frequency is assigned the configured rotation sign and its missing radial
velocity is zero. Diagnostic-only reconstruction reads the new phase directly,
never integrates eccentric `Omega_orb`. Missing new files are synthesized only
for old circular folders.

For `e=0`, the Kepler producer reduces to constant `R`, `R_dot=0`,
`theta=s*phase`, and `|theta_dot|=n g`. Existing circular inspiral, isolated-BH
geometry, Wigner/Gaunt construction, radial discretization, CN solve,
preconditioners, CAP, damping mask and `RollingPhaseAverage` are unchanged.

## Known limitations and future work

- Eccentric GW reaction, cloud backreaction, accretion-driven mass evolution,
  and evolving orbital elements are not implemented; `e>0` requires
  `evolve_orbit=false`.
- The molecular rotating-frame generator is circular-only. Atomic initial
  states translated to the initial pericenter are supported.
- The existing CN method evaluates the time-dependent Hamiltonian at the step
  endpoint. A midpoint method would be a separate numerical-method change.
- Absorption profiles are rebuilt at the instantaneous eccentric positions;
  no interpolation/cache approximation is introduced.

## Verification record

The final Release build completed with
`cmake --build build-codex-verify --config Release --parallel 4`.
`ctest --test-dir build-codex-verify --output-on-failure -C Release` reported
16 passed tests and one dependency skip: CMake's selected Python 3.12 lacks
SymPy. The skipped initial-state suite was run directly with the available
Python environment and all six checks passed, including a generated-state
solver run. The dashboard's six Python checks also passed against real solver
output. No long scientific baseline was rerun.

Measured focused errors were:

- maximum Kepler residual through `e=0.999999`: `1.77636e-15`;
- maximum specific-energy and angular-momentum errors over an `e=0.7` orbit:
  `1.04083e-16` and `1.77636e-15`;
- centered finite-difference relative error for orbital `dVdt`: `1.71708e-12`;
- centered finite-difference relative error for switching `dVdt`:
  `1.75603e-14`;
- moving eccentric sink/reference relative errors at two endpoint positions:
  `1.43803e-18` and `2.84933e-18`;
- eccentric checkpoint/restart diagnostics and checkpoint: byte-for-byte equal
  to uninterrupted evolution at the tested fixed thread count.
