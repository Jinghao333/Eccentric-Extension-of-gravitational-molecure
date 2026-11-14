# Circular-orbit implementation audit

## 1. Purpose, scope, and evidence labels

This document is a read-only audit of the circular binary implementation as it exists before an `OrbitState` migration. It traces the executable beginning at `src/main.cpp`, follows helper calls, and records numerical choices that a behavior-preserving refactor must not silently change.

Evidence labels used below:

- **Verified**: stated directly by the current source or an existing test.
- **Physical assumption**: the mathematical model inferred from the implemented equations.
- **Suspected inconsistency**: two implementation choices appear misaligned, or a quantity's name/documentation does not match its exact computation. These are not corrected here.
- **Later proposal**: an improvement deliberately deferred beyond the circular behavior-preservation stage.

Source references use the repository state audited on 2026-11-14. Line numbers refer to the production and test files before any future refactor.

## 2. Units, field normalization, and coordinates

### 2.1 Units

**Verified.** The solver uses geometrized units `G = c = 1`, with the initial primary mass set to one (`include/CrankNicolsonModes.hpp:42-44`, `src/CrankNicolsonModes.cpp:30-31`):

\[
M_1(0)=1,\qquad M_2(0)=q.
\]

Lengths and times are therefore expressed in units of the initial primary mass. The boson mass `mu` has inverse-length/inverse-time units. The Schrödinger Hamiltonian has inverse-time units.

### 2.2 Field decomposition and normalization

The field is represented as (`include/CrankNicolsonModes.hpp:6-12`)

\[
\psi(t,r,\vartheta,\varphi)
=\sum_{LM}\frac{U_{LM}(t,r)}{r}Y_{LM}(\vartheta,\varphi).
\]

The state vector is mode-major: `psi[j*Nr+i]` is \(U_{L_jM_j}(r_i)\). The continuum norm is

\[
N=\sum_{LM}\int |U_{LM}|^2\,dr,
\]

and cloud mass is \(M_c=\mu N\). The initial-state scripts normalize with a trapezoidal radial integral (`initial_states/get_initial_state_COM_fast.py:162-170`; `initial_states/get_initial_state_COM.py:495-518`). `set_initial_state_full()` then multiplies the loaded state by `sqrt(M_c/mu)` (`src/CrankNicolsonModes.cpp:130-142`). `tests/config_and_io.cpp:216-225` verifies

\[
\mu\,dr\sum_k |\psi_k|^2=M_c
\]

for a unit-discrete-norm input.

### 2.3 Center-of-mass frame and body positions

The spherical-harmonic origin is the binary center of mass. The orbit is equatorial. The primary is at

\[
\mathbf x_1=R_1(\cos\theta,\sin\theta,0),
\]

and the secondary is opposite it,

\[
\mathbf x_2=-R_2(\cos\theta,\sin\theta,0)
=R_2(\cos(\theta+\pi),\sin(\theta+\pi),0).
\]

This convention is verified independently by the absorption reference (`tests/SphericalAbsorptionReference.hpp:60-74`) and direct three-dimensional quadrature (`tests/spherical_absorption.cpp:44-65`). The centroid test at `tests/spherical_absorption.cpp:93-111` checks the primary center against \(R_1e^{i\theta}\).

For a circular COM configuration,

\[
R=R_1+R_2,\qquad
R_1=\frac{M_2}{M_1+M_2}R,\qquad
R_2=\frac{M_1}{M_1+M_2}R.
\]

Initialization uses the equivalent `q` formulas (`src/CrankNicolsonModes.cpp:24-25`).

## 3. Function-level execution graph

### 3.1 New run

```text
main                                                src/main.cpp:20-109
├─ load_config                                      include/Config.hpp:94-155
├─ CrankNicolsonModes::CrankNicolsonModes           src/CrankNicolsonModes.cpp:11-50
│  ├─ build_grid                                  src/CrankNicolsonModes.cpp:100-105
│  ├─ build_modes                                 src/CrankNicolsonModes.cpp:107-128
│  └─ build_stationary_hamiltonian                src/CrankNicolsonModes_hamiltonian.cpp:44-66
├─ load_initial_modes_binary                        src/InitialData.cpp:7-79
├─ set_initial_state_full                           src/CrankNicolsonModes.cpp:130-143
├─ set_power_average_orbits                         src/CrankNicolsonModes.cpp:52-56
├─ initialize                                       src/CrankNicolsonModes.cpp:76-98
│  ├─ build_coupling_cache                       src/CrankNicolsonModes_hamiltonian.cpp:68-104
│  ├─ build_sparse_element_caches                src/CrankNicolsonModes_hamiltonian.cpp:111-174
│  ├─ build_time_dependent_perturbation(0)
│  │  ├─ build_static_perturbation             src/CrankNicolsonModes_hamiltonian.cpp:176-222
│  │  │  ├─ F1, F2
│  │  │  └─ orbital_phase
│  │  ├─ build_spherical_absorption (optional) src/CrankNicolsonModes_absorption.cpp:153-182
│  │  │  ├─ update_absorption_profile
│  │  │  └─ build_swave_projector
│  │  └─ assemble_hamiltonian                  src/CrankNicolsonModes_hamiltonian.cpp:224-243
│  └─ build_CN_matrices                          src/CrankNicolsonModes_matrices.cpp:12-25
│     └─ update_CN_values
├─ build_damping_mask (optional)                    src/CrankNicolsonModes_outgoing.cpp:10-35
├─ build_cap (optional)                             src/CrankNicolsonModes_outgoing.cpp:37-62
└─ run                                              src/CrankNicolsonModes.cpp:236-259
   ├─ set_extraction_radius
   ├─ create_output_folder                         src/Utils.cpp:11-27
   ├─ open_diagnostics
   ├─ save_run_info
   ├─ build_time_dependent_perturbation(0)          refreshes a post-initialize CAP
   ├─ save_snapshot(psi)
   ├─ save_snapshot(schrodinger_rhs(psi))
   ├─ save_checkpoint
   └─ evolve(1, n_steps)
      └─ step(n)                                   section 9
```

### 3.2 Resume

```text
main
├─ load config from run folder
├─ construct and initialize solver from config
├─ rebuild optional mask/CAP
└─ resume                                           src/CrankNicolsonModes_restart.cpp:288-399
   ├─ set_extraction_radius
   ├─ verify modes.txt ordering
   ├─ find common complete diagnostic row count
   ├─ find latest complete snapshot with enough rows
   ├─ load_checkpoint
   │  └─ RollingPhaseAverage::load
   ├─ or reconstruct_state_from_diagnostics
   │  └─ RollingPhaseAverage::reset/update
   ├─ truncate every diagnostic to snapshot boundary
   ├─ rebuild H(t_n) from restored orbit
   └─ evolve(start + 1, n_steps)
```

## 4. Configuration, construction, grid, and modes

### 4.1 Configuration

`load_config()` is header-defined at `include/Config.hpp:94-155`. Every key is optional. Orbit-relevant inputs are:

- `[binary].q`, `R`, `counter` (`include/Config.hpp:128-130`).
- `[simulation].dt`, `ramp_periods`, `ramp_secondary_potential`, `steps_per_period`, `power_average_orbits`, `evolve_orbit`, `absorption`, and `epsilon_smear` (`include/Config.hpp:109-126`).
- Presence of `[mask]` or `[cap]` enables the respective outer treatment (`include/Config.hpp:139-153`).

`tests/config_and_io.cpp:40-127` protects defaults, explicit values, and the limited validation currently performed. There is no production validation of positive `R`, nonnegative `q`, nonnegative `ramp_periods`, mask/CAP shape parameters, or extraction radius relative to `epsilon` until later use.

### 4.2 Constructor order

The constructor initializer list computes

\[
dr=\frac{L}{N_r-1},\qquad
P_0=2\pi\sqrt{\frac{R^3}{1+q}},
\]

before entering the body (`src/CrankNicolsonModes.cpp:11-15`). The body sets `Omega=0`, `Omega_0=sqrt((1+q)/R^3)`, `P=initial_period`, COM radii, zero angles, and masses (`src/CrankNicolsonModes.cpp:20-31`). It then constructs the grid, mode list, and `H0`.

**Suspected inconsistency.** `evolve_orbit` has a public default of `true` (`include/CrankNicolsonModes.hpp:135-137`). The constructor prints the vacuum coalescence time when that default is true (`src/CrankNicolsonModes.cpp:36-40`), but `main()` applies the configured value only after construction (`src/main.cpp:76`). A fixed-orbit run may therefore print a coalescence time even though feedback is subsequently disabled.

### 4.3 Radial grid

`build_grid()` uses (`src/CrankNicolsonModes.cpp:100-105`)

\[
r_i=\epsilon+i\,dr,\qquad dr=\frac{L}{N_r-1}.
\]

`tests/config_and_io.cpp:167-185` explicitly protects this definition.

**Suspected inconsistency.** For nonzero `epsilon`, the last solver node is `L + epsilon`, although `L` is documented as the outer radius. The COM initial-state scripts instead create `linspace(0,L,N)` and replace only the first point by `epsilon` (`initial_states/get_initial_state_COM_fast.py:148-149`; `initial_states/get_initial_state_COM.py:139-140`). The grids agree when `epsilon=0`, as in the main examples, but not in general.

### 4.4 Mode list

`build_modes()` loops in increasing `l`, then increasing `m`, and retains only even `abs(l-m)` (`src/CrankNicolsonModes.cpp:107-128`). For every retained mode, \(l-m\) is even. This is the equatorial-reflection sector excited by equatorial sources. The count is

\[
n_{\rm modes}=\sum_{l=0}^{l_{\max}}(l+1)
=\frac{(l_{\max}+1)(l_{\max}+2)}2.
\]

The ordering and parity are verified by `tests/config_and_io.cpp:167-185` and later enforced on resume by `src/CrankNicolsonModes_restart.cpp:297-316`.

## 5. Initial cloud loading and placement

### 5.1 MODESBIN loading

`load_initial_modes_binary()` reads magic `MODESBIN`, version 1, radial size, and stored mode profiles (`src/InitialData.cpp:7-56`). It then maps by `(l,m)`, ignoring stored modes absent from the solver and zero-filling solver modes absent from the file (`src/InitialData.cpp:58-78`). `tests/config_and_io.cpp:130-156` protects this behavior.

The loader does not verify that the stored radial count equals the solver's `Nr`; the size mismatch is caught later by `set_initial_state_full()` because the flattened vector size differs.

### 5.2 Initial-state coordinate convention

The recommended fast script reads `q` and `R`, sets `M1` from optional Python-only `[binary].M1` with default 1, derives `M2=q*M1`, and computes `R1`, `R2`, and circular `Omega` (`initial_states/get_initial_state_COM_fast.py:54-71`). It places the primary on `+x` and the secondary on `-x` (`initial_states/get_initial_state_COM_fast.py:73-80`, `573-586`). The slower COM script documents the same convention at `initial_states/get_initial_state_COM.py:5-16` and implements it at lines 55-80.

The COM scripts generate all `(l,m)` profiles through `l_max` (`initial_states/get_initial_state_COM_fast.py:148-151`), whereas the C++ solver retains only even `l-m`. `load_initial_modes_binary()` maps only requested solver modes and ignores the rest. For the intended equatorially symmetric translated states, the discarded sector should be zero up to numerical error; `tests/python/test_initial_states.py:145-173` checks this for the fast primary initializer.

The C++ configuration loader does not read `[binary].M1`; production evolution always begins with `M1=1`. Supplying a non-unit `M1` to an initial-state script can therefore create an initial cloud inconsistent with the C++ binary.

The molecular initializer is circular by construction. It reads `M1`, `q`, `R`, and `counter` (`initial_states/get_initial_state_molecular.py:170-204`), places the bodies on the x axis at the COM positions (`lines 606-617`), and solves the stationary co-rotating problem

\[
H_{\rm rot}(R)=H_{\rm bin}(R)-\Omega(R)L_z,
\qquad
\Omega=s\sqrt{(M_1+M_2)/R^3},
\]

as documented at lines 19-29 and implemented at lines 583-603 and 2169-2267. It is therefore sensitive to any future redefinition of circular geometry even though the C++ MODESBIN loader is not.

The Gaussian initializer is orbit-independent and populates only `(0,0)` (`initial_states/get_initial_state_gaussian.py:33-59`). It uses `linspace(0,L,N)` and does not read `simulation.epsilon`, so it shares the nonzero-`epsilon` grid mismatch noted in section 4.3.

## 6. Stationary radial Hamiltonian

For each mode, `build_stationary_hamiltonian()` constructs the centered second-order finite-difference operator on interior nodes (`src/CrankNicolsonModes_hamiltonian.cpp:44-66`):

\[
(H_0U)_{LM,i}=
-\frac{U_{i-1}-2U_i+U_{i+1}}{2\mu\,dr^2}
+\frac{L(L+1)}{2\mu r_i^2}U_i.
\]

The boundary rows contain no `H0` entries. Because the union matrix includes identity diagonals, the CN boundary equations preserve the supplied boundary values. Initial-state scripts set the radial boundaries to zero, and tests use zero boundary values.

`H0` is built once. The independent dense reference at `tests/parallel_correctness.cpp:49-106` verifies its coefficients.

## 7. Multipole expansion and angular coupling

### 7.1 Implemented Newtonian potential

**Physical assumption.** The scalar experiences the Newtonian potential energy

\[
V(\mathbf x,t)=-\mu\left(
\frac{M_1}{|\mathbf x-\mathbf x_1(t)|}
+\frac{M_2}{|\mathbf x-\mathbf x_2(t)|}
\right).
\]

The radial equations are

\[
i\partial_tU_A=
\left[-\frac{\partial_r^2}{2\mu}
+\frac{L_A(L_A+1)}{2\mu r^2}\right]U_A
+\sum_B\chi_{AB}(r,t)U_B.
\]

### 7.2 Radial source factors `F1()` and `F2()`

The addition theorem uses

\[
g_l(r;d)=\frac{r_<^l}{r_>^{l+1}},\qquad
r_<=\min(r,d),\quad r_>=\max(r,d).
\]

`F1()` implements (`src/CrankNicolsonModes_hamiltonian.cpp:12-26`)

\[
F_1(r,l_*)=\frac{1}{1+q}g_{l_*}(r;R_1).
\]

`F2()` implements (`src/CrankNicolsonModes_hamiltonian.cpp:28-42`)

\[
F_2(r,l_*,m_*)=(-1)^{m_*}\frac{q}{1+q}g_{l_*}(r;R_2).
\]

In production, `radial_f2` is filled with `F2(..., 0)` (`src/CrankNicolsonModes_hamiltonian.cpp:183-185`) and the secondary parity is applied later (`lines 198, 207`). Thus it occurs once, not twice. Multiplication by the common `M1+M2` factor yields the physical masses because `q=M2/M1`:

\[
(M_1+M_2)\frac1{1+q}=M_1,
\qquad
(M_1+M_2)\frac q{1+q}=M_2.
\]

`F1()` and `F2()` are evaluated at every perturbation rebuild, even on a fixed orbit. Their value buffers are caches of repeated mode-pair-independent work, not time-independent caches.

### 7.3 Wigner/Gaunt cache and selection rules

For row mode \(A=(L_A,M_A)\) and column mode \(B=(L_B,M_B)\), production defines (`src/CrankNicolsonModes_hamiltonian.cpp:77-80`)

\[
m_*=M_B-M_A.
\]

It considers

\[
|L_B-L_A|\le l_*\le L_A+L_B,
\]

and rejects terms when \(|m_*|>l_*\) or \(L_A+L_B+l_*\) is odd (`lines 81-84`). The cached angular factor is (`lines 85-94`)

\[
G^{l_*}_{AB}=
(-1)^{M_B}
\sqrt{\frac{(2L_A+1)(2L_B+1)(2l_*+1)}{4\pi}}
\begin{pmatrix}l_*&L_B&L_A\\0&0&0\end{pmatrix}
\begin{pmatrix}l_*&L_B&L_A\\-m_*&M_B&-M_A\end{pmatrix}
Y_{l_*m_*}\!\left(\frac\pi2,0\right),
\]

with the negative-`m` convention reconstructed from GSL's nonnegative-`m` spherical Legendre value. The cache depends only on the retained modes and is built once.

### 7.4 Orbital phase and body parity

`orbital_phase()` returns (`src/CrankNicolsonModes_hamiltonian.cpp:106-109`)

\[
e^{im_*\theta}.
\]

The secondary lies at `theta + pi`, so its harmonic contributes

\[
e^{im_*(\theta+\pi)}=(-1)^{m_*}e^{im_*\theta}.
\]

The implemented matrix elements are therefore

\[
(V_1)_{AB}(r,t)=
e^{i(M_B-M_A)\theta}
\sum_{l_*}
-\mu M_1\frac{4\pi}{2l_*+1}G^{l_*}_{AB}g_{l_*}(r;R_1),
\]

\[
(V_2)_{AB}(r,t)=
e^{i(M_B-M_A)\theta}
\sum_{l_*}
-\mu M_2\frac{4\pi}{2l_*+1}G^{l_*}_{AB}(-1)^{M_B-M_A}g_{l_*}(r;R_2).
\]

`tests/parallel_correctness.cpp:49-106` independently reconstructs these formulas without the production caches and compares `H0`, `V1`, `V2`, and `dVdt`.

### 7.5 Sparse storage

`build_sparse_element_caches()` creates one radial-diagonal entry for each coupled mode pair and interior radius (`src/CrankNicolsonModes_hamiltonian.cpp:111-149`). `V1`, `V2`, `V`, and `dVdt` share the pattern. The union pattern of `H`, `A`, and `B` contains `H0`, every potential coupling, and every identity diagonal, including locations needed by a CAP (`lines 150-174`). Compressed offsets are cached, so later steps update values without structural insertions.

## 8. The four startup/activation schedules

The code contains four related but distinct policies. All are anchored to

\[
t_{\rm ramp}=N_{\rm ramp}P_0,
\qquad
P_0=\texttt{initial_period},
\]

where `initial_period` is constant (`include/CrankNicolsonModes.hpp:242-246`). Changes in the evolving period do not move these gates. `tests/orbital_ramp.cpp:55-99` explicitly protects this.

Define

\[
S(x)=3x^2-2x^3.
\]

### 8.1 Frequency ramp

At each Hamiltonian rebuild (`src/CrankNicolsonModes_hamiltonian.cpp:245-261`):

\[
\Omega(t)=
\begin{cases}
\Omega_0 S(t/t_{\rm ramp}),&t<t_{\rm ramp},\\
\sqrt{(M_1+M_2)/R^3},&t\ge t_{\rm ramp}.
\end{cases}
\]

This ramp applies even when `evolve_orbit=false`. It controls rotation, not whether separation and masses evolve.

### 8.2 Optional companion-potential ramp

`potential_ramp_factor` is (`src/CrankNicolsonModes_hamiltonian.cpp:262`)

\[
\lambda(t)=
\begin{cases}
S(t/t_{\rm ramp}),&\texttt{ramp_secondary_potential=true},\\
1,&\texttt{false}.
\end{cases}
\]

and

\[
V=V_1+\lambda V_2.
\]

The primary potential is never ramped. `tests/orbital_ramp.cpp:119-142` protects both potential policies while requiring the same frequency ramp in both cases.

### 8.3 Horizon-absorption activation

`build_spherical_absorption()` uses (`src/CrankNicolsonModes_absorption.cpp:158-171`)

\[
A(t)=
\begin{cases}
0,&t<t_{\rm ramp},\\
S((t-t_{\rm ramp})/P_0),&t_{\rm ramp}\le t<t_{\rm ramp}+P_0,\\
1,&t\ge t_{\rm ramp}+P_0.
\end{cases}
\]

The horizon-area weights are

\[
A_{h,a}=16\pi M_a^2 A(t).
\]

`assemble_hamiltonian()` marks sinks active only when `time >= ramp_end_time()` (`src/CrankNicolsonModes_hamiltonian.cpp:264-267`). At the exact gate, the sink is active but its strength is zero. `tests/fixed_orbit_absorption.cpp:14-80` and `tests/orbital_ramp.cpp:55-99` protect the schedule.

### 8.4 Orbital-backreaction activation

`update_time_dependent_increment()` applies the same one-`P0` post-startup smoothstep to the separation drift (`src/CrankNicolsonModes_evolution.cpp:32-50`). The function is not called until `t_i > t_ramp` (`src/CrankNicolsonModes.cpp:189-195`). If a step crosses the gate, only

\[
\Delta t_{\rm orbit}=t_{n+1}-\max(t_n,t_{\rm ramp})
\]

is used, but the smoothstep weight is its endpoint value. This is forward-Euler-like, not an integral of the ramp across the step. `tests/orbital_ramp.cpp:101-117` protects that exact crossing behavior.

### 8.5 `ramp_periods=0`

At initialization, `time < ramp_end_time()` is false because both are zero. Therefore:

- `Omega=Omega_0` immediately.
- `ramp_factor=1` and the full companion potential is present.
- `assemble_hamiltonian(absorption && time >= 0)` activates absorption immediately with full strength because the absorption ramp's second branch also has zero duration before `P0`: at `t=0`, its `x` is zero, so the strength is zero; it then ramps over one `P0`.
- Orbital feedback is eligible after the first positive-time endpoint and is smoothed over the first `P0`.

`tests/orbital_ramp.cpp:144-155` protects zero-ramp fixed motion. `tests/fixed_orbit_absorption.cpp:91-100` protects the remaining one-period absorption activation.

## 9. Exact timestep ordering

The ordering in `CrankNicolsonModes::step()` (`src/CrankNicolsonModes.cpp:145-234`) is a compatibility contract for the refactor.

Let the solver enter a step with state \(U^n\), time \(t_n\), timestep \(\Delta t_n\), and orbit variables left by the previous step.

1. `psi_last = psi` saves \(U^n\).
2. `step_start = t_i`; then `t_i += dt`, so `t_i=t_{n+1}`.
3. The coordinate angle advances with the **previous** frequency and completed-step timestep:

   \[
   \theta_{n+1}=\theta_n+s\Omega_n\Delta t_n,
   \qquad s=\begin{cases}+1&\text{co-rotating},\\-1&\text{counter-rotating}.
   \end{cases}
   \]

4. The averaging phase advances without the sign:

   \[
   \Phi_{n+1}=\Phi_n+\Omega_n\Delta t_n.
   \]

5. `build_time_dependent_perturbation(t_i)` evaluates the ramp at the endpoint, sets a new `Omega`, rebuilds potentials/absorbers, and assembles `H`.
6. `update_CN_values()` forms `A` and `B` from that endpoint Hamiltonian.
7. The CN system advances the field.
8. An optional damping mask multiplies the solved endpoint state.
9. `dpsi_dt=schrodinger_rhs(psi)` uses the post-mask endpoint state and the just-built Hamiltonian/sink. The mask is not part of this derivative.
10. If enabled and past the gate, orbital feedback updates separation and masses **after the field solve**.
11. Diagnostics are calculated. `p_bin` used in step 10 is the value carried from the previous diagnostic update; only now is new instantaneous power calculated and inserted into the rolling average.
12. Orbital scalars, field diagnostics, snapshots, and checkpoints are written.

Pseudocode matching the source exactly:

```text
psi_last = psi
step_start = t_i
t_i += dt
theta += sign(counter) * Omega * dt
phase += Omega * dt

build endpoint ramp, Omega, V1, V2, V, dVdt, sinks, H at t_i
build A and B from endpoint H and current dt
rhs = B * psi_old - optional sink RHS term
solve A_full * psi_new = rhs
apply optional post-solve damping mask
dpsi_dt = full endpoint generator applied to post-mask psi_new

if evolve_orbit and t_i > ramp_end:
    update R and masses using previous p_bin and midpoint field state
    derive q, R1, R2, Omega, P, and next-step dt

compute flux, energy, new instantaneous binary power
update rolling p_bin using t_i and phase
write orbit and field diagnostics
save snapshots/checkpoint on save interval
```

### 9.1 `theta` versus `phase`

`theta` is the signed source coordinate appearing in `exp(i*m_star*theta)` and in absorber positions. `phase` is strictly increasing bookkeeping for `RollingPhaseAverage`. For `counter=true`, `theta` decreases while `phase` increases. With a fixed zero-ramp circular orbit starting from zero, `theta=-phase`. This is protected by `tests/orbital_ramp.cpp:144-155` and `tests/fixed_orbit_absorption.cpp:60-64`.

The phase advance is a discrete left-rectangle integration of the ramped `Omega`: the old `Omega` advances the angle, and the endpoint ramp then computes the `Omega` used by the next step. It is not the analytic integral of the smoothstep. This must remain unchanged during the behavior-preserving migration.

### 9.2 Endpoint Hamiltonian

Both CN matrices use the Hamiltonian rebuilt at `t_{n+1}`, not `t_n` or `t_{n+1/2}`. `tests/parallel_correctness.cpp:139-177` independently reproduces that endpoint construction and dense CN solve. `tests/flux_time_alignment.cpp:30-118` verifies that diagnostics pair the endpoint state with the endpoint generator, including after a damping mask.

**Physical/numerical consequence.** CN is second-order for a fixed Hamiltonian, but selecting `H(t_{n+1})` for a varying Hamiltonian is not the symmetric midpoint treatment. This audit does not recommend changing it during the circular migration.

### 9.3 Fixed orbit

When `evolve_orbit=false`, `update_time_dependent_increment()` is never called. `R`, masses, `q`, COM radii, `P`, and `dt` remain fixed. The startup frequency and optional potential ramps still run, `theta` still advances, absorption can still remove field norm, and black-hole masses do not grow. `tests/fixed_orbit_absorption.cpp:14-80` protects this separation of field absorption from orbital/mass feedback.

## 10. Horizon absorption

### 10.1 Kernel geometry

Each body has a compact spherical kernel centered on the body,

\[
K_a(s)\propto
\begin{cases}
\exp[-1/(1-s^2/h_a^2)],&s<h_a,\\
0,&s\ge h_a,
\end{cases}
\qquad h_a=\epsilon_{\rm smear}M_a.
\]

`update_absorption_profile()` computes Legendre moments by angular-cap Gauss-Legendre integration (`src/CrankNicolsonModes_absorption.cpp:12-80`) and normalizes them with the discrete three-dimensional volume

\[
dr\sum_i r_i^2\int K_a\,d\Omega=1.
\]

Profiles are recomputed only when radius, width, or size changes (`lines 17-21`). A width smaller than `dr`, nonpositive width, or negative radius produces a zero profile (`line 22`).

### 10.2 S-wave projector and sink

`build_swave_projector()` reuses the potential's Gaunt cache and orbital phases (`src/CrankNicolsonModes_absorption.cpp:127-150`). It constructs a vector \(v_a\) such that

\[
dr\,v_a^\dagger U
=\int K_a(|\mathbf x-\mathbf x_a|)\psi(\mathbf x)\,d^3x.
\]

The rank-one operator is

\[
\Gamma_a=g_av_av_a^\dagger,
\qquad
g_a=16\pi M_a^2 A(t)\,dr.
\]

The full generator is

\[
H_{\rm full}=H-\frac{i}{2}\sum_a\Gamma_a.
\]

`H` itself remains sparse and excludes the dense rank-one sinks (`src/CrankNicolsonModes_hamiltonian.cpp:224-243`). `apply_sink()` and `schrodinger_rhs()` apply them separately (`src/CrankNicolsonModes_absorption.cpp:100-117`).

The independent addition-theorem implementation in `tests/SphericalAbsorptionReference.hpp:17-89`, the direct geometry tests in `tests/spherical_absorption.cpp`, and the s-wave selectivity test in `tests/swave_projector.cpp` protect this construction.

## 11. Hamiltonian and Crank–Nicolson solve

### 11.1 Hamiltonian composition

At each rebuild,

\[
V=V_1+\lambda V_2,
\qquad
H=H_0+V+\mathrm{CAP}.
\]

The CAP, if present, is a diagonal negative-imaginary potential constructed by `src/CrankNicolsonModes_outgoing.cpp:37-62` and inserted by `assemble_hamiltonian()` (`src/CrankNicolsonModes_hamiltonian.cpp:238-240`).

**Suspected documentation mismatch.** `H` is called the "conservative Hamiltonian" in `include/CrankNicolsonModes.hpp:197`, but it is non-Hermitian when a CAP is enabled. Horizon sinks are excluded; the CAP is not.

### 11.2 CN matrices

`update_CN_values()` (`src/CrankNicolsonModes_matrices.cpp:27-40`) forms

\[
A=I+\frac{i\Delta t}{2}H,
\qquad
B=I-\frac{i\Delta t}{2}H.
\]

With sinks,

\[
A_{\rm full}=A+\frac{\Delta t}{4}\Gamma,
\qquad
B_{\rm full}=B-\frac{\Delta t}{4}\Gamma.
\]

The RHS correction appears at `src/CrankNicolsonModes.cpp:170-173`, and the left-operator correction at `src/CrankNicolsonModes_matrices.cpp:105-115`.

### 11.3 Linear solution

`solve_CN()` (`src/CrankNicolsonModes_matrices.cpp:94-215`) normalizes the RHS for amplitude-independent stopping, starts from the old `psi`, and accepts only a physical relative residual at most `1e-9`.

The default preconditioner retains and exactly factorizes each mode's radial tridiagonal block, ignoring inter-mode couplings only in the preconditioner (`include/RadialBlockPreconditioner.hpp:13-19`). With sinks, it is Woodbury-corrected for the low-rank additions. If the radial preconditioner fails, full ILUT is attempted. `tests/parallel_correctness.cpp:16-47` checks the block solve, and lines 109-194 compare radial/ILUT, serial/parallel, and dense reference evolution.

## 12. Circular inspiral, cloud feedback, and accretion

`update_time_dependent_increment()` is called after the field solve only when `evolve_orbit` is true and the startup gate has been crossed (`src/CrankNicolsonModes.cpp:189-195`).

### 12.1 Midpoint cloud state

Accretion rates use

\[
U_{\rm mid}=\frac{U^{n+1}+U^n}{2}
\]

(`src/CrankNicolsonModes_evolution.cpp:19-21`). If absorption is enabled,

\[
\dot M_a
=\mu\,dr\,g_a|v_a^\dagger U_{\rm mid}|^2
\]

(`src/CrankNicolsonModes_absorption.cpp:83-98`). The sink geometry and strength are those built at the field step's endpoint before feedback.

### 12.2 Separation drift

The vacuum circular quadrupole law is (`src/CrankNicolsonModes_evolution.cpp:8-14`)

\[
\dot R_{\rm GW}
=-\frac{64}{5}\frac{M_1M_2(M_1+M_2)}{R^3}.
\]

Cloud power contributes (`src/CrankNicolsonModes_evolution.cpp:45-50`)

\[
\dot R_{\rm cloud}
=-\frac{2R^2}{M_1M_2}p_{\rm bin}.
\]

**Physical assumption.** This relation follows from the circular Newtonian orbital energy

\[
E_{\rm orb}=-\frac{M_1M_2}{2R}
\]

with a power convention in which positive `p_bin` removes orbital energy. Direct separation drift from changing masses is explicitly omitted (`src/CrankNicolsonModes_evolution.cpp:48-49`).

The endpoint feedback weight `ramp` multiplies both drift terms, and forward Euler updates

\[
R\leftarrow R+(\dot R_{\rm GW}+\dot R_{\rm cloud})A(t_{n+1})\Delta t_{\rm orbit}.
\]

`p_bin` here is the rolling value from the preceding diagnostic update, so cloud feedback is one diagnostic step lagged.

### 12.3 Masses and derived geometry

The same `orbit_dt` updates (`src/CrankNicolsonModes_evolution.cpp:52-70`)

\[
M_a\leftarrow M_a+\dot M_a\Delta t_{\rm orbit},
\qquad q\leftarrow M_2/M_1,
\]

then derives COM radii, Kepler frequency, and period:

\[
R_1=\frac{M_2}{M_1+M_2}R,
\quad
R_2=\frac{M_1}{M_1+M_2}R,
\quad
\Omega=\sqrt{\frac{M_1+M_2}{R^3}},
\quad
P=\frac{2\pi}{\Omega}.
\]

When `ramp_factor==1`, the next-step timestep becomes

\[
dt=P/\texttt{steps_per_period}.
\]

The completed step always retains its old timestamp; the new `dt` applies to the next step. `tests/orbital_ramp.cpp:101-117` and `tests/physics.cpp:301-352` protect this behavior. `tests/physics.cpp:354-382` verifies field-plus-black-hole mass balance under ordinary post-gate steps.

### 12.4 Timing caveats recorded, not corrected

**Verified mixed-time behavior.** `H`, `dVdt`, `dpsi_dt`, and field observables are based on the orbit before the feedback update. `R`, `R1`, `R2`, masses, `q`, `Omega`, and `P` are then updated, and the same diagnostic row writes those post-feedback orbital values (`src/CrankNicolsonModes.cpp:184-220`). Thus an evolving-orbit row does not represent a single fully synchronized orbit/Hamiltonian snapshot.

**Suspected crossing-step inconsistency.** If a large step crosses the startup gate, the endpoint sink can act in the CN solve for the whole step, while black-hole accretion uses only the truncated post-gate `orbit_dt`. Existing crossing tests protect separation timing but do not test mass balance across such a step.

## 13. `dVdt`, power, energy, flux, and norm

### 13.1 Implemented `dVdt`

`build_static_perturbation()` sets (`src/CrankNicolsonModes_hamiltonian.cpp:201-219`)

\[
(\partial_tV)_{AB}^{\rm code}
=i(M_B-M_A)s\Omega
\left[(V_1)_{AB}+\lambda(V_2)_{AB}\right].
\]

This is the rotational derivative for fixed radii. It is Hermitian in tested configurations (`tests/physics.cpp:77-104`).

**Suspected inconsistency.** During the companion-potential ramp, the true explicit derivative of `V1 + lambda*V2` also contains `lambda_dot*V2`. The code omits it, and the independent reference test deliberately reproduces the omission (`tests/parallel_correctness.cpp:101`). This is current behavior, not an accidental omission in this audit.

### 13.2 Binary power

`get_binary_power()` returns (`src/CrankNicolsonModes_extraction.cpp:65-68`)

\[
p_{\rm inst}=\operatorname{Re}
\left[dr\sum_{r_i\le r_{\rm ex}}U_i^\dagger(dVdt\,U)_i\right].
\]

`restricted_expectation()` includes radial indices `0..idx_ex` (`src/CrankNicolsonModes_extraction.cpp:8-29`). `p_bin_inst` records this value; `p_bin` records the rolling time average.

### 13.3 Energy and per-mode energy

`get_energy()` returns the real part of the restricted `H` expectation (`src/CrankNicolsonModes_extraction.cpp:105-108`). `get_energy_per_mode()` computes `H*psi`, then sums the real part mode by mode (`lines 110-129`).

**Suspected semantic issue.** When a CAP is active, `H` is non-Hermitian. The reported real expectation is not the expectation of a conservative Hamiltonian, although it remains a defined diagnostic.

`get_dVexpectation_dt()` computes (`src/CrankNicolsonModes_extraction.cpp:70-88`)

\[
\operatorname{Re}\langle\dot U|V|U\rangle_{r\le r_{\rm ex}}.
\]

Despite the output filename `dV_bin_dt.txt`, this is not by itself \(d\langle V\rangle/dt\), which would include twice this term plus \(\langle\partial_tV\rangle\).

### 13.4 Radial fluxes

At extraction index `idx_ex`, the code approximates

\[
r\partial_r(U/r)=\partial_rU-U/r
\]

by centered difference (`src/CrankNicolsonModes_extraction.cpp:31-45`). The per-mode mass flux is

\[
\mathcal F_M=
\operatorname{Im}\left[U^*(\partial_rU-U/r)\right].
\]

`tests/physics.cpp:253-299` verifies the corresponding discrete inner-mass continuity equation and identifies this as `mu*j`.

The energy flux is (`src/CrankNicolsonModes_extraction.cpp:47-63`)

\[
\mathcal F_E=
\frac1\mu\operatorname{Re}
\left[(\partial_tU)^*(\partial_rU-U/r)\right].
\]

`tests/flux_time_alignment.cpp` verifies use of the endpoint state and endpoint derivative.

### 13.5 Norms and cutoff conventions

`get_norms()` uses `psi.segment(off, idx_ex)` (`src/CrankNicolsonModes_extraction.cpp:90-103`), which includes indices `0..idx_ex-1`. `get_energy_per_mode()` uses the same exclusive length. `restricted_expectation()` includes `0..idx_ex`.

**Suspected inconsistency.** Norm/per-mode energy and total energy/power use different inclusion conventions at the extraction node.

**Suspected inconsistency.** `idx_ex=floor(r_ex/dr)` (`src/CrankNicolsonModes.cpp:261-269`) does not subtract `epsilon`, although the grid is `epsilon+i*dr`. For nonzero `epsilon`, the selected node is not generally the node at or immediately below the requested coordinate radius.

## 14. Orbital phase averaging

`RollingPhaseAverage` stores piecewise-linear samples `(time, phase, value)` and integrates the value in time with the trapezoidal rule (`src/RollingPhase.cpp:24-128`). Its window is selected by

\[
\Phi\ge \Phi_{\rm current}-2\pi N_{\rm orbits},
\]

not by a fixed time duration. The first segment crossing the phase cut is linearly truncated. Both time and phase must strictly increase.

In `step()`, the new instantaneous power is passed to

```cpp
power_average.update(t_i, phase, p_bin_last)
```

at `src/CrankNicolsonModes.cpp:202-205`. Counter-rotation does not violate monotonicity because `phase` is unsigned.

`tests/rolling_average.cpp` protects readiness after a phase window, exact averaging of a linear time signal, periodic-signal cancellation, changing-frequency windows, monotonicity checks, reset, and binary serialization.

## 15. Output, checkpoint, and restart

### 15.1 New-run output

`create_output_folder()` makes `Simulation_YYYY_MM_DD_HH_MM_SS`, `psi/`, and `psi_dt/` (`src/Utils.cpp:11-27`). `save_run_info()` copies the exact config and writes `r.txt` and `modes.txt` (`src/CrankNicolsonModes_utils.cpp:11-46`).

At `t=0`, `run()` writes `psi/snapshot_000000.txt`, `psi_dt/snapshot_000000.txt`, and `checkpoint.bin`, but no ordinary diagnostic row (`src/CrankNicolsonModes.cpp:236-258`). Each subsequent `step(n)` writes one row. Therefore a run with `n_steps` time samples writes `n_steps-1` diagnostic rows. `tests/command_line.cpp:85-133` protects file shapes and contents.

When `num_steps=0`, `main()` sets `n_steps=floor(num_periods*P/dt)` (`src/main.cpp:98-102`). Because sample zero is included, the number of evolved intervals and diagnostic rows is one less than that value. `tests/command_line.cpp:156-167` protects this current run-length convention.

### 15.2 Checkpoint payload

`save_checkpoint()` writes magic `GMCKPT01`, step id, grid/mode dimensions, then (`src/CrankNicolsonModes_restart.cpp:139-163`):

```text
t_i, theta, phase, R, M1, M2, q, R1, R2, Omega, P, dt, p_bin,
psi, complete RollingPhaseAverage state
```

The temporary file is renamed over the prior checkpoint only after a complete write (`lines 172-176`). `load_checkpoint()` restores the exact payload (`lines 179-232`).

### 15.3 Diagnostic reconstruction

Without a matching checkpoint, `reconstruct_state_from_diagnostics()` reads time, orbit, mass, and power columns (`src/CrankNicolsonModes_restart.cpp:238-248`). It reconstructs `phase` by the same left-endpoint rule used during evolution:

\[
\Phi_k=\Phi_{k-1}+\Omega_{k-1}(t_k-t_{k-1}),
\]

starting from zero and initial `Omega=0` if a startup ramp exists (`lines 250-262`). It replays instantaneous power into the rolling average. At the selected snapshot it derives `M2=q*M1`, `R1`, and `R2`, restores `Omega` and `p_bin`, and reconstructs post-ramp `P` and `dt` (`lines 267-285`).

The snapshot text has limited precision, so this path is accurate but not bitwise. `tests/restart.cpp:108-195` protects exact checkpoint continuation, crash recovery, approximate diagnostic reconstruction, incomplete snapshot rejection, extension, and mode compatibility.

### 15.4 Hamiltonian rebuild on resume

After restoring the state, `resume()` calls `build_time_dependent_perturbation(t_i)` to rebuild `H`, but saves and restores the checkpoint/diagnostic `Omega` around that call (`src/CrankNicolsonModes_restart.cpp:386-390`). The comment explains that the builder derives circular `Omega` from `R`; restoration preserves the exact frequency used for the next discrete phase advance.

## 16. Orbit-dependent variable catalog

The complete function-by-function dependency list is in `docs/orbit_dependency_map.md`. The following table gives the mathematical role and ownership category.

| Variable | Definition and units | Category |
|---|---|---|
| `R` | Instantaneous binary separation; length | Independent dynamic state when feedback is on; constant configured state otherwise |
| `R1` | `M2/(M1+M2)*R`; length | Derived COM geometry |
| `R2` | `M1/(M1+M2)*R`; length | Derived COM geometry |
| `M1` | Primary mass, initially 1; mass/length | Dynamic state under accretion |
| `M2` | Secondary mass, initially `q`; mass/length | Dynamic state under accretion |
| `q` | `M2/M1`; dimensionless | Configured initially, derived after accretion |
| `theta` | Signed primary azimuth; radians | Numerical orbit coordinate, explicitly integrated |
| `phase` | Unsigned accumulated angle; radians | Numerical averaging bookkeeping |
| `Omega` | Nonnegative circular angular speed used on the next phase advance; inverse time | Ramped/derived dynamic state |
| `Omega_0` | `sqrt((1+q)/R^3)` at construction; inverse time | Immutable initial circular frequency |
| `P` | Initially `2*pi/Omega_0`, later `2*pi/Omega`; time | Derived dynamic state |
| `dt` | Completed/next field timestep; time | Numerical integration control, updated after feedback |
| `t_i` | Current endpoint time; time | Numerical clock |
| `counter` | Selects sign of `theta` evolution; dimensionless boolean | Immutable orbit orientation flag |
| `ramp_factor` | Startup smoothstep used for `Omega` and optionally `V2`; dimensionless | Endpoint-derived policy value |
| `potential_ramp_factor` | Weight of `V2`; dimensionless | Endpoint-derived policy value |
| `initial_period` | Initial Kepler period; time | Immutable schedule anchor |
| `p_bin` | Rolling time average of instantaneous binary power over a phase window; energy/time | Feedback state and diagnostic |
| `dVdt` | Matrix containing the code's rotational explicit potential derivative; energy/time | Derived operator |

## 17. Existing tests that protect a future circular refactor

| Test | Behavior protected |
|---|---|
| `tests/config_and_io.cpp` | Config defaults/values, grid, modes, initial Kepler geometry, MODESBIN mapping, cloud-mass scaling |
| `tests/orbital_ramp.cpp` | Fixed initial-period gates, frequency and potential ramps, crossing-step truncation, zero-ramp behavior, fixed/counter motion |
| `tests/parallel_correctness.cpp` | Independent `H0`, Wigner/multipole potential, current `dVdt`, endpoint CN solve, sink, CAP/mask, preconditioners, thread invariance |
| `tests/fixed_orbit_absorption.cpp` | Absorption ramp, fixed masses/separation/dt under `evolve_orbit=false`, counter `theta/phase`, CN loss balance |
| `tests/spherical_absorption.cpp` | Sink geometry, normalization, source rotation, profile invalidation when radius/mass changes |
| `tests/swave_projector.cpp` | Local s-wave selectivity of the horizon sink |
| `tests/flux_time_alignment.cpp` | Endpoint state/derivative pairing, post-mask diagnostics, initial CAP refresh |
| `tests/physics.cpp` | Hermiticity without CAP, unitarity, Jacobi behavior, co/counter mirror symmetry, flux continuity, Peters inspiral, accretion mass balance, mask and CAP |
| `tests/rolling_average.cpp` | Phase-window semantics, changing `Omega`, monotonic inputs, checkpoint state |
| `tests/restart.cpp` | Exact checkpoint resume, approximate diagnostic reconstruction, crash truncation, extension |
| `tests/command_line.cpp` | End-to-end config/initial state/run output, period-based run length, resume/extension |

The Stage 2 migration should run these existing tests unchanged. New characterization tests should be added only where current behavior is not yet directly protected: mixed pre/post-feedback diagnostics, lagged `p_bin`, the first ramped phase increment, and checkpoint field order/version.

## 18. Suspected inconsistencies and unresolved risks

These findings are recorded for review; none was changed:

1. The constructor may print a coalescence time before `main()` applies `evolve_orbit=false`.
2. The C++ and Python radial grids differ when `epsilon != 0`.
3. Python initial-state scripts accept `[binary].M1`; the C++ solver fixes initial `M1=1`.
4. `dVdt` omits the explicit derivative of the companion-potential ramp.
5. An evolving-orbit diagnostic row mixes pre-feedback field/Hamiltonian quantities with post-feedback orbit scalars.
6. Cloud backreaction uses the preceding rolling `p_bin`, not the instantaneous power computed later in the same step.
7. A gate-crossing step may apply endpoint absorption to the full CN step but accrete only over truncated `orbit_dt`.
8. Total energy/power include the extraction node while per-mode energy/norm exclude it.
9. Extraction-index calculation ignores `epsilon`.
10. `H` is labeled conservative even when it contains a CAP.
11. `dV_bin_dt.txt` stores `Re<dot(psi)|V|psi>`, not the complete time derivative of the potential expectation.
12. The CAP profile denominator in `src/CrankNicolsonModes_outgoing.cpp:56` is `r-L_cap`, where `L_cap=(1-f)L`; existing tests protect sign and support but not this exact intended shape.
13. Several physical input domains (`R`, `q`, ramp duration, mask/CAP parameters) are not validated at config load.
14. Period-based `n_steps` includes the initial sample, so the evolved duration is one timestep shorter than `floor(num_periods*P/dt)*dt`.

The refactor must first preserve these behaviors unless a separately reviewed change explicitly targets one of them.
