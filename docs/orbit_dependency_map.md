# Orbit dependency map

## 1. Purpose and notation

This document maps circular-orbit state to every production function that initializes, reads, derives, modifies, reports, or serializes it. It complements `docs/circular_orbit_implementation_audit.md`; mathematical derivations and timestep chronology are kept there.

Operation codes:

- **I**: initialize from configuration or constructor arguments.
- **R**: read or consume.
- **W**: modify stored state.
- **D**: derive from other state.
- **O**: write to a human-readable diagnostic/output.
- **S**: serialize or deserialize exact checkpoint state.
- **X**: reconstruct approximately from diagnostic output.

`theta` denotes the code member `theta`; `phase` denotes the separate monotonically increasing averaging coordinate. They are intentionally not synonyms.

## 2. High-level dependency matrix

| Function | R | R1/R2 | M1/M2/q | theta | phase | Omega/Omega_0/P | dt/t_i | ramp factors/initial period | p_bin/dVdt |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `load_config()` | I |  | I(q) |  |  |  | I(dt) | I(ramp settings) | I(average window) |
| constructor | I | D/W | I/W | W | W | D/W | I(dt), W(t=0 implicit) | D/I | W(p=0 default) |
| `configure_isolated_black_hole()` | W | W | W | W | W | W |  | W(ramp duration) |  |
| `initialize()` |  |  |  |  |  |  | W(t_i=0) |  | D via perturbation |
| `F1()` / `F2()` |  | R | R(q) |  |  |  |  |  |  |
| `orbital_phase()` |  |  |  | R |  |  |  |  |  |
| `build_static_perturbation()` |  | R | R | R |  | R(Omega), R(counter) |  | R(potential factor) | W(dVdt) |
| `build_time_dependent_perturbation()` | R |  | R |  |  | W/D(Omega) | R(time) | W/D(ramp factors) | D(dVdt) |
| absorber builders |  | R | R | R |  |  | R(t_i) | R(initial period/ramp end) |  |
| `step()` pre-solve | R(print) |  |  | W | W | R(old Omega) | W(t_i) |  |  |
| `step()` post-solve | O |  | O | O | R | O | O | R(ramp end) | R/W/O |
| `update_time_dependent_increment()` | R/W | D/W | R/W/D |  |  | D/W | W(dt) / R(t_i) | R(ramp values) | R(p_bin) |
| extraction/diagnostics |  |  |  |  |  |  |  |  | R(dVdt) |
| `save_checkpoint()` / `load_checkpoint()` | S | S | S | S | S | S | S | average serialized separately | S(p_bin) |
| `reconstruct_state_from_diagnostics()` | X | D/W | X/D/W | X/W | X/W | X/D/W | X/D/W | R | X/W(p_bin) |
| `resume()` | R |  |  |  |  | save/restore around rebuild | R/W |  | rebuild dVdt |

## 3. Variable-by-variable map

### 3.1 `R`: binary separation

**Definition:** instantaneous coordinate distance between the black holes; length in initial-`M1` units.

**Status:** configured independent state for a fixed orbit; forward-Euler dynamic state for circular feedback.

| Operation | Function and evidence |
|---|---|
| I | Constructor argument/member initializer, `src/CrankNicolsonModes.cpp:11-15`; config read at `include/Config.hpp:128-130`; passed by `src/main.cpp:55-71`. |
| R | Constructor derives `initial_period`, `Omega_0`, `R1`, `R2`, and coalescence estimate, `src/CrankNicolsonModes.cpp:14,21,24-25,38`. |
| W | `configure_isolated_black_hole()` assigns `R=abs(x)`, `src/CrankNicolsonModes.cpp:58-74`. |
| R | `build_time_dependent_perturbation()` derives post-ramp circular `Omega=sqrt((M1+M2)/R^3)`, `src/CrankNicolsonModes_hamiltonian.cpp:245-267`. |
| R/W | `update_time_dependent_increment()` evaluates GW/cloud drift and applies `R += dR_dt*orbit_dt`, `src/CrankNicolsonModes_evolution.cpp:45-52`. |
| R | Same function derives `R1`, `R2`, `Omega`, `P`, `dt`, lines 59-70. |
| R | `step()` prints `R` for an evolving orbit and writes `R_orb.txt`, `src/CrankNicolsonModes.cpp:159-162,209`. |
| R | `evolve()` stops an evolving binary when `R<=6`, `src/CrankNicolsonModes.cpp:272-284`. Fixed orbits bypass this check. |
| S | Checkpoint scalar number 4 in save/load order, `src/CrankNicolsonModes_restart.cpp:156,220`. |
| X | Read from `R_orb.txt` and restored at the selected row, `src/CrankNicolsonModes_restart.cpp:241-275`. |
| Test protection | Initial geometry: `tests/config_and_io.cpp:167-185`; fixed orbit: `tests/fixed_orbit_absorption.cpp:41-43`; ramps: `tests/orbital_ramp.cpp`; Peters law: `tests/physics.cpp:301-352`; restart: `tests/restart.cpp`. |

**Downstream influence:** COM radii, gravitational radial profiles, absorber centers, circular frequency/period, feedback drift, adaptive timestep, stop condition, output, and restart.

### 3.2 `R1`: primary COM radius

**Definition:** \(R_1=M_2R/(M_1+M_2)\); nonnegative length.

**Status:** derived geometry, except the isolated-BH/testing APIs can assign it directly.

| Operation | Function and evidence |
|---|---|
| D/W | Constructor uses `q/(1+q)*R`, `src/CrankNicolsonModes.cpp:24`. |
| W | `configure_isolated_black_hole()` sets `R1=abs(x)`, `src/CrankNicolsonModes.cpp:68`. |
| R | `F1()` selects the interior/exterior multipole branch and evaluates powers of `R1`, `src/CrankNicolsonModes_hamiltonian.cpp:12-26`. |
| R | `update_absorption_profile(profile1, R1, epsilon_smear*M1)`, `src/CrankNicolsonModes_absorption.cpp:153-156`. It is also the cache-invalidation radius. |
| D/W | Recomputed after feedback, `src/CrankNicolsonModes_evolution.cpp:59-62`. |
| S | Checkpoint scalar number 8, `src/CrankNicolsonModes_restart.cpp:156,224`. |
| X | Derived from reconstructed masses and `R`, `src/CrankNicolsonModes_restart.cpp:267-275`. |
| Test protection | COM relation: `tests/config_and_io.cpp:183`; evolving relation: `tests/physics.cpp:342`; absorber geometry: `tests/spherical_absorption.cpp`; fixed value: `tests/fixed_orbit_absorption.cpp:41-43`. |

**Downstream influence:** primary multipole break radius and strength; primary horizon-kernel center; initial-state displacement in Python (`initial_states/get_initial_state_COM_fast.py:69-71,573-579`).

### 3.3 `R2`: secondary COM radius

**Definition:** \(R_2=M_1R/(M_1+M_2)\); nonnegative length.

**Status:** derived geometry.

| Operation | Function and evidence |
|---|---|
| D/W | Constructor uses `1/(1+q)*R`, `src/CrankNicolsonModes.cpp:25`. |
| W | Isolated-BH mode sets it to zero, `src/CrankNicolsonModes.cpp:69`. |
| R | `F2()` selects branch and evaluates powers of `R2`, `src/CrankNicolsonModes_hamiltonian.cpp:28-42`. |
| R | Secondary absorber profile center, `src/CrankNicolsonModes_absorption.cpp:156`. |
| D/W | Recomputed after feedback, `src/CrankNicolsonModes_evolution.cpp:62`. |
| S | Checkpoint scalar number 9, `src/CrankNicolsonModes_restart.cpp:156,225`. |
| X | Derived from reconstructed masses and `R`, `src/CrankNicolsonModes_restart.cpp:275`. |
| Test protection | Initial COM relation: `tests/config_and_io.cpp:184`; fixed value: `tests/fixed_orbit_absorption.cpp:41-43`; independent sink geometry: `tests/SphericalAbsorptionReference.hpp`. |

**Downstream influence:** secondary multipole break radius and secondary horizon-kernel center; Python secondary displacement is represented as negative `R2` on the same `+x` translation axis (`initial_states/get_initial_state_COM_fast.py:580-586`).

### 3.4 `M1`: primary mass

**Definition:** primary black-hole mass; initialized to 1; mass/length units.

**Status:** independent dynamic state under accretion.

| Operation | Function and evidence |
|---|---|
| I/W | Hard-coded to `1.0` in constructor, `src/CrankNicolsonModes.cpp:30`; isolated mode also sets 1, line 66. |
| R | Potential common factor `M1+M2`, `src/CrankNicolsonModes_hamiltonian.cpp:204-205`. Combined with `q` fractions to recover physical source masses. |
| R | Absorber width `epsilon_smear*M1` and area `16*pi*M1^2`, `src/CrankNicolsonModes_absorption.cpp:155,170`. |
| R/W | Accretion update `M1 += Mdot1*orbit_dt`, `src/CrankNicolsonModes_evolution.cpp:23-30,54`. |
| R | GW drift, cloud drift denominator, mass total, `q`, COM radii, and Kepler frequency, lines 45-64. |
| O | Written to `M_1.txt`, `src/CrankNicolsonModes.cpp:212`. |
| S | Checkpoint scalar number 5, `src/CrankNicolsonModes_restart.cpp:156,221`. |
| X | Read from `M_1.txt`, then used to derive `M2=q*M1`, `src/CrankNicolsonModes_restart.cpp:245,271-275`. |
| Test protection | Initial unit convention: `tests/config_and_io.cpp:185`; mass balance: `tests/physics.cpp:354-382`; CLI accretion: `tests/command_line.cpp:128-133`. |

**Downstream influence:** physical primary potential mass, primary absorber size/strength, GW losses, cloud backreaction conversion, COM geometry, frequency/period, and mass diagnostic.

### 3.5 `M2`: secondary mass

**Definition:** secondary black-hole mass; initialized to `q`; mass/length units.

**Status:** independent dynamic state under accretion after initially deriving it from configured `q` and `M1=1`.

| Operation | Function and evidence |
|---|---|
| I/W | Constructor `M2=q`, `src/CrankNicolsonModes.cpp:31`; isolated mode sets zero, line 67. |
| R | Potential common factor, absorber width/area, both separation drifts, COM geometry, and Kepler frequency. See `src/CrankNicolsonModes_hamiltonian.cpp:204-205`, `CrankNicolsonModes_absorption.cpp:156,171`, and `CrankNicolsonModes_evolution.cpp:45-64`. |
| W | `M2 += Mdot2*orbit_dt`, `src/CrankNicolsonModes_evolution.cpp:55`. |
| S | Checkpoint scalar number 6, `src/CrankNicolsonModes_restart.cpp:156,222`. |
| X | Derived as `q*M1`, `src/CrankNicolsonModes_restart.cpp:273`. |
| Test protection | Initial convention: `tests/config_and_io.cpp:185`; accretion balance: `tests/physics.cpp:354-382`. |

There is no separate `M_2.txt`; diagnostic reconstruction derives `M2` from `M1` and `q`.

### 3.6 `q`: mass ratio

**Definition:** \(q=M_2/M_1\); dimensionless.

**Status:** configured independent value initially; derived after accretion.

| Operation | Function and evidence |
|---|---|
| I | Config/default at `include/Config.hpp:128`; constructor member initializer `src/CrankNicolsonModes.cpp:11-15`. |
| R | `initial_period`, `Omega_0`, COM radii, and coalescence estimate, `src/CrankNicolsonModes.cpp:14,21,24-25,33-39`. |
| R | `F1` and `F2` mass fractions, `src/CrankNicolsonModes_hamiltonian.cpp:18,22,34,38`. |
| W | Isolated mode sets zero, `src/CrankNicolsonModes.cpp:67`. |
| D/W | After accretion, `q=M2/M1`, `src/CrankNicolsonModes_evolution.cpp:57`. |
| O | Written to `q.txt`, `src/CrankNicolsonModes.cpp:213`. |
| S | Checkpoint scalar number 7, `src/CrankNicolsonModes_restart.cpp:156,223`. |
| X | Read from `q.txt`; used to derive `M2`, `src/CrankNicolsonModes_restart.cpp:246,272-275`. |
| Test protection | Config/initial geometry: `tests/config_and_io.cpp`; `q=M2/M1` after accretion: `tests/physics.cpp:378`; restart: `tests/restart.cpp`. |

### 3.7 `theta`: signed orbital angle

**Definition:** primary azimuth in the COM equatorial plane; radians.

**Status:** numerical orbit coordinate, explicitly advanced.

| Operation | Function and evidence |
|---|---|
| I/W | Constructor sets zero, `src/CrankNicolsonModes.cpp:27`; isolated mode sets 0 or pi from the sign of `x`, line 70. |
| W | `step()` applies `theta += (counter ? -1 : 1)*Omega*dt` before rebuilding `Omega`, `src/CrankNicolsonModes.cpp:150-153`. |
| R | `orbital_phase(m_star)=exp(i*m_star*theta)`, `src/CrankNicolsonModes_hamiltonian.cpp:106-109`. |
| R | Potential and absorber projector assembly through `orbital_phase()`, `src/CrankNicolsonModes_hamiltonian.cpp:186-219`; `src/CrankNicolsonModes_absorption.cpp:127-150`. |
| O | Written to `theta_orb.txt`, `src/CrankNicolsonModes.cpp:210`. |
| S | Checkpoint scalar number 2, `src/CrankNicolsonModes_restart.cpp:156,218`. |
| X | Read from `theta_orb.txt`, `src/CrankNicolsonModes_restart.cpp:243,270`. |
| Test protection | Counter sign and separation from `phase`: `tests/orbital_ramp.cpp:144-155`, `tests/fixed_orbit_absorption.cpp:60-64`; mirror symmetry: `tests/physics.cpp:221-251`; absorber rotation: `tests/spherical_absorption.cpp:145-152`. |

**Downstream influence:** all complex potential phases and both horizon-projector orientations. It is not used to select the rolling-average window.

### 3.8 `phase`: monotonic averaging phase

**Definition:** unsigned accumulated circular angle; radians.

**Status:** numerical bookkeeping, not a body coordinate.

| Operation | Function and evidence |
|---|---|
| I/W | Constructor and isolated mode set zero, `src/CrankNicolsonModes.cpp:28,71`. |
| W | `step()` applies `phase += Omega*dt` with old `Omega`, `src/CrankNicolsonModes.cpp:153`. |
| R | Passed to `power_average.update(t_i, phase, p_bin_last)`, `src/CrankNicolsonModes.cpp:202-205`. |
| S | Checkpoint scalar number 3; rolling-average internal phase state is also serialized, `src/CrankNicolsonModes_restart.cpp:156,163,219,231`. |
| X | Reset to zero and rebuilt from prior `Omega`/time rows, while replaying instantaneous power, `src/CrankNicolsonModes_restart.cpp:250-262`. |
| Test protection | Strict monotonic requirement and window behavior: `tests/rolling_average.cpp`; counter behavior: `tests/orbital_ramp.cpp:152-154`; exact restart: `tests/restart.cpp`. |

`phase` has no direct influence on `V`, `H`, absorption geometry, or the CN solve. It influences future feedback indirectly through the phase-windowed `p_bin`.

### 3.9 `Omega`: current circular angular speed

**Definition:** nonnegative angular speed used for the next discrete `theta` and `phase` increment; inverse time.

**Status:** ramped during startup; Kepler-derived after startup and after feedback.

| Operation | Function and evidence |
|---|---|
| I/W | Constructor sets zero, `src/CrankNicolsonModes.cpp:20`; isolated mode sets zero, line 72. |
| R | `step()` uses the entry value to advance `theta` and `phase`, lines 150-153. |
| W/D | `build_time_dependent_perturbation()` assigns zero for isolated mode, `Omega_0*ramp_factor` during startup, or `sqrt((M1+M2)/R^3)` afterwards, `src/CrankNicolsonModes_hamiltonian.cpp:245-261`. |
| R | `build_static_perturbation()` forms rotational `dVdt`, `src/CrankNicolsonModes_hamiltonian.cpp:189-211`. |
| D/W | Feedback recomputes Kepler `Omega`, `src/CrankNicolsonModes_evolution.cpp:64`. |
| O | Written to `Omega_orb.txt`, `src/CrankNicolsonModes.cpp:211`. |
| S | Checkpoint scalar number 10, `src/CrankNicolsonModes_restart.cpp:156,226`. |
| X | Read from `Omega_orb.txt`; previous values reconstruct `phase`; final value restores current `Omega`, `src/CrankNicolsonModes_restart.cpp:244,250-260,276`. |
| Resume special case | Saved before `build_time_dependent_perturbation(t_i)` and restored after it, `src/CrankNicolsonModes_restart.cpp:386-390`. |
| Test protection | Frequency ramp: `tests/orbital_ramp.cpp:119-142`; Kepler update: `tests/physics.cpp:341`; fixed/counter phase: `tests/fixed_orbit_absorption.cpp:60-64`; restart: `tests/restart.cpp`. |

### 3.10 `Omega_0`: initial Kepler frequency

**Definition:** \(\Omega_0=\sqrt{(1+q)/R^3}\) at construction; inverse time.

**Status:** derived initial constant; not serialized because the config reconstructs it.

| Operation | Function and evidence |
|---|---|
| D/W | Constructor, `src/CrankNicolsonModes.cpp:21`. |
| R | Startup frequency ramp, `src/CrankNicolsonModes_hamiltonian.cpp:256`. |
| R | Diagnostic-only restart chooses initial previous frequency as zero or `Omega_0`, `src/CrankNicolsonModes_restart.cpp:250-252`. |
| Test protection | `tests/orbital_ramp.cpp:133-134,150`; period/geometry tests indirectly protect it. |

It is based on initial `q` and `R`; later accretion/inspiral does not change it.

### 3.11 `P`: current circular period

**Definition:** initially `initial_period`; after feedback \(2\pi/\Omega\); time units.

**Status:** derived state and run-length/adaptive-step input.

| Operation | Function and evidence |
|---|---|
| D/W | Constructor `P=initial_period`, `src/CrankNicolsonModes.cpp:22`. |
| R | `main()` computes period-based `num_steps=floor(num_periods*P/dt)`, `src/main.cpp:98-102`. At this point it is the initial period. |
| D/W | Feedback `P=2*pi/Omega`, `src/CrankNicolsonModes_evolution.cpp:66`. |
| R | Feedback derives next `dt=P/steps_per_period`, line 70. |
| S | Checkpoint scalar number 11, `src/CrankNicolsonModes_restart.cpp:156,227`. |
| X | Diagnostic reconstruction derives it from final `Omega` after the gate, lines 279-284. |
| Test protection | Initial period: `tests/config_and_io.cpp:182`; run length: `tests/command_line.cpp:156-167`; adaptive period: `tests/physics.cpp:343`; fixed period: `tests/fixed_orbit_absorption.cpp:41-43`. |

### 3.12 `dt`: field timestep and next orbit increment

**Definition:** completed/next numerical timestep; time units.

**Status:** configured integration control during startup/fixed evolution; derived after fully ramped feedback.

| Operation | Function and evidence |
|---|---|
| I | Config `simulation.dt`, `include/Config.hpp:109`; constructor member initializer, `src/CrankNicolsonModes.cpp:11-15`. |
| R | `step()` advances `t_i`, `theta`, and `phase`; builds CN matrices and sink terms, `src/CrankNicolsonModes.cpp:150-173`. |
| R | `update_CN_values()` uses `dt/2`, `src/CrankNicolsonModes_matrices.cpp:27-40`; `solve_CN()` uses `dt/4` for sinks, lines 105-155. |
| W/D | Feedback sets the next `dt=P/steps_per_period` when `ramp_factor==1`, `src/CrankNicolsonModes_evolution.cpp:68-71`. |
| R | `main()` uses initial `dt` to translate `num_periods` to sample count, `src/main.cpp:98-102`. |
| S | Checkpoint scalar number 12, `src/CrankNicolsonModes_restart.cpp:156,228`. |
| X | Reconstructed as `P/steps_per_period` only after the gate, `src/CrankNicolsonModes_restart.cpp:279-285`. |
| Test protection | Input validation: `tests/config_and_io.cpp:188-214`; adaptive timing and crossing: `tests/orbital_ramp.cpp`; Kepler stepping: `tests/physics.cpp:343`; restart: `tests/restart.cpp`. |

`dt` is numerical integration state, not physical orbit geometry. A future `OrbitState` should not own it.

### 3.13 `t_i`: current simulation time

**Definition:** time attached to the current endpoint state; time units.

**Status:** numerical clock.

| Operation | Function and evidence |
|---|---|
| W | `initialize()` and `run()` set zero, `src/CrankNicolsonModes.cpp:85-87,241-244`. |
| W | `step()` increments by current `dt` before orbit/ramp/H construction, lines 150-166. |
| R | Ramp evaluation, sink activation, feedback gate/weight, rolling average, output timestamp, and checkpoint. |
| O | Written to `t.txt`, `src/CrankNicolsonModes.cpp:214`. |
| S | Checkpoint scalar number 1, `src/CrankNicolsonModes_restart.cpp:156,217`. |
| X | Read from `t.txt`, final row restored, `src/CrankNicolsonModes_restart.cpp:241,268`. |
| Test protection | Endpoint alignment: `tests/flux_time_alignment.cpp:94-106`; crossing timestamp: `tests/orbital_ramp.cpp:116`; CLI monotonic time: `tests/command_line.cpp:128-131`; restart. |

### 3.14 `counter`: rotation-sense flag

**Definition:** `false` means increasing `theta`; `true` means decreasing `theta`; dimensionless immutable configuration.

| Operation | Function and evidence |
|---|---|
| I | Config and constructor, `include/Config.hpp:130`; `src/main.cpp:66`; `src/CrankNicolsonModes.cpp:13`. |
| R | `step()` sign of `theta`, `src/CrankNicolsonModes.cpp:152`. |
| R | `build_static_perturbation()` sign of rotational `dVdt`, `src/CrankNicolsonModes_hamiltonian.cpp:189-211`. |
| Not used | `phase` remains positive regardless of `counter`. `Omega` remains nonnegative. |
| Test protection | `tests/orbital_ramp.cpp:144-155`; `tests/fixed_orbit_absorption.cpp:60-64`; `tests/physics.cpp:221-251`; `tests/spherical_absorption.cpp`. |

`counter` is not serialized separately; resume reconstructs it from the copied config.

### 3.15 `ramp_factor`

**Definition:** endpoint startup smoothstep \(S(t/t_{ramp})\), or 1 after the gate; dimensionless.

**Status:** derived policy/bookkeeping value stored as a member.

| Operation | Function and evidence |
|---|---|
| W/D | Reset to 1, then possibly assigned smoothstep by `build_time_dependent_perturbation()`, `src/CrankNicolsonModes_hamiltonian.cpp:245-261`. |
| R | Sets ramped `Omega`; copied into `potential_ramp_factor` when configured, line 262. |
| R | `update_time_dependent_increment()` checks equality to 1 before adapting `dt`, `src/CrankNicolsonModes_evolution.cpp:68-71`. |
| Test protection | `tests/orbital_ramp.cpp`. |

It is not checkpointed. It is recomputed whenever `H` is rebuilt from `t_i`.

### 3.16 `potential_ramp_factor`

**Definition:** scalar companion-potential weight `ramp_factor` or 1; dimensionless.

**Status:** derived policy value.

| Operation | Function and evidence |
|---|---|
| W/D | `build_time_dependent_perturbation()`, `src/CrankNicolsonModes_hamiltonian.cpp:262`. |
| R | `build_static_perturbation()` sets `V=V1+factor*V2` and includes the same factor in current angular `dVdt`, lines 201-219. |
| Test protection | Both policy settings: `tests/orbital_ramp.cpp:119-142`; independent current formula: `tests/parallel_correctness.cpp:96-102,147-170`. |

It is not checkpointed and is recomputed from config plus time.

### 3.17 `initial_period`

**Definition:** immutable initial circular Kepler period

\[
P_0=2\pi\sqrt{R(0)^3/(1+q(0))}.
\]

**Status:** constructor-derived constant and schedule anchor.

| Operation | Function and evidence |
|---|---|
| D/I | Const member initialized in constructor initializer list, `src/CrankNicolsonModes.cpp:14`; declared `include/CrankNicolsonModes.hpp:242-246`. |
| R | `ramp_end_time()=ramp_periods*initial_period`. |
| R | Absorption and backreaction ramps last one additional initial period, `src/CrankNicolsonModes_absorption.cpp:160-168`; `src/CrankNicolsonModes_evolution.cpp:32-43`. |
| Test protection | Fixed gate despite changing current `P`: `tests/orbital_ramp.cpp:55-99`. |

It is not checkpointed because the copied configuration reconstructs it.

### 3.18 `p_bin`: rolling binary power

**Definition:** time average of instantaneous `get_binary_power()` over the last configured number of accumulated-phase orbits; energy/time units.

**Status:** feedback state plus diagnostic.

| Operation | Function and evidence |
|---|---|
| Initial value | In-class member initializer `0.0`, `include/CrankNicolsonModes.hpp:381`. |
| R | `update_time_dependent_increment()` computes cloud-driven `dR_dt`, `src/CrankNicolsonModes_evolution.cpp:45-50`. This occurs before the current step's power is calculated. |
| W | After feedback, `step()` assigns `power_average.update(t_i, phase, p_bin_last)`, `src/CrankNicolsonModes.cpp:202-205`. Isolated mode forces zero. |
| O | Written to `p_bin.txt`; instantaneous input written separately to `p_bin_inst.txt`, lines 218-219. |
| S | Checkpoint scalar number 13; full averaging history separately serialized, `src/CrankNicolsonModes_restart.cpp:156,163,229,231`. |
| X | Averaging history replayed from `p_bin_inst.txt`; final stored average read from `p_bin.txt`, `src/CrankNicolsonModes_restart.cpp:247-259,277`. |
| Test protection | `tests/rolling_average.cpp`; exact restart `tests/restart.cpp`; output shape `tests/command_line.cpp`. No test directly asserts the one-step feedback lag. |

### 3.19 `dVdt`: explicit potential derivative matrix

**Definition:** current implementation of the rotational explicit derivative of the active binary potential; energy/time operator.

**Status:** derived sparse operator sharing the `V` pattern.

| Operation | Function and evidence |
|---|---|
| Pattern | Allocated by `build_sparse_element_caches()`, `src/CrankNicolsonModes_hamiltonian.cpp:111-149`. |
| W/D | Refilled by `build_static_perturbation()`, lines 176-222, using `theta`, `counter`, `Omega`, source matrices, and `potential_ramp_factor`. |
| R | `get_binary_power()` passes it to `restricted_expectation()`, `src/CrankNicolsonModes_extraction.cpp:65-68`. |
| Indirect O | Its expectation is written to `p_bin_inst.txt` and averaged into `p_bin.txt`. The matrix itself is not serialized. |
| Rebuild | Initialization, every step, initial CAP refresh in `run()`, and restored-H rebuild in `resume()`. |
| Test protection | Independent exact formula: `tests/parallel_correctness.cpp:49-106,167-170`; Hermiticity: `tests/physics.cpp:77-104`. |

Current formula:

\[
dVdt=i(M_B-M_A)s\Omega(V_1+\lambda V_2).
\]

It intentionally maps current behavior, including omission of `lambda_dot*V2` during startup.

### 3.20 Orbit-control and absorption-control values

| Value | Initialized | Read/modified by | Effect and persistence |
|---|---|---|---|
| `evolve_orbit` | Public default `true`, `include/CrankNicolsonModes.hpp:137`; overwritten from config by `src/main.cpp:76` | `step()` gates `update_time_dependent_increment()`; `evolve()` gates `R<=6`; diagnostic reconstruction gates adaptive `P/dt` | Selects post-solve feedback only. Not checkpointed; reconstructed from copied config. |
| `absorption` | Constructor from config, `src/main.cpp:67`; stored at `include/CrankNicolsonModes.hpp:333` | Sparse-cache setup, perturbation builder, absorption builder, accretion calculation | Enables sink construction/evolution. It does not by itself make masses grow when `evolve_orbit=false`. Copied config, not a checkpoint scalar. |
| `ramp_periods` | Constructor/config, `src/CrankNicolsonModes.cpp:13`; `include/Config.hpp:115` | `ramp_end_time()` and all four activation schedules | Multiplies immutable `initial_period`. Isolated mode sets zero. Copied config. |
| `ramp_secondary_potential` | Public default `true`, overwritten by `src/main.cpp:77` | `build_time_dependent_perturbation()` selects `potential_ramp_factor`, `src/CrankNicolsonModes_hamiltonian.cpp:262` | Does not disable the frequency ramp. Copied config. |
| `steps_per_period` | Public default/config via `src/main.cpp:78` | Feedback derives next `dt`; diagnostic restart reconstructs `dt` | Only relevant when circular feedback is active after startup. Copied config. |
| `orbits_average` / `power_average` | Default 2; reset through `set_power_average_orbits()`, `src/CrankNicolsonModes.cpp:52-56` | `step()` updates; checkpoint saves full average; diagnostic restart replays it | The configured window is phase width `2*pi*n_orbits`, while values are time-averaged. Exact history is checkpointed. |
| `epsilon_smear` | Constructor/config | Absorber widths `epsilon_smear*M1/M2` | Copied config; profile cache also records current physical width. |
| `sink_active` | Initially false; set by `assemble_hamiltonian(include_absorption)`, `src/CrankNicolsonModes_hamiltonian.cpp:224-243` | RHS, CN operator, `schrodinger_rhs()`, public status/strength queries | True at/after startup gate when absorption configured, even if smoothstep strength is exactly zero at the gate. Not serialized; rebuilt. |
| `Ah1`, `Ah2` | Zero members; rebuilt by `build_spherical_absorption()` | Converted to sink strengths `Ah*dr` | `16*pi*M_a^2*absorption_ramp`; not serialized. |

## 4. Derived operator and subsystem dependencies

| Derived object/subsystem | Direct orbit inputs | Notes |
|---|---|---|
| `radial_f1` | `R1`, `q` | Rebuilt every perturbation call. |
| `radial_f2` | `R2`, `q` | Cached without secondary parity; parity applied per mode pair. |
| complex phase cache | `theta` | Contains `exp(i*m*theta)` for `m=-2*lmax..2*lmax`. |
| `V1` | `M1+M2`, `q`, `R1`, `theta` | Algebraically gives source mass `M1`. |
| `V2` | `M1+M2`, `q`, `R2`, `theta` | Algebraically gives source mass `M2`; includes `(-1)^m_star`. |
| `V` | `V1`, `V2`, `potential_ramp_factor` | Active field potential. |
| `dVdt` | `V1`, `V2`, `theta`, `counter`, `Omega`, potential ramp weight | Rotational term only. |
| primary absorber profile | `R1`, `M1`, `epsilon_smear` | Radial Legendre moments; invalidated by radius or width change. |
| secondary absorber profile | `R2`, `M2`, `epsilon_smear` | Same. |
| sink projectors | profile moments, `theta`, secondary parity | Rotate every rebuild even if radial profiles are cached. |
| sink strengths | `M1`, `M2`, `t_i`, `initial_period`, ramp end | `16*pi*M_a^2*ramp*dr`. |
| `H` | `H0`, `V`, CAP, absorption-active flag | Rank-one sinks remain outside sparse `H`. |
| `A`, `B` | `H`, `dt` | Updated every field step. |
| full CN operator | `A`, `B`, sinks, `dt` | Sink factors are `dt/4`. |
| GW separation drift | `R`, `M1`, `M2`, endpoint feedback ramp | Circular quadrupole law. |
| cloud separation drift | `R`, `M1`, `M2`, prior `p_bin`, endpoint feedback ramp | Uses rolling power from preceding update. |
| mass accretion | `M1`, `M2`, sink/projector, midpoint field, `orbit_dt` | Only when `evolve_orbit=true`. |
| rolling power | `t_i`, `phase`, instantaneous `dVdt` expectation | Window selected in phase; integral averaged in time. |

## 5. Serialization and output map

| Quantity | Text file | Checkpoint | Diagnostic reconstruction |
|---|---|---|---|
| `t_i` | `t.txt` | yes | direct |
| `R` | `R_orb.txt` | yes | direct |
| `R1`, `R2` | no | yes | derived from masses and `R` |
| `M1` | `M_1.txt` | yes | direct |
| `M2` | no | yes | `q*M1` |
| `q` | `q.txt` | yes | direct |
| `theta` | `theta_orb.txt` | yes | direct |
| `phase` | no | yes | replay old `Omega*delta_t` |
| `Omega` | `Omega_orb.txt` | yes | direct; history also rebuilds phase |
| `Omega_0` | no | no | constructor/config |
| `P` | no | yes | initial value or `2*pi/Omega` after gate |
| `dt` | no | yes | configured initial value or `P/steps_per_period` after gate |
| `counter` | copied config | no | copied config |
| ramp factors | no | no | recomputed from `t_i`, initial period, and config |
| `initial_period` | indirectly copied config | no | constructor/config |
| `p_bin` | `p_bin.txt` | yes | direct final value |
| instantaneous power | `p_bin_inst.txt` | no | replayed into rolling average |
| rolling-average segments | no | yes | reconstructed by replay |
| `dVdt` matrix | no | no | rebuilt from restored orbit |

Every file in `diagnostic_files()` participates in common-row crash recovery (`src/CrankNicolsonModes_restart.cpp:109-137,318-325`). Adding a new diagnostic to that list changes the minimum complete-row calculation and restart truncation behavior.

## 6. Mutation order within one evolving-orbit step

| Order | Mutated value | Inputs used | Visible consequence |
|---:|---|---|---|
| 1 | `psi_last` | old `psi` | midpoint accretion state later |
| 2 | `t_i` | old `t_i`, old `dt` | endpoint timestamp |
| 3 | `theta` | old `Omega`, old `dt`, `counter` | endpoint source angle |
| 4 | `phase` | old `Omega`, old `dt` | endpoint averaging coordinate |
| 5 | `ramp_factor`, `potential_ramp_factor`, `Omega` | endpoint `t_i`, pre-feedback `R/M` | endpoint `V`, `dVdt`, sinks, `H` |
| 6 | `psi` | endpoint `H`, old `dt`, old `psi` | solved endpoint field |
| 7 | `psi` again if mask | damping mask | post-mask diagnostic state |
| 8 | `R`, `M1`, `M2` | prior `p_bin`, midpoint field, endpoint sink, endpoint ramp, `orbit_dt` | post-feedback orbit |
| 9 | `q`, `R1`, `R2`, `Omega`, `P`, next `dt` | newly updated `R/M` | next-step orbit/control |
| 10 | `p_bin` | new instantaneous power, `t_i`, `phase` | next-step cloud feedback |
| 11 | checkpoint/output | mixed state described above | restart boundary |

This order is directly encoded in `src/CrankNicolsonModes.cpp:145-233` and must be preserved in Stage 2 unless a later, separately approved numerical change targets it.

## 7. Known dependency risks for refactoring

1. Moving `Omega` evaluation before the phase update would change the startup trajectory.
2. Making `theta` and `phase` one field would break counter-rotation and `RollingPhaseAverage` monotonicity.
3. Treating `dt` as orbit geometry would obscure that the feedback update changes only the next timestep.
4. Recomputing orbit state after feedback before field diagnostics would change mixed-time output and restart bytes.
5. Replacing `initial_period` by current `P` in any gate would reintroduce the moving-gate behavior explicitly rejected by `tests/orbital_ramp.cpp`.
6. Serializing only `OrbitState` without the rolling-average segments would break bitwise checkpoint continuation.
7. Deriving `phase` from `theta` would fail for `counter=true` and would not reproduce current diagnostic-only restart.
8. Passing post-feedback `R1/R2` into the already-completed step's sink or potential would change field evolution and mass balance.
