# OrbitState migration tracker

## Baseline and scope

- Base commit: `c746321bbbc3fbddef3fb65de729f7bd6345c581` (`Initial commit`).
- Working branch: `orbit-state-change-01`.
- Current authorized scope: Change 01 only. `CrankNicolsonModes` does not yet include or use `OrbitState`.
- Authoritative dependency reference: `docs/orbit_dependency_map.md`.

## Change history

| Change | Commit | Status | Contents |
|---|---|---|---|
| 01: introduce `OrbitState` | `309df62` | Complete | Minimal state/getters, focused unit tests, CMake test registration. |
| 02: centralize ownership | Not authorized | Pending | Replace solver orbital scalar storage; migrate updates, diagnostics, and restart. |
| 03: explicit operator inputs | Not authorized | Pending | Pass `const OrbitState&` to potential/absorption builders; separate ramped angular-rate update from assembly. |
| 04: verify and clean up | Not authorized | Pending | Dependency sweep, compatibility verification, and removal of temporary migration code. |

## State mapping

Status in this table describes the solver, not the availability of the new getter.

| Legacy value | Current consumers | Final source | Status |
|---|---|---|---|
| `M1`, `M2` | Constructor; isolated setup; potential and absorption assembly; circular feedback/accretion; diagnostics; checkpoint and reconstruction | `orbit.M1`, `orbit.M2` | Pending Change 02/03 |
| `R` | Constructor; isolated setup; ramp/post-ramp frequency; circular feedback; coalescence stop; diagnostics; checkpoint and reconstruction | `orbit.R` | Pending Change 02/03 |
| `theta` | Constructor; isolated setup; discrete phase advance; `orbital_phase()`; potential/sink assembly; diagnostics; checkpoint and reconstruction | `orbit.theta` | Pending Change 02/03 |
| `phase` | Constructor; discrete positive phase advance; rolling average; checkpoint and diagnostic reconstruction | `orbit.phase` | Pending Change 02 |
| radial velocity (previously transient `dR_dt`) | `update_time_dependent_increment()` | `orbit.R_dot`; zero for fixed circular and isolated states | Pending Change 02 |
| signed angular velocity (previously `counter` plus unsigned `Omega`) | `step()`; startup ramp; `dVdt`; diagnostics; checkpoint and reconstruction | `orbit.theta_dot`; `orbit.angular_speed()` where an unsigned value is required | Pending Change 02/03 |
| `q` | Initial geometry/coalescence estimate; `F1()`/`F2()`; feedback; diagnostics; checkpoint and reconstruction | `orbit.q()` | Getter complete; consumers pending |
| `R1`, `R2` | `F1()`/`F2()`; absorption profiles; feedback; checkpoint and reconstruction | `orbit.R1()`, `orbit.R2()` for ordinary binaries only | Getters complete; consumers pending |
| `Omega` | Old-frequency phase advance; ramp update; `dVdt`; feedback; diagnostics; checkpoint and reconstruction | `orbit.angular_speed()` with the sign retained in `orbit.theta_dot` | Getter complete; consumers pending |
| `P` | Initial run length; adaptive `dt`; checkpoint and reconstruction | `orbit.circular_kepler_period()` for an ordinary post-ramp circular state; immutable `initial_period` remains outside the state | Getter complete; consumers pending |
| `Omega_0` | Startup frequency ramp and diagnostic reconstruction | Startup/reference scalar outside `OrbitState` | Intentionally unchanged |
| `counter` | Selects rotation sign in phase advance and `dVdt` | Configuration input outside `OrbitState`; instantaneous rotation sign lives in `orbit.theta_dot` | Consumer migration pending Change 02/03 |

The displaced isolated-black-hole path is intentionally exceptional: it stores
the displacement magnitude in `orbit.R`, but source builders must use that value
directly and must not use the ordinary-binary `R1()` getter.

## Function migration status

| File/function | Planned Change | Status |
|---|---:|---|
| `include/OrbitState.hpp` | 01 | Complete |
| `tests/orbit_state.cpp` | 01 | Complete |
| `CrankNicolsonModes` constructor | 02 | Pending |
| `configure_isolated_black_hole()` | 02 | Pending |
| `step()` discrete `theta`/`phase` advance and diagnostic writes | 02 | Pending |
| `run()` / `evolve()` orbit reads | 02 | Pending |
| `update_time_dependent_increment()` | 02 | Pending |
| `save_checkpoint()` / `load_checkpoint()` | 02 | Pending |
| `reconstruct_state_from_diagnostics()` / `resume()` | 02 | Pending |
| `src/main.cpp` period-based step count | 02 | Pending |
| Direct orbital accesses in existing tests/benchmarks | 02 | Pending |
| `F1()` / `F2()` / `orbital_phase()` | 03 | Pending |
| `build_static_perturbation()` | 03 | Pending |
| `build_time_dependent_perturbation()` ramped angular-rate update | 03 | Pending |
| `build_swave_projector()` / `build_spherical_absorption()` | 03 | Pending |
| Hidden-dependency and duplicate-storage sweep | 04 | Pending |

## Verification

| Commit/state | Command | Result |
|---|---|---|
| `309df62` plus getter arithmetic review | `cmake --build build-codex-verify --target orbit_state --config Release` | Passed |
| Same | `ctest --test-dir build-codex-verify -R '^orbit_state$' --output-on-failure -C Release` | 1/1 passed |
| Same | `ctest --test-dir build-codex-verify -L unit --output-on-failure -C Release` | 3/3 passed |

The focused tests cover total mass and `q`, the `+R1`/`-R2` COM convention,
circular Kepler frequency/period, signed `theta_dot` versus unsigned angular
speed and positive `phase`, and the explicit isolated-source geometry exception.
No baseline simulation was run.

## Behavior and open risks

- Behavioral differences in Change 01: none. The solver does not consume the new type.
- Getter arithmetic follows the current solver's `(mass / total) * R` and
  `sqrt(total / pow(R, 3))` operation order to avoid gratuitous rounding changes.
- The getters assume the existing physical domain `M1 > 0`, `M1 + M2 > 0`, and
  `R > 0` for Kepler quantities; no new validation policy was introduced.
- The version-1 checkpoint stores unsigned `Omega` but not `R_dot`. Change 02
  must reconstruct signed `theta_dot` from `counter` and load old checkpoints
  without changing their field order. Persisting a meaningful future `R_dot`
  would require a later checkpoint version; it is not needed to preserve the
  present circular update.
- The isolated path needs an explicit source-radius branch in operator builders;
  using `R1()` there would move the displaced source to the origin.
- Direct public solver fields are used by tests and may be external API. Change
  02 should replace them coherently with one `OrbitState`, not retain writable
  aliases that create duplicate ownership.
