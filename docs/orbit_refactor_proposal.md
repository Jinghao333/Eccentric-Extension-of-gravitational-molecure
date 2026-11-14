# OrbitState refactor proposal

## 1. Decision status and non-goals

This is a proposal for Stage 2: introduce an `OrbitState` abstraction and migrate the existing circular implementation to it without changing observable behavior. It is not an implementation and is not an ADR; unresolved ownership and API questions remain for review.

Evidence labels:

- **Verified behavior** is directly established by current source/tests and must be preserved in Stage 2.
- **Physical assumption** describes the circular model encoded by the source.
- **Suspected inconsistency** is documented but must not be silently corrected during migration.
- **Proposed design** is new architecture for review.

Stage 2 non-goals:

- No eccentric elements or Kepler-equation solver.
- No analytic replacement for the discrete `theta`/`phase` update.
- No midpoint Hamiltonian.
- No ramp formula or gate change.
- No `dVdt` correction.
- No diagnostic time-alignment correction.
- No change to GW, cloud-feedback, or accretion formulas.
- No production change to CAP, mask, absorption, CN, or preconditioning algorithms.

## 2. Behavior contract to preserve

The migration is successful only if all of the following remain true.

### 2.1 Circular geometry and source convention

1. Initial masses are `M1=1`, `M2=q`.
2. `R1=M2/(M1+M2)*R`, `R2=M1/(M1+M2)*R`.
3. Primary position is `R1*(cos(theta),sin(theta),0)`.
4. Secondary position is at `theta+pi` with radius `R2`.
5. The potential phase remains `exp(i*(M_B-M_A)*theta)`.
6. `counter` changes only the signed coordinate motion and rotational derivative sign; `Omega` and `phase` remain nonnegative/increasing.

### 2.2 Exact discrete phase/ramp sequence

For every step, in this order:

```text
step_start = current time
completed_dt = dt
endpoint_time = step_start + completed_dt

theta += sign(counter) * current_Omega * completed_dt
phase += current_Omega * completed_dt

evaluate startup ramp at endpoint_time
set endpoint Omega from ramp or current circular Kepler state
build endpoint potential, dVdt, absorption, and H
build and solve endpoint-H CN system
apply optional mask
evaluate endpoint generator on post-mask field
apply optional post-solve circular feedback
calculate and write diagnostics
calculate new rolling p_bin only after feedback
```

The old `Omega` advances `theta` and `phase`. The endpoint ramp supplies the `Omega` for the next step. This is verified by `src/CrankNicolsonModes.cpp:145-168`, `tests/orbital_ramp.cpp`, and `tests/parallel_correctness.cpp:139-177`.

### 2.3 Ramp policy

The immutable schedule anchor remains

\[
t_{ramp}=\texttt{ramp_periods}\,P_0.
\]

Current `P` must never replace `P0` in that gate.

- Frequency: smoothstep from 0 to `Omega_0` during startup; Kepler value afterward.
- Companion potential: either the same smoothstep or full strength from the start.
- Absorption: zero through startup, then smoothstep over one additional `P0`.
- Backreaction: unavailable through startup, then smoothstep over one additional `P0`.
- A crossing step uses only post-gate `orbit_dt`, multiplied by the endpoint backreaction weight.

### 2.4 Fixed versus evolving circular modes

- `evolve_orbit=false` prevents separation, mass, `q`, COM-radius, `P`, and `dt` feedback updates.
- It does not disable frequency/potential startup ramps, source rotation, field absorption, CAP, mask, power diagnostics, or rolling averaging.
- With fixed-orbit absorption, field norm can decrease while black-hole masses remain fixed.
- `evolve_orbit=true` applies feedback after the field solve, using midpoint field accretion rates and the prior rolling `p_bin`.

### 2.5 Persistence

- Checkpoint continuation must remain bitwise identical at fixed thread count.
- Diagnostic-only restart must retain its current precision-limited behavior.
- Existing `GMCKPT01` files should remain readable if the binary field order is not changed.
- The rolling average's full internal state remains part of the exact checkpoint.

## 3. Proposed domain split

The current class mixes four distinct concepts:

1. Instantaneous binary state.
2. Immutable initial circular reference values.
3. Startup/activation policy.
4. Post-solve circular feedback policy.

The refactor should separate them without changing their order of application.

### 3.1 `OrbitState`: one instantaneous circular snapshot

**Proposed design:**

```cpp
struct OrbitState
{
    double time;

    double separation;       // current R
    double primary_radius;   // current R1
    double secondary_radius; // current R2

    double primary_mass;     // current M1
    double secondary_mass;   // current M2
    double mass_ratio;       // current q = M2/M1

    double angle;            // signed theta
    double phase;            // monotonically increasing averaging phase
    double angular_frequency;// nonnegative current Omega
    double period;           // current P
};
```

Recommended invariants for ordinary circular binary mode:

\[
R_1+R_2=R,
\qquad
q=M_2/M_1,
\qquad
R_1=M_2R/(M_1+M_2),
\qquad
R_2=M_1R/(M_1+M_2),
\qquad
P=2\pi/\Omega
\]

after startup and after feedback. During the startup frequency ramp, `period` intentionally remains the initial Kepler period even though `Omega` is ramped; therefore `P=2*pi/Omega` is not a universal invariant during startup.

The fields above include values currently serialized separately. Retaining them avoids changing floating-point recomputation and checkpoint semantics during the migration. A later design may choose to derive redundant fields, but Stage 2 should not.

`dt` should not be in `OrbitState`: it is field-integrator control and is changed for the next step only. `counter` should not be in the snapshot: it is an immutable motion convention. Ramp weights, `p_bin`, and `dVdt` are also not orbit geometry.

### 3.2 `CircularOrbitReference`: immutable construction data

**Proposed design:**

```cpp
struct CircularOrbitReference
{
    double initial_frequency; // current Omega_0
    double initial_period;    // current initial_period
    bool counter_rotating;    // current counter
};
```

This makes it explicit that startup schedules use initial, not current, orbital quantities. It also prevents an evolving `P` from accidentally moving the ramp gate.

The configured initial separation and mass ratio need not remain here after constructing `OrbitState`, provided config copying/restart still retains them. `atom.a` remains primary spin metadata and must not be renamed to an orbital quantity.

### 3.3 `StartupRampPolicy`: schedule, not orbit dynamics

**Proposed design:**

```cpp
struct StartupRampState
{
    double frequency_factor;            // current ramp_factor during startup
    double companion_potential_factor;  // current potential_ramp_factor
    double absorption_factor;           // separate post-startup ramp
    double backreaction_factor;         // separate post-startup ramp
    bool absorption_operator_active;    // current time >= ramp end
};

class StartupRampPolicy
{
public:
    StartupRampState evaluate(double endpoint_time,
                              const CircularOrbitReference&) const;
    double ramp_end_time(const CircularOrbitReference&) const;
    double post_gate_duration(double step_start,
                              double endpoint_time,
                              const CircularOrbitReference&) const;
};
```

The policy owns `ramp_periods` and `ramp_secondary_potential`. It must reproduce the existing formulas exactly. In particular:

- It should return weights, not integrate them.
- It should not return derivatives such as `potential_factor_dot`; the present `dVdt` omits that term.
- `post_gate_duration()` must reproduce `t_i-max(step_start,ramp_end)`.
- Equality and branch behavior at `ramp_end_time()` must match current `<` and `>=` conditions.

`ramp_factor` and `potential_ramp_factor` can remain cached members during the first migration, but their authoritative source should become the returned `StartupRampState`, not independent mutation scattered through the solver.

### 3.4 Circular motion and feedback policies

Separate coordinate advancement from post-solve feedback.

**Proposed design:**

```cpp
class CircularKinematics
{
public:
    void advance_angle_and_phase(OrbitState&, double completed_dt,
                                 const CircularOrbitReference&) const;

    void set_endpoint_frequency(OrbitState&,
                                const CircularOrbitReference&,
                                const StartupRampState&) const;

    void refresh_com_geometry(OrbitState&) const;
};
```

`advance_angle_and_phase()` must use the state's entry `angular_frequency`. `set_endpoint_frequency()` then applies the endpoint startup factor or the circular Kepler formula.

Post-solve evolution can be represented by two policies sharing one interface:

```cpp
struct OrbitFeedbackInput
{
    double orbit_dt;
    double backreaction_factor;
    double phase_averaged_cloud_power; // prior p_bin
    double primary_accretion_rate;
    double secondary_accretion_rate;
};

class OrbitEvolutionPolicy
{
public:
    virtual void apply(OrbitState&, const OrbitFeedbackInput&,
                       double steps_per_period, double& next_dt) const = 0;
};

class FixedCircularEvolution : public OrbitEvolutionPolicy;
class CircularInspiralEvolution : public OrbitEvolutionPolicy;
```

`FixedCircularEvolution::apply()` is a no-op. `CircularInspiralEvolution` contains the current forward-Euler GW/cloud/accretion formulas.

A virtual interface is not required if a tagged variant or two explicit functions is preferred. The important design boundary is that `evolve_orbit` selects only the post-solve feedback policy; it must not select or bypass startup kinematics.

### 3.5 Feedback bookkeeping outside `OrbitState`

Keep these in solver/control state:

- `double dt`: integrator control.
- `double p_bin`: lagged rolling feedback input.
- `RollingPhaseAverage power_average`: averaging history.
- `Eigen::VectorXcd psi_last`: field midpoint support.
- `SparseComplex dVdt`: derived field operator.
- `StartupRampState`: endpoint policy snapshot for the current field step.

This separation makes it harder to confuse physical orbit geometry with numerical or diagnostic state.

## 4. Proposed Stage 2 step orchestration

The refactored `step()` should make the existing sequence explicit without changing it:

```text
old_field = psi
step_start = orbit.time
completed_dt = dt

orbit.time += completed_dt

# Must use entry Omega, before endpoint ramp evaluation.
circular_kinematics.advance_angle_and_phase(
    orbit, completed_dt, circular_reference)

ramp = startup_policy.evaluate(orbit.time, circular_reference)
circular_kinematics.set_endpoint_frequency(
    orbit, circular_reference, ramp)

build_orbit_dependent_perturbation(orbit, ramp)
update_CN_values(completed_dt)
solve existing CN system
apply optional mask
compute endpoint dpsi_dt from post-mask field

if circular inspiral policy is selected and orbit.time > ramp_end:
    orbit_dt = startup_policy.post_gate_duration(
        step_start, orbit.time, circular_reference)
    compute accretion rates from midpoint field and endpoint sink
    feedback_input.phase_averaged_cloud_power = p_bin  # previous value
    circular_inspiral.apply(orbit, feedback_input,
                            steps_per_period, dt)

compute field diagnostics using already-built endpoint H/dVdt
compute current p_bin_inst
p_bin = rolling_average.update(orbit.time, orbit.phase, p_bin_inst)
write post-feedback orbit scalars and field diagnostics
checkpoint
```

Two states are intentionally implicit in current behavior:

- The **field-step orbit** used to build `H`, sinks, and `dVdt` before feedback.
- The **post-feedback orbit** stored in the member and written to output/checkpoint.

The migration should either preserve the field-step `OrbitState` as a local copy or guarantee that operators are fully built before mutating the main state. It must not rebuild `H` after feedback inside the completed step.

## 5. Signature and ownership changes

### 5.1 Functions that should change

| Current function/location | Proposed ownership/signature direction | Reason |
|---|---|---|
| `CrankNicolsonModes` constructor, `src/CrankNicolsonModes.cpp:11-50` | Construct `OrbitState`, `CircularOrbitReference`, ramp policy, and evolution policy from existing arguments/config. Later accept a cohesive circular-orbit config object. | Remove scattered scalar initialization while preserving the public constructor during the first increment if desired. |
| `step(int)`, lines 145-234 | Orchestrate mutations through `OrbitState` and policy objects; preserve exact order. | Central behavior boundary. |
| `build_time_dependent_perturbation(double)`, `src/CrankNicolsonModes_hamiltonian.cpp:245-268` | `build_orbit_dependent_perturbation(const OrbitState&, const StartupRampState&)`. Stop deriving schedule and geometry inside the matrix builder. | Separate orbit/ramp policy from operator construction. |
| `build_static_perturbation()`, lines 176-222 | Accept `const OrbitState&` and companion-potential factor explicitly. | Make all orbit inputs visible. Preserve current `dVdt` formula. |
| `F1()` / `F2()`, lines 12-42 | Accept `const OrbitState&`, or accept explicit source radius and mass fraction. | Remove hidden reads of mutable members. |
| `orbital_phase(int)`, lines 106-109 | `orbital_phase(int, double theta)` or `(..., const OrbitState&)`. | Explicit angle dependency; useful to potential and sink. |
| `build_spherical_absorption()`, `src/CrankNicolsonModes_absorption.cpp:153-182` | Accept orbit snapshot and absorption factor/operator-active policy. | Separate absorption schedule from geometry and operator assembly. |
| `build_swave_projector()`, lines 127-150 | Accept explicit `theta` or orbit snapshot. | Remove hidden orientation read. |
| `update_time_dependent_increment(double)`, `src/CrankNicolsonModes_evolution.cpp:17-72` | Split into accretion-rate calculation and `CircularInspiralEvolution::apply(OrbitState&, feedback, ..., next_dt)`. | Separate field-derived input from orbit mutation. |
| `assemble_hamiltonian(bool)`, `src/CrankNicolsonModes_hamiltonian.cpp:224-243` | Optionally accept `absorption_operator_active`; otherwise unchanged internally. | Activation comes from ramp policy. |
| `run()` / `evolve()`, `src/CrankNicolsonModes.cpp:236-289` | Read time/separation from `OrbitState`; retain run order and stop rule. | State ownership migration. |
| `save_checkpoint()` / `load_checkpoint()`, `src/CrankNicolsonModes_restart.cpp:139-233` | Add field-wise `write_orbit_state()` / `read_orbit_state()` in the exact existing scalar order. Do not raw-write the struct. | Encapsulation without padding/version/byte changes. |
| `reconstruct_state_from_diagnostics()`, lines 238-286 | Return/fill `OrbitState`; retain left-endpoint phase replay. | Consolidate reconstruction. |
| `resume()`, lines 288-399 | Rebuild operators from restored `OrbitState`; preserve saved/restored `Omega` behavior until builder no longer mutates it. | Prevent restart drift. |
| diagnostic writes in `step()` | Read orbit columns from `OrbitState`. | Centralized output mapping. |
| `configure_isolated_black_hole()`, `src/CrankNicolsonModes.cpp:58-74` | Populate a valid stationary `OrbitState` and stationary policy/reference. | The same solver class and persistence code serve isolated mode. |

### 5.2 Functions that should remain numerically unchanged

| Function/subsystem | Why unchanged |
|---|---|
| `load_config()` parsing behavior | Stage 2 can rearrange how values are passed, but existing circular keys/defaults remain valid. |
| `build_grid()` | No orbit dependency. |
| `build_modes()` | No dynamic orbit dependency; equatorial parity remains. |
| `set_initial_state_full()` | Cloud normalization is independent of orbit-state ownership. |
| `load_initial_modes_binary()` | File mapping is independent of orbit evolution. |
| `build_stationary_hamiltonian()` | Kinetic/centrifugal operator is orbit-independent. |
| `build_coupling_cache()` | Wigner/Gaunt coefficients and selection rules are orbit-independent. |
| `build_sparse_element_caches()` | Circular state changes values, not sparsity. |
| `update_absorption_profile()` | Already takes explicit radius and width; its numerical quadrature/cache rule remains. |
| `apply_sink()` | Consumes already-built sink vectors/strengths. |
| `schrodinger_rhs()` | Consumes already-built `H` and sinks. |
| `build_CN_matrices()` / `update_CN_values()` | Same `A=I+i dt H/2`, `B=I-i dt H/2`. A signature may accept `dt` explicitly without changing arithmetic. |
| `solve_CN()` and preconditioners | No orbit ownership; preserve tolerances, fallback, and Woodbury handling. |
| `restricted_expectation()` | Operator/state diagnostic helper. |
| `get_flux()` / `get_energy_flux()` | No direct orbit input. |
| `get_binary_power()` | Continue taking the already-built current `dVdt`. |
| `get_dVexpectation_dt()` | Preserve current diagnostic meaning. |
| `get_norms()` / `get_energy()` / `get_energy_per_mode()` | Preserve current cutoff conventions during Stage 2. |
| `build_damping_mask()` / `build_cap()` | Independent outer-boundary algorithms. |
| `save_snapshot()`, `save_run_info()`, quantity writers | File formats remain unchanged. |
| `RollingPhaseAverage` | Already accepts explicit time/phase/value and serializes itself. |

### 5.3 Temporary compatibility layer

Current tests and potential library users access public scalars such as `cn.R`, `cn.theta`, and `cn.Omega` directly. There are three migration options:

1. **Breaking cleanup:** replace them with `cn.orbit().separation`, etc., and update tests/users.
2. **Deprecated accessors:** provide `R()`, `theta()`, and related references returning fields of the authoritative `OrbitState`.
3. **Two-phase internal migration:** first add explicit capture/apply helpers while old fields remain authoritative; then replace storage after all functions consume `OrbitState`.

The proposal recommends option 3 followed quickly by option 2 or 1. Do not keep both public scalars and a nested `OrbitState` as independently writable authoritative copies; dual state would create harder-to-detect desynchronization than the current design.

A safe transitional sequence is:

```text
Phase A:
    define OrbitState
    add make_orbit_state_from_members() and apply_orbit_state_to_members()
    pass snapshots explicitly into pure consumers
    preserve old member storage and all tests

Phase B:
    make one nested OrbitState authoritative
    replace member reads/writes
    provide narrow compatibility accessors if required
    update tests to use the authoritative state API

Phase C:
    remove transitional copy/apply helpers
```

Each phase must avoid long-lived dual mutation.

## 6. Incremental migration plan

### Step 0: characterization before structural change

Add small, non-baseline tests for behavior not directly asserted today:

- First startup step advances with entry `Omega=0` before evaluating endpoint ramp.
- `phase` increases while counter `theta` decreases.
- Field operators remain pre-feedback while written orbit scalars are post-feedback.
- Feedback consumes the previous `p_bin`.
- Binary checkpoint scalar order and magic remain stable.

Do not run or create a new long baseline.

### Step 1: introduce data types only

Add `OrbitState`, `CircularOrbitReference`, and `StartupRampState` with field-wise equality/test helpers. No production path uses them yet. Unit-test mapping to/from current members.

### Step 2: centralize derived circular geometry

Add one helper that reproduces exactly:

```text
q = M2/M1
R1 = M2/(M1+M2)*R
R2 = M1/(M1+M2)*R
Omega = sqrt((M1+M2)/R^3)
P = 2*pi/Omega
```

Use it only where the current code performs the same complete block: constructor initialization as applicable, feedback update, and diagnostic restart. Do not call it during startup in a way that overwrites ramped `Omega`.

### Step 3: extract startup policy

Move smoothstep/gate computation into `StartupRampPolicy`. Keep the call site at the same endpoint position in `step()`. Compare every ramp test before and after.

Do not calculate an analytic integrated phase and do not add a potential-ramp derivative.

### Step 4: pass orbit snapshots into operator builders

Migrate `F1`, `F2`, `orbital_phase`, potential assembly, and absorber projector construction to explicit orbit inputs. Wigner caches, radial formulas, sparse offsets, and matrix arithmetic remain identical.

`tests/parallel_correctness.cpp` is the primary guard for this step.

### Step 5: extract fixed/inspiral feedback policy

Move the body of `update_time_dependent_increment()` behind the policy boundary while preserving:

- midpoint field state;
- endpoint sink geometry/strength;
- endpoint backreaction weight;
- crossing-step `orbit_dt`;
- prior `p_bin` input;
- omitted direct mass-driven `R` drift;
- forward Euler;
- next-step-only `dt` update.

Use a no-op fixed policy rather than branching inside the mathematical update.

### Step 6: make nested `OrbitState` authoritative

Replace scalar storage only after consumers and mutators use explicit state. Update public API/tests or add compatibility accessors. Run all existing unit, regression, physics, integration, and Python initial-state tests.

### Step 7: migrate persistence

Serialize `OrbitState` fields individually in the existing order:

```text
time, angle, phase, separation,
primary_mass, secondary_mass, mass_ratio,
primary_radius, secondary_radius,
angular_frequency, period
```

then existing `dt` and `p_bin`. Keep `GMCKPT01` only if the byte layout truly remains identical. Otherwise bump the version and retain a reader for version 1.

Refactor diagnostic reconstruction to return the same state values and replay phase/power identically.

### Step 8: remove transitional duplication and document the API

Delete snapshot/copy shims, retain one source of truth, and update Doxygen/comments to distinguish:

- coordinate angle from averaging phase;
- instantaneous state from initial circular reference;
- startup policy from feedback policy;
- integrator `dt` from orbital period.

Only after this stage passes review should eccentric fields/providers be designed.

## 7. Test plan for the migration

No new full baseline simulation is needed. Use the existing targeted suite:

1. `config_and_io`: construction mapping, grid/modes, normalization.
2. `orbital_ramp`: highest-priority compatibility suite for endpoint ramp order and fixed gates.
3. `parallel_correctness`: exact matrix values and endpoint CN solve.
4. `fixed_orbit_absorption`, `spherical_absorption`, `swave_projector`: geometry and activation.
5. `flux_time_alignment`: endpoint/post-mask diagnostic state.
6. `physics`: fixed-orbit Jacobi convergence, mirror symmetry, inspiral, mass balance, outer absorbers.
7. `rolling_average`: independent `phase` semantics and serialization.
8. `restart`: checkpoint byte identity and diagnostic reconstruction.
9. `command_line`: end-to-end file shapes, run duration, resume/extension.
10. Python initial-state tests: unchanged circular placement convention.

Suggested new focused tests, all short:

- `orbit_state_mapping`: old members round-trip through `OrbitState` exactly.
- `first_ramp_step_order`: records old-`Omega` coordinate advance and endpoint new `Omega`.
- `feedback_snapshot_order`: independently compares pre-feedback `H` with post-feedback reported state.
- `lagged_cloud_power`: injects distinguishable successive powers and verifies which drives `R`.
- `checkpoint_v1_layout`: loads a small committed version-1 fixture or verifies explicit field ordering.

## 8. Suspected inconsistencies deliberately left unchanged

The interface migration must carry these behaviors forward until separately reviewed:

1. `dVdt` excludes `potential_ramp_factor`'s time derivative.
2. Field/Hamiltonian diagnostics and orbit diagnostics are on opposite sides of post-solve feedback.
3. Cloud feedback consumes the preceding rolling `p_bin`.
4. Crossing-step absorption and mass accretion may cover different effective durations.
5. Norm/per-mode energy and total energy/power use different extraction-node inclusion.
6. Extraction indexing ignores `epsilon`.
7. C++ and Python grids differ for nonzero `epsilon`.
8. `H` may contain a CAP despite being documented as conservative.
9. The CAP profile denominator may not match the apparent intended layer normalization.
10. Constructor logging observes the default `evolve_orbit=true` before `main()` applies config.
11. Period-based run length counts the initial sample, producing one fewer evolved interval than the raw `floor(num_periods*P/dt)` value suggests.

These should become individual, test-backed change proposals after the `OrbitState` migration, not incidental edits within it.

## 9. Unresolved design decisions and risks

1. **Time ownership:** Should `OrbitState.time` become authoritative, or should a separate simulation clock own `t_i` and time-tag orbit snapshots? Including time improves snapshot coherence; separating it keeps orbit data purely geometric.
2. **Public API compatibility:** Are direct public fields (`R`, `theta`, etc.) considered supported library API? This determines whether accessors/deprecation are required.
3. **Policy representation:** Use virtual `OrbitEvolutionPolicy`, `std::variant`, or explicit fixed/inspiral branches. The behavior boundary matters more than runtime polymorphism.
4. **Checkpoint compatibility:** Preserve `GMCKPT01` byte-for-byte or introduce a versioned `OrbitState` block with a version-1 reader.
5. **Redundant derived fields:** Keep `q`, `R1`, `R2`, and `P` stored for exact compatibility, or derive them from minimal state later. Stage 2 should keep them stored.
6. **Field-step versus post-feedback state:** Should a named local `field_orbit` snapshot make current mixed-time behavior explicit? This is safer than relying on operator matrices to retain the pre-feedback state implicitly.
7. **Location of `p_bin`:** It is feedback controller state, not geometry. Decide whether it belongs to the solver, a feedback-state object, or the inspiral policy.
8. **Isolated-black-hole support:** The same state and policy API must represent `q=0`, `R2=0`, `Omega=0`, no phase averaging update, and immediate isolated absorption.
9. **Floating-point identity:** Reordering algebra while deriving COM radii, masses, or frequency can break bitwise restart tests even when mathematically equivalent.
10. **Future eccentric extensibility:** Do not add eccentric-only fields now. The circular migration should establish a stable snapshot/policy boundary without predicting the final eccentric anomaly representation.
