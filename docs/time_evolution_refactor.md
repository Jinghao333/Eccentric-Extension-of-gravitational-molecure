# Time evolution refactor

Reference revision: `54f529f1f3f2b00ec57e5c06d570a92d99125a59`.
This is a structural refactor; it does not improve the order of time convergence.

## Ownership and evaluation

`CrankNicolsonModes` still owns the single `OrbitState`, field, operators,
simulation clock, next timestep and rolling power history. `OrbitMode` and
`TimeStepRule` are derived from existing settings, never stored independently.

The pre-field stage advances phase with the previous driving rate, then evaluates
endpoint orbital rates. The fixed Kepler trajectory and legacy circular-feedback
equations are unchanged. Feedback remains a separate post-field stage.

For a fixed grid and physics configuration, `(OrbitState, HamiltonianEvaluation)`
supplies the time-dependent operator inputs. The evaluation context is temporary:
it carries the evaluation time, companion weight and its derivative, and the
separate absorption weight/activation. Its factory evaluates existing schedules
without advancing the orbit. Builders no longer read the global simulation time
or orbital-clock ramp. Geometry settings and the legacy circular omission of
radial work remain unchanged.

`step_dt` is captured on entry and passed through CN matrix assembly, sink
corrections, preconditioning and solve. `select_next_timestep()` is the only
evolution-time timestep selector. It retains the strict post-startup gate and
`ramp_factor == 1.0`, selecting `P_circular / steps_per_period` only for legacy
circular inspiral. Constructor initialization, public user assignment and exact
checkpoint restoration of `dt` remain supported.

## Preserved ordering

1. Save the previous field and advance time using the captured step duration.
2. Advance orbital phase and evaluate endpoint orbital rates.
3. Build the endpoint Hamiltonian and absorbers.
4. Apply the existing CN solve and then the damping mask.
5. Evaluate the Schrodinger derivative on the post-mask endpoint field.
6. Apply legacy circular feedback using the midpoint field, endpoint absorbers
   and the **previous** rolling power. A crossing step uses only its portion
   after startup activation.
7. Select the next timestep from the updated orbit.
8. Compute observables with the completed step's unchanged Hamiltonian, then
   update rolling power once. Write the updated average to `p_bin.txt` and
   instantaneous orbital power to `p_bin_inst.txt`.
9. Write diagnostics in their existing order, then snapshot, flush and checkpoint
   at the existing cadence. Run-end and signal-stop checkpoints are unchanged.

For inspiral, saved orbital quantities are post-feedback while field diagnostics
use the pre-feedback endpoint operator. This legacy distinction is intentional.
Restart restores the saved next timestep, phase, field and power history without
selecting a new timestep or advancing the orbit. Checkpoint versions are unchanged.

## Verification

The new `time_workflow` fixture was first compiled against the original production
implementation. Seven cases cover fixed circular/eccentric motion, counter
rotation, optional potential switching, zero startup, exact/crossed startup
boundaries, legacy inspiral, displaced isolated BH, absorption, CAP and masks.
Each case takes 24 steps; uninterrupted results are compared with restarts at
steps 1, 3 and 8, including active feedback and timestep changes.

Full-precision binary traces record the field, orbital state, next timestep,
Schrodinger derivative, sparse entries/patterns of H, V1, V2, V, dVdt,
dVdt_switch, A and B, sink strengths and sink action. The complete sink matrix is
recorded at steps 0, 3, 8 and 24; a fixed complex probe is used at other steps.
All diagnostic files, snapshots and final checkpoints are also compared bytewise.
Reference captures are build artifacts, not scientific simulation data.

| Stage | Build and regression result | Original-reference comparison |
| --- | --- | --- |
| A | Original suite: 17 passed, initial-state Python test initially skipped for missing SymPy; subsequently passed against original executable with build-local dependencies | Reference captured; all 21 restart splits identical |
| B | Build and 8 focused tests passed | All seven cases byte-identical |
| C | Build and 10 focused tests passed | All seven cases byte-identical |
| D | Build and full suite: 18/18 passed, including both Python tests | All seven cases / 336 files (76,640,248 bytes) byte-identical; zero numerical differences in captured state/operators |

One Stage B capture attempt failed with Windows `Permission denied` during
checkpoint replacement. An isolated retry passed without changes to checkpoint
code. This filesystem failure was not treated as a numerical discrepancy.

The build uses GCC 16.2, Release (`-O3 -DNDEBUG`), Eigen/GSL from MSYS2 UCRT64.
The bytewise workflow capture uses one OpenMP/Eigen thread; the existing parallel
regression exercises one/four threads and both preconditioners.

The fixture also checks derived mode/rule getters and requires each completed
CN matrix to use the timestep captured before feedback, even when the next
timestep differs. No comparison tolerance was relaxed. Stage logs and reference
captures are retained under `build-time-refactor`.

```powershell
cmake -S . -B build-time-refactor -G Ninja -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON -DCMAKE_CXX_COMPILER=C:/msys64/ucrt64/bin/c++.exe -DCMAKE_MAKE_PROGRAM=C:/msys64/ucrt64/bin/ninja.exe -DCMAKE_PREFIX_PATH=C:/msys64/ucrt64 -DPython3_EXECUTABLE=C:/Users/HP/AppData/Local/Programs/Python/Python312/python.exe
cmake --build build-time-refactor -j 4
# Reference was recorded BEFORE production edits; do not overwrite it.
build-time-refactor/time_workflow.exe --compare build-time-refactor/baseline
$env:PYTHONPATH = (Resolve-Path build-time-refactor/python-deps).Path
ctest --test-dir build-time-refactor --output-on-failure -j 2
```

SymPy and mpmath were installed only under `build-time-refactor/python-deps` to
enable the existing initial-state tests. No global Python environment was changed.

## Limits

No midpoint Hamiltonian, adaptive timestep, trial propagation, arbitrary-time
orbital prediction or eccentric inspiral is implemented. Endpoint CN, startup
clock discretization, explicit circular feedback and the previous-power delay
retain their existing convergence properties. Bytewise equivalence is verified
on the recorded short fixtures with the same compiler/environment, not asserted
across compilers or for every possible configuration. No long baseline simulation
or existing scientific data was modified.
