# Stationary primary 211 diagnostic

Use the production TOML and the original primary-centered MODESBIN file. The diagnostic copies both inputs into a fresh output directory. Its two controls use the same grid, retained modes, cloud-mass scaling and physical timestep. C0 places the primary at the origin and samples the analytic 211. D0 freezes the primary at its initial production position and loads the production data. Horizon absorption, boundary damping and orbital motion are disabled for both controls.

## Build and run

```sh
cmake -S . -B build-211 -DBUILD_CLOUD211_DIAGNOSTIC=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build-211 --target Cloud211Diagnostic
ctest --test-dir build-211 -R cloud211_diagnostic --output-on-failure
```

For the November 23 eccentric example, run from the repository root. Substitute the executable path printed by your generator. On Windows with Ninja it is under `build-211/diagnostics/cloud211/`.

```sh
build-211/diagnostics/cloud211/Cloud211Diagnostic runs/rerun_e04_q1e-5_60/20261123_131600/config.toml runs/rerun_e04_q1e-5_60/20261123_131600/initial_modes.bin --cases C0,D0 --periods 2 --radii 500,1000,10000 --output runs/211-pilot
python diagnostics/cloud211/analyze.py case runs/211-pilot/C0
python diagnostics/cloud211/analyze.py case runs/211-pilot/D0
python diagnostics/cloud211/analyze.py controls runs/211-pilot/C0 runs/211-pilot/D0 runs/211-pilot/comparison.json
```

Analysis requires Python 3.11+, NumPy and SciPy. The plots use Matplotlib. Pass `--no-plots` for a JSON report alone. The original period is the time unit for comparisons. Controls default to the configured physical timestep, or the startup spacing for a normalized clock. Supply `--dt` to refine it. Each radius selects the largest solver node at or below the requested value.

For a fresh checkout, use `diagnostics/cloud211/validation.toml`. Generate its input with `python initial_states/get_initial_state_COM_fast.py --config diagnostics/cloud211/validation.toml --output validation211.bin`, then run the diagnostic with that TOML and binary, `--periods 0.02 --radii 30,60` and a new output directory.

To repeat the production compatibility check, run `compatibility.py OLD_EXE BEFORE CONFIG...`, then `compatibility.py NEW_EXE AFTER --reference BEFORE CONFIG...` with the same config paths. Use fresh directories for BEFORE and AFTER. The check compares resolved settings and file hashes from fixed circular, normalized eccentric and circular-inspiral evolutions, including continuations from copies of the original checkpoints.

## Read the results

Radii and elapsed time are in primary-mass units. With the solver's Schrödinger convention, energy per particle is in inverse time units. Absolute fluxes retain the production cloud normalization; fractions in the report use its measured initial number (Q_0).

- `manifest.json` records the grid, case, signed displacement, reference period, effective timestep and mode order. `report.json` includes hashes of the exact input copies, full-domain initial measurements, bound-cluster energies, positive-energy estimates and threshold uncertainty. Frozen endpoints with appreciable weight invalidate the interior spectral interpretation. MODESBIN does not encode radii. For nonzero epsilon, verify the source generator used the solver grid.
- The centered channel admits a tridiagonal negative-spectrum projection. For large displaced runs, the analysis first resolves the negative 211 cluster, then measures the remaining state's positive-energy weight by Lanczos quadrature. A verified banded solve applies the sign-preserving spectral map `E/(E+s)` with positive definite `H+sI`. The report includes finite-order threshold brackets, the shift and its solve residual. Box and grid convergence are separate measurements.
- `flux.csv` contains every mode and a total row at every radius and time. `JQ` is outward number per unit time, `JE` is outward nonrelativistic energy per unit time, and `Jexcess = JE - E211*JQ` with signed continuum-reference `E211 = -mu^3/8` for the mass-one primary. `deltaQ_node` is signed trapezoidal transport from the node flux. The CN face current has its own radius and midpoint time.
- `balance.csv` reports full number, number through the extraction node, integrated CN face transport, and the residual of interior continuity. Compare its magnitude with the transport being interpreted. Divide number measurements by the measured initial number for fractions.
- `state_*.c128` are full-precision mode-major complex128 snapshots. `interior.coo` stores the actual conservative interior matrix as little-endian int32 row, int32 column, float64 real and float64 imaginary. `summary.png` shows the total and (1,1) node fluxes, continuity, radial disturbance and departure from a single CN phase. Times above 5% of the observed peak describe the sampled window. Tail and arrival claims require convergence and a later arrival at the outer sphere.
- The nested `Simulation_*` folder contains the solver's native intermediate files. The diagnostic's `manifest.json` describes the effective stationary run. Its nested copy of the source TOML describes the original binary.

## Convergence and follow-up

At the same source displacement, compare half the radial spacing at fixed box size, a higher angular cutoff, half the timestep, and a larger box at fixed spacing. Generate fresh translated inputs with the existing initializer for each changed spatial grid or cutoff. Use the same extraction radii. Compare signed areas, peak times, residuals and bound/positive-energy weights. The integrals measure net transport. A stationary mode current can cancel against other modes.

After both pilot reports exist, run `python diagnostics/cloud211/convergence.py PILOT EXECUTABLE VARIATION NEW_DIRECTORY`. Choose `radial`, `angular`, `timestep`, or `box`. Each invocation makes the matching control runs and writes differences relative to the pilot. Spatial variations assume the source used the repository's primary COM initializer and require its Python dependencies.

After reviewing D0, write a corrected input with `analyze.py case RUN/D0 --corrected NEW.bin`. This preserves the complex analytic polarization by projection into the numerical cluster and restores the input norm before cloud-mass scaling. Inspect the recorded cluster splitting and verify the state by repeating the driver with `--cases D1` and the new file. Compare D0 and D1 at identical grid, timestep, extraction radii and elapsed time.

For B0 and B1 run the existing `GravitationalMolecule config.toml INPUT.bin` from separate new directories, using the original and corrected inputs respectively. An archived B0 requires verified input and output provenance. Run `analyze.py binary B0_SIMULATION B1_SIMULATION NEW.csv` to compare matched histories. The reader converts the current writer's raw number flux by dividing by mu, and reverses its energy-flux sign. Repeat using the actual circular archive to investigate a circular signal.

To compare an existing B0 directly with the pilot, run `analyze.py archive RUN/D0 B0_SIMULATION NEW.csv`. It checks the config hash, integrator marker, mode order, extraction node and initial snapshot before interpolating the D0 histories to the recorded B0 times. Use `analyze.py plot RUN/D0 --mode 2 -2` to view any retained mode alongside the total. Each radius gets its own scale and replotting reads the stored measurements directly.
