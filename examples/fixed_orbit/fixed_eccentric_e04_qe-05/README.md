# Fixed eccentric binary with p = 49.994

The example uses midpoint Hamiltonian Crank-Nicolson, a prescribed co-rotating
orbit with q = 1e-5 and e = 0.4, and a primary-centred 211 cloud.
The config separation 35.71 is the pericenter. The semimajor axis is 59.5167,
the apocenter is 83.3233, and the period is 2884.9321. The primary starts
on the positive x axis. The Kepler clock ramps up over two periods. The
orbital elements and masses stay fixed while both sources move.
The mask attenuates the outer field once per step. Horizon absorption is off.

## Start a fresh run

Run from the repository root in PowerShell with the project Python environment
active. This builds the executable, keeps output in a new local directory,
and generates initial data from this example's config.

```powershell
cmake -S . -B build-time-refactor -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build-time-refactor --target GravitationalMolecule -j 4
$exe = (Resolve-Path build-time-refactor/GravitationalMolecule.exe).Path
$cfg = (Resolve-Path examples/fixed_orbit/fixed_eccentric_e04_qe-05/config.toml).Path
$gen = (Resolve-Path initial_states/get_initial_state_COM_fast.py).Path
$work = New-Item -ItemType Directory -Force -Path ("runs/e04_q1e-5_p50/" + [guid]::NewGuid())
Push-Location $work.FullName
python $gen --config $cfg --output initial_modes.bin
if ($LASTEXITCODE -ne 0) { throw "Initial-state generation failed" }
& $exe $cfg initial_modes.bin
Pop-Location
```

## Dashboard and restart

From the repository root, point the dashboard at the campaign to pick a run.
Open http://127.0.0.1:8050 in the browser.

```powershell
python dashboard/live_dashboard.py runs/e04_q1e-5_p50 --port 8050
```

To resume, pass the exact Simulation_* folder printed by the solver to
`build-time-refactor/GravitationalMolecule.exe`. Preserve its config,
integrator marker, diagnostics and checkpoint. To extend a finished run,
increase its requested run length in its own config.toml. Keep the original
physics and integration method. Historical folders produced with the legacy
endpoint method cannot continue as midpoint runs.

To check the orbit, signed power, rate and dashboard endpoints on a completed
midpoint run with at least two steps, use
`python tests/python/check_e04_run.py "path/to/Simulation_..."`.

## Geometry-controlled midpoint run

Select `config_geometry.toml` in the same launch commands for a separate fresh run.
It retains the Kepler orbit, grid and atomic state, uses 50 as the maximum
timestep, and stops after 150 initial periods of physical time. The two-period
startup ends on a step boundary. The separate example has no outer absorber:
for production flux measurements check boundary effects or calibrate a CAP.
Compare against constant-timestep controls using the same boundary treatment.
Use `dt_step.txt` and `t.txt` to inspect the accepted physical-time mesh.
