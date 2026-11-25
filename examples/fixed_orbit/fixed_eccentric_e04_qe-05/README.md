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
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
if ($LASTEXITCODE -ne 0) { throw "Configure failed" }
cmake --build build --target GravitationalMolecule -j 4
if ($LASTEXITCODE -ne 0) { throw "Build failed" }
$exe = (Resolve-Path build/GravitationalMolecule.exe).Path
$cfg = (Resolve-Path examples/fixed_orbit/fixed_eccentric_e04_qe-05/config.toml).Path
& $exe --check-config $cfg
if ($LASTEXITCODE -ne 0) { throw "Solver cannot read this config" }
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
`build/GravitationalMolecule.exe`. Preserve its config,
integrator marker, diagnostics and checkpoint. To extend a finished run,
increase its requested run length in its own config.toml. Keep the original
physics and integration method. Historical folders produced with the legacy
endpoint method cannot continue as midpoint runs.

To check the orbit, signed power, rate and dashboard endpoints on a completed
midpoint run with at least two steps, use
`python tests/python/check_e04_run.py "path/to/Simulation_..."`.

## Geometry-controlled midpoint run

Select `config_geometry.toml` in the same launch commands for a separate fresh run.
Read the mass ratio, timestep rule and boundary directly from the chosen TOML
with `--check-config`. Regenerate the translated initial state with that same file.
Use `dt_step.txt` and `t.txt` to inspect the accepted physical-time mesh.

## Geometry run with the solve density of fixed dt=10

`config_normalized_dt10.toml` is the production-sized normalized-geometry
example. For this orbit, fixed `dt=10` gives 288.493 solves per initial
period. The normalized grid requires an integer count, so the closest choice
is 288 steps per orbit. Its mean spacing is 10.01712536, only 0.171% larger,
while the accepted nodes are concentrated near rapidly varying geometry.

Run it directly, without a Python launch wrapper. From the repository root:

```powershell
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
if ($LASTEXITCODE -ne 0) { throw "Configure failed" }
cmake --build build --target GravitationalMolecule -j 8
if ($LASTEXITCODE -ne 0) { throw "Build failed" }

$repo = (Resolve-Path .).Path
$work = Join-Path $repo "runs/orbital-geometry-dt10"
New-Item -ItemType Directory -Force -Path $work | Out-Null
$cfg = (Resolve-Path "examples/fixed_orbit/fixed_eccentric_e04_qe-05/config_normalized_dt10.toml").Path
$exe = "$repo/build/GravitationalMolecule.exe"
& $exe --check-config $cfg
if ($LASTEXITCODE -ne 0) { throw "Solver cannot read this config" }

$env:OMP_NUM_THREADS = "8"
$env:OMP_DYNAMIC = "FALSE"
$env:OPENBLAS_NUM_THREADS = "1"

Push-Location $work
python "$repo/initial_states/get_initial_state_COM_fast.py" --config $cfg --output initial_modes.bin
if ($LASTEXITCODE -ne 0) { throw "Initial-state generation failed" }
& $exe $cfg initial_modes.bin
Pop-Location
```

The `Simulation_*` directory is created under `$work`. To resume, keep its
checkpoint, diagnostics, `config.toml`, and `integrator.txt` together, then
pass the simulation directory itself to the same executable:

```powershell
& "$repo/build/GravitationalMolecule.exe" `
  "$work/Simulation_YYYY_MM_DD_HH_MM_SS"
```

Monitor the campaign from another terminal with:

```powershell
python dashboard/live_dashboard.py "$work" --port 8010
```

Open `http://127.0.0.1:8010/` (without a trailing comma). The historical
fixed-`dt=10` calculation used legacy endpoint CN and a per-step mask, whereas
this variable-step run necessarily uses midpoint CN and a CAP. It is therefore
a useful production reference but not a one-variable timestep comparison. A
strict controller comparison must use the same integrator and CAP in both runs.

After the solver has stopped writing, verify the stored orbit, power/rate
normalizations, and dashboard API with:

```powershell
$run = Get-ChildItem $work -Directory -Filter "Simulation_*" |
       Sort-Object LastWriteTime -Descending | Select-Object -First 1
python tests/python/check_e04_run.py $run.FullName
```

Do not use this completed-run checker concurrently with the solver: it uses
plain NumPy text loading and can encounter the final row while that row is
still being appended. Use the dashboard for live monitoring. The accepted
physical timestep is in `dt_step.txt`; `t.txt`, `R_orb.txt`, `R_dot_orb.txt`,
`theta_orb.txt`, `phase_orb.txt`, and `Omega_orb.txt` reconstruct the orbit.
The main power files are `p_bin_inst.txt`, `p_bin.txt`, `energy_flux.txt`, and
`flux.txt`. `checkpoint.bin` is the authoritative restart state.

`config_endpoint.toml` remains available for legacy endpoint-CN compatibility
runs. It is a separate fixed-timestep method and must not be used to resume a
midpoint run.
