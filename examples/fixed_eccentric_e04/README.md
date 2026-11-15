# Fixed eccentric `e = 0.4` workflow

This is the `e=0.4` counterpart of `examples/minimal/baseline.toml`: it keeps
the baseline grid, cloud, mass ratio, pericenter, runtime, mask and numerical
settings. It is therefore a large scientific run, not a short smoke test.

The existing baseline initial state is reusable because the eccentric run also
starts at `R=60`, at pericenter, with the primary on the positive x axis. Run
these commands from the repository root in PowerShell:

```powershell
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build --parallel 4

Set-Location examples/fixed_eccentric_e04
..\..\build\GravitationalMolecule.exe `
  config.toml `
  ..\minimal\initial_q1e-5_R60_baseline.bin
```

The last command creates a new `Simulation_YYYY_MM_DD_HH_MM_SS` directory and
never overwrites an existing run. To inspect it in the dashboard, keep the run
terminal open or use a second PowerShell terminal:

```powershell
Set-Location "C:\physics\Thesis\Computational work\boson-cloud-binary-main"
python dashboard/live_dashboard.py examples/fixed_eccentric_e04 --port 8050
```

Then open `http://127.0.0.1:8050`.

For this configuration, `R=60` is the pericenter, the semimajor axis is
`a=60/(1-0.4)=100`, the apocenter is `140`, and the Kepler period is about
`6283.15`. The initial orientation has the primary at `+R1` on the x axis. The
two-period startup ramp advances a smoothly slowed Kepler clock; after startup
the orbit follows the standard prescribed ellipse. `evolve_orbit=false` fixes
the orbital elements and masses, not the source positions.

To regenerate the initial state instead of reusing the baseline binary:

```powershell
Set-Location "C:\physics\Thesis\Computational work\boson-cloud-binary-main"
python initial_states/get_initial_state_COM_fast.py `
  --config examples/fixed_eccentric_e04/config.toml `
  --output examples/fixed_eccentric_e04/initial_modes.bin
```
