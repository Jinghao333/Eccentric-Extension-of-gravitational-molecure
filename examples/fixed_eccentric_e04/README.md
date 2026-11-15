# Fixed eccentric `e = 0.4` workflow

Run these commands from the repository root in PowerShell:

```powershell
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build --parallel 4

python initial_states/get_initial_state_COM_fast.py `
  --config examples/fixed_eccentric_e04/config.toml `
  --output examples/fixed_eccentric_e04/initial_modes.bin

Set-Location examples/fixed_eccentric_e04
..\..\build\GravitationalMolecule.exe config.toml initial_modes.bin
```

The last command creates a new `Simulation_YYYY_MM_DD_HH_MM_SS` directory and
never overwrites an existing run. To inspect it in the dashboard, keep the run
terminal open or use a second PowerShell terminal:

```powershell
Set-Location "C:\physics\Thesis\Computational work\boson-cloud-binary-main"
python dashboard/live_dashboard.py examples/fixed_eccentric_e04 --port 8050
```

Then open `http://127.0.0.1:8050`.

For this configuration, `R=20` is the pericenter, the semimajor axis is
`a=20/(1-0.4)=33.3333`, the apocenter is `46.6667`, and the initial orientation
has the primary at `+R1` on the x axis. The startup ramp advances a slowed
Kepler clock; after startup the orbit follows the standard prescribed ellipse.
