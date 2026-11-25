# Xu Appendix A — Scientific Analysis Dashboard

**Scope:** outgoing scalar ionization, fixed prescribed eccentric Newtonian orbit, small-`q` Appendix-A response, saved-v4 schema. The dashboard does **not** rerun the production solver. Its endpoint is one-point physics interpretation, **not** a new physics stage.

## Launch

Install from the standalone project root:

```bash
python -m pip install -r requirements.txt
python -m analysis_dashboard.server --runs runs --port 8765
```

Open **http://127.0.0.1:8765** in a local web browser. `--runs` may be a parent directory with one folder per run, or one particular completed run folder. Each completed run must contain **only** the authoritative `config.json` and `physics.h5`. The program **reads** these files, never mutates them, and creates no CSV/HDF5 summaries automatically.

This package is **self-contained for both computation and analysis**: the full driver and its nine local physics modules live in the same project. The nine local modules are independently reconstructed implementations; cross-checking against the original remote source is an additional scientific validation task.

## Architecture

```
config.json + physics.h5              # only canonical saved simulation inputs/results
        |
        v
primitive_store.py                    # indexed unique Fourier and outgoing-channel primitives
        |
        v
derive_quantities.py                 # verified outgoing-flux reconstruction
        |
        v
analysis_dashboard/core.py            # analytic background, summaries, channel sources/overlaps
        |
        v
analysis_dashboard/server.py          # read-only localhost JSON API
        |
        v
static/{index.html,app.js,style.css}  # responsive interactive Plotly dashboard
```

`Plotly.js` is served from the installed Python Plotly package at `/plotly.min.js`, **without a remote CDN**. It supports interactive zoom, autoscaling, pan, tooltips, legend toggles and vector-figure export. The existing live-dashboard project provided design and interaction inspiration, not the time-domain physics definitions.

## Dashboard pages

1. **Physical system** — automatic physical/orbital parameter report, `T=2π/Ω`, unperturbed equatorial cloud mass-density heatmap with moving *companion marker only*, fixed Kepler ellipse, peak radii, radial mass shell distribution and enclosed mass fraction; total orbital ionization power, Fig. 8 normalization, leading-order fractional mass-loss rate.
2. **Channel atlas** — faithful Fig. 11 angular grouping `sum_g P_(ell_f,m_f,g)` as four marker families `m_f=ell_f, ell_f−2, ell_f−4, ell_f−6`; top-10 **individual** `(ell_f,m_f,g)` ranks and cumulative fraction; angular heatmap; eccentric harmonic spectrum, optionally filtered by the chosen angular mode. Click a nonzero channel to open the third page.
3. **Channel physics** — for selected nonzero open channel: the saved `u_(l,m,g)(r_j)`, reconstructed partial `S_(l,c)(r_j)`, saved coherent `S_c(r_j)`, analytic regular Whittaker `Z_-(r_j)`, `Z_-S`, grid-trapezoid *diagnostic* cumulative overlap, saved **adaptive** partial overlaps `I_(l,c)`, interference cross terms, analytic Wronskian and outgoing amplitude `A_c=Z_c^infinity`. All expensive orbit Fourier/production Green integrals remain untouched.

## Essential physics equations

The **production** field radial function obeys

\[
\int_0^\infty r^2 |R_i(r)|^2\,dr=\frac{M_b}{2\mu^2}.
\]

The unperturbed *mass* density plotted in the `xy` plane is therefore

\[
\rho_M(\mathbf r)=2\mu^2|R_i(r)Y_{\ell_i m_i}(\theta,\phi)|^2,
\qquad \int\rho_M\,d^3x=M_b.
\]

Do **not** call `|R_i Y|²` the cloud's physical mass density without the `2μ²` factor. This distinction corrects a normalization error in the previous draft dashboard contract. It is independently checked numerically against the frozen production radial profile in tests.

For each `c=(ell_f,m_f,g)`, with `m=m_f-m_i`,

\[
\begin{aligned}
S_{l,c}(r)&=2\mu^2rR_i(r)C(l,m;\ell_i,m_i;\ell_f,m_f)u_{lmg}^{(t)}(r),\\
I_{l,c}&=\int Z_{-,\ell_fg}(r)S_{l,c}(r)\,dr,\qquad I_c=\sum_l I_{l,c},\\
A_c&=e^{\pi\eta_g/2}I_c/W_{\ell_fg},\\
\dot E_c^{s,\infty}&=2g\Omega_{\rm orb}k_{{\rm KG},g}|A_c|^2.
\end{aligned}
\]

The dashboard computes the total power from stored **complex** outgoing amplitudes and configuration. The normalized Fig. 8 ordinate is `P_total / [q²(M_b/M)]`. The outgoing positive current in the Stage-VI convention is `J_c=2*k_KG,g*|A_c|²`. A **leading-order fractional** cloud depletion rate is shown as `Γ_ion = μ*sum_c(J_c)/M_b`, assuming the cloud rest mass is `μ Q_b`; this normalization remains a model-convention check before claiming agreement with a published ionization rate.

`S_total` is checked against a coherent sum of reconstructed `S_l`, and `A_c` is independently checked against `exp(πη/2)*sum_l(I_l)/W` by the data reader. Failure is explicit: the dashboard will not silently pretend inconsistent data are valid.

For multipoles `(l,l')`, the interference contribution is

\[
2\operatorname{Re}(I_{l,c}I_{l',c}^*),\quad
\mathcal C_c=\frac{|\sum_l I_{l,c}|}{\sum_l|I_{l,c}|}.
\]

The `I_c(r_j)` curve is **not** the production adaptive overlap: it is a cumulative **trapezoidal** diagnostic sampled on the saved grid, with the actual saved production `I_c` shown as a reference line. No convergence sweep is performed.

## Python API (same functions used by dashboard)

```python
from analysis_dashboard.core import (
    load_run, run_overview, channel_atlas, channel_detail,
    orbital_parameters, radial_field, density_slice,
)

run = load_run('runs/<my_completed_run>')
overview = run_overview(run)
atlas = channel_atlas(run)

# Inspect the most important *individual* outgoing channel.
lead = atlas['top'][0]
selected = channel_detail(run, lead['ell'], lead['m'], lead['g'])
print(selected['response'])
```

This is suitable for importing into a subsequent physics-learning notebook. Such a notebook should call these functions rather than maintaining duplicate implementations of the equations.

## Run tests

```bash
python -m pytest -q tests
```

The packaged tests include analytic cloud/Kepler checks, a manufactured v4 run, complex source reassembly, partial-overlap amplitude closure, and read-only HTTP endpoint tests. The built-in browser UI was also rendered and exercised in an isolated manufactured-data harness. **No full real-data point was executed** and none of these tests establish numerical convergence, publication agreement or a physical resonance.

## Limitations and errors to avoid

- The completed run's `config.json` and `physics.h5` are authoritative; no new simulation state is saved by the dashboard.
- **Only the unperturbed** cloud density is visualized. The orbital-phase slider changes the companion's marker, not the cloud density.
- Modes with `Delta_g≤0` have no outgoing continuum record; zero and forbidden cells are not faked with small positive numbers on logarithmic axes.
- `m_f` and temporal `g` are not identified in eccentric motion.
- A bright channel at a single point does **not** prove resonance. Source strength, coherent radial cancellation, Coulomb response and phase-space kinematics must be distinguished. A parameter scan is a separate scientific study.
- Saved forcing and source profiles exist only at the original sampled radii. The analytic regular solution can be regenerated cheaply by comparison with production Green integrals, but **grid cumulative integrals are for interpretation only**.
- Run-to-run comparison and full convergence studies are separate, deliberately deferred work.
