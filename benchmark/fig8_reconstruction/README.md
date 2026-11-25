# Fig. 8 reconstruction — standalone Xu Appendix-A research package

**Purpose.** Compute one point in Xu et al. (arXiv:2606.21439), Appendix A, in the prescribed Newtonian eccentric small-`q` model; save minimal open-channel primitives; analyze the physical system and outgoing-channel structure in an interactive dashboard.

**Status.** Independently bundled and offline-tested implementation candidate, **not** a convergence-validated or paper-matched Fig. 8 point. The nine physics modules in `physics/` are local **interface-compatible reconstructions** of the inspected original equations, *not* verbatim vendored copies or a demonstrated line-by-line equivalence with the GitHub production modules. The full `(ell_f_max,g_max)=(10,10)` physics job has not been run here. See `VALIDATION_REPORT.md`.

## Fully standalone source layout

```text
fig8_reconstruction/
├── physics/
│   ├── 01_atom.py                    # hydrogenic state, field normalization
│   ├── 02_kepler_orbit.py            # prescribed eccentric Kepler orbit
│   ├── 03_tidal_field.py             # direct and indirect tidal multipoles
│   ├── 04_fourier_forcing.py         # branch-aware Fourier decomposition
│   ├── 05_angular_coupling.py        # Gaunt and selection rules
│   ├── 06_radial_source.py           # coherent scalar radial source
│   ├── 07_green_response.py         # Coulomb/Whittaker Green solver
│   ├── 08_flux.py                   # physical outgoing amplitude/flux
│   └── 09_response_channel_catalogue.py # exact channel pruning
├── one_point_fast.py                # one-run scheduler, batched or reference modes
├── batched_quadrature.py            # harmonic-batched, grouped radial quadrature
├── primitive_store.py               # two-file run schema and validation
├── derive_quantities.py             # all derived power and channel diagnostics
├── radial_diagnostics.py            # optional on-demand partial radial profiles
├── compare_runs.py                  # compare equivalent-run outputs
├── analysis_dashboard/              # Python analysis API + local web dashboard
├── configs/
│   ├── fig8_point.json              # complete paper-point configuration
│   ├── smoke_real.json              # ONE channel, deliberately coarse smoke test
│   ├── multi_channel_smoke.json      # 12-channel, deliberately coarse smoke test
│   └── quick_regression.json         # smaller reference configuration
├── tests/                            # unit, orchestration, real solver and API tests
├── requirements.txt                  # only external runtime: standard PyPI libraries
├── OUTPUT_DICTIONARY.md              # exact scientific storage contract
└── VALIDATION_REPORT.md
```

No local directory outside this package is needed to run or analyze the physics. As usual, install the declared NumPy/SciPy/SymPy/mpmath/h5py/Plotly dependencies with `pip` (these are **not** vendored). Plotly.js is served locally from installed Plotly, not a CDN.

## Setup and tests

Python **3.11+**:

```bash
cd fig8_reconstruction
python -m pip install -r requirements.txt
python -m pytest -q tests
```

**One quick real calculation** (not a Fig. 8 reproduction or convergence test):

```bash
python one_point_fast.py --config configs/smoke_real.json --output-root runs
```

**Full requested figure point** (may be computationally expensive; not yet tested in full):

```bash
python one_point_fast.py --config configs/fig8_point.json --output-root runs
```

**Optimized-versus-reference comparison** on the same chosen configuration:

```bash
python one_point_fast.py --config configs/smoke_real.json \
  --forcing adaptive_cached --radial reference_per_channel --output-root runs
```

List structural channel counts without a physical solve:

```bash
python one_point_fast.py --config configs/fig8_point.json --catalogue-only
```

`--module-dir` is optional and defaults to the **bundled** `physics/`; it is retained solely for explicit regression testing against an alternate implementation. You never need to provide it to use this package.

## Output contract

Each completed invocation creates its own new directory under `runs/` and writes **exactly two authoritative files**:

```text
runs/<point_name>_<UTC_timestamp>/
    config.json        # effective inputs, registered before computation
    physics.h5         # unique Fourier u_lmg(r_j), S_c(r_j), A_out,c, partial I_l,c
```

The saved arrays are restricted to computed *open continuum* channels and their required forcing coefficients. The radial grid, cloud background, Gaunt coefficients, Wronskians, mode powers, normalized total, Fig. 11 aggregates and dashboard charts are **derived** from these inputs. `physics.h5` contains only small integrity/provenance attributes in addition to primary results. A failure never writes a completed `physics.h5`; failures inside the channel-evaluation loop also leave `FAILED.txt`.

Read and analyze without rerunning Fourier or Green integrals:

```bash
python derive_quantities.py runs/<your_run_folder>
python -m analysis_dashboard.server --runs runs --port 8765
```

Open `http://127.0.0.1:8765`. Pages: physical system with cloud/orbit scales; Fig. 11 / ranked channels / eccentric harmonic decomposition; selected-channel Fourier, source, Green-overlap and interference explorer. For an optional standalone radial diagnostic:

```bash
python radial_diagnostics.py runs/<your_run_folder> --top 5
```

## Physics conventions and what has not been verified

- `alpha=M*mu`, `q=m_p/M`, `E_i=-mu*alpha**2/(2*N_i**2)`.
- Initial `|211>` atom, cloud radial-field normalization `∫r²|R|² dr = M_b/(2mu²)`.
- Primary-centered prescribed Kepler orbit, `M_a=Ω(t-t_p)`, `φ_*=varpi+ν`; signed/atan2 true-anomaly branch.
- Fourier convention `u_lmg^(M)=(2π)^-1∫dM_a u_lm exp(+ig M_a)`; time-basis conversion `exp(+ig Ω t_p)`.
- Angular coupling `m=m_f-m_i`, parity/triangle/Wigner-3j, and circular-only `g=m` locking.
- `Delta_g=E_i+gΩ`; open continuum `Delta_g>0`; no epsilon shift at exact threshold.
- Coherent `l` sum inside each `(ell_f,m_f,g)`; noncoherent flux sum between channels.
- NR Coulomb radial response with retarded Whittaker branch; `A_out=exp(pi eta/2)*(∑I_l)/W`; Xu-compatible outgoing KG momentum in flux.
- No bound-sector radiation, no relativistic eccentric orbit, no finite-`q` dynamics.

**Scientific validation stop:** tests verify internal coherence and limited numerical pathway agreement. They do **not** establish equivalence to the original GitHub source or the numerical accuracy of the full paper figure. Full point agreement and convergence remain separately scheduled work.

Original inspected implementation: https://github.com/Jinghao333/Gravitational_Molecules/tree/main/fig8_reconstruction .
