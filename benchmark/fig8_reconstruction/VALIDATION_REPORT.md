# Standalone conversion — implementation and validation report

**Research object:** Xu Appendix-A small-`q` outgoing scalar ionization at one prescribed Newtonian eccentric orbit. **Scope boundary:** reorganize and test this existing benchmark; do not begin finite-`q` or relativistic dynamics, and defer numerical convergence.

## Dependency audit

The optimized package initially required nine external `01_atom.py` ... `09_response_channel_catalogue.py` physics modules under an adjacent `fig8_reconstruction/` directory. That broke clean-room reuse and real execution. This release supplies local source implementations in `physics/`, preserves the existing API used by the fast driver, includes all solver/dashboard code, configuration files and tests, and changes `--module-dir` default to that local folder.

**All nine local physics modules are independently authored interfaces matching the inspected equations, not verbatim copies of the external code.** The standalone-implementation property has been tested. Exhaustive line-by-line equivalence to the remote code has not been tested and is not claimed.

External packages (declared in `requirements.txt`): numpy, scipy, sympy, mpmath, h5py, plotly, pytest (tests only). All source code and frontend HTML/JS/CSS is included in the ZIP; the frontend Plotly runtime is provided from the installed plotly package, without CDN calls. No runtime imports are required from any external repository.

## Actual checks completed

- 22 pytest tests passed, including *six new tests with the real bundled physics modules*.
- Tests cover import paths for all 9 modules; radial-field cloud mass normalization; Newtonian bound energy and Kepler scales; signed second-half true anomaly; exact Gaunt selection; interior dipole cancellation; circular harmonic locking; eccentric Fourier conjugation; numerical-versus-analytic Coulomb Wronskian.
- A real **one-channel** computation ran through Fourier -> scalar source -> radial Green overlap -> outgoing amplitude -> energy flux -> two-file HDF5 serialization -> derived quantities -> dashboard analysis.
- For that deliberately small finite grid (`configs/smoke_real.json`), the batched and original-*method* adaptive paths produced normalized values `1.4405062424608437e-08` and `1.4405062424806146e-08`, respectively; relative difference about `1.37e-11`. **This only checks internal numerical route agreement, not agreement with the remote production source.**
- A separate real **multi-channel** calculation with `ell_f_max=2`, `g=[1,2]` completed 12 structurally eligible open channels across six `(ell_f,g)` response groups. The saved data reconstructed consistent channel rankings and Fig. 11 angular sums.
- Catalogue-only check for the configured full figure point (`ell_f_max=10`, ten positive harmonics): 1,210 candidate channels, 660 structurally eligible open channels.
- Dashboard API was exercised against actual one- and multiple-channel run outputs; previous frontend smoke checks were performed against manufactured data.

## Not established

- Full 660-channel (`ell_f_max=10`, `g=1..10`) Fig. 8 ordinate has **not** been executed here.
- Agreement with the original GitHub production implementation, printed Figure 8, or Figure 11 has **not** been established.
- Spatial grid, radial tail, harmonic, angular, quadrature and precision convergence have **not** been established; intentionally deferred.
- Overall performance speedup versus the original implementation has **not** been measured on an equal-environment full run.
- Leading-order cloud fractional ionization rate uses the dashboard's `Q_b=M_b/mu` assumption, which should be independently reconciled with paper's global charge/mass conventions before reporting an absolute rate.

## Source provenance

- User's prior v4 optimized driver and v1 interactive analysis dashboard files were used as the orchestration/dashboard starting point.
- Physics formulas and class interface behavior were reconstructed using inspection of the public repository `Jinghao333/Gravitational_Molecules`, `fig8_reconstruction/`.
- User-supplied Fig. 11 image used to specify angular grouping only, not as a numerical validation dataset.

## Stop boundary

Do not claim paper reproduction until separately validating the frozen model and numerical convergence. The next practical step is to run a medium-sized known config on the user's machine and compare local-vs-original production values at the same parameters and grid. A future convergence campaign remains distinct.
