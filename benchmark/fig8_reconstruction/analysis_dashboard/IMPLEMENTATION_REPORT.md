# Analysis dashboard — integrated standalone implementation note

The dashboard is a read-only interface for the two canonical simulation outputs `config.json` and `physics.h5`. It calls `analysis_dashboard/core.py`, `derive_quantities.py` and `primitive_store.py` from this **same project**. It does not need another checkout, and it does not save redundant derived quantities in the run directory.

Views: physical system and equatorial cloud density/orbit; angular Fig. 11 aggregation, top individual channels, harmonic spectrum and heatmap; selected-channel complex Fourier forcing, source decomposition, Green-response overlap and interference.

The Python frontend server is in `analysis_dashboard/server.py`. Plotly.js is loaded from the locally installed Plotly package, not from an external CDN. Backend/unit/API tests are included under `tests/`.

The root `VALIDATION_REPORT.md` is the source of truth for the **current** test results, scope and unverified conditions. It replaces historical v1 dashboard reports produced before the standalone conversion.
