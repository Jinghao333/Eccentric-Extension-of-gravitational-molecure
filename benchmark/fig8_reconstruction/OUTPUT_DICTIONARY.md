# Xu Fig. 8 — open-channel output contract (v4)

**Status:** implementation candidate. Offline orchestration and manufactured-data tests passed; comparison with the live original physics modules and convergence checks have **not** been run.

## 1. Run directory

Every production invocation creates a **new, non-existing** output directory. The effective configuration, including optional CLI engine overrides, is written **before** numerical computation. A successful run contains exactly:

```text
runs/
  xu_fig8_newtonian_alpha020_e04_p50_20261116T...Z/
    config.json                # sole saved input configuration
    physics.h5                # open-channel numerical results only
```

A failed execution preserves `config.json` and may write `FAILED.txt`, but will never present an incomplete HDF5 as a complete run. `--profile` explicitly adds `profile.pstats`; derived CSV exports are opt-in and should go **outside** the authoritative run directory.

## 2. Physics and notation

The input state is `|N_i,ell_i,m_i>` with `alpha=M mu`, cloud mass `M_b` and mass ratio `q=m_p/M`.

```math
E_i=-\frac{\mu\alpha^2}{2N_i^2},\quad
\Omega=\frac{(1-e^2)^{3/2}}{Mp^{3/2}},\quad
\Delta_g=E_i+g\Omega,\quad
\omega_g=\mu+\Delta_g.
```

An output channel is `c=(ell_f,m_f,g)`, with tidal forcing magnetic number `m=m_f-m_i`. Only **structurally admissible, spectrally open, nonzero-sampled** channels have saved source and outgoing-amplitude records. Forcing coefficients are saved only when they actually contribute to at least one such channel. Open channels with identically zero sampled source have small **keys-only** markers; closed, threshold and symmetry-forbidden channels are absent and can be reconstructed from config + selection rules.

The Fourier coefficient in the **coordinate-time basis** is

```math
u_{lmg}^{(M)}(r)=\frac1{2\pi}\int_0^{2\pi}
 u^T_{lm}(r,E)\,e^{ig(E-e\sin E)}(1-e\cos E)\,dE,
\qquad u_{lmg}^{(t)}=e^{ig\Omega t_p}u_{lmg}^{(M)}.
```

The scalar forcing source is

```math
S_{l,c}(r)=2\mu^2rR_i(r)\,
 C(l,m;\ell_i,m_i;\ell_f,m_f)\,u_{lmg}^{(t)}(r),
\qquad S_c(r)=\sum_l S_{l,c}(r).
```

The reduced Green-function overlap and **physical outgoing-wave coefficient** are

```math
I_{l,c}=\int_{r_{\min}}^{r_{\max}}
 Z_{-,\ell_fg}(r)S_{l,c}(r)\,dr,
\qquad I_c=\sum_l I_{l,c},
\qquad Z_c^{\infty}:=A_{\mathrm{out},c}
=e^{\pi\eta_g/2}\frac{I_c}{W_{\ell_fg}},
\quad \eta_g=\frac{M\mu^2}{\sqrt{2\mu\Delta_g}}.
```

> **Naming resolution:** this package defines the stored `Z_infinity` as the **complex physical outgoing amplitude** `A_out` of the existing Stage VI code (unit-Coulomb-phase normalization). It is a **scalar coefficient per open channel**, not the radial homogeneous outer solution `Z_+(r)` or the regular inner solution `Z_-(r)`. Both analytic radial solutions can be generated on demand; do not interchange them with `Z_infinity`.

The original source uses nonrelativistic Coulomb/Whittaker radial response and the Xu Klein–Gordon momentum in flux kinematics:

```math
k_{{\rm NR},g}=\sqrt{2\mu\Delta_g},\quad
k_{{\rm KG},g}=\sqrt{(\mu+\Delta_g)^2-\mu^2},\quad
\dot E^{s,\infty}_c
=2g\Omega\,k_{{\rm KG},g}\,\left|Z_c^\infty\right|^2,
\quad \mathcal F=\frac{\sum_c\dot E_c^{s,\infty}}{q^2(M_b/M)}.
```

## 3. `config.json`

Saved once per run. Contains all physical, numerical, truncation and Fourier-basis setup (`M,alpha,q,M_b_over_M,N_i,ell_i,m_i,p,e,t_p,varpi,ell_f_max,g_values,r_min,r_max,log_points,orbit_interval_points,forcing_*,response_*,mpmath_dps,forcing_engine,radial_engine`). **No derived orbital properties, grid array, cloud profile, normalization constants, or total flux are duplicated.**

CLI overrides are **resolved first**, then this *effective* config is registered. The stored config is sufficient to recreate the radial sample locations using the frozen grid recipe below.

## 4. `physics.h5`

| HDF5 field | Shape, type | Physics meaning | Why persisted |
|---|---|---|---|
| `/forcing/keys_lmg` | `(N_u,3)` int32 | Unique Fourier mode labels `(l,m,g)` | Index of expensive Fourier coefficients |
| `/forcing/u_lmg` | `(N_u,N_r)` complex128 | `u_lmg^(t)(r_j)` | Complex forcing phase and radial dependence; evaluated once for every *needed* `(l,m,g)` |
| `/channels/keys_ell_m_g` | `(N_c,3)` int32 | Computed open-channel labels `(ell_f,m_f,g)` | Direct indexing and matching of source, amplitude and overlaps |
| `/channels/S_total` | `(N_c,N_r)` complex128 | `S_c(r_j)=sum_l S_l,c(r_j)` | Requested coherent radial source; useful for plots and physical interpretation without regeneration |
| `/channels/Z_infinity` | `(N_c,)` complex128 | `Z_c^∞=A_out,c` | Requested physical coefficient entering ionization power; retains phase, not just squared modulus |
| `/channels/sampled_zero_keys` | `(N_0,3)` int32 | Structurally possible open channel with sampled `S=0` | Tiny integrity marker; distinguishes genuine sampled zero from an accidentally omitted result |
| `/overlap/keys_ell_m_g_l` | `(N_I,4)` int32 | `I_l,c` labels `(ell_f,m_f,g,l)` | Index of computed partial Green-function overlaps |
| `/overlap/I_l` | `(N_I,)` complex128 | `I_l,c=int Z_- S_l,c dr` | Required for **multipole interference** without repeating expensive Green integrations |

All mode-indexed arrays are sorted lexicographically by their corresponding key table. Each forcing coefficient is saved **once across all channels sharing `(l,m,g)`**. The `S_total`, `Z_infinity` and `I_l` arrays contain data only for *nonzero-sampled, open* channels. **`I_l` is the only additional expensive primitive**: the coherent total `I_c` would lose individual multipole phase contributions.

### Small HDF5 attributes (not per-sample output arrays)

* `schema`, `complete`, `grid_recipe`, `grid_sha256`, `n_radius`, `fourier_basis`, `Fourier_sign`, `spherical_harmonics`, `Z_infinity_definition`, `Z_infinity_kinematics`.
* `provenance_json`: original-module and pipeline code SHA-256 hashes, plus library versions. **These are not recoverable from numerical inputs**, and allow future users to detect code drift.
* `quality_json`: only `fourier_all_accepted`, maximum reported Fourier absolute error, maximum radial *group* error norm, count of adaptive fallbacks, and elapsed seconds. These are **diagnostic indicators, not certified uncertainty bounds**. No per-radius quadrature error/evaluation/status arrays or grouped-integration arrays are kept.

The radius grid is **not stored**. Reconstruct it from `config.json` as:

```python
r = np.unique(np.concatenate([
    np.geomspace(r_min, r_max, log_points),
    np.linspace(max(r_min, M*p/(1+e)),
                min(r_max, M*p/(1-e)), orbit_interval_points)
    # only if orbit_interval_points>0 AND clipped periapsis < apoapsis
]))
```

The package uses the shared `primitive_store.radial_grid(config)` implementation, tags recipe `geomspace_plus_inclusive_peri_apo_v1`, and saves only the SHA-256 hash of the numeric float64 grid to detect changed recipes or float rounding.

## 5. Why NOT separate output files per `(l,m,g)`?

`(l,m,g)` identifies an **orbital Fourier forcing**, not an outgoing channel: outgoing channels are `(ell_f,m_f,g)`, and several outgoing channels reuse one forcing coefficient. Storing one file per coefficient would create unnecessary repeated metadata and slow filesystem access for hundreds of small objects; storing the source and outgoing amplitude in forcing-named files would also hide their **different index spaces**.

Instead, use one HDF5 per run containing four indexed numerical tables, with small readers:

```python
from primitive_store import load_primitives, get_forcing, get_channel

run = load_primitives('runs/my_run')
r, u311 = get_forcing(run, 3, 1, 1)
ch = get_channel(run, ell_f=2, m_f=2, g=1)
r, S, Z_inf, I_l = ch['r'], ch['S'], ch['Z_infinity'], ch['I_l']
```

These functions return NumPy arrays directly; analysis exports should be generated **on demand**, outside the canonical run directory.

## 6. Derivable outputs (do not store them)

- `r_j`, `r_p`, `r_a`, `a`, `Omega`, `E_i`, atomic profile `R_i(r_j)` and the Kepler trajectory: directly from config and frozen model formulas.
- Partial sources `S_l,c(r_j)`: from unique `u_lmg`, `R_i`, and exact Gaunt coefficients; no partial-source arrays needed. `S_total` is retained **by explicit request** as a convenient computed result even though it can be regenerated from these ingredients.
- Inner/outer radial homogeneous functions `Z_-(r), Z_+(r)`, the analytic Wronskian `W_(ell_f,g)`, total overlap `I_c=sum_l I_l,c`, and `k_NR`, `k_KG`, `eta`, phases and spectral threshold distances: regenerate without new orbital Fourier quadrature.
- Field-energy, signed charge and orbital-energy fluxes, normalised channel fluxes, Fig. 8 ordinate, harmonic/angular/cumulative sum tables: compute using `Z_infinity` from `physics.h5` and `config.json`.
- Multipole interference `2Re(I_l,c I_l',c*)`, coherence `|sum_l I_l|/sum_l|I_l|`, and relative phases: from saved partial complex `I_l` values.

## 7. Completeness and limits

The reader rejects incompatible radial grid hashes, duplicated keys, inconsistent arrays, nonfinite values, and partial overlaps assigned to missing channels. The derived analysis also checks that **every** structurally eligible, open channel is either computed or explicitly marked as sampled zero; and that `Z_infinity` is consistent with the saved `I_l` values and the analytic Wronskian. The optional spatial diagnostic reconstructs the source from Fourier/Gaunt data and compares it with the stored `S_total`.

**Not included:** closed/bound response, distributed time-dependent scalar fields, a source away from the configured sampled domain (without new evaluation), or convergence proof. The code retains the Appendix-A small-`q`, Newtonian eccentric benchmark—not the future finite-`q` molecule.

**Note on phase:** `Z_infinity` is phase-convention-dependent. Always compare complex amplitudes only after confirming the same time origin, periapsis angle, harmonic sign and spherical-harmonic normalization.
