# Local Xu Appendix-A physics modules

The nine numbered modules provide the API required by `one_point_fast.py` and the optional `radial_diagnostics.py`, with **no imports from another repository**. The labels follow the original Appendix-A stages: unperturbed atom, orbital trajectory and tidal forcing, angular projection, radial response, outgoing amplitude and flux, and candidate channel selection.

These modules are **local reconstructions following the inspected physics formulas**. They are not copied byte-for-byte from the external GitHub project. Therefore the code is self-contained, and its internal formulas have analytic/unit/small-run regression checks, but cross-implementation equality to the original source is not a proven property of this release.

## Model/approximation ledger

| Scheme | Control / inputs | Retained | Omitted | Validity |
|---|---|---|---|---|
| Companion small-`q` | `q=m_p/M`, `M,alpha,M_b,p,e` fixed in a formal q->0 comparison | leading scalar source O(q), outgoing power O(q²) | nonlinear cloud response O(q²), power corrections beyond q² | q << 1, prescribed orbit on orbital timescales |
| Nonrelativistic atom | `alpha=M mu` small, `N_i,ell_i,m_i` fixed | Coulomb hydrogenic bound profile E_i=-mu alpha²/(2 N_i²) | fine structure, relativistic gravitational corrections | alpha <<1 (formal expansion; alpha=.2 is a benchmark parameter) |
| Kepler orbit | `p,e,M` fixed, e<1 | Newtonian ellipse with Ω=(1-e²)^(3/2)/(Mp^(3/2)) | periapse precession, PN frequency splitting, radiation reaction | weak field/slow motion; separately from small-q |
| Linear response | initial state fixed, forcing O(q) | Green/Whittaker inhomogeneous scalar response | backreaction, self-interactions, finite-q two-center Hamiltonian | leading perturbative cloud response |
| Output flux kinematics | selected open Delta_g>0 | Xu KG momentum in asymptotic current | relativistic correction to radial operator | mixed imported prescription; not pure NR flux |
| Numerical cutoffs | config ell_f_max, g_values, radial grid, quadrature | finite sums and finite domain | all excluded modes/tails | requires later convergence validation |

## Principal equations

`02_kepler_orbit.py`: `r_*(E)=a(1-e cos E)`, `M_a=E-e sin E`, signed true anomaly through `atan2`; `a=pM/(1-e²)`.

`03_tidal_field.py`: primary-centered direct multipoles plus indirect `l=1` term. Interior dipole cancels exactly; exterior retains both contributions.

`04_fourier_forcing.py`: `u_lmg^(M)(r)=(2π)^-1∫_0^(2π) dE (1-e cos E)u_lm(r,E) exp[+ig(E-e sin E)]`; `u^(t)=exp(+igΩt_p)u^(M)`.

`05_angular_coupling.py`: `C=∫Y_lm Y_li,mi Y_lf,mf* dΩ`; exact m, triangle and parity filters with Wigner 3j.

`06_radial_source.py`: `S_l,c(r)=2mu²r R_i(r) C_l,c u_lmg^(t)(r)`; sum complex amplitudes coherently.

`07_green_response.py`: Coulomb regular `M_{κ,lf+1/2}(λr)`, retarded exterior `W_{κ,lf+1/2}(λr)`, `λ=-2ik_NR` in the open sector, `κ=2Mmu²/λ`, `W_r=-λ Γ(2lf+2)/Γ(lf+1-κ)`, `I_l,c=∫dr Z_-(r)S_l,c(r)`.

`08_flux.py`: `A_out,c=exp(π η_g/2)∑_l I_l,c/W`; `k_KG=sqrt[(mu+Delta_g)²-mu²]`; `dot E^s_c=2g Ω k_KG |A_out,c|²`.

`09_response_channel_catalogue.py`: structurally prune exact angular/equatorial zeros, circular temporal zeros, missing degrees and spectral closed/threshold modes before numerical source construction.

For rigor, no physical resonance claim should be drawn from a bright mode at one simulation point. Selection rules and analytic normalization are tested, while full paper figure agreement remains open.
