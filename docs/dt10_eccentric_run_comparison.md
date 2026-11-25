# Fixed and geometry-controlled timesteps for the eccentric \(p\simeq50\) run

The historical calculation uses a fixed timestep of 10. The newer calculation distributes 288 solves over each complete post-startup orbit. Their mean spacings differ by 0.171%. The new method also integrates the startup clock analytically, samples the Hamiltonian at the step midpoint, and applies outer damping per unit physical time. This report distinguishes the measured phase improvement, controlled field-refinement results, and agreement between the two production runs.

| Run | Time integrator | Outer boundary |
|---|---|---|
| [Simulation_2026_11_15_16_53_10](../Simulation_2026_11_15_16_53_10/config.toml) | Legacy endpoint CN, fixed \(h=10\) | Per-step mask |
| [Simulation_2026_11_20_00_06_35](../runs/orbital-geometry-dt10-comparison/Simulation_2026_11_20_00_06_35/config.toml) | Midpoint CN, normalized geometry | Complex absorbing potential (CAP) |

## 1. The timestep function

Both simulations start at pericenter \(R_p=35.71\) on the same prescribed ellipse. With \(M_1=1,\ M_2=10^{-5},\ e=0.4\), the semimajor axis is \(a=R_p/(1-e)\), the mean motion is \(n=\sqrt{(M_1+M_2)/a^3}\), and \(P_0=2\pi/n=2884.9321039\). Here \(\theta\) is the primary's unwrapped azimuth and \(R\) is the instantaneous binary separation.

**During startup** (\(0\le t< T_r=2P_0\)), the new solver uses uniform physical-time nodes at multiples of \(P_0/N_s\), clipped at the ramp and final boundaries. Here \(N_s=288\), so startup takes 576 steps of \(10.01712536\).

**After startup**, let \(N_o=288\). The most direct statement of the implemented rule is that each complete step \([t_i,t_i+h_i]\) satisfies

$$
\boxed{
\frac12\left[
\frac{h_i}{P_0}
+
\frac{2|\theta(t_i+h_i)-\theta(t_i)|
+3\displaystyle\int_{t_i}^{t_i+h_i}\frac{|\dot R(s)|}{R(s)}\,ds}
{4\pi+6\ln[(1+e)/(1-e)]}
\right]=\frac1{N_o}.}
$$

The first term measures elapsed time. The second measures changes in the potential geometry. The angular weight 2 represents a quadrupole phase, and the radial weight 3 represents the \(R^{-3}\) dependence of an exterior quadrupolar tide. These are fixed sensitivity weights. The denominator is the geometry accumulated over one revolution: \(2|\Delta\theta|=4\pi\), and the outward and inward legs each contribute \(\ln(R_a/R_p)\) to the total variation of \(\ln R\). The absolute value inside the integral prevents cancellation across a turning point. A larger local rate \(2|\dot\theta|+3|\dot R|/R\) produces shorter steps.

**Implementation.** For a cycle beginning at \(t_c=T_r+kP_0\), the code sets \(u=(t-t_c)/P_0\) and evaluates the boxed left-hand side cumulatively from \(t_c\) to \(t\), obtaining \(C(u)\). It finds \(u_j=C^{-1}(j/N_o)\) by bisection and takes \(h_j=P_0(u_{j+1}-u_j)\). For the radial integral it splits the orbit at each apsis and sums \(|\ln(R_{\rm end}/R_{\rm start})|\) over the monotonic pieces. It anchors the endpoints to absolute cycle times and clips the final step to the requested physical end time. The two resolution knobs are **startup_steps_per_period** and **orbit_steps_per_period**. Holding the first fixed while increasing the second refines only the post-startup orbit. The earlier tolerance-based **orbital_geometry** controller is a separate option; this run selected **normalized_geometry**. [Scheduler implementation](../src/CrankNicolsonModes_evolution.cpp)

**Measured example.** In the first full orbit after startup, the step leaving pericenter is \(7.08532\). One of the shortest steps is \(6.56232\), ending at \(R=37.4546\) on approach to pericenter. The step reaching apocenter is \(15.00149\), ending at \(R=83.3233\). The radial-speed term explains why the smallest step can occur slightly away from pericenter. All 288 accepted durations sum to \(2884.9321039\). The historical grid averages \(P_0/10=288.4932\) steps per orbit. These counts compare linear solves; midpoint propagation also requires an endpoint operator rebuild.

## 2. The startup phase and Hamiltonian evaluation

The prescribed motion obeys \(\dot{\mathcal M}=ng(t)\), with \(g=3x^2-2x^3\) and \(x=\mathrm{clamp}(t/T_r,0,1)\). Both runs solve \(E-e\sin E=\mathcal M\) and use \(R=a(1-e\cos E)\) and the Kepler true anomaly for \(\theta\). Both radial and angular velocities carry the factor \(g\). The companion's potential has a separate weight \(\lambda=g\) for these configurations, giving \(V=V_1+\lambda V_2\).

The historical run uses \(\mathcal M_{i+1}=\mathcal M_i+nh_i g(t_i)\). The new run uses \(\mathcal M(s)=\mathcal M_i+n[F(s)-F(t_i)]\), where

$$
F(t)=\begin{cases}
T_r[(t/T_r)^3-\tfrac12(t/T_r)^4],&0\le t\le T_r,\\
t-T_r/2,&t\ge T_r.
\end{cases}
$$

The code evaluates the interval integral in factored form and handles a step crossing the ramp endpoint. At \(T_r\) the exact clock has accumulated one full revolution. In saved post-startup diagnostics through \(95P_0\), the old mean anomaly lags the analytic clock by \(0.01088966\) rad, equivalent to about five time units. The newer run's maximum absolute difference is about \(6.0\times10^{-11}\) rad. This is a direct production measurement of improved orbital timing. [Clock implementation](../src/CrankNicolsonModes_evolution.cpp) · [Kepler geometry](../src/FixedKeplerOrbit.cpp)

Both field solvers use

$$
(I+i h_iH_*/2)\psi_{i+1}=(I-i h_iH_*/2)\psi_i.
$$

The historical run takes \(H_*=H(t_i+h_i)\). The new run evaluates the midpoint and endpoint orbits independently from the same starting state. It propagates with \(H_*=H(t_i+h_i/2)\), including the midpoint companion weight, then rebuilds the endpoint Hamiltonian for the field derivative, energy, flux and binary-power diagnostics. It uses the captured \(h_i\) throughout the solve. Endpoint sampling generally gives first-order time accuracy for varying \(H\); midpoint sampling gives second-order accuracy for sufficiently smooth \(H\). Radial multipoles can change analytic branch when a source crosses a radial grid node. The achieved order must therefore be checked on the quantity and mesh of interest. [Timestep workflow](../src/CrankNicolsonModes.cpp) · [CN solve](../src/CrankNicolsonModes_matrices.cpp)

## 3. Outer damping, physical controls and reproducibility

- **Outer damping.** Both layers start at \(r_b=12750\), outside \(r_{\rm ex}=10000\). Write \(x=(r-r_b)/2250\) inside the layer. The old mask multiplies the field by \(\exp(-x^6)\) after every step of 10. The new Hamiltonian includes \(-i\gamma\), with \(\gamma=0.1x^6\). Its exact local attenuation over ten time units is \(\exp(-10\gamma)=\exp(-x^6)\). This matches the continuous amplitude-decay rate. For pure damping the discrete CN multiplier is \((1-h\gamma/2)/(1+h\gamma/2)\); transport through the layer also affects reflection. Equality of reflected waves has not been established. [Boundary code](../src/CrankNicolsonModes_outgoing.cpp)
- **Physical controls.** Both prescribe moving fixed-mass orbits with feedback and horizon absorption disabled. Both configure a \(211\) cloud with \(\mu=0.2,\ M_c=1\), 10,000 radial nodes over \(L=15000\), and 28 angular modes through \(\ell=6\). They share the multipole coupling, radial discretization and sparse CN solve. The new configuration requests eight threads; the old one inherits the environment. The eccentric orbital work term \(\dot R\,\partial_R V+\dot\theta\,\partial_\theta V\) and separate switching work \(\dot\lambda V_2\) were already present in the old model. [Hamiltonian construction](../src/CrankNicolsonModes_hamiltonian.cpp)
- **Timing and restart.** The old driver converts the requested periods to a fixed sample count using a floor. The new driver targets \(200P_0\) directly, records each accepted duration in **dt_step.txt**, and reconstructs the anchored grid on resume. Checkpoint validation includes integrator and physics settings. Both runs request snapshots every 100 solves, so the physical time between new snapshots varies after startup. [Driver](../src/main.cpp) · [Restart implementation](../src/CrankNicolsonModes_restart.cpp)

## 4. What the convergence comparisons establish

For complex-field refinement at common final time, define \(Q=\|\psi_h-\psi_{h/2}\|_2/\|\psi_{h/2}-\psi_{h/4}\|_2\). Second-order behavior gives \(Q\to4\). The archived midpoint fixture gives \(Q=3.9698\) for \(e=0\) and \(Q=3.96932\) for \(e=0.4\). The separate normalized-grid fixture, with startup held fixed and orbital counts 60, 120 and 240, gives successive field differences \(0.0101243\) and \(0.00637388\), or \(Q=1.58841\). It demonstrates improvement on refinement in that test; the measured ratio falls below four. [Midpoint test](../tests/midpoint_cn.cpp) · [Normalized mesh test](../tests/geometry_timestep.cpp)

For the production observable, average the signed dashboard ionization power at the extraction sphere with \(\langle S\rangle=(b-a)^{-1}\int_a^b S(t)\,dt\). Use physical-time trapezoidal weights and interpolate at the two window boundaries. An unweighted sample average would give excessive weight to short steps. Negative power denotes outgoing ionization in this convention.

| Mean signed power | Fixed step + mask | Normalized geometry + CAP | Change in magnitude |
|---|---:|---:|---:|
| Total, \(0\!-\!80P_0\) | \(-6.502515\times10^{-6}\) | \(-6.521995\times10^{-6}\) | +0.2996% |
| Total, \(90\!-\!95P_0\) | \(-2.823045\times10^{-6}\) | \(-2.819010\times10^{-6}\) | −0.1429% |
| \((1,1)\), \(90\!-\!95P_0\) | \(-6.989280\times10^{-7}\) | \(-6.943440\times10^{-7}\) | −0.6559% |
| All other modes, \(90\!-\!95P_0\) | \(-2.124117\times10^{-6}\) | \(-2.124666\times10^{-6}\) | +0.0258% |

The strongest production result is the removal of the measured startup phase bias. The two power windows show agreement at nearly equal solve density, with the largest late-window difference in the \((1,1)\) channel. Determining which timestep distribution yields smaller production power error requires refinement at common physical times with the same midpoint propagator, startup treatment and CAP on both grids.
