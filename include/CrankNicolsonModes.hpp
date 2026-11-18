/**
 * @file CrankNicolsonModes.hpp
 * @brief Crank–Nicolson evolution of a scalar cloud in a binary black-hole
 *        potential, expanded in spherical-harmonic modes about the centre of mass.
 *
 * The field is written as
 * @f[
 *   \psi(t, r, \theta, \phi) = \sum_{\ell m} \frac{U_{\ell m}(t, r)}{r}\, Y_{\ell m}(\theta, \phi),
 * @f]
 * and the radial functions @f$U_{\ell m}@f$ of all retained modes are stacked
 * into one vector of length @c Nr * n_modes (mode-major: entry @c j*Nr+i is
 * mode @c j at radial node @c i). The Hamiltonian is the free radial operator
 * plus the multipole expansion of the two Newtonian potentials, which couples
 * the angular modes. Optional ingredients are an s-wave absorber on each black
 * hole, a complex absorbing potential (CAP) or damping mask at the outer
 * boundary, and back-reaction of the cloud on the orbit.
 */
#ifndef CRANK_NICOLSON_MODES_HPP
#define CRANK_NICOLSON_MODES_HPP

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <csignal>
#include <complex>
#include <cstddef>
#include <fstream>
#include <string>
#include <utility>
#include <vector>
#include "Definitions.hpp"
#include "FixedKeplerOrbit.hpp"
#include "OrbitState.hpp"
#include "RollingPhase.hpp"
#include "RadialBlockPreconditioner.hpp"
#include "Utils.hpp"
#include "Config.hpp"

/// Row-major complex sparse matrix, so Eigen can parallelise sparse products.
using SparseComplex = Eigen::SparseMatrix<dcomplex, Eigen::RowMajor>;

/**
 * @brief Mode-decomposed Schrödinger solver for a gravitational atom in a binary.
 *
 * Units are @f$G = c = 1@f$ with the primary mass @f$M_1 = 1@f$ at start. The
 * primary sits at distance @c R1 and the companion (mass @c q) at @c R2 from
 * the centre of mass, separated by @c R.
 *
 * Typical use:
 * @code
 *   CrankNicolsonModes cn(N, L, dt, epsilon, l_max, mu, M_c, q, a, R, counter,
 *                         absorption, save_every, ramp_periods, epsilon_smear,
 *                         "config.toml");
 *   cn.set_initial_state_full(load_initial_modes_binary("initial_modes.bin", cn.modes));
 *   cn.initialize();
 *   cn.run(num_steps, r_ex);
 * @endcode
 *
 * Each step solves @f$(1 + i\,dt\,H_{\rm full}/2)\,\psi^{n+1} = (1 - i\,dt\,H_{\rm full}/2)\,\psi^n@f$
 * with BiCGSTAB, where @f$H_{\rm full} = H - i\Gamma/2@f$ and @f$\Gamma@f$ is the
 * rank-one s-wave sink of each black hole (kept outside the sparse @ref H).
 */
class CrankNicolsonModes
{
public:
    enum class OrbitMode { Isolated, FixedCircular, FixedEccentric, CircularInspiralLegacy };
    enum class TimeStepRule { FixedDt, LegacyCircularInspiral };

    /// Derived from existing settings; neither mode is independently mutable.
    OrbitMode orbit_mode() const noexcept;
    TimeStepRule timestep_rule() const noexcept;

    /**
     * @brief Build the grid, the mode list and the free Hamiltonian.
     *
     * @param Nx            Number of radial nodes (>= 3).
     * @param L             Outer radius of the radial domain.
     * @param dt            Initial timestep (kept until the startup ramp ends).
     * @param epsilon       Offset of the first radial node, @f$r_i = \epsilon + i\,dr@f$.
     * @param l_max         Largest retained multipole @f$\ell@f$.
     * @param mu            Boson mass @f$\mu@f$.
     * @param M_c           Cloud mass; the initial state is scaled to it.
     * @param q             Mass ratio @f$M_2/M_1@f$.
     * @param a             Primary spin (stored; not used by the evolution).
     * @param R             Initial binary separation.
     * @param counter       @c true for a counter-rotating orbit.
     * @param absorption    Enable the s-wave black-hole absorbers.
     * @param save_every    Snapshot / checkpoint interval in steps.
     * @param ramp_periods  Length of the startup ramp in initial orbital periods.
     * @param epsilon_smear Absorber support radius in units of each body's mass.
     * @param config_file   Configuration copied into the output folder.
     * @throws std::invalid_argument for an invalid grid or timestep.
     */
    CrankNicolsonModes(int Nx, double L, double dt, double epsilon, int l_max,
                       double mu, double M_c, double q, double a, double R,
                       bool counter, bool absorption, int save_every,
                       double ramp_periods, double epsilon_smear,
                       std::string config_file, double eccentricity = 0.0);

    /// @name Grid, modes and physical parameters
    /// @{
    int Nr;          ///< Number of radial nodes.
    int save_every;  ///< Snapshot interval in steps.
    int l_max;       ///< Largest retained multipole.
    int n_modes;     ///< Number of retained (l, m) modes.
    double epsilon;  ///< Radial offset of the first node.
    double mu;       ///< Boson mass.
    OrbitState orbit; ///< One authoritative instantaneous orbital state.
    double a;        ///< Primary spin (informational).
    double L;        ///< Outer radius.
    double dr;       ///< Radial spacing.
    double dt;       ///< Timestep for the next call to step(); restored verbatim on resume.
    bool counter;    ///< Counter-rotating orbit.
    double M_c;      ///< Cloud mass.
    /// @}

    /**
     * @brief Advance the state by one timestep and append all diagnostics.
     * @param step_id Index of the time sample produced by this step.
     */
    void step(int step_id);

    /**
     * @brief Start a new run in a fresh @c Simulation_<date> folder.
     * @param n_steps Number of time samples, including @f$t = 0@f$.
     * @param r_ex    Extraction radius for fluxes and inner-region integrals.
     */
    void run(int n_steps, double r_ex);

    /**
     * @brief Continue the run stored in @p folder from @c checkpoint.bin.
     *
     * The binary checkpoint is the authoritative physical state. Diagnostics and
     * snapshots beyond its completed step are removed before output is appended.
     */
    void resume(const std::string &folder, int n_steps, double r_ex);

    /**
     * @brief Poll an async-signal-safe stop flag between completed timesteps.
     *
     * The owner must keep @p flag alive for the duration of run() or resume().
     * A signal handler may safely request a stop by assigning 1 to the flag.
     */
    void set_stop_flag(const volatile std::sig_atomic_t *flag) noexcept
    {
        stop_flag = flag;
    }

    /// Feed the binary power and accretion back onto the orbit. Absorption of
    /// the field is controlled separately by the constructor's @c absorption.
    bool evolve_orbit = true;

    /// The startup ramp always slows the orbital clock from zero to its full
    /// rate. When true it also ramps the companion's potential; when false the
    /// companion is present at full strength from @f$t = 0@f$.
    bool ramp_secondary_potential = true;

    /// With orbital evolution, the post-ramp timestep is @c P / steps_per_period.
    double steps_per_period = 50.0;

    /**
     * @brief Set the number of orbits in the rolling phase average of @c p_bin.
     *
     * Call before run(); a resumed run keeps the window stored in its checkpoint.
     * @throws std::invalid_argument unless @p n_orbits > 0.
     */
    void set_power_average_orbits(double n_orbits);

    /// Use the radial block preconditioner (default) instead of a full ILUT.
    bool use_block_preconditioner = true;
    double last_solver_residual = 0.0;  ///< Relative residual of the last CN solve.
    int last_solver_iterations = 0;     ///< BiCGSTAB iterations of the last CN solve.

    Eigen::VectorXd t;    ///< Time of each sample.
    Eigen::VectorXd r;    ///< Radial nodes.
    Eigen::VectorXcd psi; ///< Current state, stacked radial functions U_lm.
    bool mask = false;    ///< A damping mask is applied after every step.
    bool cap = false;     ///< A complex absorbing potential is part of H.

    /// Fixed eccentricity (zero for every legacy circular and isolated run).
    double eccentricity() const noexcept { return fixed_kepler.eccentricity(); }
    /// Kepler period of the prescribed ellipse, or current circular period for inspiral.
    double orbital_period() const noexcept;
    /// Initial/reference period used by startup, absorption and run-length schedules.
    double initial_orbital_period() const noexcept { return initial_period; }
    /// Mean motion of the fixed Kepler trajectory (not instantaneous angular speed).
    double mean_motion() const noexcept { return fixed_kepler.mean_motion(); }

    /**
     * @brief Set the initial state from all modes, scaling it to the cloud mass.
     * @param psi_0 Unit-norm stacked state of size @c Nr * n_modes, ordered as @ref modes.
     * @throws std::invalid_argument if the size does not match.
     */
    void set_initial_state_full(Eigen::VectorXcd psi_0);

    /// Retained modes; only @f$\ell - m@f$ even (equatorial symmetry) are kept.
    std::vector<mode> modes;

    /**
     * @brief Build the coupling caches, the initial potential and the CN matrices.
     *
     * Call after the grid and modes exist and before step() or run(). The sparse
     * matrices must not be structurally modified afterwards, because their
     * compressed offsets are cached.
     */
    void initialize();

    /**
     * @brief Replace the binary by one stationary mass-1 black hole at (x, 0, 0).
     *
     * Call before initialize(). There is no rotation or companion, and the
     * constructor's @c epsilon_smear is then the physical absorber radius.
     * @throws std::invalid_argument if the hole lies outside the domain.
     * @throws std::logic_error if called after initialize().
     */
    void configure_isolated_black_hole(double x);

    /// @name Operators
    /// @{
    SparseComplex H;    ///< Conservative Hamiltonian H0 + V (+ CAP).
    SparseComplex H0;   ///< Free radial Hamiltonian, including the centrifugal term.
    SparseComplex A;    ///< CN left operator 1 + i dt H / 2.
    SparseComplex B;    ///< CN right operator 1 - i dt H / 2.
    SparseComplex V1;   ///< Potential of the primary.
    SparseComplex V2;   ///< Potential of the companion.
    SparseComplex V;    ///< Active potential V1 + ramp * V2.
    SparseComplex dVdt; ///< Orbital derivative Rdot d_R V + theta_dot d_theta V.
    SparseComplex dVdt_switch; ///< Artificial startup switching derivative lambda_dot V2.
    SparseComplex Id;   ///< Identity, used to fix the union sparsity pattern.
    SparseComplex CAP;  ///< Diagonal complex absorbing potential.
    /// @}

    /// Apply the absorption operator @f$\Gamma = \sum_a g_a |v_a\rangle\langle v_a|@f$ to @p state.
    Eigen::VectorXcd apply_sink(const Eigen::VectorXcd &state) const;

    /// Return @f$-i H_{\rm full}\,\psi = -i H \psi - \Gamma\psi/2@f$, the Schrödinger right-hand side.
    Eigen::VectorXcd schrodinger_rhs(const Eigen::VectorXcd &state) const;

    /// Largest eigenvalue @f$g_a |v_a|^2@f$ of one body's sink (@p body = 1 or 2); zero when inactive.
    double sink_strength(int body) const;

    /// True once the absorbers are part of the evolution operator.
    bool absorption_sink_active() const { return sink_active; }

    /**
     * @brief Build a multiplicative damping mask for the outer region.
     *
     * For @f$r \ge r_m = f L@f$ every mode is multiplied after each step by
     * @f$\exp[-W ((r - r_m)/((1-f) L))^n]@f$.
     * @param mask_L_factor Fraction @f$f@f$ of the domain where the mask starts.
     * @param mask_W        Strength @f$W@f$.
     * @param mask_n        Power @f$n@f$.
     */
    void build_damping_mask(double mask_L_factor, double mask_W, int mask_n);

    /**
     * @brief Build a diagonal complex absorbing potential for the outer region.
     * @param cap_L_factor Fraction of the domain where the CAP starts.
     * @param cap_W        Strength.
     * @param cap_n        Power of the profile.
     */
    void build_cap(double cap_L_factor, double cap_W, int cap_n);

private:
    bool isolated_black_hole = false;
    double ramp_periods;
    /// Initial Kepler period. Startup schedules never use instantaneous angular speed.
    double initial_period;
    FixedKeplerOrbit fixed_kepler;
    /// Time at which the startup ramp ends and the orbit may evolve.
    double ramp_end_time() const { return ramp_periods * initial_period; }

    Eigen::VectorXcd psi_last; ///< State at the start of the current step.
    Eigen::VectorXcd damping_mask;

    std::string config_file;

    double ramp_factor = 0.0; ///< Orbital clock's previous endpoint factor, used by advance_orbit.

    /// Temporary operator inputs, not persistent orbit/clock state. Together
    /// with OrbitState these specify the evaluation for the fixed grid/physics
    /// configuration. In particular, absorption has its own later activation.
    struct HamiltonianEvaluation {
        double time;
        double companion_weight;
        double companion_weight_derivative;
        double absorption_weight;
        bool absorption_active;
    };
    HamiltonianEvaluation hamiltonian_evaluation(double time) const;

    double r_ex;    ///< Extraction radius.
    int idx_ex;     ///< Radial index of the extraction radius.

    std::string output_folder;

    double t_i = 0.0; ///< Current simulation time.
    const volatile std::sig_atomic_t *stop_flag = nullptr;

    Eigen::BiCGSTAB<SparseComplex, RadialBlockPreconditioner> solver;
    Eigen::BiCGSTAB<SparseComplex, Eigen::IncompleteLUT<dcomplex>> ilut_solver;

    /**
     * @brief Solve the CN system @f$A_{\rm full}\psi = \mathrm{rhs}@f$ into @c psi.
     *
     * Uses the radial block preconditioner and falls back to full ILUT. When
     * the sinks are active, the preconditioner is corrected for them with the
     * Woodbury identity. Every accepted solution satisfies a relative
     * residual of at most 1e-9.
     * @throws std::runtime_error if no preconditioner reaches the tolerance.
     */
    void solve_CN(const Eigen::VectorXcd &rhs, double step_dt);
    /// CN assembly, solve and discrete mask using only the captured step duration.
    void propagate_field(double step_dt);

    double orbits_average = 2.0;
    RollingPhaseAverage power_average{orbits_average};

    /// Radial factor of the primary's multipole @p l_star, @f$r_<^{\ell}/r_>^{\ell+1}@f$ with mass weight.
    double F1(const OrbitState &state, double r, int l_star, int m_star) const;
    /// Radial factor of the companion's multipole @p l_star, including its @f$(-1)^{m}@f$ parity.
    double F2(const OrbitState &state, double r, int l_star, int m_star) const;
    double dF1_dR(const OrbitState &state, double r, int l_star) const;
    double dF2_dR(const OrbitState &state, double r, int l_star) const;

    void build_modes(); ///< Fill @ref modes with all (l, m), l <= l_max, l - m even.
    void build_grid();  ///< Fill @ref r with the radial nodes.

    /// Build @ref H0 with second-order finite differences plus the centrifugal term.
    void build_stationary_hamiltonian();

    /// Precompute the Gaunt coefficients of every coupled mode pair.
    void build_coupling_cache();
    /// Evaluate V1, V2, V and dVdt for the current orbit.
    void build_static_perturbation(const OrbitState &state, const HamiltonianEvaluation &evaluation);
    /// Build potentials, absorbers and H for the already-updated endpoint state.
    void build_time_dependent_perturbation(const OrbitState &state, const HamiltonianEvaluation &evaluation);
    /// Advance the orbital clock using the old rate and evaluate endpoint state/rates.
    void advance_orbit(double step_dt, double endpoint_time);
    /// Evaluate startup weights and endpoint velocities without assembling operators.
    void update_orbit_and_ramp(double time, bool evaluate_geometry);
    void validate_orbit_configuration() const;

    /// One multipole contribution to a mode pair.
    struct AngularTerm
    {
        int l_star;       ///< Multipole order of the potential.
        dcomplex angular; ///< Angular (Gaunt) coefficient at the equator.
    };

    /// All multipoles coupling mode @c row_mode to mode @c col_mode.
    struct ModePair
    {
        int row_mode;
        int col_mode;
        int m_star; ///< Azimuthal number of the coupling, @c m_col - m_row.
        std::vector<AngularTerm> terms;
        /// Positions in the shared, compressed perturbation sparsity pattern.
        std::vector<Eigen::Index> offsets;
    };

    std::vector<ModePair> mode_pairs;
    std::vector<double> radial_f1, radial_f2;    ///< Cached F1, F2 indexed l*Nr+i.
    std::vector<double> radial_df1, radial_df2;  ///< Their derivatives with separation R.
    std::vector<dcomplex> stationary_values;     ///< H0 values in H's sparsity pattern.
    std::vector<Eigen::Index> potential_offsets; ///< Map from H slots to V slots (-1 if none).
    /// Refill H in place from H0, V and CAP; @p include_absorption activates the sinks.
    void assemble_hamiltonian(bool include_absorption);
    /// Refill A and B in place from H and the captured step duration.
    void update_CN_values(double step_dt);

    /// Build the fixed sparsity patterns of V, H, A, B and their offset caches.
    void build_sparse_element_caches();

    /// @f$e^{i m \theta}@f$ for the current orbital angle.
    dcomplex orbital_phase(int m_star, const OrbitState &state) const;

    /// @name Absorption
    /// @{
    bool absorption;
    /// Legendre moments of one body's spherical kernel about the origin.
    struct AbsorptionProfile {
        double radius = -1.0, width = -1.0;
        /// @f$2\pi\int K(s)\,P_\ell(\cos\gamma)\,d\cos\gamma@f$, indexed l*Nr+i.
        std::vector<double> multipoles;
    };
    AbsorptionProfile absorption_profile1, absorption_profile2;
    /**
     * @brief Recompute a body's kernel moments if its position or width changed.
     *
     * The kernel is the compact bump @f$\exp[-1/(1 - s^2/h^2)]@f$ in the
     * distance @f$s@f$ from the body, normalised to unit volume on the grid.
     */
    void update_absorption_profile(AbsorptionProfile &profile, double radius, double width);
    double Ah1 = 0.0, Ah2 = 0.0; ///< Horizon-area factors 16 pi M_a^2 (times the absorption ramp).
    double epsilon_smear;
    /**
     * Rank-one s-wave sink of one body: @f$\Gamma = g\,v v^\dagger@f$, where
     * @f$dr\,v^\dagger U = \int K(|x - x_{BH}|)\,\psi\,d^3x@f$ and @f$g = A_h\,dr@f$.
     * It only sees the @f$\ell = 0@f$ component of @f$\psi@f$ about the black hole.
     */
    struct SwaveSink {
        Eigen::VectorXcd v;
        double strength = 0.0;
    };
    SwaveSink sink1, sink2;
    bool sink_active = false;
    int mode00 = -1; ///< Index of the (0, 0) mode, or -1.
    /// Build the projector vector @p v of one body from its kernel moments.
    void build_swave_projector(const AbsorptionProfile &profile, bool secondary,
                               const OrbitState &state, Eigen::VectorXcd &v) const;
    /// Update both absorbers for the current orbit and absorption ramp.
    void build_spherical_absorption(const OrbitState &state, const HamiltonianEvaluation &evaluation);
    /// Mass accretion rate of the primary for @p state.
    double get_absorption_rate_1(const Eigen::VectorXcd &state) const;
    /// Mass accretion rate of the companion for @p state.
    double get_absorption_rate_2(const Eigen::VectorXcd &state) const;
    /// @}

    /**
     * @brief Evolve the orbit over @p orbit_dt.
     *
     * The separation shrinks by gravitational-wave emission and by the
     * phase-averaged cloud power @c p_bin; the masses grow by accretion. The
     * circular angular speed follows the new separation. Timestep selection is separate.
     */
    void apply_legacy_circular_feedback(double orbit_dt);
    /// Post-field stage, including the exact legacy activation/crossing gate.
    void apply_orbital_feedback(double step_start);
    /// Select only the next step's dt, after feedback; never called during resume.
    double select_next_timestep() const;

    double p_bin = 0.0; ///< Phase-averaged binary power.

    /// Temporary physical measurements, with no orbital state or rolling history.
    struct StepObservables {
        std::vector<double> fluxes, energy_flux, norms, energies;
        double energy, binary_power, switching_power, dV_bin_dt;
    };
    StepObservables compute_step_observables(const Eigen::VectorXcd &dpsi_dt);
    /// Exactly once per completed step, after feedback has used the previous p_bin.
    void update_rolling_power(double instantaneous_power);
    void report_step_progress(int step_id) const;
    void write_step_output(int step_id, const StepObservables &observables);

    /// Allocate A and B and configure both iterative solvers.
    void build_CN_matrices();

    /// Write @p psi to @c output/name/snapshot_NNNNNN.txt as "re im" lines.
    void save_snapshot(int step_id, const Eigen::VectorXcd &psi, std::string name);
    /// Write config.toml, r.txt and modes.txt to the output folder.
    void save_run_info();
    /// Flush every diagnostic stream.
    void flush_diagnostics();

    /// @name Restart support
    /// @{
    /// Set @c r_ex and @c idx_ex, validating them against the grid.
    void set_extraction_radius(double r_ex_, int n_steps);
    /// Open every diagnostic file in the output folder with @p mode.
    void open_diagnostics(std::ios::openmode mode);
    /// Run steps [@p first_step, @p n_steps), then flush and checkpoint the last completed step.
    void evolve(int first_step, int n_steps);
    /// Every diagnostic stream with its file name.
    std::vector<std::pair<std::ofstream *, std::string>> diagnostic_files();
    /// Write the exact binary state to checkpoint.bin, replaced atomically.
    void save_checkpoint(int step_id) const;
    /// Restore the exact state from checkpoint.bin and return its completed step, or -1.
    int load_checkpoint();
    /// @}

    /// @name Diagnostic files
    /// @{
    std::ofstream flux_file;
    std::ofstream energy_flux_file;
    std::ofstream energy_file;
    std::ofstream norms_file;
    std::ofstream time_file;
    std::ofstream energy_modes_file;
    std::ofstream R_file;
    std::ofstream theta_file;
    std::ofstream Omega_file;
    std::ofstream p_bin_file;
    std::ofstream p_bin_inst_file;
    std::ofstream d_V_bin_dt_file;
    std::ofstream phase_file;
    std::ofstream R_dot_file;
    std::ofstream p_ramp_inst_file;
    std::ofstream M_1_file;
    std::ofstream q_file;
    /// @}

    /// Probability flux of each mode through the extraction radius.
    std::vector<double> get_flux(const Eigen::VectorXcd &psi);
    /// Energy flux of each mode through the extraction radius. @p dpsi_dt is the
    /// actual derivative @f$-iH(t)\psi(t)@f$, at the same time as @p psi.
    std::vector<double> get_energy_flux(const Eigen::VectorXcd &psi, const Eigen::VectorXcd &dpsi_dt);

    /// Real part of @f$\langle\psi|A|\psi\rangle@f$ restricted to nodes 0..@p i_cut.
    double restricted_expectation(const Eigen::VectorXcd &psi, const SparseComplex &A, int i_cut) const;

    /// @f$\mathrm{Re}\langle\dot\psi|V|\psi\rangle@f$ inside the extraction radius.
    double get_dVexpectation_dt(const Eigen::VectorXcd &psi, const Eigen::VectorXcd &dpsi_dt);

    /// Instantaneous orbital-motion power, excluding companion switching work.
    double get_binary_power(const Eigen::VectorXcd &psi);
    /// Artificial work from switching the companion potential during startup.
    double get_switching_power(const Eigen::VectorXcd &psi);

    /// Norm of each mode inside the extraction radius.
    std::vector<double> get_norms(const Eigen::VectorXcd &psi);
    /// Energy @f$\langle\psi|H|\psi\rangle@f$ inside the extraction radius.
    double get_energy(const Eigen::VectorXcd &psi);
    /// Energy of each mode inside the extraction radius.
    std::vector<double> get_energy_per_mode(const Eigen::VectorXcd &psi);

    /// Append one value as a line of @p quantity_file.
    void write_quantity_line(const double &quantity, std::ofstream &quantity_file);
    /// Append one value per mode as a space-separated line of @p quantity_file.
    void write_quantity_line_modes(const std::vector<double> &quantity, std::ofstream &quantity_file);
};

#endif
