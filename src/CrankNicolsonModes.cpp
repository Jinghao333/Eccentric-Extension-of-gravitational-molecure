/**
 * @file CrankNicolsonModes.cpp
 * @brief Construction, initial state, single timestep and run driver.
 */
#include "CrankNicolsonModes.hpp"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <iostream>
#include <vector>

CrankNicolsonModes::CrankNicolsonModes(
    int Nr, double L, double dt_i, double epsilon, int l_max, double mu,
    double M_c, double q, double a, double R, bool counter, bool absorption,
    int save_every, double ramp_periods, double epsilon_smear,
    std::string config_file, double eccentricity, FieldIntegrator integrator, GeometryTimestep geometry)
    : Nr(Nr), save_every(save_every), l_max(l_max), epsilon(epsilon), mu(mu),
      a(a), L(L), dr(L / (Nr - 1)), dt(dt_i), counter(counter), M_c(M_c),
      integrator(integrator), geometry(geometry), maximum_dt(dt_i), ramp_periods(ramp_periods), initial_period(0.0), config_file(config_file),
      absorption(absorption), epsilon_smear(epsilon_smear)
{

    if (Nr < 3 || l_max < 0 || L <= 0.0 || dt <= 0.0 || save_every <= 0)
        throw std::invalid_argument("Require Nr >= 3, l_max >= 0, and positive L, dt, save_every");
    orbit.M1 = 1.0;
    orbit.M2 = q;
    orbit.R = R;
    orbit.theta = 0.0;
    orbit.R_dot = 0.0;
    orbit.theta_dot = 0.0;
    orbit.phase = 0.0;
    fixed_kepler.configure(orbit.total_mass(), R, eccentricity, counter);
    initial_period = fixed_kepler.period();

    if (q > 0.0)
        std::cout << "Orbital period = " << initial_period << std::endl;

    if (evolve_orbit && q > 0.0 && eccentricity == 0.0)
    {
        const double t_c = 5.0 / 256.0 * pow(R, 4) / (q * (1 + q));
        std::cout << "Vacuum coalescence time = " << t_c << std::endl;
    }

    std::cout << "Building grid..." << std::endl;
    build_grid();

    std::cout << "Building modes list for l_max = " << l_max << "..." << std::endl;
    build_modes();

    std::cout << "Building stationary hamiltonian..." << std::endl;
    build_stationary_hamiltonian();
}


void CrankNicolsonModes::set_power_average_orbits(double n_orbits)
{
    power_average = RollingPhaseAverage(n_orbits); // throws unless n_orbits > 0
    orbits_average = n_orbits;
}

void CrankNicolsonModes::configure_isolated_black_hole(double x)
{
    if (!std::isfinite(x) || std::abs(x) >= L)
        throw std::invalid_argument("The isolated BH must lie inside the radial domain");
    if (A.rows() != 0)
        throw std::logic_error("Configure the isolated BH before initialize()");
    isolated_black_hole = true;
    evolve_orbit = false;
    orbit.M1 = 1.0;
    orbit.M2 = 0.0;
    orbit.R = std::abs(x);
    orbit.theta = x < 0.0 ? constants::pi : 0.0;
    orbit.R_dot = 0.0;
    orbit.theta_dot = 0.0;
    orbit.phase = 0.0;
    ramp_periods = 0.0;
}

void CrankNicolsonModes::initialize()
{
    validate_orbit_configuration();

    std::cout << "Building cache of mode couplings..." << std::endl;
    build_coupling_cache();

    std::cout << "Building sparse element caches..." << std::endl;
    build_sparse_element_caches();

    std::cout << "Building initial perturbation..." << std::endl;
    t_i = 0.0;
    update_orbit_and_ramp(0.0, !evolve_orbit);
    build_time_dependent_perturbation(orbit, hamiltonian_evaluation(0.0));

    std::cout << "Building sparse matrices for the Crank-Nicolson Method..." << std::endl;
    build_CN_matrices();

    if (absorption)
    {
        std::cout << "Spherical absorption activated." << std::endl;
    }

    std::cout << "Solver ready to run!" << std::endl;
}

void CrankNicolsonModes::build_grid()
{
    r.resize(Nr);
    for (int i = 0; i < Nr; i++)
        r(i) = epsilon + i * dr;
}

void CrankNicolsonModes::build_modes()
{
    n_modes = 0;
    std::cout << "Modes list: " << std::endl;

    for (int l = 0; l <= l_max; ++l)
    {
        for (int m = -l; m <= l; ++m)
        {
            // Equatorial symmetry: modes with odd l - m are never excited.
            if (abs(l - m) % 2 == 1)
            {
                continue;
            }
            modes.push_back({l, m});
            std::cout << "(" << l << ", " << m << ") ";
            n_modes += 1;
        }
    }

    std::cout << std::endl;
}

void CrankNicolsonModes::set_initial_state_full(Eigen::VectorXcd psi_0)
{
    if (psi_0.size() != Nr * n_modes)
        throw std::invalid_argument("Initial state size does not match the radial grid and modes");
    std::cout << "Assigning initial state..." << std::endl;
    psi.resize(Nr * n_modes);

    for (int i = 0; i < psi.size(); i++)
        psi[i] = psi_0[i];

    psi *= std::sqrt(M_c / mu);

    psi_last = Eigen::VectorXcd::Zero(Nr * n_modes);
}

void CrankNicolsonModes::step(int step_id)
{
    if (geometry.enabled && !(t_i < target_time))
        throw std::logic_error("No remaining geometric step; set a later physical-time target");
    // Keep the previous state for the orbit's midpoint accretion rates.
    psi_last = psi;

    const double step_start = t_i;
    const double step_dt = dt;
    const double endpoint = step_start + step_dt;
    const bool midpoint = integrator == FieldIntegrator::MidpointCN;
    OrbitState propagation_orbit;
    double propagation_time = endpoint;
    if (midpoint)
    {
        const OrbitState start = orbit;
        propagation_time = step_start + 0.5 * step_dt;
        propagation_orbit = evaluate_prescribed_orbit_at(start, step_start, propagation_time);
        orbit = evaluate_prescribed_orbit_at(start, step_start, endpoint);
        update_orbit_and_ramp(endpoint, false); // Commit the endpoint clock factor.
    }
    t_i = endpoint;
    if (!midpoint)
        advance_orbit(step_dt, t_i);

    report_step_progress(step_id);

    t(step_id) = t_i;

    build_time_dependent_perturbation(midpoint ? propagation_orbit : orbit,
                                      hamiltonian_evaluation(propagation_time));
    propagate_field(step_dt);
    // A/B retain the propagation operator. All field measurements below use
    // the endpoint H, potential derivatives and horizon projectors.
    if (midpoint)
        build_time_dependent_perturbation(orbit, hamiltonian_evaluation(t_i));

    // The flux and snapshots use the post-mask endpoint state. Evaluate the
    // Schrodinger RHS on that same state, with the completed step's H(t_i).
    // The mask itself is a discrete damping operation, not part of this RHS.
    const Eigen::VectorXcd dpsi_dt = schrodinger_rhs(psi);

    apply_orbital_feedback(step_start);
    dt = select_next_timestep();

    // Keep the completed step's pre-feedback Hamiltonian for diagnostics.
    const auto observables = compute_step_observables(dpsi_dt);
    update_rolling_power(observables.binary_power);
    if (geometry.enabled)
        dt_step_file << std::setprecision(17) << step_dt << '\n';
    write_step_output(step_id, observables);
}

void CrankNicolsonModes::set_time_target(double final_time)
{
    if (!geometry.enabled || !std::isfinite(final_time) || final_time <= 0)
        throw std::invalid_argument("Physical-time driver requires orbital_geometry and a finite positive target");
    validate_orbit_configuration();
    target_time = final_time;
}

void CrankNicolsonModes::run_until(double final_time, double extraction)
{
    set_time_target(final_time);
    run_impl(1, extraction);
}

void CrankNicolsonModes::run(int n_steps, double extraction)
{
    if (geometry.enabled) throw std::invalid_argument("Use run_until for orbital_geometry");
    run_impl(n_steps, extraction);
}

void CrankNicolsonModes::run_impl(int n_steps, double r_ex_)
{

    std::cout << "Running simulation..." << std::endl;

    set_extraction_radius(r_ex_, n_steps);
    t_i = 0.0;
    t.resize(n_steps);
    t(0) = t_i;

    initial_physics = describe_initial_physics();
    output_folder = create_output_folder();
    std::cout << "Output folder: " << output_folder << std::endl;

    open_diagnostics(std::ios::out);

    save_run_info();
    // Refresh the initial operator, including a CAP enabled after initialize().
    update_orbit_and_ramp(t_i, !evolve_orbit);
    build_time_dependent_perturbation(orbit, hamiltonian_evaluation(t_i));
    save_snapshot(0, psi, "psi");
    save_snapshot(0, schrodinger_rhs(psi), "psi_dt");
    if (geometry.enabled) dt = select_next_timestep();
    save_checkpoint(0);

    evolve(1, n_steps);
}

void CrankNicolsonModes::set_extraction_radius(double r_ex_, int n_steps)
{
    r_ex = r_ex_;
    idx_ex = (int)std::floor(r_ex / dr);

    std::cout << "Computing fluxes at: " << r_ex_ << " [Index: " << idx_ex << "]" << std::endl;

    if (n_steps < 1 || idx_ex < 1 || idx_ex >= Nr - 1)
        throw std::invalid_argument("Require at least one time sample and an interior extraction radius");
}

void CrankNicolsonModes::evolve(int first_step, int n_steps)
{
    int last_completed = first_step - 1;
    bool stopped = false;

    for (int n = first_step; geometry.enabled ? t_i < target_time : n < n_steps; ++n)
    {
        if (stop_flag && *stop_flag)
        {
            stopped = true;
            break;
        }
        if (!isolated_black_hole && evolve_orbit)
        {
            if (orbit.R <= 6.0)
            {
                break;
            }
        }
        if (geometry.enabled) {
            if (n >= std::numeric_limits<int>::max()-1)
                throw std::overflow_error("Too many geometric steps for checkpoint step index");
            if (n >= t.size())
                t.conservativeResize(std::clamp<long long>(2LL*t.size(),64,
                                                           std::numeric_limits<int>::max()));
        }
        step(n);
        last_completed = n;
        if (stop_flag && *stop_flag)
        {
            stopped = true;
            break;
        }
    }

    std::cout << std::endl;
    if (geometry.enabled) t.conservativeResize(last_completed+1);
    flush_diagnostics();
    save_checkpoint(last_completed);
    if (stopped)
        std::cout << "Stop requested; checkpoint saved after step " << last_completed << "." << std::endl;
    else
        std::cout << "Simulation completed; checkpoint saved after step " << last_completed << "." << std::endl;
}
