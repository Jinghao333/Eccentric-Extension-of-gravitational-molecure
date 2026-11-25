/**
 * @file main.cpp
 * @brief Command-line driver for binary simulations.
 *
 * Usage:
 * @code
 *   GravitationalMolecule config.toml initial_modes.bin   # new run
 *   GravitationalMolecule Simulation_YYYY_MM_DD_HH_MM_SS   # resume a run
 * @endcode
 */
#include <cmath>
#include <csignal>
#include <iostream>
#include <filesystem>
#include <iomanip>
#include <sstream>
#include <Eigen/Core>
#include <omp.h>
#include "Config.hpp"
#include "CrankNicolsonModes.hpp"
#include "InitialData.hpp"

namespace
{
volatile std::sig_atomic_t stop_requested = 0;

extern "C" void request_stop(int)
{
    stop_requested = 1;
}

void report_configuration(const Config &cfg, const char *exe, const std::string &file)
{
    const FixedKeplerOrbit trajectory(1.0 + cfg.q, cfg.R, cfg.e, cfg.counter);
    std::ostringstream out;
    out << std::setprecision(12) << std::boolalpha
        << "Executable: " << std::filesystem::absolute(exe).string()
        << "\nConfiguration: " << std::filesystem::absolute(file).string()
        << "\nq=" << cfg.q << " e=" << cfg.e << " R_peri=" << cfg.R
        << " P0=" << trajectory.period() << " evolve_orbit=" << cfg.evolve_orbit
        << " counter=" << cfg.counter
        << "\nintegrator=" << (cfg.integrator == FieldIntegrator::MidpointCN ? "midpoint_cn" : "legacy_endpoint")
        << " timestep_rule=";
    if (cfg.geometry.normalized)
        out << "normalized_geometry startup_steps_per_period=" << cfg.geometry.startup_steps
            << " orbit_steps_per_period=" << cfg.geometry.orbit_steps
            << " startup_spacing=" << trajectory.period() / cfg.geometry.startup_steps;
    else if (cfg.geometry.enabled)
        out << "orbital_geometry maximum_dt=" << cfg.dt << " geometry_tolerance=" << cfg.geometry.tolerance
            << " effective_azimuthal_order=" << cfg.geometry.azimuthal_order
            << " effective_radial_order=" << cfg.geometry.radial_order
            << " ramp_weight_tolerance=" << cfg.geometry.ramp_tolerance;
    else if (cfg.evolve_orbit)
        out << "legacy_circular_inspiral initial_dt=" << cfg.dt << " steps_per_period=" << cfg.steps_per_period;
    else out << "fixed dt=" << cfg.dt;
    out << "\nramp_periods=" << cfg.ramp_periods << " ramp_secondary_potential=" << cfg.ramp_secondary_potential
        << " absorption=" << cfg.absorption << " N=" << cfg.N << " L=" << cfg.L
        << " l_max=" << cfg.l_max << " r_ex=" << cfg.r_ex
        << "\nmask=" << cfg.mask;
    if (cfg.mask) out << " L_factor=" << cfg.mask_L_factor << " W=" << cfg.mask_W << " n=" << cfg.mask_n;
    out << " cap=" << cfg.cap;
    if (cfg.cap) out << " L_factor=" << cfg.cap_L_factor << " W_edge=" << cfg.cap_W << " n=" << cfg.cap_n;
    std::cout << out.str() << '\n';
}
}

int main(int argc, char **argv) try
{
    std::signal(SIGINT, request_stop);
    std::signal(SIGTERM, request_stop);

    // A single directory argument resumes the run stored there, from its
    // config.toml and binary checkpoint.
    const bool resume = argc == 2 && std::filesystem::is_directory(argv[1]);
    const bool check = argc == 3 && std::string(argv[1]) == "--check-config";
    if (argc != 3 && !resume)
    {
        std::cerr << "Usage: " << argv[0] << " config.toml initial_modes.bin\n"
                  << "       " << argv[0] << " --check-config config.toml\n"
                  << "       " << argv[0] << " simulation_folder   (resume)\n";
        return 1;
    }
    std::string config_file = "config.toml";
    std::string resume_folder;

    if (resume)
    {
        resume_folder = std::filesystem::path(argv[1]).lexically_normal().string();
        while (resume_folder.size() > 1 && resume_folder.back() == '/')
            resume_folder.pop_back();
        config_file = resume_folder + "/config.toml";
    }
    else
        config_file = argv[check ? 2 : 1];

    Config cfg = load_config(config_file);
    if (cfg.N < 3 || cfg.l_max < 0 || cfg.L <= 0 || cfg.dt <= 0 || cfg.save_every <= 0)
        throw std::invalid_argument("Require N >= 3, l_max >= 0, and positive L, dt, save_every");
    const double dr = cfg.L / (cfg.N - 1);
    if (cfg.r_ex < dr || std::floor(cfg.r_ex / dr) >= cfg.N - 1)
        throw std::invalid_argument("Require an interior extraction radius");
    if (cfg.num_steps < 0 || (cfg.num_steps == 0 && cfg.num_periods <= 0))
        throw std::invalid_argument("Set positive simulation.num_steps or simulation.num_periods");
    report_configuration(cfg, argv[0], config_file);
    if (check) return 0;

    // A positive configuration value overrides OMP_NUM_THREADS; otherwise use
    // the OpenMP runtime's environment/default for both Eigen and our kernels.
    if (cfg.threads > 0)
        omp_set_num_threads(cfg.threads);
    Eigen::setNbThreads(omp_get_max_threads());
    std::cout << "Using up to " << omp_get_max_threads() << " OpenMP threads.\n";

    CrankNicolsonModes solver(
        cfg.N,
        cfg.L,
        cfg.dt,
        cfg.epsilon,
        cfg.l_max,
        cfg.mu,
        cfg.M_c,
        cfg.q,
        cfg.a,
        cfg.R,
        cfg.counter,
        cfg.absorption,
        cfg.save_every,
        cfg.ramp_periods,
        cfg.epsilon_smear,
        config_file,
        cfg.e, cfg.integrator, cfg.geometry);

    if (!resume)
        solver.set_initial_state_full(load_initial_modes_binary(argv[2], solver.modes));

    solver.evolve_orbit = cfg.evolve_orbit;
    solver.ramp_secondary_potential = cfg.ramp_secondary_potential;
    solver.steps_per_period = cfg.steps_per_period;
    solver.set_power_average_orbits(cfg.power_average_orbits);
    std::cout << "Startup ramp over " << cfg.ramp_periods << " initial periods: "
              << (cfg.ramp_secondary_potential ? "frequency and companion potential"
                                               : "frequency only (companion potential at full strength)")
              << "\n";

    solver.use_block_preconditioner = cfg.preconditioner == "radial";
    solver.set_stop_flag(&stop_requested);
    solver.initialize();

    if (cfg.mask)
    {
        solver.build_damping_mask(cfg.mask_L_factor, cfg.mask_W, cfg.mask_n);
    }

    if (cfg.cap)
    {
        solver.build_cap(cfg.cap_L_factor, cfg.cap_W, cfg.cap_n);
    }

    if (cfg.geometry.enabled) {
        const double final_time = cfg.num_periods * solver.initial_orbital_period();
        if (resume) solver.resume_until(resume_folder, final_time, cfg.r_ex);
        else solver.run_until(final_time, cfg.r_ex);
        return 0;
    }

    // num_steps takes precedence; otherwise run for num_periods initial periods.
    const int num_steps = cfg.num_steps
        ? cfg.num_steps
        : static_cast<int>(std::floor(cfg.num_periods * solver.initial_orbital_period() / solver.dt));

    if (resume)
        solver.resume(resume_folder, num_steps, cfg.r_ex);
    else
        solver.run(num_steps, cfg.r_ex);

    return 0;
}
catch (const std::exception &error)
{
    std::cerr << "GravitationalMolecule: " << error.what() << '\n';
    return 1;
}
