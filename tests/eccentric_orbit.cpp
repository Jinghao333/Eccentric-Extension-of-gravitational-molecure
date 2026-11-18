#include "CrankNicolsonModes.hpp"
#include "FixedKeplerOrbit.hpp"
#include "SphericalAbsorptionReference.hpp"
#include "TestUtils.hpp"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <cstdint>
#include <cstring>

namespace fs = std::filesystem;
using test::close;
using test::require;

namespace {

std::unique_ptr<CrankNicolsonModes> make_solver(double e, double ramp_periods,
                                                bool absorb = false,
                                                bool counter = false,
                                                int save_every = 5)
{
    std::ofstream("config.toml") << "# eccentric test fixture\n";
    auto solver = std::make_unique<CrankNicolsonModes>(
        49, 72.0, 0.4, 0.13, 3, 0.2, 0.02, 0.3, 0.0, 18.0,
        counter, absorb, save_every, ramp_periods, 5.0, "config.toml", e);
    solver->evolve_orbit = false;
    solver->set_power_average_orbits(0.5);
    return solver;
}

Eigen::VectorXcd initial_state(const CrankNicolsonModes &solver)
{
    Eigen::VectorXcd state = Eigen::VectorXcd::Zero(solver.Nr * solver.n_modes);
    for (int j = 0; j < solver.n_modes; ++j)
        for (int i = 1; i < solver.Nr - 1; ++i)
            state[j * solver.Nr + i] =
                std::polar(std::exp(-std::pow((solver.r[i] - 10.0 - j) / 5.0, 2)),
                           0.17 * j);
    return state / std::sqrt(state.squaredNorm() * solver.dr);
}

struct Operators
{
    Eigen::MatrixXcd V, dVdt, switching;
    OrbitState orbit;
    std::vector<mode> modes;
    int Nr;
    Eigen::Index V_nonzeros, dVdt_nonzeros;
};

Operators operators_at(double mean_anomaly, double eccentricity, bool counter = false)
{
    auto solver = make_solver(eccentricity, 0.0, false, counter);
    solver->orbit.phase = mean_anomaly;
    solver->set_initial_state_full(Eigen::VectorXcd::Zero(solver->Nr * solver->n_modes));
    solver->initialize();
    return {Eigen::MatrixXcd(solver->V), Eigen::MatrixXcd(solver->dVdt),
            Eigen::MatrixXcd(solver->dVdt_switch), solver->orbit,
            solver->modes, solver->Nr, solver->V.nonZeros(), solver->dVdt.nonZeros()};
}

void eccentric_orbital_clock_ramp()
{
    test::ScratchDir dir("eccentric_clock");
    auto solver = make_solver(0.3, 1.0);
    solver->set_initial_state_full(initial_state(*solver));
    solver->initialize();
    const double period = solver->initial_orbital_period();
    solver->dt = 0.25 * period;
    solver->run(1, 36.0);
    solver->t.resize(4);

    solver->step(1);
    close(solver->orbit.phase, 0.0, 0.0,
          "First ramp step did not use the old zero clock rate");
    close(solver->orbit.R, 18.0, 2e-14, "Orbit left pericenter on the first ramp step");
    const double g_quarter = 3.0 * 0.25 * 0.25 - 2.0 * 0.25 * 0.25 * 0.25;
    close(solver->orbit.theta_dot,
          g_quarter * solver->mean_motion() * std::sqrt(1.0 - 0.3 * 0.3) /
              std::pow(1.0 - 0.3, 2),
          3e-14, "Endpoint angular rate does not use the ramped Kepler clock");

    solver->step(2);
    const double expected_phase = 2.0 * constants::pi * 0.25 * g_quarter;
    close(solver->orbit.phase, expected_phase, 3e-14,
          "Eccentric mean anomaly did not advance with the previous ramp rate");
    FixedKeplerOrbit trajectory(1.3, 18.0, 0.3, false);
    OrbitState expected = solver->orbit;
    expected.phase = expected_phase;
    const double g_half = 0.5;
    trajectory.evaluate(expected, g_half);
    close(solver->orbit.R, expected.R, 3e-14, "Ramped eccentric radius is wrong");
    close(solver->orbit.theta, expected.theta, 3e-14, "Ramped eccentric angle is wrong");
    close(solver->orbit.R_dot, expected.R_dot, 3e-14, "Ramped radial velocity is wrong");
    close(solver->orbit.theta_dot, expected.theta_dot, 3e-14,
          "Ramped angular velocity is wrong");
}

void orbital_derivative_matches_finite_difference()
{
    constexpr double e = 0.3, mean_anomaly = 1.17, delta_t = 2e-5;
    const Operators center = operators_at(mean_anomaly, e);
    FixedKeplerOrbit trajectory(1.3, 18.0, e, false);
    const double delta_M = trajectory.mean_motion() * delta_t;
    const Operators minus = operators_at(mean_anomaly - delta_M, e);
    const Operators plus = operators_at(mean_anomaly + delta_M, e);
    const Eigen::MatrixXcd finite_difference = (plus.V - minus.V) / (2.0 * delta_t);
    const double scale = std::max(1.0, center.dVdt.norm());
    const double derivative_error = (finite_difference - center.dVdt).norm() / scale;
    require(derivative_error < 2e-8,
            "Eccentric orbital dVdt disagrees with a centered finite difference");
    require((center.V - center.V.adjoint()).norm() < 2e-13 * std::max(1.0, center.V.norm()),
            "Eccentric gravitational potential is not Hermitian");
    require((center.dVdt - center.dVdt.adjoint()).norm() <
                2e-13 * std::max(1.0, center.dVdt.norm()),
            "Eccentric orbital derivative is not Hermitian");

    Eigen::MatrixXcd angular = Eigen::MatrixXcd::Zero(center.V.rows(), center.V.cols());
    for (int row = 0; row < center.V.rows(); ++row)
        for (int col = 0; col < center.V.cols(); ++col)
        {
            const int m_row = center.modes[row / center.Nr].m;
            const int m_col = center.modes[col / center.Nr].m;
            angular(row, col) = I * static_cast<double>(m_col - m_row) *
                                center.orbit.theta_dot * center.V(row, col);
        }
    require(center.orbit.R_dot != 0.0 && (center.dVdt - angular).norm() > 1e-8,
            "Eccentric dVdt has no radial contribution");
    require(center.dVdt_nonzeros == center.V_nonzeros,
            "Orbital derivative changed the potential sparsity");
    std::cout << "  relative orbital dVdt error = " << derivative_error << '\n';
}

Eigen::MatrixXcd switching_at(double time_fraction, bool derivative)
{
    test::ScratchDir dir("switching_operator");
    auto solver = make_solver(0.3, 1.0);
    solver->set_initial_state_full(initial_state(*solver));
    solver->initialize();
    solver->run(1, 36.0);
    solver->t.resize(2);
    solver->dt = time_fraction * solver->initial_orbital_period();
    solver->step(1); // old ramp rate is zero, so geometry remains at pericenter
    return derivative ? Eigen::MatrixXcd(solver->dVdt_switch)
                      : Eigen::MatrixXcd(solver->V);
}

void switching_derivative_is_separate()
{
    constexpr double x = 0.37, delta_x = 1e-5;
    const Eigen::MatrixXcd center = switching_at(x, true);
    const double period = FixedKeplerOrbit(1.3, 18.0, 0.3, false).period();
    const Eigen::MatrixXcd finite_difference =
        (switching_at(x + delta_x, false) - switching_at(x - delta_x, false)) /
        (2.0 * delta_x * period);
    const double switching_error = (finite_difference - center).norm() /
                                   std::max(1.0, center.norm());
    require(switching_error < 2e-9,
            "Companion switching derivative disagrees with finite difference");

    const Operators no_ramp = operators_at(1.0, 0.3);
    require(no_ramp.switching.norm() == 0.0,
            "Switching work leaked into an un-ramped orbit");
    std::cout << "  relative switching dVdt error = " << switching_error << '\n';
}

void absorption_follows_moving_sources()
{
    test::ScratchDir dir("eccentric_absorption_geometry");
    auto solver = make_solver(0.3, 0.0, true, false);
    const Eigen::VectorXcd probe = initial_state(*solver);
    solver->set_initial_state_full(probe);
    solver->initialize();
    solver->run(1, 36.0);
    solver->t.resize(3);
    solver->dt = 0.2 * solver->initial_orbital_period();

    solver->step(1);
    const double ramp1 = 3.0 * 0.2 * 0.2 - 2.0 * 0.2 * 0.2 * 0.2;
    const Eigen::VectorXcd sink1 = solver->apply_sink(probe);
    const Eigen::VectorXcd reference1 = sink_reference(*solver, 5.0, ramp1) * probe;
    const double error1 = (sink1 - reference1).norm() / std::max(1.0, reference1.norm());
    require(error1 < 5e-10,
            "Eccentric absorber does not follow its first endpoint positions");
    require(probe.dot(sink1).real() >= -1e-13,
            "Moving eccentric sink is not positive semidefinite");
    const double radius1 = solver->orbit.R;

    solver->step(2);
    const double ramp2 = 3.0 * 0.4 * 0.4 - 2.0 * 0.4 * 0.4 * 0.4;
    const Eigen::VectorXcd sink2 = solver->apply_sink(probe);
    const Eigen::VectorXcd reference2 = sink_reference(*solver, 5.0, ramp2) * probe;
    const double error2 = (sink2 - reference2).norm() / std::max(1.0, reference2.norm());
    require(error2 < 5e-10,
            "Eccentric absorber does not follow its later endpoint positions");
    require(solver->orbit.R != radius1 && (sink2 - sink1).norm() > 1e-8,
            "Absorption profile remained fixed while the eccentric sources moved");
    std::cout << "  moving-sink reference errors = " << error1 << ", " << error2 << '\n';
}

fs::path run_eccentric(const fs::path &root, int samples)
{
    fs::create_directories(root);
    const fs::path previous = fs::current_path();
    fs::current_path(root);
    auto solver = make_solver(0.4, 0.25, true, true);
    solver->set_initial_state_full(initial_state(*solver));
    solver->initialize();
    solver->run(samples, 36.0);
    fs::current_path(previous);
    return test::simulation_folder(root);
}

std::string contents(const fs::path &path)
{
    std::ifstream in(path, std::ios::binary);
    return std::string(std::istreambuf_iterator<char>(in), {});
}

void convert_checkpoint_to_legacy_v1(const fs::path &path)
{
    std::ifstream in(path, std::ios::binary);
    char magic[8];
    std::int32_t header[3];
    double v2[10];
    in.read(magic, sizeof(magic));
    in.read(reinterpret_cast<char *>(header), sizeof(header));
    in.read(reinterpret_cast<char *>(v2), sizeof(v2));
    require(in && std::memcmp(magic, "GMCKPT02", 8) == 0,
            "Eccentric fixture is not a version-2 checkpoint");
    const std::string tail{std::istreambuf_iterator<char>(in), {}};

    const double total = v2[6] + v2[7];
    const double omega = std::abs(v2[5]);
    const double v1[13] = {
        v2[0], v2[1], v2[2], v2[3], v2[6], v2[7], v2[7] / v2[6],
        v2[7] / total * v2[3], v2[6] / total * v2[3], omega,
        2.0 * constants::pi / omega, v2[8], v2[9]};

    std::ofstream out(path, std::ios::binary | std::ios::trunc);
    out.write("GMCKPT01", 8);
    out.write(reinterpret_cast<const char *>(header), sizeof(header));
    out.write(reinterpret_cast<const char *>(v1), sizeof(v1));
    out.write(tail.data(), static_cast<std::streamsize>(tail.size()));
    require(bool(out), "Could not create legacy eccentric checkpoint fixture");
}

void eccentric_evolution_and_restart()
{
    test::ScratchDir dir("eccentric_restart");
    // Continue long enough after the checkpoint for any restart-generated
    // high-frequency packet to cross the small test grid and reach r_ex.
    constexpr int final_step = 121;
    const fs::path reference = run_eccentric(fs::current_path() / "reference", final_step);
    const fs::path interrupted = run_eccentric(fs::current_path() / "interrupted", 32);

    const fs::path previous = fs::current_path();
    fs::current_path(interrupted.parent_path());
    auto resumed = make_solver(0.4, 0.25, true, true);
    resumed->initialize();
    resumed->resume(interrupted.string(), final_step, 36.0);
    resumed.reset(); // close diagnostic streams before bytewise comparison on Windows
    fs::current_path(previous);

    // Restart must preserve not only the orbit, but also the scalar state and
    // every flux used by the ionization dashboard.  A rounded text-snapshot
    // fallback can leave the orbit correct while launching a delayed wave in
    // the dominant cloud mode, so compare the field diagnostics byte-for-byte.
    for (const char *name : {"t.txt", "R_orb.txt", "R_dot_orb.txt", "theta_orb.txt",
                             "phase_orb.txt", "Omega_orb.txt", "M_1.txt", "q.txt",
                             "energy.txt", "p_bin.txt", "p_bin_inst.txt",
                             "p_ramp_inst.txt", "dV_bin_dt.txt", "norms.txt",
                             "flux.txt", "energy_flux.txt", "energy_modes.txt",
                             "checkpoint.bin"})
        require(contents(reference / name) == contents(interrupted / name),
                std::string(name) + " differs after eccentric checkpoint restart");

    for (int step = 0; step < final_step; step += 5)
    {
        std::ostringstream name;
        name << "snapshot_" << std::setw(6) << std::setfill('0') << step << ".txt";
        require(contents(reference / "psi" / name.str()) ==
                    contents(interrupted / "psi" / name.str()),
                name.str() + " differs after eccentric checkpoint restart");
    }

    const auto radii = test::read_rows(reference / "R_orb.txt");
    const auto phases = test::read_rows(reference / "phase_orb.txt");
    const auto masses = test::read_rows(reference / "M_1.txt");
    double min_radius = radii.front()[0], max_radius = radii.front()[0];
    for (std::size_t i = 0; i < radii.size(); ++i)
    {
        require(std::isfinite(radii[i][0]) && std::isfinite(phases[i][0]),
                "Eccentric integration produced nonfinite diagnostics");
        min_radius = std::min(min_radius, radii[i][0]);
        max_radius = std::max(max_radius, radii[i][0]);
        close(masses[i][0], 1.0, 0.0, "Absorption fed back onto a fixed eccentric orbit");
        if (i > 0)
            require(phases[i][0] > phases[i - 1][0],
                    "Eccentric mean-anomaly diagnostic is not increasing");
    }
    require(max_radius > min_radius, "Fixed eccentric bodies did not move around their ellipse");
}

void legacy_eccentric_checkpoint_is_rejected()
{
    test::ScratchDir dir("eccentric_legacy_restart");
    const fs::path interrupted = run_eccentric(fs::current_path() / "interrupted", 12);
    convert_checkpoint_to_legacy_v1(interrupted / "checkpoint.bin");

    const fs::path previous = fs::current_path();
    fs::current_path(interrupted.parent_path());
    bool threw = false;
    try
    {
        auto resumed = make_solver(0.4, 0.25, true, true);
        resumed->initialize();
        resumed->resume(interrupted.string(), 21, 36.0);
    }
    catch (const std::runtime_error &)
    {
        threw = true;
    }
    fs::current_path(previous);
    require(threw, "Legacy eccentric checkpoint was accepted without complete orbit state");
    require(test::read_rows(interrupted / "t.txt").size() == 11,
            "Rejected legacy eccentric restart modified diagnostics");
}

} // namespace

int main()
{
    return test::run_all(
        std::pair{"eccentric orbital-clock ramp", &eccentric_orbital_clock_ramp},
        std::pair{"orbital derivative finite difference", &orbital_derivative_matches_finite_difference},
        std::pair{"separate switching derivative", &switching_derivative_is_separate},
        std::pair{"absorption follows moving sources", &absorption_follows_moving_sources},
        std::pair{"eccentric evolution and restart", &eccentric_evolution_and_restart},
        std::pair{"legacy eccentric checkpoint is rejected",
                  &legacy_eccentric_checkpoint_is_rejected});
}
