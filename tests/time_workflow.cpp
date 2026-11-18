// Deterministic workflow/restart fixture. --record DIR captures reference data;
// --compare DIR requires byte identity with a capture from the original solver.
#include "CrankNicolsonModes.hpp"
#include "TestUtils.hpp"
#include <omp.h>
#include <algorithm>
#include <memory>

namespace fs = std::filesystem;
using test::require;
namespace {
constexpr int last_step = 24;

std::string bytes(const fs::path &path)
{
    std::ifstream in(path, std::ios::binary);
    require(bool(in), "Cannot read " + path.string());
    return {std::istreambuf_iterator<char>(in), {}};
}

void compare_tree(const fs::path &reference, const fs::path &actual)
{
    std::size_t count = 0, actual_count = 0;
    for (const auto &entry : fs::recursive_directory_iterator(reference)) {
        if (!entry.is_regular_file()) continue;
        const auto relative = fs::relative(entry.path(), reference);
        require(bytes(entry.path()) == bytes(actual / relative),
                "Byte mismatch: " + relative.string());
        ++count;
    }
    for (const auto &entry : fs::recursive_directory_iterator(actual))
        actual_count += entry.is_regular_file();
    require(count == actual_count && count > 20, "Output file set changed");
}

template<class T> void raw(std::ostream &out, const T &value)
{
    out.write(reinterpret_cast<const char *>(&value), sizeof(value));
}

void capture(std::ostream &out, const CrankNicolsonModes &cn, int step)
{
    for (double value : {cn.t[step], cn.dt, cn.orbit.M1, cn.orbit.M2,
                        cn.orbit.R, cn.orbit.theta, cn.orbit.R_dot,
                        cn.orbit.theta_dot, cn.orbit.phase,
                        cn.sink_strength(1), cn.sink_strength(2)}) raw(out, value);
    auto vector = [&](const Eigen::VectorXcd &v) {
        for (Eigen::Index i = 0; i < v.size(); ++i) {
            raw(out, v[i].real()); raw(out, v[i].imag());
        }
    };
    vector(cn.psi);
    vector(cn.schrodinger_rhs(cn.psi));
    for (const auto *matrix : {&cn.H, &cn.V1, &cn.V2, &cn.V,
                              &cn.dVdt, &cn.dVdt_switch, &cn.A, &cn.B}) {
        raw(out, static_cast<std::int64_t>(matrix->nonZeros()));
        for (int row = 0; row < matrix->rows(); ++row)
            for (SparseComplex::InnerIterator it(*matrix, row); it; ++it) {
                raw(out, static_cast<std::int64_t>(it.row()));
                raw(out, static_cast<std::int64_t>(it.col()));
                raw(out, it.value().real()); raw(out, it.value().imag());
            }
    }
    // Full sink matrix at representative endpoints, probes at every other step.
    Eigen::VectorXcd probe = Eigen::VectorXcd::Zero(cn.psi.size());
    if (step == 0 || step == 3 || step == 8 || step == last_step) {
        for (int j = 0; j < probe.size(); ++j) {
            probe[j] = 1.0;
            vector(cn.apply_sink(probe));
            probe[j] = 0.0;
        }
    } else {
        for (int j = 0; j < probe.size(); ++j)
            probe[j] = std::polar(1.0 / (1 + j), 0.13 * j);
        vector(cn.apply_sink(probe));
    }
    require(bool(out), "Trace write failed");
}

std::unique_ptr<CrankNicolsonModes> make_solver(int scenario)
{
    const bool eccentric = scenario == 2 || scenario == 3;
    auto cn = std::make_unique<CrankNicolsonModes>(33, 64.0, 2.0, 0.1, 2,
        0.2, 0.001, 0.25, 0.0, 16.0, scenario % 2 != 0, scenario != 0,
        1, scenario == 0 || scenario == 5 ? 0.0 :
           (scenario == 3 ? 2.5 / 16 : 2.0 / 16), 4.0, "fixture.toml",
        eccentric ? 0.4 : 0.0);
    cn->evolve_orbit = scenario == 4 || scenario == 5;
    cn->ramp_secondary_potential = scenario != 3;
    cn->steps_per_period = 24;
    cn->set_power_average_orbits(0.5);
    if (scenario == 6) cn->configure_isolated_black_hole(-8.0);
    cn->dt = cn->initial_orbital_period() / 16;
    using Mode = CrankNicolsonModes::OrbitMode;
    using Rule = CrankNicolsonModes::TimeStepRule;
    const Mode expected = scenario == 6 ? Mode::Isolated :
        cn->evolve_orbit ? Mode::CircularInspiralLegacy :
        eccentric ? Mode::FixedEccentric : Mode::FixedCircular;
    require(cn->orbit_mode() == expected, "Incorrect derived orbit mode");
    require(cn->timestep_rule() == (cn->evolve_orbit ? Rule::LegacyCircularInspiral : Rule::FixedDt),
            "Incorrect derived timestep rule");
    if (cn->evolve_orbit) {
        cn->evolve_orbit = false;
        require(cn->orbit_mode() == Mode::FixedCircular && cn->timestep_rule() == Rule::FixedDt,
                "Modes are stale copies of configuration");
        cn->evolve_orbit = true;
    }
    return cn;
}

fs::path run_case(int scenario, int restart_step, const fs::path &root)
{
    fs::create_directories(root);
    const auto previous = fs::current_path();
    fs::current_path(root);
    std::ofstream("fixture.toml") << "# time workflow fixture\n";
    auto cn = make_solver(scenario);
    Eigen::VectorXcd initial(cn->Nr * cn->n_modes);
    for (int j = 0; j < initial.size(); ++j)
        initial[j] = std::polar(std::sin(constants::pi * (j % cn->Nr) /
                                      (cn->Nr - 1)) / (1 + j / cn->Nr), 0.13 * j);
    initial /= std::sqrt(initial.squaredNorm() * cn->dr);
    cn->set_initial_state_full(initial);
    cn->initialize();
    cn->build_cap(0.8, 0.002, 3);
    if (scenario == 3 || scenario == 6) cn->build_damping_mask(0.8, 0.01, 2);
    cn->run(1, 40.0);
    cn->t.conservativeResize(last_step + 1);
    const auto folder = test::simulation_folder(root);
    std::ofstream trace(folder / "workflow.bin", std::ios::binary);
    capture(trace, *cn, 0);
    for (int step = 1; step <= last_step; ++step) {
        const double step_dt = cn->dt;
        cn->step(step);
        // The stored matrices describe the completed solve, even when feedback
        // has already selected a different cn.dt for the next step.
        const dcomplex scale = I * (step_dt * 0.5);
        for (int row = 0; row < cn->H.rows(); ++row)
            for (SparseComplex::InnerIterator it(cn->H, row); it; ++it) {
                const double diagonal = it.row() == it.col() ? 1.0 : 0.0;
                require(cn->A.coeff(it.row(), it.col()) == diagonal + scale * it.value() &&
                        cn->B.coeff(it.row(), it.col()) == diagonal - scale * it.value(),
                        "CN used the next timestep instead of the captured duration");
            }
        capture(trace, *cn, step);
        if (step == restart_step) {
            cn.reset();
            cn = make_solver(scenario);
            cn->initialize();
            cn->build_cap(0.8, 0.002, 3);
            if (scenario == 3 || scenario == 6) cn->build_damping_mask(0.8, 0.01, 2);
            cn->resume(folder.string(), step + 1, 40.0);
            cn->t.conservativeResize(last_step + 1);
        }
    }
    cn.reset();
    trace.close();
    fs::current_path(previous);
    return folder;
}
}

int main(int argc, char **argv)
{
    try {
        omp_set_dynamic(0);
        omp_set_num_threads(1);
        Eigen::setNbThreads(1);
        require(argc == 1 || (argc == 3 &&
            (std::string(argv[1]) == "--record" || std::string(argv[1]) == "--compare")),
            "Usage: time_workflow [--record|--compare DIRECTORY]");
        const auto reference = argc == 3 ? fs::absolute(argv[2]) : fs::path{};
        if (argc == 3 && std::string(argv[1]) == "--record")
            require(!fs::exists(reference), "Refusing to overwrite baseline");
        test::ScratchDir scratch("workflow");
        for (int scenario = 0; scenario < 7; ++scenario) {
            const auto full = run_case(scenario, -1, scratch.path() / ("full" + std::to_string(scenario)));
            for (int split : {1, 3, 8}) {
                const auto resumed = run_case(scenario, split, scratch.path() /
                    ("resume" + std::to_string(scenario) + "_" + std::to_string(split)));
                compare_tree(full, resumed);
            }
            if (argc == 3) {
                const auto target = reference / std::to_string(scenario);
                if (std::string(argv[1]) == "--record") {
                    fs::create_directories(target);
                    fs::copy(full, target, fs::copy_options::recursive);
                } else compare_tree(target, full);
            }
        }
        std::cout << "\nSeven scenarios: exact state/operator/output/restart comparisons passed.\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
