// Resuming runs: exact continuation from checkpoint.bin, cleanup of output
// beyond it, legacy circular compatibility, and run extension.
#include "CrankNicolsonModes.hpp"
#include "TestUtils.hpp"
#include <omp.h>
#include <filesystem>
#include <fstream>
#include <cstdint>
#include <cstring>
#include <memory>

namespace fs = std::filesystem;
using test::require;

namespace {

constexpr double r_ex = 40.0;

std::unique_ptr<CrankNicolsonModes> make_solver()
{
    std::ofstream("config.toml") << "# restart fixture\n";
    auto cn = std::make_unique<CrankNicolsonModes>(65, 64.0, 2.0, 0.0, 2, 0.2, 0.05, 0.5, 0.0, 20.0, false, true,
                                                   10, 0.5, 4.0, "config.toml");
    cn->evolve_orbit = true;
    cn->steps_per_period = 32;
    cn->set_power_average_orbits(1.0);
    return cn;
}

Eigen::VectorXcd initial_state(const CrankNicolsonModes &cn)
{
    Eigen::VectorXcd psi = Eigen::VectorXcd::Zero(cn.Nr * cn.n_modes);
    for (int j = 0; j < cn.n_modes; ++j)
        for (int i = 1; i < cn.Nr - 1; ++i)
            psi[j * cn.Nr + i] = std::polar(std::exp(-std::pow((cn.r[i] - 12.0 - 3 * j) / 6.0, 2)), 0.3 * j);
    return psi / std::sqrt(psi.squaredNorm() * cn.dr);
}

/// Start a run of @p n_steps samples in a fresh directory and return its folder.
fs::path new_run(const std::string &name, int n_steps)
{
    const fs::path root = fs::current_path() / name;
    fs::create_directories(root);
    const fs::path previous = fs::current_path();
    fs::current_path(root);
    auto cn = make_solver();
    cn->set_initial_state_full(initial_state(*cn));
    cn->initialize();
    cn->run(n_steps, r_ex);
    fs::current_path(previous);
    return test::simulation_folder(root);
}

/// Start a run with an already-raised signal flag; no new timestep may begin.
fs::path stopped_run(const std::string &name, int n_steps)
{
    const fs::path root = fs::current_path() / name;
    fs::create_directories(root);
    const fs::path previous = fs::current_path();
    fs::current_path(root);
    auto cn = make_solver();
    cn->set_initial_state_full(initial_state(*cn));
    cn->initialize();
    volatile std::sig_atomic_t stop = 1;
    cn->set_stop_flag(&stop);
    cn->run(n_steps, r_ex);
    fs::current_path(previous);
    return test::simulation_folder(root);
}

void resume_run(const fs::path &folder, int n_steps)
{
    auto cn = make_solver();
    cn->initialize();
    cn->resume(folder.string(), n_steps, r_ex);
}

fs::path copy_run(const fs::path &folder, const std::string &name)
{
    const fs::path target = fs::current_path() / name / folder.filename();
    fs::create_directories(target.parent_path());
    fs::copy(folder, target, fs::copy_options::recursive);
    return target;
}

std::string contents(const fs::path &path)
{
    std::ifstream in(path, std::ios::binary);
    return std::string(std::istreambuf_iterator<char>(in), {});
}

/// Every file of @p expected exists in @p actual with identical bytes.
void require_identical(const fs::path &expected, const fs::path &actual)
{
    int compared = 0;
    for (const auto &entry : fs::recursive_directory_iterator(expected)) {
        if (!entry.is_regular_file()) continue;
        const fs::path relative = fs::relative(entry.path(), expected);
        require(fs::exists(actual / relative), "Missing " + relative.string() + " after resume");
        require(contents(entry.path()) == contents(actual / relative),
                relative.string() + " differs from the uninterrupted run");
        ++compared;
    }
    require(compared > 20, "Too few files compared");
}

const char *diagnostics[] = {"t.txt", "R_orb.txt", "theta_orb.txt", "Omega_orb.txt", "M_1.txt", "q.txt",
                             "energy.txt", "p_bin.txt", "norms.txt", "flux.txt", "energy_flux.txt"};

int checkpoint_step(const fs::path &path)
{
    std::ifstream in(path, std::ios::binary);
    char magic[8]{};
    std::int32_t header[3] = {-1, -1, -1};
    in.read(magic, sizeof(magic));
    in.read(reinterpret_cast<char *>(header), sizeof(header));
    require(in && (std::memcmp(magic, "GMCKPT01", 8) == 0 ||
                   std::memcmp(magic, "GMCKPT02", 8) == 0),
            "Invalid checkpoint header");
    return header[0];
}

void convert_checkpoint_to_v1(const fs::path &path)
{
    std::ifstream in(path, std::ios::binary);
    char magic[8];
    std::int32_t header[3];
    double v2[10];
    in.read(magic, sizeof(magic));
    in.read(reinterpret_cast<char *>(header), sizeof(header));
    in.read(reinterpret_cast<char *>(v2), sizeof(v2));
    require(in && std::memcmp(magic, "GMCKPT02", 8) == 0,
            "Fixture is not a version-2 checkpoint");
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
    require(bool(out), "Could not write version-1 checkpoint fixture");
}

void changed_physics_rejected()
{
    test::ScratchDir dir("changed_physics");
    const auto run=new_run("run",25), original=copy_run(run,"original");
    for (int k=0;k<5;++k) {
        auto cn=make_solver();
        if (k==0) cn->mu*=1.1;
        if (k==1) cn->counter=true;
        if (k==2) cn->dt*=0.5;
        if (k==3) cn->steps_per_period++;
        cn->initialize();
        if (k==4) cn->build_cap(0.7,1,2);
        try { cn->resume(run.string(),41,r_ex); }
        catch (const std::runtime_error &error) {
            require(std::string(error.what()).find("physics")!=std::string::npos,"Wrong resume error");
            require_identical(original,run);
            continue;
        }
        throw std::runtime_error("Changed physics accepted");
    }
}

void exact_resume_from_checkpoint()
{
    test::ScratchDir dir("exact_resume");
    const fs::path reference = new_run("reference", 41);
    // An orderly stop after step 24 is between regular snapshots 20 and 30.
    // Restart must use checkpoint 24 directly and continue with step 25.
    const fs::path interrupted = new_run("interrupted", 25);
    require(test::read_rows(interrupted / "t.txt").size() == 24, "Interrupted run has the wrong row count");
    require(checkpoint_step(interrupted / "checkpoint.bin") == 24,
            "Orderly stop did not checkpoint its final completed step");
    require(!fs::exists(interrupted / "psi" / "snapshot_000024.txt"),
            "The test interruption unexpectedly coincides with a text snapshot");
    resume_run(interrupted, 41);
    require_identical(reference, interrupted);
}

void signal_stop_and_resume()
{
    test::ScratchDir dir("signal_stop");
    const fs::path reference = new_run("reference", 41);
    const fs::path interrupted = stopped_run("interrupted", 41);
    require(checkpoint_step(interrupted / "checkpoint.bin") == 0,
            "Raised stop flag did not preserve the last completed step");
    require(test::read_rows(interrupted / "t.txt").empty(),
            "Solver began another timestep after the stop request");
    resume_run(interrupted, 41);
    require_identical(reference, interrupted);
}

void resume_after_crash_mid_write()
{
    test::ScratchDir dir("crash");
    const fs::path reference = new_run("reference", 41);
    const fs::path crashed = new_run("crashed", 25);
    // A crash can leave a partial row and a snapshot written ahead of the rows.
    std::ofstream(crashed / "flux.txt", std::ios::app) << "0.123 0.4";
    fs::copy_file(crashed / "psi" / "snapshot_000020.txt", crashed / "psi" / "snapshot_000030.txt");
    resume_run(crashed, 41);
    require_identical(reference, crashed);
}

void rejects_missing_checkpoint()
{
    test::ScratchDir dir("missing_checkpoint");
    const fs::path interrupted = new_run("interrupted", 25);
    fs::remove(interrupted / "checkpoint.bin");
    bool threw = false;
    try { resume_run(interrupted, 41); }
    catch (const std::runtime_error &) { threw = true; }
    require(threw, "Resume reconstructed physical state without checkpoint.bin");
}

void old_checkpoint_and_missing_diagnostics()
{
    test::ScratchDir dir("old_checkpoint");
    const fs::path reference = new_run("reference", 31);
    const fs::path legacy = new_run("legacy", 16);
    convert_checkpoint_to_v1(legacy / "checkpoint.bin");
    resume_run(legacy, 31);
    require_identical(reference, legacy);

    const fs::path missing = copy_run(reference, "missing_columns");
    convert_checkpoint_to_v1(missing / "checkpoint.bin");
    fs::remove(missing / "phase_orb.txt");
    fs::remove(missing / "R_dot_orb.txt");
    fs::remove(missing / "p_ramp_inst.txt");
    resume_run(missing, 41);
    require(test::read_rows(missing / "phase_orb.txt").size() == 40 &&
                test::read_rows(missing / "R_dot_orb.txt").size() == 40 &&
                test::read_rows(missing / "p_ramp_inst.txt").size() == 40,
            "Resume did not backfill diagnostics absent from an old circular run");
}

void checkpoint_does_not_require_snapshot()
{
    test::ScratchDir dir("no_snapshot");
    const fs::path reference = new_run("reference", 41);
    const fs::path interrupted = new_run("interrupted", 31);
    for (const auto &entry : fs::directory_iterator(interrupted / "psi"))
        fs::remove(entry.path());
    resume_run(interrupted, 41);
    for (const char *name : diagnostics)
        require(contents(reference / name) == contents(interrupted / name),
                std::string(name) + " differs after checkpoint-only restart");
    require(contents(reference / "checkpoint.bin") == contents(interrupted / "checkpoint.bin"),
            "Final checkpoint differs after restart without text snapshots");
}

void extend_finished_run()
{
    test::ScratchDir dir("extend");
    const fs::path longer = new_run("long", 61);
    const fs::path shorter = new_run("short", 41);
    resume_run(shorter, 41);
    require(test::read_rows(shorter / "t.txt").size() == 40, "Resuming a finished run changed it");
    resume_run(shorter, 61);
    require_identical(longer, shorter);
}

void rejects_mismatched_modes()
{
    test::ScratchDir dir("mismatch");
    const fs::path folder = new_run("run", 12);
    std::ofstream(folder / "modes.txt", std::ios::trunc) << "0,0\n1,1\n";
    bool threw = false;
    try {
        resume_run(folder, 20);
    } catch (const std::runtime_error &) {
        threw = true;
    }
    require(threw, "Resume accepted a folder with a different mode list");

    const fs::path empty = new_run("empty", 5);
    fs::remove(empty / "checkpoint.bin");
    threw = false;
    try {
        resume_run(empty, 20);
    } catch (const std::runtime_error &) {
        threw = true;
    }
    require(threw, "Resume without a checkpoint did not fail");
}

} // namespace

int main()
{
    omp_set_dynamic(0);
    omp_set_num_threads(2);
    Eigen::setNbThreads(2);
    return test::run_all(std::pair{"reject changed physics before writes",&changed_physics_rejected},std::pair{"exact resume from checkpoint", &exact_resume_from_checkpoint},
                         std::pair{"signal stop and resume", &signal_stop_and_resume},
                         std::pair{"resume after a crash mid-write", &resume_after_crash_mid_write},
                         std::pair{"reject missing checkpoint", &rejects_missing_checkpoint},
                         std::pair{"old checkpoint and missing diagnostics", &old_checkpoint_and_missing_diagnostics},
                         std::pair{"checkpoint does not require snapshot", &checkpoint_does_not_require_snapshot},
                         std::pair{"extend a finished run", &extend_finished_run},
                         std::pair{"reject mismatched folders", &rejects_mismatched_modes});
}
