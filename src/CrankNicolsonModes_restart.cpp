/**
 * @file CrankNicolsonModes_restart.cpp
 * @brief Checkpoints and resuming a run from its output folder.
 */
#include "CrankNicolsonModes.hpp"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#endif

namespace
{
constexpr char checkpoint_magic_v1[8] = {'G', 'M', 'C', 'K', 'P', 'T', '0', '1'};
constexpr char checkpoint_magic_v2[8] = {'G', 'M', 'C', 'K', 'P', 'T', '0', '2'};

template <class T>
void write_pod(std::ostream &out, const T &value)
{
    out.write(reinterpret_cast<const char *>(&value), sizeof(T));
}

template <class T>
void read_pod(std::istream &in, T &value)
{
    in.read(reinterpret_cast<char *>(&value), sizeof(T));
}

// Byte offsets just past each '\n'; a trailing partial line is not counted.
std::vector<std::uintmax_t> line_ends(const std::string &path)
{
    std::ifstream in(path, std::ios::binary);
    if (!in)
        throw std::runtime_error("Cannot open " + path);

    std::vector<std::uintmax_t> ends;
    std::vector<char> buffer(1 << 20);
    std::uintmax_t offset = 0;
    while (in)
    {
        in.read(buffer.data(), static_cast<std::streamsize>(buffer.size()));
        const std::streamsize got = in.gcount();
        for (std::streamsize i = 0; i < got; ++i)
            if (buffer[static_cast<std::size_t>(i)] == '\n')
                ends.push_back(offset + static_cast<std::uintmax_t>(i) + 1);
        offset += static_cast<std::uintmax_t>(got);
    }
    return ends;
}

// First n_rows values of a one-column diagnostic file.
std::vector<double> read_column(const std::string &path, int n_rows)
{
    std::ifstream in(path);
    if (!in)
        throw std::runtime_error("Cannot open " + path);

    std::vector<double> values;
    values.reserve(static_cast<std::size_t>(n_rows));
    std::string line;
    while (static_cast<int>(values.size()) < n_rows && std::getline(in, line))
        values.push_back(std::strtod(line.c_str(), nullptr));

    if (static_cast<int>(values.size()) != n_rows)
        throw std::runtime_error("Too few rows in " + path);
    return values;
}

void write_column(const std::string &path, const std::vector<double> &values)
{
    std::ofstream out(path, std::ios::trunc);
    if (!out)
        throw std::runtime_error("Cannot create " + path);
    out << std::setprecision(12);
    for (double value : values)
        out << value << '\n';
}

// Old circular run folders predate the explicit phase, radial-velocity and
// switching-work diagnostics. Backfill only those missing columns so their
// row alignment and diagnostic-only restart behavior remain usable.
void backfill_orbit_diagnostics(const std::string &folder, bool isolated,
                                double initial_angular_speed,
                                double ramp_end)
{
    namespace fs = std::filesystem;
    const std::string prefix = folder + "/";
    const int rows = static_cast<int>(line_ends(prefix + "t.txt").size());
    if (!fs::exists(prefix + "R_dot_orb.txt"))
        write_column(prefix + "R_dot_orb.txt", std::vector<double>(rows, 0.0));
    if (!fs::exists(prefix + "p_ramp_inst.txt"))
        write_column(prefix + "p_ramp_inst.txt", std::vector<double>(rows, 0.0));
    if (!fs::exists(prefix + "phase_orb.txt"))
    {
        const std::vector<double> times = read_column(prefix + "t.txt", rows);
        const std::vector<double> omegas = read_column(prefix + "Omega_orb.txt", rows);
        std::vector<double> phases(rows, 0.0);
        double phase = 0.0;
        double previous_time = 0.0;
        double previous_omega = isolated ? 0.0 : (ramp_end > 0.0 ? 0.0 : initial_angular_speed);
        for (int k = 0; k < rows; ++k)
        {
            phase += previous_omega * (times[k] - previous_time);
            phases[k] = phase;
            previous_time = times[k];
            previous_omega = omegas[k];
        }
        write_column(prefix + "phase_orb.txt", phases);
    }
}

}

std::vector<std::pair<std::ofstream *, std::string>> CrankNicolsonModes::diagnostic_files()
{
    std::vector<std::pair<std::ofstream *, std::string>> files = {
        {&flux_file, "flux.txt"},
        {&energy_flux_file, "energy_flux.txt"},
        {&energy_file, "energy.txt"},
        {&energy_modes_file, "energy_modes.txt"},
        {&norms_file, "norms.txt"},
        {&time_file, "t.txt"},
        {&p_bin_file, "p_bin.txt"},
        {&p_bin_inst_file, "p_bin_inst.txt"},
        {&d_V_bin_dt_file, "dV_bin_dt.txt"},
        {&p_ramp_inst_file, "p_ramp_inst.txt"},
        {&R_file, "R_orb.txt"},
        {&R_dot_file, "R_dot_orb.txt"},
        {&Omega_file, "Omega_orb.txt"},
        {&M_1_file, "M_1.txt"},
        {&q_file, "q.txt"},
        {&theta_file, "theta_orb.txt"},
        {&phase_file, "phase_orb.txt"},
    };
    if (geometry.enabled) files.emplace_back(&dt_step_file,"dt_step.txt");
    return files;
}

void CrankNicolsonModes::open_diagnostics(std::ios::openmode mode)
{
    for (auto &[stream, name] : diagnostic_files())
    {
        stream->open(output_folder + "/" + name, mode);
        if (!*stream)
            throw std::runtime_error("Cannot open " + output_folder + "/" + name);
    }
}

// Captured from the effective solver settings before either run or restore.
std::string CrankNicolsonModes::describe_initial_physics() const
{
    std::ostringstream out;
    out.imbue(std::locale::classic());
    out << "GM-PHYSICS-1" << std::setprecision(17);
    for (double x : {double(integrator), double(orbit_mode()), double(Nr), L, epsilon,
                     double(l_max), mu, M_c, a, orbit.M1, orbit.M2, orbit.R, orbit.theta,
                     eccentricity(), double(counter), geometry.normalized ? 0. : (geometry.enabled ? maximum_dt : dt), ramp_periods, double(ramp_secondary_potential),
                     steps_per_period, orbits_average, double(absorption), epsilon_smear,
                     double(use_block_preconditioner), r_ex, double(mask), double(cap)})
        out << ' ' << x;
    for (const auto &layer : {mask_settings,cap_settings})
        out << ' ' << layer.factor << ' ' << layer.strength << ' ' << layer.exponent;
    if (geometry.normalized)
        out << " normalized_geometry_v1 " << geometry.startup_steps << ' ' << geometry.orbit_steps << " 2 3";
    else if (geometry.enabled)
        out << " orbital_geometry_v1 " << geometry.tolerance << ' ' << geometry.azimuthal_order
            << ' ' << geometry.radial_order << ' ' << geometry.ramp_tolerance;
    return out.str();
}

void CrankNicolsonModes::save_checkpoint(int step_id) const
{
    namespace fs = std::filesystem;
    const std::string final_path = output_folder + "/checkpoint.bin";
    const std::string tmp_path = final_path + ".tmp";

    {
        std::ofstream out(tmp_path, std::ios::binary | std::ios::trunc);
        if (!out)
            throw std::runtime_error("Could not write checkpoint " + tmp_path);

        out.write(checkpoint_magic_v2, sizeof(checkpoint_magic_v2));
        const std::int32_t header[3] = {step_id, Nr, n_modes};
        out.write(reinterpret_cast<const char *>(header), sizeof(header));
        for (const double value : {t_i, orbit.theta, orbit.phase, orbit.R,
                                   orbit.R_dot, orbit.theta_dot, orbit.M1,
                                   orbit.M2, dt, p_bin})
            write_pod(out, value);

        const std::uint64_t size = static_cast<std::uint64_t>(psi.size());
        write_pod(out, size);
        out.write(reinterpret_cast<const char *>(psi.data()),
                  static_cast<std::streamsize>(size * sizeof(dcomplex)));
        power_average.save(out);
        out << initial_physics;
        out.close();

        if (!out)
            throw std::runtime_error("Could not write checkpoint " + tmp_path);
    }

    // A crash mid-write leaves the previous checkpoint intact.
    std::error_code error;
#ifdef _WIN32
    const auto src = fs::path(tmp_path).wstring(), dst = fs::path(final_path).wstring();
    for (int attempt=0; ; ++attempt) {
        if (MoveFileExW(src.c_str(),dst.c_str(),MOVEFILE_REPLACE_EXISTING | MOVEFILE_WRITE_THROUGH)) break;
        const auto code = GetLastError();
        if (attempt == 20 || (code != ERROR_SHARING_VIOLATION && code != ERROR_ACCESS_DENIED)) {
            error.assign(code,std::system_category());
            break;
        }
        Sleep(25);
    }
#else
    fs::rename(tmp_path, final_path, error);
#endif
    if (error)
        throw std::runtime_error("Could not replace checkpoint: " + error.message());
}

int CrankNicolsonModes::load_checkpoint()
{
    std::ifstream in(output_folder + "/checkpoint.bin", std::ios::binary);
    if (!in)
        return -1;

    char magic[sizeof(checkpoint_magic_v2)];
    in.read(magic, sizeof(magic));
    std::int32_t header[3] = {-1, -1, -1};
    in.read(reinterpret_cast<char *>(header), sizeof(header));
    const bool version1 = std::memcmp(magic, checkpoint_magic_v1, sizeof(magic)) == 0;
    const bool version2 = std::memcmp(magic, checkpoint_magic_v2, sizeof(magic)) == 0;
    if (!in || (!version1 && !version2) || header[0] < 0 ||
        header[1] != Nr || header[2] != n_modes)
        return -1;

    // GMCKPT01 predates eccentric trajectories: it has no radial velocity and
    // labels its angular rate/period as circular quantities.  Using it for an
    // eccentric run can silently freeze the binary at the saved separation.
    // Keep legacy circular compatibility, but require a complete v2 checkpoint
    // for eccentric motion rather than reconstructing physical state from text.
    if (version1 && eccentricity() > 0.0)
    {
        std::cerr << "Legacy GMCKPT01 cannot exactly restore an eccentric run; "
                     "a GMCKPT02 checkpoint is required."
                  << std::endl;
        return -1;
    }

    std::vector<double> values(version1 ? 13 : 10);
    for (double &value : values)
        read_pod(in, value);
    std::uint64_t size = 0;
    read_pod(in, size);
    if (!in || size != static_cast<std::uint64_t>(Nr) * static_cast<std::uint64_t>(n_modes))
        return -1;

    Eigen::VectorXcd state(static_cast<Eigen::Index>(size));
    in.read(reinterpret_cast<char *>(state.data()),
            static_cast<std::streamsize>(size * sizeof(dcomplex)));
    if (!in)
        return -1;

    RollingPhaseAverage average{orbits_average};
    try
    {
        average.load(in);
    }
    catch (const std::exception &)
    {
        return -1;
    }

    std::string recorded;
    std::getline(in,recorded);
    if (recorded.empty() && in.eof()) {
        // Only an independently archived original input can attest settings
        // (e.g. mu) absent from an old checkpoint. Never use editable config.toml.
        const auto original = std::filesystem::path(output_folder)/"trusted_original_config.toml";
        if (!std::filesystem::is_regular_file(original))
            throw std::runtime_error("Unverified legacy checkpoint: supply the archived original input as trusted_original_config.toml");
        const auto c = load_config(original.string());
        if (geometry.enabled || c.geometry.enabled || isolated_black_hole || c.cap || c.mask || cap || mask ||
            c.integrator != integrator || c.N != Nr || c.L != L || c.epsilon != epsilon ||
            c.l_max != l_max || c.mu != mu || c.M_c != M_c || c.a != a ||
            c.q != orbit.M2 || c.e != eccentricity() ||
            c.R/(1-c.e) != fixed_kepler.semimajor_axis() ||
            c.counter != counter || c.dt != dt || c.evolve_orbit != evolve_orbit ||
            c.ramp_periods != ramp_periods || c.ramp_secondary_potential != ramp_secondary_potential ||
            c.steps_per_period != steps_per_period || c.power_average_orbits != orbits_average ||
            c.absorption != absorption || c.epsilon_smear != epsilon_smear ||
            (c.preconditioner == "radial") != use_block_preconditioner || c.r_ex != r_ex)
            throw std::runtime_error("Legacy checkpoint settings differ or its boundary model cannot be reproduced");
    } else if (recorded != initial_physics || !in.eof())
        throw std::runtime_error("Checkpoint physics settings disagree with this solver or the record is damaged");

    t_i = values[0];
    orbit.theta = values[1];
    orbit.phase = values[2];
    orbit.R = values[3];
    if (version1)
    {
        orbit.R_dot = 0.0;
        orbit.M1 = values[4];
        orbit.M2 = values[5];
        orbit.theta_dot = (counter ? -1.0 : 1.0) * values[9];
        dt = values[11];
        p_bin = values[12];
    }
    else
    {
        orbit.R_dot = values[4];
        orbit.theta_dot = values[5];
        orbit.M1 = values[6];
        orbit.M2 = values[7];
        dt = values[8];
        p_bin = values[9];
    }
    psi = state;
    power_average = average;
    return header[0];
}

void CrankNicolsonModes::resume_until(const std::string &folder, double final_time, double extraction)
{
    set_time_target(final_time);
    resume_impl(folder, 1, extraction);
}

void CrankNicolsonModes::resume(const std::string &folder, int n_steps, double extraction)
{
    if (geometry.enabled) throw std::invalid_argument("Use resume_until for orbital_geometry");
    resume_impl(folder, n_steps, extraction);
}

void CrankNicolsonModes::resume_impl(const std::string &folder, int n_steps, double r_ex_)
{
    namespace fs = std::filesystem;

    std::cout << "Resuming simulation in " << folder << std::endl;

    const auto marker_path = fs::path(folder) / "integrator.txt";
    const bool midpoint_run = fs::exists(marker_path);
    if (midpoint_run)
    {
        std::ifstream marker(marker_path);
        std::string name;
        if (!std::getline(marker, name) || name != "midpoint_cn")
            throw std::runtime_error("Unrecognized integrator in " + marker_path.string());
    }
    if (midpoint_run != (integrator == FieldIntegrator::MidpointCN))
        throw std::runtime_error("Cannot change simulation.integrator when resuming a run");

    set_extraction_radius(r_ex_, n_steps);
    initial_physics = describe_initial_physics();
    output_folder = folder;

    // The stored mode ordering must match the one rebuilt from config.toml.
    {
        std::ifstream modes_in(folder + "/modes.txt");
        std::string line;
        int count = 0;
        while (std::getline(modes_in, line))
        {
            if (line.empty())
                continue;
            int l = 0, m = 0;
            char comma = 0;
            std::istringstream row(line);
            row >> l >> comma >> m;
            if (count >= n_modes || modes[count].l != l || modes[count].m != m)
                throw std::runtime_error("modes.txt does not match the modes built from config.toml");
            ++count;
        }
        if (modes_in.is_open() && count != n_modes)
            throw std::runtime_error("modes.txt does not match the modes built from config.toml");
    }

    // The binary checkpoint is the sole authority for physical restart state.
    const int start = load_checkpoint();
    if (start < 0)
        throw std::runtime_error("No usable checkpoint.bin in " + folder);
    std::cout << "Exact state restored from checkpoint.bin at step " << start << "." << std::endl;

    // Restore the endpoint ramp factor without advancing the orbit, then make
    // sure an eccentric checkpoint state belongs to the trajectory
    // defined by this run's config.toml.  This catches mixed executables,
    // changed orbital parameters and circularized legacy restart state before
    // any output is truncated.
    const double restored_R_dot = orbit.R_dot;
    const double restored_theta_dot = orbit.theta_dot;
    update_orbit_and_ramp(t_i, false);
    if (!isolated_black_hole && !evolve_orbit && eccentricity() > 0.0)
    {
        OrbitState expected = orbit;
        fixed_kepler.evaluate(expected, ramp_factor);
        const auto inconsistent = [](double saved, double reference) {
            const double scale = std::max({1.0, std::abs(saved), std::abs(reference)});
            return std::abs(saved - reference) > 5e-9 * scale;
        };
        if (inconsistent(orbit.R, expected.R) ||
            inconsistent(orbit.theta, expected.theta) ||
            inconsistent(orbit.R_dot, expected.R_dot) ||
            inconsistent(orbit.theta_dot, expected.theta_dot))
            throw std::runtime_error(
                "Restored eccentric orbit is inconsistent with config.toml; "
                "refusing to continue with a mixed or circularized trajectory");
    }

    // Complete diagnostic rows available in every file. Missing columns from
    // old circular runs can still be backfilled without reconstructing state.
    if (!midpoint_run)
        backfill_orbit_diagnostics(folder, isolated_black_hole,
                                   fixed_kepler.mean_motion(), ramp_end_time());
    std::size_t complete_rows = std::numeric_limits<std::size_t>::max();
    std::size_t maximum_rows = 0;
    std::vector<std::vector<std::uintmax_t>> ends;
    std::vector<std::string> diagnostic_names;
    for (auto &[stream, name] : diagnostic_files())
    {
        ends.push_back(line_ends(folder + "/" + name));
        diagnostic_names.push_back(name);
        complete_rows = std::min(complete_rows, ends.back().size());
        maximum_rows = std::max(maximum_rows, ends.back().size());
    }
    if (complete_rows != maximum_rows)
    {
        std::cerr << "Warning: diagnostic row counts disagree; truncating every file "
                     "to checkpoint step " << start << ":\n";
        for (std::size_t i = 0; i < ends.size(); ++i)
            if (ends[i].size() != complete_rows)
                std::cerr << "  " << diagnostic_names[i] << ": " << ends[i].size() << " rows\n";
    }

    if (complete_rows < static_cast<std::size_t>(start))
        throw std::runtime_error("Diagnostics end before checkpoint step " + std::to_string(start));

    if (geometry.enabled) {
        if (target_time < t_i || start >= std::numeric_limits<int>::max()-1)
            throw std::runtime_error("Invalid geometric restart target or step index");
        const auto durations = read_column(folder+"/dt_step.txt",start);
        const auto timestamps = geometry.normalized ? read_column(folder+"/t.txt",start)
                                                    : std::vector<double>{};
        t.resize(start+1);
        t(0)=0.;
        for (int i=0;i<start;++i) {
            t(i+1)=geometry.normalized ? timestamps[i] : t(i)+durations[i];
            if (geometry.normalized && t(i+1)-t(i)!=durations[i])
                throw std::runtime_error("Normalized node times disagree with accepted steps");
            if (!(t(i+1)>t(i))) throw std::runtime_error("Invalid accepted step history");
        }
        if (t(start)!=t_i) throw std::runtime_error("Accepted steps disagree with checkpoint time");
        dt=select_next_timestep();
    }

    // Drop diagnostics and snapshots written after the checkpoint, then append.
    int file_index = 0;
    for (auto &[stream, name] : diagnostic_files())
    {
        const auto &file_ends = ends[file_index++];
        const std::uintmax_t keep = start == 0 ? 0 : file_ends[static_cast<std::size_t>(start - 1)];
        fs::resize_file(folder + "/" + name, keep);
    }
    const std::string snapshot_prefix = "snapshot_";
    const std::string snapshot_suffix = ".txt";
    for (const auto &entry : fs::directory_iterator(folder + "/psi"))
    {
        const std::string name = entry.path().filename().string();
        if (name.rfind(snapshot_prefix, 0) != 0 ||
            name.size() <= snapshot_prefix.size() + snapshot_suffix.size() ||
            name.compare(name.size() - snapshot_suffix.size(), snapshot_suffix.size(), snapshot_suffix) != 0)
            continue;
        const std::string digits = name.substr(
            snapshot_prefix.size(), name.size() - snapshot_prefix.size() - snapshot_suffix.size());
        if (digits.find_first_not_of("0123456789") == std::string::npos &&
            std::stoi(digits) > start)
            fs::remove(entry.path());
    }
    fs::remove(folder + "/checkpoint.bin.tmp");
    open_diagnostics(std::ios::out | std::ios::app);

    if (!geometry.enabled) {
        t.resize(n_steps);
        t.setZero();
        if (start < n_steps) t(start) = t_i;
    }
    psi_last = Eigen::VectorXcd::Zero(psi.size());

    // Restore schedule weights without advancing the orbit, then rebuild the
    // endpoint operator from the authoritative checkpoint/diagnostic state.
    orbit.R_dot = restored_R_dot;
    orbit.theta_dot = restored_theta_dot;
    build_time_dependent_perturbation(orbit, hamiltonian_evaluation(t_i));

    if (geometry.enabled ? t_i == target_time : start + 1 >= n_steps)
    {
        std::cout << "The run had already reached its last step ("
                  << (geometry.enabled ? start : n_steps - 1) << ")." << std::endl;
        return;
    }

    evolve(start + 1, n_steps);
}
