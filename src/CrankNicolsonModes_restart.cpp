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
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

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

bool load_snapshot(const std::string &path, Eigen::VectorXcd &state, Eigen::Index size)
{
    std::ifstream in(path);
    if (!in)
        return false;

    Eigen::VectorXcd values(size);
    std::string line;
    Eigen::Index k = 0;
    while (std::getline(in, line))
    {
        if (line.empty())
            continue;
        if (k == size)
            return false;
        char *end = nullptr;
        const double re = std::strtod(line.c_str(), &end);
        char *end_im = nullptr;
        const double im = std::strtod(end, &end_im);
        if (end == line.c_str() || end_im == end)
            return false;
        values[k++] = dcomplex(re, im);
    }

    if (k != size)
        return false;
    state = values;
    return true;
}
}

std::vector<std::pair<std::ofstream *, std::string>> CrankNicolsonModes::diagnostic_files()
{
    return {
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

void CrankNicolsonModes::save_checkpoint(int step_id) const
{
    namespace fs = std::filesystem;
    const std::string final_path = output_folder + "/checkpoint.bin";
    const std::string tmp_path = final_path + ".tmp";

    {
        std::ofstream out(tmp_path, std::ios::binary | std::ios::trunc);
        if (!out)
        {
            std::cerr << "Could not write checkpoint " << tmp_path << "\n";
            return;
        }

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

        if (!out)
        {
            std::cerr << "Could not write checkpoint " << tmp_path << "\n";
            return;
        }
    }

    // A crash mid-write leaves the previous checkpoint intact.
    std::error_code error;
    fs::rename(tmp_path, final_path, error);
    if (error)
        std::cerr << "Could not replace checkpoint: " << error.message() << "\n";
}

bool CrankNicolsonModes::load_checkpoint(int step_id)
{
    std::ifstream in(output_folder + "/checkpoint.bin", std::ios::binary);
    if (!in)
        return false;

    char magic[sizeof(checkpoint_magic_v2)];
    in.read(magic, sizeof(magic));
    std::int32_t header[3] = {-1, -1, -1};
    in.read(reinterpret_cast<char *>(header), sizeof(header));
    const bool version1 = std::memcmp(magic, checkpoint_magic_v1, sizeof(magic)) == 0;
    const bool version2 = std::memcmp(magic, checkpoint_magic_v2, sizeof(magic)) == 0;
    if (!in || (!version1 && !version2) || header[0] != step_id ||
        header[1] != Nr || header[2] != n_modes)
        return false;

    std::vector<double> values(version1 ? 13 : 10);
    for (double &value : values)
        read_pod(in, value);
    std::uint64_t size = 0;
    read_pod(in, size);
    if (!in || size != static_cast<std::uint64_t>(Nr) * static_cast<std::uint64_t>(n_modes))
        return false;

    Eigen::VectorXcd state(static_cast<Eigen::Index>(size));
    in.read(reinterpret_cast<char *>(state.data()),
            static_cast<std::streamsize>(size * sizeof(dcomplex)));
    if (!in)
        return false;

    RollingPhaseAverage average{orbits_average};
    try
    {
        average.load(in);
    }
    catch (const std::exception &)
    {
        return false;
    }

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
    return true;
}

// For runs without a matching checkpoint: rebuild the evolving scalars from
// the diagnostic rows 0..step_id-1 (row k is written by step k+1). These are
// stored with 12 significant digits, so the restart is accurate, not bitwise.
void CrankNicolsonModes::reconstruct_state_from_diagnostics(int step_id)
{
    const std::string f = output_folder + "/";
    const std::vector<double> times = read_column(f + "t.txt", step_id);
    const std::vector<double> Rs = read_column(f + "R_orb.txt", step_id);
    const std::vector<double> thetas = read_column(f + "theta_orb.txt", step_id);
    const std::vector<double> Omegas = read_column(f + "Omega_orb.txt", step_id);
    const std::vector<double> M1s = read_column(f + "M_1.txt", step_id);
    const std::vector<double> qs = read_column(f + "q.txt", step_id);
    const std::vector<double> phases = read_column(f + "phase_orb.txt", step_id);
    const std::vector<double> R_dots = read_column(f + "R_dot_orb.txt", step_id);
    const std::vector<double> p_bins = read_column(f + "p_bin.txt", step_id);
    const std::vector<double> p_bin_insts = read_column(f + "p_bin_inst.txt", step_id);

    power_average.reset();
    for (int k = 0; k < step_id; ++k)
    {
        if (!isolated_black_hole)
            power_average.update(times[k], phases[k], p_bin_insts[k]);
    }

    if (step_id == 0)
        return; // The constructor's state is the initial state.

    const int last = step_id - 1;
    t_i = times[last];
    orbit.R = Rs[last];
    orbit.R_dot = R_dots[last];
    orbit.theta = thetas[last];
    orbit.phase = phases[last];
    orbit.M1 = M1s[last];
    orbit.M2 = qs[last] * orbit.M1;
    orbit.theta_dot = (counter ? -1.0 : 1.0) * Omegas[last];
    p_bin = p_bins[last];

    // Circular feedback switches to period-based dt after the startup gate.
    if (!isolated_black_hole && evolve_orbit && t_i > ramp_end_time())
        dt = orbit.circular_kepler_period() / steps_per_period;
}

void CrankNicolsonModes::resume(const std::string &folder, int n_steps, double r_ex_)
{
    namespace fs = std::filesystem;

    std::cout << "Resuming simulation in " << folder << std::endl;

    set_extraction_radius(r_ex_, n_steps);
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

    // Complete diagnostic rows available in every file.
    backfill_orbit_diagnostics(folder, isolated_black_hole,
                               fixed_kepler.mean_motion(), ramp_end_time());
    std::size_t complete_rows = std::numeric_limits<std::size_t>::max();
    std::vector<std::vector<std::uintmax_t>> ends;
    for (auto &[stream, name] : diagnostic_files())
    {
        ends.push_back(line_ends(folder + "/" + name));
        complete_rows = std::min(complete_rows, ends.back().size());
    }

    // Latest snapshot that is complete and has all of its diagnostic rows.
    std::vector<int> snapshot_steps;
    const std::regex pattern(R"(snapshot_(\d+)\.txt)");
    for (const auto &entry : fs::directory_iterator(folder + "/psi"))
    {
        std::smatch match;
        const std::string name = entry.path().filename().string();
        if (std::regex_match(name, match, pattern))
            snapshot_steps.push_back(std::stoi(match[1].str()));
    }
    std::sort(snapshot_steps.rbegin(), snapshot_steps.rend());

    int start = -1;
    for (const int candidate : snapshot_steps)
    {
        if (static_cast<std::size_t>(candidate) > complete_rows)
            continue;
        std::ostringstream path;
        path << folder << "/psi/snapshot_" << std::setw(6) << std::setfill('0') << candidate << ".txt";
        if (load_snapshot(path.str(), psi, static_cast<Eigen::Index>(Nr) * n_modes))
        {
            start = candidate;
            break;
        }
        std::cout << "Skipping incomplete snapshot " << path.str() << std::endl;
    }
    if (start < 0)
        throw std::runtime_error("No usable snapshot found in " + folder + "/psi");

    std::cout << "Restarting from snapshot " << start << std::endl;

    if (load_checkpoint(start))
    {
        std::cout << "Exact state restored from checkpoint.bin." << std::endl;
    }
    else
    {
        std::cout << "No checkpoint for this snapshot: rebuilding the orbit state from the "
                     "diagnostics (snapshot and diagnostic precision limit the match)."
                  << std::endl;
        reconstruct_state_from_diagnostics(start);
    }

    // Drop everything written after the snapshot, then append.
    int file_index = 0;
    for (auto &[stream, name] : diagnostic_files())
    {
        const auto &file_ends = ends[file_index++];
        const std::uintmax_t keep = start == 0 ? 0 : file_ends[static_cast<std::size_t>(start - 1)];
        fs::resize_file(folder + "/" + name, keep);
    }
    open_diagnostics(std::ios::out | std::ios::app);

    t.resize(n_steps);
    t.setZero();
    if (start < n_steps)
        t(start) = t_i;
    psi_last = Eigen::VectorXcd::Zero(psi.size());

    // Restore schedule weights without advancing the orbit, then rebuild the
    // endpoint operator from the authoritative checkpoint/diagnostic state.
    const double restored_R_dot = orbit.R_dot;
    const double restored_theta_dot = orbit.theta_dot;
    update_orbit_and_ramp(t_i, false);
    orbit.R_dot = restored_R_dot;
    orbit.theta_dot = restored_theta_dot;
    build_time_dependent_perturbation(orbit, t_i);

    if (start + 1 >= n_steps)
    {
        std::cout << "The run had already reached its last step (" << n_steps - 1 << ")." << std::endl;
        return;
    }

    evolve(start + 1, n_steps);
}
