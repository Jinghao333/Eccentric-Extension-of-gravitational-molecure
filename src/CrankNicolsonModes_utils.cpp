/**
 * @file CrankNicolsonModes_utils.cpp
 * @brief Output: run metadata, snapshots and diagnostic lines.
 */
#include "CrankNicolsonModes.hpp"
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <sstream>

void CrankNicolsonModes::report_step_progress(int step_id) const
{
    if (!evolve_orbit)
    {
        std::cout << "\rStep: " << step_id << "/" << t.size() << " | Time " << t_i << std::flush;
    }
    else
    {
        std::cout << "\rStep: " << step_id << "/" << t.size() << " | Time " << t_i << " | R = " << orbit.R << std::flush;
    }
}

void CrankNicolsonModes::write_step_output(int step_id, const StepObservables &observables)
{
    write_quantity_line(orbit.R, R_file);
    write_quantity_line(orbit.theta, theta_file);
    write_quantity_line(orbit.angular_speed(), Omega_file);
    write_quantity_line(orbit.M1, M_1_file);
    write_quantity_line(orbit.q(), q_file);
    write_quantity_line(orbit.phase, phase_file);
    write_quantity_line(orbit.R_dot, R_dot_file);
    write_quantity_line(t_i, time_file);
    write_quantity_line_modes(observables.energy_flux, energy_flux_file);
    write_quantity_line_modes(observables.fluxes, flux_file);
    write_quantity_line(observables.energy, energy_file);
    write_quantity_line(p_bin, p_bin_file);
    write_quantity_line(observables.binary_power, p_bin_inst_file);
    write_quantity_line(observables.switching_power, p_ramp_inst_file);
    write_quantity_line(observables.dV_bin_dt, d_V_bin_dt_file);

    write_quantity_line_modes(observables.norms, norms_file);

    write_quantity_line_modes(observables.energies, energy_modes_file);

    if (step_id % save_every == 0)
    {
        save_snapshot(step_id, psi, "psi");
        flush_diagnostics();
        save_checkpoint(step_id);
    }
}

void CrankNicolsonModes::save_run_info()
{
    namespace fs = std::filesystem;

    fs::copy_file(config_file, output_folder + "/config.toml", fs::copy_options::overwrite_existing);
    // Unmarked folders use the original endpoint integrator.
    if (integrator == FieldIntegrator::MidpointCN)
    {
        std::ofstream marker(output_folder + "/integrator.txt");
        marker << "midpoint_cn\n";
        marker.close();
        if (!marker)
            throw std::runtime_error("Cannot record midpoint integrator in " + output_folder);
    }

    std::stringstream filename;
    filename << output_folder << "/r.txt";

    std::ofstream out(filename.str());
    if (!out)
    {
        std::cerr << "Could not open file: " << filename.str() << "\n";
        return;
    }

    for (int i = 0; i < Nr; ++i)
    {
        out << r[i] << "\n";
    }

    std::stringstream modes_filename;
    modes_filename << output_folder << "/modes.txt";

    std::ofstream modes_out(modes_filename.str());
    if (!modes_out)
    {
        std::cerr << "Could not open file: " << modes_filename.str() << "\n";
        return;
    }

    for (int n = 0; n < n_modes; n++)
    {
        modes_out << modes[n].l << "," << modes[n].m << std::endl;
    }
}

void CrankNicolsonModes::save_snapshot(int step_id, const Eigen::VectorXcd &psi, std::string name)
{
    namespace fs = std::filesystem;

    std::stringstream filename;
    filename << output_folder << "/" << name << "/snapshot_" << std::setw(6)
             << std::setfill('0') << step_id << ".txt";

    std::ofstream out(filename.str());
    if (!out)
    {
        std::cerr << "Could not open file: " << filename.str() << "\n";
        return;
    }

    for (int i = 0; i < psi.size(); i++)
    {
        out << psi[i].real() << " " << psi[i].imag() << "\n";
    }
}

void CrankNicolsonModes::write_quantity_line_modes(const std::vector<double> &quantity, std::ofstream &quantity_file)
{
    for (size_t i = 0; i < quantity.size(); i++)
    {
        quantity_file << std::scientific << std::setprecision(12) << quantity[i];
        if (i + 1 < quantity.size())
            quantity_file << " ";
    }
    quantity_file << '\n';
}

void CrankNicolsonModes::write_quantity_line(const double &quantity, std::ofstream &quantity_file)
{
    quantity_file << std::scientific << std::setprecision(12) << quantity << '\n';
}

void CrankNicolsonModes::flush_diagnostics()
{
    for (auto &[stream, name] : diagnostic_files())
    {
        (void)name;
        stream->flush();
    }
}
