/**
 * @file CrankNicolsonModes_evolution.cpp
 * @brief Orbital evolution: gravitational-wave decay, cloud back-reaction and accretion.
 */
#include "CrankNicolsonModes.hpp"
#include <algorithm>
#include <cmath>

CrankNicolsonModes::OrbitMode CrankNicolsonModes::orbit_mode() const noexcept
{
    if (isolated_black_hole) return OrbitMode::Isolated;
    if (evolve_orbit) return OrbitMode::CircularInspiralLegacy;
    return eccentricity() > 0.0 ? OrbitMode::FixedEccentric : OrbitMode::FixedCircular;
}

CrankNicolsonModes::TimeStepRule CrankNicolsonModes::timestep_rule() const noexcept
{
    return orbit_mode() == OrbitMode::CircularInspiralLegacy
        ? TimeStepRule::LegacyCircularInspiral : TimeStepRule::FixedDt;
}

double CrankNicolsonModes::orbital_period() const noexcept
{
    if (!evolve_orbit || eccentricity() > 0.0)
        return fixed_kepler.period();
    return orbit.circular_kepler_period();
}

void CrankNicolsonModes::validate_orbit_configuration() const
{
    if (eccentricity() > 0.0 && evolve_orbit)
        throw std::invalid_argument(
            "Eccentric inspiral is not implemented: set evolve_orbit=false when e>0");
}

void CrankNicolsonModes::update_orbit_and_ramp(double time, bool evaluate_geometry)
{
    ramp_factor = 1.0;
    const double duration = ramp_end_time();
    if (!isolated_black_hole && duration > 0.0 && time < duration)
    {
        const double x = std::clamp(time / duration, 0.0, 1.0);
        ramp_factor = 3.0 * x * x - 2.0 * x * x * x;
    }

    if (orbit_mode() == OrbitMode::Isolated)
    {
        orbit.R_dot = 0.0;
        orbit.theta_dot = 0.0;
        return;
    }

    if (orbit_mode() != OrbitMode::CircularInspiralLegacy)
    {
        if (evaluate_geometry)
            fixed_kepler.evaluate(orbit, ramp_factor);
        return;
    }

    // Legacy circular inspiral: startup changes only the angular rate.  The
    // post-field feedback update remains responsible for R, masses and R_dot.
    orbit.theta_dot = (counter ? -1.0 : 1.0) *
                      (time < duration ? fixed_kepler.mean_motion() * ramp_factor
                                       : orbit.circular_kepler_frequency());
}

void CrankNicolsonModes::advance_orbit(double step_dt, double endpoint_time)
{
    if (orbit_mode() == OrbitMode::Isolated)
    {
        update_orbit_and_ramp(endpoint_time, false);
        return;
    }

    if (orbit_mode() != OrbitMode::CircularInspiralLegacy)
    {
        // The old endpoint convention is deliberate: phase uses the previous
        // step's clock factor, then endpoint geometry/rates use the new factor.
        fixed_kepler.advance_phase(orbit, step_dt, ramp_factor);
        update_orbit_and_ramp(endpoint_time, true);
        return;
    }

    orbit.theta += orbit.theta_dot * step_dt;
    orbit.phase += orbit.angular_speed() * step_dt;
    update_orbit_and_ramp(endpoint_time, false);
}

namespace
{
/// Quadrupole gravitational-wave decay rate of a circular orbit, dR/dt.
double dR_dt_GW(double R, double M1, double M2)
{
    return -(64.0 / 5.0) * M1 * M2 * (M1 + M2) / std::pow(R, 3);
}
} // namespace

void CrankNicolsonModes::apply_orbital_feedback(double step_start)
{
    if (orbit_mode() == OrbitMode::CircularInspiralLegacy && t_i > ramp_end_time())
    {
        // Only the part of a crossing step after activation contributes.
        const double orbit_dt = t_i - std::max(step_start, ramp_end_time());
        apply_legacy_circular_feedback(orbit_dt);
    }
}

double CrankNicolsonModes::select_next_timestep() const
{
    // Keep both the original caller's strict gate and the original dt gate.
    if (timestep_rule() == TimeStepRule::LegacyCircularInspiral &&
        t_i > ramp_end_time() && ramp_factor == 1.0)
        return orbit.circular_kepler_period() / steps_per_period;
    return dt;
}

void CrankNicolsonModes::apply_legacy_circular_feedback(double orbit_dt)
{
    // Accretion rates are evaluated on the midpoint of the step.
    Eigen::VectorXcd psi_mid = (psi + psi_last);
    psi_mid *= 0.5;

    double Mdot1 = 0.0;
    double Mdot2 = 0.0;

    if (absorption)
    {
        Mdot1 = get_absorption_rate_1(psi_mid);
        Mdot2 = get_absorption_rate_2(psi_mid);
    }

    // Orbital back-reaction turns on smoothly over one initial period after the ramp.
    double ramp = 1.0;

    if (t_i < ramp_end_time())
    {
        ramp = 0.0;
    }
    else if (t_i < ramp_end_time() + initial_period)
    {
        const double x = (t_i - ramp_end_time()) / initial_period;
        ramp = 3.0 * x * x - 2.0 * x * x * x;
    }

    const double dR_dt_gw = dR_dt_GW(orbit.R, orbit.M1, orbit.M2);
    const double dR_dt_cloud = -2.0 * orbit.R * orbit.R * p_bin /
                               (orbit.M1 * orbit.M2);

    // Accretion changes the masses below; its direct drift of the separation,
    // R * (Mdot1 / M1 + Mdot2 / M2), is not included.
    orbit.R_dot = (dR_dt_gw + dR_dt_cloud) * ramp;
    orbit.R += orbit.R_dot * orbit_dt;

    orbit.M1 += Mdot1 * orbit_dt;
    orbit.M2 += Mdot2 * orbit_dt;

    orbit.theta_dot = (counter ? -1.0 : 1.0) *
                      orbit.circular_kepler_frequency();
}
