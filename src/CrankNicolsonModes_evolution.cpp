/**
 * @file CrankNicolsonModes_evolution.cpp
 * @brief Orbital evolution: gravitational-wave decay, cloud back-reaction and accretion.
 */
#include "CrankNicolsonModes.hpp"
#include <algorithm>
#include <cmath>

namespace
{
double orbital_clock_factor(double time, double duration)
{
    if (time >= duration) return 1.0;
    const double x = time / duration;
    return std::clamp(3.0 * x * x - 2.0 * x * x * x, 0.0, 1.0);
}

// F(t) = T*(x^3 - x^4/2) inside startup, t - T/2 thereafter.
// Evaluate F(end)-F(begin) directly so nearby large times retain their difference.
double integrated_clock_interval(double begin, double end, double duration)
{
    if (begin >= duration) return end - begin; // includes T=0
    const double ramp_stop = std::min(end, duration);
    const double x = begin / duration;
    const double y = ramp_stop / duration;
    // Divide the differences of cubes and fourth powers by (y-x) analytically.
    const double mean_factor = x*x + x*y + y*y - 0.5*(x+y)*(x*x + y*y);
    return (ramp_stop - begin) * mean_factor + (end - ramp_stop);
}
} // namespace

OrbitState CrankNicolsonModes::evaluate_prescribed_orbit_at(
    const OrbitState &start, double t_start, double t_eval) const
{
    const auto mode = orbit_mode();
    if (mode != OrbitMode::FixedCircular && mode != OrbitMode::FixedEccentric)
        throw std::logic_error("Prescribed orbit evaluation requires a fixed circular or eccentric binary");
    if (!std::isfinite(t_start) || !std::isfinite(t_eval) || t_start < 0.0 || t_eval < t_start)
        throw std::invalid_argument("Prescribed orbit times must be finite with 0 <= t_start <= t_eval");

    const double duration = ramp_end_time();
    OrbitState result = start;
    result.phase += fixed_kepler.mean_motion() * integrated_clock_interval(t_start, t_eval, duration);
    fixed_kepler.evaluate(result, orbital_clock_factor(t_eval, duration));
    return result;
}

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
    if (integrator == FieldIntegrator::MidpointCN && (isolated_black_hole || evolve_orbit))
        throw std::invalid_argument("midpoint_cn requires a fixed circular or fixed eccentric binary");
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
