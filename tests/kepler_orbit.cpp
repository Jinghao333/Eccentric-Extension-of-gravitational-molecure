#include "FixedKeplerOrbit.hpp"
#include "TestUtils.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>

using test::close;
using test::require;

namespace {

OrbitState state_at(const FixedKeplerOrbit &trajectory, double mean_anomaly,
                    double clock_factor = 1.0)
{
    OrbitState state;
    state.M1 = 1.0;
    state.M2 = 0.4;
    state.phase = mean_anomaly;
    trajectory.evaluate(state, clock_factor);
    return state;
}

void circular_limit()
{
    FixedKeplerOrbit orbit(1.4, 12.0, 0.0, false);
    const OrbitState state = state_at(orbit, 7.25);
    close(state.R, 12.0, 0.0, "Circular radius changed");
    close(state.theta, 7.25, 2e-15, "Circular true anomaly differs from mean anomaly");
    close(state.R_dot, 0.0, 0.0, "Circular orbit has radial velocity");
    close(state.theta_dot, orbit.mean_motion(), 2e-15, "Circular angular velocity is wrong");
}

void kepler_residuals()
{
    double max_residual = 0.0;
    for (double e : std::array<double, 5>{0.01, 0.3, 0.8, 0.95, 0.999999})
    {
        FixedKeplerOrbit orbit(1.4, 3.0, e, false);
        for (double M : std::array<double, 7>{0.0, 1e-8, 0.3, 3.0, 5.9, 2 * constants::pi, 19.2})
        {
            const double E = orbit.eccentric_anomaly(M);
            max_residual = std::max(max_residual, std::abs(E - e * std::sin(E) - M));
            close(E - e * std::sin(E), M, 2e-14, "Kepler residual is too large");
        }
    }
    std::cout << "  max |E-e sin(E)-M| = " << max_residual << '\n';
}

void apsides_and_velocities()
{
    constexpr double total_mass = 1.4, pericenter = 8.0, e = 0.3;
    FixedKeplerOrbit trajectory(total_mass, pericenter, e, false);
    const double a = pericenter / (1.0 - e);
    const OrbitState peri = state_at(trajectory, 0.0);
    const OrbitState apo = state_at(trajectory, constants::pi);

    close(peri.R, a * (1.0 - e), 2e-15, "Pericenter is wrong");
    close(apo.R, a * (1.0 + e), 2e-15, "Apocenter is wrong");
    close(peri.R_dot, 0.0, 2e-15, "Radial velocity is nonzero at pericenter");
    close(apo.R_dot, 0.0, 2e-15, "Radial velocity is nonzero at apocenter");
    close(peri.theta_dot,
          trajectory.mean_motion() * std::sqrt(1.0 - e * e) / std::pow(1.0 - e, 2),
          2e-15, "Pericenter angular velocity is wrong");

    const OrbitState quarter = state_at(trajectory, 0.5 * constants::pi);
    const double E = trajectory.eccentric_anomaly(0.5 * constants::pi);
    const double denominator = 1.0 - e * std::cos(E);
    close(quarter.R_dot, a * trajectory.mean_motion() * e * std::sin(E) / denominator,
          2e-15, "Radial velocity is wrong");
}

void invariants_and_full_orbit()
{
    constexpr double total_mass = 1.4, pericenter = 7.0, e = 0.7;
    FixedKeplerOrbit trajectory(total_mass, pericenter, e, false);
    const double a = trajectory.semimajor_axis();
    const double expected_energy = -total_mass / (2.0 * a);
    const double expected_h = std::sqrt(total_mass * a * (1.0 - e * e));
    double max_energy_error = 0.0, max_h_error = 0.0;
    for (int k = 0; k <= 100; ++k)
    {
        const OrbitState state = state_at(trajectory, 2.0 * constants::pi * k / 100.0);
        const double tangential = state.R * state.theta_dot;
        const double energy = 0.5 * (state.R_dot * state.R_dot + tangential * tangential) -
                              total_mass / state.R;
        const double h = state.R * state.R * state.theta_dot;
        max_energy_error = std::max(max_energy_error, std::abs(energy - expected_energy));
        max_h_error = std::max(max_h_error, std::abs(h - expected_h));
        close(energy, expected_energy, 2e-13, "Kepler energy is not conserved");
        close(h, expected_h, 2e-13,
              "Kepler angular momentum is not conserved");
    }

    const OrbitState start = state_at(trajectory, 0.0);
    const OrbitState end = state_at(trajectory, 2.0 * constants::pi);
    close(end.R, start.R, 2e-14, "Radius does not return after one orbit");
    close(end.theta, 2.0 * constants::pi, 2e-14, "Angle is not unwrapped after one orbit");
    std::cout << "  max Kepler energy error = " << max_energy_error
              << ", max angular-momentum error = " << max_h_error << '\n';
}

void counter_rotation_and_ramped_clock()
{
    FixedKeplerOrbit co(1.4, 9.0, 0.4, false);
    FixedKeplerOrbit counter(1.4, 9.0, 0.4, true);
    const OrbitState a = state_at(co, 8.0, 0.35);
    const OrbitState b = state_at(counter, 8.0, 0.35);
    close(b.R, a.R, 0.0, "Counter rotation changed the radius");
    close(b.R_dot, a.R_dot, 0.0, "Counter rotation changed radial motion");
    close(b.theta, -a.theta, 2e-15, "Counter angle is not mirrored");
    close(b.theta_dot, -a.theta_dot, 2e-15, "Counter angular velocity is not mirrored");

    OrbitState state;
    state.phase = 2.0;
    co.advance_phase(state, 0.7, 0.25);
    close(state.phase, 2.0 + co.mean_motion() * 0.25 * 0.7, 2e-15,
          "Ramped mean-anomaly clock advanced incorrectly");
    require(state.phase > 2.0, "Unwrapped phase did not increase");
}

} // namespace

int main()
{
    return test::run_all(
        std::pair{"circular limit", &circular_limit},
        std::pair{"Kepler residuals", &kepler_residuals},
        std::pair{"apsides and velocities", &apsides_and_velocities},
        std::pair{"invariants and full orbit", &invariants_and_full_orbit},
        std::pair{"counter rotation and ramped clock", &counter_rotation_and_ramped_clock});
}
