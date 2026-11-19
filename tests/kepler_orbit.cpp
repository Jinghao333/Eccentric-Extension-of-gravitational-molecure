#include "FixedKeplerOrbit.hpp"
#include "CrankNicolsonModes.hpp"
#include "TestUtils.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>

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

namespace {
CrankNicolsonModes fixed_solver(double e, bool counter, double ramp_periods)
{
    // Guaranteed copy elision; callers select fixed motion on the solver.
    return CrankNicolsonModes(16, 64, 1, 0.1, 1, 0.2, 0.01, 0.4, 0,
                              8, counter, false, 1, ramp_periods, 2, "", e);
}

void prescribed_geometry()
{
    for (double e : {0.0, 1e-12, 0.01, 0.4, 0.9}) {
        for (bool counter : {false, true}) {
            auto solver = fixed_solver(e, counter, 0.0);
            solver.evolve_orbit = false;
            const double a = 8.0 / (1 - e), n = std::sqrt(1.4 / (a*a*a));
            const double sign = counter ? -1 : 1;
            // Specify E first, giving an independent Kepler solution with known M.
            for (double E : {0.0, 1.2, constants::pi, 4.0}) {
                const double turns = 6 * constants::pi;
                const double M = turns + E - e*std::sin(E);
                OrbitState start = solver.orbit;
                start.phase = turns;
                start.M1 = 1.1;
                start.M2 = 0.3;
                const double end = 1 + (M - turns) / n;
                const auto s = solver.evaluate_prescribed_orbit_at(start, 1, end);
                const double d = 1 - e*std::cos(E), b = std::sqrt(1-e*e);
                const double f = 2*std::atan2(std::sqrt(1+e)*std::sin(E/2),
                                              std::sqrt(1-e)*std::cos(E/2));
                close(s.phase, M, 1e-14, "Zero-ramp phase");
                close(s.R, a*d, 2e-13, "Separation at specified anomaly");
                close(s.theta, sign*(turns+f), 2e-13, "Signed unwrapped angle");
                close(s.R_dot, a*n*e*std::sin(E)/d, 2e-13, "Radial speed");
                close(s.theta_dot, sign*n*b/(d*d), 2e-13, "Angular speed");
                require(start.phase == turns && s.M1 == start.M1 && s.M2 == start.M2,
                        "Starting phase or supplied masses changed");
                if (e <= 1e-12) {
                    close(s.theta, sign*M, 3e-12, "Near-circular angle");
                    close(s.R, 8.0, 3e-12, "Near-circular separation");
                }
            }
        }
    }
}

void integrated_startup()
{
    for (double e : {0.0, 0.4}) {
        auto solver = fixed_solver(e, true, 1.0);
        solver.evolve_orbit = false;
        const double T = solver.initial_orbital_period(), n = solver.mean_motion();
        OrbitState start = solver.orbit;
        // E=pi/2 gives nonzero radial speed on the eccentric trajectory.
        const double target = 4*constants::pi + constants::pi/2 - e;
        for (bool switch_potential : {false, true}) {
            solver.ramp_secondary_potential = switch_potential;
            for (const auto interval : {std::pair{0.,0.}, {0.,0.5}, {0.25,1.},
                                        {0.75,1.125}, {1.25,1.5}}) {
                const auto [u,v] = interval;
                // 1.125 is the midpoint of the step [0.75,1.5] in units of T.
                const auto F = [](double x) { return x < 1 ? std::pow(x,3)-0.5*std::pow(x,4) : x-0.5; };
                const double g = v < 1 ? v*v*(3-2*v) : 1;
                start.phase = target - n*T*(F(v)-F(u));
                const auto end = solver.evaluate_prescribed_orbit_at(start, u*T, v*T);
                close(end.phase, target, 2e-15, "Integrated startup phase");
                close(end.R_dot, g*8/(1-e)*n*e, 2e-14, "Ramped radial velocity");
                close(end.theta_dot, -g*n*std::sqrt(1-e*e), 2e-14, "Ramped angular velocity");
                const double mid = (u+v)*T/2;
                const auto half = solver.evaluate_prescribed_orbit_at(start, u*T, mid);
                const auto split = solver.evaluate_prescribed_orbit_at(half, mid, v*T);
                close(split.phase, end.phase, 2e-15, "Split interval phase");
                close(split.R, end.R, 2e-14, "Split interval radius");
                close(split.theta, end.theta, 2e-14, "Split interval angle");
            }
        }
    }
}

void evaluation_errors()
{
    auto solver = fixed_solver(0, false, 1);
    solver.evolve_orbit = false;
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    for (const auto [begin,end] : {std::pair{-1.,0.}, {1.,0.}, {nan,1.},
                                   {0.,nan}, {inf,inf}, {0.,inf}}) {
        try { (void)solver.evaluate_prescribed_orbit_at(solver.orbit, begin, end); }
        catch (const std::invalid_argument &) { continue; }
        throw std::runtime_error("Invalid evaluation times accepted");
    }
    solver.evolve_orbit = true;
    for (bool isolated : {false, true}) {
        if (isolated) solver.configure_isolated_black_hole(-4);
        try { (void)solver.evaluate_prescribed_orbit_at(solver.orbit, 0, 1); }
        catch (const std::logic_error &) { continue; }
        throw std::runtime_error("Unsupported orbit mode accepted");
    }
}
} // namespace

int main()
{
    return test::run_all(
        std::pair{"prescribed geometry", &prescribed_geometry},
        std::pair{"integrated startup", &integrated_startup},
        std::pair{"evaluation errors", &evaluation_errors},
        std::pair{"circular limit", &circular_limit},
        std::pair{"Kepler residuals", &kepler_residuals},
        std::pair{"apsides and velocities", &apsides_and_velocities},
        std::pair{"invariants and full orbit", &invariants_and_full_orbit},
        std::pair{"counter rotation and ramped clock", &counter_rotation_and_ramped_clock});
}
