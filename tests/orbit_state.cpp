// Focused unit tests for the minimal authoritative orbital state.
#include "OrbitState.hpp"
#include "TestUtils.hpp"

#include <cmath>

using test::close;
using test::require;

namespace {

void masses_and_com_geometry()
{
    OrbitState orbit;
    orbit.M1 = 2.0;
    orbit.M2 = 0.5;
    orbit.R = 10.0;

    close(orbit.total_mass(), 2.5, 0.0, "Incorrect total mass");
    close(orbit.q(), 0.25, 0.0, "Incorrect mass ratio");
    close(orbit.R1(), 2.0, 0.0, "Incorrect primary COM radius");
    close(orbit.R2(), 8.0, 0.0, "Incorrect secondary COM radius");

    const double x1 = orbit.R1();
    const double x2 = -orbit.R2();
    close(x1 - x2, orbit.R, 0.0, "COM positions do not have separation R");
    close(orbit.M1 * x1 + orbit.M2 * x2, 0.0, 0.0,
          "Primary-at-positive-R1 convention does not place the COM at zero");
}

void circular_kepler_quantities()
{
    OrbitState orbit;
    orbit.M1 = 1.0;
    orbit.M2 = 0.2;
    orbit.R = 12.0;

    const double expected_frequency = std::sqrt(1.2 / (12.0 * 12.0 * 12.0));
    close(orbit.circular_kepler_frequency(), expected_frequency, 1e-15,
          "Incorrect circular Kepler frequency");
    close(orbit.circular_kepler_period(), 2.0 * constants::pi / expected_frequency, 1e-15,
          "Incorrect circular Kepler period");
}

void signed_angle_and_positive_phase_are_distinct()
{
    OrbitState orbit;
    orbit.theta = -0.4;
    orbit.theta_dot = -0.125;
    orbit.phase = 2.75;

    require(orbit.theta_dot < 0.0, "Counter-rotating angular velocity lost its sign");
    close(orbit.angular_speed(), 0.125, 0.0, "Angular speed is not |theta_dot|");
    require(orbit.phase > 0.0, "Rolling-average phase is not independently positive");
}

void isolated_source_requires_explicit_geometry()
{
    OrbitState orbit;
    orbit.M1 = 1.0;
    orbit.M2 = 0.0;
    orbit.R = 7.5; // displaced isolated source location, not a binary separation

    close(orbit.q(), 0.0, 0.0, "Isolated state has nonzero mass ratio");
    close(orbit.R1(), 0.0, 0.0,
          "Ordinary COM getter should place the sole primary at the COM");
    close(orbit.R2(), orbit.R, 0.0,
          "Degenerate COM relation changed for M2 = 0");
    require(orbit.R != orbit.R1(),
            "Displaced isolated source was silently treated as ordinary COM geometry");
}

} // namespace

int main()
{
    return test::run_all(
        std::pair{"masses and COM geometry", &masses_and_com_geometry},
        std::pair{"circular Kepler quantities", &circular_kepler_quantities},
        std::pair{"signed angle and positive phase are distinct", &signed_angle_and_positive_phase_are_distinct},
        std::pair{"isolated source requires explicit geometry", &isolated_source_requires_explicit_geometry});
}
