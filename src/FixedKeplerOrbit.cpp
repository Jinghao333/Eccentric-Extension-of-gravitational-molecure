/**
 * @file FixedKeplerOrbit.cpp
 * @brief Robust solution and evaluation of a prescribed Kepler ellipse.
 */
#include "FixedKeplerOrbit.hpp"

#include "Definitions.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

FixedKeplerOrbit::FixedKeplerOrbit(double total_mass, double pericenter,
                                   double eccentricity, bool counter_rotating)
{
    configure(total_mass, pericenter, eccentricity, counter_rotating);
}

void FixedKeplerOrbit::configure(double total_mass, double pericenter,
                                 double eccentricity, bool counter_rotating)
{
    if (!(total_mass > 0.0) || !(pericenter > 0.0) ||
        !(eccentricity >= 0.0 && eccentricity < 1.0) ||
        !std::isfinite(total_mass) || !std::isfinite(pericenter) ||
        !std::isfinite(eccentricity))
        throw std::invalid_argument(
            "FixedKeplerOrbit requires finite M>0, pericenter>0, and 0<=e<1");

    eccentricity_ = eccentricity;
    semimajor_axis_ = pericenter / (1.0 - eccentricity);
    mean_motion_ = std::sqrt(total_mass / std::pow(semimajor_axis_, 3));
    direction_ = counter_rotating ? -1 : 1;
}

double FixedKeplerOrbit::period() const noexcept
{
    return 2.0 * constants::pi / mean_motion_;
}

double FixedKeplerOrbit::eccentric_anomaly(double mean_anomaly) const
{
    if (!std::isfinite(mean_anomaly) || mean_anomaly < 0.0)
        throw std::invalid_argument("Mean anomaly must be finite and nonnegative");

    const double two_pi = 2.0 * constants::pi;
    const double revolutions = std::floor(mean_anomaly / two_pi);
    const double reduced_mean = mean_anomaly - revolutions * two_pi;
    if (eccentricity_ == 0.0)
        return mean_anomaly;

    double lower = 0.0;
    double upper = two_pi;
    double E = eccentricity_ < 0.8
                   ? reduced_mean
                   : (reduced_mean < constants::pi ? 0.5 * constants::pi
                                                   : 1.5 * constants::pi);
    E = std::clamp(E, lower, upper);

    const double tolerance = 8.0 * std::numeric_limits<double>::epsilon() *
                             (1.0 + std::abs(reduced_mean));
    for (int iteration = 0; iteration < 80; ++iteration)
    {
        const double residual = E - eccentricity_ * std::sin(E) - reduced_mean;
        if (std::abs(residual) <= tolerance)
            return E + revolutions * two_pi;

        if (residual > 0.0)
            upper = E;
        else
            lower = E;

        const double derivative = 1.0 - eccentricity_ * std::cos(E);
        const double newton = E - residual / derivative;
        E = (std::isfinite(newton) && newton > lower && newton < upper)
                ? newton
                : 0.5 * (lower + upper);
    }

    throw std::runtime_error("Kepler equation did not converge");
}

void FixedKeplerOrbit::evaluate(OrbitState &state, double clock_factor) const
{
    if (!(clock_factor >= 0.0 && clock_factor <= 1.0) ||
        !std::isfinite(clock_factor))
        throw std::invalid_argument("Kepler clock factor must lie in [0,1]");

    const double E_unwrapped = eccentric_anomaly(state.phase);
    const double two_pi = 2.0 * constants::pi;
    const double revolutions = std::floor(E_unwrapped / two_pi);
    const double E = E_unwrapped - revolutions * two_pi;
    const double denominator = 1.0 - eccentricity_ * std::cos(E);
    const double root = std::sqrt((1.0 - eccentricity_) * (1.0 + eccentricity_));
    double true_anomaly = std::atan2(root * std::sin(E),
                                    std::cos(E) - eccentricity_);
    if (true_anomaly < 0.0)
        true_anomaly += two_pi;
    true_anomaly += revolutions * two_pi;

    state.R = semimajor_axis_ * denominator;
    state.theta = static_cast<double>(direction_) * true_anomaly;
    state.R_dot = clock_factor * semimajor_axis_ * mean_motion_ * eccentricity_ *
                  std::sin(E) / denominator;
    state.theta_dot = clock_factor * static_cast<double>(direction_) * mean_motion_ *
                      root / (denominator * denominator);
}

void FixedKeplerOrbit::advance_phase(OrbitState &state, double dt,
                                     double clock_factor) const
{
    if (!(dt >= 0.0) || !std::isfinite(dt))
        throw std::invalid_argument("Kepler phase increment requires finite dt>=0");
    if (!(clock_factor >= 0.0 && clock_factor <= 1.0) ||
        !std::isfinite(clock_factor))
        throw std::invalid_argument("Kepler clock factor must lie in [0,1]");
    state.phase += mean_motion_ * clock_factor * dt;
}
