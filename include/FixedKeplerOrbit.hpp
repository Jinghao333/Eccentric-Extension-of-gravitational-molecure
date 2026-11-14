/**
 * @file FixedKeplerOrbit.hpp
 * @brief Prescribed bound Kepler trajectory for a fixed binary.
 */
#pragma once

#include "OrbitState.hpp"

/**
 * @brief Fixed planar Kepler ellipse parameterized by unwrapped mean anomaly.
 *
 * This object owns only immutable trajectory parameters.  The instantaneous
 * physical state is always written to an external OrbitState.
 */
class FixedKeplerOrbit
{
public:
    FixedKeplerOrbit() = default;
    FixedKeplerOrbit(double total_mass, double pericenter, double eccentricity,
                     bool counter_rotating);

    void configure(double total_mass, double pericenter, double eccentricity,
                   bool counter_rotating);

    [[nodiscard]] double eccentricity() const noexcept { return eccentricity_; }
    [[nodiscard]] double semimajor_axis() const noexcept { return semimajor_axis_; }
    [[nodiscard]] double mean_motion() const noexcept { return mean_motion_; }
    [[nodiscard]] double period() const noexcept;
    [[nodiscard]] int direction() const noexcept { return direction_; }

    /** Solve E - e sin(E) = M, preserving complete revolutions of M. */
    [[nodiscard]] double eccentric_anomaly(double mean_anomaly) const;

    /** Evaluate geometry and velocities at state.phase (the mean anomaly). */
    void evaluate(OrbitState &state, double clock_factor = 1.0) const;

    /** Advance mean anomaly with the supplied previous-step clock factor. */
    void advance_phase(OrbitState &state, double dt, double clock_factor) const;

private:
    double eccentricity_ = 0.0;
    double semimajor_axis_ = 1.0;
    double mean_motion_ = 1.0;
    int direction_ = 1;
};
