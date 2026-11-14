/**
 * @file OrbitState.hpp
 * @brief Minimal instantaneous state of the binary orbit.
 */
#pragma once

#include "Definitions.hpp"

#include <cmath>

/**
 * @brief Authoritative physical state of an ordinary planar binary orbit.
 *
 * Units are G = c = 1.  For the ordinary binary, @c R is the separation,
 * @c theta is the azimuth of the primary, and the primary and secondary lie
 * at +R1 and -R2 respectively along that direction.  The signed angular
 * velocity is @c theta_dot, while @c phase is the separately accumulated,
 * non-decreasing phase used by the rolling power average.
 *
 * The displaced isolated-black-hole configuration is a compatibility case:
 * it stores the source displacement in @c R, but must not use the COM-radius
 * getters to locate that source.
 */
struct OrbitState
{
    double M1 = 1.0;       ///< Primary black-hole mass.
    double M2 = 0.0;       ///< Secondary black-hole mass.
    double R = 0.0;        ///< Binary separation (isolated case: source displacement magnitude).
    double theta = 0.0;    ///< Azimuth of the primary in the binary COM convention.
    double R_dot = 0.0;    ///< Time derivative of the binary separation.
    double theta_dot = 0.0; ///< Signed angular velocity of the primary azimuth.
    double phase = 0.0;    ///< Non-decreasing accumulated phase for rolling averages.

    [[nodiscard]] double total_mass() const noexcept { return M1 + M2; }

    /** @brief Mass ratio q = M2/M1.  The physical domain requires M1 > 0. */
    [[nodiscard]] double q() const noexcept { return M2 / M1; }

    /** @brief Primary distance R1 = M2 R / (M1 + M2) from the COM. */
    [[nodiscard]] double R1() const noexcept
    {
        return M2 / total_mass() * R;
    }

    /** @brief Secondary distance R2 = M1 R / (M1 + M2) from the COM. */
    [[nodiscard]] double R2() const noexcept
    {
        return M1 / total_mass() * R;
    }

    /** @brief Unsigned angular speed corresponding to the legacy Omega value. */
    [[nodiscard]] double angular_speed() const noexcept { return std::abs(theta_dot); }

    /** @brief Circular Kepler frequency sqrt((M1 + M2) / R^3). */
    [[nodiscard]] double circular_kepler_frequency() const noexcept
    {
        return std::sqrt(total_mass() / std::pow(R, 3));
    }

    /** @brief Circular Kepler period 2 pi / Omega_K. */
    [[nodiscard]] double circular_kepler_period() const noexcept
    {
        return 2.0 * constants::pi / circular_kepler_frequency();
    }
};
