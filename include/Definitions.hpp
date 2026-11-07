/**
 * @file Definitions.hpp
 * @brief Common types: complex scalars and angular modes.
 */
#pragma once
#include <complex>

namespace constants
{
/// Pi, provided here because the non-standard M_PI macro is not portable.
inline constexpr double pi = 3.141592653589793238462643383279502884;
}

/// Double-precision complex number.
using dcomplex = std::complex<double>;
/// The imaginary unit.
extern const dcomplex I;

/// A spherical-harmonic mode (l, m).
struct mode
{
    int l;
    int m;

    bool operator==(const mode &other) const
    {
        return l == other.l && m == other.m;
    }
};
