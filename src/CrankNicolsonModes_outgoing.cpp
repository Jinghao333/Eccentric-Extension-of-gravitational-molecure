/**
 * @file CrankNicolsonModes_outgoing.cpp
 * @brief Outer-boundary treatments: damping mask and complex absorbing potential.
 */
#include "CrankNicolsonModes.hpp"
#include <cmath>
#include <iostream>
#include <vector>

namespace {
void validate_layer(double factor, double W, int n, const Eigen::VectorXd &r, double L)
{
    if (!(factor >= 0 && factor < 1) || !std::isfinite(W) || W < 0 || n < 1 ||
        !r.allFinite() || r[r.size()-2] <= factor*L || !(r[r.size()-1] > r[r.size()-2]))
        throw std::invalid_argument("Outer absorber requires 0 <= L_factor < 1, finite W >= 0, n >= 1 and at least one damped interior radial node");
}
}

void CrankNicolsonModes::build_damping_mask(double factor, double W, int n)
{
    if (geometry.enabled) throw std::invalid_argument("orbital_geometry cannot use a per-step mask; use CAP or no outer absorber");
    validate_layer(factor,W,n,r,L);
    damping_mask.resize(Nr*n_modes);
    for (int i=0; i<Nr; ++i) {
        // Preserve the original nominal mask width; clip before the onset.
        const double x = std::max(0.0,(r[i]-factor*L)/((1-factor)*L));
        const double weight = W == 0 ? 1.0 : std::exp(-W*std::pow(x,n));
        for (int j=0; j<n_modes; ++j) damping_mask[j*Nr+i]=weight;
    }
    mask_settings = {factor,W,n};
    mask = true;
}

void CrankNicolsonModes::build_cap(double factor, double W, int n)
{
    validate_layer(factor,W,n,r,L);
    const double onset = factor*L;
    const double width = r[Nr-1]-onset;
    std::vector<Eigen::Triplet<dcomplex>> entries;
    for (int i=0; i<Nr; ++i) {
        if (r[i] <= onset) continue;
        const double strength = W*std::pow((r[i]-onset)/width,n);
        for (int j=0; j<n_modes; ++j) entries.emplace_back(j*Nr+i,j*Nr+i,-I*strength);
    }
    CAP.resize(Nr*n_modes,Nr*n_modes);
    CAP.setFromTriplets(entries.begin(),entries.end());
    cap_settings = {factor,W,n};
    cap = true;
}
