/**
 * @file Config.hpp
 * @brief Simulation configuration and its TOML loader.
 */
#ifndef CONFIG_HPP
#define CONFIG_HPP

#include <cmath>
#include <iostream>
#include <string>
#include <stdexcept>
#include <toml++/toml.hpp>

enum class FieldIntegrator { LegacyEndpoint, MidpointCN };

struct GeometryTimestep
{
    bool enabled = false;
    double tolerance = 0.25;
    int azimuthal_order = 2;
    int radial_order = 3;
    double ramp_tolerance = 0.02;

    bool normalized = false;
    int startup_steps = 120;
    int orbit_steps = 120;

    void validate(int lmax) const
    {
        if (normalized) {
            if (!enabled || startup_steps <= 0 || orbit_steps <= 0)
                throw std::invalid_argument("normalized_geometry requires positive integer step counts");
            return;
        }
        if (!std::isfinite(tolerance) || tolerance <= 0 ||
            !std::isfinite(ramp_tolerance) || ramp_tolerance <= 0 ||
            azimuthal_order < 0 || azimuthal_order > 2LL*lmax || radial_order < 1)
            throw std::invalid_argument("Geometry control needs positive finite tolerances, radial order >=1 and 0 <= azimuthal order <= 2*l_max");
    }
};

/**
 * @brief Run parameters read from a TOML configuration file.
 *
 * Every key is optional; defaults are those of load_config(). The tables
 * @c [mask] and @c [cap] switch on the outer damping mask and the complex
 * absorbing potential simply by being present.
 */
struct Config
{
    /// @name [grid]
    /// @{
    int N;    ///< Number of radial nodes (default 2048).
    double L; ///< Outer radius (default 5000).
    /// @}

    /// @name [simulation]
    /// @{
    double dt;                          ///< Initial timestep (default 10).
    GeometryTimestep geometry;
    FieldIntegrator integrator = FieldIntegrator::LegacyEndpoint;
    double num_periods;                 ///< Run length in initial orbital periods (may be fractional), used if num_steps = 0.
    int save_every;                     ///< Snapshot interval in steps (default 50).
    double ramp_periods;                ///< Startup ramp length in initial periods (default 0).
    bool ramp_secondary_potential = true; ///< Ramp the companion's potential with the orbital clock.
    double steps_per_period = 50.0;     ///< Post-ramp timestep is P / steps_per_period.
    double power_average_orbits = 2.0;  ///< Orbits in the rolling average of the binary power.
    bool evolve_orbit;                  ///< Feed the cloud power and accretion back on the orbit (default false).
    bool absorption;                    ///< Enable the black-hole absorbers (default false).
    int l_max;                          ///< Largest retained multipole (default 5).
    int num_steps;                      ///< Number of time samples; 0 uses num_periods.
    double epsilon;                     ///< Offset of the first radial node (default 0.001).
    double epsilon_smear;               ///< Absorber radius in units of each body's mass (default 2).
    double r_ex;                        ///< Extraction radius for fluxes (default 1000).
    /// @}

    /// @name [binary]
    /// @{
    double q;     ///< Mass ratio M2/M1 (default 0.001).
    double R;     ///< Initial separation; pericenter distance when e>0 (default 50).
    double e;     ///< Orbital eccentricity, 0 <= e < 1 (default 0).
    bool counter; ///< Counter-rotating orbit (default false).
    /// @}

    /// @name [atom]
    /// @{
    double mu;  ///< Boson mass (default 0.2).
    double a;   ///< Primary spin (default 0).
    double M_c; ///< Cloud mass (default 0).
    int n;      ///< Principal quantum number of the initial state (informational).
    int l;      ///< Orbital quantum number of the initial state (informational).
    int m;      ///< Magnetic quantum number of the initial state (informational).
    /// @}

    /// @name [parallel] and [solver]
    /// @{
    int threads = 0;                     ///< OpenMP threads; 0 uses OMP_NUM_THREADS.
    std::string preconditioner = "radial"; ///< "radial" (default) or "ilut".
    /// @}

    /// @name [mask]: multiplicative damping near the outer boundary
    /// @{
    bool mask = false;
    double mask_L_factor; ///< Fraction of L where the mask starts (required when enabled).
    double mask_W;        ///< Strength (default 0.1).
    int mask_n;           ///< Power of the profile (default 3).
    /// @}

    /// @name [cap]: complex absorbing potential near the outer boundary
    /// @{
    bool cap = false;
    double cap_L_factor; ///< Fraction of L where the CAP starts (required when enabled).
    double cap_W;        ///< Strength (default 0.1).
    int cap_n;           ///< Power of the profile (default 3).
    /// @}
};

/**
 * @brief Parse a TOML configuration file.
 * @param filename Path to the file.
 * @return The configuration, with defaults for missing keys.
 * @throws std::invalid_argument for invalid thread counts, preconditioners,
 *         steps_per_period or power_average_orbits.
 * @throws toml::parse_error if the file cannot be parsed.
 */
inline Config load_config(const std::string &filename)
{
    Config c;
    auto tbl = toml::parse_file(filename);

    c.threads = tbl["parallel"]["threads"].value_or(0);
    c.preconditioner = tbl["solver"]["preconditioner"].value_or(std::string("radial"));
    if (c.threads < 0)
        throw std::invalid_argument("parallel.threads must be nonnegative");
    if (c.preconditioner != "radial" && c.preconditioner != "ilut")
        throw std::invalid_argument("solver.preconditioner must be radial or ilut");

    c.N = tbl["grid"]["N"].value_or(2048);
    c.L = tbl["grid"]["L"].value_or(5000.0);

    c.dt = tbl["simulation"]["dt"].value_or(10.0);
    if (auto key = tbl["simulation"]["integrator"]; key)
    {
        const auto name = key.value<std::string>();
        if (name == "midpoint_cn")
            c.integrator = FieldIntegrator::MidpointCN;
        else if (name != "legacy_endpoint")
            throw std::invalid_argument("simulation.integrator must be legacy_endpoint or midpoint_cn");
    }
    c.num_periods = tbl["simulation"]["num_periods"].value_or(0.0);
    c.num_steps = tbl["simulation"]["num_steps"].value_or(0);
    c.save_every = tbl["simulation"]["save_every"].value_or(50);
    c.l_max = tbl["simulation"]["l_max"].value_or(5);
    c.r_ex = tbl["simulation"]["r_ex"].value_or(1000.0);
    c.ramp_periods = tbl["simulation"]["ramp_periods"].value_or(0.0);
    c.ramp_secondary_potential = tbl["simulation"]["ramp_secondary_potential"].value_or(true);
    c.steps_per_period = tbl["simulation"]["steps_per_period"].value_or(50.0);
    if (!(c.steps_per_period > 0.0))
        throw std::invalid_argument("simulation.steps_per_period must be positive");
    c.power_average_orbits = tbl["simulation"]["power_average_orbits"].value_or(2.0);
    if (!(c.power_average_orbits > 0.0))
        throw std::invalid_argument("simulation.power_average_orbits must be positive");
    c.epsilon = tbl["simulation"]["epsilon"].value_or(0.001);
    c.epsilon_smear = tbl["simulation"]["epsilon_smear"].value_or(2.0);
    c.evolve_orbit = tbl["simulation"]["evolve_orbit"].value_or(false);
    if (c.integrator == FieldIntegrator::MidpointCN && c.evolve_orbit)
        throw std::invalid_argument("midpoint_cn requires a prescribed binary with evolve_orbit=false");
    c.absorption = tbl["simulation"]["absorption"].value_or(false);

    c.q = tbl["binary"]["q"].value_or(0.001);
    c.R = tbl["binary"]["R"].value_or(50.0);
    c.e = tbl["binary"]["e"].value_or(0.0);
    c.counter = tbl["binary"]["counter"].value_or(false);
    if (!(c.q >= 0.0) || !(c.R > 0.0) || !(c.e >= 0.0 && c.e < 1.0) ||
        !std::isfinite(c.q) || !std::isfinite(c.R) || !std::isfinite(c.e))
        throw std::invalid_argument("binary requires finite q>=0, R>0, and 0<=e<1");
    if (c.e > 0.0 && c.evolve_orbit)
        throw std::invalid_argument(
            "Eccentric inspiral is not implemented: require binary.e=0 when simulation.evolve_orbit=true");

    c.mu = tbl["atom"]["mu"].value_or(0.2);
    c.n = tbl["atom"]["n"].value_or(2);
    c.l = tbl["atom"]["l"].value_or(1);
    c.m = tbl["atom"]["m"].value_or(1);
    c.a = tbl["atom"]["a"].value_or(0.0);
    c.M_c = tbl["atom"]["M_c"].value_or(0.0);

    auto read_layer = [&](const char *name, bool &enabled, double &factor, double &strength, int &exponent) {
        const auto layer = tbl[name];
        if (!layer) return;
        if (!layer.is_table()) throw std::invalid_argument(std::string(name)+" must be a table");
        enabled = true;
        factor = layer["L_factor"].value<double>().value_or(-1.);
        strength = layer["W"] ? layer["W"].value<double>().value_or(-1.) : 0.1;
        exponent = layer["n"] ? layer["n"].value<int>().value_or(0) : 3;
        if (!(factor >= 0 && factor < 1) || !std::isfinite(strength) || strength < 0 ||
            exponent < 1 || c.N < 3 || !std::isfinite(c.L) || !(c.L>0) ||
            !(c.epsilon+(c.N-2)*(c.L/(c.N-1)) > factor*c.L))
            throw std::invalid_argument(std::string(name)+" requires a resolved layer with explicit 0<=L_factor<1, W>=0 and integer n>=1");
    };
    read_layer("mask",c.mask,c.mask_L_factor,c.mask_W,c.mask_n);
    read_layer("cap",c.cap,c.cap_L_factor,c.cap_W,c.cap_n);

    auto sim = tbl["simulation"];
    auto setting = [&](const char *key, auto fallback) {
        if (!sim[key]) return fallback;
        const auto v = sim[key].value<decltype(fallback)>();
        if (!v) throw std::invalid_argument(std::string("Invalid simulation.")+key);
        return *v;
    };
    const auto rule = setting("timestep_rule", std::string("fixed"));
    if (rule != "fixed" && rule != "orbital_geometry" && rule != "normalized_geometry")
        throw std::invalid_argument("simulation.timestep_rule must be fixed, orbital_geometry or normalized_geometry");
    auto &g = c.geometry;
    g.normalized = rule == "normalized_geometry";
    g.enabled = rule != "fixed";
    if (g.normalized) {
        for (auto key : {"startup_steps_per_period", "orbit_steps_per_period"})
            if (sim[key] && !sim[key].is_integer())
                throw std::invalid_argument(std::string("Expected a positive integer for simulation.") + key);
        g.startup_steps = setting("startup_steps_per_period", 120);
        g.orbit_steps = setting("orbit_steps_per_period", 120);
    }
    g.tolerance = setting("geometry_tolerance", g.tolerance);
    g.ramp_tolerance = setting("ramp_weight_tolerance", g.ramp_tolerance);
    g.azimuthal_order = setting("effective_azimuthal_order", c.l_max == 0 ? 0 : 2);
    g.radial_order = setting("effective_radial_order", g.radial_order);
    g.validate(c.l_max);
    if (g.enabled) {
        c.dt = g.normalized ? 1.0 : setting("dt",50.0); // Initialization placeholder until first proposal.
        if (c.mask) throw std::invalid_argument("Variable timesteps cannot use a per-step mask; use CAP or no outer absorber");
        if (c.integrator != FieldIntegrator::MidpointCN || c.evolve_orbit ||
            setting("num_steps",0) != 0 || !std::isfinite(c.num_periods) || c.num_periods <= 0 ||
            !std::isfinite(c.dt) || c.dt <= 0 || !std::isfinite(c.ramp_periods) || c.ramp_periods < 0)
            throw std::invalid_argument("Variable timesteps require prescribed midpoint_cn, finite positive dt and num_periods, num_steps=0 and nonnegative ramp_periods");
    }
    return c;
}

#endif
