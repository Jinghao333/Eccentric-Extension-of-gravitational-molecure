"""Numerical batch kernels for Xu Appendix-A one-point reconstruction.

No physical approximation is introduced here: the first routine evaluates the
same Kepler-anomaly Fourier integral with a branch-split Gauss–Legendre rule;
the second computes the same source/regular-Coulomb overlaps with a shared
adaptive vector quadrature. Both operate on original floating-point inputs.

These routines are independently unit-testable without the gravitational-atom
implementation being present.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from typing import Callable
import math

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad_vec
from scipy.interpolate import CubicSpline


@lru_cache(maxsize=16)
def _gauss_rule(n: int) -> tuple[np.ndarray, np.ndarray]:
    nodes, weights = leggauss(n)
    return nodes, weights


@dataclass(frozen=True)
class FourierBatch:
    values: np.ndarray       # shape (n_g,)
    errors: np.ndarray       # shape (n_g,); difference between consecutive orders
    evaluations: int         # orbital nodes, shared across all g
    accepted: bool
    order: int


def fourier_at_radius(
    orbit,
    tidal_multipole: Callable,
    radial_branch_structure: Callable,
    ell: int, m: int,
    g_values: tuple[int, ...],
    radius: float,
    q: float,
    *, initial_order: int = 16, max_order: int = 128,
    epsabs: float = 1e-10, epsrel: float = 1e-8,
) -> FourierBatch:
    """Batch all g at one radius, splitting exactly at orbital radius crossings.

    u_lmg(r) = (2pi)^-1 integral_0^{2pi} dE
      (1-e*cos E) u_lm[r; E] exp(i*g*(E-e*sin E)).

    The second-half anomaly is calculated directly from the orbit state,
    NEVER folded via |tan(E/2)|. The integration does not truncate harmonics
    beyond the configured g list; it evaluates those exact Fourier integrals.

    A two-order difference is a quadrature diagnostic, not a proof of accuracy.
    If it cannot meet epsabs/epsrel, the caller MUST fall back to the original
    adaptive coefficient calculation rather than accepting an unconverged value.
    """
    if initial_order < 4 or max_order < initial_order * 2:
        raise ValueError("Require max_order >= 2*initial_order >= 8")
    if not g_values:
        raise ValueError("g_values must be nonempty")
    if not (0 <= orbit.eccentricity < 1):
        raise ValueError("elliptic orbit required")
    e = orbit.eccentricity
    g_arr = np.asarray(g_values, dtype=float)
    structure = radial_branch_structure(orbit, float(radius))
    # At exact periapsis the branch finder may report E=0 as a tangency.
    # That is an endpoint, not a positive-width quadrature interval.
    interior = sorted(set(float(v) for v in structure.eccentric_anomaly_crossings
                          if 0.0 < float(v) < 2.0*math.pi))
    boundaries = (0.0, *interior, 2.0*math.pi)
    if not np.all(np.diff(boundaries) > 0.0):
        raise ValueError("invalid orbital splitting points")

    def evaluate(order: int) -> tuple[np.ndarray, int]:
        nodes, weights = _gauss_rule(order)
        accum = np.zeros(len(g_values), dtype=np.complex128)
        count = 0
        for a, b in zip(boundaries[:-1], boundaries[1:]):
            t = 0.5*(b-a)*nodes + 0.5*(b+a)
            state = orbit.state_from_eccentric_anomaly(t)
            tidal = np.asarray(tidal_multipole(state, radius, ell, m, q), dtype=np.complex128)
            if tidal.shape != t.shape:
                raise ValueError("tidal integrand must broadcast to orbital-node vector")
            anomaly = t - e*np.sin(t)
            weighted = tidal * (1.0 - e*np.cos(t)) * weights * (0.5*(b-a)/(2.0*math.pi))
            accum += weighted @ np.exp(1j * anomaly[:, None] * g_arr[None, :])
            count += order
        return accum, count

    old, n_eval = evaluate(initial_order)
    order = initial_order
    err = np.full(len(g_values), np.inf)
    while order*2 <= max_order:
        order *= 2
        new, n = evaluate(order)
        n_eval += n
        err = np.abs(new - old)
        if np.all(err <= epsabs + epsrel * np.abs(new)):
            return FourierBatch(new, err, n_eval, True, order)
        old = new
    return FourierBatch(old, err, n_eval, False, order)


@dataclass(frozen=True)
class RadialBatch:
    partial_integrals: np.ndarray
    error_norm: float
    evaluations: int
    success: bool
    message: str


def radial_overlap_batch(
    radii: np.ndarray,
    partial_sources: np.ndarray,
    inner_solution: Callable[[float], complex],
    *, epsabs: float = 1e-10, epsrel: float = 1e-8,
    limit: int = 300, breakpoints: tuple[float, ...] = (),
) -> RadialBatch:
    """Integrate all partial source overlaps in a group with ONE inner Whittaker.

    partial_sources has shape (n_partials, n_radii); each column is interpolated
    by the SAME CubicSpline scheme as original SampledRadialSource. No outgoing
    W function is evaluated: its asymptotic coefficient uses analytic W.
    """
    radii = np.asarray(radii, dtype=float)
    partials = np.asarray(partial_sources, dtype=np.complex128)
    if radii.ndim != 1 or partials.ndim != 2 or partials.shape[1] != len(radii):
        raise ValueError("expected partial_sources shape (n_partials, n_radii)")
    if partials.shape[0] == 0:
        return RadialBatch(np.empty(0, dtype=complex), 0.0, 0, True, "no sources")
    if len(radii) < 4 or np.any(np.diff(radii) <= 0):
        raise ValueError("radial grid must be strictly increasing, with >=4 samples")
    spline = CubicSpline(radii, partials.T, axis=0)

    def integrand(radius: float) -> np.ndarray:
        return complex(inner_solution(float(radius))) * np.asarray(spline(radius), dtype=np.complex128)

    points = tuple(float(x) for x in breakpoints if radii[0] < x < radii[-1])
    values, error, info = quad_vec(
        integrand, float(radii[0]), float(radii[-1]),
        epsabs=epsabs, epsrel=epsrel, limit=limit,
        points=points or None, full_output=True, norm="max",
    )
    return RadialBatch(
        np.asarray(values, dtype=np.complex128), float(error),
        int(getattr(info, "neval", -1)),
        bool(getattr(info, "success", True)),
        str(getattr(info, "message", "")),
    )
