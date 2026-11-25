#!/usr/bin/env python3
"""Compare orbit-averaged ionization power from two simulation folders.

Usage:
    python dashboard/compare_ionization.py SIMULATION_1 SIMULATION_2 START END

START and END are in units of each run's initial Kepler period.  The reported
power uses the dashboard's signed convention: outgoing ionization is negative.
"""

import argparse
import math
from pathlib import Path
import sys
import tomllib

import numpy as np


REQUESTED_MODES = ((1, 1), (3, 3), (2, 2), (4, 4))


def read_complete_table(path: Path, columns: int, limit: int | None = None) -> np.ndarray:
    """Read complete whitespace-delimited rows, ignoring a live partial tail."""
    rows = []
    with path.open("r", encoding="utf8") as stream:
        for line_number, line in enumerate(stream, 1):
            if limit is not None and len(rows) >= limit:
                break
            values = np.fromstring(line, sep=" ")
            if values.size == 0 and not line.strip():
                continue
            if values.size != columns:
                if line.endswith("\n"):
                    raise ValueError(
                        f"{path}:{line_number} has {values.size} columns; expected {columns}"
                    )
                break
            if not np.isfinite(values).all():
                raise ValueError(f"{path}:{line_number} contains a non-finite value")
            rows.append(values)
    if not rows:
        raise ValueError(f"{path} contains no complete rows")
    return np.asarray(rows)


def initial_period(config: dict) -> float:
    binary = config["binary"]
    total_mass = 1.0 + float(binary["q"])
    semimajor_axis = float(binary["R"]) / (1.0 - float(binary.get("e", 0.0)))
    return 2.0 * math.pi * math.sqrt(semimajor_axis**3 / total_mass)


def time_average(times: np.ndarray, values: np.ndarray, start: float, end: float) -> np.ndarray:
    if times[-1] < end:
        raise ValueError(
            f"requested endpoint {end:.15g} exceeds the last complete diagnostic time "
            f"{times[-1]:.15g}"
        )
    inside = (times > start) & (times < end)
    grid = np.r_[start, times[inside], end]
    sampled = np.column_stack(
        [np.interp(grid, times, values[:, column]) for column in range(values.shape[1])]
    )
    return np.trapezoid(sampled, grid, axis=0) / (end - start)


def load_average(folder: Path, orbit_start: float, orbit_end: float) -> dict:
    folder = folder.resolve()
    with (folder / "config.toml").open("rb") as stream:
        config = tomllib.load(stream)
    modes = np.loadtxt(folder / "modes.txt", delimiter=",", dtype=int, ndmin=2)
    number_of_modes = len(modes)
    period = initial_period(config)
    start = orbit_start * period
    end = orbit_end * period

    times = read_complete_table(folder / "t.txt", 1)[:, 0]
    if not np.all(np.diff(times) > 0):
        raise ValueError(f"{folder / 't.txt'} is not strictly increasing")
    required_rows = int(np.searchsorted(times, end, side="left")) + 1
    if required_rows > len(times):
        raise ValueError(
            f"{folder.name} reaches only t/P={times[-1] / period:.9g}; "
            f"the requested interval ends at {orbit_end:g}"
        )

    def table(name: str, columns: int) -> np.ndarray:
        values = read_complete_table(folder / f"{name}.txt", columns, required_rows)
        if len(values) < required_rows:
            raise ValueError(
                f"{folder / (name + '.txt')} has only {len(values)} complete rows; "
                f"{required_rows} are required"
            )
        return values[:required_rows]

    times = times[:required_rows]
    norms = table("norms", number_of_modes)
    energy = table("energy", 1)[:, 0]
    flux = table("flux", number_of_modes)
    energy_flux = table("energy_flux", number_of_modes)

    mu = float(config["atom"]["mu"])
    q = float(config["binary"]["q"])
    if mu <= 0 or q == 0:
        raise ValueError(f"{folder.name} requires mu > 0 and q != 0 for this normalization")
    particle_number = norms.sum(axis=1)
    if np.any(particle_number <= 0):
        raise ValueError(f"{folder.name} has non-positive cloud particle number")

    epsilon = energy / particle_number
    power_by_mode = (
        energy_flux + epsilon[:, None] * flux / mu
    ) / (mu * particle_number[:, None] * q * q)
    averaged_modes = time_average(times, power_by_mode, start, end)
    mode_index = {tuple(mode): index for index, mode in enumerate(modes.tolist())}
    total = float(averaged_modes.sum())
    mode11 = float(averaged_modes[mode_index[(1, 1)]]) if (1, 1) in mode_index else None

    contributions = {
        "Total": total,
        "(1,1)": mode11,
        "Remaining (all except 1,1)": total - mode11 if mode11 is not None else total,
    }
    for mode in REQUESTED_MODES[1:]:
        contributions[f"({mode[0]},{mode[1]})"] = (
            float(averaged_modes[mode_index[mode]]) if mode in mode_index else None
        )
    return {
        "folder": folder,
        "period": period,
        "time_start": start,
        "time_end": end,
        "contributions": contributions,
    }


def number(value: float | None) -> str:
    return "N/A" if value is None else f"{value:.9e}"


def percentage(value: float | None, total: float) -> str:
    return "N/A" if value is None or total == 0 else f"{100.0 * value / total:.6f}%"


def magnitude_change(first: float | None, second: float | None) -> str:
    if first is None or second is None or first == 0:
        return "N/A"
    return f"{100.0 * (abs(second) / abs(first) - 1.0):+.6f}%"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("simulation_1", type=Path)
    parser.add_argument("simulation_2", type=Path)
    parser.add_argument("orbit_start", type=float)
    parser.add_argument("orbit_end", type=float)
    args = parser.parse_args()
    if not math.isfinite(args.orbit_start) or not math.isfinite(args.orbit_end):
        parser.error("orbit bounds must be finite")
    if args.orbit_start < 0 or args.orbit_end <= args.orbit_start:
        parser.error("require 0 <= orbit_start < orbit_end")

    try:
        first = load_average(args.simulation_1, args.orbit_start, args.orbit_end)
        second = load_average(args.simulation_2, args.orbit_start, args.orbit_end)
    except (OSError, ValueError, KeyError, tomllib.TOMLDecodeError) as error:
        parser.exit(2, f"error: {error}\n")

    print(f"Interval: {args.orbit_start:g} <= t/P0 <= {args.orbit_end:g}")
    print(f"Simulation 1: {first['folder']}  (P0={first['period']:.12g})")
    print(f"Simulation 2: {second['folder']}  (P0={second['period']:.12g})")
    if not math.isclose(first["period"], second["period"], rel_tol=1e-12, abs_tol=1e-12):
        print("Note: periods differ; each interval uses its own run's P0.")
    print("Signed convention: outgoing ionization power is negative.\n")
    print("| Contribution | Simulation 1 | % of total 1 | Simulation 2 | % of total 2 | Magnitude change 2 vs 1 |")
    print("|---|---:|---:|---:|---:|---:|")
    total1 = first["contributions"]["Total"]
    total2 = second["contributions"]["Total"]
    for label, value1 in first["contributions"].items():
        value2 = second["contributions"].get(label)
        print(
            f"| {label} | {number(value1)} | {percentage(value1, total1)} | "
            f"{number(value2)} | {percentage(value2, total2)} | "
            f"{magnitude_change(value1, value2)} |"
        )
    print("\nRemaining includes every evolved mode except (1,1), so the explicitly listed")
    print("(2,2), (3,3), and (4,4) rows are subsets of Remaining, not additional terms.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
