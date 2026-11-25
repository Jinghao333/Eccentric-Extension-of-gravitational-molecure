"""Repeat one matched pilot with exactly one resolution or boundary change."""
import argparse
import json
from pathlib import Path
import re
import subprocess
import sys
import tomllib

from analyze import analyze, compare_controls


def replace_setting(text, table, key, value):
    """Edit a scratch copy, keeping all other production settings verbatim."""
    lines = text.splitlines(keepends=True)
    section = None
    for i,line in enumerate(lines):
        header = re.match(r"^\s*\[([\w]+)\]\s*(?:#.*)?$", line)
        if header:
            section = header.group(1)
        elif section == table and re.match(rf"^\s*{key}\s*=", line):
            lines[i] = f"{key} = {value}\n"
            result = "".join(lines)
            if tomllib.loads(result)[table][key] != value:
                raise ValueError(f"Failed to set {table}.{key}")
            return result
    raise ValueError(f"Expected explicit source setting {table}.{key}")


def repeat(pilot, solver, output, variation):
    base = json.loads((pilot/"D0"/"report.json").read_text())
    m = base["manifest"]
    output.mkdir(parents=True,exist_ok=False)
    original = (pilot/"source.toml").read_text(encoding="utf-8-sig")
    cfg = tomllib.loads(original)
    config, state = pilot/"source.toml", pilot/"source.bin"
    dt = m["dt"]
    if variation == "timestep":
        dt /= 2
    else:
        variant = original
        if variation in ("radial","box"):
            variant = replace_setting(variant,"grid","N",2*(m["N"]-1)+1)
        if variation == "box":
            variant = replace_setting(variant,"grid","L",2*cfg["grid"]["L"])
        if variation == "angular":
            variant = replace_setting(variant,"simulation","l_max",cfg["simulation"]["l_max"]+2)
        config, state = output/"variant.toml", output/"variant.bin"
        config.write_text(variant, encoding="utf-8")
        generator = Path(__file__).resolve().parents[2]/"initial_states"/"get_initial_state_COM_fast.py"
        subprocess.run([sys.executable,str(generator),"--config",str(config),"--output",str(state)],
                       check=True,cwd=output)
    cases = ["D0"] if variation=="angular" else ["C0","D0"]
    subprocess.run([str(solver),str(config),str(state),"--cases",",".join(cases),
                    "--radii",",".join(map(str,m["radii"])),"--periods",str(m["duration"]/m["period"]),
                    "--dt",str(dt),"--output",str(output/"controls")],check=True)
    for name in cases:
        case = output/"controls"/name
        analyze(case,plots=False)
        compare_controls(pilot/name,case,output/f"{name}_difference.json")


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pilot",type=Path)
    parser.add_argument("solver",type=Path)
    parser.add_argument("variation",choices=["radial","angular","timestep","box"])
    parser.add_argument("output",type=Path)
    args=parser.parse_args()
    repeat(args.pilot.resolve(),args.solver.resolve(),args.output.resolve(),args.variation)
