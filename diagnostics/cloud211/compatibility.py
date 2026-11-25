"""Record or verify production CLI runs in a fresh scratch directory."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import struct
import subprocess

import numpy as np


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def files(folder):
    return {p.relative_to(folder).as_posix():sha(p) for p in folder.rglob("*") if p.is_file()}


def execute(exe, cwd, *args):
    result = subprocess.run([str(exe), *map(str,args)], cwd=cwd, capture_output=True, text=True)
    if result.returncode:
        raise RuntimeError(f"{result.stdout}\n{result.stderr}")
    return result.stdout


def capture(exe, root, configs, reference):
    root.mkdir(parents=True, exist_ok=False)
    summary = {"configs":{}, "runs":{}}
    for config in configs:
        text = execute(exe, root, "--check-config", config)
        summary["configs"][str(config)] = [sha(config),
            [line for line in text.splitlines() if not line.startswith("Executable:")]]
    nr = 65
    r = np.arange(nr)*2
    u = r*r*np.exp(-r/15)
    u[[0,-1]]=0
    u /= np.linalg.norm(u)*np.sqrt(2)
    data = struct.pack("<8siqqii", b"MODESBIN",1,nr,1,1,1)+u.astype("<f8").tobytes()+bytes(8*nr)
    for model in ("fixed_circular","fixed_eccentric","circular_inspiral"):
        directory=root/model
        directory.mkdir()
        ecc = model=="fixed_eccentric"
        clock = ('integrator="midpoint_cn"\ntimestep_rule="normalized_geometry"\n'
                 'startup_steps_per_period=16\norbit_steps_per_period=16\nnum_periods=0.5'
                 if ecc else 'dt=2.0\nnum_steps=9')
        config = f"""[grid]
N=65
L=128.0
[binary]
q=0.1
R=16.0
e={0.4 if ecc else 0.0}
[atom]
n=2
l=1
m=1
mu=0.2
M_c=0.01
[simulation]
{clock}
save_every=4
r_ex=60.0
l_max=2
epsilon=0.0
evolve_orbit={str(model=="circular_inspiral").lower()}
[parallel]
threads=1
"""
        (directory/"config.toml").write_text(config)
        (directory/"initial.bin").write_bytes(data)
        execute(exe,directory,"config.toml","initial.bin")
        folder = next(directory.glob("Simulation_*"))
        summary["runs"][model]=files(folder)
        if reference:
            checkpoint = next((reference/model).glob("Simulation_*"))
        else:
            checkpoint = folder
        resumed=directory/"continued"
        shutil.copytree(checkpoint,resumed)
        cfg=resumed/"config.toml"
        cfg.write_text(config.replace("num_steps=9","num_steps=17").replace("num_periods=0.5","num_periods=1.0"))
        execute(exe,directory,resumed)
        summary["runs"][model+"_resume"]=files(resumed)
    (root/"compatibility.json").write_text(json.dumps(summary,indent=2))
    if reference and summary != json.loads((reference/"compatibility.json").read_text()):
        raise RuntimeError("Production comparison differs; inspect compatibility.json files")
    print(f"{len(configs)} configuration checks, three runs and three checkpoint continuations completed")
    if reference:
        print("All fields, time grids, orbit histories, diagnostics and checkpoints byte-identical")


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("solver",type=Path)
    p.add_argument("output",type=Path)
    p.add_argument("--reference",type=Path)
    p.add_argument("configs",nargs="*",type=Path)
    args=p.parse_intermixed_args()
    capture(args.solver.resolve(),args.output.resolve(),[c.resolve() for c in args.configs],
            args.reference.resolve() if args.reference else None)
