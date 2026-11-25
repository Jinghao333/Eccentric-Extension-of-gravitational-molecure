"""Analyze one completed companion-off case, or compare production archives."""
import argparse
import csv
import hashlib
import json
import struct
import tomllib
from pathlib import Path

import numpy as np
from spectrum import characterize


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def write_json(path, obj):
    path.write_text(json.dumps(obj, indent=2, allow_nan=False)+"\n", encoding="utf-8")


def read_samples(path):
    return np.genfromtxt(path, delimiter=",", names=True, encoding="utf-8", dtype=None, ndmin=1)


def plot_flux(case,mode):
    """View any stored mode with independent scales at each extraction radius."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    report=json.loads((case/"report.json").read_text())
    meta=report["manifest"]
    if list(mode) not in meta["modes"]:
        raise ValueError(f"Mode {mode} is absent")
    flux=read_samples(case/"flux.csv")
    fig,axes=plt.subplots(len(report["surfaces"]),2,squeeze=False,figsize=(11,3*len(report["surfaces"])),
                          sharex=True,layout="constrained")
    for row,surface in zip(axes,report["surfaces"],strict=True):
        data=flux[flux["r_node"]==surface["radius"]]
        total=data[data["mode"]=="total"]
        selected=data[(data["l"]==mode[0]) & (data["m"]==mode[1])]
        for ax,field,quantity in zip(row,("JQ","JE"),("Number flux","Energy flux"),strict=True):
            ax.plot(total["t"]/meta["period"],total[field],label="Total")
            ax.plot(selected["t"]/meta["period"],selected[field],"--",label=f"({mode[0]}, {mode[1]})")
            ax.set_ylabel(quantity)
            title=f'r = {surface["radius"]:.6g}'
            if surface.get("transport_above_balance_floor") is False:
                title+="   |   Transport below continuity floor"
            ax.set_title(title,fontsize=10)
            ax.legend(fontsize=8)
    for ax in axes[-1]:
        ax.set_xlabel("Time / original period")
    fig.savefig(case/f"flux_by_radius_{mode[0]}_{mode[1]}.png",dpi=140)
    plt.close(fig)


def write_modes(path, modes, state):
    with path.open("xb") as out:
        out.write(struct.pack("<8siqq", b"MODESBIN", 1, state.shape[1], len(modes)))
        for (l, m), row in zip(modes, state, strict=True):
            out.write(struct.pack("<ii", l, m))
            out.write(row.real.astype("<f8").tobytes())
            out.write(row.imag.astype("<f8").tobytes())


def input_labels(path):
    with path.open("rb") as stream:
        magic,version,nr,count = struct.unpack("<8siqq",stream.read(28))
        if magic!=b"MODESBIN" or version!=1:
            raise ValueError("Unknown input format")
        labels=[]
        for _ in range(count):
            labels.append(list(struct.unpack("<ii",stream.read(8))))
            stream.seek(16*nr,1)
    return labels


def read_production_fluxes(folder, until=None):
    cfg = tomllib.loads((folder/"config.toml").read_text(encoding="utf-8-sig"))
    modes = np.loadtxt(folder/"modes.txt", delimiter=",", dtype=int, ndmin=2)
    rows = None
    if until is None:
        t = np.loadtxt(folder/"t.txt", ndmin=1)
    else:
        # A live archive can be appending new rows. Only read the completed
        # prefix needed for the requested time window.
        times=[]
        with (folder/"t.txt").open() as stream:
            for line in stream:
                value=float(line)
                if value>until+8*np.spacing(until):
                    break
                times.append(value)
        rows=len(times)
        if not rows:
            raise ValueError("No production samples in the requested interval")
        t=np.asarray(times)
    mu = cfg.get("atom",{}).get("mu",0.2)
    number = np.loadtxt(folder/"flux.txt", ndmin=2, max_rows=rows)/mu
    energy = -np.loadtxt(folder/"energy_flux.txt", ndmin=2, max_rows=rows)
    if (number.shape != (len(t), len(modes)) or energy.shape != number.shape
            or not np.all(np.diff(t)>0) or not np.isfinite(number).all()
            or not np.isfinite(energy).all() or not np.isfinite(t).all()):
        raise ValueError("Invalid or incomplete production flux history")
    return t, modes, number, energy


def compare_archive(control, archive, out):
    report=json.loads((control/"report.json").read_text())
    meta=report["manifest"]
    if meta["case"]!="D0" or digest(archive/"config.toml")!=meta["source_config_sha256"]:
        raise ValueError("Archive config must match the D0 source")
    cfg=tomllib.loads((archive/"config.toml").read_text(encoding="utf-8-sig"))
    simulation=cfg.get("simulation",{})
    expected=simulation.get("integrator","legacy_endpoint")
    marker=archive/"integrator.txt"
    actual=marker.read_text().strip() if marker.exists() else "legacy_endpoint"
    if actual!=expected:
        raise ValueError("Archive integrator marker disagrees with its config")
    archived=np.loadtxt(archive/"psi"/"snapshot_000000.txt",ndmin=2)
    full=np.fromfile(control/"state_0.c128",dtype="<c16")
    if archived.shape!=(len(full),2):
        raise ValueError("Archive initial snapshot has a different dimension")
    rel_error=np.linalg.norm(archived[:,0]+1j*archived[:,1]-full)/np.linalg.norm(full)
    if not rel_error <= 5e-6:
        raise ValueError("Archive initial snapshot disagrees with D0 beyond its text precision")
    t,lm,number,energy=read_production_fluxes(archive,meta["duration"])
    if not np.array_equal(lm,meta["modes"]):
        raise ValueError("Archive mode order differs from D0")
    k=int(np.floor(simulation.get("r_ex",1000.)/meta["dr"]))
    radius=meta["epsilon"]+k*meta["dr"]
    data=read_samples(control/"flux.csv")
    data=data[np.isclose(data["r_node"],radius,rtol=0,atol=4*np.spacing(radius))]
    total=data[data["mode"]=="total"]
    if not len(total):
        raise ValueError("D0 did not extract at the archived solver node")
    selected=data[data["mode"]=="mode"]
    aligned=[]
    for field in ("JQ","JE"):
        values=selected[field].reshape(len(total),len(lm))
        aligned.append(np.column_stack([np.interp(t,total["t"],column) for column in values.T]))
    with out.open("x",newline="") as stream:
        writer=csv.writer(stream)
        writer.writerow(["t","r_node","mode","l","m","JQ_D0","JQ_B0","JE_D0","JE_B0"])
        for i,time in enumerate(t):
            for j,(l,m) in enumerate(lm):
                writer.writerow([time,radius,"mode",l,m,aligned[0][i,j],number[i,j],aligned[1][i,j],energy[i,j]])
            writer.writerow([time,radius,"total","","",aligned[0][i].sum(),number[i].sum(),aligned[1][i].sum(),energy[i].sum()])
    return float(rel_error)


def compare_binary(b0, b1, out):
    if digest(b0/"config.toml") != digest(b1/"config.toml"):
        raise ValueError("B0/B1 require identical production configs")
    for name in ("integrator.txt", "dt_step.txt", "R_orb.txt", "theta_orb.txt",
                 "phase_orb.txt", "Omega_orb.txt"):
        left, right = b0/name, b1/name
        if left.exists() != right.exists() or (left.exists() and digest(left) != digest(right)):
            raise ValueError(f"Binary evolution mismatch in {name}")
    t, modes, q0, e0 = read_production_fluxes(b0)
    other_t, other_modes, q1, e1 = read_production_fluxes(b1)
    if not np.array_equal(t, other_t) or not np.array_equal(modes, other_modes):
        raise ValueError("Different sample times or mode order")
    with out.open("x", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["t","mode","l","m","JQ_B0","JQ_B1","JE_B0","JE_B1"])
        for i, time in enumerate(t):
            for j,(l,m) in enumerate(modes):
                writer.writerow([time,"mode",l,m,q0[i,j],q1[i,j],e0[i,j],e1[i,j]])
            writer.writerow([time,"total","","",q0[i].sum(),q1[i].sum(),e0[i].sum(),e1[i].sum()])


def compare_controls(first, second, out):
    reports = [json.loads((case/"report.json").read_text()) for case in (first,second)]
    a,b = reports
    for key in ("mu","M_c","source_x","period"):
        if a["manifest"][key] != b["manifest"][key]:
            raise ValueError(f"Different source physics: {key}")
    if not np.isclose(a["manifest"]["duration"], b["manifest"]["duration"],rtol=1e-13,atol=0):
        raise ValueError("Compare runs at the same elapsed time")
    def spectrum_value(report, key):
        return report["spectrum"].get(key)
    result = dict(first=str(first), second=str(second), differences={})
    for key in ("mean_energy","variance","rms_residual","analytical_residual"):
        result["differences"][key]=b["initial"][key]-a["initial"][key]
    for key in ("positive_estimate","cluster_weight","fidelity"):
        values=[spectrum_value(report,key) for report in reports]
        if None not in values:
            result["differences"][key]=values[1]-values[0]
    if a["manifest"]["radii"] != b["manifest"]["radii"]:
        raise ValueError("Different requested extraction radii")
    histories = [read_samples(case/"flux.csv") for case in (first,second)]
    result["surfaces"] = []
    for left,right in zip(a["surfaces"],b["surfaces"],strict=True):
        surface = dict(first_node=left["radius"],second_node=right["radius"],
            differences={key:right[key]-left[key] for key in
                         ("peak_outward","peak_time","last_flux","net_face_fraction","net_node_fraction")},
            max_balance_fraction=max(left["max_balance_fraction"],right["max_balance_fraction"]))
        for label in ("total","mode_1_1"):
            series=[]
            for data,row in zip(histories,(left,right)):
                selection=(data["mode"]=="total") if label=="total" else ((data["l"]==1)&(data["m"]==1))
                series.append(data[selection & (data["r_node"]==row["radius"])])
            x,y=series
            if label=="total":
                surface["transport_resolved_in_both"]=all(
                    max(abs(data["deltaQ_face"]))/report["initial"]["Q0"] >
                    max(row["max_balance_fraction"],16*np.finfo(float).eps)
                    for data,report,row in zip(series,reports,(left,right),strict=True))
            metrics={}
            for key in ("JQ","JE","deltaQ_face"):
                aligned=np.interp(x["t"],y["t"],y[key])
                difference=float(np.max(np.abs(aligned-x[key])))
                scale=float(np.max(np.abs(x[key])))
                metrics[key]=dict(max_abs_difference=difference,reference_peak_abs=scale,
                                  difference_over_peak=difference/scale if scale and surface["transport_resolved_in_both"] else None)
            surface[label]=metrics
        result["surfaces"].append(surface)
    with out.open("x") as stream:
        json.dump(result,stream,indent=2,allow_nan=False)


def analyze(case, corrected=None, plots=True, mode=(1,1)):
    if not (case/"completed.txt").is_file():
        raise ValueError("Case did not finish")
    meta = json.loads((case/"manifest.json").read_text())
    if list(mode) not in meta["modes"]:
        raise ValueError(f"Requested mode {mode} is absent from the solver basis")
    initial = json.loads((case/"initial.json").read_text())
    manifest = dict(meta, source_config_sha256=digest(case.parent/"source.toml"),
                    source_input_sha256=digest(case.parent/"source.bin"))
    labels = input_labels(case.parent/"source.bin")
    manifest["input_modes"] = labels
    manifest["retained_modes_zero_filled_by_reader"] = [lm for lm in meta["modes"] if lm not in labels]
    shape = (len(meta["modes"]), meta["N"])
    state = np.fromfile(case/"state_0.c128", dtype="<c16").reshape(shape)
    spectral, candidate = characterize(case, meta, state, initial["E211"])
    if corrected is not None:
        if meta["case"] != "D0" or candidate is None:
            raise ValueError("A resolved D0 reference is required for a corrected input")
        write_modes(corrected, meta["modes"], candidate)
    flux = read_samples(case/"flux.csv")
    balance = read_samples(case/"balance.csv")
    manifest["extraction"] = [dict(requested=float(row["r_requested"]),
                                 node=float(row["r_node"]),face=float(row["r_face"]))
                              for row in flux[(flux["step"]==0) & (flux["mode"]=="total")]]
    write_json(case/"manifest.json",manifest)
    summaries = []
    for radius in np.unique(balance["r_node"]):
        b = balance[balance["r_node"] == radius]
        f = flux[(flux["r_node"] == radius) & (flux["mode"] == "total")]
        peak = int(np.argmax(f["JQ"]))
        threshold = max(0., f["JQ"][peak])*0.05
        balance_floor = float(max(abs(b["continuity"]))/initial["Q0"])
        transport = float(max(abs(b["deltaQ_face"]))/initial["Q0"])
        resolved = transport > max(balance_floor,16*np.finfo(float).eps)
        above = np.flatnonzero(f["JQ"] > threshold) if resolved else []
        summaries.append(dict(radius=float(radius), peak_outward=float(f["JQ"][peak]),
                              peak_time=float(f["t"][peak]),
                              max_abs_face_fraction=transport,transport_above_balance_floor=bool(resolved),
                              first_above_5pct_of_peak=float(f["t"][above[0]]) if len(above) else None,
                              last_above_5pct_of_peak=float(f["t"][above[-1]]) if len(above) else None,
                              last_flux=float(f["JQ"][-1]), net_node_fraction=float(f["deltaQ_node"][-1]/initial["Q0"]),
                              net_face_fraction=float(b["deltaQ_face"][-1]/initial["Q0"]),
                              max_balance_fraction=balance_floor))
    snapshots = sorted(case.glob("state_*.c128"), key=lambda p:int(p.stem.split("_")[1]))
    profiles, times, populations, departures = [], [], [], []
    # The energy-expectation phase provides a single-frequency stationarity test.
    # Use the actual CN phase including the shortened final interval.
    sample_times = balance[balance["r_node"]==balance["r_node"][0]]["t"]
    phase = np.r_[0., np.cumsum(-2*np.arctan(np.diff(sample_times)*initial["mean_energy"]/2))]
    for path in snapshots:
        step = int(path.stem.split("_")[1])
        u = np.fromfile(path, dtype="<c16").reshape(shape)
        profiles.append(np.sum(abs(u)**2-abs(state)**2, axis=0))
        populations.append((meta["dr"]*np.sum(abs(u)**2, axis=1)).tolist())
        times.append(float(sample_times[step]))
        departures.append(float(np.linalg.norm(u-state*np.exp(1j*phase[step]))/np.linalg.norm(state)))
    report = dict(manifest=manifest, initial=initial, spectrum=spectral, surfaces=summaries,
                  max_norm_drift=float(max(abs(balance["Q_full"]-initial["Q0"]))),
                  snapshot_times=times, snapshot_populations=populations,
                  relative_departure_from_single_CN_phase=departures)
    if meta["epsilon"] != 0:
        report["grid_provenance"] = ("MODESBIN has no radii. Verify generator coordinates: "
            "the COM script offsets only the first node, the solver offsets all nodes.")
    write_json(case/"report.json", report)
    if plots:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(3, 2, figsize=(12,10), layout="constrained")
        for s in summaries:
            data = flux[flux["r_node"]==s["radius"]]
            total = data[data["mode"]=="total"]
            main = data[(data["l"]==mode[0]) & (data["m"]==mode[1])]
            label = f'r={s["radius"]:g}'
            for series, style in [(total,"-"), (main,"--")]:
                for ax, field in [(axes[0,0],"JQ"),(axes[0,1],"JE"),(axes[1,0],"deltaQ_node")]:
                    ax.plot(series["t"]/meta["period"],series[field],style,
                            label=label+(" total" if style=="-" else f" ({mode[0]},{mode[1]})"))
            b = balance[balance["r_node"]==s["radius"]]
            axes[1,1].plot(b["t"]/meta["period"], b["continuity"]/initial["Q0"],label=label)
        for ax,title in zip(axes.flat,["Outward number flux","Outward energy flux",
                             "Net outward number at node","Interior balance / initial number"]):
            ax.set(xlabel="Time / original period",ylabel=title)
            ax.legend(fontsize=7)
        r = meta["epsilon"]+np.arange(meta["N"])*meta["dr"]
        from matplotlib.colors import SymLogNorm
        density_scale = max(np.max(np.abs(profiles)),np.finfo(float).tiny)
        colors = SymLogNorm(linthresh=density_scale*1e-9, vmin=-density_scale, vmax=density_scale)
        mesh = axes[2,0].pcolormesh(r,np.array(times)/meta["period"],profiles,
                                    shading="auto",cmap="RdBu_r",norm=colors)
        axes[2,0].set(xlabel="Radius",ylabel="Time / original period")
        fig.colorbar(mesh,ax=axes[2,0],label="Change in radial number density")
        axes[2,1].plot(np.array(times)/meta["period"],departures)
        axes[2,1].set(xlabel="Time / original period",ylabel="Relative departure from a single CN phase")
        filename = "summary.png" if tuple(mode)==(1,1) else f"summary_{mode[0]}_{mode[1]}.png"
        fig.savefig(case/filename, dpi=150)
        plt.close(fig)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("case")
    p.add_argument("directory", type=Path)
    p.add_argument("--corrected", type=Path)
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--mode",type=int,nargs=2,metavar=("L","M"),default=(1,1))
    p = sub.add_parser("binary")
    p.add_argument("B0", type=Path)
    p.add_argument("B1", type=Path)
    p.add_argument("output_csv", type=Path)
    view = sub.add_parser("plot")
    view.add_argument("directory",type=Path)
    view.add_argument("--mode",type=int,nargs=2,metavar=("L","M"),default=(1,1))
    p = sub.add_parser("controls")
    p.add_argument("first", type=Path)
    p.add_argument("second", type=Path)
    p.add_argument("output_json", type=Path)
    p = sub.add_parser("archive")
    p.add_argument("D0", type=Path)
    p.add_argument("B0", type=Path)
    p.add_argument("output_csv", type=Path)
    args = parser.parse_args()
    if args.command == "case":
        analyze(args.directory, args.corrected, not args.no_plots,args.mode)
    elif args.command == "binary":
        compare_binary(args.B0, args.B1, args.output_csv)
    elif args.command == "archive":
        print("Archive initial-state relative difference:",compare_archive(args.D0,args.B0,args.output_csv))
    elif args.command == "plot":
        plot_flux(args.directory,args.mode)
    else:
        compare_controls(args.first,args.second,args.output_json)
