"""End-to-end isolated diagnostic regression using independent small dense solves."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile

try:
    import numpy as np
    import scipy
except ImportError:
    sys.exit(77)

from analyze import analyze, digest, read_samples, read_production_fluxes, write_modes
from spectrum import analytical_target, lanczos_weights, load_operator, sign_transform


def test(exe, root):
    cfg = root/"config.toml"
    cfg.write_text("""[grid]
N=101
L=200.0
[binary]
q=0.1
R=16.0
e=0.4
[atom]
n=2
l=1
m=1
mu=0.5
M_c=0.25
[simulation]
integrator="midpoint_cn"
timestep_rule="normalized_geometry"
startup_steps_per_period=60
orbit_steps_per_period=60
num_periods=1
l_max=1
save_every=5
epsilon=0.0
r_ex=30.0
[parallel]
threads=1
[cap]
L_factor=0.8
W=0.1
n=4
""")
    modes = [[0,0], [1,-1], [1,1]]
    settings = dict(N=101, dr=2., epsilon=0., mu=0.5, modes=modes, x=16./11)
    state = analytical_target(settings,32)
    source = root/"source.bin"
    write_modes(source, modes, state)
    hashes = [digest(cfg),digest(source)]
    def run(file, destination, cases, periods=".001"):
        cmd = [str(exe), str(cfg), str(file), "--output", str(destination),
               "--cases",cases,"--periods",periods,"--dt","0.25","--radii","30,60"]
        subprocess.run(cmd,check=True,stdout=subprocess.DEVNULL,cwd=root)
        return cmd
    output = root/"controls"
    cmd = run(source,output,"C0,D0")
    assert subprocess.run(cmd,cwd=root,capture_output=True).returncode != 0
    corrected = root/"corrected.bin"
    for name in ["C0","D0"]:
        case = output/name
        report = analyze(case, corrected if name=="D0" else None, plots=False)
        meta, spec = report["manifest"], report["spectrum"]
        np.testing.assert_allclose(meta["x"], 0. if name=="C0" else settings["x"])
        np.testing.assert_allclose(report["initial"]["Q0"],0.5*meta["input_norm"],rtol=1e-14)
        assert spec["fidelity"] > 1-1e-12 and spec["endpoint_fraction"] < 1e-14
        assert report["max_norm_drift"] < 1e-9*report["initial"]["Q0"]
        assert max(s["max_balance_fraction"] for s in report["surfaces"]) < 1e-9
        flux = read_samples(case/"flux.csv")
        sphere = flux[(flux["step"]>0) & (flux["r_node"]==30) & (flux["mode"]=="mode")]
        t, lm, q, e = read_production_fluxes(next(case.glob("Simulation_*")))
        np.testing.assert_array_equal(lm,modes)
        np.testing.assert_allclose(sphere["JQ"].reshape(q.shape),q,rtol=1e-11,atol=1e-18)
        np.testing.assert_allclose(sphere["JE"].reshape(e.shape),e,rtol=1e-11,atol=1e-18)
        u = np.fromfile(case/"state_0.c128", "<c16").reshape(3,101)
        h,_ = load_operator(case,297)
        eye = np.eye(h.shape[0])
        # Include the exact forcing from the two frozen endpoint samples.
        ref = u[:,1:-1].flatten()
        forcing = np.zeros_like(u[:,1:-1])
        forcing[:,[0,-1]] = -u[:,[0,-1]]/(2*meta["mu"]*meta["dr"]**2)
        for dt in np.diff(np.r_[0.,t]):
            ref = np.linalg.solve(eye+0.5j*dt*h.toarray(),
                                  (eye-0.5j*dt*h.toarray())@ref-1j*dt*forcing.ravel())
        final = np.fromfile(case/f'state_{meta["steps"]}.c128',"<c16").reshape(3,101)
        np.testing.assert_allclose(final[:,1:-1].ravel(),ref,atol=2e-10,rtol=1e-9)
    run(corrected,root/"corrected_run","D1")
    corrected_report = analyze(root/"corrected_run"/"D1", plots=False)
    assert corrected_report["spectrum"]["cluster_weight"] > 1-1e-12
    np.testing.assert_allclose(corrected_report["manifest"]["input_norm"],meta["input_norm"],rtol=1e-13)

    # Known single eigenstate: CN changes only its complex phase.
    vals, vecs = np.linalg.eigh(h.toarray())
    chosen = np.argmin(abs(vals+0.5**3/8))
    stationary = np.zeros_like(u)
    stationary[:,1:-1] = vecs[:,chosen].reshape(3,-1)
    eigenfile = root/"eigen.bin"
    write_modes(eigenfile,modes,stationary)
    run(eigenfile,root/"stationary","D1")
    eigen_report = analyze(root/"stationary"/"D1", plots=False)
    assert max(eigen_report["relative_departure_from_single_CN_phase"]) < 1e-8

    # An outward packet of positive energy gives positive number and energy flux.
    packet = np.zeros_like(u)
    radii = 2*np.arange(101)
    packet[2] = np.exp(-((radii-30)/8)**2 + 0.5j*radii)
    packet[:,[0,-1]]=0
    packetfile=root/"packet.bin"
    write_modes(packetfile,modes,packet)
    run(packetfile,root/"packet","D0")
    packet_flux=read_samples(root/"packet"/"D0"/"flux.csv")
    first=packet_flux[(packet_flux["step"]==0)&(packet_flux["mode"]=="total")][0]
    assert first["JQ"]>0 and first["JE"]>0
    weights=lanczos_weights(np.diag([-2.,-1.,1.,2.]),np.array([1.,2.,3.,4.]),4)
    np.testing.assert_allclose(weights["estimate"],25./30)
    assert weights["lower"] <= 25./30 <= weights["upper"]
    transform,shift,error=sign_transform(h,99,-0.5**3/8)
    random=np.random.default_rng(17)
    probe=random.normal(size=297)+1j*random.normal(size=297)
    as_radial=lambda v:v.reshape(3,99).T.flatten()
    exact=probe-shift*np.linalg.solve(h.toarray()+shift*np.eye(297),probe)
    np.testing.assert_allclose(transform@as_radial(probe),as_radial(exact),atol=1e-11)
    assert error<1e-10
    assert hashes==[digest(cfg),digest(source)]


if __name__=="__main__":
    with tempfile.TemporaryDirectory() as temp:
        test(Path(sys.argv[1]).resolve(),Path(temp))
