"""Read-only check of a completed midpoint e=0.4, q=1e-5, p~50 run and its dashboard."""
import argparse
import json
import sys
import threading
import tomllib
from http.server import ThreadingHTTPServer
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "dashboard"))
import live_dashboard as dashboard


def check(folder):
    cfg = tomllib.loads((folder / "config.toml").read_text())
    sim, binary = cfg["simulation"], cfg["binary"]
    e, q, rp = binary["e"], binary["q"], binary["R"]
    assert e == .4 and q == 1e-5 and abs(rp*(1+e)-50) < .01
    assert sim["integrator"] == (folder/"integrator.txt").read_text().strip() == "midpoint_cn"
    assert not sim["evolve_orbit"] and not binary["counter"]

    def column(name):
        return np.loadtxt(folder/(name+".txt"), ndmin=1)

    t = column("t")
    n = np.sqrt((1+q)/(rp/(1-e))**3)
    period = 2*np.pi/n
    ramp = sim.get("ramp_periods", 0.)*period
    phase = n*t
    if ramp > 0:
        x = np.minimum(t/ramp, 1.)
        phase = n*(ramp*x**3*(1-x/2)+np.maximum(0., t-ramp))
    np.testing.assert_allclose(column("phase_orb"), phase, rtol=3e-12, atol=1e-13)
    anomaly = phase.copy()
    for _ in range(10):
        anomaly -= (anomaly-e*np.sin(anomaly)-phase)/(1-e*np.cos(anomaly))
    separation = rp/(1-e)*(1-e*np.cos(anomaly))
    np.testing.assert_allclose(column("R_orb"), separation, rtol=3e-12)
    np.testing.assert_array_equal(column("M_1"), np.ones_like(t))
    np.testing.assert_array_equal(column("q"), np.full_like(t,q))

    dashboard.Handler.root = dashboard.Handler.root_real = str(folder)
    dashboard.Handler.max_points = len(t)
    with ThreadingHTTPServer(("127.0.0.1", 0), dashboard.Handler) as server:
        threading.Thread(target=server.serve_forever, daemon=True).start()
        def get(path):
            with urlopen(f"http://127.0.0.1:{server.server_port}"+path,timeout=60) as response:
                return response.read()
        try:
            assert b"<html" in get("/").lower()
            assert json.loads(get("/api/runs"))[0]["params"]["q"] == q
            query = urlencode({"run": str(folder)})
            data = json.loads(get("/api/data?"+query))
            assert data["ready"] and data["theory"] is None
            assert data["info"]["rows"] == len(t)
            assert data["others"] == list(range(len(data["info"]["modes"])))
            assert [1,1] in data["info"]["modes"]
            assert "midpoint_cn" in json.dumps(data["info"]["setup"])
            step = max(data["info"]["snapshots"])
            assert json.loads(get(f"/api/snapshot?{query}&step={step}&grid=64"))["step"] == step
        finally:
            server.shutdown()
    series = data["series"]
    np.testing.assert_allclose(series["t_over_P"], t/period, rtol=6e-9)
    number = column("norms").sum(axis=1)
    mu = cfg["atom"]["mu"]
    mass = mu*number
    eps = column("energy")/number
    current = column("flux").sum(axis=1)/mu
    power = (column("energy_flux").sum(axis=1)+eps*current)/(q*q*mass)
    rate = (mu+eps)*current/(q*q*mass)
    np.testing.assert_allclose(series["total_energy_flux"],power,rtol=6e-7,atol=1e-280)
    np.testing.assert_allclose(series["total_mass_flux"],rate,rtol=6e-7,atol=1e-280)
    print(f"{folder.name}: {len(t)} steps, t={t[-1]}, P={period:.9f}, "
          f"maximum separation error={max(abs(column('R_orb')-separation)):.3g}; dashboard verified")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folders", nargs="+", type=Path)
    for path in parser.parse_args().folders:
        check(path.resolve())
