#!/usr/bin/env python3
"""
N3 / S2 , the parameter search + the ★ TBM GATE (TBM angles from the closed-loop field dynamics).

S1 (n3_mass_matrix.py) showed the loop overlaps give a mu-tau matrix whose MAGIC residual (x+y)-(z+w)
crosses zero as the mirror tilt alpha varies, and theta12 passes through the trimaximal 35.264 deg right
at the crossing (theta23=45, theta13=0 hold by the mirror symmetry). This file:

  (1) ROOT-FINDS the magic crossing alpha*(geometry) and confirms the 3 TBM angles emerge there from the
      DYNAMICS (the gate, 10a "the 3 angles emerge from the simulation within tolerance");
  (2) MAPS alpha* and the gate quality across the other geometry knobs (R_loop, core_vox, kappa, q) and
      grid resolution n, in parallel over 16 cores, to test whether the TBM gate is ROBUST (a real locus)
      or a fine-tuned coincidence , and whether alpha* is a special angle;
  (3) writes the gate locus + a verdict for the findings.

Honest FAIL clause (10a gate): if no geometry gives the 3 TBM angles within tolerance, that is the
reportable result (closed-loop SO(3) dynamics do not yield TBM). Convention index-0. Headless f64.
LOCAL (#236 N-program, HELD). Run: python3 n3_search.py
"""

import json
import os
import numpy as np
from concurrent.futures import ProcessPoolExecutor
from scipy.optimize import brentq

from n3_mass_matrix import flavour_mass_matrix
from n3_derisk import TH12_TBM

HERE = os.path.dirname(os.path.abspath(__file__))
TBM_TOL = 0.05            # deg; "the 3 angles emerge within tolerance" = the gate pass bar


def magic_residual(alpha, geom):
    r = flavour_mass_matrix(alpha=alpha, **geom)
    return r["magic_residual"]


def find_alpha_star(geom, bracket=(0.45, 1.15)):
    """Root-find the magic crossing alpha* in `bracket`; return (alpha*, full result at alpha*)."""
    lo, hi = bracket
    flo, fhi = magic_residual(lo, geom), magic_residual(hi, geom)
    if flo * fhi > 0:
        # widen once
        for b in (0.3, 1.3, 0.2, 1.4):
            fb = magic_residual(b, geom)
            if flo * fb < 0:
                hi, fhi = b, fb; break
            if fhi * fb < 0:
                lo, flo = b, fb; break
        if flo * fhi > 0:
            return None, None
    astar = brentq(lambda a: magic_residual(a, geom), lo, hi, xtol=1e-6, rtol=1e-10)
    return astar, flavour_mass_matrix(alpha=astar, **geom)


def _worker(geom):
    """Find the gate for one geometry; return a compact record."""
    try:
        astar, r = find_alpha_star(geom)
    except Exception as e:                                   # noqa: BLE001 (record, don't crash the map)
        return {"geom": geom, "error": repr(e)}
    if astar is None:
        return {"geom": geom, "alpha_star": None, "gate_pass": False, "reason": "no magic crossing"}
    a = r["angles"]
    return {
        "geom": geom, "alpha_star": float(astar),
        "angles": a, "tbm_err": r["tbm_err"],
        "eigenvalues": r["eigenvalues"], "magic_residual": r["magic_residual"],
        "mutau_violation": r["mutau_violation"],
        "gate_pass": bool(r["tbm_err"] < TBM_TOL),
    }


def main():
    print("=" * 80)
    print("N3 / S2 , parameter search + the TBM GATE (TBM from closed-loop dynamics)")
    print("=" * 80)

    # ---- (1) the reference gate: alpha* at the default geometry, all 3 angles ----
    geom0 = {"n": 48, "R_loop": 9.0, "delta": 1.0e-10, "q": 0.5, "core_vox": 2.0, "kappa": 0.0}
    astar0, r0 = find_alpha_star(geom0)
    a0 = r0["angles"]
    print(f"\n[GATE] reference geometry {geom0}")
    print(f"[GATE] alpha* = {astar0:.6f} rad = {np.degrees(astar0):.4f} deg   (magic crossing)")
    print(f"[GATE] angles AT alpha*: th12={a0['theta12']:.4f}  th23={a0['theta23']:.4f}  "
          f"th13={a0['theta13']:.2e}")
    print(f"[GATE] vs TBM (35.2644/45/0): max err = {r0['tbm_err']:.4f} deg  -> "
          f"{'PASS' if r0['tbm_err'] < TBM_TOL else 'FAIL'}")
    print(f"[GATE] eigenvalues (mass spectrum) = {[f'{v:.4e}' for v in r0['eigenvalues']]}")

    # ---- (2) robustness map: alpha* + gate quality across geometry knobs (parallel) ----
    grid = []
    for n in (40, 48, 56):
        for R_loop in (7.0, 9.0, 11.0):
            for core_vox in (1.5, 2.0, 2.5):
                for kappa in (0.0, 0.25, 1.0):
                    grid.append({"n": n, "R_loop": R_loop, "delta": 1.0e-10, "q": 0.5,
                                 "core_vox": core_vox, "kappa": kappa})
    print(f"\n[MAP] scanning {len(grid)} geometries over 16 cores for the gate locus ...")
    with ProcessPoolExecutor(max_workers=16) as ex:
        records = list(ex.map(_worker, grid))

    found = [r for r in records if r.get("alpha_star") is not None]
    passes = [r for r in records if r.get("gate_pass")]
    astars = np.array([r["alpha_star"] for r in found])
    tbm_errs = np.array([r["tbm_err"] for r in found])
    print(f"[MAP] magic crossing found in {len(found)}/{len(grid)} geometries; "
          f"gate PASS (<{TBM_TOL} deg) in {len(passes)}/{len(grid)}")
    if len(found):
        print(f"[MAP] alpha* range = [{astars.min():.4f}, {astars.max():.4f}] rad "
              f"= [{np.degrees(astars.min()):.2f}, {np.degrees(astars.max()):.2f}] deg")
        print(f"[MAP] TBM err range = [{tbm_errs.min():.4f}, {tbm_errs.max():.4f}] deg "
              f"(median {np.median(tbm_errs):.4f})")

    # is alpha* a SPECIAL angle? compare to candidate magic angles
    candidates = {"arctan(sqrt2)=54.7deg": np.arctan(np.sqrt(2)),
                  "arccos(1/sqrt3)=54.7deg": np.arccos(1/np.sqrt(3)),
                  "pi/4=45deg": np.pi/4, "arctan(1/sqrt2)=35.26deg": np.arctan(1/np.sqrt(2))}
    amed = float(np.median(astars)) if len(found) else float("nan")
    nearest = min(candidates.items(), key=lambda kv: abs(kv[1] - amed)) if len(found) else (None, None)
    print(f"[MAP] median alpha* = {amed:.4f} rad ({np.degrees(amed):.2f} deg); "
          f"nearest candidate magic angle: {nearest[0]} (delta {abs(nearest[1]-amed):.3f} rad)")

    verdict = ("GATE PASS , the 3 TBM angles emerge from the closed-loop field dynamics at the magic "
               "crossing alpha*, robustly across geometry") if len(passes) >= 0.8 * len(grid) else (
               "GATE PARTIAL , TBM emerges at the magic crossing but quality varies with geometry") \
        if len(passes) else "GATE FAIL , no geometry yields TBM within tolerance (reportable)"
    print("\n" + "=" * 80)
    print(f"N3 / S2 TBM GATE: {verdict}")
    print("=" * 80)

    summary = {
        "gate_tolerance_deg": TBM_TOL,
        "reference": {"geom": geom0, "alpha_star_rad": float(astar0),
                      "alpha_star_deg": float(np.degrees(astar0)), "angles": a0,
                      "tbm_err_deg": r0["tbm_err"], "eigenvalues": r0["eigenvalues"]},
        "robustness_map": {
            "n_geometries": len(grid), "magic_crossing_found": len(found),
            "gate_pass": len(passes),
            "alpha_star_rad_range": [float(astars.min()), float(astars.max())] if len(found) else None,
            "alpha_star_deg_median": np.degrees(amed) if len(found) else None,
            "tbm_err_deg_range": [float(tbm_errs.min()), float(tbm_errs.max())] if len(found) else None,
            "nearest_candidate_magic_angle": {"name": nearest[0],
                                              "delta_rad": float(abs(nearest[1]-amed))} if len(found) else None,
        },
        "records": records,
        "verdict": verdict,
    }
    with open(os.path.join(HERE, "n3_search_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    print(f"summary -> {os.path.join(HERE, 'n3_search_summary.json')}")
    return len(passes) > 0


if __name__ == "__main__":
    raise SystemExit(0 if main() else 1)
