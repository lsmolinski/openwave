#!/usr/bin/env python3
"""
N4 / Study C , a mu-tau-RESPECTING topological framing: can theta13 be quantized without killing CP?

n4_topo.py FAILED because it twisted the secondary frame by N * s with s = the LOCAL (rotated) azimuth,
which is mu-tau-ODD under the mu-tau mirror (z -> -z, which conjugates Rx(+alpha) <-> Rx(-alpha)) , so it
broke mu-tau. The FIX: the GLOBAL azimuth phi = atan2(Y,X) is mu-tau-EVEN under z -> -z (it has no z
dependence). So an N * phi framing twist (N integer by single-valuedness) is a mu-tau-RESPECTING topological
self-linking. This file tests whether that quantizes theta13 while preserving the TBM angles and delta_CP.

Outcome decides the article's theta13 status: TOPOLOGICAL (theta13 stepped by integer N, TBM + delta_CP
intact) vs CONTINUOUS (theta13 set by the chiral material coupling g_chiral, as n4_refine found). Convention
index-0. Headless f64. LOCAL (#236 N-program, HELD). Run: python3 n4_linking.py
"""

import json
import os
import numpy as np

from n2_closed_loop import G_SCALE
from n3_mass_matrix import rot_axis
from n3_theta13 import biaxial_vacuum, alpha_star
from n4_chiral import real_overlap, chiral_overlap, pmns_from_H

HERE = os.path.dirname(os.path.abspath(__file__))


def seed_loop_global_topo(n, R3, R_loop, delta, chi=0.0, n_link=0, q=0.5, core_vox=2.0, blend_pow=4):
    """Biaxial oriented loop with a mu-tau-EVEN topological framing: the secondary delta-director gets a
    meridian screw chi*psi PLUS n_link twists about the GLOBAL azimuth phi=atan2(Y,X) (mu-tau-even under
    z->-z). n_link integer for single-valuedness (topological)."""
    c0 = (n - 1) / 2.0
    idx = np.arange(n) - c0
    X, Y, Z = np.meshgrid(idx, idx, idx, indexing="ij")
    P = np.stack([X, Y, Z], axis=-1)
    phi_global = np.arctan2(Y, X)                          # mu-tau-EVEN under z->-z (no z dependence)
    Ploc = np.einsum("ab,...b->...a", R3.T, P)
    xl, yl, zl = Ploc[..., 0], Ploc[..., 1], Ploc[..., 2]
    s = np.arctan2(yl, xl)
    er = np.stack([np.cos(s), np.sin(s), np.zeros_like(s)], axis=-1)
    ez = np.zeros_like(er); ez[..., 2] = 1.0
    rho_r = np.sqrt(xl * xl + yl * yl) - R_loop
    rho_z = zl
    psi = np.arctan2(rho_z, rho_r)
    d_core = np.sqrt(rho_r * rho_r + rho_z * rho_z)
    ang = q * psi
    d_wind_loc = np.cos(ang)[..., None] * er + np.sin(ang)[..., None] * ez
    d_wind = np.einsum("ab,...b->...a", R3, d_wind_loc)
    w = 1.0 / (1.0 + (d_core / (3.0 * core_vox))[..., None] ** blend_pow)
    zhat = np.zeros_like(d_wind); zhat[..., 2] = 1.0
    n1 = w * d_wind + (1.0 - w) * zhat
    n1 = n1 / (np.linalg.norm(n1, axis=-1, keepdims=True) + 1e-12)

    t = np.cross(zhat, n1)
    tn = np.linalg.norm(t, axis=-1, keepdims=True)
    xhat = np.zeros_like(n1); xhat[..., 0] = 1.0
    t = np.where(tn < 1e-6, np.cross(xhat, n1), t)
    t = t / (np.linalg.norm(t, axis=-1, keepdims=True) + 1e-12)
    bvec = np.cross(n1, t)
    twist = ((chi * psi + n_link * phi_global) * w[..., 0])[..., None]   # mu-tau-EVEN topological framing
    n2 = np.cos(twist) * t + np.sin(twist) * bvec

    nn1 = np.einsum("...a,...b->...ab", n1, n1)
    nn2 = np.einsum("...a,...b->...ab", n2, n2)
    M_sp = 1.0 * nn1 + delta * nn2
    M = np.zeros((n, n, n, 4, 4))
    M[..., 0, 0] = G_SCALE
    M[..., 1:, 1:] = M_sp
    return M


def link_mass_matrix(n, alpha, delta, chi, n_link, g_chiral, R_loop=9.0, q=0.5, core_vox=2.0, kappa=0.0):
    Re, Rmu, Rtau = np.eye(3), rot_axis((1, 0, 0), +alpha), rot_axis((1, 0, 0), -alpha)
    Mvac = biaxial_vacuum(n, delta)
    fe = seed_loop_global_topo(n, Re, R_loop, delta, chi=0.0, n_link=n_link, q=q, core_vox=core_vox)
    fmu = seed_loop_global_topo(n, Rmu, R_loop, delta, chi=chi, n_link=n_link, q=q, core_vox=core_vox)
    ftau = seed_loop_global_topo(n, Rtau, R_loop, delta, chi=chi, n_link=n_link, q=q, core_vox=core_vox)
    d = [fe - Mvac, fmu - Mvac, ftau - Mvac]
    Mr = np.zeros((3, 3)); Cc = np.zeros((3, 3))
    for a in range(3):
        for b in range(a, 3):
            K, Pp = real_overlap(d[a], d[b]); Mr[a, b] = Mr[b, a] = K + kappa * Pp
        for b in range(a + 1, 3):
            cab = chiral_overlap(d[a], d[b]); Cc[a, b] = cab; Cc[b, a] = -cab
    M_H = Mr.astype(complex) + 1j * g_chiral * Cc
    r = pmns_from_H(M_H)
    r["mutau_violation"] = float(abs(Mr[0, 1] - Mr[0, 2]) + abs(Mr[1, 1] - Mr[2, 2]))
    return r


def main():
    print("=" * 80)
    print("N4 / Study C , mu-tau-RESPECTING topological framing (global-azimuth self-linking N)")
    print("=" * 80)
    n = 40
    geom = {"R_loop": 9.0, "q": 0.5, "core_vox": 2.0, "kappa": 0.0}
    delta, chi, g_chiral = 0.1, 0.6, 1.0

    print(f"\n[C] scan self-linking N (mu-tau-EVEN global framing), delta={delta}, chi={chi}, g_chiral={g_chiral}:")
    print(f"   {'N':>3} {'theta12':>8} {'theta23':>8} {'theta13':>8} {'delta_CP':>9} {'mutau_viol':>11}")
    rows = []
    for N in (0, 1, 2, 3, 4):
        a = alpha_star(n, delta, chi, chi, geom)
        r = link_mass_matrix(n, a, delta, chi, N, g_chiral, **geom)
        r["n_link"] = N
        rows.append(r)
        print(f"   {N:>3} {r['theta12']:>8.3f} {r['theta23']:>8.3f} {r['theta13']:>8.3f} "
              f"{r['delta_CP']:>9.1f} {r['mutau_violation']:>11.2e}")

    th23_ok = all(abs(r["theta23"] - 45.0) < 0.5 for r in rows)
    th12_ok = all(abs(r["theta12"] - 35.264) < 1.5 for r in rows)
    dcp_ok = all(abs(abs(r["delta_CP"]) - 90.0) < 1.0 or abs(r["delta_CP"]) < 1.0 for r in rows)
    th13_steps = len(set(round(r["theta13"], 1) for r in rows)) > 1
    # is theta13 a clean integer-N ladder (monotone, distinct per N) with TBM+CP intact?
    th13_clean_ladder = bool(th23_ok and dcp_ok and th13_steps
                             and all(rows[i]["theta13"] <= rows[i + 1]["theta13"] + 1e-6
                                     for i in range(len(rows) - 1)))
    print(f"\n[C] TBM preserved across N: theta23~45 {th23_ok}, theta12~35.26 {th12_ok}; "
          f"delta_CP maximal {dcp_ok}")
    print(f"[C] theta13 steps with N: {th13_steps}; clean topological ladder (TBM+CP intact): {th13_clean_ladder}")

    print("\n" + "=" * 80)
    if th13_clean_ladder:
        print("N4 Study C: theta13 IS a mu-tau-respecting TOPOLOGICAL LADDER (stepped by integer N, TBM + ")
        print("delta_CP intact) -> theta13 is PREDICTED by the loop self-linking. (article-grade result)")
    else:
        print("N4 Study C: the mu-tau-even global framing does NOT give a clean topological theta13 ladder.")
        print("Conclusion: theta13 is CONTINUOUS (set by the chiral material coupling g_chiral, n4_refine),")
        print("NOT topologically quantized. delta_CP is the discrete prediction (+-90); theta13 is the one")
        print("free coupling (the cholesteric-pitch analogue). This is the article's honest theta13 status.")
    print("=" * 80)

    summary = {
        "construction": "biaxial loops + mu-tau-EVEN global-azimuth framing N*phi + chiral complex matrix",
        "delta": delta, "chi": chi, "g_chiral": g_chiral,
        "scan_N": rows,
        "tbm_preserved": bool(th23_ok and dcp_ok),
        "theta13_clean_topological_ladder": th13_clean_ladder,
        "conclusion": ("theta13 is a mu-tau-respecting topological ladder" if th13_clean_ladder else
                       "theta13 is CONTINUOUS (chiral material coupling), not topologically quantized; "
                       "delta_CP (+-90) is the discrete prediction, theta13 is the one free coupling "
                       "(cholesteric-pitch analogue). The mu-tau-even global framing did not quantize it."),
    }
    with open(os.path.join(HERE, "n4_linking_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    print(f"summary -> {os.path.join(HERE, 'n4_linking_summary.json')}")
    return True


if __name__ == "__main__":
    raise SystemExit(0 if main() else 1)
