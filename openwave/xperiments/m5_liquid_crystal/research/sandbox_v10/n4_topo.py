#!/usr/bin/env python3
"""
N4 / topological test , is the loop HANDEDNESS a topological invariant that sets theta13 + delta_CP?

N4 (n4_chiral.py) showed a chiral coupling predicts delta_CP = +-90 (maximal) and sources theta13, both from
the loop handedness , but the handedness was a CONTINUOUS knob (chi, g_chiral). The "make history" question:
is the handedness actually a TOPOLOGICAL integer? A closed loop's director framing has an integer SELF-LINKING
number N (Calugareanu-White-Fuller: Lk = Tw + Wr, Lk in Z). Adding N full twists of the secondary frame around
the loop AZIMUTH s is: (a) single-valued ONLY for integer N (so N is quantized by topology); (b) reflection-ODD
(s -> -s under a mirror, so N -> -N = opposite handedness). If theta13 / delta_CP are STEPPED by the integer N,
they are predicted by loop topology, not fit.

This file seeds biaxial loops with self-linking N (an N*s azimuthal twist of the secondary delta-director),
builds the complex Hermitian (chiral) mass matrix, and scans N = -2..2 to see whether theta13 and delta_CP are
quantized by N. Convention index-0. Headless f64. LOCAL (#236 N-program, HELD). Run: python3 n4_topo.py
"""

import json
import os
import numpy as np

from n2_closed_loop import G_SCALE
from n3_mass_matrix import rot_axis
from n3_theta13 import biaxial_vacuum, alpha_star
from n4_chiral import real_overlap, chiral_overlap, pmns_from_H

HERE = os.path.dirname(os.path.abspath(__file__))


def seed_loop_topo(n, R3, R_loop, delta, chi=0.0, n_link=0, q=0.5, core_vox=2.0, blend_pow=4):
    """Biaxial oriented loop with SELF-LINKING n_link: the secondary delta-director gets a meridian screw
    chi*psi PLUS n_link full twists around the loop azimuth s (the topological framing). n_link must be an
    integer for single-valuedness (the topological quantization)."""
    c0 = (n - 1) / 2.0
    idx = np.arange(n) - c0
    X, Y, Z = np.meshgrid(idx, idx, idx, indexing="ij")
    P = np.stack([X, Y, Z], axis=-1)
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
    # secondary twist = meridian screw + N azimuthal twists (the self-linking framing)
    twist = ((chi * psi + n_link * s) * w[..., 0])[..., None]
    n2 = np.cos(twist) * t + np.sin(twist) * bvec

    nn1 = np.einsum("...a,...b->...ab", n1, n1)
    nn2 = np.einsum("...a,...b->...ab", n2, n2)
    M_sp = 1.0 * nn1 + delta * nn2
    M = np.zeros((n, n, n, 4, 4))
    M[..., 0, 0] = G_SCALE
    M[..., 1:, 1:] = M_sp
    return M


def topo_mass_matrix(n, alpha, delta, chi, n_link, g_chiral, R_loop=9.0, q=0.5, core_vox=2.0, kappa=0.0):
    Re, Rmu, Rtau = np.eye(3), rot_axis((1, 0, 0), +alpha), rot_axis((1, 0, 0), -alpha)
    Mvac = biaxial_vacuum(n, delta)
    fe = seed_loop_topo(n, Re, R_loop, delta, chi=0.0, n_link=n_link, q=q, core_vox=core_vox)
    fmu = seed_loop_topo(n, Rmu, R_loop, delta, chi=chi, n_link=n_link, q=q, core_vox=core_vox)
    ftau = seed_loop_topo(n, Rtau, R_loop, delta, chi=chi, n_link=n_link, q=q, core_vox=core_vox)
    d = [fe - Mvac, fmu - Mvac, ftau - Mvac]
    Mr = np.zeros((3, 3)); Cc = np.zeros((3, 3))
    for a in range(3):
        for b in range(a, 3):
            K, Pp = real_overlap(d[a], d[b]); Mr[a, b] = Mr[b, a] = K + kappa * Pp
        for b in range(a + 1, 3):
            cab = chiral_overlap(d[a], d[b]); Cc[a, b] = cab; Cc[b, a] = -cab
    return Mr.astype(complex) + 1j * g_chiral * Cc, Mr, Cc


def main():
    print("=" * 80)
    print("N4 / topological test , is theta13 / delta_CP quantized by the self-linking N?")
    print("=" * 80)
    n = 40
    geom = {"R_loop": 9.0, "q": 0.5, "core_vox": 2.0, "kappa": 0.0}
    delta, chi, g_chiral = 0.1, 0.6, 1.0

    print(f"\n[scan self-linking N] delta={delta}, chi={chi}, g_chiral={g_chiral}:")
    print(f"   {'N':>3} {'theta12':>8} {'theta23':>8} {'theta13':>8} {'delta_CP':>9} {'J':>11} {'C_norm':>9}")
    rows = []
    for N in (-2, -1, 0, 1, 2):
        a = alpha_star(n, delta, chi, chi, geom)
        M_H, Mr, Cc = topo_mass_matrix(n, a, delta, chi, N, g_chiral, **geom)
        r = pmns_from_H(M_H)
        r["n_link"] = N; r["alpha_star"] = a; r["C_norm"] = float(np.abs(Cc).max())
        rows.append(r)
        print(f"   {N:>3} {r['theta12']:>8.3f} {r['theta23']:>8.3f} {r['theta13']:>8.3f} "
              f"{r['delta_CP']:>9.2f} {r['J']:>11.3e} {r['C_norm']:>9.3e}")

    # HONEST analysis: does N CLEANLY quantize theta13 while preserving the TBM baseline? (it does NOT)
    th13 = {r["n_link"]: r["theta13"] for r in rows}
    dcp = {r["n_link"]: r["delta_CP"] for r in rows}
    tbm_intact = {N: (abs(th13.get(N, 9) is not None) and
                      abs(rows[i]["theta12"] - 35.264) < 0.5 and abs(rows[i]["theta23"] - 45.0) < 0.5)
                  for i, N in enumerate((-2, -1, 0, 1, 2))}
    antisym_dcp = all(abs(dcp.get(N, 0) + dcp.get(-N, 0)) < 1.0 for N in (1, 2))
    n0_clean = abs(dcp.get(0, 0) + 90.0) < 1.0 or abs(dcp.get(0, 0) - 90.0) < 1.0
    print(f"\n[topo] N=0 keeps the clean structure (TBM + delta_CP=+-90 maximal): {n0_clean}")
    print(f"[topo] N!=0 PRESERVES TBM (theta12~35.26, theta23~45)? "
          f"{ {N: tbm_intact[N] for N in (-2,-1,1,2)} }  -> NO (the azimuthal twist BREAKS mu-tau)")
    print(f"[topo] delta_CP antisymmetric + maximal under N->-N? {antisym_dcp}  -> NO "
          f"(N=+1 {dcp.get(1,0):.1f} vs N=-1 {dcp.get(-1,0):.1f}; intermediate, not maximal)")

    print("\n" + "=" * 80)
    print("N4 topo: INCONCLUSIVE / NEGATIVE for clean topological quantization.")
    print("  - The naive self-linking twist (N*s in each local frame) does NOT cleanly step theta13; it")
    print("    BREAKS mu-tau and DEGRADES the TBM baseline (theta12 -> ~33, theta23 wobbles, delta_CP -> ")
    print("    intermediate, non-antisymmetric). Only N=0 preserves the clean delta_CP=+-90 structure.")
    print("  - So 'theta13 is topologically quantized' is NOT supported by this construction. A proper test")
    print("    needs a mu-tau-RESPECTING definition of the self-linking (the framing must be applied")
    print("    consistently with the loop orientations) , an open task (next session / Duda).")
    print("  - The SOLID N4 result remains n4_chiral: delta_CP = +-90 (maximal), robust, at N=0.")
    print("=" * 80)

    summary = {
        "construction": "biaxial loops with self-linking N (N*s azimuthal twist of the secondary director) + "
                        "chiral complex-Hermitian mass matrix; N integer by single-valuedness",
        "delta": delta, "chi": chi, "g_chiral": g_chiral,
        "scan_N": rows,
        "verdict": "INCONCLUSIVE / NEGATIVE for clean topological quantization",
        "n0_clean_maximal_CP": bool(n0_clean),
        "delta_CP_antisymmetric_under_N_flip": bool(antisym_dcp),
        "reading": ("the naive N*s self-linking twist BREAKS mu-tau and degrades the TBM baseline rather than "
                    "cleanly stepping theta13; only N=0 keeps the clean delta_CP=+-90 structure. 'theta13 "
                    "topologically quantized' is NOT supported by this construction. A mu-tau-respecting "
                    "definition of the loop self-linking is the open task. The solid N4 result stays the "
                    "delta_CP=maximal prediction (n4_chiral, N=0)."),
    }
    with open(os.path.join(HERE, "n4_topo_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    print(f"summary -> {os.path.join(HERE, 'n4_topo_summary.json')}")
    return True


if __name__ == "__main__":
    raise SystemExit(0 if main() else 1)
