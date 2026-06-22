#!/usr/bin/env python3
"""
N4 , theta13 + delta_CP from a CHIRAL loop coupling (turn the theta13 FIT into a prediction + add CP).

N3 reached the TBM gate with a REAL symmetric mass matrix -> U real -> delta_CP locked to {0,180}, and
theta13 = a free O(0.1) mu-tau asymmetry (chiral in origin). N4 adds the CP-odd (Lifshitz / cholesteric)
CHIRAL coupling the real matrix was missing, making the flavour mass matrix COMPLEX HERMITIAN:

    M_H = M_real_sym  +  i * g_chiral * C ,
    C_ab = INT [ <dx dM_a, dy dM_b> - <dy dM_a, dx dM_b> + (y->z) + (z->x) ]_s     (real, antisymmetric)

C is the curl-like chiral overlap (reflection-ODD = handedness); i*C is Hermitian. U = complex eigenvectors
-> theta13 = |U_e3|^2 AND a genuine delta_CP via the (rephasing-invariant) Jarlskog invariant
J = Im(U_e1 U_mu2 U*_e2 U*_mu1). The chirality is sourced by the secondary delta-twist SCREW with a GLOBAL
handedness (same sign on the loops). g_chiral -> 0 recovers N3 (real, CP-conserving).

Tests: (1) does the chiral term give delta_CP != {0,180}? (2) does it also source theta13? (3) can ONE
handedness give theta13 ~ 8.5 AND delta_CP ~ 212 (NuFIT 6.0 NO) together? Convention index-0. Headless f64,
16-core. LOCAL (#236 N-program, HELD). Run: python3 n4_chiral.py
"""

import json
import os
import numpy as np
from concurrent.futures import ProcessPoolExecutor

from n1_precision_method import _grads, _sdot
from n3_mass_matrix import rot_axis
from n3_theta13 import seed_loop_biaxial, biaxial_vacuum, alpha_star
from n3_derisk import TH12_TBM, NUFIT

HERE = os.path.dirname(os.path.abspath(__file__))


# ----------------------------------------------------------------------------------
# complex-U PMNS extraction (theta's + delta_CP via Jarlskog)
# ----------------------------------------------------------------------------------
def canonical_U_complex(V):
    """Column (mass) assignment by electron content |U_e| (mass3 = smallest), as in n3_derisk."""
    e = np.abs(V[0, :])
    i3 = int(np.argmin(e))
    rest = sorted([c for c in range(3) if c != i3], key=lambda c: -e[c])
    return V[:, [rest[0], rest[1], i3]]


def pmns_from_H(M_H):
    """Complex Hermitian flavour matrix -> (theta12, theta23, theta13, delta_CP) via standard PDG + Jarlskog."""
    w, V = np.linalg.eigh(M_H)                     # Hermitian: real eigenvalues, complex orthonormal cols
    U = canonical_U_complex(V)
    s13_sq = float(abs(U[0, 2]) ** 2)
    th13 = np.degrees(np.arcsin(np.sqrt(min(s13_sq, 1.0))))
    denom = max(1.0 - s13_sq, 1e-15)
    th12 = np.degrees(np.arcsin(np.sqrt(min(abs(U[0, 1]) ** 2 / denom, 1.0))))
    th23 = np.degrees(np.arcsin(np.sqrt(min(abs(U[1, 2]) ** 2 / denom, 1.0))))
    # Jarlskog (rephasing-invariant): J = Im(U_e1 U_mu2 U*_e2 U*_mu1)
    J = float(np.imag(U[0, 0] * U[1, 1] * np.conj(U[0, 1]) * np.conj(U[1, 0])))
    s12, c12 = np.sin(np.radians(th12)), np.cos(np.radians(th12))
    s23, c23 = np.sin(np.radians(th23)), np.cos(np.radians(th23))
    s13, c13 = np.sin(np.radians(th13)), np.cos(np.radians(th13))
    Jmax = s12 * c12 * s23 * c23 * s13 * c13 * c13
    sin_dcp = float(np.clip(J / Jmax, -1.0, 1.0)) if abs(Jmax) > 1e-18 else 0.0
    delta_cp = float(np.degrees(np.arcsin(sin_dcp)))    # principal value; sign from J
    return {"theta12": th12, "theta23": th23, "theta13": th13,
            "delta_CP": delta_cp, "sin_delta_CP": sin_dcp, "J": J,
            "eigenvalues": w.tolist()}


# ----------------------------------------------------------------------------------
# chiral overlap (the antisymmetric, reflection-odd curl-like coupling)
# ----------------------------------------------------------------------------------
def chiral_overlap(dA, dB):
    Ax, Ay, Az = _grads(dA)
    Bx, By, Bz = _grads(dB)
    c = (_sdot(Ax, By) - _sdot(Ay, Bx)
         + _sdot(Ay, Bz) - _sdot(Az, By)
         + _sdot(Az, Bx) - _sdot(Ax, Bz))
    return float(np.sum(c))


def real_overlap(dA, dB):
    Ax, Ay, Az = _grads(dA)
    Bx, By, Bz = _grads(dB)
    K = float(np.sum(_sdot(Ax, Bx) + _sdot(Ay, By) + _sdot(Az, Bz)))
    P = float(np.sum(_sdot(dA[1:-1, 1:-1, 1:-1], dB[1:-1, 1:-1, 1:-1])))
    return K, P


def chiral_mass_matrix(n, alpha, delta, chi, g_chiral, R_loop=9.0, q=0.5, core_vox=2.0, kappa=0.0,
                       tilt_axis=(1.0, 0.0, 0.0)):
    """3 biaxial loops with a GLOBAL secondary-screw handedness chi (same sign on mu & tau; e = ref).
    Real part = grad+field overlap (the N3 matrix); chiral part = i*g_chiral*C (antisymmetric)."""
    Re, Rmu, Rtau = np.eye(3), rot_axis(tilt_axis, +alpha), rot_axis(tilt_axis, -alpha)
    Mvac = biaxial_vacuum(n, delta)
    fe = seed_loop_biaxial(n, Re, R_loop, delta, q=q, core_vox=core_vox, chi=0.0)
    fmu = seed_loop_biaxial(n, Rmu, R_loop, delta, q=q, core_vox=core_vox, chi=+chi)
    ftau = seed_loop_biaxial(n, Rtau, R_loop, delta, q=q, core_vox=core_vox, chi=+chi)
    d = [fe - Mvac, fmu - Mvac, ftau - Mvac]
    Mr = np.zeros((3, 3)); Cc = np.zeros((3, 3))
    for a in range(3):
        for b in range(a, 3):
            K, P = real_overlap(d[a], d[b])
            Mr[a, b] = Mr[b, a] = K + kappa * P
        for b in range(a + 1, 3):
            cab = chiral_overlap(d[a], d[b])
            Cc[a, b] = cab; Cc[b, a] = -cab
    M_H = Mr.astype(complex) + 1j * g_chiral * Cc
    return M_H, Mr, Cc


def _scan_point(args):
    n, delta, chi, g_chiral, geom = args
    try:
        a = alpha_star(n, delta, chi, chi, geom)        # magic crossing (mu,tau same chi -> real mu-tau-sym)
        M_H, Mr, Cc = chiral_mass_matrix(n, a, delta, chi, g_chiral, **geom)
        r = pmns_from_H(M_H)
        r.update({"delta": delta, "chi": chi, "g_chiral": g_chiral, "alpha_star": a,
                  "C_norm": float(np.abs(Cc).max())})
        return r
    except Exception as e:                               # noqa: BLE001
        return {"delta": delta, "chi": chi, "g_chiral": g_chiral, "error": repr(e),
                "theta13": float("nan"), "delta_CP": float("nan")}


def main():
    print("=" * 80)
    print("N4 , theta13 + delta_CP from a CHIRAL loop coupling")
    print("=" * 80)
    n = 40
    geom = {"R_loop": 9.0, "q": 0.5, "core_vox": 2.0, "kappa": 0.0}

    # ---- sanity: g_chiral = 0 recovers the real N3 result (delta_CP in {0,180}, theta13 from sym = 0) ----
    base = _scan_point((n, 0.1, 0.6, 0.0, geom))
    print(f"\n[check] g_chiral=0 (real N3): theta13={base['theta13']:.3e}  delta_CP={base['delta_CP']:.1f}  "
          f"J={base['J']:.2e}  (expect theta13~0, J~0)")

    # ---- does the chiral term turn on theta13 AND delta_CP? scan g_chiral at fixed delta, chi ----
    gs = np.linspace(0.0, 2.0, 11)
    with ProcessPoolExecutor(max_workers=16) as ex:
        scan_g = list(ex.map(_scan_point, [(n, 0.1, 0.6, float(g), geom) for g in gs]))
    print(f"\n[scan g_chiral] delta=0.1, chi=0.6:")
    print(f"   {'g_chiral':>9} {'theta12':>8} {'theta23':>8} {'theta13':>8} {'delta_CP':>9} {'J':>11}")
    for r in scan_g:
        print(f"   {r['g_chiral']:>9.2f} {r['theta12']:>8.3f} {r['theta23']:>8.3f} {r['theta13']:>8.3f} "
              f"{r['delta_CP']:>9.2f} {r['J']:>11.3e}")

    # ---- 2D scan (chi, g_chiral) at delta=0.1: find a simultaneous theta13~8.5 AND delta_CP~212 ----
    chis = np.linspace(0.2, 1.2, 9)
    gs2 = np.linspace(0.2, 2.0, 9)
    grid = [(n, 0.1, float(c), float(g), geom) for c in chis for g in gs2]
    with ProcessPoolExecutor(max_workers=16) as ex:
        scan2 = list(ex.map(_scan_point, grid))
    # score vs NuFIT (theta13, delta_CP)
    def score(r):
        if not np.isfinite(r.get("theta13", np.nan)):
            return 1e9
        return abs(r["theta13"] - NUFIT["theta13"]) + 0.1 * abs(abs(r["delta_CP"]) - (180 - (212 - 180)))
    best = min(scan2, key=score)
    print(f"\n[best simultaneous] chi={best['chi']:.2f} g_chiral={best['g_chiral']:.2f}: "
          f"theta13={best['theta13']:.2f} (tgt 8.56)  delta_CP={best['delta_CP']:.1f} "
          f"(tgt 212 -> |dCP| near 32 from CP-cons or ~180-)  J={best['J']:.3e}")
    th13_range = [min(r['theta13'] for r in scan2 if np.isfinite(r['theta13'])),
                  max(r['theta13'] for r in scan2 if np.isfinite(r['theta13']))]
    dcp_range = [min(r['delta_CP'] for r in scan2 if np.isfinite(r['delta_CP'])),
                 max(r['delta_CP'] for r in scan2 if np.isfinite(r['delta_CP']))]
    print(f"[ranges] theta13 in [{th13_range[0]:.2f},{th13_range[1]:.2f}]  "
          f"delta_CP in [{dcp_range[0]:.1f},{dcp_range[1]:.1f}]")

    chiral_gives_cp = bool(max(abs(r['J']) for r in scan_g if np.isfinite(r['J'])) > 1e-6)
    chiral_gives_th13 = bool(max(r['theta13'] for r in scan_g if np.isfinite(r['theta13'])) > 1.0)
    print("\n" + "=" * 80)
    print(f"N4 (first pass): chiral coupling gives delta_CP != 0 = {chiral_gives_cp}; "
          f"gives theta13 = {chiral_gives_th13}")
    print("=" * 80)

    summary = {
        "mechanism": "complex Hermitian M_H = M_real + i*g_chiral*C (Lifshitz/cholesteric chiral overlap); "
                     "global secondary-screw handedness chi; delta_CP via Jarlskog",
        "check_g0": base, "scan_g_chiral": scan_g,
        "best_simultaneous": best, "theta13_range": th13_range, "delta_CP_range": dcp_range,
        "chiral_gives_cp": chiral_gives_cp, "chiral_gives_theta13": chiral_gives_th13,
        "nufit_target": {"theta13": NUFIT["theta13"], "delta_CP": NUFIT["delta_CP"]},
    }
    with open(os.path.join(HERE, "n4_chiral_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    print(f"summary -> {os.path.join(HERE, 'n4_chiral_summary.json')}")
    return True


if __name__ == "__main__":
    raise SystemExit(0 if main() else 1)
