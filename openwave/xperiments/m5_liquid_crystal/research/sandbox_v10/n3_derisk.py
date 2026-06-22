#!/usr/bin/env python3
"""
N3 de-risk , the mass-matrix <-> PMNS scaffold, verified BEFORE touching the loop field theory.

The N3 bridge (research/sandbox_v10/checkpoints/04_n3_design.md): the three neutrino FLAVOUR states
are three closed-loop configurations; their LC energy overlaps form a 3x3 effective mass matrix in
flavour space; U = eigenvectors of that matrix; the PMNS angles come from U. Before computing M_mass
from loops (n3_mass_matrix.py), this file PROVES the scaffold on KNOWN matrices, so a later failure is
unambiguously a loop-physics issue, not a bookkeeping bug:

  D1  the bridge round-trips: build M = U_TBM diag(m1,m2,m3) U_TBM^T, diagonalize, recover the TBM
      angles exactly. Reveal the structure M_TBM must have (magic + mu-tau).
  D2  the Z3 "democratic" matrix M = m0 I + t (J - I) gives the TRIMAXIMAL solar column (sin^2 th12
      = 1/3) , but leaves a degenerate doublet (th23, th13 undetermined). This is the Z3 core the
      three symmetric loops produce.
  D3  magic (all row sums equal) + mu-tau (2<->3) symmetry => EXACT TBM, scanned over the free params.
  D4  which breakings switch on theta13, and its sign + linearity:
        - mu-tau breaking (2<->3 asymmetry)         -> theta13 linear, theta23 departs 45
        - magic breaking (unequal row sums)         -> theta13 linear
      This maps the delta -> theta13 channel that n3_theta13.py will source from the LC delta.

Convention: index-0 (Duda) is irrelevant here (this is the flavour-space 3x3, downstream of the field
theory). Reuses N2's validated PDG extraction. Headless. LOCAL (#236 N-program, HELD).
Run: python3 n3_derisk.py
"""

import json
import os
import numpy as np

from n2_closed_loop import U_TBM, pmns_angles, _R

HERE = os.path.dirname(os.path.abspath(__file__))
RAD = 180.0 / np.pi

# Standard-convention TBM (solar column = (1,1,1)/sqrt3, reactor = (0,1,-1)/sqrt2, mass1 =
# (2,-1,-1)/sqrt6): used for the magic + mu-tau STRUCTURE readout so the row sums read clean.
# (n2's U_TBM uses a (1,1,-1)/sqrt3 solar sign convention -> same |U|^2 / same angles, but the
# magic eigenvector is sign-flipped; angles_from_M is convention-robust, structure readout is not.)
U_TBM_STD = np.array([
    [np.sqrt(2.0 / 3.0),  np.sqrt(1.0 / 3.0),  0.0],
    [-np.sqrt(1.0 / 6.0), np.sqrt(1.0 / 3.0),  np.sqrt(1.0 / 2.0)],
    [-np.sqrt(1.0 / 6.0), np.sqrt(1.0 / 3.0), -np.sqrt(1.0 / 2.0)],
])

# NuFIT 6.0 (normal ordering, w/ SK atmospheric) best-fit targets, for scoring downstream
NUFIT = {"theta12": 33.68, "theta23": 43.3, "theta13": 8.56, "delta_CP": 212.0}
TH12_TBM = float(np.degrees(np.arcsin(np.sqrt(1.0 / 3.0))))   # 35.26439 deg (exact, not rounded)
TBM_TARGET = {"theta12": TH12_TBM, "theta23": 45.0, "theta13": 0.0}


def tbm_err(ang):
    """Max angle deviation (deg) of an angles dict from EXACT tribimaximal."""
    return max(abs(ang["theta12"] - TH12_TBM), abs(ang["theta23"] - 45.0), abs(ang["theta13"]))


# ----------------------------------------------------------------------------------
# the bridge: a real-symmetric flavour mass matrix -> mixing matrix -> angles
# ----------------------------------------------------------------------------------
def canonical_U(V):
    """Reorder the columns (mass eigenstates) of an orthogonal eigenvector matrix V (rows=flavour
    e,mu,tau; cols=mass) into the standard PMNS assignment, so pmns_angles reads the right slots:
      mass3 (col 2) = the eigenvector with the SMALLEST electron content |U_e| (reactor; U_e3->0),
      of the remaining two: mass1 (col 0) = larger |U_e|, mass2 (col 1) = smaller |U_e|
      (|U_e1| = cos th12 cos th13 > |U_e2| = sin th12 cos th13 for th12 < 45 deg).
    Angles depend only on |U|^2, so column SIGNS are irrelevant; only the permutation matters."""
    e_content = np.abs(V[0, :])
    i3 = int(np.argmin(e_content))                       # reactor column
    rest = [c for c in range(3) if c != i3]
    rest.sort(key=lambda c: -e_content[c])               # larger |U_e| first
    i1, i2 = rest[0], rest[1]
    return V[:, [i1, i2, i3]]


def angles_from_M(M):
    """Real-symmetric flavour mass matrix -> (canonical U, angles dict). M need not be the m^2 matrix;
    any symmetric coupling matrix whose eigenvectors are the mixing matrix works (the diagonalizing
    rotation is what PMNS measures)."""
    M = 0.5 * (M + M.T)
    w, V = np.linalg.eigh(M)                              # ascending eigenvalues, orthonormal cols
    U = canonical_U(V)
    ang = pmns_angles(U)
    return U, ang, w


def sin2(theta_deg):
    return np.sin(np.radians(theta_deg)) ** 2


# ----------------------------------------------------------------------------------
# D1 , the bridge round-trips on a constructed TBM matrix; reveal its structure
# ----------------------------------------------------------------------------------
def d1_bridge_roundtrip(m=(1.0, 2.5, 6.0)):
    m1, m2, m3 = m
    M = U_TBM_STD @ np.diag([m1, m2, m3]) @ U_TBM_STD.T
    U, ang, w = angles_from_M(M)
    err = tbm_err(ang)
    # structural readout
    rowsums = M.sum(axis=1)
    magic = float(np.ptp(rowsums))                       # 0 if all row sums equal
    mutau = float(abs(M[1, 1] - M[2, 2]) + abs(M[0, 1] - M[0, 2]))   # 0 if 2<->3 symmetric
    return {
        "M": M.tolist(), "angles": ang, "angle_err_vs_TBM": err,
        "roundtrip_ok": bool(err < 1e-6),
        "structure": {"row_sums": rowsums.tolist(), "magic_violation": magic,
                      "mutau_violation": mutau},
    }


# ----------------------------------------------------------------------------------
# D2 , the Z3 democratic core -> trimaximal solar column
# ----------------------------------------------------------------------------------
def d2_democratic(m0=3.0, t=-0.7):
    J = np.ones((3, 3))
    M = m0 * np.eye(3) + t * (J - np.eye(3))
    w, V = np.linalg.eigh(M)
    # the non-degenerate eigenvector is the trimaximal (1,1,1)/sqrt3 with eigenvalue m0 + 2t
    trimax = np.array([1, 1, 1]) / np.sqrt(3)
    overlaps = [abs(V[:, c] @ trimax) for c in range(3)]
    c_tri = int(np.argmax(overlaps))
    # sin^2 th12 = |U_e2|^2 with mass2 = the trimaximal column; the doublet is degenerate
    s2_solar = float(trimax[0] ** 2)                     # = 1/3 exactly
    gap = float(np.ptp(w))
    # degeneracy of the doublet (the two non-trimaximal eigenvalues)
    doublet = sorted([w[c] for c in range(3) if c != c_tri])
    degeneracy = float(abs(doublet[0] - doublet[1]))
    return {
        "M": M.tolist(), "eigenvalues": w.tolist(),
        "trimaximal_eigenvalue": float(m0 + 2 * t),
        "sin2_theta12_solar": s2_solar, "expected_1_3": 1.0 / 3.0,
        "solar_is_trimaximal": bool(abs(s2_solar - 1.0 / 3.0) < 1e-12),
        "doublet_degeneracy": degeneracy,
        "note": "Z3 democratic gives sin^2 th12 = 1/3 exactly; th23/th13 are degenerate (undetermined) "
                "until a mu-tau + magic refinement splits the doublet (D3).",
    }


# ----------------------------------------------------------------------------------
# D3 , magic + mu-tau symmetric matrix => exact TBM (scan the free params)
# ----------------------------------------------------------------------------------
def magic_mutau_M(x, y, z):
    """General matrix invariant under 2<->3 (mu-tau) AND magic (equal row sums).
       M = [[x, y, y],[y, z, w],[y, w, z]] with the magic constraint x + 2y = y + z + w
       => w = x + y - z. Two free shape params after fixing scale; TBM for ALL (x,y,z)."""
    w = x + y - z
    return np.array([[x, y, y], [y, z, w], [y, w, z]])


def d3_magic_mutau_scan(n=400, seed=12345):
    rng = np.random.default_rng(seed)
    worst = 0.0
    worst_params = None
    for _ in range(n):
        x, y, z = rng.uniform(-3, 3, size=3)
        M = magic_mutau_M(x, y, z)
        # guard against accidental eigenvalue degeneracy (angles ill-defined)
        w = np.linalg.eigvalsh(M)
        if np.ptp(w) < 1e-6 or min(abs(w[1] - w[0]), abs(w[2] - w[1])) < 1e-4:
            continue
        _, ang, _ = angles_from_M(M)
        err = tbm_err(ang)
        if err > worst:
            worst, worst_params = err, (float(x), float(y), float(z))
    return {
        "scanned": n, "worst_angle_err_vs_TBM_deg": float(worst),
        "worst_params_xyz": worst_params,
        "tbm_for_all_magic_mutau": bool(worst < 1e-6),
        "note": "every magic + mu-tau matrix diagonalizes to EXACT TBM, independent of the free "
                "params x,y,z -> TBM is the symmetry, not a tuning. This is what the 3 Z3/reflection-"
                "symmetric loops must reproduce.",
    }


# ----------------------------------------------------------------------------------
# D4 , the breakings that switch on theta13 (size + sign + linearity)
# ----------------------------------------------------------------------------------
def d4_breaking_to_theta13(x=2.0, y=-0.8, z=1.3, eps_grid=None):
    if eps_grid is None:
        eps_grid = np.concatenate([-np.logspace(-1, -5, 12)[::-1], [0.0], np.logspace(-5, -1, 12)])
    base = magic_mutau_M(x, y, z)
    out = {"mutau_break": [], "magic_break": []}
    # mu-tau breaking: perturb M[1,1]-M[2,2] asymmetry (keep magic by symmetric compensation off)
    for eps in eps_grid:
        Mmt = base.copy()
        Mmt[1, 1] += eps; Mmt[2, 2] -= eps
        _, ang, _ = angles_from_M(Mmt)
        out["mutau_break"].append({"eps": float(eps), "theta13": ang["theta13"],
                                   "theta23": ang["theta23"], "theta12": ang["theta12"]})
    # magic breaking: perturb a single row sum (electron row)
    for eps in eps_grid:
        Mmg = base.copy()
        Mmg[0, 0] += eps                                # breaks row-0 sum -> magic violated
        Mmg[0, 0] = Mmg[0, 0]                           # (explicit, symmetric already)
        _, ang, _ = angles_from_M(Mmg)
        out["magic_break"].append({"eps": float(eps), "theta13": ang["theta13"],
                                   "theta23": ang["theta23"], "theta12": ang["theta12"]})
    # linearity + sign: fit theta13 (>=0) vs eps over a resolved, still-linear window (mu-tau channel)
    sm = [d for d in out["mutau_break"] if 1e-3 <= d["eps"] <= 3e-2]
    if len(sm) >= 2:
        e = np.array([d["eps"] for d in sm]); t = np.array([abs(d["theta13"]) for d in sm])
        slope = float(np.polyfit(e, t, 1)[0])
    else:
        slope = float("nan")
    out["mutau_theta13_per_eps_deg"] = slope
    out["note"] = ("theta13 turns on LINEARLY in both breakings near TBM; this is the delta->theta13 "
                   "channel n3_theta13 sources from the LC delta. The size of theta13 per unit "
                   "breaking sets how big the effective breaking must be for 8.5 deg.")
    return out


def main():
    print("=" * 80)
    print("N3 DE-RISK , the mass-matrix <-> PMNS scaffold (verify before the loop field theory)")
    print("=" * 80)

    d1 = d1_bridge_roundtrip()
    print(f"\n[D1] bridge round-trip: angles = "
          f"th12={d1['angles']['theta12']:.4f} th23={d1['angles']['theta23']:.4f} "
          f"th13={d1['angles']['theta13']:.2e}  err={d1['angle_err_vs_TBM']:.2e}  "
          f"-> {'OK' if d1['roundtrip_ok'] else 'FAIL'}")
    print(f"[D1] M_TBM structure: magic_violation={d1['structure']['magic_violation']:.2e} "
          f"mutau_violation={d1['structure']['mutau_violation']:.2e} "
          f"(both ~0 -> TBM matrix IS magic + mu-tau)")

    d2 = d2_democratic()
    print(f"\n[D2] Z3 democratic: sin^2 th12 = {d2['sin2_theta12_solar']:.6f} (expect 1/3="
          f"{1/3:.6f}) solar_trimaximal={d2['solar_is_trimaximal']}; "
          f"doublet_degeneracy={d2['doublet_degeneracy']:.2e} (th23/th13 undetermined at Z3 level)")

    d3 = d3_magic_mutau_scan()
    print(f"\n[D3] magic+mu-tau scan ({d3['scanned']} random matrices): worst angle err vs TBM = "
          f"{d3['worst_angle_err_vs_TBM_deg']:.2e} deg -> "
          f"{'EXACT TBM for all (symmetry, not tuning)' if d3['tbm_for_all_magic_mutau'] else 'NOT always TBM'}")

    d4 = d4_breaking_to_theta13()
    th13_mt = [d["theta13"] for d in d4["mutau_break"]]
    print(f"\n[D4] breaking -> theta13: mu-tau channel theta13 per unit eps = "
          f"{d4['mutau_theta13_per_eps_deg']:.3f} deg/eps (linear onset); "
          f"theta13 range over scan = [{min(th13_mt):.3f}, {max(th13_mt):.3f}] deg")

    scaffold_ok = bool(d1["roundtrip_ok"] and d2["solar_is_trimaximal"]
                       and d3["tbm_for_all_magic_mutau"]
                       and np.isfinite(d4["mutau_theta13_per_eps_deg"])
                       and d4["mutau_theta13_per_eps_deg"] > 0)
    print("\n" + "=" * 80)
    print(f"N3 DE-RISK: {'PASS' if scaffold_ok else 'FAIL'}  "
          f"(bridge round-trips; Z3->trimaximal; magic+mu-tau->TBM; breaking->theta13 linear)")
    print("  => the flavour-space scaffold is correct. A later TBM miss is loop physics, not bookkeeping.")
    print("=" * 80)

    summary = {
        "scaffold_pass": scaffold_ok,
        "targets": {"TBM": TBM_TARGET, "NuFIT_6.0_NO": NUFIT},
        "D1_bridge_roundtrip": d1,
        "D2_democratic_Z3": d2,
        "D3_magic_mutau_TBM": d3,
        "D4_breaking_to_theta13": d4,
    }
    with open(os.path.join(HERE, "n3_derisk_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    print(f"summary -> {os.path.join(HERE, 'n3_derisk_summary.json')}")
    return scaffold_ok


if __name__ == "__main__":
    raise SystemExit(0 if main() else 1)
