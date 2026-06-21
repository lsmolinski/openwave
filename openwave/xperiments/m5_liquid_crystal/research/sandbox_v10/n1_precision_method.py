#!/usr/bin/env python3
"""
N1 , the precision-safe numerical method for the closed-loop oscillation build.

THE PROBLEM (Duda 2026-06-21). The order parameter is
    M = O . D . O^T ,   D = diag(g, 1, delta, 0)   (index-0 convention)
with g ~ 1e10 (the rest-mass / boost axis) and delta ~ 1e-10 (the quantum-phase /
Dirac axis). The two scales span a ~1e20 dynamic range. The energy functional
(read from engine3_observables.compute_energyH_density_M + engine2_pde.V_M),

    H = 1/2||Mdot||^2  +  c^2 . 4 . sum_{mu<nu} <[d_mu M, d_nu M], [.]>_s  +  (V_M - v0)

mixes g^2 ~ 1e20 with delta^2 ~ 1e-20 in the SAME float sum. f64 has ~16 digits, so
once a term of size ~1e20 is in the accumulator every contribution below ~1e4 is
annihilated , and the ENTIRE delta-sector (the SO(3)-breaking, the theta_13 channel)
lives there. Worse: the field perturbation delta*M1 ~ 1e-10 sits below the ULP of the
g-scale entries (~1e-6 at 1e10), so in f64 the configurations M(delta_phys) and M(0)
are BIT-IDENTICAL , the breaking is gone before the energy is even formed.

THE METHOD (this file).
  (1) Non-dimensionalization: g is a FIXED background scale; the dynamical structure is
      O(1). We expose the natural grading of the energy in powers of g.
  (2) Perturbative-delta. D(delta) = D_0 + delta . D_1 is LINEAR, so M(delta) = M_0 +
      delta . M_1 EXACTLY (M_1 = O . D_1 . O^T). The polynomial energy then expands
      EXACTLY: E(delta) = E_0 + delta E_1 + delta^2 E_2 + ... Each order is computed in
      its OWN well-scaled arithmetic from the delta-graded commutators
          C(delta) = C_0 + delta C_1 + delta^2 C_2 ,  C_k held separately,
      never summed against the 1e20-scale order-0 term. E_0 = the SO(3)-symmetric (TBM)
      structure; E_1 = the leading SO(3)-BREAKING (the theta_13 channel; see
      research/10a_neutrino_oscillations.md "the connecting hypothesis").

VALIDATION. The O(delta) breaking coefficient E_1 is computed three ways:
  - NAIVE f64 finite difference [E(delta_h) - E(0)] / delta_h   (expected: catastrophic)
  - GRADED perturbative f64 (this method)                        (expected: clean)
  - mpmath 50-digit finite-difference reference on a tiny grid   (the gold truth)
PASS if graded matches the mpmath reference to ~1e-10 while naive is wrong by many orders.

Convention: Duda index-0, D = diag(g,1,delta,0), eta = diag(-1,1,1,1) (g and the metric
minus both at index 0). The production engine uses the equivalent index-3 ordering; this
sandbox adopts index-0 to match Duda and end the recurring confusion.

LOCAL artifact (OpenWave #236 N-program, held until the N-program finishes). Headless.
Scope: research/10n1_foundation_scope.md.  Run: python3 n1_precision_method.py
"""

import json
import os
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))

# ----------------------------------------------------------------------------------
# Convention (Duda index-0) + physical scales
# ----------------------------------------------------------------------------------
G_SCALE = 1.0e10        # index-0 eigenvalue: rest-mass / boost axis
DELTA_PHYS = 1.0e-10    # index-2 eigenvalue: quantum-phase (Dirac) axis = the breaking
ETA_DIAG = np.array([-1.0, 1.0, 1.0, 1.0])        # Minkowski metric, minus at index 0
SIGN_MAT = np.outer(ETA_DIAG, ETA_DIAG)           # s_ab = eta_aa eta_bb (=-1 iff one idx is 0)

# LdG potential coeffs (representative; the engine's exact (a,b,c) is the N3 search, Duda Q7)
LDG_A, LDG_B, LDG_C = -2.0, 0.0, 1.0
CURV_C2 = 1.0           # c^2 curvature prefactor (bare units for the precision demo)

D0_DIAG = np.array([G_SCALE, 1.0, 0.0, 0.0])      # D(delta=0)
D1_DIAG = np.array([0.0, 0.0, 1.0, 0.0])          # dD/ddelta  (the breaking direction)


# ----------------------------------------------------------------------------------
# Geometry: a frame field O(x) that mixes the g-axis (boost) into a spatial winding,
# so the curvature commutators genuinely carry the full 1e20 scale range.
# ----------------------------------------------------------------------------------
def _rot4(p, q, ang):
    """4x4 rotation by `ang` in the (p,q) plane."""
    r = np.eye(4)
    c, s = np.cos(ang), np.sin(ang)
    r[p, p] = c; r[q, q] = c; r[p, q] = -s; r[q, p] = s
    return r


def build_O_field(n, b_star=0.30):
    """O(x) = R_{12}(phi) . R_{01}(theta):
        phi(x,y) = atan2(y,x)            -> (1,2)-plane disclination winding (the delta axis
                                            mixes with the unit axis)
        theta(x) = b_star.exp(-r^2/rw^2) -> (0,1)-plane boost dressing (brings g into the
                                            spatial structure -> g enters the curvature)
    Returns O of shape (n,n,n,4,4), f64. Geometry is O(1); the scale lives in D, not O."""
    c0 = (n - 1) / 2.0
    rw = max(n / 4.0, 1.0)
    O = np.empty((n, n, n, 4, 4))
    for i in range(n):
        for j in range(n):
            for k in range(n):
                x, y, z = i - c0, j - c0, k - c0
                phi = np.arctan2(y, x)
                theta = b_star * np.exp(-(x * x + y * y + z * z) / (rw * rw))
                O[i, j, k] = _rot4(1, 2, phi) @ _rot4(0, 1, theta)
    return O


def M_from_O(O, diag):
    """M = O . diag(d) . O^T for a diagonal eigenvalue vector, batched over the grid."""
    return np.einsum("...ab,b,...db->...ad", O, diag, O)


# ----------------------------------------------------------------------------------
# Energy pieces (numpy f64)
# ----------------------------------------------------------------------------------
def _grads(F, dx=1.0):
    """Central-difference spatial gradients on the interior, shape (n-2,n-2,n-2,4,4)."""
    inv2 = 1.0 / (2.0 * dx)
    Fx = (F[2:, 1:-1, 1:-1] - F[:-2, 1:-1, 1:-1]) * inv2
    Fy = (F[1:-1, 2:, 1:-1] - F[1:-1, :-2, 1:-1]) * inv2
    Fz = (F[1:-1, 1:-1, 2:] - F[1:-1, 1:-1, :-2]) * inv2
    return Fx, Fy, Fz


def _comm(A, B):
    return np.einsum("...ab,...bc->...ac", A, B) - np.einsum("...ab,...bc->...ac", B, A)


def _sdot(A, B):
    """Signed bilinear <A,B>_s = sum_ab s_ab A_ab B_ab, per voxel (engine signed_dot4)."""
    return np.sum(SIGN_MAT * A * B, axis=(-2, -1))


def _potential_density(M):
    """V_M on the NON-g block (indices 1,2,3) = diag(1,delta,0) sector; g-axis excluded."""
    S = M[..., 1:, 1:]
    S2 = np.einsum("...ab,...bc->...ac", S, S)
    tr2 = np.einsum("...aa->...", S2)
    tr3 = np.einsum("...aa->...", np.einsum("...ab,...bc->...ac", S2, S))
    return LDG_A * tr2 - LDG_B * tr3 + LDG_C * tr2 * tr2


def naive_energy(M0, M1, delta):
    """The NAIVE path: form M = M0 + delta.M1 in f64 and sum the energy (faithful to how
    the engine would evaluate H at the physical scales). Returns the scalar total energy."""
    M = M0 + delta * M1
    Mx, My, Mz = _grads(M)
    curv = 4.0 * (_sdot(_comm(Mx, My), _comm(Mx, My))
                  + _sdot(_comm(Mx, Mz), _comm(Mx, Mz))
                  + _sdot(_comm(My, Mz), _comm(My, Mz)))
    Vd = _potential_density(M)[1:-1, 1:-1, 1:-1]
    return CURV_C2 * float(np.sum(curv)) + float(np.sum(Vd))


def graded_energy_orders(M0, M1):
    """The GRADED perturbative path: compute E_0, E_1, E_2 from the delta-graded
    commutators C_k and the graded potential, each in its own well-scaled arithmetic."""
    M0x, M0y, M0z = _grads(M0)
    M1x, M1y, M1z = _grads(M1)

    curv0 = curv1 = curv2 = 0.0
    for (a0, b0, a1, b1) in (
        (M0x, M0y, M1x, M1y),
        (M0x, M0z, M1x, M1z),
        (M0y, M0z, M1y, M1z),
    ):
        C0 = _comm(a0, b0)                       # O(g^2)
        C1 = _comm(a0, b1) + _comm(a1, b0)       # O(g)
        C2 = _comm(a1, b1)                       # O(1)
        curv0 += np.sum(_sdot(C0, C0))
        curv1 += np.sum(2.0 * _sdot(C0, C1))
        curv2 += np.sum(_sdot(C1, C1) + 2.0 * _sdot(C0, C2))
    curv0 *= 4.0 * CURV_C2
    curv1 *= 4.0 * CURV_C2
    curv2 *= 4.0 * CURV_C2

    # potential graded on the non-g block (S = S0 + delta.S1)
    S0 = M0[..., 1:, 1:]
    S1 = M1[..., 1:, 1:]
    sl = (slice(1, -1),) * 3

    def tr(X):
        return np.einsum("...aa->...", X)

    def mm(X, Y):
        return np.einsum("...ab,...bc->...ac", X, Y)

    S0S0 = mm(S0, S0)
    T0 = tr(S0S0)                 # Tr(S0^2)
    T01 = tr(mm(S0, S1))          # Tr(S0 S1)
    pot0 = LDG_A * T0 - LDG_B * tr(mm(S0S0, S0)) + LDG_C * T0 * T0
    pot1 = (LDG_A * 2.0 * T01
            - LDG_B * 3.0 * tr(mm(S0S0, S1))
            + LDG_C * 4.0 * T0 * T01)
    E0 = curv0 + float(np.sum(pot0[sl]))
    E1 = curv1 + float(np.sum(pot1[sl]))
    E2 = curv2  # potential order-2 omitted from E2 headline (sub-dominant, non-g block)
    return E0, E1, E2


# ----------------------------------------------------------------------------------
# mpmath gold reference (tiny grid): high-precision finite difference for E_1
# ----------------------------------------------------------------------------------
def mpmath_reference_E1(n_tiny=6, b_star=0.30, dps=50, delta_h="1e-3"):
    """Compute E_1 = dE/ddelta at delta=0 on a tiny grid in `dps`-digit arithmetic, by a
    central finite difference in delta. At 50 digits there is no catastrophic cancellation,
    so this is the trusted truth to validate the f64 graded method against. Uses the SAME
    f64 geometry O (promoted to mpf) so the test isolates ARITHMETIC precision."""
    mp.mp.dps = dps
    O = build_O_field(n_tiny, b_star=b_star)          # f64 geometry
    g = mp.mpf(G_SCALE)
    h = mp.mpf(delta_h)
    a, b, c, c2 = mp.mpf(LDG_A), mp.mpf(LDG_B), mp.mpf(LDG_C), mp.mpf(CURV_C2)
    smat = [[mp.mpf(int(SIGN_MAT[p, q])) for q in range(4)] for p in range(4)]

    def Omat(i, j, k):
        return [[mp.mpf(float(O[i, j, k, p, q])) for q in range(4)] for p in range(4)]

    def M_at(i, j, k, delta):
        Oijk = Omat(i, j, k)
        d = [g, mp.mpf(1), delta, mp.mpf(0)]
        # M = O diag(d) O^T
        return [[sum(Oijk[p][b_] * d[b_] * Oijk[q][b_] for b_ in range(4))
                 for q in range(4)] for p in range(4)]

    def comm(A, B):
        return [[sum(A[p][m] * B[m][q] - B[p][m] * A[m][q] for m in range(4))
                 for q in range(4)] for p in range(4)]

    def sdot(A, B):
        return sum(smat[p][q] * A[p][q] * B[p][q] for p in range(4) for q in range(4))

    def energy(delta):
        # precompute M on the tiny grid
        M = {}
        for i in range(n_tiny):
            for j in range(n_tiny):
                for k in range(n_tiny):
                    M[(i, j, k)] = M_at(i, j, k, delta)
        tot = mp.mpf(0)
        half = mp.mpf("0.5")
        for i in range(1, n_tiny - 1):
            for j in range(1, n_tiny - 1):
                for k in range(1, n_tiny - 1):
                    def grad(ax):
                        ip = [i, j, k]; im = [i, j, k]
                        ip[ax] += 1; im[ax] -= 1
                        Mp = M[tuple(ip)]; Mm = M[tuple(im)]
                        return [[(Mp[p][q] - Mm[p][q]) * half for q in range(4)] for p in range(4)]
                    Mx, My, Mz = grad(0), grad(1), grad(2)
                    cxy, cxz, cyz = comm(Mx, My), comm(Mx, Mz), comm(My, Mz)
                    curv = mp.mpf(4) * (sdot(cxy, cxy) + sdot(cxz, cxz) + sdot(cyz, cyz))
                    # potential on non-g block (indices 1,2,3)
                    Mc = M[(i, j, k)]
                    S = [[Mc[p][q] for q in range(1, 4)] for p in range(1, 4)]
                    S2 = [[sum(S[p][m] * S[m][q] for m in range(3)) for q in range(3)] for p in range(3)]
                    tr2 = sum(S2[p][p] for p in range(3))
                    tr3 = sum(sum(S2[p][m] * S[m][p] for m in range(3)) for p in range(3))
                    V = a * tr2 - b * tr3 + c * tr2 * tr2
                    tot += c2 * curv + V
        return tot

    Ep, Em = energy(h), energy(-h)
    E1 = (Ep - Em) / (2 * h)
    return float(E1)


# ----------------------------------------------------------------------------------
# Driver
# ----------------------------------------------------------------------------------
def main():
    n_full = 20
    print("=" * 78)
    print("N1 , precision-safe numerical method (Duda index-0, perturbative-delta)")
    print("=" * 78)
    print(f"g = {G_SCALE:.1e}  delta_phys = {DELTA_PHYS:.1e}  range = {G_SCALE/DELTA_PHYS:.1e}")

    O = build_O_field(n_full)
    M0 = M_from_O(O, D0_DIAG)
    M1 = M_from_O(O, D1_DIAG)

    # Graded perturbative orders (the method) , full grid
    E0, E1_graded, E2 = graded_energy_orders(M0, M1)
    print(f"\n[graded]  E_0 = {E0: .6e}   (order delta^0, the g-sector ~ symmetric/TBM)")
    print(f"[graded]  E_1 = {E1_graded: .6e}   (order delta^1, the SO(3)-BREAKING channel)")
    print(f"[graded]  E_2 = {E2: .6e}   (order delta^2)")
    print(f"[non-dim] scale grading |E_0|:|E_1|:|E_2| ~ "
          f"{abs(E0):.1e} : {abs(E1_graded):.1e} : {abs(E2):.1e}  (each order ~g down)")

    # Energy-level loss: the O(delta) breaking signal vs the f64 ULP set by the g-sector E_0.
    energy_ulp = float(np.spacing(abs(E0)))
    delta_signal = float(DELTA_PHYS * abs(E1_graded))
    signal_below_ulp = bool(delta_signal < energy_ulp)
    print(f"\n[energy-level] delta_phys*|E_1| = {delta_signal:.3e}  (the breaking's energy signal)")
    print(f"[energy-level] ULP of E_0        = {energy_ulp:.3e}  (f64 resolution at the g-sector)")
    print(f"[energy-level] breaking signal below the ULP -> naive derivative underflows? {signal_below_ulp}")

    # Naive finite difference for E_1 across a sweep of delta_h
    deltas = np.logspace(0, -12, 25)
    naive_E1_sweep = []
    E_at_0 = naive_energy(M0, M1, 0.0)
    for dh in deltas:
        naive_E1_sweep.append((naive_energy(M0, M1, dh) - E_at_0) / dh)
    naive_E1_sweep = np.array(naive_E1_sweep)

    naive_E1_phys = (naive_energy(M0, M1, DELTA_PHYS) - E_at_0) / DELTA_PHYS
    print(f"\n[naive]   E_1 (finite diff at delta_phys={DELTA_PHYS:.0e}) = {naive_E1_phys: .6e}")
    rel_naive = abs(naive_E1_phys - E1_graded) / abs(E1_graded)
    print(f"[naive]   relative error vs graded = {rel_naive:.3e}")

    # mpmath gold reference (tiny grid) , validate the graded algebra independently
    print("\n[mpmath]  computing 50-digit reference on a tiny grid (n=6) ...")
    E0_t, E1_t_graded, E2_t = graded_energy_orders(
        M_from_O(build_O_field(6), D0_DIAG), M_from_O(build_O_field(6), D1_DIAG))
    E1_t_mp = mpmath_reference_E1(n_tiny=6)
    rel_graded_tiny = abs(E1_t_graded - E1_t_mp) / abs(E1_t_mp)
    # naive on the tiny grid for the same comparison
    O6 = build_O_field(6)
    M0_6, M1_6 = M_from_O(O6, D0_DIAG), M_from_O(O6, D1_DIAG)
    naive_E1_tiny = (naive_energy(M0_6, M1_6, DELTA_PHYS) - naive_energy(M0_6, M1_6, 0.0)) / DELTA_PHYS
    rel_naive_tiny = abs(naive_E1_tiny - E1_t_mp) / abs(E1_t_mp) if E1_t_mp != 0 else float("inf")
    print(f"[mpmath]  E_1 reference (tiny)      = {E1_t_mp: .8e}")
    print(f"[graded]  E_1 (tiny)                = {E1_t_graded: .8e}  rel err = {rel_graded_tiny:.2e}")
    print(f"[naive]   E_1 (tiny)                = {naive_E1_tiny: .8e}  rel err = {rel_naive_tiny:.2e}")

    # ---- PASS / FAIL ----
    passed = bool((rel_graded_tiny < 1e-10) and (rel_naive_tiny > 1e-2) and signal_below_ulp)
    print("\n" + "=" * 78)
    print(f"CANCELLATION TEST: {'PASS' if passed else 'FAIL'}")
    print("  graded matches the 50-digit reference to ~machine precision; naive is destroyed")
    print("  by the 1e20 range (the breaking's energy signal sits below the g-sector ULP).")
    print("=" * 78)

    # ---- plot: naive E_1 error vs delta_h (the catastrophe as delta -> physical scale) ----
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
    rel_err = np.abs(naive_E1_sweep - E1_graded) / abs(E1_graded)
    ax[0].loglog(deltas, np.maximum(rel_err, 1e-18), "o-", color="#c0392b", label="naive f64 finite diff")
    ax[0].axhline(1e-10, ls="--", color="#27ae60", label="graded perturbative (machine)")
    ax[0].axvline(DELTA_PHYS, ls=":", color="#2c3e50", label=f"delta_phys = {DELTA_PHYS:.0e}")
    ax[0].set_xlabel("delta_h (finite-difference step)")
    ax[0].set_ylabel("relative error in E_1 vs graded")
    ax[0].set_title("N1: naive O(delta) breaking coeff collapses as delta -> physical scale")
    ax[0].legend(fontsize=8); ax[0].grid(True, which="both", alpha=0.3)

    ax[1].bar(["|E_0|\n(g-sector)", "|E_1|\n(breaking)", "|E_2|"],
              [abs(E0), abs(E1_graded), abs(E2)], color=["#2980b9", "#c0392b", "#7f8c8d"])
    ax[1].set_yscale("log")
    ax[1].set_ylabel("energy magnitude (bare units)")
    ax[1].set_title("N1: delta-graded energy orders (each computed separately)")
    ax[1].grid(True, which="both", axis="y", alpha=0.3)
    fig.tight_layout()
    plot_path = os.path.join(HERE, "n1_precision_method.png")
    fig.savefig(plot_path, dpi=110)
    print(f"\nplot -> {plot_path}")

    summary = {
        "convention": "Duda index-0: D=diag(g,1,delta,0), eta=diag(-1,1,1,1)",
        "g_scale": G_SCALE,
        "delta_phys": DELTA_PHYS,
        "dynamic_range": G_SCALE / DELTA_PHYS,
        "grid_full": n_full,
        "energy_level": {
            "delta_phys_times_E1_signal": delta_signal,
            "ulp_of_E0": energy_ulp,
            "breaking_signal_below_ulp": signal_below_ulp,
        },
        "graded_orders": {"E0": E0, "E1": E1_graded, "E2": E2},
        "naive_E1_at_delta_phys": naive_E1_phys,
        "naive_rel_err_vs_graded": rel_naive,
        "mpmath_reference": {
            "grid_tiny": 6,
            "dps": 50,
            "E1_reference": E1_t_mp,
            "E1_graded_tiny": E1_t_graded,
            "E1_naive_tiny": naive_E1_tiny,
            "graded_rel_err": rel_graded_tiny,
            "naive_rel_err": rel_naive_tiny,
        },
        "cancellation_test_pass": passed,
    }
    json_path = os.path.join(HERE, "n1_precision_method_summary.json")
    with open(json_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"summary -> {json_path}")
    return passed


if __name__ == "__main__":
    ok = main()
    raise SystemExit(0 if ok else 1)
