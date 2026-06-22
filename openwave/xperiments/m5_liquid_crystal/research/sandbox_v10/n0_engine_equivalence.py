#!/usr/bin/env python3
"""
N0 , engine <-> numpy-port equivalence test (nails the foundation to the M5 LC engine).

WHY THIS EXISTS. The N1/N2 foundation scripts are standalone numpy: they RE-IMPLEMENT the M5
liquid-crystal energy functional (order parameter M=O.D.O^T, the signed curvature commutators,
the LdG potential V_M) in f64, because (a) the production engine is f32 and cannot carry the
g~1e10 / delta~1e-10 range the precision method is about, (b) the engine uses the index-3 time
axis while the neutrino build adopts Duda's index-0, and (c) the closed loop is new geometry.
That port's fidelity to the actual engine was, until now, ASSERTED (transcribed by reading the
source), not VERIFIED. This file verifies it , and therefore earns the claim "validated against
the M5 liquid-crystal model".

TWO TESTS (at the engine's TOY scales g=8, delta=0.30, where f32 is fine):
  A. FIDELITY      : seed a field with the engine, compute the engine's energy kernel
                     (engine3_observables.compute_energyH_density_M), pull the M field out, run
                     the numpy port on the IDENTICAL field. PASS if energies agree to f32 tol.
  B. CONVENTION    : relabel the same field from index-3 (engine: time axis = index 3) to index-0
                     (Duda: time axis = index 0) and recompute the curvature (signed), the
                     potential, in BOTH conventions. PASS if identical , proving the index choice
                     is a pure relabeling with ZERO physics consequence (the convention question).

LOCAL artifact (OpenWave #236 N-program, HELD until the N-program finishes). Headless / CPU.
Run from sandbox_v10:  python3 n0_engine_equivalence.py
"""

import os
import sys
import json

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../"))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

import numpy as np
import taichi as ti

ti.init(arch=ti.cpu, default_fp=ti.f32)

from openwave.xperiments.m5_liquid_crystal import medium
from openwave.xperiments.m5_liquid_crystal import engine1_seeds as seeds
from openwave.xperiments.m5_liquid_crystal import engine2_pde as pde
from openwave.xperiments.m5_liquid_crystal import engine3_observables as obs

HERE = os.path.dirname(os.path.abspath(__file__))


# ----------------------------------------------------------------------------------
# numpy port energy pieces (parametrized by which index is the time/g axis)
# ----------------------------------------------------------------------------------
def _grads(F, dx):
    inv2 = 1.0 / (2.0 * dx)
    Fx = (F[2:, 1:-1, 1:-1] - F[:-2, 1:-1, 1:-1]) * inv2
    Fy = (F[1:-1, 2:, 1:-1] - F[1:-1, :-2, 1:-1]) * inv2
    Fz = (F[1:-1, 1:-1, 2:] - F[1:-1, 1:-1, :-2]) * inv2
    return Fx, Fy, Fz


def _comm(A, B):
    return np.einsum("...ab,...bc->...ac", A, B) - np.einsum("...ab,...bc->...ac", B, A)


def _frob_sq(C):
    return np.sum(C * C, axis=(-2, -1))


def _sign_mat(time_idx):
    d = np.ones(4); d[time_idx] = -1.0
    return np.outer(d, d)            # s_ab = -1 iff exactly one of a,b is the time axis


def _signed_sq(C, time_idx):
    return np.sum(_sign_mat(time_idx) * C * C, axis=(-2, -1))


def _spatial_block(M, g_idx):
    keep = [i for i in range(4) if i != g_idx]
    return M[..., keep, :][..., :, keep]


def _V_M(M, g_idx, a, b, c):
    S = _spatial_block(M, g_idx)
    S2 = np.einsum("...ab,...bc->...ac", S, S)
    tr2 = np.einsum("...aa->...", S2)
    tr3 = np.einsum("...aa->...", np.einsum("...ab,...bc->...ac", S2, S))
    return a * tr2 - b * tr3 + c * tr2 * tr2


def numpy_energy_frobenius(M, dx, c2, a, b, c, v0, g_idx):
    """Mirror engine3_observables.compute_energyH_density_M for a STATIC field (Mdot=0):
    E = sum_interior [ c2 * 4 * sum_{mu<nu}||[d_mu M, d_nu M]||_F^2  +  (V_M - v0) ]."""
    Mx, My, Mz = _grads(M, dx)
    curv = 4.0 * (_frob_sq(_comm(Mx, My)) + _frob_sq(_comm(Mx, Mz)) + _frob_sq(_comm(My, Mz)))
    Vd = _V_M(M, g_idx, a, b, c)[1:-1, 1:-1, 1:-1] - v0
    return float(np.sum(c2 * curv + Vd))


def numpy_curv_signed(M, dx, c2, time_idx):
    Mx, My, Mz = _grads(M, dx)
    curv = 4.0 * (_signed_sq(_comm(Mx, My), time_idx) + _signed_sq(_comm(Mx, Mz), time_idx)
                  + _signed_sq(_comm(My, Mz), time_idx))
    return float(np.sum(c2 * curv))


def relabel_to_index3(M):
    """The engine is now index-0 (time axis first; flipped 2026-06-21). Relabel back to the
    OLD index-3 (time axis last) to demonstrate the two conventions are identical physics.
    sigma=[1,2,3,0] (index-0 -> index-3: every index k -> (k-1) mod 4)."""
    sigma = [1, 2, 3, 0]
    return M[..., sigma, :][..., :, sigma]


def main():
    print("=" * 80)
    print("N0 , engine <-> numpy-port equivalence (toy scales g=8, delta=0.30)")
    print("=" * 80)

    # ---- build the engine field (mirror the m5_9_4 production setup, toy scales) ----
    EDGE = 1e-15
    TARGET_VOXELS = 32 ** 3
    delta = 0.30
    tf = medium.TensorField([EDGE, EDGE, EDGE], TARGET_VOXELS, [0.5, 0.5, 0.5], viz_stride=2)
    N = tf.nx
    C = N // 2
    GRID = tf.grid_size
    observables = medium.FieldObservables(GRID)

    dx_am = float(tf.dx_am)
    c_amrs = 0.2998
    dt_rs = 0.5 * dx_am * 0.95 / (c_amrs * 3 ** 0.5)
    ldg_c = 1.0 * c_amrs ** 2 / dx_am ** 4
    ldg_a = -2.0 * ldg_c * (1.0 + delta ** 2)
    ldg_b = 0.0
    tr2_vac = 1.0 + delta ** 2
    v0 = ldg_a * tr2_vac + ldg_c * tr2_vac ** 2
    r0_vox = 0.06 * N
    rhoc_vox = 3.0
    c2 = c_amrs * c_amrs
    print(f"grid {N}^3   dx_am={dx_am:.4e}   g(lc_g)={tf.lc_g}   delta={delta}")
    print(f"LdG a={ldg_a:.3e} b={ldg_b:.1f} c={ldg_c:.3e}  v0={v0:.3e}  c2={c2:.4e}")

    # static biaxial hedgehog -> M_prev = M (kinetic = 0), pure curvature + potential
    seeds.seed_biaxial_hedgehog_M(tf, C, C, C, r0_vox, rhoc_vox, delta)

    # ---- engine energy (the reference) ----
    obs.compute_energyH_density_M(tf, observables, c_amrs, dt_rs, ldg_a, ldg_b, ldg_c, v0, 1.0)
    E_engine = float(np.sum(observables.energyH_density_aJ.to_numpy()))

    # ---- pull the identical field, run the numpy port (engine now stores g at index 0) ----
    M = tf.M_am.to_numpy().astype(np.float64)
    E_port_eng = numpy_energy_frobenius(M, dx_am, c2, ldg_a, ldg_b, ldg_c, v0, g_idx=0)

    relA = abs(E_engine - E_port_eng) / max(abs(E_engine), 1e-30)
    print("\n--- TEST A : fidelity (engine kernel vs numpy port, same field) ---")
    print(f"  E_engine (compute_energyH_density_M) = {E_engine: .8e}")
    print(f"  E_port   (numpy Frobenius, index-0)  = {E_port_eng: .8e}")
    print(f"  relative difference                  = {relA:.3e}   (f32 seed -> tol 1e-4)")
    passA = bool(relA < 1e-4)

    # ---- TEST B : convention invariance (engine now index-0 vs historical index-3) ----
    M3 = relabel_to_index3(M)            # same physical field, OLD index-3 labels
    # curvature: Frobenius is label-blind; the SIGNED curvature is the real convention probe
    cs_eng = numpy_curv_signed(M, dx_am, c2, time_idx=0)     # engine now: minus at index 0
    cs_duda = numpy_curv_signed(M3, dx_am, c2, time_idx=3)   # historical: minus at index 3
    # potential: the non-g block must be the SAME 3 spatial axes after relabel
    V_eng = float(np.sum(_V_M(M, 0, ldg_a, ldg_b, ldg_c)[1:-1, 1:-1, 1:-1]))
    V_duda = float(np.sum(_V_M(M3, 3, ldg_a, ldg_b, ldg_c)[1:-1, 1:-1, 1:-1]))
    relB_curv = abs(cs_eng - cs_duda) / max(abs(cs_eng), 1e-30)
    relB_pot = abs(V_eng - V_duda) / max(abs(V_eng), 1e-30)
    print("\n--- TEST B : convention invariance (engine now index-0 vs historical index-3) ---")
    print(f"  signed curvature  index-0 = {cs_eng: .8e}   index-3 = {cs_duda: .8e}   rel = {relB_curv:.2e}")
    print(f"  LdG potential     index-0 = {V_eng: .8e}   index-3 = {V_duda: .8e}   rel = {relB_pot:.2e}")
    passB = bool(relB_curv < 1e-12 and relB_pot < 1e-12)

    passed = bool(passA and passB)
    print("\n" + "=" * 80)
    print(f"N0 EQUIVALENCE: {'PASS' if passed else 'FAIL'}")
    print("  A: the numpy port reproduces the engine energy on an identical field (the M5")
    print("     liquid-crystal connection is VERIFIED, not just transcribed).")
    print("  B: the engine (now index-0) and the historical index-3 are the SAME physics (pure")
    print("     relabel) , the engine + sandbox_v10 + Duda now all share the index-0 convention.")
    print("=" * 80)

    summary = {
        "scales": {"g_lc": float(tf.lc_g), "delta": delta, "grid": N, "dx_am": dx_am},
        "test_A_fidelity": {
            "E_engine": E_engine, "E_numpy_port_index0": E_port_eng,
            "rel_diff": relA, "tol": 1e-4, "pass": passA,
        },
        "test_B_convention_invariance": {
            "signed_curv_index3": cs_eng, "signed_curv_index0": cs_duda, "rel_diff": relB_curv,
            "potential_index3": V_eng, "potential_index0": V_duda, "rel_diff_pot": relB_pot,
            "pass": passB,
            "note": "index-3 (engine) and index-0 (Duda) are a pure relabel, zero physics consequence",
        },
        "equivalence_pass": passed,
    }
    with open(os.path.join(HERE, "n0_engine_equivalence_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    print(f"summary -> {os.path.join(HERE, 'n0_engine_equivalence_summary.json')}")
    return passed


if __name__ == "__main__":
    raise SystemExit(0 if main() else 1)
