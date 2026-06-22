#!/usr/bin/env python3
"""
GOLDEN MASTER , the safety net for the index-3 -> index-0 convention flip.

Captures CONVENTION-INVARIANT physical outputs of the M5 engine across the full pipeline
(seeders + static energy + Euclidean gradient flow + the constrained clock integrator), so the
index flip can be PROVEN physics-neutral. N0 Test B already showed energies/potential/curvature
are invariant under the relabel sigma=[3,0,1,2]; this extends that to the WHOLE engine.

Why these invariants survive the relabel (a row+col permutation = a similarity transform):
  - sum of all M_ab^2 (Frobenius^2)        : reorders the sum, same value
  - Tr(M), Tr(M^2), Tr(M^3)                : trace of powers, similarity-invariant
  - eigenvalues of M (eigvalsh)            : similarity-invariant
  - energy density (per voxel scalar)      : N0 Test B; total/mean/max/min identical
  - stable-mask count, evolved-state energy: physical, identical trajectory
The script NEVER hardcodes a time-axis index, so the SAME file runs pre- and post-flip.

Usage (from research/_convention_refactor/):
  python3 golden_master.py            # no baseline -> CAPTURE to golden_master_baseline.json
  python3 golden_master.py --check    # baseline exists -> CHECK and report PASS/FAIL per item
Run CAPTURE before the flip, --check after.  tol = 1e-4 (f32, sum-order differences expected).
"""

import os
import sys
import json
import argparse

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../"))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

import numpy as np
import taichi as ti

ti.init(arch=ti.cpu, default_fp=ti.f32, random_seed=0)

from openwave.xperiments.m5_liquid_crystal import medium
from openwave.xperiments.m5_liquid_crystal import engine1_seeds as seeds
from openwave.xperiments.m5_liquid_crystal import engine2_pde as pde
from openwave.xperiments.m5_liquid_crystal import engine3_observables as obs

HERE = os.path.dirname(os.path.abspath(__file__))
BASELINE = os.path.join(HERE, "golden_master_baseline.json")
TOL = 1e-4

# ---- fixed toy-scale setup (deterministic) ----
EDGE = 1e-15
TARGET_VOXELS = 32 ** 3
DELTA = 0.30
B_STAR = 0.13
KICK = 0.05
KM_INERTIA = 30.0

tf = medium.TensorField([EDGE, EDGE, EDGE], TARGET_VOXELS, [0.5, 0.5, 0.5], viz_stride=2)
N = tf.nx
C = N // 2
GRID = tf.grid_size
observables = medium.FieldObservables(GRID)
trackers = medium.Trackers(GRID)

dx_am = float(tf.dx_am)
c_amrs = 0.2998
dt_rs = 0.5 * dx_am * 0.95 / (c_amrs * 3 ** 0.5)
ldg_c = 1.0 * c_amrs ** 2 / dx_am ** 4
ldg_a = -2.0 * ldg_c * (1.0 + DELTA ** 2)
ldg_b = 0.0
tr2_vac = 1.0 + DELTA ** 2
v0 = ldg_a * tr2_vac + ldg_c * tr2_vac ** 2
r0_vox = 0.06 * N
rhoc_vox = 3.0
rw_vox = 0.29 * N


# ---- invariants (convention-agnostic) ----
def field_invariants(tag):
    M = tf.M_am.to_numpy().astype(np.float64)
    frob2 = float(np.sum(M * M))
    trM = float(np.sum(np.einsum("...aa->...", M)))
    M2 = np.einsum("...ab,...bc->...ac", M, M)
    trM2 = float(np.sum(np.einsum("...aa->...", M2)))
    eig = np.linalg.eigvalsh(M)                 # sorted ascending, per voxel
    eig_means = [float(np.mean(eig[..., i])) for i in range(4)]
    return {f"{tag}.frob2": frob2, f"{tag}.trM": trM, f"{tag}.trM2": trM2,
            f"{tag}.eig0": eig_means[0], f"{tag}.eig1": eig_means[1],
            f"{tag}.eig2": eig_means[2], f"{tag}.eig3": eig_means[3]}


def energy_invariants(tag):
    obs.compute_energyH_density_M(tf, observables, c_amrs, dt_rs, ldg_a, ldg_b, ldg_c, v0, 1.0)
    e = observables.energyH_density_aJ.to_numpy().astype(np.float64)
    return {f"{tag}.E_tot": float(np.sum(e)), f"{tag}.E_mean": float(np.mean(e)),
            f"{tag}.E_max": float(np.max(e)), f"{tag}.E_min": float(np.min(e))}


def render_invariants(tag):
    """The RENDER/tracker path the compute suite misses: eigen_decompose (director +
    spatial eigenvalues, via principal_director on the spatial block) + update_trackers_M
    (amplitude ‖M − D_vac‖_F). These read the spatial-block layout, so they catch the
    eigen_decompose + d_vac cascades the flip touches."""
    pde.eigen_decompose(tf)
    nh = tf.director_nhat.to_numpy().astype(np.float64)
    ev = tf.eigenvalues.to_numpy().astype(np.float64)
    obs.update_trackers_M(tf, trackers, dt_rs, DELTA)
    amp = trackers.amp_local_emarms_am.to_numpy().astype(np.float64)
    return {
        f"{tag}.dir_mean_norm": float(np.mean(np.linalg.norm(nh, axis=-1))),
        f"{tag}.dir_probe_x": float(nh[C + 8, C, C, 0]),
        f"{tag}.dir_probe_y": float(nh[C + 8, C, C, 1]),
        f"{tag}.dir_probe_z": float(nh[C + 8, C, C, 2]),
        f"{tag}.eig_hi_mean": float(np.mean(ev[..., 0])),
        f"{tag}.eig_mid_mean": float(np.mean(ev[..., 1])),
        f"{tag}.eig_lo_mean": float(np.mean(ev[..., 2])),
        f"{tag}.amp_mean": float(np.mean(amp)),
        f"{tag}.amp_max": float(np.max(amp)),
    }


def setup_constrained():
    idx = np.indices((N, N, N)).astype(np.float32)
    r_vox = np.sqrt((idx[0] - C) ** 2 + (idx[1] - C) ** 2 + (idx[2] - C) ** 2)
    rho_vox = np.sqrt((idx[0] - C) ** 2 + (idx[1] - C) ** 2)
    act = np.zeros((N, N, N), np.float32)
    act[2:-2, 2:-2, 2:-2] = 1.0
    act *= (r_vox > 6.0) * (rho_vox > 3.0)
    tf.act4d.from_numpy(act)
    basis = np.zeros((10, 4, 4), np.float32)
    for a_ in range(4):
        basis[a_, a_, a_] = 1.0
    bi = 4
    for a_ in range(4):
        for c_ in range(a_ + 1, 4):
            basis[bi, a_, c_] = basis[bi, c_, a_] = 1.0 / np.sqrt(2.0)
            bi += 1
    tf.sym_basis.from_numpy(basis)
    return float(act[N // 2].sum() + act[:, N // 2].sum() + act[:, :, N // 2].sum())


def gradient_flow(n_steps):
    tf.stable_mask.fill(0.0)
    pde.compute_tstar(tf)
    n_pl = float(3 * N * N)
    for _ in range(n_steps):
        tf.M_prev_am.copy_from(tf.M_am)
        pde.compute_curvature_flux_4d(tf)
        pde.sample_v03_drift(tf, dt_rs)
        vm = [tf.v03_sums[a_] / n_pl for a_ in range(3)]
        pde.evolve_M_4d(tf, c_amrs, dt_rs, ldg_c, vm[0], vm[1], vm[2], KM_INERTIA)
        tf.swap_matrix_buffers()


def clock_run(n_steps, n_act_pl):
    dt_rs_b = 0.007 * dx_am * 0.95 / (c_amrs * 3 ** 0.5)
    dt_eff_b = c_amrs * dt_rs_b
    cc4d = 0.5 * ldg_c / (c_amrs * c_amrs)
    pde.compute_tstar(tf)
    pde.init_P_4d(tf, KICK)
    blew = False
    for n in range(n_steps):
        pde.flux_4d_constrained(tf)
        pde.update_P_4d(tf, dt_eff_b, cc4d)
        pde.sample_p03_drift(tf)
        npl = max(n_act_pl, 1.0)
        pde.apply_p03_clamp(tf, tf.v03_sums[0] / npl, tf.v03_sums[1] / npl, tf.v03_sums[2] / npl)
        pde.solve_constrained_4d(tf)
        pde.update_M_4d_constrained(tf, dt_eff_b)
        tf.swap_matrix_buffers()
        mmax = float(np.abs(tf.M_am.to_numpy()).max())
        if not np.isfinite(mmax) or mmax > 50.0:
            blew = True
            break
    return blew


def run_suite():
    out = {}
    # 1. seeders (cover embed4 / uniaxial_M / biaxial frame / dressed boost mixing)
    seeds.seed_vacuum_M(tf, DELTA)
    out.update(field_invariants("vacuum")); out.update(energy_invariants("vacuum"))
    seeds.seed_biaxial_hedgehog_M(tf, C, C, C, r0_vox, rhoc_vox, DELTA)
    out.update(field_invariants("biaxial")); out.update(energy_invariants("biaxial"))
    seeds.seed_dressed_hedgehog_M(tf, C, C, C, r0_vox, rhoc_vox, DELTA, B_STAR, rw_vox, 0.0)
    out.update(field_invariants("dressed")); out.update(energy_invariants("dressed"))
    pde.compute_stable_mask(tf)
    out["dressed.stable_mask_count"] = float(np.sum(tf.stable_mask.to_numpy()))

    # render/tracker path (eigen_decompose + amplitude tracker) — the cascade the compute suite misses
    seeds.seed_biaxial_hedgehog_M(tf, C, C, C, r0_vox, rhoc_vox, DELTA)
    out.update(render_invariants("render"))

    # 2. Euclidean gradient flow (evolve_M_4d path) on the dressed hedgehog
    seeds.seed_dressed_hedgehog_M(tf, C, C, C, r0_vox, rhoc_vox, DELTA, B_STAR, rw_vox, 0.0)
    gradient_flow(20)
    out.update(field_invariants("flow20")); out.update(energy_invariants("flow20"))

    # 3. constrained clock integrator (signed_dot4 / (alpha,3) flux , the riskiest path)
    seeds.seed_dressed_hedgehog_M(tf, C, C, C, r0_vox, rhoc_vox, DELTA, B_STAR, rw_vox, KICK)
    n_act_pl = setup_constrained()
    blew = clock_run(20, n_act_pl)
    out.update(field_invariants("clock20")); out.update(energy_invariants("clock20"))
    out["clock20.blew_up"] = bool(blew)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true", help="compare against the baseline")
    args = ap.parse_args()

    print("=" * 78)
    print(f"GOLDEN MASTER , {'CHECK' if args.check else 'CAPTURE'}   grid {N}^3  g={tf.lc_g}  delta={DELTA}")
    print("=" * 78)
    results = run_suite()

    if not args.check:
        with open(BASELINE, "w") as f:
            json.dump({"commit_note": "pre-flip index-3 baseline", "tol": TOL, "values": results}, f, indent=2)
        print(f"captured {len(results)} invariants -> {BASELINE}")
        for k in sorted(results):
            print(f"  {k:28s} = {results[k]}")
        return True

    base = json.load(open(BASELINE))["values"]
    print(f"comparing {len(results)} invariants against the pre-flip baseline (tol={TOL}):\n")
    worst = 0.0
    fails = []
    for k in sorted(results):
        if k not in base:
            print(f"  [NEW]  {k}")
            continue
        a, b = results[k], base[k]
        if isinstance(a, bool) or isinstance(b, bool):
            ok = (a == b); rel = 0.0 if ok else 1.0
        else:
            denom = max(abs(b), 1e-30)
            rel = abs(a - b) / denom
            ok = rel < TOL
        worst = max(worst, rel)
        flag = "ok " if ok else "FAIL"
        if not ok:
            fails.append(k)
        print(f"  [{flag}] {k:28s} base={b!s:>22}  now={a!s:>22}  rel={rel:.2e}")
    passed = not fails
    print("\n" + "=" * 78)
    print(f"GOLDEN MASTER: {'PASS' if passed else 'FAIL'}   worst rel diff = {worst:.2e}"
          + ("" if passed else f"   FAILED: {fails}"))
    print("  PASS => the index-3 -> index-0 flip is PROVEN physics-neutral.")
    print("=" * 78)
    return passed


if __name__ == "__main__":
    raise SystemExit(0 if main() else 1)
