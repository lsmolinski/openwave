# ⚠️ ARCHIVE (M5.8.1 4×4 promotion, 2026-06-04): written for the 3×3 matrix substrate.
# It builds 3×3 M arrays (M_am.from_numpy of diag(1,δ,0)), now a shape mismatch vs the
# 4×4 field block-diag(spatial, g). Kept as a historical record of the M5.4–M5.6 milestone;
# would need a 4×4 migration to run. Live 4×4 regression: sandbox_v8/m5_8_1_headless_check.py.
"""
M5.4 — 1/d Coulomb on the MATRIX substrate (headless gate)

THE M5.4 EXIT CRITERION. Reproduces the M5.1 Coulomb result on the new order-
parameter substrate M = O·D·O^T, confirming the static-topology physics carries
over after the Vector(3)→matrix migration (decided 2026-05-26: director-equivalent
path — isolates the substrate swap from any physics change).

Director-equivalent pipeline (per the confirmed M5.4 gate design):
    1. seed_hedgehog_M  → M_am (matrix) + director_nhat            [matrix seeder]
    2. eigen_decompose  → director_nhat ← principal eigenvector(M) [THE bridge]
    3. flow the extracted director into the validated M5.1 working buffer (psi_am)
    4. relax_director_step (UNCHANGED M5.1 kernel) → ground state
    5. compute_energyF_density (UNCHANGED M5.1 kernel) → E_total(d)
    6. fit E(d) = a + b/d, R²

Because steps 4–5 are the byte-identical M5.1 relaxation + Frank kernels operating
on the eigenvector extracted from M, the Coulomb result MUST carry over — this gate
verifies the matrix seeder + eigen_decompose pipeline feeds them a faithful director.

Config mirrors research/sandbox_v2/m5_1_coulomb.py (N=65, N_RELAX=2000, same
separations / blend / τ) so the numbers are directly comparable to M5.1's R²=0.978.

PASS: opposite-charge R² > 0.95 AND b < 0 (attractive).

USAGE:
    python -m openwave.xperiments.m5_liquid_crystal.research.sandbox_v4.m5_4_coulomb_matrix
"""

import sys
from pathlib import Path

import numpy as np
import taichi as ti

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[4]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from openwave.xperiments.m5_liquid_crystal import medium  # noqa: E402
from openwave.xperiments.m5_liquid_crystal import engine1_seeds as seeds  # noqa: E402
from openwave.xperiments.m5_liquid_crystal import engine2_pde as pde  # noqa: E402
from openwave.xperiments.m5_liquid_crystal import engine3_observables as observables  # noqa: E402


# ================================================================
# CONFIG (mirrors m5_1_coulomb.py for direct comparability)
# ================================================================
N = 65
UNIVERSE_EDGE_M = 1e-15
TARGET_VOXELS = N**3
N_RELAX = 2000
SEPARATIONS = [8, 10, 12, 14, 16, 18, 20]
DOMAIN_QUARTER_FRACTION = 0.20
TAU_FRAC = 0.4
R2_PASS_THRESHOLD = 0.95

PLOT_DIR = HERE / "plots"
PLOT_PATH = PLOT_DIR / "m5_4_coulomb_matrix.png"


# ================================================================
# HELPERS
# ================================================================


def setup_defect_pair(wf, d_vox, sign_pair):
    """Build centers + signs ti.fields for a hedgehog pair at offset ±d/2 along x."""
    centers = ti.field(dtype=ti.i32, shape=(2, 3))
    signs = ti.field(dtype=ti.i32, shape=(2,))
    half = d_vox // 2
    centers[0, 0] = wf.nx // 2 - half
    centers[0, 1] = wf.ny // 2
    centers[0, 2] = wf.nz // 2
    centers[1, 0] = wf.nx // 2 + half
    centers[1, 1] = wf.ny // 2
    centers[1, 2] = wf.nz // 2
    signs[0] = sign_pair[0]
    signs[1] = sign_pair[1]
    return centers, signs


def relax_and_measure_matrix(wf, obs, centers, signs, n_defects, n_steps, tau, dq, delta):
    """Matrix-substrate director-equivalent relaxation → total Frank energy.

    Seeds M, extracts the director from M via eigen_decompose, then relaxes the
    director directly — relax_director_step + compute_energyF_density operate on
    director_nhat (= principal eigenvector of M), the same gradient-descent math
    as M5.1, so the Coulomb result carries over.
    """
    # 1. seed the matrix order parameter (writes M_am + director_nhat)
    seeds.seed_hedgehog_M(wf, centers, signs, dq, n_defects, delta)
    # 2. extract the director from M (the bridge: director_nhat ← eigenvector(M_am))
    pde.eigen_decompose(wf)
    # 3. relax the director directly (M5.4 relax gradient-descends director_nhat)
    for _ in range(n_steps):
        pde.relax_director_step(wf, tau, centers, signs, n_defects)
        wf.director_nhat.copy_from(wf.director_nhat_new)
    # 4. Frank elastic energy of the relaxed director
    observables.compute_energyF_density(wf, obs, observables.K_FRANK)
    return float(obs.energyF_density_aJ.to_numpy().sum())


def fit_coulomb(d, E):
    """Fit E(d) = a + b/d via linear lstsq; return (a, b, R²)."""
    d, E = np.asarray(d, dtype=float), np.asarray(E, dtype=float)
    X = np.column_stack([np.ones_like(d), 1.0 / d])
    coef, *_ = np.linalg.lstsq(X, E, rcond=None)
    a, b = coef
    E_pred = a + b / d
    ss_res = ((E - E_pred) ** 2).sum()
    ss_tot = ((E - E.mean()) ** 2).sum()
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 1e-12 else 1.0
    return float(a), float(b), float(r2)


# ================================================================
# MAIN
# ================================================================


def main():
    ti.init(arch=ti.gpu)
    print("=" * 70)
    print("M5.4 — 1/d Coulomb on the MATRIX substrate (M = O·D·O^T)")
    print("=" * 70)

    wf = medium.WaveField([UNIVERSE_EDGE_M] * 3, TARGET_VOXELS)
    obs = medium.FieldObservables(wf.grid_size)
    delta = wf.lc_delta

    cfl_bound = (wf.dx_am**2) / 6.0
    tau = TAU_FRAC * cfl_bound
    dq = float(DOMAIN_QUARTER_FRACTION * max(wf.nx, wf.ny, wf.nz))

    print(f"Grid: {wf.nx}³  dx={wf.dx_am:.3f} am   δ={delta}")
    print(f"Relax steps per separation: {N_RELAX}  τ={tau:.3f}")
    print(f"Separations (voxels): {SEPARATIONS}   D/4 blend: {dq:.1f}")
    print()

    # --- Opposite-charge sweep (the headline test) ---
    print("[Opposite-charge pair (+1, −1) — expect ATTRACTIVE, b < 0]")
    E_opp = []
    for d in SEPARATIONS:
        centers, signs = setup_defect_pair(wf, d, (+1, -1))
        E_total = relax_and_measure_matrix(wf, obs, centers, signs, 2, N_RELAX, tau, dq, delta)
        E_opp.append(E_total)
        print(f"  d = {d:3d}  E = {E_total:.4e}")
    a_opp, b_opp, r2_opp = fit_coulomb(SEPARATIONS, E_opp)
    pass_opp = r2_opp > R2_PASS_THRESHOLD
    interaction_opp = "ATTRACTIVE" if b_opp < 0 else "REPULSIVE"
    print(f"\n  Fit: E(d) = {a_opp:.4e} + ({b_opp:.4e})/d")
    print(f"  R² = {r2_opp:.4f}  → {interaction_opp}")
    print(f"  Pass R² > {R2_PASS_THRESHOLD} → {'PASS' if pass_opp else 'FAIL'}")
    print()

    # --- Same-charge sweep (informational sanity check) ---
    print("[Same-charge pair (+1, +1) — expect REPULSIVE, b > 0]")
    E_same = []
    for d in SEPARATIONS:
        centers, signs = setup_defect_pair(wf, d, (+1, +1))
        E_total = relax_and_measure_matrix(wf, obs, centers, signs, 2, N_RELAX, tau, dq, delta)
        E_same.append(E_total)
        print(f"  d = {d:3d}  E = {E_total:.4e}")
    a_same, b_same, r2_same = fit_coulomb(SEPARATIONS, E_same)
    interaction_same = "ATTRACTIVE" if b_same < 0 else "REPULSIVE"
    print(f"\n  Fit: E(d) = {a_same:.4e} + ({b_same:.4e})/d   R² = {r2_same:.4f} → {interaction_same}")
    print(f"  (informational only — pinned same-sign defects don't reach clean equilibrium)")
    print()

    # --- Plot ---
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8))
    d_arr = np.array(SEPARATIONS, dtype=float)
    d_plot = np.linspace(d_arr.min(), d_arr.max(), 200)

    ax = axes[0]
    ax.plot(d_arr, E_opp, "o", ms=8, label="opp-charge measured")
    ax.plot(d_plot, a_opp + b_opp / d_plot, "--", label=f"fit R²={r2_opp:.3f}")
    ax.plot(d_arr, E_same, "s", ms=6, alpha=0.6, label="same-charge measured")
    ax.plot(d_plot, a_same + b_same / d_plot, ":", label=f"fit R²={r2_same:.3f}")
    ax.set_xlabel("separation d (voxels)"); ax.set_ylabel("Frank energy E")
    ax.set_title(f"(a) E(d) on matrix substrate — opp: {interaction_opp}")
    ax.legend(); ax.grid(alpha=0.3)

    ax = axes[1]
    inv_d = 1.0 / d_arr
    ax.plot(inv_d, E_opp, "o", ms=8, label="opp-charge")
    ax.plot(1.0 / d_plot, a_opp + b_opp / d_plot, "--", label=f"linear (b={b_opp:.2e})")
    ax.plot(inv_d, E_same, "s", ms=6, alpha=0.6, label="same-charge")
    ax.set_xlabel("1/d"); ax.set_ylabel("Frank energy E")
    ax.set_title("(b) E vs 1/d — straight line ⇒ Coulomb")
    ax.legend(); ax.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(PLOT_PATH, dpi=110)
    plt.close()
    print(f"Plot saved to {PLOT_PATH}")

    sign_opp_ok = b_opp < 0
    all_pass = pass_opp and sign_opp_ok
    print("\n" + "=" * 70)
    print(f"  OPPOSITE CHARGES (+1, −1) — MATRIX-SUBSTRATE GATE:")
    print(f"    R² = {r2_opp:.4f}  (threshold > {R2_PASS_THRESHOLD})  →  {'PASS' if pass_opp else 'FAIL'}")
    print(f"    b  = {b_opp:+.3e}  ({interaction_opp})  →  {'PASS (correct sign)' if sign_opp_ok else 'FAIL (wrong sign)'}")
    print(f"  M5.1 reference (Vector(3) substrate): R² = 0.978, b < 0")
    print("=" * 70)
    print(f"OVERALL: {'PASS' if all_pass else 'FAIL'}  (M5.4 gate: Coulomb reproduced on M = O·D·O^T)")
    print("=" * 70)
    return 0 if all_pass else 1


if __name__ == "__main__":
    sys.exit(main())
