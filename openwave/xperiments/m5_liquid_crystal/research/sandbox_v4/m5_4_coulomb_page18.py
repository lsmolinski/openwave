# ✅ MIGRATED to the 4×4 substrate (2026-06-05) — RUNNABLE key-finding reproduction.
# The spatial 3×3 M = n⊗n now embeds as block-diag(M_spatial, g) per the M5.8.1 promotion;
# the constant-g time block adds 0 to gradients/commutators, so the page-18 physics check
# is UNCHANGED. (TensorField was WaveField pre-M5.8.)
# Re-validated on 4×4 (2026-06-05): R²=0.9959, b<0 ATTRACTIVE — PASS.
"""
M5.4 — Page-18 Mathematica Coulomb cross-validation (headless)

Companion to m5_4_coulomb_matrix.py. Where that gate relaxes a seeded director and
measures Frank energy, THIS ports Duda's analytical page-18 Mathematica check (4a §11)
— a DIFFERENT, complementary energy functional on the ANALYTICAL (closed-form, un-
relaxed) director:

    cos = 1 + (z−d)/√((z−d)²+r²) − (z+d)/√((z+d)²+r²)      (pair on z-axis at ±d)
    n   = {√(1−cos²)·x/r, √(1−cos²)·y/r, cos}
    M   = n⊗n                                               (pure projector, δ=0)
    H   = Σ_{i<j} ‖[∂_iM, ∂_jM]‖²_F                          (commutator curvature)
    Es  = ∫ H d³r
    Fit:  V(d) ≈ 1589.56 − 25.16/d                          (Duda's analytical Coulomb)

This validates the production matrix operators end-to-end: M from a director, central-
difference ∂_iM, and pde.commutator — the same [·,·] that builds Eq.19/20. The 1/d FORM
(R² high, b < 0 attractive) is the cross-val target; the absolute constant 25.16 is
units/domain-dependent (Duda's NIntegrate normalization), so we report the fit and the
1/d quality, not an exact-constant match.

USAGE:
    python -m openwave.xperiments.m5_liquid_crystal.research.sandbox_v4.m5_4_coulomb_page18
"""

import sys
from pathlib import Path

import numpy as np
import taichi as ti

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[4]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

ti.init(arch=ti.gpu, default_fp=ti.f32)

from openwave.xperiments.m5_liquid_crystal import medium  # noqa: E402
from openwave.xperiments.m5_liquid_crystal import engine2_pde as pde  # noqa: E402

N = 65
UNIVERSE_EDGE_M = 1e-15
HALF_SEPARATIONS = [6, 8, 10, 12, 14, 16]  # defects on z-axis at center ± d
PLOT_DIR = HERE / "plots"
PLOT_PATH = PLOT_DIR / "m5_4_coulomb_page18.png"


@ti.kernel
def compute_curvature_energy(tf: ti.template(), out: ti.template()):  # type: ignore
    """Per-voxel page-18 Hamiltonian density H = Σ_{i<j} ‖[∂_iM, ∂_jM]‖²_F.

    Central-differences M_am along each axis, commutes each pair via the production
    pde.commutator, sums the squared Frobenius norms. 1-cell halo (interior only).
    """
    nx, ny, nz = tf.nx, tf.ny, tf.nz
    inv_2dx = 1.0 / (2.0 * tf.dx_am)
    for i, j, k in ti.ndrange((1, nx - 1), (1, ny - 1), (1, nz - 1)):
        dMx = (tf.M_am[i + 1, j, k] - tf.M_am[i - 1, j, k]) * inv_2dx
        dMy = (tf.M_am[i, j + 1, k] - tf.M_am[i, j - 1, k]) * inv_2dx
        dMz = (tf.M_am[i, j, k + 1] - tf.M_am[i, j, k - 1]) * inv_2dx
        cxy = pde.commutator(dMx, dMy)
        cxz = pde.commutator(dMx, dMz)
        cyz = pde.commutator(dMy, dMz)
        out[i, j, k] = cxy.norm_sqr() + cxz.norm_sqr() + cyz.norm_sqr()


def analytical_M(nx, ny, nz, d):
    """Build M = n⊗n from Duda's analytical pair director (page-18 cos profile)."""
    c = np.array([nx / 2.0, ny / 2.0, nz / 2.0])
    ii, jj, kk = np.meshgrid(np.arange(nx), np.arange(ny), np.arange(nz), indexing="ij")
    x, y, z = ii - c[0], jj - c[1], kk - c[2]
    r = np.sqrt(x**2 + y**2) + 1e-6  # cylindrical radius (soft core on z-axis)
    cos = 1.0 + (z - d) / np.sqrt((z - d) ** 2 + r**2) - (z + d) / np.sqrt((z + d) ** 2 + r**2)
    cos = np.clip(cos, -1.0, 1.0)
    sin = np.sqrt(np.maximum(1.0 - cos**2, 0.0))
    n = np.stack([sin * x / r, sin * y / r, cos], axis=-1)
    n /= (np.linalg.norm(n, axis=-1, keepdims=True) + 1e-12)
    M = n[..., :, None] * n[..., None, :]  # n⊗n (δ=0 projector, per page-18)
    # M5.8.1 — embed in the 4×4 substrate: block-diag(spatial, g). The constant-g
    # time block has zero gradient, so the commutator curvature is unchanged.
    M4 = np.zeros(M.shape[:-2] + (4, 4), np.float32)
    M4[..., :3, :3] = M
    M4[..., 3, 3] = medium.LC_G
    return M4


def fit_coulomb(d, E):
    d, E = np.asarray(d, float), np.asarray(E, float)
    X = np.column_stack([np.ones_like(d), 1.0 / d])
    (a, b), *_ = np.linalg.lstsq(X, E, rcond=None)
    Ep = a + b / d
    r2 = 1.0 - ((E - Ep) ** 2).sum() / max(((E - E.mean()) ** 2).sum(), 1e-12)
    return float(a), float(b), float(r2)


def main():
    print("=" * 70)
    print("M5.4 — Page-18 analytical Coulomb cross-validation (commutator curvature)")
    print("=" * 70)
    tf = medium.TensorField([UNIVERSE_EDGE_M] * 3, N**3)
    out = ti.field(ti.f32, shape=tf.grid_size)
    print(f"Grid: {tf.nx}³  dx={tf.dx_am:.3f} am   H = Σ_(i<j) ‖[∂_iM,∂_jM]‖²")
    print(f"Half-separations d (defects at center±d): {HALF_SEPARATIONS}\n")

    E = []
    for d in HALF_SEPARATIONS:
        tf.M_am.from_numpy(analytical_M(tf.nx, tf.ny, tf.nz, d))
        out.fill(0)
        compute_curvature_energy(tf, out)
        Es = float(out.to_numpy().sum())
        E.append(Es)
        print(f"  d = {d:3d}   Es = {Es:.4e}")

    a, b, r2 = fit_coulomb(HALF_SEPARATIONS, E)
    interaction = "ATTRACTIVE" if b < 0 else "REPULSIVE"
    print(f"\n  Fit: V(d) = {a:.4e} + ({b:.4e})/d")
    print(f"  R² = {r2:.4f}  → {interaction}")
    print(f"  Duda page-18 analytical reference: V(d) ≈ 1589.56 − 25.16/d (clean 1/d)")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    d_arr = np.array(HALF_SEPARATIONS, float)
    dp = np.linspace(d_arr.min(), d_arr.max(), 200)
    _, ax = plt.subplots(1, 2, figsize=(13, 4.8))
    ax[0].plot(d_arr, E, "o", ms=8, label="commutator energy Es")
    ax[0].plot(dp, a + b / dp, "--", label=f"fit a+b/d  R²={r2:.3f}")
    ax[0].set_xlabel("half-separation d (voxels)"); ax[0].set_ylabel("Es = ∫ H")
    ax[0].set_title(f"(a) Es(d) — {interaction}"); ax[0].legend(); ax[0].grid(alpha=0.3)
    ax[1].plot(1.0 / d_arr, E, "o", ms=8); ax[1].plot(1.0 / dp, a + b / dp, "--")
    ax[1].set_xlabel("1/d"); ax[1].set_ylabel("Es"); ax[1].set_title("(b) Es vs 1/d ⇒ Coulomb")
    ax[1].grid(alpha=0.3)
    plt.tight_layout(); plt.savefig(PLOT_PATH, dpi=110); plt.close()
    print(f"  Plot saved to {PLOT_PATH}")

    ok = r2 > 0.95 and b < 0
    print("\n" + "=" * 70)
    print(f"  CROSS-VAL: 1/d form R² = {r2:.4f} (>0.95), b = {b:+.3e} ({interaction})")
    print(f"  {'PASS' if ok else 'FAIL'} — analytical commutator-curvature Coulomb reproduced on M")
    print("=" * 70)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
