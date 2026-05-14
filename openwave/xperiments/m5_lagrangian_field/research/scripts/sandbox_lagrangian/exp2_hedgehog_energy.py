"""
Experiment 2: Hedgehog Energy vs Distance

Numerical Experiment 2 in the Lagrangian Framework Sub-Project.

Hypothesis:
Two hedgehog defects in a 3D director field produce a clean 1/r Coulomb potential
without sinc oscillation.

Setup:
- 3D grid with unit-vector director field n(x) = (nx, ny, nz), |n|=1
- Seed vacuum: n = ẑ everywhere
- Seed hedgehog at c1 (Q=+1): n(x) = (x-c1)/|x-c1| within a blending region
- Seed anti-hedgehog at c2 (Q=-1): n(x) = -(x-c2)/|x-c2|
- Relax with gradient descent on Frank elastic energy: H = (K/2)·∫|∇n|² d³r
- Keep defect cores fixed (Dirichlet-like BC on small ball around each center)
- Measure E(d) for separations d ∈ [2, 8]

Expected: E(d) ≈ const + C/d (see Duda arXiv:2108.07896 Fig. 2).

Spec:    ../3_LAGRANGIAN_FRAMEWORK.md  (Experiment 2)
Results: ../3b_lagrangian_experiments.md  (Experiment 2)
"""

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# PARAMETERS
# ---------------------------------------------------------------------------

N = 48                         # grid points per axis (48³ ≈ 110k voxels)
DOMAIN = 16.0                  # physical extent [-8, +8] each axis
DX = DOMAIN / (N - 1)
K_FRANK = 1.0                  # one-constant Frank elastic constant
CORE_RADIUS = 0.8              # hedgehog core size (only center point pinned)
# Stability: τ < dx²/(2·dim·K) — for dx=0.34 this is ≈ 0.019
TAU = 0.008                    # gradient-descent step size (conservative)
N_RELAX = 60                   # gentle relaxation (keep topology stable)
SEPARATIONS = [2.0, 3.0, 4.0, 5.0, 6.0, 8.0]

RESULTS_DIR = Path(__file__).parent / "exp2_results"


# ---------------------------------------------------------------------------
# GRID AND FIELDS
# ---------------------------------------------------------------------------

_X = _Y = _Z = None  # cached coordinate grids


def make_grid():
    """Return cached 3D meshgrid coordinates (X, Y, Z)."""
    global _X, _Y, _Z
    if _X is None:
        x = np.linspace(-DOMAIN / 2, DOMAIN / 2, N)
        _X, _Y, _Z = np.meshgrid(x, x, x, indexing="ij")
    return _X, _Y, _Z


def seed_vacuum():
    """Uniform director n(x) = ẑ everywhere (ground state)."""
    n = np.zeros((N, N, N, 3), dtype=np.float64)
    n[..., 2] = 1.0
    return n


def seed_hedgehog_pair(center1, sign1, center2, sign2):
    """Two hedgehogs via weighted superposition + unit-normalization.

    Returns (n, r1, r2) where r1, r2 are distance fields to each defect.
    """
    X, Y, Z = make_grid()
    dx1, dy1, dz1 = X - center1[0], Y - center1[1], Z - center1[2]
    dx2, dy2, dz2 = X - center2[0], Y - center2[1], Z - center2[2]
    r1 = np.sqrt(dx1**2 + dy1**2 + dz1**2)
    r2 = np.sqrt(dx2**2 + dy2**2 + dz2**2)
    r1_safe = np.maximum(r1, 0.2 * DX)
    r2_safe = np.maximum(r2, 0.2 * DX)

    # Each hedgehog's ideal director field
    n1 = np.stack([sign1 * dx1 / r1_safe,
                   sign1 * dy1 / r1_safe,
                   sign1 * dz1 / r1_safe], axis=-1)
    n2 = np.stack([sign2 * dx2 / r2_safe,
                   sign2 * dy2 / r2_safe,
                   sign2 * dz2 / r2_safe], axis=-1)

    # Proximity weighting — the closer defect dominates
    w1 = 1.0 / (r1 + 0.5)
    w2 = 1.0 / (r2 + 0.5)
    n = w1[..., None] * n1 + w2[..., None] * n2

    # Also blend in vacuum (ẑ) so far field returns to ground state
    # Scale of "far" relative to defects
    dist_to_nearest = np.minimum(r1, r2)
    w_vac = 1.0 / (1.0 + (dist_to_nearest / (DOMAIN / 4)) ** 4)  # ~1 near defects, →0 far
    vac = np.zeros_like(n)
    vac[..., 2] = 1.0
    n = w_vac[..., None] * n + (1.0 - w_vac[..., None]) * vac

    # Renormalize to unit length
    norm = np.linalg.norm(n, axis=-1, keepdims=True)
    n = n / np.maximum(norm, 1e-12)
    return n, r1, r2


def enforce_cores(n, center1, sign1, center2, sign2, r1, r2, core_radius):
    """Pin only the closest voxel to each defect center (soft Dirichlet).

    Prevents the hedgehog from numerically dissolving into vacuum during
    relaxation (since on a discrete grid, topology is not strictly preserved
    without pinning). Only one voxel per defect is pinned — this avoids the
    sharp-discontinuity-at-core-boundary problem a ball-shaped Dirichlet
    would create.
    """
    X, Y, Z = make_grid()
    # Single closest voxel per defect
    i1 = np.unravel_index(np.argmin(r1), r1.shape)
    i2 = np.unravel_index(np.argmin(r2), r2.shape)

    for (i, center, sign) in [(i1, center1, sign1), (i2, center2, sign2)]:
        # At the exact closest voxel, set n to a canonical direction
        # (for a hedgehog, which direction doesn't matter — the topology
        # is carried by the surrounding field; we just keep it from drifting)
        n[i[0], i[1], i[2], :] = np.array([0.0, 0.0, sign])  # ±ẑ
    return n


# ---------------------------------------------------------------------------
# ENERGY + RELAXATION
# ---------------------------------------------------------------------------

def frank_energy(n, dx):
    """H = (K/2) · ∫ |∇n|² d³r  (one-constant Frank elastic approximation)."""
    grad2 = 0.0
    for comp in range(3):
        for axis in range(3):
            g = np.gradient(n[..., comp], dx, axis=axis)
            grad2 = grad2 + g**2
    return 0.5 * K_FRANK * grad2.sum() * dx**3


def laplacian_vec(n, dx):
    """Vector Laplacian via 6-point stencil on each component (periodic BC via roll)."""
    lap = np.zeros_like(n)
    for comp in range(3):
        f = n[..., comp]
        lap[..., comp] = (
            np.roll(f, 1, axis=0) + np.roll(f, -1, axis=0)
            + np.roll(f, 1, axis=1) + np.roll(f, -1, axis=1)
            + np.roll(f, 1, axis=2) + np.roll(f, -1, axis=2)
            - 6.0 * f
        ) / dx**2
    return lap


def relax(n, center1, sign1, center2, sign2, n_steps=N_RELAX, tau=TAU,
          core_radius=CORE_RADIUS, verbose=False):
    """Gradient descent on Frank energy: ∂n/∂τ = ∇²n − (n·∇²n)n, renormalize.

    Defect cores held fixed as Dirichlet BC.
    """
    _, _, _ = make_grid()
    X, Y, Z = make_grid()
    dx1, dy1, dz1 = X - center1[0], Y - center1[1], Z - center1[2]
    dx2, dy2, dz2 = X - center2[0], Y - center2[1], Z - center2[2]
    r1 = np.sqrt(dx1**2 + dy1**2 + dz1**2)
    r2 = np.sqrt(dx2**2 + dy2**2 + dz2**2)

    energies = []
    for step in range(n_steps):
        lap = laplacian_vec(n, DX)
        ndot_lap = (n * lap).sum(axis=-1, keepdims=True)
        dn = lap - ndot_lap * n   # tangent projection
        n = n + tau * dn
        # Unit-length constraint
        norm = np.linalg.norm(n, axis=-1, keepdims=True)
        n = n / np.maximum(norm, 1e-12)
        # Dirichlet BC on cores
        n = enforce_cores(n, center1, sign1, center2, sign2, r1, r2, core_radius)
        if step % 20 == 0 or step == n_steps - 1:
            E = frank_energy(n, DX)
            energies.append((step, E))
            if verbose:
                print(f"  step {step:4d}  E = {E:.4f}")
    return n, energies


# ---------------------------------------------------------------------------
# FIT
# ---------------------------------------------------------------------------

def fit_coulomb(d, E):
    """Fit E(d) = a + b/d via linear least squares."""
    d, E = np.asarray(d, dtype=float), np.asarray(E, dtype=float)
    X = np.column_stack([np.ones_like(d), 1.0 / d])
    coef, *_ = np.linalg.lstsq(X, E, rcond=None)
    a, b = coef
    E_pred = a + b / d
    ss_res = ((E - E_pred) ** 2).sum()
    ss_tot = ((E - E.mean()) ** 2).sum()
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 1e-12 else 1.0
    return float(a), float(b), float(r2)


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------

def main():
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Grid: {N}³,  domain=[{-DOMAIN/2}, {DOMAIN/2}],  dx={DX:.4f}")
    print(f"Frank K={K_FRANK},  core={CORE_RADIUS},  τ={TAU},  N_relax={N_RELAX}")
    print(f"Separations: {SEPARATIONS}\n")

    # Sanity: single-hedgehog energy scaling
    n = seed_vacuum()
    X, Y, Z = make_grid()
    r_single = np.sqrt(X**2 + Y**2 + Z**2)
    r_safe = np.maximum(r_single, 0.2 * DX)
    n_single = np.stack([X / r_safe, Y / r_safe, Z / r_safe], axis=-1)
    E_single_raw = frank_energy(n_single, DX)
    print(f"Single hedgehog (no vacuum blend) — pre-relax E = {E_single_raw:.3f}")
    print(f"  analytical scale  8π·K·R  with R=domain/2 = {8*np.pi*K_FRANK*DOMAIN/2:.3f}\n")

    # Main run: E(d)
    prerelax_E = []
    measured_E = []
    relax_histories = []
    for d in SEPARATIONS:
        print(f"d = {d:.2f}")
        c1 = np.array([-d / 2.0, 0.0, 0.0])
        c2 = np.array([+d / 2.0, 0.0, 0.0])
        n, r1, r2 = seed_hedgehog_pair(c1, +1, c2, -1)
        E0 = frank_energy(n, DX)
        prerelax_E.append(E0)
        print(f"  pre-relax  E = {E0:.4f}")
        n_relaxed, hist = relax(n, c1, +1, c2, -1)
        E_final = frank_energy(n_relaxed, DX)
        print(f"  post-relax E = {E_final:.4f}")
        measured_E.append(E_final)
        relax_histories.append(hist)

    d_arr = np.array(SEPARATIONS)
    E_arr = np.array(measured_E)
    E_pre = np.array(prerelax_E)

    # Report both fits
    a_pre, b_pre, r2_pre = fit_coulomb(d_arr, E_pre)
    print(f"\nPre-relax fit:   E(d) = {a_pre:.4f} + ({b_pre:.4f})/d    R²={r2_pre:.4f}")
    print(f"  → interaction coefficient b = {b_pre:.2f}  "
          f"({'ATTRACTIVE' if b_pre < 0 else 'REPULSIVE'})")

    a, b, r2 = fit_coulomb(d_arr, E_arr)
    print(f"\nPost-relax fit:  E(d) = {a:.4f} + ({b:.4f})/d    R²={r2:.4f}")
    # Sign convention: E(d) = a + b/d   (with b ∈ ℝ)
    #   b > 0 → E is HIGH near, decays to `a` far → REPULSIVE (like-charges)
    #   b < 0 → E is LOW near,  rises   to `a` far → ATTRACTIVE (opposite charges)
    print(f"  → interaction coefficient b = {b:.2f}  "
          f"({'ATTRACTIVE' if b < 0 else 'REPULSIVE'})")

    # ---------------- PLOTS ----------------
    d_plot = np.linspace(d_arr.min(), d_arr.max(), 200)
    E_fit = a + b / d_plot

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 4.5))

    E_pre_fit = a_pre + b_pre / d_plot
    ax1.plot(d_arr, E_pre, "o", ms=8, label="pre-relax E(d)")
    ax1.plot(d_plot, E_pre_fit, "--",
             label=f"fit: {a_pre:.1f} + ({b_pre:.1f})/d  (R²={r2_pre:.3f})")
    ax1.plot(d_arr, E_arr, "s", ms=6, alpha=0.6, label="post-relax E(d)")
    ax1.set_xlabel("separation d")
    ax1.set_ylabel("Frank elastic energy E")
    ax1.set_title("Exp 2 — Hedgehog pair energy vs distance")
    ax1.legend()
    ax1.grid(alpha=0.3)

    ax2.plot(1.0 / d_arr, E_pre, "o", ms=8, label="measured (pre-relax)")
    ax2.plot(1.0 / d_plot, E_pre_fit, "--", label="linear fit")
    ax2.set_xlabel("1/d")
    ax2.set_ylabel("E")
    ax2.set_title("Exp 2 — E vs 1/d (straight line ⇒ Coulomb)")
    ax2.legend()
    ax2.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(RESULTS_DIR / "energy_vs_distance.png", dpi=110)
    plt.close()

    # Convergence plot
    plt.figure(figsize=(8, 4.5))
    for d, hist in zip(SEPARATIONS, relax_histories):
        steps, Es = zip(*hist)
        plt.plot(steps, Es, "-", label=f"d={d}")
    plt.xlabel("relaxation step")
    plt.ylabel("E")
    plt.title("Exp 2 — Gradient-descent convergence")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(RESULTS_DIR / "relax_convergence.png", dpi=110)
    plt.close()

    print(f"\nPlots saved to: {RESULTS_DIR}")
    print("\n=== SUMMARY ===")
    print(f"{'d':>6}  {'E_pre':>10}  {'E_post':>10}")
    for d, ep, E in zip(d_arr, E_pre, E_arr):
        print(f"{d:>6.2f}  {ep:>10.4f}  {E:>10.4f}")
    print(f"\nPre-relax:   E(d) = {a_pre:.4f} + ({b_pre:.4f})/d    R²={r2_pre:.4f}")
    print(f"Post-relax:  E(d) = {a:.4f} + ({b:.4f})/d    R²={r2:.4f}")


if __name__ == "__main__":
    main()
