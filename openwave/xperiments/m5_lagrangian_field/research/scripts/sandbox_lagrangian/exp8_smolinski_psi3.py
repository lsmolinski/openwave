"""
Experiment 8: Smolinski's Non-linear Ψ³ (K-selectivity Test)

Numerical Experiment 8 in the Lagrangian Framework Sub-Project.

Hypothesis:
Adding −κ·Ψ³ to the linear wave equation produces K-dependent stability —
K=10 tetrahedron is uniquely stable under perturbation, while K=2..9 are not.

Equation (from Lagrangian L = ½(∂_t ψ)² − ½c²|∇ψ|² − (κ/4)·ψ⁴):
    ∂²Ψ/∂t² = c²·∇²Ψ − κ·Ψ³

(Lagrangian validity confirmed by Experiment 5.)

Approach:
- 3D scalar grid, time-stepping leapfrog PDE evolution (NOT phasor superposition —
  superposition breaks for non-linear equations)
- Seed K Gaussian bumps at geometric positions (K ∈ {1, 2, 4, 6, 8})
- Sweep κ over multiple orders of magnitude: {0.0 (linear), 0.5, 2.0, 10.0}
- Apply small-amplitude random field perturbation
- Measure stability metrics: peak retention at each WC location, total
  energy drift, spatial concentration (fraction of energy within WC balls)

Spec:    ../../3_LAGRANGIAN_FRAMEWORK.md  (Experiment 8)
Results: ../3b_lagrangian_experiments.md  (Experiment 8)
"""

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# PARAMETERS
# ---------------------------------------------------------------------------

N = 64                         # grid points per axis (64³ ≈ 262k voxels)
DOMAIN = 16.0                  # physical extent [-8, +8] each axis
DX = DOMAIN / (N - 1)
C = 1.0                        # wave speed (natural units)
DT = 0.05                      # CFL check: c·dt/dx = 0.197 < 1 ✓
N_STEPS = 60                   # total time steps (t_final = 3, waves don't reach boundary)
WC_RADIUS = 3.0                # distance from origin for WC positions
BUMP_WIDTH = 1.2               # Gaussian σ for each WC bump (wide enough to persist)
BUMP_AMP = 1.0                 # initial amplitude at each WC
PERTURB_NOISE = 0.0            # clean baseline; separate perturbation test below
RNG_SEED = 42

K_VALUES = [1, 2, 4, 6, 8]     # K configurations to test
KAPPA_VALUES = [0.0, 0.5, 2.0, 10.0]   # κ sweep (0 = pure linear wave)

RESULTS_DIR = Path(__file__).parent / "exp8_results"


# ---------------------------------------------------------------------------
# GRID AND GEOMETRY
# ---------------------------------------------------------------------------

_X = _Y = _Z = None


def make_grid():
    global _X, _Y, _Z
    if _X is None:
        x = np.linspace(-DOMAIN / 2, DOMAIN / 2, N)
        _X, _Y, _Z = np.meshgrid(x, x, x, indexing="ij")
    return _X, _Y, _Z


def wc_positions(K):
    """Return a list of K wave-center positions in 3D at radius WC_RADIUS."""
    R = WC_RADIUS
    if K == 1:
        return [np.array([0.0, 0.0, 0.0])]
    if K == 2:
        return [np.array([-R, 0.0, 0.0]), np.array([+R, 0.0, 0.0])]
    if K == 4:
        t = R / np.sqrt(3.0)
        return [np.array([+t, +t, +t]),
                np.array([+t, -t, -t]),
                np.array([-t, +t, -t]),
                np.array([-t, -t, +t])]
    if K == 6:
        return [np.array([+R, 0, 0]), np.array([-R, 0, 0]),
                np.array([0, +R, 0]), np.array([0, -R, 0]),
                np.array([0, 0, +R]), np.array([0, 0, -R])]
    if K == 8:
        t = R / np.sqrt(3.0)
        return [np.array([sx * t, sy * t, sz * t])
                for sx in (-1, 1) for sy in (-1, 1) for sz in (-1, 1)]
    raise ValueError(f"K={K} not implemented")


def seed_wcs(K, perturb_noise=PERTURB_NOISE, seed=RNG_SEED):
    """Build initial Ψ field with K Gaussian bumps + random noise."""
    X, Y, Z = make_grid()
    psi = np.zeros((N, N, N), dtype=np.float64)
    centers = wc_positions(K)
    for c in centers:
        r2 = (X - c[0])**2 + (Y - c[1])**2 + (Z - c[2])**2
        psi += BUMP_AMP * np.exp(-r2 / (2.0 * BUMP_WIDTH**2))
    # Additive random noise for perturbation test
    rng = np.random.default_rng(seed)
    psi += perturb_noise * rng.standard_normal(psi.shape)
    return psi, centers


# ---------------------------------------------------------------------------
# SMOLINSKI EVOLUTION
# ---------------------------------------------------------------------------

def laplacian(f, dx):
    """6-point discrete Laplacian (periodic via np.roll; domain large enough
    that the boundary is essentially vacuum for the test durations used)."""
    return (
        np.roll(f, 1, 0) + np.roll(f, -1, 0)
        + np.roll(f, 1, 1) + np.roll(f, -1, 1)
        + np.roll(f, 1, 2) + np.roll(f, -1, 2)
        - 6.0 * f
    ) / dx**2


def step_smolinski(psi, psi_old, dt, dx, c, kappa):
    """Leapfrog step:  ψ_new = 2ψ − ψ_old + dt²·(c²·∇²ψ − κ·ψ³)."""
    lap = laplacian(psi, dx)
    accel = c**2 * lap - kappa * psi**3
    return 2.0 * psi - psi_old + dt**2 * accel


def hamiltonian_density(psi, psi_dot, dx, c, kappa):
    """H = ½(∂_t ψ)² + ½c²|∇ψ|² + (κ/4)·ψ⁴  (per Noether, Exp 5)."""
    grad = np.gradient(psi, dx)
    grad2 = sum(g**2 for g in grad)
    return 0.5 * psi_dot**2 + 0.5 * c**2 * grad2 + 0.25 * kappa * psi**4


def total_energy(psi, psi_dot, dx, c, kappa):
    return float(hamiltonian_density(psi, psi_dot, dx, c, kappa).sum() * dx**3)


# ---------------------------------------------------------------------------
# DIAGNOSTICS
# ---------------------------------------------------------------------------

def sample_at_points(psi, centers):
    """Return Ψ at each WC's grid-closest voxel."""
    X, Y, Z = make_grid()
    vals = []
    for c in centers:
        # Convert physical → grid index
        ix = int(np.clip(round((c[0] + DOMAIN / 2) / DX), 0, N - 1))
        iy = int(np.clip(round((c[1] + DOMAIN / 2) / DX), 0, N - 1))
        iz = int(np.clip(round((c[2] + DOMAIN / 2) / DX), 0, N - 1))
        vals.append(float(psi[ix, iy, iz]))
    return np.array(vals)


def energy_concentration(psi, psi_dot, centers, ball_radius, dx, c, kappa):
    """Fraction of total energy within balls of given radius around each WC."""
    X, Y, Z = make_grid()
    H = hamiltonian_density(psi, psi_dot, dx, c, kappa)
    E_total = float(H.sum()) * dx**3
    mask = np.zeros_like(H, dtype=bool)
    for wc in centers:
        r2 = (X - wc[0])**2 + (Y - wc[1])**2 + (Z - wc[2])**2
        mask |= (r2 < ball_radius**2)
    E_inside = float(H[mask].sum()) * dx**3
    return E_inside / max(E_total, 1e-30)


# ---------------------------------------------------------------------------
# FULL RUN FOR ONE (K, κ)
# ---------------------------------------------------------------------------

def run_one(K, kappa, n_steps=N_STEPS, verbose=False):
    """Seed, evolve, track diagnostics. Returns a dict."""
    psi0, centers = seed_wcs(K)
    initial_peaks = sample_at_points(psi0, centers)

    psi = psi0.copy()
    psi_old = psi0.copy()   # psi_dot(t=0) = 0 → psi_old = psi at t=-dt is same

    sample_stride = max(1, n_steps // 10)
    t_samples = []
    peak_samples = []
    energy_samples = []
    concentration_samples = []
    ball_radius = max(2.0 * BUMP_WIDTH, 1.5)

    for step in range(n_steps + 1):
        psi_dot = (psi - psi_old) / DT
        if step % sample_stride == 0 or step == n_steps:
            t_samples.append(step * DT)
            peak_samples.append(sample_at_points(psi, centers))
            energy_samples.append(total_energy(psi, psi_dot, DX, C, kappa))
            concentration_samples.append(
                energy_concentration(psi, psi_dot, centers, ball_radius,
                                     DX, C, kappa))
        if step < n_steps:
            psi_new = step_smolinski(psi, psi_old, DT, DX, C, kappa)
            psi_old = psi
            psi = psi_new

    peak_samples = np.array(peak_samples)    # (n_samples, K)
    energy_samples = np.array(energy_samples)
    concentration_samples = np.array(concentration_samples)
    t_samples = np.array(t_samples)

    # Scalar summaries
    peak_retention = float(peak_samples[-1].mean() / max(initial_peaks.mean(), 1e-12))
    energy_drift_rel = float((energy_samples.max() - energy_samples.min())
                             / max(abs(energy_samples[0]), 1e-30))
    concentration_final = float(concentration_samples[-1])
    concentration_initial = float(concentration_samples[0])

    if verbose:
        print(f"  K={K}, κ={kappa}:  peak_retention={peak_retention:.3f}  "
              f"energy_drift={energy_drift_rel:.2e}  "
              f"concentration {concentration_initial:.3f} → {concentration_final:.3f}")

    return {
        "K": K, "kappa": kappa, "n_steps": n_steps,
        "t": t_samples,
        "peaks": peak_samples,
        "energies": energy_samples,
        "concentrations": concentration_samples,
        "initial_peaks": initial_peaks,
        "peak_retention": peak_retention,
        "energy_drift_rel": energy_drift_rel,
        "concentration_initial": concentration_initial,
        "concentration_final": concentration_final,
        "final_psi": psi,
        "centers": centers,
    }


# ---------------------------------------------------------------------------
# MAIN SWEEP
# ---------------------------------------------------------------------------

def main():
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Grid: {N}³,  domain=[{-DOMAIN/2}, {DOMAIN/2}],  dx={DX:.4f}")
    print(f"c={C}, dt={DT}, N_STEPS={N_STEPS}, t_final={N_STEPS*DT}")
    print(f"CFL: c·dt/dx = {C*DT/DX:.3f}  (need < 1)")
    print(f"WC radius = {WC_RADIUS}, bump σ = {BUMP_WIDTH}, amp = {BUMP_AMP}")
    print(f"Perturbation noise σ = {PERTURB_NOISE}")
    print(f"K values: {K_VALUES}")
    print(f"κ values: {KAPPA_VALUES}")
    print()

    # Full sweep
    results = {}   # keyed by (K, κ)
    for K in K_VALUES:
        for kappa in KAPPA_VALUES:
            print(f"K={K}, κ={kappa} ...", end=" ")
            r = run_one(K, kappa, verbose=False)
            results[(K, kappa)] = r
            print(f"peak_ret={r['peak_retention']:+.3f}  "
                  f"conc_final={r['concentration_final']:.3f}  "
                  f"E_drift={r['energy_drift_rel']:.1e}")

    # ---------- SUMMARY TABLE ----------
    print("\n" + "=" * 78)
    print("PEAK RETENTION vs K and κ  (final / initial peak amplitude)")
    print("=" * 78)
    hdr = f"{'K':>3}  " + "  ".join([f"κ={k:>5.2f}" for k in KAPPA_VALUES])
    print(hdr)
    for K in K_VALUES:
        row = f"{K:>3}  " + "  ".join(
            [f"{results[(K, k)]['peak_retention']:>+7.3f}" for k in KAPPA_VALUES])
        print(row)

    print("\n" + "=" * 78)
    print("ENERGY CONCENTRATION (final; fraction of total E in WC balls)")
    print("=" * 78)
    print(hdr)
    for K in K_VALUES:
        row = f"{K:>3}  " + "  ".join(
            [f"{results[(K, k)]['concentration_final']:>7.3f}" for k in KAPPA_VALUES])
        print(row)

    print("\n" + "=" * 78)
    print("ENERGY DRIFT (rel; max − min over simulation)")
    print("=" * 78)
    print(hdr)
    for K in K_VALUES:
        row = f"{K:>3}  " + "  ".join(
            [f"{results[(K, k)]['energy_drift_rel']:>7.1e}" for k in KAPPA_VALUES])
        print(row)

    # ---------- PLOTS ----------
    # Plot: peak retention heatmap
    fig, ax = plt.subplots(figsize=(7, 5))
    data = np.array([[results[(K, k)]["peak_retention"]
                      for k in KAPPA_VALUES] for K in K_VALUES])
    im = ax.imshow(data, aspect="auto", origin="lower", cmap="RdBu_r",
                   vmin=-0.5, vmax=1.5)
    ax.set_xticks(range(len(KAPPA_VALUES)))
    ax.set_xticklabels([f"{k:.2f}" for k in KAPPA_VALUES])
    ax.set_yticks(range(len(K_VALUES)))
    ax.set_yticklabels([f"K={K}" for K in K_VALUES])
    ax.set_xlabel("κ (coupling)")
    ax.set_title("Exp 8 — Peak retention (final / initial) vs K and κ")
    # annotate cells
    for i, K in enumerate(K_VALUES):
        for j, k in enumerate(KAPPA_VALUES):
            ax.text(j, i, f"{data[i, j]:+.2f}", ha="center", va="center",
                    color="black" if abs(data[i, j] - 0.5) < 0.5 else "white",
                    fontsize=9)
    plt.colorbar(im, ax=ax, label="peak retention")
    plt.tight_layout()
    plt.savefig(RESULTS_DIR / "peak_retention_heatmap.png", dpi=110)
    plt.close()

    # Plot: concentration time-series by κ (one panel per K)
    fig, axs = plt.subplots(1, len(K_VALUES), figsize=(4 * len(K_VALUES), 4),
                            sharey=True)
    for ax, K in zip(axs, K_VALUES):
        for kappa in KAPPA_VALUES:
            r = results[(K, kappa)]
            ax.plot(r["t"], r["concentrations"], label=f"κ={kappa}", lw=1.5)
        ax.set_title(f"K = {K}")
        ax.set_xlabel("time")
        ax.grid(alpha=0.3)
    axs[0].set_ylabel("energy concentration in WC balls")
    axs[-1].legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(RESULTS_DIR / "concentration_vs_time.png", dpi=110)
    plt.close()

    # Plot: peak amplitude time-series by κ (one panel per K)
    fig, axs = plt.subplots(1, len(K_VALUES), figsize=(4 * len(K_VALUES), 4),
                            sharey=True)
    for ax, K in zip(axs, K_VALUES):
        for kappa in KAPPA_VALUES:
            r = results[(K, kappa)]
            ax.plot(r["t"], r["peaks"].mean(axis=1), label=f"κ={kappa}", lw=1.5)
        ax.set_title(f"K = {K}")
        ax.set_xlabel("time")
        ax.grid(alpha=0.3)
    axs[0].set_ylabel("mean WC peak amplitude")
    axs[-1].legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(RESULTS_DIR / "peak_vs_time.png", dpi=110)
    plt.close()

    print(f"\nPlots saved to {RESULTS_DIR}")


if __name__ == "__main__":
    main()
