"""
Experiment 1: Sine-Gordon 1D Solitons

Numerical Experiment 1 in the Lagrangian Framework Sub-Project.

Hypothesis:
The Sine-Gordon equation produces stable kink solitons that exhibit:
  - Topological stability (kinks cannot dissipate)
  - Pair creation/annihilation (kink + anti-kink collisions)
  - Lorentz contraction (relativistic kinematics from wave equation)

Equation:
    ∂²φ/∂t² - c²·∂²φ/∂x² + (m²c²/ℏ²) · sin(φ) = 0

This comes from the Lagrangian:
    L = ½(∂φ/∂t)² - ½c²(∂φ/∂x)² - (m²c⁴/ℏ²)·(1 - cos(φ))

Spec: ../../3_LAGRANGIAN_FRAMEWORK.md  (Experiment 1 section)
Results: ../3b_lagrangian_experiments.md  (Experiment 1 section)
"""

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# PARAMETERS
# ---------------------------------------------------------------------------

# Grid
N_X = 1024            # number of spatial points
DX = 0.1              # spatial step
DOMAIN = (N_X - 1) * DX   # ensure np.linspace step == DX exactly

# Physics (natural units: c = ℏ = m = 1)
C = 1.0
M = 1.0
HBAR = 1.0

# Time
DT = 0.05             # CFL: c·dt/dx = 0.5 < 1 ✓
N_STEPS = 2000        # total steps → t_final = 100

# Natural kink width
KINK_WIDTH = HBAR / (M * C)   # L = ℏ/(mc) = 1

# Output
RESULTS_DIR = Path(__file__).parent / "exp1_results"


# ---------------------------------------------------------------------------
# INITIAL CONDITION: kink / anti-kink ansatz (Lorentz-boosted)
#   kink:     φ(x, 0) = 4·arctan(exp( γ·(x - x₀) / L))    goes 0 → 2π
#   antikink: φ(x, 0) = -4·arctan(exp( γ·(x - x₀) / L))   goes 0 → -2π
# ---------------------------------------------------------------------------

def kink_profile(x, x0, v, L, c, sign=+1):
    """Static or moving Sine-Gordon kink (sign=+1) or anti-kink (sign=-1)."""
    gamma = 1.0 / np.sqrt(1.0 - (v / c) ** 2)
    return sign * 4.0 * np.arctan(np.exp(gamma * (x - x0) / L))


def kink_antikink_pair(x, xL, xR, v, L, c):
    """Kink at xL moving +v, anti-kink at xR moving -v (approaching)."""
    gammaL = 1.0 / np.sqrt(1.0 - (v / c) ** 2)
    gammaR = gammaL
    # kink (goes 0→2π) + (2π - kink_at_xR) to get a bump 0 → 2π → 0
    kink_left = 4.0 * np.arctan(np.exp(gammaL * (x - xL) / L))
    kink_right = 4.0 * np.arctan(np.exp(gammaR * (x - xR) / L))
    return kink_left - kink_right  # bump: goes 0 → 2π → 0


# ---------------------------------------------------------------------------
# TIME EVOLUTION: leapfrog
#   φ_{n+1} = 2·φ_n - φ_{n-1} + dt²·[c²·∂²φ/∂x² - (m²c⁴/ℏ²)·sin(φ)]
# ---------------------------------------------------------------------------

def laplacian_1d(phi, dx):
    """Second spatial derivative via 3-point stencil. Endpoints held fixed."""
    lap = np.zeros_like(phi)
    lap[1:-1] = (phi[2:] - 2.0 * phi[1:-1] + phi[:-2]) / dx ** 2
    return lap


def step_sine_gordon(phi, phi_prev, dt, dx, c, m, hbar):
    """One leapfrog step."""
    lap = laplacian_1d(phi, dx)
    # EoM: ∂²ₜφ = c²∂²ₓφ - (m²c⁴/ℏ²)·sin(φ)
    mass_coef = (m * c ** 2 / hbar) ** 2
    accel = c ** 2 * lap - mass_coef * np.sin(phi)
    phi_next = 2.0 * phi - phi_prev + dt ** 2 * accel
    return phi_next


# ---------------------------------------------------------------------------
# DIAGNOSTICS
# ---------------------------------------------------------------------------

def kink_energy(phi, phi_dot, dx, c, m, hbar):
    """Total Sine-Gordon Hamiltonian energy.
    H = ∫ [½(∂ₜφ)² + ½c²(∂ₓφ)² + (m²c⁴/ℏ²)·(1-cos(φ))] dx
    """
    kinetic = 0.5 * phi_dot ** 2
    grad = np.gradient(phi, dx)
    gradient = 0.5 * c ** 2 * grad ** 2
    potential = (m * c ** 2 / hbar) ** 2 * (1.0 - np.cos(phi))
    return float(np.sum(kinetic + gradient + potential) * dx)


def kink_center(x, phi, target=np.pi):
    """Locate kink center by linear interpolation of phi = target crossing.
    Assumes a monotonically increasing kink (0 → 2π)."""
    # find the first index where phi crosses target from below to above
    below = phi < target
    # transitions from True to False
    idx = int(np.argmax(np.diff(below.astype(np.int8)) == -1))
    if idx == 0 and not (below[0] and not below[1]):
        # fallback: argmin of |phi - target|
        idx = int(np.argmin(np.abs(phi - target)))
    # linear interp between phi[idx] (below) and phi[idx+1] (above or equal)
    i = max(0, min(idx, len(phi) - 2))
    y0, y1 = phi[i], phi[i + 1]
    if abs(y1 - y0) > 1e-12:
        frac = (target - y0) / (y1 - y0)
    else:
        frac = 0.0
    return float(x[i] + frac * (x[i + 1] - x[i]))


def kink_width(x, phi):
    """Measure kink width by fitting the slope at center.
    For ideal kink: dφ/dx|_{x=x₀} = 2/L  →  L = 2 / slope_at_center.
    """
    grad = np.gradient(phi, x[1] - x[0])
    peak_idx = int(np.argmax(grad))
    peak_slope = grad[peak_idx]
    if peak_slope > 1e-9:
        return 2.0 / peak_slope
    return float("nan")


# ---------------------------------------------------------------------------
# EVOLVE WITH DIAGNOSTICS
# ---------------------------------------------------------------------------

def evolve(phi0, phi_prev0, n_steps, dt, dx, c, m, hbar, track_center=True,
           n_snapshots=11):
    phi = phi0.copy()
    phi_prev = phi_prev0.copy()
    x = np.linspace(-DOMAIN / 2, DOMAIN / 2, len(phi))

    energies = np.empty(n_steps + 1)
    positions = np.empty(n_steps + 1) if track_center else None

    phi_dot0 = (phi - phi_prev) / dt
    energies[0] = kink_energy(phi, phi_dot0, dx, c, m, hbar)
    if track_center:
        positions[0] = kink_center(x, phi)

    snapshot_stride = max(1, n_steps // (n_snapshots - 1))
    snapshots = [(0.0, phi.copy())]

    for n in range(1, n_steps + 1):
        phi_next = step_sine_gordon(phi, phi_prev, dt, dx, c, m, hbar)
        phi_dot = (phi_next - phi_prev) / (2.0 * dt)
        energies[n] = kink_energy(phi, phi_dot, dx, c, m, hbar)
        if track_center:
            positions[n] = kink_center(x, phi)
        if n % snapshot_stride == 0:
            snapshots.append((n * dt, phi.copy()))
        phi_prev = phi
        phi = phi_next

    return {
        "phi_final": phi,
        "x": x,
        "energies": energies,
        "positions": positions,
        "snapshots": snapshots,
    }


# ---------------------------------------------------------------------------
# TESTS
# ---------------------------------------------------------------------------

def test_static_kink():
    """Test 1: static kink holds shape, energy conserved at 8·m·c² = 8."""
    print("\n" + "=" * 72)
    print("TEST 1 — Static kink stability")
    print("=" * 72)
    x = np.linspace(-DOMAIN / 2, DOMAIN / 2, N_X)
    phi0 = kink_profile(x, x0=0.0, v=0.0, L=KINK_WIDTH, c=C)
    phi_prev0 = phi0.copy()  # v=0 → phi(t=-dt) = phi(t=0)
    res = evolve(phi0, phi_prev0, N_STEPS, DT, DX, C, M, HBAR)

    E = res["energies"]
    E_predicted = 8.0 * M * C ** 2
    print(f"Predicted kink rest energy  E = 8·m·c²          = {E_predicted:.6f}")
    print(f"Measured initial energy      E₀                   = {E[0]:.6f}")
    print(f"Measured final energy        E_T                   = {E[-1]:.6f}")
    print(f"Max relative drift           max|ΔE|/E₀            = "
          f"{np.max(np.abs(E - E[0])) / E[0]:.2e}")
    w = kink_width(res["x"], res["phi_final"])
    print(f"Measured final kink width    L                     = {w:.4f}  "
          f"(predicted ℏ/mc = {KINK_WIDTH:.4f})")

    # Plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.2))
    for t, p in res["snapshots"]:
        ax1.plot(res["x"], p, lw=0.9, alpha=0.75, label=f"t={t:.1f}")
    ax1.set_xlim(-20, 20)
    ax1.set_title("Test 1 — Static kink (snapshots)")
    ax1.set_xlabel("x"); ax1.set_ylabel("φ(x,t)")
    ax1.legend(fontsize=7, ncol=2)

    ax2.plot(np.arange(N_STEPS + 1) * DT, E)
    ax2.axhline(E_predicted, color="red", ls="--", label="8·m·c²")
    ax2.set_title("Test 1 — Energy vs time")
    ax2.set_xlabel("t"); ax2.set_ylabel("E")
    ax2.legend()
    plt.tight_layout()
    plt.savefig(RESULTS_DIR / "test1_static_kink.png", dpi=110)
    plt.close()

    return {
        "E_predicted": E_predicted,
        "E_initial": E[0],
        "E_final": E[-1],
        "E_drift_rel": float(np.max(np.abs(E - E[0])) / E[0]),
        "L_measured": w,
    }


def test_moving_kink(v=0.5):
    """Test 2: moving kink at v, measure velocity + Lorentz-contracted width."""
    print("\n" + "=" * 72)
    print(f"TEST 2 — Moving kink (v = {v}·c), Lorentz contraction")
    print("=" * 72)
    x = np.linspace(-DOMAIN / 2, DOMAIN / 2, N_X)
    x0 = -20.0
    gamma = 1.0 / np.sqrt(1.0 - v ** 2)
    L_predicted = KINK_WIDTH / gamma

    phi0 = kink_profile(x, x0=x0, v=v, L=KINK_WIDTH, c=C)
    phi_prev0 = kink_profile(x, x0=x0 - v * DT, v=v, L=KINK_WIDTH, c=C)
    res = evolve(phi0, phi_prev0, N_STEPS, DT, DX, C, M, HBAR)

    # fit velocity from position over time (use central portion to avoid edge effects)
    t = np.arange(N_STEPS + 1) * DT
    pos = res["positions"]
    mask = (t > 5) & (t < N_STEPS * DT * 0.9)
    coef = np.polyfit(t[mask], pos[mask], 1)
    v_measured = coef[0]

    # predicted kink rest energy (moving): E = γ·8·m·c²
    E_predicted = gamma * 8.0 * M * C ** 2
    print(f"Lorentz factor γ                                   = {gamma:.6f}")
    print(f"Predicted velocity (input)                          = {v:.4f}·c")
    print(f"Measured velocity (slope of x(t))                   = {v_measured:.4f}·c")
    print(f"Predicted kink width L/γ                            = {L_predicted:.6f}")
    w = kink_width(res["x"], res["phi_final"])
    print(f"Measured final kink width                           = {w:.6f}")
    print(f"Predicted moving kink energy γ·8·m·c²               = {E_predicted:.6f}")
    print(f"Measured initial energy                              = {res['energies'][0]:.6f}")
    print(f"Max relative drift                                   = "
          f"{np.max(np.abs(res['energies'] - res['energies'][0])) / res['energies'][0]:.2e}")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.2))
    for tt, p in res["snapshots"]:
        ax1.plot(res["x"], p, lw=0.9, alpha=0.8, label=f"t={tt:.1f}")
    ax1.set_title(f"Test 2 — Moving kink (v={v}c)")
    ax1.set_xlabel("x"); ax1.set_ylabel("φ(x,t)")
    ax1.legend(fontsize=7, ncol=2)

    ax2.plot(t, pos, label="measured x(t)")
    ax2.plot(t, x0 + v * t, "--", color="red", label=f"predicted x₀+v·t")
    ax2.set_title("Test 2 — Kink position vs time")
    ax2.set_xlabel("t"); ax2.set_ylabel("x_kink")
    ax2.legend()
    plt.tight_layout()
    plt.savefig(RESULTS_DIR / "test2_moving_kink.png", dpi=110)
    plt.close()

    return {
        "v_input": v,
        "v_measured": float(v_measured),
        "gamma": float(gamma),
        "L_predicted": float(L_predicted),
        "L_measured": float(w),
        "E_predicted": float(E_predicted),
        "E_initial": float(res['energies'][0]),
        "E_drift_rel": float(np.max(np.abs(res['energies'] - res['energies'][0])) / res['energies'][0]),
    }


def test_kink_antikink_collision(v=0.5):
    """Test 3: kink + anti-kink with inward velocities → collision behavior."""
    print("\n" + "=" * 72)
    print(f"TEST 3 — Kink + anti-kink collision (v = ±{v}·c)")
    print("=" * 72)
    x = np.linspace(-DOMAIN / 2, DOMAIN / 2, N_X)
    xL, xR = -15.0, +15.0

    # kink at xL moves +v; anti-kink at xR moves -v (computed via initial phi_prev)
    gamma = 1.0 / np.sqrt(1.0 - v ** 2)

    # Compose the pair profile (0 → 2π → 0 bump)
    phi0 = kink_antikink_pair(x, xL, xR, v, KINK_WIDTH, C)
    # at t=-dt, kink was at xL - v·dt; anti-kink was at xR + v·dt
    phi_prev0 = kink_antikink_pair(x, xL - v * DT, xR + v * DT, v, KINK_WIDTH, C)

    res = evolve(phi0, phi_prev0, N_STEPS, DT, DX, C, M, HBAR, track_center=False)

    E = res["energies"]
    E_pair_predicted = 2.0 * gamma * 8.0 * M * C ** 2
    print(f"Predicted pair energy 2·γ·8·m·c²                    = {E_pair_predicted:.6f}")
    print(f"Measured initial energy                              = {E[0]:.6f}")
    print(f"Measured final energy                                = {E[-1]:.6f}")
    print(f"Max relative drift                                   = "
          f"{np.max(np.abs(E - E[0])) / E[0]:.2e}")
    # estimate peak amplitude over time (bump amplitude tells us about breather vs annihilation)
    amps = np.array([np.max(np.abs(s[1])) for s in res["snapshots"]])
    print(f"Peak |φ| across snapshots:                          "
          f"{np.array2string(amps, precision=3)}")

    fig, axs = plt.subplots(2, 1, figsize=(11, 7))
    for tt, p in res["snapshots"]:
        axs[0].plot(res["x"], p, lw=0.9, alpha=0.7, label=f"t={tt:.1f}")
    axs[0].set_title("Test 3 — Kink + anti-kink evolution")
    axs[0].set_xlabel("x"); axs[0].set_ylabel("φ(x,t)")
    axs[0].legend(fontsize=7, ncol=3)
    axs[0].set_xlim(-40, 40)

    axs[1].plot(np.arange(N_STEPS + 1) * DT, E)
    axs[1].axhline(E_pair_predicted, color="red", ls="--",
                   label="2·γ·8·m·c² (predicted)")
    axs[1].set_title("Test 3 — Total energy vs time")
    axs[1].set_xlabel("t"); axs[1].set_ylabel("E")
    axs[1].legend()
    plt.tight_layout()
    plt.savefig(RESULTS_DIR / "test3_kink_antikink.png", dpi=110)
    plt.close()

    return {
        "E_pair_predicted": float(E_pair_predicted),
        "E_initial": float(E[0]),
        "E_final": float(E[-1]),
        "E_drift_rel": float(np.max(np.abs(E - E[0])) / E[0]),
    }


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------

def main():
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Grid: N={N_X}, dx={DX}, domain=[{-DOMAIN/2:.1f}, {DOMAIN/2:.1f}]")
    print(f"Physics: c={C}, m={M}, ℏ={HBAR}  (natural units)")
    print(f"Time:  dt={DT}, N_STEPS={N_STEPS}, t_final={N_STEPS * DT}")
    print(f"CFL:   c·dt/dx = {C * DT / DX:.3f}  (need < 1)")
    print(f"Natural kink width L = ℏ/(mc) = {KINK_WIDTH}")

    r1 = test_static_kink()
    r2 = test_moving_kink(v=0.5)
    r3 = test_kink_antikink_collision(v=0.5)

    print("\n" + "=" * 72)
    print("SUMMARY")
    print("=" * 72)
    print(f"Plots saved to: {RESULTS_DIR}")
    print(f"\nTest 1 (static):  E/E_pred = {r1['E_initial']/r1['E_predicted']:.4f}  "
          f"drift = {r1['E_drift_rel']:.2e}")
    print(f"Test 2 (moving):  v_meas = {r2['v_measured']:.4f}c (input {r2['v_input']}c),  "
          f"γ = {r2['gamma']:.4f},  L/L_pred = {r2['L_measured']/r2['L_predicted']:.3f}")
    print(f"Test 3 (pair):    E/E_pred = {r3['E_initial']/r3['E_pair_predicted']:.4f}  "
          f"drift = {r3['E_drift_rel']:.2e}")


if __name__ == "__main__":
    main()
