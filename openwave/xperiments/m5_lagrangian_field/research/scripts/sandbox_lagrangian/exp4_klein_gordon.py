"""
Experiment 4: Klein-Gordon from Twist Dynamics

Numerical Experiment 4 in the Lagrangian Framework Sub-Project.

Hypothesis:
Perturbations of a uniaxial director vacuum produce massive Klein-Gordon
dispersion: ω² = c²·k² + m². More generally, any scalar field with the
Lagrangian
    L = ½(∂_t ψ)² − ½c²·|∇ψ|² − ½m²·ψ²
obeys the Klein-Gordon PDE
    ∂²ψ/∂t² − c²·∇²ψ + m²·ψ = 0
whose plane-wave solutions satisfy the relativistic dispersion relation
    ω² = c²·k² + m².

Physics connection:
In the Landau-de Gennes uniaxial limit, the twist degree of freedom of
the director field reduces (to leading order) to a scalar with rotational
inertia χ and twist elastic constant K. With a mass term from the LdG
potential V(M), the EoM for small twist perturbations is Klein-Gordon.
Testing the scalar Klein-Gordon dispersion here validates the core
mechanism; a full 3D director simulation (Exp 6) will test the biaxial
generalization later.

Method:
1. 1D periodic grid.
2. For each wavenumber k in a sweep, seed ψ(x, 0) = cos(k·x), ψ_dot = 0.
3. Evolve with leapfrog: ψ_new = 2ψ − ψ_old + dt²·(c²·∇²ψ − m²·ψ).
4. Sample ψ at a fixed point over time → extract dominant angular
   frequency ω via FFT.
5. Fit ω² = c²·k² + m² across the k sweep.

Spec:    ../../3_LAGRANGIAN_FRAMEWORK.md  (Experiment 4)
Results: ../3b_lagrangian_experiments.md (Experiment 4)
"""

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# PARAMETERS
# ---------------------------------------------------------------------------

N = 1024                       # spatial grid points (periodic)
DOMAIN = 40.0                  # spatial extent
DX = DOMAIN / N                # periodic: N cells, width DX
C = 1.0                        # wave speed
N_STEPS = 3000                 # temporal samples per run
DT = 0.02                      # CFL: c·dt/dx ≈ 0.512 < 1 ✓
# Mass parameter for the primary test (m² > 0 ⇒ massive dispersion)
M = 0.7
# Wavenumber sweep (rad/length). Choose k·L to be integer multiples of 2π/DOMAIN
# so each plane wave fits a whole number of periods into the periodic box.
K_MODES = [1, 2, 3, 5, 7, 10, 14, 20, 28]   # mode numbers n (k = 2π·n / DOMAIN)
# Control: re-run at m = 0 on a couple of k values to confirm the
# light-cone dispersion ω = c·k (zero-mass limit).
CONTROL_MASSES = [0.0, 0.7]

RESULTS_DIR = Path(__file__).parent / "exp4_results"


# ---------------------------------------------------------------------------
# PDE EVOLUTION
# ---------------------------------------------------------------------------

def laplacian_1d(psi, dx):
    """1D periodic Laplacian via 3-point stencil (np.roll wraps boundaries)."""
    return (np.roll(psi, 1) - 2.0 * psi + np.roll(psi, -1)) / dx**2


def step_klein_gordon(psi, psi_old, dt, dx, c, m):
    """Leapfrog step: ψ_new = 2·ψ − ψ_old + dt²·(c²·∇²ψ − m²·ψ)."""
    lap = laplacian_1d(psi, dx)
    accel = c**2 * lap - m**2 * psi
    return 2.0 * psi - psi_old + dt**2 * accel


def evolve(k_wave, mass, n_steps=N_STEPS, dt=DT, dx=DX, c=C):
    """Seed ψ = cos(k·x), ψ_dot = 0 (periodic); evolve, sample at x=0."""
    x = np.linspace(-DOMAIN / 2, DOMAIN / 2, N, endpoint=False)
    psi_curr = np.cos(k_wave * x)
    psi_old = psi_curr.copy()     # zero initial velocity
    sample_idx = N // 2           # x=0 in the centred grid

    times = np.empty(n_steps + 1, dtype=np.float64)
    psi_samples = np.empty(n_steps + 1, dtype=np.float64)
    times[0] = 0.0
    psi_samples[0] = psi_curr[sample_idx]

    for step in range(1, n_steps + 1):
        psi_new = step_klein_gordon(psi_curr, psi_old, dt, dx, c, mass)
        psi_old = psi_curr
        psi_curr = psi_new
        times[step] = step * dt
        psi_samples[step] = psi_curr[sample_idx]

    return times, psi_samples


# ---------------------------------------------------------------------------
# FFT-BASED FREQUENCY EXTRACTION
# ---------------------------------------------------------------------------

def measure_omega(times, psi_samples):
    """Extract the dominant positive angular frequency from a time series."""
    n = len(times)
    dt = times[1] - times[0]
    # Remove DC, apply Hann window to reduce spectral leakage
    sig = psi_samples - psi_samples.mean()
    window = np.hanning(n)
    fft = np.fft.rfft(sig * window)
    freqs = np.fft.rfftfreq(n, dt)  # cycles per unit time (Hz)
    power = np.abs(fft)**2
    # Skip DC bin
    peak = int(np.argmax(power[1:]) + 1)
    # Parabolic refinement for sub-bin accuracy
    if 0 < peak < len(power) - 1:
        y0, y1, y2 = power[peak - 1], power[peak], power[peak + 1]
        denom = y0 - 2.0 * y1 + y2
        if abs(denom) > 1e-12:
            shift = 0.5 * (y0 - y2) / denom
        else:
            shift = 0.0
    else:
        shift = 0.0
    f_peak = (peak + shift) * (freqs[1] - freqs[0])
    omega = 2.0 * np.pi * f_peak
    return float(omega)


# ---------------------------------------------------------------------------
# FITS
# ---------------------------------------------------------------------------

def fit_dispersion(k_arr, omega_arr):
    """Fit ω² = a·k² + b  (expect a = c², b = m²)."""
    k_arr = np.asarray(k_arr)
    omega_arr = np.asarray(omega_arr)
    X = np.column_stack([k_arr**2, np.ones_like(k_arr)])
    coef, *_ = np.linalg.lstsq(X, omega_arr**2, rcond=None)
    a, b = float(coef[0]), float(coef[1])
    pred = a * k_arr**2 + b
    ss_res = float(((omega_arr**2 - pred) ** 2).sum())
    ss_tot = float(((omega_arr**2 - omega_arr.mean()**2) ** 2).sum() + 1e-30)
    r2 = 1.0 - ss_res / ss_tot
    return a, b, r2


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------

def main():
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Grid:  N={N}, dx={DX:.4f}, domain=[{-DOMAIN/2}, {DOMAIN/2}]")
    print(f"Time:  dt={DT}, N_STEPS={N_STEPS}, t_final={N_STEPS*DT}")
    print(f"CFL:   c·dt/dx = {C*DT/DX:.3f}  (need < 1)")
    print(f"Mass sweep: {CONTROL_MASSES}")
    print(f"Mode sweep: n={K_MODES} → k = 2π·n/DOMAIN ∈ "
          f"{[round(2*np.pi*n/DOMAIN, 3) for n in K_MODES]}")
    print()

    all_results = {}   # keyed by mass
    for mass in CONTROL_MASSES:
        print(f"\n{'='*78}\nTest run — mass m = {mass}\n{'='*78}")
        k_values = []
        omega_measured = []
        omega_expected = []
        for mode_n in K_MODES:
            k = 2.0 * np.pi * mode_n / DOMAIN
            times, psi = evolve(k, mass)
            omega = measure_omega(times, psi)
            w_predicted = np.sqrt(C**2 * k**2 + mass**2)
            print(f"  n={mode_n:>3}  k={k:>6.3f}  "
                  f"ω_measured={omega:>7.4f}  "
                  f"ω_predicted={w_predicted:>7.4f}  "
                  f"ratio={omega / w_predicted:.4f}")
            k_values.append(k)
            omega_measured.append(omega)
            omega_expected.append(w_predicted)
        k_arr = np.array(k_values)
        om_arr = np.array(omega_measured)
        expected_arr = np.array(omega_expected)

        a, b, r2 = fit_dispersion(k_arr, om_arr)
        print(f"\n  Fit ω² = a·k² + b:   a = {a:.4f}  (expect c² = {C**2:.4f})")
        print(f"                       b = {b:.4f}  (expect m² = {mass**2:.4f})")
        print(f"                       R² = {r2:.6f}")

        all_results[mass] = {
            "k": k_arr, "omega_measured": om_arr, "omega_expected": expected_arr,
            "a_fit": a, "b_fit": b, "r2": r2,
        }

    # --------------------------- PLOTS ---------------------------
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

    for mass, r in all_results.items():
        ax1.plot(r["k"], r["omega_measured"], "o", ms=8, label=f"m={mass} measured")
        k_cont = np.linspace(r["k"].min() * 0.95, r["k"].max() * 1.05, 200)
        w_curve = np.sqrt(C**2 * k_cont**2 + mass**2)
        ax1.plot(k_cont, w_curve, "--", label=f"m={mass} predicted")
    ax1.set_xlabel("k")
    ax1.set_ylabel("ω")
    ax1.set_title("Exp 4 — Dispersion ω vs k")
    ax1.legend()
    ax1.grid(alpha=0.3)

    for mass, r in all_results.items():
        ax2.plot(r["k"]**2, r["omega_measured"]**2, "o", ms=8,
                 label=f"m={mass} measured")
        k2_line = np.linspace(0, r["k"].max()**2 * 1.05, 50)
        ax2.plot(k2_line, C**2 * k2_line + mass**2, "--",
                 label=f"m={mass}: ω²=c²k²+m²")
    ax2.set_xlabel("k²")
    ax2.set_ylabel("ω²")
    ax2.set_title("Exp 4 — Dispersion ω² vs k² (linear fit check)")
    ax2.legend()
    ax2.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(RESULTS_DIR / "dispersion.png", dpi=110)
    plt.close()

    # Summary
    print("\n" + "=" * 78)
    print("SUMMARY")
    print("=" * 78)
    for mass, r in all_results.items():
        print(f"  m = {mass}:   a={r['a_fit']:.4f} (expect {C**2:.4f}),  "
              f"b={r['b_fit']:.4f} (expect {mass**2:.4f}),  R²={r['r2']:.6f}")
    print(f"\nPlots saved to {RESULTS_DIR}")


if __name__ == "__main__":
    main()
