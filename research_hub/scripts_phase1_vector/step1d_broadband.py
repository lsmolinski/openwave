#!/usr/bin/env python3
"""
Phase 1d: Broadband WC Out-Wave from Yee & Hauger Shells
=========================================================
Hypothesis: if each shell of the particle emits at a DIFFERENT wavelength,
the far-field out-wave is broadband (multiple k values). The interaction:

  E_int ∝ Σ_k |Z₁(k)| · |Z₂(k)| · cos(k·Δr + Δφ)

At each k, cos(k·Δr) oscillates at a different rate. When summed, the
oscillations partially cancel (destructive interference of the cosines
themselves). The net might converge to a definite direction.

Key: K=1 has 1 shell → monochromatic → sinc persists → neutral
     K=10 has 10 shells at different λ → broadband → sinc smooths → Coulomb?

Each shell n has wavelength λ_n = 2(K-n)λ₀ and contributes an out-wave
at wavenumber k_n = 2π/λ_n. The WC out-wave is the sum of these K components.

Usage: python step1d_broadband.py
"""

import numpy as np
import time

# ============================================================
# CONSTANTS
# ============================================================
A0   = 0.92
LAM0 = 28.5
F0   = 0.0105
RHO0 = 38.6
C    = 0.3
K0   = 2 * np.pi / LAM0

# ============================================================
# CONFIGS
# ============================================================
K_VALUES = [1, 2, 3, 5, 10]
N_SEPS = 50

PHASE_CONFIGS = [
    (0.0,   0.0,   "same"),
    (0.0,   np.pi, "opp"),
]


# ============================================================
# BROADBAND WC OUT-WAVE
# ============================================================

def shell_wavelengths(K):
    """
    Yee & Hauger shell wavelengths: λ_n = 2(K-n)λ₀ for n=1..K
    Returns array of K wavelengths (excluding n=K which gives λ=0).
    """
    n = np.arange(1, K + 1)
    lam_n = 2 * (K - n) * LAM0
    # Exclude zero-wavelength shells and very small ones
    valid = lam_n > LAM0 * 0.5
    return lam_n[valid]


def broadband_outwave_1d(x, wc_x, phase, K, amplitude):
    """
    Broadband WC out-wave: sum of spherical waves at each shell wavelength.

    Each shell n contributes: A_n · sin(k_n·r - ωt - φ) / (k_n·r)
    with amplitude A_n proportional to shell width (more energy in wider shells).

    In phasor form: Z_n = -i·A_n·exp(+i(k_n·r - φ)) / (k_n·r)
    """
    r = np.abs(x - wc_x)
    r_safe = np.maximum(r, LAM0 * 0.01)

    lambdas = shell_wavelengths(K)
    if len(lambdas) == 0:
        # K=1 edge case: only one shell with λ=0 → use fundamental
        lambdas = np.array([LAM0])

    # Amplitude per shell: proportional to shell width (wider = more energy)
    # Normalize so total energy matches monochromatic case
    weights = lambdas / lambdas.sum()
    A_shells = amplitude * np.sqrt(weights)  # sqrt because E ∝ A²

    Z_total = np.zeros_like(x, dtype=complex)

    for lam_n, A_n in zip(lambdas, A_shells):
        k_n = 2 * np.pi / lam_n
        kr_n = k_n * r_safe
        # Out-wave phasor at this wavelength
        Z_n = -1j * A_n * np.exp(1j * (kr_n - phase)) / kr_n
        Z_total += Z_n

    return Z_total


def monochromatic_outwave_1d(x, wc_x, phase, amplitude):
    """Single-frequency out-wave for comparison."""
    r = np.abs(x - wc_x)
    r_safe = np.maximum(r, LAM0 * 0.01)
    kr = K0 * r_safe
    return -1j * amplitude * np.exp(1j * (kr - phase)) / kr


# ============================================================
# FORCE COMPUTATION
# ============================================================

def interaction_force_1d(x, Z1, Z2):
    """Compute interaction force from two phasor fields."""
    dx = x[1] - x[0]
    E_12 = RHO0 * dx * F0**2 * np.abs(Z1 + Z2)**2 / 2
    E_1 = RHO0 * dx * F0**2 * np.abs(Z1)**2 / 2
    E_2 = RHO0 * dx * F0**2 * np.abs(Z2)**2 / 2
    E_int = E_12 - E_1 - E_2
    return -np.gradient(E_int, x)


# ============================================================
# MAIN
# ============================================================

def main():
    print("=" * 70)
    print("  PHASE 1d :  Broadband WC from Yee & Hauger Shells")
    print("=" * 70)
    print(f"  K values: {K_VALUES}")
    print()

    # Show shell wavelengths for each K
    for K in K_VALUES:
        lams = shell_wavelengths(K)
        ks = 2 * np.pi / lams if len(lams) > 0 else np.array([K0])
        print(f"  K={K:2d}: {len(lams)} shells, λ/λ₀ = [{', '.join(f'{l/LAM0:.1f}' for l in lams)}]")

    for phase1, phase2, charge_label in PHASE_CONFIGS:
        print(f"\n{'#' * 70}")
        print(f"  {charge_label.upper()} CHARGE (φ₁={phase1:.1f}, φ₂={phase2:.1f})")
        print(f"  Coulomb expects: {'REP' if charge_label == 'same' else 'ATT'}")
        print(f"{'#' * 70}")

        for K in K_VALUES:
            r_particle = K**2 * LAM0
            sep_min = max(2.0 * r_particle, 5 * LAM0)
            sep_max = max(4.0 * r_particle, 15 * LAM0)
            seps = np.linspace(sep_min, sep_max, N_SEPS)

            # Grid
            margin = 20 * LAM0
            n_pts = 50000
            x = np.linspace(-margin, seps.max() + margin, n_pts)

            broad_dirs = []
            mono_dirs = []

            for sep in seps:
                # Broadband WC out-waves
                Z1_broad = broadband_outwave_1d(x, 0.0, phase1, K, A0)
                Z2_broad = broadband_outwave_1d(x, sep, phase2, K, A0)
                F_broad = interaction_force_1d(x, Z1_broad, Z2_broad)

                # Monochromatic reference
                Z1_mono = monochromatic_outwave_1d(x, 0.0, phase1, A0)
                Z2_mono = monochromatic_outwave_1d(x, sep, phase2, A0)
                F_mono = interaction_force_1d(x, Z1_mono, Z2_mono)

                # Sample at WC1
                idx = np.argmin(np.abs(x))
                broad_dirs.append("ATT" if F_broad[idx] > 0 else "REP")
                mono_dirs.append("ATT" if F_mono[idx] > 0 else "REP")

            # Count
            ba = broad_dirs.count("ATT")
            br = broad_dirs.count("REP")
            ma = mono_dirs.count("ATT")
            mr = mono_dirs.count("REP")
            bc = "CONSISTENT" if (ba == N_SEPS or br == N_SEPS) else "MIXED"
            mc = "CONSISTENT" if (ma == N_SEPS or mr == N_SEPS) else "MIXED"
            bd = "ATT" if ba > br else "REP"
            md = "ATT" if ma > mr else "REP"

            coulomb_dir = "REP" if charge_label == "same" else "ATT"
            b_match = "COULOMB" if (bc == "CONSISTENT" and bd == coulomb_dir) else ""
            if bc == "CONSISTENT" and bd != coulomb_dir:
                b_match = "INVERTED"

            print(f"\n  K={K:2d}  seps={seps[0]/LAM0:.0f}-{seps[-1]/LAM0:.0f}λ")
            print(f"    Mono:  {ma:2d}/{N_SEPS}A {mr:2d}/{N_SEPS}R  {mc:>10s}  dom={md}")
            print(f"    Broad: {ba:2d}/{N_SEPS}A {br:2d}/{N_SEPS}R  {bc:>10s}  dom={bd}  {b_match}")

    print(f"\n{'=' * 70}")
    print(f"  Does broadband break the sinc at higher K?")
    print(f"  K=1 mono=broad (only 1 shell)")
    print(f"  K=10: if broad=CONSISTENT and mono=MIXED → broadband IS the mechanism")
    print(f"{'=' * 70}")


if __name__ == "__main__":
    main()
