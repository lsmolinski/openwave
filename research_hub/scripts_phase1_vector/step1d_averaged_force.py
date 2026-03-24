#!/usr/bin/env python3
"""
Phase 1d: Statistical Averaging of Force Over Particle Radius
=============================================================
Hypothesis: the Coulomb force is the sinc interaction force AVERAGED over
the particle's standing wave volume (radius K²λ).

At K=1 (neutrino): averaging window = 1λ → covers ~2 sinc half-cycles → no smoothing → neutral
At K=10 (electron): averaging window = 100λ → covers ~200 sinc half-cycles → smooths out → Coulomb

This directly explains why neutrinos are neutral (K=1, no averaging) and
electrons are charged (K=10, averaging over many oscillations).

Method:
  For each separation d between two WCs, compute the interaction force
  F_int(x) = -dE_int/dx along the axis, then average F_int over a window
  of size K²λ centered on WC1. Compare raw vs averaged force.

Usage: python step1d_averaged_force.py
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
K_VALUES = [1, 2, 3, 5, 10]  # Test multiple K to see transition

# Fine separation sweep — many points to see the sinc and averaging clearly
N_SEPS = 50

PHASE_CONFIGS = [
    (0.0,   0.0,   "same"),
    (0.0,   np.pi, "opp"),
]


# ============================================================
# 1D COMPUTATION (high resolution along axis)
# ============================================================

def compute_interaction_force_1d(sep, phase1, phase2, n_points=50000):
    """
    Compute the 1D interaction force between two WCs along the x-axis.
    Returns: x array, F_interaction(x) array
    """
    # Grid: extends well beyond both WCs
    margin = 20 * LAM0
    x = np.linspace(-margin, sep + margin, n_points)
    dx = x[1] - x[0]

    # WC positions
    wc1_x = 0.0
    wc2_x = sep

    # Distance from each WC
    r1 = np.abs(x - wc1_x)
    r2 = np.abs(x - wc2_x)

    # Phasors (out-wave only: -i·exp(+i(kr-φ))/kr)
    kr1 = K0 * np.maximum(r1, LAM0 * 0.01)
    kr2 = K0 * np.maximum(r2, LAM0 * 0.01)

    Z1 = -1j * A0 * np.exp(1j * (kr1 - phase1)) / kr1
    Z2 = -1j * A0 * np.exp(1j * (kr2 - phase2)) / kr2

    # Interaction energy: E(1+2) - E(1) - E(2)
    E_combined = RHO0 * dx * F0**2 * np.abs(Z1 + Z2)**2 / 2
    E_1 = RHO0 * dx * F0**2 * np.abs(Z1)**2 / 2
    E_2 = RHO0 * dx * F0**2 * np.abs(Z2)**2 / 2
    E_int = E_combined - E_1 - E_2

    # Force: F = -dE/dx
    F_int = -np.gradient(E_int, x)

    return x, F_int


def average_force_over_window(x, F, center, window_radius):
    """Average the force over a window centered on `center`."""
    mask = np.abs(x - center) <= window_radius
    if mask.sum() == 0:
        return 0.0
    return F[mask].mean()


# ============================================================
# MAIN
# ============================================================

def main():
    print("=" * 70)
    print("  PHASE 1d :  Statistical Averaging of Force Over Particle Radius")
    print("=" * 70)
    print(f"  K values: {K_VALUES}")
    print(f"  Hypothesis: K=1 → sinc (neutral), K=10 → smooth (Coulomb)")
    print()

    for phase1, phase2, charge_label in PHASE_CONFIGS:
        print(f"\n{'#' * 70}")
        print(f"  CHARGE CONFIG: {charge_label} (φ₁={phase1:.1f}, φ₂={phase2:.1f})")
        print(f"  Coulomb expects: {'REPULSION (same)' if charge_label == 'same' else 'ATTRACTION (opposite)'}")
        print(f"{'#' * 70}")

        # Results by K
        for K in K_VALUES:
            r_particle = K**2 * LAM0
            window = r_particle  # Average over particle radius

            # Separation range: from 2× particle radius to 4× (far-field)
            sep_min = 2.5 * r_particle
            sep_max = 4.0 * r_particle
            seps = np.linspace(sep_min, sep_max, N_SEPS)

            raw_forces = []
            avg_forces = []

            for sep in seps:
                x, F_int = compute_interaction_force_1d(sep, phase1, phase2)
                # Raw force at WC1 center
                idx_wc1 = np.argmin(np.abs(x))
                f_raw = F_int[idx_wc1]
                # Averaged force over particle radius
                f_avg = average_force_over_window(x, F_int, 0.0, window)

                raw_forces.append(f_raw)
                avg_forces.append(f_avg)

            raw_forces = np.array(raw_forces)
            avg_forces = np.array(avg_forces)
            seps_lam = seps / LAM0

            # Count direction flips
            raw_att = np.sum(raw_forces > 0)
            raw_rep = np.sum(raw_forces < 0)
            avg_att = np.sum(avg_forces > 0)
            avg_rep = np.sum(avg_forces < 0)

            raw_cons = "CONSISTENT" if (raw_att == N_SEPS or raw_rep == N_SEPS) else "MIXED"
            avg_cons = "CONSISTENT" if (avg_att == N_SEPS or avg_rep == N_SEPS) else "MIXED"

            raw_dom = "ATT" if raw_att > raw_rep else "REP"
            avg_dom = "ATT" if avg_att > avg_rep else "REP"

            # Check if averaged force matches Coulomb expectation
            coulomb_dir = "REP" if charge_label == "same" else "ATT"
            avg_correct = "COULOMB" if (avg_cons == "CONSISTENT" and avg_dom == coulomb_dir) else ""
            if avg_cons == "CONSISTENT" and avg_dom != coulomb_dir:
                avg_correct = "INVERTED"

            print(f"\n  K={K:2d}  window={K**2:4d}λ  seps={seps_lam[0]:.0f}-{seps_lam[-1]:.0f}λ")
            print(f"    Raw:  {raw_att:2d}/{N_SEPS}ATT  {raw_rep:2d}/{N_SEPS}REP  {raw_cons:>10s}  dom={raw_dom}")
            print(f"    Avg:  {avg_att:2d}/{N_SEPS}ATT  {avg_rep:2d}/{N_SEPS}REP  {avg_cons:>10s}  dom={avg_dom}  {avg_correct}")

            # Print a few sample values
            for i in [0, N_SEPS//4, N_SEPS//2, 3*N_SEPS//4, N_SEPS-1]:
                d_raw = "ATT" if raw_forces[i] > 0 else "REP"
                d_avg = "ATT" if avg_forces[i] > 0 else "REP"
                print(f"      sep={seps_lam[i]:7.1f}λ  raw={raw_forces[i]:+10.3e}({d_raw})  "
                      f"avg={avg_forces[i]:+10.3e}({d_avg})")

    # ---- Summary ----
    print(f"\n{'=' * 70}")
    print(f"  SUMMARY: Does averaging over K²λ produce Coulomb?")
    print(f"{'=' * 70}")
    print(f"  Coulomb = CONSISTENT direction matching same→REP, opp→ATT")
    print(f"  If K=1 MIXED and K=10 CONSISTENT → K determines charge emergence")
    print(f"{'=' * 70}")


if __name__ == "__main__":
    main()
