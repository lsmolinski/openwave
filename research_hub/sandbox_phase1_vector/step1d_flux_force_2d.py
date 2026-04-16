#!/usr/bin/env python3
"""
Phase 1d: 2D Flux-Based Force (Radiation Pressure)
===================================================
LaFreniere's animation shows 2D wave interference with clear amplitude
asymmetry: lower between opposite-phase WCs (destructive), higher outside.
The 2D off-axis paths smooth the sinc — different angles have different
effective separations (d·cos(θ)), so the force flips cancel when integrated.

Compute radiation pressure force in 2D:
  1. Full time-domain wave field ψ(x,y,t) on 2D grid
  2. Energy flux S = -c²·ψ·∇ψ (2D vector)
  3. Time-average over one full period
  4. Force on WC1 = net flux through a circle surrounding WC1

This is what LaFreniere's animation actually computes — 2D, not 1D.

Usage: python step1d_flux_force_2d.py
"""

import numpy as np
import time

# ============================================================
# CONSTANTS
# ============================================================
A0    = 0.92
LAM0  = 28.5
F0    = 0.0105
RHO0  = 38.6
C     = 0.3
K0    = 2 * np.pi / LAM0
OMEGA = 2 * np.pi * F0

# ============================================================
# GRID — large enough for K=1, good resolution
# ============================================================
N_GRID   = 512
GRID_LAM = 16.0
L        = GRID_LAM * LAM0
DX       = L / N_GRID

# Weight function (same as M3)
TRANSITION_LAM = 1.25
WEIGHT_POWER   = 8

# Time averaging
N_TIME_SAMPLES = 32

# Separations
SEPS_LAM = np.arange(2.0, 7.5, 0.25)  # Fine sweep every quarter-λ

PHASE_CONFIGS = [
    (0.0,   0.0,   "same"),
    (0.0,   np.pi, "opp"),
]


# ============================================================
# 2D GRID
# ============================================================

def build_grid_2d(n, dx):
    ax = (np.arange(n) - n / 2.0 + 0.5) * dx
    gx, gy = np.meshgrid(ax, ax, indexing='ij')
    return gx, gy


# ============================================================
# TIME-DOMAIN WAVE FIELD (2D)
# ============================================================

def wave_field_2d(gx, gy, wc_x, wc_y, phase, amplitude, t):
    """
    Weighted partial standing wave in 2D at time t.
    ψ = A · [w·sin(kr + ωt + φ) + sin(kr - ωt - φ)] / kr
    Returns ψ(x,y,t).
    """
    rx = gx - wc_x
    ry = gy - wc_y
    r = np.sqrt(rx**2 + ry**2)
    r_safe = np.maximum(r, LAM0 * 0.01)
    kr = K0 * r_safe

    w = 1.0 / (1.0 + (r_safe / (TRANSITION_LAM * LAM0)) ** WEIGHT_POWER)

    in_wave = w * np.sin(kr + OMEGA * t + phase)
    out_wave = np.sin(kr - OMEGA * t - phase)

    psi = amplitude * (in_wave + out_wave) / kr
    return psi


# ============================================================
# FLUX AND RADIATION PRESSURE (2D)
# ============================================================

def compute_flux_2d(psi, dx):
    """
    Energy flux vector: S = -c²·ψ·∇ψ
    Returns Sx, Sy on 2D grid.
    """
    dpsi_dx = np.zeros_like(psi)
    dpsi_dy = np.zeros_like(psi)
    dpsi_dx[1:-1, :] = (psi[2:, :] - psi[:-2, :]) / (2 * dx)
    dpsi_dy[:, 1:-1] = (psi[:, 2:] - psi[:, :-2]) / (2 * dx)

    Sx = -C**2 * psi * dpsi_dx
    Sy = -C**2 * psi * dpsi_dy
    return Sx, Sy


def radiation_force_on_wc(Sx, Sy, gx, gy, wc_x, wc_y, dx, sample_radius):
    """
    Net radiation pressure force on WC by integrating flux through a circle.

    Integrate the OUTWARD flux through a circle of radius sample_radius
    centered on the WC. Net outward flux = energy leaving → no force.
    Asymmetric flux = net force.

    Force = -∮ S · n̂ dℓ (inward momentum = push)
    In discrete: sum flux·normal over ring of grid points.
    """
    rx = gx - wc_x
    ry = gy - wc_y
    r = np.sqrt(rx**2 + ry**2)

    # Select ring of grid points at ~sample_radius
    ring_width = dx * 1.5
    ring = (r >= sample_radius - ring_width) & (r <= sample_radius + ring_width)

    if ring.sum() == 0:
        return 0.0, 0.0

    # Unit outward normal at each ring point
    r_ring = r[ring]
    nx = rx[ring] / r_ring
    ny = ry[ring] / r_ring

    # Outward flux component at each ring point
    S_outward = Sx[ring] * nx + Sy[ring] * ny

    # Force = NEGATIVE of net outward flux (inward momentum transfer)
    # If more flux exits to the right → WC pushed LEFT
    # So Fx = -Σ(S_outward · nx), Fy = -Σ(S_outward · ny)
    Fx = -np.sum(S_outward * nx) * dx  # approximate line integral
    Fy = -np.sum(S_outward * ny) * dx

    return Fx, Fy


# ============================================================
# MAIN
# ============================================================

def main():
    print("=" * 70)
    print("  PHASE 1d :  2D Flux-Based Force (Radiation Pressure)")
    print("=" * 70)
    print(f"  Grid: {N_GRID}², {GRID_LAM}λ, {LAM0/DX:.1f} vox/λ")
    print(f"  Time samples: {N_TIME_SAMPLES} per period")
    print(f"  Separations: {SEPS_LAM[0]:.1f} to {SEPS_LAM[-1]:.1f}λ (step {SEPS_LAM[1]-SEPS_LAM[0]:.2f}λ)")
    print(f"  Sample radius for flux ring: 1.0λ from WC center")
    print()

    gx, gy = build_grid_2d(N_GRID, DX)
    period = 1.0 / F0
    times = np.linspace(0, period, N_TIME_SAMPLES, endpoint=False)
    sample_r = 1.0 * LAM0  # flux ring at 1λ from WC

    for phase1, phase2, label in PHASE_CONFIGS:
        coulomb_dir = "REP" if label == "same" else "ATT"

        print(f"\n{'#' * 70}")
        print(f"  {label.upper()} charge (expect {coulomb_dir})")
        print(f"{'#' * 70}")

        flux_dirs = []
        grad_dirs = []
        flux_forces = []

        t0 = time.time()

        for sep_lam in SEPS_LAM:
            sep = sep_lam * LAM0
            wc1_x, wc1_y = -sep / 2, 0.0
            wc2_x, wc2_y = +sep / 2, 0.0

            # Time-averaged flux force
            Fx_avg = 0.0

            for t in times:
                psi1 = wave_field_2d(gx, gy, wc1_x, wc1_y, phase1, A0, t)
                psi2 = wave_field_2d(gx, gy, wc2_x, wc2_y, phase2, A0, t)
                psi_total = psi1 + psi2

                Sx, Sy = compute_flux_2d(psi_total, DX)
                Fx, Fy = radiation_force_on_wc(Sx, Sy, gx, gy, wc1_x, wc1_y, DX, sample_r)
                Fx_avg += Fx

            Fx_avg /= N_TIME_SAMPLES

            # Direction: positive Fx = toward WC2 = ATT
            d_flux = "ATT" if Fx_avg > 0 else "REP"
            flux_dirs.append(d_flux)
            flux_forces.append(Fx_avg)

        dt = time.time() - t0

        # Summary
        n = len(SEPS_LAM)
        fa = flux_dirs.count("ATT")
        fr = flux_dirs.count("REP")
        fc = "CONSISTENT" if (fa == n or fr == n) else "MIXED"
        fd = "ATT" if fa > fr else "REP"
        match = "COULOMB" if (fc == "CONSISTENT" and fd == coulomb_dir) else ""
        if fc == "CONSISTENT" and fd != coulomb_dir:
            match = "INVERTED"

        print(f"  Computed in {dt:.1f}s")
        print(f"  Flux 2D: {fa:2d}/{n}A {fr:2d}/{n}R  {fc:>10s}  dom={fd}  {match}")
        print()

        # Print per-separation details
        print(f"  {'sep':>5s}  {'Fx_flux':>12s}  {'dir':>4s}  {'correct':>7s}")
        print(f"  {'─'*5}  {'─'*12}  {'─'*4}  {'─'*7}")
        for i, sep_lam in enumerate(SEPS_LAM):
            correct = "YES" if flux_dirs[i] == coulomb_dir else "no"
            print(f"  {sep_lam:5.2f}  {flux_forces[i]:12.4e}  {flux_dirs[i]:>4s}  {correct:>7s}")

    print(f"\n{'=' * 70}")
    print(f"  2D FLUX RESULTS")
    print(f"{'=' * 70}")
    print(f"  Does 2D radiation pressure give Coulomb?")
    print(f"  same→REP at ALL separations = electric repulsion")
    print(f"  opp→ATT at ALL separations = electric attraction")
    print(f"{'=' * 70}")


if __name__ == "__main__":
    main()
