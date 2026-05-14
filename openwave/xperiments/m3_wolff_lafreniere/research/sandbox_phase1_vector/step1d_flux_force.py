#!/usr/bin/env python3
"""
Phase 1d: Flux-Based Force (Radiation Pressure / Poynting-like)
===============================================================
LaFreniere's Coulomb mechanism: force from wave MOMENTUM FLUX,
not from energy density gradient (-∇E).

The key difference:
  -∇E: force from WHERE energy is stored (spatial gradient of scalar)
  Flux: force from WHERE energy is FLOWING (directional momentum transfer)

For waves: energy flux ∝ ψ · ∂ψ/∂t (displacement × velocity)
In phasor: flux ∝ Im(P* · ∇P) — the Poynting-like vector for scalar waves.
This is inherently directional — it knows which way the wave is traveling.

Standing waves: zero net flux (energy sloshes back and forth)
Traveling waves: nonzero net flux (energy flows outward)

The force on a WC is the net momentum flux through a surface surrounding it.
If more flux arrives from one side (due to the other WC's traveling wave),
the imbalance pushes the WC — this is radiation pressure.

For same phase: flux between WCs cancels (destructive) → less outward
  pressure between → particles pushed TOGETHER? Or apart?
For opposite phase: flux between WCs reinforces → more pressure between
  → particles pushed APART? Need to check.

Usage: python step1d_flux_force.py
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
OMEGA = 2 * np.pi * F0

# ============================================================
# CONFIGS
# ============================================================
N_SEPS = 50
N_PTS  = 50000

PHASE_CONFIGS = [
    (0.0,   0.0,   "same"),
    (0.0,   np.pi, "opp"),
]

K_VALUES = [1, 10]


# ============================================================
# WAVE FUNCTIONS (time-domain for flux computation)
# ============================================================

def weighted_partial_standing_wave(x, wc_x, phase, amplitude, t):
    """
    Full time-domain weighted partial standing wave from M3:
      ψ = A · [w·sin(kr + ωt + φ) + sin(kr - ωt - φ)] / kr

    Returns displacement ψ(x,t) and velocity ∂ψ/∂t(x,t).
    """
    r = np.abs(x - wc_x)
    r_safe = np.maximum(r, LAM0 * 0.01)
    kr = K0 * r_safe

    # Weight function
    transition = 1.25 * LAM0
    w = 1.0 / (1.0 + (r_safe / transition) ** 8)

    # Direction: sign of (x - wc_x) for 1D gradient
    # Not needed for scalar displacement, but needed for flux direction

    # Displacement
    in_wave = w * np.sin(kr + OMEGA * t + phase)
    out_wave = np.sin(kr - OMEGA * t - phase)
    psi = amplitude * (in_wave + out_wave) / kr

    # Velocity: ∂ψ/∂t
    d_in_wave = w * OMEGA * np.cos(kr + OMEGA * t + phase)
    d_out_wave = -OMEGA * np.cos(kr - OMEGA * t - phase)
    dpsi_dt = amplitude * (d_in_wave + d_out_wave) / kr

    return psi, dpsi_dt


def compute_energy_flux_1d(x, wc1_x, wc2_x, phase1, phase2, amplitude, t):
    """
    Compute energy flux (Poynting-like) along x-axis.

    For 1D scalar waves: flux S = -c² · ψ · (∂ψ/∂x)
    This is the acoustic intensity — energy flow per unit area.
    Positive = rightward flow, negative = leftward.

    Also compute: flux S = ψ · (∂ψ/∂t) — displacement × velocity
    (equivalent up to constants for propagating waves).

    We compute BOTH individual WC fluxes and the combined flux.
    The radiation pressure on WC1 is the NET flux imbalance at WC1's position.
    """
    dx = x[1] - x[0]

    # Combined field
    psi1, dpsi1_dt = weighted_partial_standing_wave(x, wc1_x, phase1, amplitude, t)
    psi2, dpsi2_dt = weighted_partial_standing_wave(x, wc2_x, phase2, amplitude, t)

    psi_total = psi1 + psi2
    dpsi_total_dt = dpsi1_dt + dpsi2_dt

    # Spatial gradient of total displacement
    dpsi_dx = np.gradient(psi_total, x)

    # Energy flux: S = -c² · ψ · ∂ψ/∂x (acoustic intensity)
    S_spatial = -C**2 * psi_total * dpsi_dx

    # Alternative: S = ψ · ∂ψ/∂t (displacement × velocity)
    S_temporal = psi_total * dpsi_total_dt

    return S_spatial, S_temporal, psi_total


def radiation_pressure_force(x, S, wc_x, core_radius):
    """
    Radiation pressure force on a WC:
    Net momentum flux through the WC = flux from right minus flux from left.

    F = S(x_wc + ε) - S(x_wc - ε)

    If more flux arrives from the right → pushed left (negative force).
    If more flux arrives from the left → pushed right (positive force).
    Positive force = toward +x.
    """
    idx = np.argmin(np.abs(x - wc_x))
    # Sample flux just outside the WC on both sides
    offset = max(3, int(core_radius / (x[1] - x[0])))

    idx_right = min(idx + offset, len(S) - 1)
    idx_left = max(idx - offset, 0)

    # Net flux: flux arriving from left (S at left, pointing right = positive)
    # minus flux arriving from right (S at right, pointing left = negative)
    # Radiation pressure: F ∝ S_left - S_right (net inward flux)
    # If S_left > S_right → net push to the right (toward +x)
    F_rad = S[idx_left] - S[idx_right]

    return F_rad


# ============================================================
# TIME-AVERAGED FLUX
# ============================================================

def time_averaged_flux_force(x, wc1_x, wc2_x, phase1, phase2, amplitude, n_samples=32):
    """
    Average the flux-based force over one full wave period.
    The flux oscillates at 2ω; time averaging extracts the DC component.
    """
    period = 1.0 / F0  # one full period in rs
    times = np.linspace(0, period, n_samples, endpoint=False)
    core_r = LAM0

    F_spatial_sum = 0.0
    F_temporal_sum = 0.0

    for t in times:
        S_s, S_t, _ = compute_energy_flux_1d(x, wc1_x, wc2_x, phase1, phase2, amplitude, t)
        F_spatial_sum += radiation_pressure_force(x, S_s, wc1_x, core_r)
        F_temporal_sum += radiation_pressure_force(x, S_t, wc1_x, core_r)

    return F_spatial_sum / n_samples, F_temporal_sum / n_samples


# ============================================================
# MAIN
# ============================================================

def main():
    print("=" * 70)
    print("  PHASE 1d :  Flux-Based Force (Radiation Pressure)")
    print("=" * 70)
    print(f"  Force = net energy flux imbalance at WC position")
    print(f"  Time-averaged over one full period (32 samples)")
    print(f"  Two flux definitions: S_spatial = -c²·ψ·∂ψ/∂x")
    print(f"                        S_temporal = ψ·∂ψ/∂t")
    print()

    for K in K_VALUES:
        r_particle = K**2 * LAM0
        if K == 1:
            seps_lam = np.linspace(2.0, 7.0, N_SEPS)
        else:
            seps_lam = np.linspace(105.0, 200.0, N_SEPS)

        print(f"\n{'#' * 70}")
        print(f"  K = {K}  seps = {seps_lam[0]:.0f}-{seps_lam[-1]:.0f}λ")
        print(f"{'#' * 70}")

        for phase1, phase2, label in PHASE_CONFIGS:
            coulomb = "REP" if label == "same" else "ATT"

            spatial_dirs = []
            temporal_dirs = []
            gradient_dirs = []

            x = np.linspace(-20 * LAM0, (seps_lam[-1] + 20) * LAM0, N_PTS)

            for sep_lam in seps_lam:
                sep = sep_lam * LAM0

                # Flux-based force (time-averaged)
                F_s, F_t = time_averaged_flux_force(x, 0.0, sep, phase1, phase2, A0)
                spatial_dirs.append("ATT" if F_s > 0 else "REP")
                temporal_dirs.append("ATT" if F_t > 0 else "REP")

                # Gradient-based force (for comparison, from interaction energy)
                r1 = np.abs(x)
                r2 = np.abs(x - sep)
                kr1 = K0 * np.maximum(r1, LAM0 * 0.01)
                kr2 = K0 * np.maximum(r2, LAM0 * 0.01)
                Z1 = -1j * A0 * np.exp(1j * (kr1 - phase1)) / kr1
                Z2 = -1j * A0 * np.exp(1j * (kr2 - phase2)) / kr2
                dx = x[1] - x[0]
                E_int = RHO0 * dx * F0**2 * (np.abs(Z1+Z2)**2 - np.abs(Z1)**2 - np.abs(Z2)**2) / 2
                F_grad = -np.gradient(E_int, x)
                idx = np.argmin(np.abs(x))
                gradient_dirs.append("ATT" if F_grad[idx] > 0 else "REP")

            # Count
            sa = spatial_dirs.count("ATT")
            sr = spatial_dirs.count("REP")
            ta = temporal_dirs.count("ATT")
            tr = temporal_dirs.count("REP")
            ga = gradient_dirs.count("ATT")
            gr = gradient_dirs.count("REP")

            sc = "CONSISTENT" if (sa == N_SEPS or sr == N_SEPS) else "MIXED"
            tc = "CONSISTENT" if (ta == N_SEPS or tr == N_SEPS) else "MIXED"
            gc = "CONSISTENT" if (ga == N_SEPS or gr == N_SEPS) else "MIXED"

            sd = "ATT" if sa > sr else "REP"
            td = "ATT" if ta > tr else "REP"
            gd = "ATT" if ga > gr else "REP"

            s_match = "COULOMB" if (sc == "CONSISTENT" and sd == coulomb) else ("INVERTED" if sc == "CONSISTENT" else "")
            t_match = "COULOMB" if (tc == "CONSISTENT" and td == coulomb) else ("INVERTED" if tc == "CONSISTENT" else "")

            print(f"\n  {label.upper()} charge (expect {coulomb}):")
            print(f"    Flux(spatial):  {sa:2d}/{N_SEPS}A {sr:2d}/{N_SEPS}R  {sc:>10s}  dom={sd}  {s_match}")
            print(f"    Flux(temporal): {ta:2d}/{N_SEPS}A {tr:2d}/{N_SEPS}R  {tc:>10s}  dom={td}  {t_match}")
            print(f"    Gradient(-∇E):  {ga:2d}/{N_SEPS}A {gr:2d}/{N_SEPS}R  {gc:>10s}  dom={gd}")

    print(f"\n{'=' * 70}")
    print(f"  KEY QUESTION: Does flux-based force give consistent Coulomb?")
    print(f"  same→REP + opp→ATT = Coulomb correct")
    print(f"  CONSISTENT flux where gradient is MIXED → flux IS the mechanism")
    print(f"{'=' * 70}")


if __name__ == "__main__":
    main()
