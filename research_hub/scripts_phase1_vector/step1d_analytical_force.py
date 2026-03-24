#!/usr/bin/env python3
"""
Phase 1d: Analytical Interaction Force Between Two WCs
======================================================
Isolates the WC-WC interaction energy from the base wave noise.

The total energy decomposes as:
  E_total = E_base + E_wc1 + E_wc2 + E_base×wc1 + E_base×wc2 + E_wc1×wc2

Only E_wc1×wc2 (the cross term) carries charge-dependent information.
F = -∇E_total is dominated by E_base×wc terms. To reveal the Coulomb force,
we compute F_int = -∇E_interaction where E_int = E_total - E_base_only - E_wc1_only - E_wc2_only.

Tests:
  1. Does the isolated interaction force depend on charge (source_offset)?
  2. Is the direction consistent across separations (no sinc flips)?
  3. Does 2D integration help vs 1D?

Usage: python step1d_analytical_force.py
"""

import numpy as np
from scipy.integrate import cumulative_trapezoid
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
K_VALUES = [1, 10]
N_GRID   = 512    # 2D grid

PHASE_CONFIGS = [
    (0.0,    0.0,    "same (0,0)"),
    (np.pi,  np.pi,  "same (π,π)"),
    (0.0,    np.pi,  "opp (0,π)"),
    (np.pi,  0.0,    "opp (π,0)"),
]


# ============================================================
# λ(r) AND WKB (from step1d)
# ============================================================

def lambda_profile(r, K):
    r = np.asarray(r, dtype=float)
    r_core = K * LAM0
    r_particle = K**2 * LAM0
    lam = np.full_like(r, LAM0)
    lam[r < r_core] = K * LAM0
    shell = (r >= r_core) & (r <= r_particle)
    if shell.any():
        rs = r[shell]
        b = -(2*K - 1)
        disc = np.maximum(b**2 - 4*(rs/LAM0 - K), 0.0)
        x = np.clip((-b - np.sqrt(disc)) / 2.0, 0, K)
        lam[shell] = np.maximum(2*(K - x)*LAM0, LAM0*0.1)
    return lam


def build_wkb_lookup(K, r_max, n=10000):
    r_f = np.linspace(0, r_max, n)
    k_f = 2*np.pi / lambda_profile(r_f, K)
    phi_f = np.zeros(n)
    phi_f[1:] = cumulative_trapezoid(k_f, r_f)
    return r_f, phi_f


# ============================================================
# 2D INFRASTRUCTURE
# ============================================================

def build_grid_2d(n, dx):
    ax = (np.arange(n) - n/2.0 + 0.5) * dx
    gx, gy = np.meshgrid(ax, ax, indexing='ij')
    return gx, gy


def wc_phasor_2d(gx, gy, wx, wy, phase, K, amp, r_table, phi_table):
    """WC out-wave phasor on 2D grid."""
    r = np.sqrt((gx - wx)**2 + (gy - wy)**2)
    r_safe = np.maximum(r, LAM0*0.01)
    phi = np.interp(r_safe, r_table, phi_table)
    phi_safe = np.where(phi > 1e-10, phi, 1e-10)
    env = amp * np.sin(phi_safe) / phi_safe
    env = np.where(phi < 1e-10, amp, env)
    return env * np.exp(1j * (phi + phase))


def energy_2d(P, V):
    """E = ρ·V·f₀²·A_rms² (constant-λ for interaction isolation)."""
    return RHO0 * V * F0**2 * np.abs(P)**2 / 2.0


def force_x_2d(E, dx):
    """dE/dx via central differences, return Fx field."""
    Fx = np.zeros_like(E)
    Fx[1:-1, :] = -(E[2:, :] - E[:-2, :]) / (2*dx)
    return Fx


def sample_force_x(Fx, gx, dx, point_x, point_y, core_r):
    """Average Fx over core region around a point."""
    n = gx.shape[0]
    origin = gx[0, 0]
    ix = int(round((point_x - origin) / dx))
    iy = int(round((point_y - origin) / dx))
    ix = np.clip(ix, 3, n-4)
    iy = np.clip(iy, 3, n-4)
    hw = max(1, min(int(core_r / dx), 5))
    sl = (slice(ix-hw, ix+hw+1), slice(iy-hw, iy+hw+1))
    return Fx[sl].mean()


# ============================================================
# MAIN
# ============================================================

def main():
    print("=" * 70)
    print("  PHASE 1d :  Analytical Interaction Force (Isolated Cross Term)")
    print("=" * 70)
    print(f"  Method: E_int = E(wc1+wc2) - E(wc1_alone) - E(wc2_alone)")
    print(f"  This isolates the charge-dependent WC×WC cross term")
    print(f"  from the charge-blind self-energy terms.")
    print()

    for K in K_VALUES:
        r_particle = K**2 * LAM0

        if K == 1:
            grid_lam = 16.0
            seps_lam = [2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0]
        else:
            grid_lam = 256.0
            seps_lam = [105.0, 110.0, 115.0, 120.0, 130.0, 140.0, 150.0, 180.0, 200.0]

        L = grid_lam * LAM0
        DX = L / N_GRID
        V = DX**2 * DX

        print(f"\n{'#' * 70}")
        print(f"  K = {K}  ({'neutrino' if K == 1 else 'electron'})")
        print(f"  Grid: {N_GRID}², {grid_lam}λ, {LAM0/DX:.1f} vox/λ")
        print(f"{'#' * 70}")

        gx, gy = build_grid_2d(N_GRID, DX)
        r_max_wkb = grid_lam * LAM0 * 1.5
        r_table, phi_table = build_wkb_lookup(K, r_max_wkb)
        r_core = K * LAM0

        print(f"\n  {'sep':>6s}  {'config':>14s}  {'Fx_int':>10s}  {'Fx_total':>10s}  "
              f"{'dir_int':>7s}  {'dir_tot':>7s}  {'charge':>8s}")
        print(f"  {'─'*6}  {'─'*14}  {'─'*10}  {'─'*10}  {'─'*7}  {'─'*7}  {'─'*8}")

        results_int = {}

        for sep_lam in seps_lam:
            sep = sep_lam * LAM0
            wc1_x, wc1_y = -sep/2, 0.0
            wc2_x, wc2_y = +sep/2, 0.0

            for phase1, phase2, label in PHASE_CONFIGS:
                # Compute individual and combined phasors (NO base wave)
                P_wc1 = wc_phasor_2d(gx, gy, wc1_x, wc1_y, phase1, K, A0,
                                     r_table, phi_table)
                P_wc2 = wc_phasor_2d(gx, gy, wc2_x, wc2_y, phase2, K, A0,
                                     r_table, phi_table)

                # Energies
                E_combined = energy_2d(P_wc1 + P_wc2, V)
                E_wc1_only = energy_2d(P_wc1, V)
                E_wc2_only = energy_2d(P_wc2, V)

                # Interaction energy: E_int = E(1+2) - E(1) - E(2)
                # = ρVf²·Re(P_wc1*·P_wc2) (the cross term)
                E_interaction = E_combined - E_wc1_only - E_wc2_only

                # Also compute total (with both) for comparison
                E_total = E_combined

                # Force from interaction energy only
                Fx_int = force_x_2d(E_interaction, DX)
                Fx_tot = force_x_2d(E_total, DX)

                # Sample at WC1
                fx_i = sample_force_x(Fx_int, gx, DX, wc1_x, wc1_y, r_core)
                fx_t = sample_force_x(Fx_tot, gx, DX, wc1_x, wc1_y, r_core)

                dir_i = "ATT" if fx_i > 0 else "REP"
                dir_t = "ATT" if fx_t > 0 else "REP"

                if label not in results_int:
                    results_int[label] = []
                results_int[label].append(dir_i)

                print(f"  {sep_lam:6.1f}  {label:>14s}  {fx_i:10.3e}  {fx_t:10.3e}  "
                      f"{dir_i:>7s}  {dir_t:>7s}")

        # ---- Charge sensitivity from interaction term ----
        print(f"\n  {'─' * 66}")
        print(f"  CHARGE SENSITIVITY FROM INTERACTION ENERGY (K={K})")
        print(f"  {'─' * 66}")

        for i, sep_lam in enumerate(seps_lam):
            d_same = results_int["same (0,0)"][i]
            d_opp = results_int["opp (0,π)"][i]
            charge = "CHARGE-DEP" if d_same != d_opp else "blind"
            print(f"    sep={sep_lam:6.1f}λ  same={d_same}  opp={d_opp}  → {charge}")

        # Direction consistency
        print(f"\n  DIRECTION CONSISTENCY (interaction force, K={K}):")
        for label in results_int:
            d = results_int[label]
            n = len(d)
            a, r = d.count("ATT"), d.count("REP")
            c = "CONSISTENT" if (a == n or r == n) else "MIXED"
            print(f"    {label:>14s}  {a}/{n}A {r}/{n}R  {c}")

    print(f"\n{'=' * 70}")
    print(f"  ANALYTICAL FORCE RESULTS")
    print(f"{'=' * 70}")
    print(f"  The interaction energy E_int = E(1+2) - E(1) - E(2)")
    print(f"  isolates the cross term Re(P₁*·P₂) that carries charge info.")
    print(f"  If CHARGE-DEP → the force IS there, just buried in E_total.")
    print(f"  If still blind → the sinc in the cross term itself is the blocker.")
    print(f"{'=' * 70}")


if __name__ == "__main__":
    main()
