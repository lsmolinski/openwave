#!/usr/bin/env python3
"""
Phase 1d: Variable λ(r) — 2D Cross-Section Force Test
======================================================
2D slice through the midplane of two WCs. Captures off-axis interference
that 1D misses, without the cubic memory explosion of 3D.

Grid: 1024² at ~4 vox/λ covers 256λ — fits K=10 particles (radius 100λ)
with near-field separations and some far-field.

Tests K=1 and K=10 with variable-λ energy equation.

Usage: python step1d_variable_lambda_2d.py
"""

import numpy as np
from scipy.integrate import cumulative_trapezoid
import time

# ============================================================
# PHYSICAL CONSTANTS
# ============================================================
A0   = 0.92
LAM0 = 28.5
F0   = 0.0105
RHO0 = 38.6
C    = 0.3
K0   = 2 * np.pi / LAM0

# ============================================================
# GRID
# ============================================================
N_GRID   = 1024   # 2D: 1024² = 1M points (manageable)
RNG_SEED = 42

# ============================================================
# CONFIGS
# ============================================================
K_VALUES = [1, 10]

PHASE_CONFIGS = [
    (0.0,    0.0,    "same (0,0)"),
    (np.pi,  np.pi,  "same (π,π)"),
    (0.0,    np.pi,  "opp (0,π)"),
    (np.pi,  0.0,    "opp (π,0)"),
]


# ============================================================
# λ(r) AND WKB
# ============================================================

def lambda_profile(r, K):
    """Yee & Hauger λ(r) shell model."""
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
    """Precompute WKB phase table."""
    r_f = np.linspace(0, r_max, n)
    k_f = 2*np.pi / lambda_profile(r_f, K)
    phi_f = np.zeros(n)
    phi_f[1:] = cumulative_trapezoid(k_f, r_f)
    return r_f, phi_f


# ============================================================
# 2D GRID AND BASE WAVE
# ============================================================

def build_grid_2d(n, dx):
    """2D grid on the xy-plane (z=0), centered at origin."""
    ax = (np.arange(n) - n/2.0 + 0.5) * dx
    gx, gy = np.meshgrid(ax, ax, indexing='ij')
    return gx, gy


def compute_base_phasor_2d(gx, gy, n_sources, seed, a0, k0):
    """
    Base wave on 2D plane: project N 3D plane waves onto z=0.
    Each plane wave k_hat contributes: amp * exp(-i*(k0*(kx*x + ky*y) + phi))
    """
    rng = np.random.default_rng(seed)
    # Fibonacci sphere for 3D directions
    idx = np.arange(n_sources, dtype=float)
    phi_s = np.arccos(1.0 - 2.0*(idx + 0.5)/n_sources)
    theta_s = np.pi*(1.0 + np.sqrt(5.0))*idx
    kx = np.sin(phi_s)*np.cos(theta_s)
    ky = np.sin(phi_s)*np.sin(theta_s)
    phases = rng.uniform(0, 2*np.pi, n_sources)

    amp = a0 / np.sqrt(n_sources)
    P = np.zeros_like(gx, dtype=np.complex128)
    for n in range(n_sources):
        kr = k0 * (kx[n]*gx + ky[n]*gy)
        P += amp * np.exp(-1j * (kr + phases[n]))
    return P


def compute_wc_outwave_2d(gx, gy, wc_x, wc_y, wc_phase, K, amplitude,
                          r_table, phi_table):
    """WC spherical out-wave on 2D plane (scalar, WKB phase)."""
    rx = gx - wc_x
    ry = gy - wc_y
    r = np.sqrt(rx**2 + ry**2)
    r_safe = np.maximum(r, LAM0*0.01)

    phi = np.interp(r_safe, r_table, phi_table)
    phi_safe = np.where(phi > 1e-10, phi, 1e-10)
    envelope = amplitude * np.sin(phi_safe) / phi_safe
    envelope = np.where(phi < 1e-10, amplitude, envelope)

    return envelope * np.exp(1j * (phi + wc_phase))


# ============================================================
# ENERGY AND FORCE (2D)
# ============================================================

def energy_var_2d(P, V, lam_local):
    A_rms_sq = np.abs(P)**2 / 2.0
    return RHO0 * V * C**2 * A_rms_sq / lam_local**2


def energy_const_2d(P, V):
    A_rms_sq = np.abs(P)**2 / 2.0
    return RHO0 * V * F0**2 * A_rms_sq


def force_2d(E, dx):
    """F = -∇E on 2D grid, returns (Nx, Ny, 2)."""
    Fx = np.zeros_like(E)
    Fy = np.zeros_like(E)
    Fx[1:-1, :] = -(E[2:, :] - E[:-2, :]) / (2*dx)
    Fy[:, 1:-1] = -(E[:, 2:] - E[:, :-2]) / (2*dx)
    return np.stack([Fx, Fy], axis=-1)


def sample_force_2d(F, gx, dx, point_x, point_y, core_r):
    """Average force over core-sized region around a point."""
    n = gx.shape[0]
    origin = gx[0, 0]
    ix = int(round((point_x - origin) / dx))
    iy = int(round((point_y - origin) / dx))
    ix = np.clip(ix, 3, n-4)
    iy = np.clip(iy, 3, n-4)
    hw = max(1, min(int(core_r / dx), 5))
    sl = (slice(ix-hw, ix+hw+1), slice(iy-hw, iy+hw+1))
    return F[sl].mean(axis=(0, 1))


# ============================================================
# MAIN
# ============================================================

def main():
    print("=" * 70)
    print("  PHASE 1d :  Variable λ(r) — 2D Cross-Section Force Test")
    print("=" * 70)

    for K in K_VALUES:
        r_particle = K**2 * LAM0

        # Scale grid to K
        if K == 1:
            grid_lam = 16.0       # 16λ, gives 4 vox/λ at 1024
            seps_lam = [2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0]
        else:
            grid_lam = 256.0      # 256λ, gives 4 vox/λ at 1024
            seps_lam = [105.0, 110.0, 115.0, 120.0, 130.0, 140.0, 150.0, 160.0, 180.0, 200.0]

        L = grid_lam * LAM0
        DX = L / N_GRID
        V_VOXEL = DX**2 * DX  # 2D "slab" volume

        print(f"\n{'#' * 70}")
        print(f"  K = {K}  ({'neutrino' if K == 1 else 'electron'})")
        print(f"  Particle radius : {r_particle:.0f} am  ({K**2}λ)")
        print(f"  Grid            : {N_GRID}², {grid_lam}λ, dx={DX:.2f} am ({LAM0/DX:.1f} vox/λ)")
        print(f"  Separations     : {seps_lam} λ")
        print(f"{'#' * 70}")

        gx, gy = build_grid_2d(N_GRID, DX)

        t0 = time.time()
        P_base = compute_base_phasor_2d(gx, gy, 200, RNG_SEED, A0, K0)
        print(f"  Base wave: {time.time()-t0:.1f}s")

        r_max_wkb = grid_lam * LAM0 * 1.5
        r_table, phi_table = build_wkb_lookup(K, r_max_wkb)

        print(f"\n  {'sep':>6s}  {'config':>14s}  "
              f"{'Fx_var':>10s}  {'Fx_const':>10s}  {'dir_v':>5s}  {'dir_c':>5s}  {'match':>5s}")
        print(f"  {'─'*6}  {'─'*14}  {'─'*10}  {'─'*10}  {'─'*5}  {'─'*5}  {'─'*5}")

        results_var = {}
        results_const = {}

        for sep_lam in seps_lam:
            sep = sep_lam * LAM0
            wc1_x, wc1_y = -sep/2, 0.0
            wc2_x, wc2_y = +sep/2, 0.0

            r1 = np.sqrt((gx - wc1_x)**2 + (gy - wc1_y)**2)
            r2 = np.sqrt((gx - wc2_x)**2 + (gy - wc2_y)**2)

            lam1 = lambda_profile(r1, K)
            lam2 = lambda_profile(r2, K)
            lam_local = np.where(r1 <= r2, lam1, lam2)

            for phase1, phase2, label in PHASE_CONFIGS:
                P_wc1 = compute_wc_outwave_2d(gx, gy, wc1_x, wc1_y, phase1, K, A0,
                                              r_table, phi_table)
                P_wc2 = compute_wc_outwave_2d(gx, gy, wc2_x, wc2_y, phase2, K, A0,
                                              r_table, phi_table)

                P_total = P_base + P_wc1 + P_wc2

                E_var = energy_var_2d(P_total, V_VOXEL, lam_local)
                E_const = energy_const_2d(P_total, V_VOXEL)

                F_var = force_2d(E_var, DX)
                F_const = force_2d(E_const, DX)

                r_core = K * LAM0
                Fv = sample_force_2d(F_var, gx, DX, wc1_x, wc1_y, r_core)
                Fc = sample_force_2d(F_const, gx, DX, wc1_x, wc1_y, r_core)

                Fxv = Fv[0]
                Fxc = Fc[0]

                dir_v = "ATT" if Fxv > 0 else "REP"
                dir_c = "ATT" if Fxc > 0 else "REP"
                match = "YES" if dir_v == dir_c else " NO"

                print(f"  {sep_lam:6.1f}  {label:>14s}  "
                      f"{Fxv:10.3e}  {Fxc:10.3e}  {dir_v:>5s}  {dir_c:>5s}  {match:>5s}")

                if label not in results_var:
                    results_var[label] = []
                    results_const[label] = []
                results_var[label].append(dir_v)
                results_const[label].append(dir_c)

        # Summary
        print(f"\n  {'─' * 66}")
        print(f"  DIRECTION CONSISTENCY (K={K}, 2D)")
        print(f"  {'─' * 66}")
        for label in results_var:
            dv = results_var[label]
            dc = results_const[label]
            n = len(dv)
            av, rv = dv.count("ATT"), dv.count("REP")
            ac, rc = dc.count("ATT"), dc.count("REP")
            cv = "CONSISTENT" if (av == n or rv == n) else "MIXED"
            cc = "CONSISTENT" if (ac == n or rc == n) else "MIXED"
            print(f"  {label:>14s}  var-λ: {av}/{n}A {rv}/{n}R {cv:>10s}  "
                  f"const-λ: {ac}/{n}A {rc}/{n}R {cc:>10s}")

        # Charge sensitivity: do same and opposite DIFFER?
        print(f"\n  CHARGE SENSITIVITY (K={K}):")
        for i, sep_lam in enumerate(seps_lam):
            dirs = {}
            for label in results_var:
                dirs[label] = results_var[label][i]
            same_dir = dirs["same (0,0)"]
            opp_dir = dirs["opp (0,π)"]
            charge_sens = "CHARGE-DEPENDENT" if same_dir != opp_dir else "charge-blind"
            print(f"    sep={sep_lam:6.1f}λ  same={same_dir}  opp={opp_dir}  → {charge_sens}")

    print(f"\n{'=' * 70}")
    print(f"  2D CROSS-SECTION RESULTS")
    print(f"{'=' * 70}")
    print(f"  Does 2D off-axis integration break sinc and reveal charge?")
    print(f"  CHARGE-DEPENDENT = same-phase and opposite-phase get different directions")
    print(f"{'=' * 70}")


if __name__ == "__main__":
    main()
