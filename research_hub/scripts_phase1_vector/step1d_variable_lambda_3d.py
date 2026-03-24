#!/usr/bin/env python3
"""
Phase 1d: Variable λ(r) — 3D Force Test
========================================
Tests whether 3D spherical integration of energy gradients breaks the sinc
oscillation that dominates the 1D on-axis force.

Hypothesis: on-axis, sinc nodes flip force every λ/2. Off-axis, the path
difference is d·cos(θ) → different oscillation period. Integrating -∇E
over all directions averages out the sinc flips, producing smooth net force.

Setup:
  - Base wave: 200 Fibonacci sphere plane waves (from Phase 1c Step 1)
  - Two WCs with WKB out-waves (variable λ from Yee & Hauger shells)
  - Energy: E = ρ·V·(c·A/λ_local)² — variable-λ equation
  - Force: F = -∇E via 3D central differences (captures ALL directions)
  - Compare K=1 (neutrino) vs K=10 (electron)

Usage: python step1d_variable_lambda_3d.py
"""

import numpy as np
from scipy.integrate import cumulative_trapezoid
import time

# ============================================================
# PHYSICAL CONSTANTS (sim units: am, rs, qg)
# ============================================================
A0   = 0.92
LAM0 = 28.5
F0   = 0.0105
RHO0 = 38.6
C    = 0.3
K0   = 2 * np.pi / LAM0

# ============================================================
# GRID PARAMETERS — scaled by K
# ============================================================
N_GRID   = 64
RNG_SEED = 42
N_SOURCES = 200

# ============================================================
# K VALUES AND SEPARATIONS
# ============================================================
K_VALUES = [1, 10]

# Separations in multiples of particle radius (K²λ)
# For K=1: r_particle = 1λ, seps = 3λ to 6λ (far-field)
# For K=10: r_particle = 100λ, grid can't fit far-field
#   → use near-field seps (within/near particle radius)
#   → grid extent = 8λ*scale, need seps < grid/2

PHASE_CONFIGS = [
    (0.0,    0.0,    "same (0,0)"),
    (np.pi,  np.pi,  "same (π,π)"),
    (0.0,    np.pi,  "opp (0,π)"),
    (np.pi,  0.0,    "opp (π,0)"),
]


# ============================================================
# HELPER FUNCTIONS
# ============================================================

def fibonacci_sphere(n):
    idx = np.arange(n, dtype=float)
    phi = np.arccos(1.0 - 2.0 * (idx + 0.5) / n)
    theta = np.pi * (1.0 + np.sqrt(5.0)) * idx
    return np.column_stack([
        np.sin(phi) * np.cos(theta),
        np.sin(phi) * np.sin(theta),
        np.cos(phi),
    ])


def build_grid(n_grid, dx):
    ax = (np.arange(n_grid) - n_grid / 2.0 + 0.5) * dx
    gx, gy, gz = np.meshgrid(ax, ax, ax, indexing='ij')
    return np.stack([gx, gy, gz], axis=-1)


def compute_base_phasor(coords, k_dirs, phases, a0, k0):
    """Base wave: N isotropic plane in-waves (scalar for speed)."""
    n_src = len(k_dirs)
    amp = a0 / np.sqrt(n_src)
    P = np.zeros(coords.shape[:-1], dtype=np.complex128)
    for n in range(n_src):
        kr = k0 * np.einsum('...j,j', coords, k_dirs[n])
        P += amp * np.exp(-1j * (kr + phases[n]))
    return P


# ============================================================
# λ(r) PROFILE AND WKB PHASE
# ============================================================

def lambda_profile(r, K):
    """Local wavelength from Yee & Hauger shell model."""
    r = np.asarray(r, dtype=float)
    r_core = K * LAM0
    r_particle = K**2 * LAM0
    lam = np.full_like(r, LAM0)

    core_mask = r < r_core
    lam[core_mask] = K * LAM0

    shell_mask = (r >= r_core) & (r <= r_particle)
    if shell_mask.any():
        r_s = r[shell_mask]
        b = -(2 * K - 1)
        c_coeff = r_s / LAM0 - K
        disc = b**2 - 4 * c_coeff
        disc = np.maximum(disc, 0.0)
        x = (-b - np.sqrt(disc)) / 2.0
        x = np.clip(x, 0, K)
        lam[shell_mask] = np.maximum(2 * (K - x) * LAM0, LAM0 * 0.1)

    return lam


def build_wkb_lookup(K, r_max, n_points=5000):
    """Precompute WKB phase lookup table for a given K."""
    r_fine = np.linspace(0, r_max, n_points)
    lam_fine = lambda_profile(r_fine, K)
    k_fine = 2 * np.pi / lam_fine
    phi_fine = np.zeros(n_points)
    phi_fine[1:] = cumulative_trapezoid(k_fine, r_fine)
    return r_fine, phi_fine


def wkb_phase_lookup(r, r_table, phi_table):
    """Interpolate WKB phase from precomputed lookup."""
    return np.interp(r, r_table, phi_table)


# ============================================================
# WC OUT-WAVE (3D SCALAR WITH WKB PHASE)
# ============================================================

def compute_wc_outwave_3d(coords, wc_pos, wc_phase, K, amplitude,
                          r_table, phi_table):
    """
    WC spherical out-wave with WKB phase on 3D grid.
    Scalar phasor (no L/T split — testing variable-λ effect on energy).
    """
    r_vec = coords - wc_pos
    r_mag = np.linalg.norm(r_vec, axis=-1)
    r_safe = np.maximum(r_mag, LAM0 * 0.01)

    # WKB phase from lookup
    phi = wkb_phase_lookup(r_safe, r_table, phi_table)

    # Modified sinc: sin(φ_wkb)/φ_wkb
    phi_safe = np.where(phi > 1e-10, phi, 1e-10)
    envelope = amplitude * np.sin(phi_safe) / phi_safe
    envelope = np.where(phi < 1e-10, amplitude, envelope)

    # Out-wave phasor
    return envelope * np.exp(1j * (phi + wc_phase))


# ============================================================
# ENERGY AND FORCE (3D)
# ============================================================

def energy_variable_lambda_3d(P, V, lam_local):
    """E = ρ·V·(c·A_rms/λ)²"""
    A_rms_sq = np.abs(P)**2 / 2.0
    return RHO0 * V * C**2 * A_rms_sq / lam_local**2


def energy_constant_lambda_3d(P, V):
    """E = ρ·V·(f₀·A_rms)²"""
    A_rms_sq = np.abs(P)**2 / 2.0
    return RHO0 * V * F0**2 * A_rms_sq


def force_3d(E, dx):
    """F = -∇E via central differences, returns (Nx,Ny,Nz,3)."""
    F = np.zeros(E.shape + (3,))
    F[1:-1, :, :, 0] = -(E[2:, :, :] - E[:-2, :, :]) / (2.0 * dx)
    F[:, 1:-1, :, 1] = -(E[:, 2:, :] - E[:, :-2, :]) / (2.0 * dx)
    F[:, :, 1:-1, 2] = -(E[:, :, 2:] - E[:, :, :-2]) / (2.0 * dx)
    return F


def sample_force(F, coords, point, dx, n_grid, core_radius):
    """Sample force at grid point closest to position, averaged over core."""
    idx = np.round((point - coords[0, 0, 0]) / dx).astype(int)
    idx = np.clip(idx, 3, n_grid - 4)
    hw = max(1, int(core_radius / dx))
    hw = min(hw, 3)  # cap window to avoid edge issues
    sl = (slice(idx[0]-hw, idx[0]+hw+1),
          slice(idx[1]-hw, idx[1]+hw+1),
          slice(idx[2]-hw, idx[2]+hw+1))
    return F[sl].mean(axis=(0, 1, 2))


# ============================================================
# MAIN
# ============================================================

def main():
    print("=" * 70)
    print("  PHASE 1d :  Variable λ(r) — 3D Force Test")
    print("=" * 70)

    for K in K_VALUES:
        r_particle = K**2 * LAM0

        # Scale grid to particle size
        # K=1: 8λ grid (plenty of room)
        # K=10: need larger grid. Use scale factor to fit.
        if K == 1:
            grid_lam = 8.0
            seps_lam = [2.0, 3.0, 4.0, 5.0, 6.0]
        else:
            # K=10: particle radius = 100λ. Grid must be >> particle.
            # Use coarser resolution to fit: 256λ grid at ~4 vox/λ
            grid_lam = 256.0
            seps_lam = [120.0, 150.0, 180.0, 210.0, 240.0]

        L = grid_lam * LAM0
        DX = L / N_GRID
        V_VOXEL = DX**3

        print(f"\n{'#' * 70}")
        print(f"  K = {K}  ({'neutrino' if K == 1 else 'electron'})")
        print(f"  Particle radius : {r_particle:.0f} am  ({K**2} λ)")
        print(f"  Grid            : {N_GRID}³, {grid_lam}λ, dx={DX:.2f} am ({LAM0/DX:.1f} vox/λ)")
        print(f"  Separations     : {seps_lam} λ")
        print(f"{'#' * 70}")

        # ---- Build grid ----
        coords = build_grid(N_GRID, DX)

        # ---- Base wave ----
        k_dirs = fibonacci_sphere(N_SOURCES)
        rng = np.random.default_rng(RNG_SEED)
        phases = rng.uniform(0, 2 * np.pi, N_SOURCES)

        t0 = time.time()
        P_base = compute_base_phasor(coords, k_dirs, phases, A0, K0)
        print(f"  Base wave: {time.time()-t0:.1f}s")

        # ---- WKB lookup table ----
        r_max_wkb = grid_lam * LAM0 * 1.5
        r_table, phi_table = build_wkb_lookup(K, r_max_wkb)

        # ---- Force test ----
        print(f"\n  {'sep':>6s}  {'config':>14s}  "
              f"{'Fx_var':>10s}  {'Fx_const':>10s}  {'dir_v':>5s}  {'dir_c':>5s}  {'match':>5s}")
        print(f"  {'─'*6}  {'─'*14}  {'─'*10}  {'─'*10}  {'─'*5}  {'─'*5}  {'─'*5}")

        results_var = {}
        results_const = {}

        for sep_lam in seps_lam:
            sep = sep_lam * LAM0

            # WC positions symmetric on x-axis
            wc1_pos = np.array([-sep / 2, 0.0, 0.0])
            wc2_pos = np.array([+sep / 2, 0.0, 0.0])

            # Distance from each WC to every grid point
            r1_vec = coords - wc1_pos
            r2_vec = coords - wc2_pos
            r1 = np.linalg.norm(r1_vec, axis=-1)
            r2 = np.linalg.norm(r2_vec, axis=-1)

            # λ(r) field: nearest WC determines local λ
            lam1 = lambda_profile(r1, K)
            lam2 = lambda_profile(r2, K)
            lam_local = np.where(r1 <= r2, lam1, lam2)

            for phase1, phase2, label in PHASE_CONFIGS:
                # WC out-waves with WKB phase
                P_wc1 = compute_wc_outwave_3d(coords, wc1_pos, phase1, K, A0,
                                              r_table, phi_table)
                P_wc2 = compute_wc_outwave_3d(coords, wc2_pos, phase2, K, A0,
                                              r_table, phi_table)

                # Total scalar phasor
                P_total = P_base + P_wc1 + P_wc2

                # Energy: variable-λ and constant-λ
                E_var = energy_variable_lambda_3d(P_total, V_VOXEL, lam_local)
                E_const = energy_constant_lambda_3d(P_total, V_VOXEL)

                # 3D force
                F_var = force_3d(E_var, DX)
                F_const = force_3d(E_const, DX)

                # Sample force at WC1 (x-component: positive = toward WC2 = attraction)
                r_core = K * LAM0
                Fv = sample_force(F_var, coords, wc1_pos, DX, N_GRID, r_core)
                Fc = sample_force(F_const, coords, wc1_pos, DX, N_GRID, r_core)

                Fxv = Fv[0]
                Fxc = Fc[0]

                dir_v = "ATT" if Fxv > 0 else "REP"
                dir_c = "ATT" if Fxc > 0 else "REP"
                match = "YES" if dir_v == dir_c else " NO"

                print(f"  {sep_lam:6.0f}  {label:>14s}  "
                      f"{Fxv:10.3e}  {Fxc:10.3e}  {dir_v:>5s}  {dir_c:>5s}  {match:>5s}")

                if label not in results_var:
                    results_var[label] = []
                    results_const[label] = []
                results_var[label].append(dir_v)
                results_const[label].append(dir_c)

        # ---- Summary ----
        print(f"\n  {'─' * 66}")
        print(f"  DIRECTION CONSISTENCY (K={K}, 3D)")
        print(f"  {'─' * 66}")
        for label in results_var:
            dv = results_var[label]
            dc = results_const[label]
            n = len(dv)
            av = dv.count("ATT")
            rv = dv.count("REP")
            ac = dc.count("ATT")
            rc = dc.count("REP")
            cv = "CONSISTENT" if (av == n or rv == n) else "MIXED"
            cc = "CONSISTENT" if (ac == n or rc == n) else "MIXED"
            print(f"  {label:>14s}  var-λ: {av}/{n}A {rv}/{n}R {cv:>10s}  "
                  f"const-λ: {ac}/{n}A {rc}/{n}R {cc:>10s}")

    print(f"\n{'=' * 70}")
    print(f"  3D vs 1D COMPARISON")
    print(f"{'=' * 70}")
    print(f"  If 3D shows CONSISTENT where 1D showed MIXED → 3D averaging")
    print(f"  breaks the sinc oscillation. The Coulomb force is a 3D phenomenon.")
    print(f"  If still MIXED → need different λ(r) profile or additional mechanism.")
    print(f"{'=' * 70}")


if __name__ == "__main__":
    main()
