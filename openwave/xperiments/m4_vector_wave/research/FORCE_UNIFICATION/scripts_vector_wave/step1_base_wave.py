#!/usr/bin/env python3
"""
Phase 1c Step 1: 3D Vector Base Wave
=====================================
Validates that an isotropic base wave (N plane waves from all directions)
produces:
  1. Mean energy matches theory: E = rho*V*f^2*A0^2/2
  2. Energy CV = 1/sqrt(3) ~ 0.577 (chi-squared with 6 DOF -- speckle)
  3. E_L/E_total = 1/3, E_T/E_total = 2/3 (isotropic L/T split)
  4. E_L + E_T = E_total exactly at every point
  5. Force is speckle noise only (no large-scale gradient)

Physics:
  Each plane wave is longitudinal — displacement along propagation direction k_hat.
  Isotropic superposition from all directions creates vector displacement at every
  point. L/T decomposition relative to any reference point yields 1/3 longitudinal
  (one radial direction) and 2/3 transverse (two perpendicular directions).

  Energy statistics: the phasor components P_x, P_y, P_z converge to independent
  complex Gaussians (CLT). |P|^2 = sum of 6 squared Gaussians = chi-squared(6).
  This gives CV = 1/sqrt(3) ~ 0.577 INDEPENDENT of N_sources — the "speckle"
  pattern is fundamental to coherent wave interference, not a numerical artifact.
  At particle scales (electron = 100 lam), these fluctuations average out.

Usage: python step1_base_wave.py

Research: Phase 1c Vector Wave Force (FORCE_UNIFICATION/01c_vector_wave.md)
Strategy: math-only numpy, no visualization.
"""

import numpy as np
import time

# ============================================================
# PHYSICAL CONSTANTS (sim units: am, rs, qg)
# ============================================================
A0   = 0.92       # Fundamental amplitude [am]
LAM0 = 28.5       # Fundamental wavelength [am]
F0   = 0.0105     # Fundamental frequency [rHz]
RHO0 = 38.6       # Medium density [qg/am^3]
C    = 0.3        # Wave speed [am/rs]
K0   = 2 * np.pi / LAM0  # Wavenumber [1/am]

# ============================================================
# GRID PARAMETERS
# ============================================================
N_GRID   = 64      # Points per axis (64^3 = 262,144 voxels)
GRID_LAM = 8.0     # Grid extent in wavelengths
L        = GRID_LAM * LAM0   # Physical size [am]
DX       = L / N_GRID        # Grid spacing [am]
V_VOXEL  = DX ** 3           # Voxel volume [am^3]

# ============================================================
# BASE WAVE PARAMETERS
# ============================================================
N_SOURCES = 200    # Number of isotropic plane wave sources
RNG_SEED  = 42     # For reproducibility


# ============================================================
# HELPER FUNCTIONS
# ============================================================

def fibonacci_sphere(n):
    """Generate n approximately uniform points on unit sphere."""
    idx = np.arange(n, dtype=float)
    # Distribute points uniformly in cos(theta)
    phi = np.arccos(1.0 - 2.0 * (idx + 0.5) / n)
    # Golden angle spacing in azimuth
    theta = np.pi * (1.0 + np.sqrt(5.0)) * idx
    return np.column_stack([
        np.sin(phi) * np.cos(theta),
        np.sin(phi) * np.sin(theta),
        np.cos(phi),
    ])


def build_grid(n_grid, dx):
    """Build 3D coordinate grid centered at origin. Returns (N,N,N,3)."""
    ax = (np.arange(n_grid) - n_grid / 2.0 + 0.5) * dx
    gx, gy, gz = np.meshgrid(ax, ax, ax, indexing='ij')
    return np.stack([gx, gy, gz], axis=-1)


# ============================================================
# CORE PHYSICS
# ============================================================

def compute_phasor_field(coords, k_dirs, phases, a0, k0):
    """
    Compute complex vector phasor from isotropic longitudinal plane waves.

    Each source n:
      psi_n(r,t) = (a0/sqrt(N)) * cos(k0*k_hat_n.r + phi_n - wt) * k_hat_n

    Complex phasor (dropping shared -wt):
      Z_n(r) = (a0/sqrt(N)) * exp(i*(k0*k_hat_n.r + phi_n)) * k_hat_n

    Returns: P(r) = sum_n Z_n(r),  shape (*grid, 3) complex128
    """
    n_src = len(k_dirs)
    amp = a0 / np.sqrt(n_src)
    P = np.zeros(coords.shape, dtype=np.complex128)

    for n in range(n_src):
        # k_hat_n . r  at every grid point
        kr = k0 * np.einsum('...j,j', coords, k_dirs[n])
        # Complex amplitude (scalar per grid point)
        z = amp * np.exp(1j * (kr + phases[n]))
        # Displacement along k_hat_n (longitudinal)
        P += z[..., np.newaxis] * k_dirs[n]

    return P


def energy_from_phasor(P):
    """
    Energy per voxel from vector phasor.
    E(r) = rho * V * f^2 * |P|^2 / 2
    The /2 converts peak amplitude to RMS: A_rms = |P|/sqrt(2).
    """
    amp_sq = np.sum(np.abs(P) ** 2, axis=-1)   # |Px|^2 + |Py|^2 + |Pz|^2
    return RHO0 * V_VOXEL * F0 ** 2 * amp_sq / 2.0


def force_from_energy(E, dx):
    """F = -grad(E) via central differences. Boundary left at zero."""
    F = np.zeros(E.shape + (3,))
    F[1:-1, :, :, 0] = -(E[2:, :, :] - E[:-2, :, :]) / (2.0 * dx)
    F[:, 1:-1, :, 1] = -(E[:, 2:, :] - E[:, :-2, :]) / (2.0 * dx)
    F[:, :, 1:-1, 2] = -(E[:, :, 2:] - E[:, :, :-2]) / (2.0 * dx)
    return F


def lt_decomposition(P, ref_point, coords):
    """
    Decompose vector phasor into longitudinal and transverse components
    relative to ref_point.

    r_hat = (r - r_ref) / |r - r_ref|
    A_L = |P . r_hat|           (projection onto radial)
    A_T = |P - (P.r_hat)r_hat|  (rejection from radial)

    Returns: A_L_sq, A_T_sq, valid_mask  (excludes ref_point vicinity)
    """
    dr = coords - ref_point
    dist = np.linalg.norm(dr, axis=-1, keepdims=True)

    # Avoid division by zero at reference point
    safe = np.where(dist > 1e-10, dist, 1.0)
    r_hat = dr / safe
    valid = dist[..., 0] > DX * 0.5

    # Complex dot product P . r_hat
    P_dot_r = np.sum(P * r_hat, axis=-1)      # complex scalar per point
    A_L_sq = np.abs(P_dot_r) ** 2              # |P.r_hat|^2
    A_T_sq = np.sum(np.abs(P) ** 2, axis=-1) - A_L_sq   # |P|^2 - A_L^2
    A_T_sq = np.maximum(A_T_sq, 0.0)           # clamp float noise

    return A_L_sq, A_T_sq, valid


# ============================================================
# MAIN
# ============================================================

def main():
    print("=" * 65)
    print("  PHASE 1c  STEP 1 :  3D Vector Base Wave Validation")
    print("=" * 65)
    print(f"  Grid         : {N_GRID}^3 = {N_GRID**3:,} voxels")
    print(f"  Extent       : {GRID_LAM} lam  ({L:.1f} am)")
    print(f"  Resolution   : dx = {DX:.2f} am  ({LAM0/DX:.1f} vox/lam)")
    print(f"  Sources      : {N_SOURCES}  (Fibonacci sphere)")
    print(f"  Constants    : A0={A0}  lam0={LAM0}  f0={F0}  rho0={RHO0}")
    print()

    # ---- Build grid ----
    coords = build_grid(N_GRID, DX)

    # ---- Isotropic sources ----
    k_dirs = fibonacci_sphere(N_SOURCES)
    rng = np.random.default_rng(RNG_SEED)
    phases = rng.uniform(0, 2 * np.pi, N_SOURCES)

    dir_bias = np.abs(k_dirs.mean(axis=0)).max()
    print(f"  Source isotropy: max|<k_hat>| = {dir_bias:.5f}  (->0 for perfect)")

    # ---- Phasor field ----
    t0 = time.time()
    P = compute_phasor_field(coords, k_dirs, phases, A0, K0)
    dt = time.time() - t0
    print(f"  Phasor field computed in {dt:.1f}s")

    # ================================================================
    # TEST 1:  Energy uniformity
    # ================================================================
    E = energy_from_phasor(P)
    E_mean = E.mean()
    E_std = E.std()
    CV = E_std / E_mean
    E_theory = RHO0 * V_VOXEL * F0 ** 2 * A0 ** 2 / 2.0

    print(f"\n{'─' * 65}")
    print(f"  TEST 1: ENERGY DENSITY (expect uniform)")
    print(f"{'─' * 65}")
    print(f"  Theory (isotropic mean) : {E_theory:.6e}")
    print(f"  Computed mean           : {E_mean:.6e}  (ratio {E_mean / E_theory:.4f})")
    print(f"  Std / Min / Max         : {E_std:.4e}  /  {E.min():.4e}  /  {E.max():.4e}")
    CV_chi2 = 1.0 / np.sqrt(3)  # chi-squared(6) prediction for 3-component complex field
    print(f"  CV (std/mean)           : {CV:.4f}   (theory: {CV_chi2:.4f} = 1/sqrt(3), chi2(6))")

    # ================================================================
    # TEST 2:  Force = speckle noise only (no large-scale gradient)
    # ================================================================
    F = force_from_energy(E, DX)
    # Skip 2 boundary layers (gradient artifacts)
    interior = np.linalg.norm(F[2:-2, 2:-2, 2:-2], axis=-1)

    # Expected noise: |F| ~ CV * E_mean * k0 (gradient of speckle pattern)
    F_noise_theory = CV * E_mean * K0

    print(f"\n{'─' * 65}")
    print(f"  TEST 2: FORCE (speckle noise only, no large-scale gradient)")
    print(f"{'─' * 65}")
    print(f"  |F| mean                : {interior.mean():.6e}")
    print(f"  |F| max                 : {interior.max():.6e}")
    print(f"  Speckle noise estimate  : {F_noise_theory:.6e}  (CV * E * k0)")
    print(f"  |F|_mean / noise_est    : {interior.mean() / F_noise_theory:.4f}")

    # ================================================================
    # TEST 3:  L/T decomposition  (expect 1/3 L, 2/3 T for isotropic)
    # ================================================================
    ref = np.zeros(3)
    A_L_sq, A_T_sq, valid = lt_decomposition(P, ref, coords)
    amp_sq_total = np.sum(np.abs(P) ** 2, axis=-1)

    frac_L = (A_L_sq[valid] / amp_sq_total[valid]).mean()
    frac_T = (A_T_sq[valid] / amp_sq_total[valid]).mean()

    print(f"\n{'─' * 65}")
    print(f"  TEST 3: L/T DECOMPOSITION (ref = grid center)")
    print(f"{'─' * 65}")
    print(f"  E_L / E_total           : {frac_L:.4f}   (theory 0.3333)")
    print(f"  E_T / E_total           : {frac_T:.4f}   (theory 0.6667)")

    # ================================================================
    # TEST 4:  Energy conservation  E_L + E_T = E_total
    # ================================================================
    E_L = RHO0 * V_VOXEL * F0 ** 2 * A_L_sq / 2.0
    E_T = RHO0 * V_VOXEL * F0 ** 2 * A_T_sq / 2.0
    decomp_err = np.abs(E_L + E_T - E).max() / E_mean

    print(f"\n{'─' * 65}")
    print(f"  TEST 4: ENERGY CONSERVATION (E_L + E_T = E)")
    print(f"{'─' * 65}")
    print(f"  max |E_L+E_T - E| / E  : {decomp_err:.2e}")

    # ================================================================
    # TEST 5:  L/T at multiple reference points (isotropy check)
    # ================================================================
    print(f"\n{'─' * 65}")
    print(f"  TEST 5: L/T AT MULTIPLE REFERENCE POINTS")
    print(f"{'─' * 65}")

    offsets = [
        ("center",     np.array([0.0, 0.0, 0.0])),
        ("+2 lam x",   np.array([2 * LAM0, 0.0, 0.0])),
        ("-2 lam y",   np.array([0.0, -2 * LAM0, 0.0])),
        ("diagonal",   np.array([LAM0, LAM0, LAM0])),
    ]
    for label, rp in offsets:
        al, at, v = lt_decomposition(P, rp, coords)
        aq = amp_sq_total
        fl = (al[v] / aq[v]).mean()
        ft = (at[v] / aq[v]).mean()
        print(f"  {label:12s}  E_L/E={fl:.4f}  E_T/E={ft:.4f}")

    # ================================================================
    # CONVERGENCE STUDY:  CV vs N_sources (32^3 grid for speed)
    # ================================================================
    print(f"\n{'─' * 65}")
    print(f"  CONVERGENCE: CV vs N_sources (32^3 grid)")
    print(f"{'─' * 65}")

    dx_s = L / 32
    coords_s = build_grid(32, dx_s)
    chi2_cv = 1.0 / np.sqrt(3)
    print(f"  Theory: CV = 1/sqrt(3) = {chi2_cv:.4f} (independent of N)")
    print()
    for n_src in [50, 100, 200, 500, 1000]:
        kd = fibonacci_sphere(n_src)
        ph = rng.uniform(0, 2 * np.pi, n_src)
        P_s = compute_phasor_field(coords_s, kd, ph, A0, K0)
        amp_sq_s = np.sum(np.abs(P_s) ** 2, axis=-1)
        cv_s = amp_sq_s.std() / amp_sq_s.mean()
        print(f"  N={n_src:5d}  CV={cv_s:.4f}  theory={chi2_cv:.4f}  err={abs(cv_s - chi2_cv):.4f}")

    # ================================================================
    # VALIDATION SUMMARY
    # ================================================================
    print(f"\n{'=' * 65}")
    print(f"  VALIDATION SUMMARY")
    print(f"{'=' * 65}")

    checks = [
        ("Mean energy matches theory (+/- 5%)",
         abs(E_mean / E_theory - 1.0) < 0.05,
         f"ratio = {E_mean / E_theory:.4f}"),
        ("CV matches chi2(6) = 1/sqrt(3) (+/- 0.05)",
         abs(CV - 1.0 / np.sqrt(3)) < 0.05,
         f"CV = {CV:.4f}, theory = {1/np.sqrt(3):.4f}"),
        ("Force = speckle noise (within 2x estimate)",
         interior.mean() / F_noise_theory < 2.0,
         f"|F|/noise = {interior.mean() / F_noise_theory:.4f}"),
        ("E_L/E = 1/3 +/- 0.05",
         abs(frac_L - 1.0 / 3) < 0.05,
         f"E_L/E = {frac_L:.4f}"),
        ("E_T/E = 2/3 +/- 0.05",
         abs(frac_T - 2.0 / 3) < 0.05,
         f"E_T/E = {frac_T:.4f}"),
        ("E_L + E_T = E (err < 1e-10)",
         decomp_err < 1e-10,
         f"err = {decomp_err:.2e}"),
    ]

    all_pass = True
    for label, passed, detail in checks:
        tag = "PASS" if passed else "FAIL"
        print(f"  [{tag}]  {label}  --  {detail}")
        if not passed:
            all_pass = False

    print()
    if all_pass:
        print("  >>> ALL CHECKS PASSED -- base wave validated for Phase 1c")
    else:
        print("  >>> SOME CHECKS FAILED -- investigate above")
    print(f"{'=' * 65}")

    return all_pass


if __name__ == "__main__":
    main()
