#!/usr/bin/env python3
"""
Phase 1d: Variable λ(r) for Coulomb Force
==========================================
Tests whether variable wavelength near WC cores produces charge-dependent
Coulomb force that constant-λ models cannot.

Key hypothesis:
  K=1 (neutrino): particle radius = 1λ, minimal λ variation → sinc persists → neutral
  K=10 (electron): particle radius = 100λ, 10 shells → λ variation breaks sinc → Coulomb

Physics:
  - λ(r) from Yee & Hauger shells: r_wavelength(n) = 2(K-n)λ₀
  - WKB phase: φ(r) = ∫k(r')dr' replaces simple kr
  - Energy: E = ρ·V·(c·A/λ(r))² — the ∇λ term creates force from wavelength gradients
  - c is constant — only λ(r) varies (not wave speed)

1D analysis along axis connecting two WCs (high resolution, fast).

Usage: python step1d_variable_lambda.py
"""

import numpy as np
from scipy.integrate import cumulative_trapezoid
import time

# ============================================================
# PHYSICAL CONSTANTS (sim units: am, rs, qg)
# ============================================================
A0   = 0.92       # Fundamental amplitude [am]
LAM0 = 28.5       # Fundamental wavelength [am]
F0   = 0.0105     # Fundamental frequency [rHz]
RHO0 = 38.6       # Medium density [qg/am³]
C    = 0.3        # Wave speed [am/rs] — CONSTANT everywhere
K0   = 2 * np.pi / LAM0  # Fundamental wavenumber [1/am]

# ============================================================
# SIMULATION PARAMETERS
# ============================================================
N_POINTS = 10000   # 1D grid resolution

# K values to test
K_VALUES = [1, 10]

# Separations to test (in multiples of particle radius K²λ)
# Near-field: inside particle radius, Far-field: beyond particle radius
SEP_FACTORS = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]

# Phase configurations: (phase1, phase2, label)
PHASE_CONFIGS = [
    (0.0,    0.0,    "same phase (0,0)"),
    (np.pi,  np.pi,  "same phase (π,π)"),
    (0.0,    np.pi,  "opposite phase (0,π)"),
    (np.pi,  0.0,    "opposite phase (π,0)"),
]


# ============================================================
# λ(r) PROFILE — Yee & Hauger Shell Model
# ============================================================

def lambda_profile(r, K):
    """
    Local wavelength λ(r) from Yee & Hauger shell model.

    Regions:
      r < K·λ₀  (core):          λ = K·λ₀ (core wavelength, proportional to K)
      K·λ₀ ≤ r ≤ K²·λ₀ (shells): λ = 2(K-x(r))·λ₀, decreasing outward
      r > K²·λ₀ (beyond):        λ = λ₀ (fundamental)

    Shell inversion: r_x = (K + 2Kx - x² - x)·λ₀ → solve for x(r)
    Local shell wavelength: λ(x) = 2(K-x)·λ₀

    Args:
        r: distance from WC [am] (scalar or array)
        K: wave center count (1=neutrino, 10=electron)

    Returns:
        local wavelength [am]
    """
    r = np.asarray(r, dtype=float)
    r_core = K * LAM0
    r_particle = K**2 * LAM0

    # Default: fundamental wavelength
    lam = np.full_like(r, LAM0)

    # Core region: wavelength proportional to K (larger core = longer λ)
    core_mask = r < r_core
    lam[core_mask] = K * LAM0

    # Shell region: solve quadratic for continuous shell number x
    # r_x = (K + (2K-1)x - x²)·λ₀  →  x² - (2K-1)x + (r/λ₀ - K) = 0
    shell_mask = (r >= r_core) & (r <= r_particle)
    if shell_mask.any():
        r_s = r[shell_mask]
        a = 1.0
        b = -(2 * K - 1)
        c_coeff = r_s / LAM0 - K
        discriminant = b**2 - 4 * a * c_coeff
        discriminant = np.maximum(discriminant, 0.0)  # numerical safety
        x = (-b - np.sqrt(discriminant)) / (2 * a)
        x = np.clip(x, 0, K)
        lam[shell_mask] = 2 * (K - x) * LAM0
        # Clamp minimum λ to avoid division by zero
        lam[shell_mask] = np.maximum(lam[shell_mask], LAM0 * 0.1)

    return lam


def wkb_phase(r, K):
    """
    WKB phase integral: φ(r) = ∫₀ʳ k(r')dr' where k(r) = 2π/λ(r).

    Computed numerically via cumulative trapezoid integration on a fine radial grid.

    Args:
        r: distance(s) from WC [am] (array)
        K: wave center count

    Returns:
        accumulated phase [radians] at each distance
    """
    # Fine radial grid for integration
    r_max = np.max(r) + LAM0
    n_int = max(5000, int(r_max / LAM0 * 100))  # ~100 points per λ
    r_fine = np.linspace(0, r_max, n_int)

    # k(r) on fine grid
    lam_fine = lambda_profile(r_fine, K)
    k_fine = 2 * np.pi / lam_fine

    # Cumulative integral
    phi_fine = np.zeros_like(r_fine)
    phi_fine[1:] = cumulative_trapezoid(k_fine, r_fine)

    # Interpolate to requested r values
    return np.interp(r, r_fine, phi_fine)


# ============================================================
# WC OUT-WAVE WITH WKB PHASE
# ============================================================

def compute_wc_outwave_1d(x, wc_pos, wc_phase, K, amplitude):
    """
    Compute WC out-wave phasor along 1D axis using WKB phase.

    Out-wave: P(x) = A · sinc_wkb(r) · exp(+i·(φ_wkb(r) + source_offset))

    where sinc_wkb = sin(φ_wkb) / φ_wkb uses WKB phase instead of kr.
    This is a modified sinc with non-uniform node spacing from variable λ(r).

    Args:
        x: 1D grid positions [am]
        wc_pos: WC position [am]
        wc_phase: source_offset [radians]
        K: wave center count
        amplitude: out-wave amplitude [am]
    """
    r = np.abs(x - wc_pos)
    r_safe = np.maximum(r, LAM0 * 0.01)  # avoid singularity at WC

    # WKB phase: ∫k(r')dr' — non-uniform node spacing
    phi = wkb_phase(r_safe, K)

    # Modified sinc envelope: sin(φ)/φ instead of sin(kr)/kr
    # At r→0: φ→0, sin(φ)/φ→1
    phi_safe = np.where(phi > 1e-10, phi, 1e-10)
    envelope = amplitude * np.sin(phi_safe) / phi_safe
    envelope = np.where(phi < 1e-10, amplitude, envelope)

    # Out-wave phasor: exp(+i·(φ_wkb + source_offset))
    phasor = envelope * np.exp(1j * (phi + wc_phase))

    return phasor


# ============================================================
# ENERGY EQUATIONS (constant-λ vs variable-λ)
# ============================================================

def energy_constant_lambda(P, dx):
    """Standard energy: E = ρ·dx·(f₀·A_rms)² — constant λ, can't see λ variation."""
    A_rms = np.abs(P) / np.sqrt(2)
    return RHO0 * dx * (F0 * A_rms)**2


def energy_variable_lambda(P, dx, lam_local):
    """Variable-λ energy: E = ρ·dx·(c·A_rms/λ(r))² — the ∇λ term is active."""
    A_rms = np.abs(P) / np.sqrt(2)
    return RHO0 * dx * (C * A_rms / lam_local)**2


def force_from_energy_1d(E, x):
    """F = -dE/dx via numpy gradient."""
    return -np.gradient(E, x)


# ============================================================
# MAIN
# ============================================================

def main():
    print("=" * 70)
    print("  PHASE 1d :  Variable λ(r) for Coulomb Force")
    print("=" * 70)
    print(f"  K values     : {K_VALUES}")
    print(f"  Grid points  : {N_POINTS}")
    print(f"  Sep factors  : {SEP_FACTORS} × K²λ (particle radius)")
    print()

    for K in K_VALUES:
        r_core = K * LAM0
        r_particle = K**2 * LAM0

        print(f"\n{'#' * 70}")
        print(f"  K = {K}  ({'neutrino' if K == 1 else 'electron'})")
        print(f"  Core radius    : {r_core:.1f} am  ({K} λ)")
        print(f"  Particle radius: {r_particle:.1f} am  ({K**2} λ)")
        print(f"{'#' * 70}")

        # ---- λ(r) profile diagnostic ----
        r_diag = np.linspace(0, r_particle * 1.5, 200)
        lam_diag = lambda_profile(r_diag, K)
        print(f"\n  λ(r) profile (K={K}):")
        print(f"  {'r/λ₀':>8s}  {'λ(r)/λ₀':>8s}  {'k(r)/k₀':>8s}  {'region':>12s}")
        for i in range(0, len(r_diag), 20):
            rr = r_diag[i]
            ll = lam_diag[i]
            region = "core" if rr < r_core else ("shells" if rr < r_particle else "far-field")
            print(f"  {rr/LAM0:8.2f}  {ll/LAM0:8.3f}  {LAM0/ll:8.3f}  {region:>12s}")

        # ---- Phase accumulation diagnostic ----
        phi_diag = wkb_phase(r_diag, K)
        phi_const = K0 * r_diag  # constant-λ reference
        print(f"\n  WKB phase vs constant-k phase:")
        print(f"  {'r/λ₀':>8s}  {'φ_wkb':>10s}  {'φ_const':>10s}  {'Δφ':>10s}  {'Δφ/π':>8s}")
        for i in range(0, len(r_diag), 20):
            rr = r_diag[i]
            pw = phi_diag[i]
            pc = phi_const[i]
            dp = pw - pc
            print(f"  {rr/LAM0:8.2f}  {pw:10.3f}  {pc:10.3f}  {dp:10.3f}  {dp/np.pi:8.3f}")

        # ---- Two-WC force test ----
        print(f"\n  {'─' * 66}")
        print(f"  TWO-WC FORCE TEST (K={K})")
        print(f"  {'─' * 66}")

        # Results table header
        print(f"\n  {'sep':>6s}  {'sep/r_p':>7s}  {'config':>25s}  "
              f"{'F_var':>10s}  {'F_const':>10s}  {'dir_v':>5s}  {'dir_c':>5s}  {'match':>5s}")
        print(f"  {'─'*6}  {'─'*7}  {'─'*25}  {'─'*10}  {'─'*10}  {'─'*5}  {'─'*5}  {'─'*5}")

        # Track results
        results_var = {}
        results_const = {}

        for sep_factor in SEP_FACTORS:
            # Separation in physical units
            d = sep_factor * r_particle + r_particle  # distance from center to center
            d_lam = d / LAM0

            # 1D grid: extends well beyond both particles
            margin = 2 * r_particle
            x = np.linspace(-margin, d + margin, N_POINTS)
            dx = x[1] - x[0]

            # WC positions
            wc1_pos = 0.0
            wc2_pos = d

            # Distance from each WC
            r1 = np.abs(x - wc1_pos)
            r2 = np.abs(x - wc2_pos)

            # λ(r) fields from each WC — use NEAREST WC's λ
            # The local medium is modified by the closest WC's standing wave structure
            # Near WC1: use WC1's λ(r1). Near WC2: use WC2's λ(r2).
            # At the midpoint: use whichever WC is closer.
            lam1 = lambda_profile(r1, K)
            lam2 = lambda_profile(r2, K)
            lam_local = np.where(r1 <= r2, lam1, lam2)

            for phase1, phase2, label in PHASE_CONFIGS:
                # WC out-waves with WKB phase
                P_wc1 = compute_wc_outwave_1d(x, wc1_pos, phase1, K, A0)
                P_wc2 = compute_wc_outwave_1d(x, wc2_pos, phase2, K, A0)

                # Base wave: standing wave (simplified 1D)
                P_base = A0 * np.cos(K0 * x)

                # Total field
                P_total = P_base + P_wc1 + P_wc2

                # Energy: variable-λ vs constant-λ
                E_var = energy_variable_lambda(P_total, dx, lam_local)
                E_const = energy_constant_lambda(P_total, dx)

                # Force at WC positions (midpoint of WC1)
                # Sample force near WC1: average over a few points around WC1
                idx_wc1 = np.argmin(np.abs(x - wc1_pos))
                idx_wc2 = np.argmin(np.abs(x - wc2_pos))

                F_var = force_from_energy_1d(E_var, x)
                F_const = force_from_energy_1d(E_const, x)

                # Force on WC1 (x-component, positive = toward WC2 = attraction)
                # Average over small region around WC1 to reduce noise
                hw = max(3, int(r_core / dx))  # half-window = core radius
                sl1 = slice(max(idx_wc1 - hw, 0), min(idx_wc1 + hw + 1, N_POINTS))
                F_var_wc1 = F_var[sl1].mean()
                F_const_wc1 = F_const[sl1].mean()

                # Direction: positive = toward WC2 (attraction), negative = away (repulsion)
                dir_v = "ATT" if F_var_wc1 > 0 else "REP"
                dir_c = "ATT" if F_const_wc1 > 0 else "REP"
                match = "YES" if dir_v == dir_c else " NO"

                print(f"  {d_lam:6.1f}  {sep_factor+1:7.2f}  {label:>25s}  "
                      f"{F_var_wc1:10.3e}  {F_const_wc1:10.3e}  {dir_v:>5s}  {dir_c:>5s}  {match:>5s}")

                # Track for summary
                key = label
                if key not in results_var:
                    results_var[key] = []
                    results_const[key] = []
                results_var[key].append(dir_v)
                results_const[key].append(dir_c)

        # ---- Direction consistency summary ----
        print(f"\n  {'─' * 66}")
        print(f"  DIRECTION CONSISTENCY (K={K})")
        print(f"  {'─' * 66}")
        for label in results_var:
            dirs_v = results_var[label]
            dirs_c = results_const[label]
            att_v = dirs_v.count("ATT")
            rep_v = dirs_v.count("REP")
            att_c = dirs_c.count("ATT")
            rep_c = dirs_c.count("REP")
            n = len(dirs_v)
            cons_v = "CONSISTENT" if (att_v == n or rep_v == n) else "MIXED"
            cons_c = "CONSISTENT" if (att_c == n or rep_c == n) else "MIXED"
            dom_v = "ATT" if att_v > rep_v else "REP"
            dom_c = "ATT" if att_c > rep_c else "REP"

            print(f"  {label:>25s}  var-λ: {att_v}/{n}ATT {rep_v}/{n}REP {cons_v:>10s} ({dom_v})"
                  f"  const-λ: {att_c}/{n}ATT {rep_c}/{n}REP {cons_c:>10s} ({dom_c})")

    # ---- Final verdict ----
    print(f"\n{'=' * 70}")
    print(f"  PHASE 1d VERDICT")
    print(f"{'=' * 70}")
    print(f"  Does variable λ(r) produce consistent charge-dependent force?")
    print(f"  K=1 (neutrino): expect MIXED (sinc/lock-in, neutral particle)")
    print(f"  K=10 (electron): expect CONSISTENT (Coulomb, charged particle)")
    print(f"  If var-λ is CONSISTENT where const-λ is MIXED → λ(r) IS the mechanism")
    print(f"{'=' * 70}")


if __name__ == "__main__":
    main()
