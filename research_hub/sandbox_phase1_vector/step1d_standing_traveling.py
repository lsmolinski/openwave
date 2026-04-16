#!/usr/bin/env python3
"""
Phase 1d: Standing vs Traveling Wave Force Decomposition
========================================================
LaFreniere's key insight: standing waves → strong force (lock-in),
traveling waves → electrostatic force (Coulomb).

The weighted partial standing wave:
  ψ = A · [w·sin(kr+ωt+φ) + sin(kr-ωt-φ)] / kr

Decomposes into phasor components:
  C_n = (w+1)·sin(kr)/kr   ← STANDING (dominates near WC, has sinc nodes)
  S_n = (w-1)·cos(kr)/kr   ← TRAVELING (dominates far, smooth 1/kr decay)

Test: compute force from traveling component ONLY for two WCs.
If traveling gives correct Coulomb (same repels, opposite attracts)
while standing gives lock-in (capture) — the decomposition works.

Usage: python step1d_standing_traveling.py
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
# GRID
# ============================================================
N_GRID = 512     # 2D
GRID_LAM = 16.0  # For K=1
L = GRID_LAM * LAM0
DX = L / N_GRID
V = DX**2 * DX

# Weight transition (same as M3 wave engine)
TRANSITION_LAM = 1.25  # wavelengths
WEIGHT_POWER = 8       # sharp rolloff

# Separations to test
SEPS_LAM = [2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0]

PHASE_CONFIGS = [
    (0.0,    0.0,    "same (0,0)"),
    (np.pi,  np.pi,  "same (π,π)"),
    (0.0,    np.pi,  "opp (0,π)"),
    (np.pi,  0.0,    "opp (π,0)"),
]


# ============================================================
# 2D GRID
# ============================================================

def build_grid_2d(n, dx):
    ax = (np.arange(n) - n/2.0 + 0.5) * dx
    gx, gy = np.meshgrid(ax, ax, indexing='ij')
    return gx, gy


# ============================================================
# WEIGHTED PARTIAL STANDING WAVE — DECOMPOSED
# ============================================================

def compute_wc_phasors_decomposed(gx, gy, wc_x, wc_y, source_offset, amplitude):
    """
    Decompose WC wave into IN-WAVE and OUT-WAVE phasors.

    From M3 wave engine:
      ψ = A · [w·sin(kr+ωt+φ) + sin(kr-ωt-φ)] / kr
           ↑ in-wave (inward)     ↑ out-wave (outward)

    Complex phasors (ψ = Re[Z·e^{-iωt}]):
      Z_in  =  w · A · i · exp(-i(kr+φ)) / kr    ← e^{-ikr} (inward)
      Z_out = -A · i · exp(+i(kr-φ)) / kr         ← e^{+ikr} (outward)

    The in-wave creates standing wave structure (lock-in).
    The out-wave carries energy outward (Coulomb / radiation pressure).

    Returns: Z_in, Z_out, Z_full (complex phasors on 2D grid)
    """
    rx = gx - wc_x
    ry = gy - wc_y
    r = np.sqrt(rx**2 + ry**2)
    kr = K0 * r
    kr_safe = np.maximum(kr, 1e-10)

    # Weight function (same as M3)
    transition_r = TRANSITION_LAM * LAM0
    weight = 1.0 / (1.0 + (r / transition_r) ** WEIGHT_POWER)

    # 1/kr envelope (sinc-like, clamped at center)
    inv_kr = np.where(kr < 1e-10, 1.0, 1.0 / kr_safe)

    # In-wave phasor: w · A · i · exp(-i(kr+φ)) / kr
    Z_in = weight * amplitude * 1j * np.exp(-1j * (kr + source_offset)) * inv_kr

    # Out-wave phasor: -A · i · exp(+i(kr-φ)) / kr
    Z_out = -amplitude * 1j * np.exp(1j * (kr - source_offset)) * inv_kr

    # Full = in + out
    Z_full = Z_in + Z_out

    return Z_in, Z_out, Z_full


# ============================================================
# ENERGY AND FORCE
# ============================================================

def energy_2d(P, V):
    return RHO0 * V * F0**2 * np.abs(P)**2 / 2.0


def force_x_2d(E, dx):
    Fx = np.zeros_like(E)
    Fx[1:-1, :] = -(E[2:, :] - E[:-2, :]) / (2*dx)
    return Fx


def sample_fx(Fx, gx, dx, px, py, core_r):
    n = gx.shape[0]
    origin = gx[0, 0]
    ix = int(round((px - origin) / dx))
    iy = int(round((py - origin) / dx))
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
    print("  PHASE 1d :  Standing vs Traveling Wave Force Decomposition")
    print("=" * 70)
    print(f"  Grid: {N_GRID}², {GRID_LAM}λ, {LAM0/DX:.1f} vox/λ")
    print(f"  Weight transition: {TRANSITION_LAM}λ, power={WEIGHT_POWER}")
    print(f"  Hypothesis: in-wave (standing) → strong force, out-wave (traveling) → Coulomb")
    print()

    gx, gy = build_grid_2d(N_GRID, DX)
    r_core = LAM0  # K=1 core radius

    # Three force components to track
    components = ["in-wave", "out-wave", "full"]

    print(f"  {'sep':>5s}  {'config':>14s}  "
          f"{'Fx_in':>10s}  {'Fx_out':>10s}  {'Fx_full':>10s}  "
          f"{'d_in':>5s}  {'d_out':>6s}  {'d_full':>6s}")
    print(f"  {'─'*5}  {'─'*14}  {'─'*10}  {'─'*10}  {'─'*10}  {'─'*5}  {'─'*6}  {'─'*6}")

    results = {c: {} for c in components}

    for sep_lam in SEPS_LAM:
        sep = sep_lam * LAM0
        wc1_x, wc1_y = -sep/2, 0.0
        wc2_x, wc2_y = +sep/2, 0.0

        for phase1, phase2, label in PHASE_CONFIGS:
            # Decomposed phasors for each WC
            S1, T1, F1 = compute_wc_phasors_decomposed(gx, gy, wc1_x, wc1_y, phase1, A0)
            S2, T2, F2 = compute_wc_phasors_decomposed(gx, gy, wc2_x, wc2_y, phase2, A0)

            forces = {}
            dirs = {}

            for comp_name, P1, P2 in [("in-wave", S1, S2),
                                       ("out-wave", T1, T2),
                                       ("full", F1, F2)]:
                # Interaction energy: E(1+2) - E(1) - E(2)
                E_int = energy_2d(P1 + P2, V) - energy_2d(P1, V) - energy_2d(P2, V)
                Fx = force_x_2d(E_int, DX)
                fx = sample_fx(Fx, gx, DX, wc1_x, wc1_y, r_core)
                forces[comp_name] = fx
                dirs[comp_name] = "ATT" if fx > 0 else "REP"

                if label not in results[comp_name]:
                    results[comp_name][label] = []
                results[comp_name][label].append(dirs[comp_name])

            print(f"  {sep_lam:5.1f}  {label:>14s}  "
                  f"{forces['in-wave']:10.3e}  {forces['out-wave']:10.3e}  {forces['full']:10.3e}  "
                  f"{dirs['in-wave']:>5s}  {dirs['out-wave']:>6s}  {dirs['full']:>6s}")

    # ---- Charge sensitivity per component ----
    for comp in components:
        print(f"\n  {'─' * 66}")
        print(f"  CHARGE SENSITIVITY — {comp.upper()} component")
        print(f"  {'─' * 66}")
        for i, sep_lam in enumerate(SEPS_LAM):
            d_same = results[comp]["same (0,0)"][i]
            d_opp = results[comp]["opp (0,π)"][i]
            charge = "CHARGE-DEP" if d_same != d_opp else "blind"
            coulomb = ""
            if d_same != d_opp:
                # Coulomb: same should REP, opposite should ATT
                if d_same == "REP" and d_opp == "ATT":
                    coulomb = "  ← COULOMB CORRECT"
                elif d_same == "ATT" and d_opp == "REP":
                    coulomb = "  ← INVERTED (strong force?)"
            print(f"    sep={sep_lam:5.1f}λ  same={d_same}  opp={d_opp}  {charge}{coulomb}")

    # ---- Consistency per component ----
    print(f"\n  {'=' * 66}")
    print(f"  DIRECTION CONSISTENCY BY COMPONENT")
    print(f"  {'=' * 66}")
    for comp in components:
        print(f"\n  {comp.upper()}:")
        for label in results[comp]:
            d = results[comp][label]
            n = len(d)
            a, r = d.count("ATT"), d.count("REP")
            c = "CONSISTENT" if (a == n or r == n) else "MIXED"
            print(f"    {label:>14s}  {a}/{n}A {r}/{n}R  {c}")

    print(f"\n{'=' * 70}")
    print(f"  KEY QUESTION: Does OUT-WAVE give correct Coulomb sign?")
    print(f"  same→REP + opp→ATT = Coulomb correct")
    print(f"  same→ATT + opp→REP = strong force / lock-in (inverted)")
    print(f"{'=' * 70}")


if __name__ == "__main__":
    main()
