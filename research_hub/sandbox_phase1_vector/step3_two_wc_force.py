#!/usr/bin/env python3
"""
Phase 1c Step 3: Two-WC Force Test with Spin
=============================================
The critical test: does L->T spin conversion produce charge-dependent
force direction between two wave centers?

Setup:
  - Base wave: N=200 isotropic plane in-waves (from Step 1)
  - Two WCs on x-axis, separated by variable distance d
  - Each WC emits spherical out-wave with L + T components (from Step 2)
  - Four spin configurations: CW-CW, CCW-CCW, CW-CCW, CCW-CW

Tests:
  1. Does force direction depend on spin combination? (charge emergence)
  2. Is same-spin repulsive and opposite-spin attractive? (or vice versa)
  3. Is the direction consistent across separations? (no sinc flips)
  4. Newton's 3rd law: F_on_WC1 = -F_on_WC2?
  5. Force magnitude scaling with distance

Physics:
  Same spin (CW-CW): T patterns from both WCs rotate same direction
    -> constructive T interference between them -> energy hill -> repulsion?
  Opposite spin (CW-CCW): T patterns rotate opposite directions
    -> destructive T interference between them -> energy well -> attraction?

Usage: python step3_two_wc_force.py
"""

import numpy as np
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
# GRID PARAMETERS
# ============================================================
N_GRID   = 64
GRID_LAM = 8.0
L        = GRID_LAM * LAM0
DX       = L / N_GRID
V_VOXEL  = DX ** 3

# ============================================================
# BASE WAVE
# ============================================================
N_SOURCES = 200
RNG_SEED  = 42

# ============================================================
# WAVE CENTER PARAMETERS
# ============================================================
WC_AMPLITUDE  = A0
SPIN_AXIS     = np.array([0.0, 0.0, 1.0])
ALPHA         = 7.2974e-3   # Fine structure constant ~ 1/137

# Separations to sweep (in wavelengths)
SEPARATIONS_LAM = [2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0]

# Spin configurations: (spin1, spin2, phase1, phase2, label)
CONFIGS = [
    (+1, +1, np.pi, np.pi, "CW-CW   (same spin, same phase)"),
    (-1, -1, np.pi, np.pi, "CCW-CCW  (same spin, same phase)"),
    (+1, -1, np.pi, 0.0,   "CW-CCW   (opp spin, opp phase)"),
    (-1, +1, 0.0,   np.pi, "CCW-CW   (opp spin, opp phase)"),
]

# Eta values to test
ETA_VALUES = [ALPHA, 0.1, 0.5]


# ============================================================
# SHARED FUNCTIONS (from Steps 1-2)
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
    n_src = len(k_dirs)
    amp = a0 / np.sqrt(n_src)
    P = np.zeros(coords.shape, dtype=np.complex128)
    for n in range(n_src):
        kr = k0 * np.einsum('...j,j', coords, k_dirs[n])
        z = amp * np.exp(-1j * (kr + phases[n]))
        P += z[..., np.newaxis] * k_dirs[n]
    return P


def compute_wc_outwave(coords, wc_pos, wc_phase, wc_amp, k0,
                       eta, spin_sign, spin_axis):
    """Spherical out-wave with L + T components from WC."""
    r_vec = coords - wc_pos
    r_mag = np.linalg.norm(r_vec, axis=-1)
    r_safe = np.maximum(r_mag, 1e-10)
    r_hat = r_vec / r_safe[..., np.newaxis]

    # sinc envelope: sin(kr)/(kr)
    kr = k0 * r_mag
    kr_safe = np.where(kr > 1e-10, kr, 1e-10)
    envelope = wc_amp * np.sin(kr_safe) / kr_safe
    envelope = np.where(kr < 1e-10, wc_amp, envelope)

    # Out-wave phasor: exp(+i*(kr + source_offset))
    phase = kr + wc_phase
    z = envelope * np.exp(1j * phase)

    # L component (radial)
    P_L = np.sqrt(1.0 - eta) * z[..., np.newaxis] * r_hat

    # T component (azimuthal from spin)
    t_vec = np.cross(spin_axis, r_hat)
    P_T = np.sqrt(eta) * spin_sign * z[..., np.newaxis] * t_vec

    return P_L + P_T


def energy_from_phasor(P):
    return RHO0 * V_VOXEL * F0 ** 2 * np.sum(np.abs(P) ** 2, axis=-1) / 2.0


def force_from_energy(E, dx):
    F = np.zeros(E.shape + (3,))
    F[1:-1, :, :, 0] = -(E[2:, :, :] - E[:-2, :, :]) / (2.0 * dx)
    F[:, 1:-1, :, 1] = -(E[:, 2:, :] - E[:, :-2, :]) / (2.0 * dx)
    F[:, :, 1:-1, 2] = -(E[:, :, 2:] - E[:, :, :-2]) / (2.0 * dx)
    return F


def force_at_point(F_field, coords, point):
    """Sample force at the grid point closest to the given position."""
    # Find closest grid point
    offset = coords[0, 0, 0]  # grid origin
    idx = np.round((point - offset) / DX).astype(int)
    # Clamp to valid range
    idx = np.clip(idx, 2, N_GRID - 3)
    return F_field[idx[0], idx[1], idx[2]]


# ============================================================
# MAIN
# ============================================================

def main():
    print("=" * 70)
    print("  PHASE 1c  STEP 3 :  Two-WC Force Test with Spin")
    print("=" * 70)
    print(f"  Grid         : {N_GRID}^3, {GRID_LAM} lam")
    print(f"  Separations  : {SEPARATIONS_LAM} lam")
    print(f"  Alpha (1/137): {ALPHA:.4e}")
    print(f"  Eta values   : {[f'{e:.4f}' for e in ETA_VALUES]}")
    print()

    # ---- Build grid ----
    coords = build_grid(N_GRID, DX)
    grid_origin = coords[0, 0, 0]

    # ---- Base wave ----
    k_dirs = fibonacci_sphere(N_SOURCES)
    rng = np.random.default_rng(RNG_SEED)
    phases = rng.uniform(0, 2 * np.pi, N_SOURCES)

    t0 = time.time()
    P_base = compute_base_phasor(coords, k_dirs, phases, A0, K0)
    print(f"  Base wave computed in {time.time() - t0:.1f}s")

    # ---- Sweep eta values ----
    for eta in ETA_VALUES:
        eta_label = f"alpha" if eta == ALPHA else f"{eta:.2f}"
        print(f"\n{'=' * 70}")
        print(f"  ETA = {eta_label} ({eta:.4e})")
        print(f"  L out = {np.sqrt(1-eta)*100:.1f}%,  T out = {np.sqrt(eta)*100:.1f}%")
        print(f"{'=' * 70}")

        # ---- Results table header ----
        print(f"\n  {'sep':>5s}  {'config':>38s}  {'Fx_WC1':>10s}  {'Fx_WC2':>10s}  "
              f"{'dir':>5s}  {'N3rd':>6s}  {'|F|':>10s}")
        print(f"  {'─' * 5}  {'─' * 38}  {'─' * 10}  {'─' * 10}  {'─' * 5}  {'─' * 6}  {'─' * 10}")

        # Track results for summary
        results = {label: [] for _, _, _, _, label in CONFIGS}

        for sep_lam in SEPARATIONS_LAM:
            sep_am = sep_lam * LAM0  # Physical separation

            # WC positions symmetric about origin on x-axis
            wc1_pos = np.array([-sep_am / 2, 0.0, 0.0])
            wc2_pos = np.array([+sep_am / 2, 0.0, 0.0])

            for spin1, spin2, phase1, phase2, label in CONFIGS:
                # Compute WC out-waves
                P_wc1 = compute_wc_outwave(coords, wc1_pos, phase1, WC_AMPLITUDE,
                                           K0, eta, spin1, SPIN_AXIS)
                P_wc2 = compute_wc_outwave(coords, wc2_pos, phase2, WC_AMPLITUDE,
                                           K0, eta, spin2, SPIN_AXIS)

                # Total field
                P_total = P_base + P_wc1 + P_wc2

                # Energy and force
                E_total = energy_from_phasor(P_total)
                F = force_from_energy(E_total, DX)

                # Sample force at each WC position
                F_wc1 = force_at_point(F, coords, wc1_pos)
                F_wc2 = force_at_point(F, coords, wc2_pos)

                # Force along x-axis (connecting the two WCs)
                Fx1 = F_wc1[0]   # x-component at WC1
                Fx2 = F_wc2[0]   # x-component at WC2

                # Direction: WC1 is at -x, WC2 at +x
                # Attraction: WC1 pushed toward +x (Fx1 > 0), WC2 toward -x (Fx2 < 0)
                # Repulsion: WC1 pushed toward -x (Fx1 < 0), WC2 toward +x (Fx2 > 0)
                if Fx1 > 0 and Fx2 < 0:
                    direction = "ATT"
                elif Fx1 < 0 and Fx2 > 0:
                    direction = "REP"
                else:
                    direction = "???"

                # Newton's 3rd law check: F1x should be ~ -F2x
                if abs(Fx1) > 1e-15 and abs(Fx2) > 1e-15:
                    n3rd_ratio = abs(Fx1 / Fx2)
                    n3rd = f"{n3rd_ratio:.3f}"
                else:
                    n3rd = "N/A"

                F_mag = (abs(Fx1) + abs(Fx2)) / 2

                print(f"  {sep_lam:5.1f}  {label:>38s}  {Fx1:10.3e}  {Fx2:10.3e}  "
                      f"{direction:>5s}  {n3rd:>6s}  {F_mag:10.3e}")

                results[label].append((sep_lam, direction, Fx1, Fx2, F_mag))

        # ---- Summary for this eta ----
        print(f"\n  {'─' * 70}")
        print(f"  DIRECTION SUMMARY (eta = {eta_label})")
        print(f"  {'─' * 70}")

        for label in results:
            directions = [r[1] for r in results[label]]
            att_count = directions.count("ATT")
            rep_count = directions.count("REP")
            unk_count = directions.count("???")
            total = len(directions)
            consistency = "CONSISTENT" if (att_count == total or rep_count == total) else "MIXED"
            dominant = "ATT" if att_count > rep_count else "REP" if rep_count > att_count else "???"

            print(f"  {label:>38s}  ATT={att_count}/{total}  REP={rep_count}/{total}  "
                  f"???={unk_count}/{total}  {consistency:>10s}  {dominant}")

    # ---- Final validation ----
    print(f"\n{'=' * 70}")
    print(f"  CHARGE EMERGENCE TEST")
    print(f"{'=' * 70}")
    print(f"  Question: does spin sign determine force direction?")
    print(f"  If same-spin ALWAYS repels and opposite-spin ALWAYS attracts")
    print(f"  (or vice versa), then charge emerges from spin.")
    print(f"  If mixed — spin alone doesn't determine charge.")
    print(f"{'=' * 70}")


if __name__ == "__main__":
    main()
