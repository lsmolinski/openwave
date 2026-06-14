"""
M5.8.2u — CLOCK AS ENERGY MINIMUM (Duda 2026-06-11).

Duda's statement (Sulich ZBW thread, models-of-particles, 2026-06-11):
  "Clock propulsion ... preferring oscillations ... these oscillations need to
   allow to reduce energy, which minimum are the preferred frequencies"
  "You can only increase its energy - above 511keV for energy minimizing:
   de Broglie frequency."
i.e. the de Broglie clock is the ENERGY-MINIMIZING state; a clock-stopped
defect has HIGHER energy. Toymodel reference: arXiv:2501.04036 (energy vs
oscillation frequency has a minimum at the de Broglie value).

This script tests the DIRECTION of that claim at the seed level, reusing the
validated rest-energy primitives from m5_8_2q_delta_scaling.py. The clock is
dressed by the time-axis boost b_star (the negative-GEM "clock-fuel" block, the
2026-06-09 boost decomposition); b_star = 0 is the clock-stopped seed. We sweep
b_star finely and locate the minimum of H_static(b_star).

PASS (direction reproduced): there is a b_star > 0 with H_static < H_static(0),
i.e. activating the clock-fuel sector LOWERS the rest energy.

Caveat: this is the SEED-LEVEL static energy with the boost dressing as the
clock-activation proxy, not the full dynamical de-Broglie-frequency sweep; the
dynamical confirmation (a settled-clock vs frozen-seed energy comparison at a
converged timestep) is the named follow-up.

Headless, numpy-only (no Taichi evolution). Writes a small npz + a plot.
"""

import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

# Reuse the validated seed-level rest-energy primitives (no edits to that module).
from m5_8_2q_delta_scaling import H_of_b, sectors_of, B_STAR  # noqa: E402

DELTA, G_TIME = 0.3, 8.0          # the calibrated gate point (m5_8_2q)
OUT_NPZ = HERE / "data" / "_m5_8_2u_clock_energy_minimum.npz"
OUT_PNG = HERE / "plots" / "_m5_8_2u_clock_energy_minimum.png"


def main():
    print("=" * 78)
    print("M5.8.2u — CLOCK AS ENERGY MINIMUM (Duda 2026-06-11)")
    print(f"  seed gate: delta={DELTA}, g={G_TIME}; clock dressed by boost b_star")
    print("=" * 78)

    # Fine sweep across the small-boost regime where the clock-fuel acts, plus a
    # few larger points to show the curvature cost taking over (the minimum).
    b_grid = np.concatenate([
        np.linspace(0.0, 0.16, 33),   # fine, brackets the production B_STAR=0.13
        np.array([0.20, 0.25, 0.30]),  # coarse tail (curvature cost dominates)
    ])

    H = np.array([H_of_b(DELTA, G_TIME, float(b))[0] for b in b_grid])
    GEM = np.array([sectors_of(DELTA, G_TIME, float(b))[1] for b in b_grid])

    H0 = H[0]                              # clock stopped (b_star = 0)
    i_min = int(np.argmin(H))
    b_min, H_min = float(b_grid[i_min]), float(H[i_min])
    drop = H0 - H_min
    drop_pct = 100.0 * drop / H0

    # H at the production dressing B_STAR (interpolated onto the grid label)
    H_prod = float(H_of_b(DELTA, G_TIME, B_STAR)[0])

    print("\n  | b_star | H_static | GEM (clock-fuel) | dH vs b=0 |")
    print("  | --- | --- | --- | --- |")
    for b, h, g in zip(b_grid, H, GEM):
        mark = "  <- MIN" if abs(b - b_min) < 1e-9 else ""
        if b <= 0.16 + 1e-9 or mark:
            print(f"  | {b:6.3f} | {h:9.4f} | {g:+10.4f} | {h - H0:+8.4f} |{mark}")

    print("\n  RESULT:")
    print(f"    H_static(clock stopped, b=0) = {H0:.4f}")
    print(f"    H_static(minimum)            = {H_min:.4f}  at b_star = {b_min:.4f}")
    print(f"    energy reduction by activating the clock-fuel = {drop:.4f} ({drop_pct:.2f}%)")
    print(f"    H_static(production B_STAR={B_STAR}) = {H_prod:.4f} "
          f"({'past' if B_STAR > b_min else 'at/before'} the minimum)")
    verdict = "PASS (direction reproduced)" if H_min < H0 else "FAIL"
    print(f"    VERDICT: {verdict} — clock activation "
          f"{'LOWERS' if H_min < H0 else 'does not lower'} the rest energy")
    print("  CAVEAT: seed-level static energy, boost-dressing as the clock proxy;")
    print("  the dynamical de-Broglie-frequency sweep is the named follow-up.")
    print("=" * 78)

    np.savez(
        OUT_NPZ,
        b_grid=b_grid, H=H, GEM=GEM,
        H0=H0, b_min=b_min, H_min=H_min, drop=drop, drop_pct=drop_pct,
        H_prod=H_prod, B_STAR=B_STAR, delta=DELTA, g_time=G_TIME,
    )

    fig, ax = plt.subplots(figsize=(7, 4.5))
    m = b_grid <= 0.16
    ax.plot(b_grid[m], H[m], "o-", color="#1f77b4", label="H_static(b_star)")
    ax.axhline(H0, ls="--", color="gray", lw=1, label=f"clock stopped (b=0): {H0:.3f}")
    ax.plot([b_min], [H_min], "*", color="crimson", ms=16,
            label=f"minimum: {H_min:.3f} @ b={b_min:.3f}")
    ax.axvline(B_STAR, ls=":", color="green", lw=1, label=f"production B_STAR={B_STAR}")
    ax.set_xlabel("boost dressing b_star  (clock-fuel knob; b=0 = clock stopped)")
    ax.set_ylabel("H_static (rest energy, seed level)")
    ax.set_title("M5.8.2u clock as energy minimum: activating the clock lowers H")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=110)
    print(f"  plot -> {OUT_PNG.name}")
    print(f"  data -> {OUT_NPZ.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
