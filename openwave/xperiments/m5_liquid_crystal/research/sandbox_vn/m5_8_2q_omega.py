"""
M5.8.2q (Phase B) — the CLOCK frequency ω₁(δ), the half of the Duda calibration
that needs the evolution. Reads δ, g from env (M58_DELTA / M58_G — the 2026-06-08
overrides) so seed_M builds the δ-seed; the evolution kernels are δ-agnostic, so
this reuses the VALIDATED constrained integrator (d2.run_beta) unchanged. No cache
clobbering — builds its seed in-process, prints the result, writes nothing.

GATE OUTCOME (2026-06-08): at δ=0.3 this gives ω≈0.60, NOT the N-1/N-3 record
ω₁≈1.10 — because run_beta from a FRESH kicked seed is still settling onto the
attractor (H decays 43→0→−3→+8 over 36τ), so the s-probe rides the decaying
envelope, not the saturated clock. The canonical ω₁ is measured on the SETTLED
state via 2h.run_dense (needs _m5_8_2g_settled.npz). ⇒ THIS DRIVER IS NOT VALID
for the absolute clock ω; the clean ω(δ) needs the 2h settled path (deferred).
The δ-conclusion does NOT depend on it — see m5_8_2q_delta_scaling.py + the ZBW
law (N-6a: ω ∝ H_rest), which with the measured δ-flat H_rest already pins the
ratio scaling ω·δ/(2H) ∝ δ. Kept as the scaffold for the proper 2h-path version.

USAGE (driven by a δ-loop in bash):
    M58_DELTA=0.3 M58_G=8 python m5_8_2q_omega.py [steps]
"""
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[4]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

# env is already set by the launcher → seed_M/D4 pick up δ, g
from openwave.xperiments.m5_liquid_crystal.research.sandbox_v8.m5_8_2cb_taichi_constrained import (  # noqa: E402
    build_grid, seed_M, make_masks, B_STAR, DT,
)
from openwave.xperiments.m5_liquid_crystal.research.sandbox_v8 import (  # noqa: E402
    m5_8_2d_quartic_saturation as d2,
)
from openwave.xperiments.m5_liquid_crystal.research.sandbox_vn.m5_8_2h_omega_attractor import (  # noqa: E402
    spec, peaks,
)

DELTA = float(os.environ.get("M58_DELTA", "0.3"))
G = float(os.environ.get("M58_G", "8.0"))
STEPS = int(sys.argv[1]) if len(sys.argv) > 1 else 15000
PROBE = d2.PROBE                  # s(t) sampled every PROBE steps
BETA = 1.558


def omega1(s, dsamp):
    """N-1 extractor: detrend dial over T∈(2.5,5,10), median of the top fundamental."""
    s = np.asarray(s, float)
    s = s[len(s) // 10:]                                  # drop the initial transient
    funda, tops = [], []
    for T in (2.5, 5.0, 10.0):
        n = max(int(T / dsamp) | 1, 3)
        pk = peaks(*spec(s, dsamp, detrend_n=n), 0.5, 3.0, n=4)
        if not pk:
            continue
        tops.append(pk[0][0])                             # strongest peak (may be 2ω₁)
        low = [w for w, _ in pk if w < 1.7]               # fundamental band
        if low:
            funda.append(low[0])
    f = float(np.median(funda)) if funda else float("nan")
    t = float(np.median(tops)) if tops else float("nan")
    return f, t


def main():
    g = build_grid()
    M0, Mth = seed_M(g, B_STAR)
    act, core = make_masks(g)
    h = g["h"]
    seed = (M0, Mth, act, core, g["rhat"], h)

    # H_static of this δ-seed (matches the 2q energy sweep)
    from openwave.xperiments.m5_liquid_crystal.research.sandbox_v8.m5_8_2c1_full_evolution import central, tw
    Mi = [central(M0, ax, h) for ax in range(3)]
    u = 0.0
    for i in range(3):
        for j in range(i + 1, 3):
            F = d2.np_commf(Mi[i], Mi[j])
            u = u + 2.0 * np.einsum("...ab,...ab->...", F, tw(F))
    H_static = float((u + BETA * u * u)[act > 0.5].sum()) * h ** 3

    dt = DT
    print(f"[2q-ω] δ={DELTA}  g={G}  steps={STEPS}  dt={dt}  H_static={H_static:.4f}")
    rec = d2.run_beta(f"d{DELTA}", BETA, seed, dt, STEPS)
    s = rec["s"]
    dsamp = PROBE * dt
    f, t = omega1(s, dsamp)
    print(f"OMEGA_RESULT delta={DELTA} g={G} omega1={f:.4f} omega_top={t:.4f} "
          f"H_static={H_static:.4f} ratio_hbarEQdelta={DELTA * f / (2 * H_static):.5f} "
          f"onset={rec['onset']} nsamp={len(s)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
