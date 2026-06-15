"""
M5.8.2z — ACCURATE E(ω): rest energy vs the SETTLED clock frequency.

The honest version of 2w. 2w used a short per-config probe whose ω rode the
still-settling transient (scattered, one nan). Here each boost-dressing point
runs LONG enough to saturate onto the molten-clock attractor, and ω₁ is measured
on the SETTLED TAIL (last 40% of s), not the transient. H_static is the seed
rest energy at that dressing (same as 2u/2w).

Physics to expect (from 2o / N-6a): the saturated clock's ω₁ is RIGID
(state-dependent mode composition, ~1.1), NOT a smooth function of E. So this is
expected to show the energy MINIMUM (the 2u 21% dip) sitting at a roughly-rigid
ω, which is a real DIFFERENCE from the ϕ⁴-toymodel's smooth E(ω) minimum, the
clock frequency does not track energy here. Reported either way.

Long run (~5-6 min/point, ~45 min total). Headless. Writes data/ + plots/.
"""
import os
import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

os.environ.setdefault("M58_DELTA", "0.3")
os.environ.setdefault("M58_G", "8.0")

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[4]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from openwave.xperiments.m5_liquid_crystal.research.sandbox_v8.m5_8_2cb_taichi_constrained import (  # noqa: E402
    build_grid, seed_M, make_masks, DT,
)
from openwave.xperiments.m5_liquid_crystal.research.sandbox_v8 import (  # noqa: E402
    m5_8_2d_quartic_saturation as d2,
)
from openwave.xperiments.m5_liquid_crystal.research.sandbox_v8.m5_8_2c1_full_evolution import (  # noqa: E402
    central, tw,
)
from openwave.xperiments.m5_liquid_crystal.research.sandbox_vn.m5_8_2q_omega import omega1  # noqa: E402

DELTA = float(os.environ["M58_DELTA"])
G = float(os.environ["M58_G"])
BETA = 1.558
STEPS = int(sys.argv[1]) if len(sys.argv) > 1 else 50000
TAIL = 0.40                       # measure omega on the last 40% (settled)
PROBE = d2.PROBE
OUT_NPZ = HERE / "data" / "_m5_8_2z_settled_energy_vs_frequency.npz"
OUT_PNG = HERE / "plots" / "_m5_8_2z_settled_energy_vs_frequency.png"

B_SWEEP = [0.0, 0.02, 0.04, 0.06, 0.085, 0.10, 0.115, 0.13]


def H_static_of(M0, act, h):
    Mi = [central(M0, ax, h) for ax in range(3)]
    u = 0.0
    for i in range(3):
        for j in range(i + 1, 3):
            F = d2.np_commf(Mi[i], Mi[j])
            u = u + 2.0 * np.einsum("...ab,...ab->...", F, tw(F))
    return float((u + BETA * u * u)[act > 0.5].sum()) * h ** 3


def settled_omega(s, dsamp):
    s = np.asarray(s, float)
    tail = s[int((1.0 - TAIL) * len(s)):]          # the settled portion
    return omega1(tail, dsamp)


def main():
    print("=" * 78)
    print(f"M5.8.2z — ACCURATE E(ω) (settled): δ={DELTA}, g={G}, {STEPS} steps/point, "
          f"ω on last {int(TAIL*100)}%")
    print("=" * 78)
    g = build_grid()
    act, core = make_masks(g)
    h = g["h"]
    rhat = g["rhat"]
    dsamp = PROBE * DT

    bs, Hs, ws, wtops = [], [], [], []
    for b in B_SWEEP:
        M0, Mth = seed_M(g, b)
        H_static = H_static_of(M0, act, h)
        seed = (M0, Mth, act, core, rhat, h)
        rec = d2.run_beta(f"b{b:g}", BETA, seed, DT, STEPS)
        w, wtop = settled_omega(rec["s"], dsamp)
        bs.append(b); Hs.append(H_static); ws.append(w); wtops.append(wtop)
        print(f"  b_star={b:6.3f}  H_static={H_static:9.4f}  omega1_settled={w:7.4f}  "
              f"omega_top={wtop:7.4f}", flush=True)

    bs, Hs, ws = np.array(bs), np.array(Hs), np.array(ws)
    i_min = int(np.nanargmin(Hs))
    wvalid = ws[~np.isnan(ws)]
    print(f"\n  energy minimum: H={Hs[i_min]:.4f} at b*={bs[i_min]:.3f}, "
          f"omega_settled={ws[i_min]:.4f}")
    print(f"  settled omega spread: {np.nanmin(ws):.4f} - {np.nanmax(ws):.4f} "
          f"(std {np.nanstd(wvalid):.4f}) over a {Hs.max()-Hs.min():.2f} energy range")
    print("  READ: if omega is ~rigid while H varies, the M5.8 clock frequency does NOT")
    print("  track energy (N-6a rigidity) - a real difference from the toymodel's E(ω).")

    np.savez(OUT_NPZ, b=bs, H_static=Hs, omega1=ws, omega_top=np.array(wtops),
             delta=DELTA, g=G, steps=STEPS, tail=TAIL)

    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    m = ~np.isnan(ws)
    ax.plot(ws[m], Hs[m], "o", ms=9, color="#1f77b4")
    ax.plot([ws[i_min]], [Hs[i_min]], "*", color="crimson", ms=18,
            label=f"energy min H={Hs[i_min]:.2f} (b*={bs[i_min]:.3f})")
    for b, w, H in zip(bs, ws, Hs):
        if not np.isnan(w):
            ax.annotate(f"b={b:g}", (w, H), fontsize=7, xytext=(4, 4), textcoords="offset points")
    ax.set_xlabel("settled clock frequency ω₁  (saturated attractor)")
    ax.set_ylabel("rest energy H_static")
    ax.set_title("M5.8.2z  accurate E(ω): energy minimum vs the SETTLED clock frequency")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=110)
    print(f"  plot -> {OUT_PNG.name}\n  data -> {OUT_NPZ.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
