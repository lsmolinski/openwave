"""
M5.8.2w — ENERGY vs CLOCK FREQUENCY, E(ω).

The direct analog of the ϕ⁴-toymodel E(ω) minimum (arXiv:2501.04036): rest energy
as a function of the clock frequency. The 2u boost-sweep showed H_static(b_star)
has a minimum (activating the clock-fuel lowers the rest energy ~21%). Here we add
the FREQUENCY axis: at each boost dressing b_star we measure the clock frequency ω₁
from the dynamics (the 2q machinery) and the seed rest energy H_static, then plot
H_static vs ω₁, so the energy minimum is shown against ω, not the boost proxy.

ω-AXIS CAVEAT (honest): ω₁ here is the per-config probe frequency from a fresh-seed
constrained run (`d2.run_beta`), which rides the still-settling envelope rather than
the fully saturated clock (see m5_8_2q_omega.py: this reads low vs the settled
ω₁≈1.1). It is measured the SAME way at every b, so the SHAPE / ordering of E(ω) is
meaningful; the ABSOLUTE ω calibration needs the settled-clock 2h path (the named
refinement). Read this as the energy-vs-frequency MINIMUM structure, not a calibrated
ω axis.

Headless. Each b runs the dynamics (~80-100 s); writes data/ + plots/.
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
STEPS = int(sys.argv[1]) if len(sys.argv) > 1 else 12000
PROBE = d2.PROBE
OUT_NPZ = HERE / "data" / "_m5_8_2w_energy_vs_frequency.npz"
OUT_PNG = HERE / "plots" / "_m5_8_2w_energy_vs_frequency.png"

# boost-dressing sweep: fine through the 2u minimum (b*≈0.085), out to the production B_STAR
B_SWEEP = [0.0, 0.02, 0.04, 0.06, 0.085, 0.10, 0.115, 0.13]


def H_static_of(M0, act, h):
    Mi = [central(M0, ax, h) for ax in range(3)]
    u = 0.0
    for i in range(3):
        for j in range(i + 1, 3):
            F = d2.np_commf(Mi[i], Mi[j])
            u = u + 2.0 * np.einsum("...ab,...ab->...", F, tw(F))
    return float((u + BETA * u * u)[act > 0.5].sum()) * h ** 3


def main():
    print("=" * 78)
    print(f"M5.8.2w — E(ω): rest energy vs clock frequency  [δ={DELTA}, g={G}, {STEPS} steps]")
    print("=" * 78)
    g = build_grid()
    act, core = make_masks(g)
    h = g["h"]
    rhat = g["rhat"]
    dt = DT
    dsamp = PROBE * dt

    bs, Hs, ws = [], [], []
    for b in B_SWEEP:
        M0, Mth = seed_M(g, b)
        H_static = H_static_of(M0, act, h)
        seed = (M0, Mth, act, core, rhat, h)
        rec = d2.run_beta(f"b{b:g}", BETA, seed, dt, STEPS)
        w, wtop = omega1(rec["s"], dsamp)
        bs.append(b); Hs.append(H_static); ws.append(w)
        print(f"  b_star={b:6.3f}  H_static={H_static:9.4f}  omega1={w:7.4f}  omega_top={wtop:7.4f}")

    bs, Hs, ws = np.array(bs), np.array(Hs), np.array(ws)
    i_min = int(np.nanargmin(Hs))
    print(f"\n  energy minimum: H={Hs[i_min]:.4f} at b*={bs[i_min]:.3f}, omega1={ws[i_min]:.4f}")
    print(f"  H(b=0, clock stopped)={Hs[0]:.4f}; reduction at min = {Hs[0]-Hs[i_min]:.4f} "
          f"({100*(Hs[0]-Hs[i_min])/Hs[0]:.1f}%)")
    print("  CAVEAT: omega axis = per-config probe (unsettled envelope), consistent across the")
    print("  sweep so the E(omega) MINIMUM SHAPE is meaningful; absolute omega needs the 2h path.")

    np.savez(OUT_NPZ, b=bs, H_static=Hs, omega1=ws, delta=DELTA, g=G, steps=STEPS)

    order = np.argsort(ws)
    fig, ax = plt.subplots(figsize=(7, 4.6))
    ax.plot(ws[order], Hs[order], "o-", color="#1f77b4")
    ax.plot([ws[i_min]], [Hs[i_min]], "*", color="crimson", ms=16,
            label=f"min H={Hs[i_min]:.2f} (b*={bs[i_min]:.3f})")
    for b, w, H in zip(bs, ws, Hs):
        ax.annotate(f"b={b:g}", (w, H), fontsize=7, xytext=(3, 3), textcoords="offset points")
    ax.set_xlabel("clock frequency ω₁  (per-config probe; relative, see caveat)")
    ax.set_ylabel("rest energy H_static")
    ax.set_title("M5.8.2w  E(ω): the clock lowers rest energy to a minimum")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=110)
    print(f"  plot -> {OUT_PNG.name}\n  data -> {OUT_NPZ.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
