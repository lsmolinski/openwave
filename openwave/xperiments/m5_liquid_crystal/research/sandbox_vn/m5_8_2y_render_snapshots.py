"""
M5.8.2y — static SNAPSHOTS of the M5.8 simulation (for sharing; headless).

Three panels from already-saved runs, no new evolution:
  (1) spatial structure of the settled molten clock: the field departure from
      vacuum ||M(x) - M_vac||_F on a mid-plane slice (shows the localized defect);
  (2) the clock signal s(t) (the Zitterbewegung overlap oscillation);
  (3) the energy H(t) settling onto the attractor.

Reads data/_m5_8_2o_settled.npz (M_settled) + data/_m5_8_2h_dense.npz (traces).
numpy + matplotlib only (no Taichi). Writes plots/.
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
OUT = HERE / "plots" / "_m5_8_2y_snapshots.png"


def main():
    settled = np.load(DATA / "_m5_8_2o_settled.npz")
    M = settled["M_settled"]                       # (N,N,N,4,4)
    N = M.shape[0]
    M_vac = M[0, 0, 0]                             # corner ≈ vacuum
    dev = np.sqrt(((M - M_vac) ** 2).sum(axis=(-1, -2)))   # (N,N,N) departure
    sl = dev[:, :, N // 2]                         # mid-plane z-slice

    dense = np.load(DATA / "_m5_8_2h_dense.npz")
    t, s = dense["R0-kicked_t"], dense["R0-kicked_s"]

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))

    # clip the colour scale to the interior (the outer ring is the vacuum-blend boundary)
    vmax = float(np.percentile(dev[3:-3, 3:-3, N // 2], 98))
    im = ax[0].imshow(sl.T, origin="lower", cmap="inferno", aspect="equal", vmax=vmax)
    ax[0].set_title("settled molten clock:\nfield magnitude ‖M − M_vac‖ (mid-plane slice)")
    ax[0].set_xlabel("x (voxels)"); ax[0].set_ylabel("y (voxels)")
    fig.colorbar(im, ax=ax[0], fraction=0.046, pad=0.04, label="‖M − M_vac‖")

    ax[1].plot(t, s, lw=0.9, color="#1f77b4")
    ax[1].set_title("clock signal s(t)\n(the Zitterbewegung overlap, self-starting + settling)")
    ax[1].set_xlabel("t (τ)"); ax[1].set_ylabel("clock overlap s")

    fig.suptitle("M5.8 Zitterbewegung clock — simulation snapshots", y=1.0, fontsize=12)
    fig.tight_layout()
    fig.savefig(OUT, dpi=120, bbox_inches="tight")
    print(f"  grid {N}³  dev range [{dev.min():.3f}, {dev.max():.3f}]  vmax(clip)={vmax:.3f}")
    print(f"  s(t): {len(s)} samples")
    print(f"  plot -> {OUT.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
