"""
M5.0i — Side-by-side profile vs M2 (laplace_propagation)

Times M2's `propagate_wave` (single fused kernel: leapfrog + trackers +
swap, all in one ndrange) against M5's split kernels at matching grid
sizes. Lets us see how much the M5 architectural cleanliness (separate
trackers / energy / swap) costs vs M2's all-in-one fusion.

Both use the same sample_avg_trackers 3-plane pattern.

NOTE: M2 stores L+T as two scalar fields (psiL_am, psiT_am — 24 bytes/voxel
total for the 6 buffers). M5 stores ψ as one Vector(3) field (36 bytes/voxel
across 3 buffers). M2 is structurally smaller per voxel.

USAGE:
    python -m openwave.xperiments.m5_lagrangian_field.research.m5_0i_profile_m2_compare
"""

import sys
import time
from pathlib import Path

import numpy as np
import taichi as ti

_THIS_DIR = Path(__file__).resolve().parent
_PROJECT_ROOT = _THIS_DIR.parents[3]
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from openwave.common import constants  # noqa: E402

# M2
from openwave.xperiments.m2_laplace_propagation import medium as m2_medium  # noqa: E402
from openwave.xperiments.m2_laplace_propagation import wave_engine as m2_engine  # noqa: E402

# M5
from openwave.xperiments.m5_lagrangian_field import lagrangian_engine as lagrange  # noqa: E402
from openwave.xperiments.m5_lagrangian_field import medium as m5_medium  # noqa: E402

UNIVERSE_EDGE = 1e-15  # m
GRID_SIZES = [128, 256, 384]
WARMUP_STEPS = 20
N_PROFILE_STEPS = 100


def stats(values):
    arr = np.asarray(values, dtype=np.float64)
    return {
        "mean": float(arr.mean()),
        "p50": float(np.median(arr)),
        "p95": float(np.percentile(arr, 95)),
    }


def profile_m2(target_voxels_per_axis):
    """Profile M2's per-step path at the given grid size."""
    target_voxels = target_voxels_per_axis**3
    wave_field = m2_medium.WaveField([UNIVERSE_EDGE, UNIVERSE_EDGE, UNIVERSE_EDGE], target_voxels)
    trackers = m2_medium.Trackers(wave_field.grid_size, scale_factor=1.0)
    nx = wave_field.nx

    cfl_safety = 0.95
    c_amrs = constants.WAVE_SPEED / constants.ATTOMETER * constants.RONTOSECOND
    dt_rs = wave_field.dx_am * cfl_safety / (c_amrs * np.sqrt(3.0))

    # Seed something non-trivial so the kernel does real work
    m2_engine.charge_gaussian(wave_field, 1.0)
    ti.sync()

    # Warmup
    for _ in range(WARMUP_STEPS):
        m2_engine.propagate_wave(wave_field, trackers, c_amrs, dt_rs, 0.0, 1.0)
    ti.sync()

    # Per-step
    step_ms = []
    for s in range(N_PROFILE_STEPS):
        t0 = time.perf_counter()
        m2_engine.propagate_wave(wave_field, trackers, c_amrs, dt_rs, float(s) * dt_rs, 1.0)
        ti.sync()
        step_ms.append((time.perf_counter() - t0) * 1000.0)

    return {"grid": (nx, nx, nx), "voxel_count": nx**3, "propagate_wave_ms": step_ms}


def profile_m5(target_voxels_per_axis):
    """Profile M5's per-step path (4 kernels + swap) at the given grid size."""
    target_voxels = target_voxels_per_axis**3
    wave_field = m5_medium.WaveField([UNIVERSE_EDGE, UNIVERSE_EDGE, UNIVERSE_EDGE], target_voxels)
    trackers = m5_medium.Trackers(wave_field.grid_size)
    nx = wave_field.nx

    cfl_safety = 0.95
    c_amrs = constants.WAVE_SPEED / constants.ATTOMETER * constants.RONTOSECOND
    dt_rs = wave_field.dx_am * cfl_safety / (c_amrs * np.sqrt(3.0))

    lagrange.seed_gaussian(wave_field, c_amrs, dt_rs, 5.0, 24.0, ti.Vector([0.0, 1.0, 0.0]), 0)
    ti.sync()

    for _ in range(WARMUP_STEPS):
        lagrange.evolve_psi(wave_field, c_amrs, dt_rs)
        wave_field.swap_buffers()
        lagrange.update_trackers(wave_field, trackers, dt_rs, 0.0)
        lagrange.compute_energyH_density(wave_field, trackers, c_amrs, dt_rs, 0.0)
    ti.sync()

    step_ms = []
    for s in range(N_PROFILE_STEPS):
        t0 = time.perf_counter()
        lagrange.evolve_psi(wave_field, c_amrs, dt_rs)
        wave_field.swap_buffers()
        lagrange.update_trackers(wave_field, trackers, dt_rs, float(s) * dt_rs)
        lagrange.compute_energyH_density(wave_field, trackers, c_amrs, dt_rs, 0.0)
        ti.sync()
        step_ms.append((time.perf_counter() - t0) * 1000.0)

    return {"grid": (nx, nx, nx), "voxel_count": nx**3, "step_ms": step_ms}


def main():
    ti.init(arch=ti.gpu)
    print(f"[M2 vs M5 compare]  grid sizes: {GRID_SIZES}  measure={N_PROFILE_STEPS}")

    rows = []
    for g in GRID_SIZES:
        print(f"\n--- {g}³ ---")
        try:
            m2 = profile_m2(g)
            m2s = stats(m2["propagate_wave_ms"])
            print(
                f"  M2 propagate_wave (fused):  mean={m2s['mean']:.2f}ms  "
                f"p95={m2s['p95']:.2f}ms"
            )
        except Exception as e:
            print(f"  M2 failed: {type(e).__name__}: {e}")
            m2s = None

        try:
            m5 = profile_m5(g)
            m5s = stats(m5["step_ms"])
            print(
                f"  M5 step (4 kernels+swap):   mean={m5s['mean']:.2f}ms  "
                f"p95={m5s['p95']:.2f}ms"
            )
        except Exception as e:
            print(f"  M5 failed: {type(e).__name__}: {e}")
            m5s = None

        if m2s and m5s:
            ratio = m5s["mean"] / m2s["mean"]
            print(f"  → M5 is {ratio:.2f}× M2 step time at {g}³")
            rows.append((g, m2s, m5s, ratio))

    print()
    print("=" * 72)
    print(f"{'grid':<8} {'M2 ms':>10} {'M5 ms':>10} {'M5/M2':>8} {'M2 fps':>10} {'M5 fps':>10}")
    print("-" * 72)
    for g, m2s, m5s, ratio in rows:
        m2_fps = 1000.0 / m2s["mean"]
        m5_fps = 1000.0 / m5s["mean"]
        print(
            f"{g}³{'':<3} {m2s['mean']:>10.2f} {m5s['mean']:>10.2f} "
            f"{ratio:>7.2f}× {m2_fps:>10.1f} {m5_fps:>10.1f}"
        )

    return 0


if __name__ == "__main__":
    sys.exit(main())
