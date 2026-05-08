"""
M5.0i — Per-step kernel profile

Times each per-step kernel of `compute_oscillation` separately at production
grid sizes to identify bottlenecks. Drives the Tier 2 optimization order.

KERNELS PROFILED (mirrors _launcher.compute_oscillation):
    propagate_psi              — leapfrog/Verlet step
    swap_buffers               — Python triple-buffer cycle (no kernel work)
    update_trackers_psi        — per-voxel amp / freq EMA
    compute_energy_density_H   — ½|ψ̇|² + ½c²|∇ψ|² + V(ψ)
    sample_avg_trackers        — 3-plane mean of amp / freq / energy (every 60 frames)

MEASUREMENT PROTOCOL:
    1. Build WaveField + Trackers at target grid size
    2. Seed Gaussian wave packet (realistic per-step load)
    3. WARMUP_STEPS to compile kernels + warm GPU clocks
    4. PER-KERNEL pass (N_PROFILE_STEPS): ti.sync() before/after each kernel,
       record ms per kernel. Slow but accurate breakdown.
    5. END-TO-END pass (N_PROFILE_STEPS): one ti.sync() per step, measure
       realistic step time without per-kernel sync overhead.
    6. Print table; flag any per-step total >50ms (=20fps target floor)
    7. Save bar chart per kernel × grid size

USAGE:
    python -m openwave.xperiments.m5_lagrangian_field.research.m5_0i_profile

Headless. Output: console table + research/plots/m5_0i_profile.png.
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
from openwave.xperiments.m5_lagrangian_field import lagrangian_engine as lagrange  # noqa: E402
from openwave.xperiments.m5_lagrangian_field import medium  # noqa: E402

# ================================================================
# CONFIG
# ================================================================
UNIVERSE_EDGE = 1e-15  # m, same as smoke test
GRID_SIZES = [128, 256, 384]  # nx=ny=nz target — actual lands on nearest odd
WARMUP_STEPS = 20
N_PROFILE_STEPS = 100  # enough for stable mean; runs ~30s at 256³ if step ~300ms
TARGET_STEP_MS = 50.0  # =20 fps floor; flag bottlenecks above this

# Seed parameters (Gaussian packet, parallel to _test_smoke)
SEED_AMPLITUDE_AM = 5.0
SEED_VOXELS_PER_LAMBDA = 24.0
SEED_POLARIZATION = ti.Vector([0.0, 1.0, 0.0])
SEED_DIRECTION = 0  # x-axis


# ================================================================
# PROFILER
# ================================================================
def profile_at_grid(target_voxels_per_axis):
    """Profile one grid size end-to-end. Returns dict of measurements."""
    target_voxels = target_voxels_per_axis**3
    print(f"\n{'=' * 72}")
    print(f"Grid target: {target_voxels_per_axis}³  (~{target_voxels / 1e6:.1f}M voxels)")
    print(f"{'=' * 72}")

    # ── Build grid + trackers ────────────────────────────────────
    wave_field = medium.WaveField(
        [UNIVERSE_EDGE, UNIVERSE_EDGE, UNIVERSE_EDGE], target_voxels
    )
    trackers = medium.Trackers(wave_field.grid_size)
    nx, ny, nz = wave_field.nx, wave_field.ny, wave_field.nz
    print(f"  actual grid: ({nx}, {ny}, {nz})  dx={wave_field.dx_am:.3f} am")

    # ── CFL-stable timestep ──────────────────────────────────────
    cfl_safety = 0.95
    c_amrs = constants.WAVE_SPEED / constants.ATTOMETER * constants.RONTOSECOND
    dt_rs = wave_field.dx_am * cfl_safety / (c_amrs * np.sqrt(3.0))
    cfl_factor = (c_amrs * dt_rs / wave_field.dx_am) ** 2
    print(f"  c={c_amrs:.4f} am/rs  dt={dt_rs:.4f} rs  cfl={cfl_factor:.4f}")

    # ── Seed Gaussian packet (realistic per-step load) ──────────
    lagrange.seed_gaussian(
        wave_field,
        c_amrs,
        dt_rs,
        SEED_AMPLITUDE_AM,
        SEED_VOXELS_PER_LAMBDA,
        SEED_POLARIZATION,
        SEED_DIRECTION,
    )
    ti.sync()

    # ── Warmup (kernel compile + GPU clock ramp) ────────────────
    print(f"  warmup ({WARMUP_STEPS} steps)...", flush=True)
    for _ in range(WARMUP_STEPS):
        lagrange.propagate_psi(wave_field, c_amrs, dt_rs)
        wave_field.swap_buffers()
        lagrange.update_trackers_psi(wave_field, trackers, dt_rs, 0.0)
        lagrange.compute_energy_density_H(wave_field, trackers, c_amrs, dt_rs)
    ti.sync()

    # ── PER-KERNEL pass: sync between every kernel ─────────────
    print(f"  per-kernel pass ({N_PROFILE_STEPS} steps)...", flush=True)
    times = {
        "propagate_psi": [],
        "swap_buffers": [],
        "update_trackers_psi": [],
        "compute_energy_density_H": [],
    }

    for step in range(N_PROFILE_STEPS):
        ti.sync()
        t0 = time.perf_counter()
        lagrange.propagate_psi(wave_field, c_amrs, dt_rs)
        ti.sync()
        t1 = time.perf_counter()
        wave_field.swap_buffers()
        ti.sync()
        t2 = time.perf_counter()
        lagrange.update_trackers_psi(wave_field, trackers, dt_rs, float(step) * dt_rs)
        ti.sync()
        t3 = time.perf_counter()
        lagrange.compute_energy_density_H(wave_field, trackers, c_amrs, dt_rs)
        ti.sync()
        t4 = time.perf_counter()

        times["propagate_psi"].append((t1 - t0) * 1000.0)
        times["swap_buffers"].append((t2 - t1) * 1000.0)
        times["update_trackers_psi"].append((t3 - t2) * 1000.0)
        times["compute_energy_density_H"].append((t4 - t3) * 1000.0)

    # ── sample_avg_trackers — measured separately (every-60-frames cadence) ──
    sample_times_ms = []
    for _ in range(10):  # 10 samples is plenty for this kernel
        ti.sync()
        t0 = time.perf_counter()
        lagrange.sample_avg_trackers(wave_field, trackers)
        ti.sync()
        sample_times_ms.append((time.perf_counter() - t0) * 1000.0)
    times["sample_avg_trackers (per call)"] = sample_times_ms

    # ── END-TO-END pass: one sync per step ──────────────────────
    print(f"  end-to-end pass ({N_PROFILE_STEPS} steps, single sync/step)...", flush=True)
    end_to_end_ms = []
    for step in range(N_PROFILE_STEPS):
        t0 = time.perf_counter()
        lagrange.propagate_psi(wave_field, c_amrs, dt_rs)
        wave_field.swap_buffers()
        lagrange.update_trackers_psi(wave_field, trackers, dt_rs, float(step) * dt_rs)
        lagrange.compute_energy_density_H(wave_field, trackers, c_amrs, dt_rs)
        ti.sync()
        end_to_end_ms.append((time.perf_counter() - t0) * 1000.0)

    return {
        "grid": (nx, ny, nz),
        "voxel_count": nx * ny * nz,
        "per_kernel": times,
        "end_to_end_ms": end_to_end_ms,
    }


def stats(values):
    arr = np.asarray(values, dtype=np.float64)
    return {
        "mean": float(arr.mean()),
        "p50": float(np.median(arr)),
        "p95": float(np.percentile(arr, 95)),
    }


def print_report(result):
    nx, ny, nz = result["grid"]
    voxel_count = result["voxel_count"]
    per_kernel = result["per_kernel"]
    end_to_end = result["end_to_end_ms"]

    e2e_stats = stats(end_to_end)
    fps_mean = 1000.0 / e2e_stats["mean"] if e2e_stats["mean"] > 0 else 0.0
    fps_p95 = 1000.0 / e2e_stats["p95"] if e2e_stats["p95"] > 0 else 0.0
    sum_per_kernel_mean = sum(stats(v)["mean"] for k, v in per_kernel.items() if k != "sample_avg_trackers (per call)")

    print()
    print(f"--- {nx}×{ny}×{nz}  ({voxel_count / 1e6:.1f}M voxels) ---")
    print(f"{'kernel':<32} {'mean ms':>10} {'p50 ms':>10} {'p95 ms':>10} {'% of step':>10}")
    print("-" * 76)
    for kname, samples in per_kernel.items():
        s = stats(samples)
        if kname == "sample_avg_trackers (per call)":
            pct = ""  # not part of every step
        else:
            pct = f"{100.0 * s['mean'] / sum_per_kernel_mean:.1f}%"
        print(f"{kname:<32} {s['mean']:>10.3f} {s['p50']:>10.3f} {s['p95']:>10.3f} {pct:>10}")
    print("-" * 76)
    print(f"{'sum (synced per kernel)':<32} {sum_per_kernel_mean:>10.3f}")
    print(
        f"{'end-to-end (1 sync/step)':<32} {e2e_stats['mean']:>10.3f} "
        f"{e2e_stats['p50']:>10.3f} {e2e_stats['p95']:>10.3f}"
    )
    print(f"  → fps mean={fps_mean:.1f}  fps@p95={fps_p95:.1f}")
    if e2e_stats["mean"] > TARGET_STEP_MS:
        over = e2e_stats["mean"] - TARGET_STEP_MS
        print(
            f"  ⚠ over 20fps budget by {over:.1f} ms "
            f"({100.0 * over / TARGET_STEP_MS:.0f}% over)"
        )
    else:
        print(f"  ✓ under 20fps budget ({TARGET_STEP_MS:.0f} ms/step)")


# ================================================================
# PLOT
# ================================================================
def plot_results(results):
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("\n[warn] matplotlib not available — skipping plot")
        return

    plot_path = _THIS_DIR / "plots" / "m5_0i_profile.png"
    plot_path.parent.mkdir(parents=True, exist_ok=True)

    kernel_order = [
        "propagate_psi",
        "swap_buffers",
        "update_trackers_psi",
        "compute_energy_density_H",
    ]
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"]
    grid_labels = [f"{r['grid'][0]}³" for r in results]
    n_grids = len(results)

    fig, ax = plt.subplots(figsize=(9, 5.5))
    bottom = np.zeros(n_grids)
    for kname, color in zip(kernel_order, colors):
        means = np.array([stats(r["per_kernel"][kname])["mean"] for r in results])
        ax.bar(grid_labels, means, bottom=bottom, label=kname, color=color)
        bottom += means

    # End-to-end overlay (red dashed at top)
    e2e_means = np.array([stats(r["end_to_end_ms"])["mean"] for r in results])
    ax.scatter(
        grid_labels,
        e2e_means,
        marker="_",
        s=400,
        color="black",
        zorder=10,
        label="end-to-end mean (1 sync/step)",
    )

    ax.axhline(
        TARGET_STEP_MS,
        color="red",
        linestyle="--",
        alpha=0.6,
        label=f"20 fps budget ({TARGET_STEP_MS:.0f} ms)",
    )
    ax.set_ylabel("ms per step")
    ax.set_title("M5.0i — Per-step kernel profile (synced per kernel, stacked)")
    ax.legend(loc="upper left", fontsize=9)
    ax.grid(alpha=0.3, axis="y")

    fig.tight_layout()
    fig.savefig(str(plot_path), dpi=120)
    print(f"\nplot saved: {plot_path}")


# ================================================================
# MAIN
# ================================================================
def main():
    ti.init(arch=ti.gpu)
    print(f"[M5.0i] grid sizes: {GRID_SIZES}")
    print(f"[M5.0i] warmup={WARMUP_STEPS}  measure={N_PROFILE_STEPS}")
    print(f"[M5.0i] target: ≤{TARGET_STEP_MS:.0f} ms/step (20 fps)")

    results = []
    for grid_size in GRID_SIZES:
        try:
            r = profile_at_grid(grid_size)
            print_report(r)
            results.append(r)
        except Exception as e:
            print(f"\n[ERROR] grid {grid_size}³ failed: {type(e).__name__}: {e}")
            print("  (likely out-of-memory; skipping larger grids)")
            break

    if results:
        plot_results(results)

    return 0


if __name__ == "__main__":
    sys.exit(main())
