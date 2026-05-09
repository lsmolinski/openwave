"""
M5.0h — Dispersion-Relation Regression Test

Validates that propagate_psi reproduces the discrete leapfrog dispersion

    ω² = (c²/dx²) · 2 · Σ_α (1 − cos(k_α · dx))

within 0.5% per-mode c² recovery, gating progression to M5.1.

METHOD:
    1. Seed 5 Dirichlet standing-wave eigenmodes simultaneously (linear
       superposition — V=0 means modes don't interact).
    2. Propagate for N_STEPS via leapfrog. Each step, project ψ onto each
       known eigenmode shape → 5 time series.
    3. Per-mode temporal FFT → ω_measured. Parabolic peak interpolation
       gives sub-bin precision.
    4. Solve discrete-dispersion equation for c²: per mode,
           c²_recovered = ω²_measured · dx² / [2 Σ_α (1−cos(k_α·dx))]
       Compare against c²_target = c_amrs² (constants.WAVE_SPEED_AMRS²).
    5. Pass iff max |c²_recovered/c²_target − 1| < 0.5% over all modes.

OUTPUT:
    - Console fit table + PASS/FAIL
    - PNG: ω(k) measured vs theoretical curve + per-mode residual bar chart
      → ./plots/m5_0h_dispersion.png

EXIT CODE:
    0 on PASS, 1 on FAIL (CI-friendly).

USAGE:
    python -m openwave.xperiments.m5_lagrangian_field.research.m5_0h_dispersion

Headless: no GUI window (would slow the run by ~70× and serves no purpose
for a regression gate). PNG is the visual artifact for inspection.
"""

import sys
from pathlib import Path

import numpy as np
import taichi as ti

# Make sibling modules importable when run as a script
_THIS_DIR = Path(__file__).resolve().parent
_PROJECT_ROOT = _THIS_DIR.parents[3]
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from openwave.common import constants  # noqa: E402
from openwave.xperiments.m5_lagrangian_field import lagrangian_engine as lagrange  # noqa: E402
from openwave.xperiments.m5_lagrangian_field import medium  # noqa: E402


# ================================================================
# POINT-SAMPLING (avoids full-grid GPU reduction)
# ================================================================
# Lesson from M2's `3-PLANE SAMPLING FOR AVERAGE TRACKERS`: full-grid
# atomic_add reductions on GPU stall on atomic contention at million-voxel
# scale. Even Taichi's auto-reduction lowers to atomics on Metal — the test
# stalled at ~25 steps/min with GPU pinned at 100% on a sum reduction.
#
# Workaround: sample ψ at a handful of strategically-chosen points (one
# antinode per mode of interest), each a single voxel read. Linear
# superposition means ψ_z(t) at any non-nodal point carries ALL mode
# frequencies — temporal FFT recovers ω_m by peak-picking in a narrow band
# around theoretical ω_m. Cost: O(N_SAMPLES) per step, not O(N³).


@ti.kernel
def record_sample_points(
    wave_field: ti.template(),  # type: ignore
    sample_points: ti.template(),  # type: ignore — shape (N_SAMPLES, 3) of i32
    polarization: ti.template(),  # type: ignore
    history: ti.template(),  # type: ignore — shape (N_STEPS, N_SAMPLES)
    step: ti.i32,  # type: ignore
):
    """Read ψ·ê at each sample point and store into history[step, s].

    Single ti.field read per sample × N_SAMPLES samples per step. No
    reduction, no atomics. The full-grid spatial structure is encoded
    in the choice of sample location (one per mode antinode).
    """
    n_samples = sample_points.shape[0]
    for s in range(n_samples):
        i = sample_points[s, 0]
        j = sample_points[s, 1]
        k = sample_points[s, 2]
        history[step, s] = wave_field.psi_am[i, j, k].dot(polarization)

# ================================================================
# TEST PARAMETERS
# ================================================================
UNIVERSE_EDGE = 1e-15  # m, same physical scale as smoke test
TARGET_VOXELS = 64**3  # ~262k voxels — small enough for ~30s run
N_STEPS = 4096  # FFT bin width / ω ~ 0.02% for slowest mode (with parabolic interp)
TOLERANCE_C2_PCT = 0.5  # pass criterion: |c²_recovered/c²_true − 1| < this %

# 5 Dirichlet standing-wave modes: vary x-axis index, transverse fixed at lowest mode (1).
# Spans ~30 → 8 voxels/wavelength on the x-axis.
MODE_INDICES = np.array(
    [
        [4, 1, 1],
        [6, 1, 1],
        [8, 1, 1],
        [12, 1, 1],
        [16, 1, 1],
    ],
    dtype=np.int32,
)
N_MODES = len(MODE_INDICES)
SEED_AMPLITUDE_AM = 1.0  # peak amplitude per mode (am)


# ================================================================
# HELPERS
# ================================================================
def discrete_K_dx2(n_triple, nm1):
    """Discrete-Laplacian eigenvalue scaled by dx²:  K·dx² = 2·Σ_α(1−cos(k_α·dx))."""
    kdx = np.pi * np.asarray(n_triple, dtype=np.float64) / nm1
    return float(2.0 * np.sum(1.0 - np.cos(kdx)))


def leapfrog_omega(n_triple, c_amrs, dx_am, dt_rs, nm1):
    """Theoretical ω for the FULL leapfrog scheme (both spatial AND temporal
    discretization).

    Spatial-only would be ω² = c²·K.  Leapfrog adds a temporal correction:
        sin(ω·dt/2) = (c·dt/2)·√K
        ω           = (2/dt)·arcsin((c·dt/2)·√K)

    For small (c·dt·√K), arcsin(x) ≈ x·(1 + x²/6), giving
        ω_leapfrog ≈ c·√K · (1 + cfl_factor·K·dx²/24)
    so the per-mode correction grows ∝ K·dx² ≈ (k·dx)² — exactly the
    systematic error we observed when fitting against the spatial-only ω.
    """
    K_dx2 = discrete_K_dx2(n_triple, nm1)  # K·dx² (dimensionless)
    arg = (c_amrs * dt_rs / (2.0 * dx_am)) * np.sqrt(K_dx2)
    return float((2.0 / dt_rs) * np.arcsin(arg))


def parabolic_peak_bin(spectrum, peak_bin):
    """Sub-bin peak refinement via parabolic fit on 3 points around peak."""
    if peak_bin <= 0 or peak_bin >= len(spectrum) - 1:
        return float(peak_bin)
    a = float(spectrum[peak_bin - 1])
    b = float(spectrum[peak_bin])
    c = float(spectrum[peak_bin + 1])
    denom = a - 2.0 * b + c
    if abs(denom) < 1e-12:
        return float(peak_bin)
    return float(peak_bin) + 0.5 * (a - c) / denom


# ================================================================
# MAIN
# ================================================================
def main():
    ti.init(arch=ti.gpu)

    # ── Build grid ───────────────────────────────────────────────
    wave_field = medium.WaveField(
        [UNIVERSE_EDGE, UNIVERSE_EDGE, UNIVERSE_EDGE], TARGET_VOXELS
    )
    # Allocate trackers (unused, but constructor keeps any shared expectations)
    _ = medium.Trackers(wave_field.grid_size)

    nx = wave_field.nx
    ny = wave_field.ny
    nz = wave_field.nz
    if not (nx == ny == nz):
        raise RuntimeError(f"Test assumes cubic grid; got ({nx},{ny},{nz})")
    N = nx
    nm1 = N - 1
    dx_am = wave_field.dx_am

    # ── Physics ──────────────────────────────────────────────────
    sim_speed = 1.0
    cfl_safety = 0.95
    c_amrs = (
        constants.WAVE_SPEED / constants.ATTOMETER * constants.RONTOSECOND * sim_speed
    )
    dt_rs = dx_am * cfl_safety / (c_amrs / sim_speed * np.sqrt(3.0))
    cfl_factor = (c_amrs * dt_rs / dx_am) ** 2

    # ── Theoretical per-mode ω (full leapfrog: space + time) ────
    omegas_theory = np.array(
        [leapfrog_omega(M, c_amrs, dx_am, dt_rs, nm1) for M in MODE_INDICES]
    )
    voxels_per_lambda_x = 2.0 * nm1 / MODE_INDICES[:, 0]

    # ── Allocate Taichi buffers for seed/projection ──────────────
    mode_indices_field = ti.field(dtype=ti.i32, shape=(N_MODES, 3))
    mode_amp_field = ti.field(dtype=ti.f32, shape=(N_MODES,))
    mode_cos_omega_dt_field = ti.field(dtype=ti.f32, shape=(N_MODES,))

    mode_indices_field.from_numpy(MODE_INDICES)
    mode_amp_field.from_numpy(np.full(N_MODES, SEED_AMPLITUDE_AM, dtype=np.float32))
    mode_cos_omega_dt_field.from_numpy(
        np.cos(omegas_theory * dt_rs).astype(np.float32)
    )

    # ── Sample-point selection (one antinode per mode) ────────────
    # Each mode m has its 1st x-antinode at i ≈ (N−1) / (2·n_x). Place sample
    # at j=k=mid for transverse antinode (sin(π·mid/(N−1)) ≈ 1).
    j_mid = (N - 1) // 2
    k_mid = (N - 1) // 2
    sample_points_np = np.zeros((N_MODES, 3), dtype=np.int32)
    for m in range(N_MODES):
        nx_m = MODE_INDICES[m, 0]
        i_anti = max(1, min(N - 2, int(round(nm1 / (2.0 * nx_m)))))
        sample_points_np[m] = [i_anti, j_mid, k_mid]
    sample_points_field = ti.field(dtype=ti.i32, shape=(N_MODES, 3))
    sample_points_field.from_numpy(sample_points_np)
    history_field = ti.field(dtype=ti.f32, shape=(N_STEPS, N_MODES))

    # ê_z polarization — decouples from any mode-axis structure
    polarization = ti.Vector([0.0, 0.0, 1.0])

    # ── Seed initial condition ───────────────────────────────────
    lagrange.seed_dispersion_modes(
        wave_field,
        mode_indices_field,
        mode_amp_field,
        mode_cos_omega_dt_field,
        polarization,
    )

    # ── Run + record per-step ψ at sample points ─────────────────
    print(
        f"[M5.0h] grid={N}³  dx={dx_am:.4f} am  dt={dt_rs:.4f} rs  "
        f"cfl_factor={cfl_factor:.4f}  steps={N_STEPS}"
    )
    print(
        f"[M5.0h] sample points: {sample_points_np.tolist()}",
        flush=True,
    )
    print(f"[M5.0h] propagating...", flush=True)
    progress_every = max(1, N_STEPS // 20)
    for step in range(N_STEPS):
        lagrange.propagate_psi(wave_field, c_amrs, dt_rs)
        wave_field.swap_buffers()
        record_sample_points(
            wave_field,
            sample_points_field,
            polarization,
            history_field,
            step,
        )
        if (step + 1) % progress_every == 0:
            ti.sync()
            print(f"  step {step + 1}/{N_STEPS}", flush=True)

    # Single GPU→CPU sync after the whole run
    ti.sync()
    sample_history = history_field.to_numpy().astype(np.float64)

    # ── Extract ω_measured per mode via FFT + windowed peak ──────
    # Each sample point's time series carries ALL mode frequencies (linear
    # superposition). For mode m, we use the time series at the sample point
    # placed at mode m's antinode (where mode m's contribution is strongest)
    # and find the FFT peak inside a ±10% window around theoretical ω_m.
    measured_omegas = np.zeros(N_MODES)
    T_total = N_STEPS * dt_rs
    bin_to_omega = 2.0 * np.pi / T_total  # ω per FFT bin

    WINDOW_PCT = 0.10  # search ±10% around theoretical ω for peak picking
    for m in range(N_MODES):
        signal = sample_history[:, m] - sample_history[:, m].mean()
        spectrum = np.abs(np.fft.rfft(signal))
        omega_lo = omegas_theory[m] * (1.0 - WINDOW_PCT)
        omega_hi = omegas_theory[m] * (1.0 + WINDOW_PCT)
        bin_lo = max(1, int(np.floor(omega_lo / bin_to_omega)))
        bin_hi = min(len(spectrum) - 1, int(np.ceil(omega_hi / bin_to_omega)))
        if bin_hi <= bin_lo:
            bin_hi = bin_lo + 1
        peak_offset = int(np.argmax(spectrum[bin_lo : bin_hi + 1]))
        peak_bin = bin_lo + peak_offset
        refined_bin = parabolic_peak_bin(spectrum, peak_bin)
        measured_omegas[m] = refined_bin * bin_to_omega

    # ── Recover c² per mode (invert full leapfrog dispersion) ────
    # sin(ω·dt/2) = (c·dt/2)·√K   ⇒  c² = (4/dt²/K)·sin²(ω_d·dt/2)
    # where K = (1/dx²)·(K·dx²) and (K·dx²) is the dimensionless Laplacian
    # eigenvalue from discrete_K_dx2.
    c2_target = c_amrs**2
    c2_recovered = np.zeros(N_MODES)
    for m in range(N_MODES):
        K_dx2 = discrete_K_dx2(MODE_INDICES[m], nm1)
        sin_half = np.sin(measured_omegas[m] * dt_rs / 2.0)
        c2_recovered[m] = 4.0 * (dx_am**2) * (sin_half**2) / (dt_rs**2 * K_dx2)

    err_pct = 100.0 * (c2_recovered / c2_target - 1.0)

    # ── Print fit table ──────────────────────────────────────────
    print()
    print(
        f"{'mode':<14} {'voxels/λ_x':>10} {'k_x·dx':>8} "
        f"{'ω_theory':>14} {'ω_measured':>14} "
        f"{'c²_rec/c²':>11} {'err %':>9}  pass"
    )
    print("-" * 96)
    all_pass = True
    for m in range(N_MODES):
        kdx_x = np.pi * MODE_INDICES[m, 0] / nm1
        passed = abs(err_pct[m]) < TOLERANCE_C2_PCT
        all_pass = all_pass and passed
        marker = "PASS" if passed else "FAIL"
        n_str = (
            f"({MODE_INDICES[m, 0]},{MODE_INDICES[m, 1]},{MODE_INDICES[m, 2]})"
        )
        print(
            f"{n_str:<14} {voxels_per_lambda_x[m]:>10.2f} {kdx_x:>8.4f} "
            f"{omegas_theory[m]:>14.6e} {measured_omegas[m]:>14.6e} "
            f"{c2_recovered[m] / c2_target:>11.6f} {err_pct[m]:>+9.4f}  {marker}"
        )
    print("-" * 96)
    print(f"c_target  = {c_amrs:.10f} am/rs")
    print(f"c²_target = {c2_target:.10e}")
    print(f"tolerance = ±{TOLERANCE_C2_PCT}% on c²")
    print()
    print(f"{'OVERALL: PASS' if all_pass else 'OVERALL: FAIL'}")

    # ── Plot ──────────────────────────────────────────────────────
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("\n[warn] matplotlib not available — skipping plot")
        return 0 if all_pass else 1

    plot_path = _THIS_DIR / "plots" / "m5_0h_dispersion.png"
    plot_path.parent.mkdir(parents=True, exist_ok=True)

    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1, figsize=(8.5, 7.0), gridspec_kw={"height_ratios": [3, 1]}
    )

    # k axis: physical wavenumber (rad/am) along x
    k_modes = np.pi * MODE_INDICES[:, 0] / (nm1 * dx_am)
    k_smooth = np.linspace(k_modes.min() * 0.9, k_modes.max() * 1.1, 200)

    # Theoretical curves — both sweep n_x with transverse fixed at 1
    n_x_smooth = k_smooth * nm1 * dx_am / np.pi
    omegas_theory_smooth = np.array(
        [leapfrog_omega([n, 1, 1], c_amrs, dx_am, dt_rs, nm1) for n in n_x_smooth]
    )
    omegas_continuum = c_amrs * k_smooth

    ax_top.plot(
        k_smooth, omegas_continuum, "k--", alpha=0.4, lw=1.0, label=r"$\omega = c k$ (continuum)"
    )
    ax_top.plot(
        k_smooth,
        omegas_theory_smooth,
        "b-",
        alpha=0.6,
        lw=1.5,
        label=r"$\omega(k)$ leapfrog (space+time)",
    )
    ax_top.scatter(
        k_modes,
        measured_omegas,
        color="red",
        s=70,
        zorder=10,
        label="measured",
        edgecolors="black",
        linewidths=0.5,
    )
    ax_top.set_xlabel(r"$k_x$ (rad / am)")
    ax_top.set_ylabel(r"$\omega$ (rad / rs)")
    title_status = "PASS" if all_pass else "FAIL"
    ax_top.set_title(
        f"M5.0h Dispersion Test  [{title_status}]   "
        f"N={N}³, {N_STEPS} steps, dt={dt_rs:.2f} rs"
    )
    ax_top.legend(loc="lower right")
    ax_top.grid(alpha=0.3)

    bar_colors = [
        "tab:green" if abs(e) < TOLERANCE_C2_PCT else "tab:red" for e in err_pct
    ]
    xs = np.arange(N_MODES)
    ax_bot.bar(xs, err_pct, color=bar_colors)
    ax_bot.axhline(y=0, color="k", linewidth=0.5)
    ax_bot.axhline(
        y=TOLERANCE_C2_PCT, color="tab:red", linestyle=":", alpha=0.6, lw=0.8
    )
    ax_bot.axhline(
        y=-TOLERANCE_C2_PCT, color="tab:red", linestyle=":", alpha=0.6, lw=0.8
    )
    ax_bot.set_xticks(xs)
    ax_bot.set_xticklabels(
        [
            f"({M[0]},{M[1]},{M[2]})\n{voxels_per_lambda_x[i]:.1f} v/λ"
            for i, M in enumerate(MODE_INDICES)
        ],
        fontsize=8,
    )
    ax_bot.set_ylabel(r"$c^2$ err (%)")
    ax_bot.grid(alpha=0.3, axis="y")

    fig.tight_layout()
    fig.savefig(str(plot_path), dpi=120)
    print(f"\nplot saved: {plot_path}")

    return 0 if all_pass else 1


if __name__ == "__main__":
    sys.exit(main())
