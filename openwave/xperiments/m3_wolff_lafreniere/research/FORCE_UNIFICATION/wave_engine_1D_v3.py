"""
1D WAVE ENGINE SANDBOX v3 — Phase 1b: Base Wave Research

Base wave candidate models for testing the fundamental energy wave field:
- Model A: Uniform oscillation (null baseline)
- Model B: Standing wave (counter-propagating)
- Model C: Stochastic multi-phase (statistical isotropy)
- Model D: Dual-phase standing wave (quadrature complementary)
- Model E: Laplacian propagation (time-stepped, reflecting BC)

WC interaction with base wave: not yet implemented (Phase 1b Step 2).

3 panels:
  1. Displacement ψ(x,t) + RMS envelope (overlay)
  2. Energy density E(x) = ρ·V·(f·A)²
  3. Force field F(x) = -∇E

Controls:
  SPACE: Pause/Resume animation
  LEFT/RIGHT arrows: Step frame when paused
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from pathlib import Path

from openwave.common import colormap, constants


# ================================================================
# Physical Constants (EWT)
# ================================================================

A0_am = constants.EWAVE_AMPLITUDE / constants.ATTOMETER  # amplitude in am
lam_am = constants.EWAVE_LENGTH / constants.ATTOMETER  # wavelength in am
f_rHz = constants.EWAVE_FREQUENCY * constants.RONTOSECOND  # frequency in rHz
rho_qgam = constants.MEDIUM_DENSITY_QGAM  # medium density in qg/am³

omega = 2.0 * np.pi * f_rHz  # angular frequency (rad/rs)
k = 2.0 * np.pi / lam_am  # wave number (rad/am)
period_rs = 2.0 * np.pi / omega  # one full period in rs
c_am_rs = lam_am * f_rHz  # wave speed in am/rs

# SI unit conversion factors (internal → SI)
ENERGY_TO_J = constants.QUECTOGRAM * constants.ATTOMETER**2 / constants.RONTOSECOND**2
FORCE_TO_N = constants.QUECTOGRAM * constants.ATTOMETER / constants.RONTOSECOND**2


# ================================================================
# Spatial Domain
# ================================================================

domain_half = 10 * lam_am  # 10λ each side
x_am = np.linspace(-domain_half, +domain_half, 4001)  # odd for exact center
dx_am = x_am[1] - x_am[0]
N_POINTS = len(x_am)


# ================================================================
# Base Wave Mode Configuration
# ================================================================
# A = Uniform Oscillation:       ψ = A₀·cos(ωt), flat energy
# B = Standing Wave:             ψ = A₀·cos(kx)·cos(ωt), nodes at λ/2
# C = Stochastic (N-source):     N random-phase standing waves, ~flat energy
# D = Dual-Phase Standing Wave:  two offset standing waves, flat energy (at 90°)
# E = Laplacian Propagation:     time-stepped wave equation, reflecting BC

BASE_WAVE_MODE = "E"

BASE_WAVE_NAMES = {
    "A": "Uniform Oscillation",
    "B": "Standing Wave",
    "C": "Stochastic (N-source)",
    "D": "Dual-Phase Standing Wave",
    "E": "Laplacian Propagation",
}

# --- Model C parameters ---
STOCHASTIC_N_SOURCES = 100
STOCHASTIC_SEED = 42

# --- Model D parameters ---
# Spatial offset between channels: π/2 → flat energy, π → nodes
DUAL_SPATIAL_OFFSET = np.pi / 2
# Temporal offset between channels: π/2 → traveling wave sum, 0 → standing wave sum
DUAL_TEMPORAL_OFFSET = np.pi / 2

# --- Model E parameters ---
LAPLACIAN_WARMUP_PERIODS = 20  # warmup in wave periods before animation
LAPLACIAN_INIT = "gaussian"  # "gaussian" pulse or "standing" analytical


# ================================================================
# Stochastic Model Precomputation (Model C)
# ================================================================

_stochastic_phases = np.random.default_rng(STOCHASTIC_SEED).uniform(
    0, 2 * np.pi, STOCHASTIC_N_SOURCES
)

# Pre-compute spatial envelope: (2A₀/√N) · Σ cos(kx + φᵢ)
_stochastic_envelope = np.zeros(N_POINTS)
for _phi in _stochastic_phases:
    _stochastic_envelope += np.cos(k * x_am + _phi)
_stochastic_envelope *= 2.0 * A0_am / np.sqrt(STOCHASTIC_N_SOURCES)


# ================================================================
# Laplacian Model State (Model E)
# ================================================================

_lap = {
    "psi": np.zeros(N_POINTS),
    "psi_old": np.zeros(N_POINTS),
    "dt": 0.0,
    "t": 0.0,
    "rms": None,
    "initialized": False,
    "steps_per_frame": 1,
}


def _lap_init():
    """Initialize Laplacian field with Gaussian pulse or analytical standing wave."""
    dt = 0.9 * dx_am / c_am_rs  # CFL condition
    _lap["dt"] = dt
    _lap["t"] = 0.0

    if LAPLACIAN_INIT == "gaussian":
        # Gaussian pulse at center — waves propagate outward, reflect off boundaries
        sigma = 2.0 * lam_am
        _lap["psi"] = A0_am * np.exp(-(x_am**2) / (2.0 * sigma**2))
        _lap["psi_old"] = _lap["psi"].copy()  # zero initial velocity
    else:
        # Analytical standing wave — should be a stable solution of the wave equation
        _lap["psi"] = A0_am * np.cos(k * x_am)
        _lap["psi_old"] = A0_am * np.cos(k * x_am) * np.cos(-omega * dt)

    # Dirichlet BC: ψ = 0 at boundaries
    _lap["psi"][0] = _lap["psi"][-1] = 0.0
    _lap["psi_old"][0] = _lap["psi_old"][-1] = 0.0

    steps_per_period = max(1, int(period_rs / dt))
    _lap["steps_per_frame"] = max(1, steps_per_period // 50)
    _lap["initialized"] = True


def _lap_step(n=1):
    """Advance Laplacian simulation by n timesteps using leap-frog."""
    if not _lap["initialized"]:
        _lap_init()

    psi = _lap["psi"]
    psi_old = _lap["psi_old"]
    dt = _lap["dt"]
    c2dt2 = (c_am_rs * dt) ** 2

    for _ in range(n):
        # 3-point Laplacian on interior
        laplacian = np.empty_like(psi)
        laplacian[0] = laplacian[-1] = 0.0
        laplacian[1:-1] = (psi[2:] - 2.0 * psi[1:-1] + psi[:-2]) / dx_am**2

        psi_new = 2.0 * psi - psi_old + c2dt2 * laplacian

        # Dirichlet BC
        psi_new[0] = 0.0
        psi_new[-1] = 0.0

        psi_old[:] = psi
        psi[:] = psi_new
        _lap["t"] += dt


def _lap_warmup():
    """Run warmup steps, then measure RMS over one full period."""
    if not _lap["initialized"]:
        _lap_init()

    dt = _lap["dt"]
    steps_per_period = max(1, int(period_rs / dt))
    warmup_steps = LAPLACIAN_WARMUP_PERIODS * steps_per_period

    print(f"Laplacian warmup: {warmup_steps} steps ({LAPLACIAN_WARMUP_PERIODS} periods)...")
    _lap_step(warmup_steps)

    # Measure RMS: track peak |ψ| at each point over one full period
    peak = np.zeros(N_POINTS)
    for _ in range(steps_per_period):
        _lap_step(1)
        np.maximum(peak, np.abs(_lap["psi"]), out=peak)

    _lap["rms"] = peak / np.sqrt(2.0)
    print(f"Laplacian warmup complete. t = {_lap['t']:.2f} rs")


# ================================================================
# Base Wave Computation
# ================================================================


def compute_base_displacement(x, t):
    """Compute base wave displacement ψ(x,t) for current mode.

    Models A–D are analytical; Model E returns current Laplacian field state.
    """
    if BASE_WAVE_MODE == "A":
        # Uniform: every point oscillates together, no spatial structure
        return A0_am * np.cos(omega * t) * np.ones_like(x)

    elif BASE_WAVE_MODE == "B":
        # Standing wave: counter-propagating left + right waves
        return A0_am * np.cos(k * x) * np.cos(omega * t)

    elif BASE_WAVE_MODE == "C":
        # Stochastic: pre-computed spatial envelope × cos(ωt)
        return _stochastic_envelope * np.cos(omega * t)

    elif BASE_WAVE_MODE == "D":
        # Dual-phase: two standing wave channels with spatial + temporal offsets
        # Channel 1: A₀ · cos(kx) · cos(ωt)
        # Channel 2: A₀ · cos(kx + δs) · cos(ωt + δt)
        # Displayed displacement is the field sum (what would be measured)
        ch1 = np.cos(k * x) * np.cos(omega * t)
        ch2 = np.cos(k * x + DUAL_SPATIAL_OFFSET) * np.cos(omega * t + DUAL_TEMPORAL_OFFSET)
        return A0_am * (ch1 + ch2)

    elif BASE_WAVE_MODE == "E":
        # Laplacian: return current simulation state
        if not _lap["initialized"]:
            _lap_warmup()
        return _lap["psi"].copy()

    return np.zeros_like(x)


def compute_base_rms(x):
    """Compute base wave RMS envelope for current mode.

    For dual-phase (D): per-channel energy sum gives combined RMS.
    Energy is independent per channel — RMS_combined = √(RMS₁² + RMS₂²).
    For Laplacian (E): measured empirically from simulation.
    """
    if BASE_WAVE_MODE == "A":
        # Uniform: peak = A₀ everywhere, RMS = A₀/√2
        return (A0_am / np.sqrt(2.0)) * np.ones_like(x)

    elif BASE_WAVE_MODE == "B":
        # Standing: peak at x = A₀|cos(kx)|, RMS = peak/√2
        return (A0_am / np.sqrt(2.0)) * np.abs(np.cos(k * x))

    elif BASE_WAVE_MODE == "C":
        # Stochastic: spatial envelope varies, RMS = |envelope|/√2
        return np.abs(_stochastic_envelope) / np.sqrt(2.0)

    elif BASE_WAVE_MODE == "D":
        # Per-channel RMS², summed (independent energy channels):
        # RMS₁² = (A₀²/2)·cos²(kx),  RMS₂² = (A₀²/2)·cos²(kx + δs)
        # RMS_combined = √(RMS₁² + RMS₂²)
        # For δs = π/2: cos² + sin² = 1 → RMS = A₀/√2 (FLAT)
        # For δs = π:   2cos² → A₀|cos(kx)| (NODES, like B but √2 larger)
        rms_sq = (A0_am**2 / 2.0) * (np.cos(k * x) ** 2 + np.cos(k * x + DUAL_SPATIAL_OFFSET) ** 2)
        return np.sqrt(rms_sq)

    elif BASE_WAVE_MODE == "E":
        # Laplacian: measured RMS from warmup
        if _lap["rms"] is None:
            _lap_warmup()
        return _lap["rms"].copy()

    return np.zeros_like(x)


# ================================================================
# Energy Density and Force
# ================================================================
# E = ρ · V · (f · A)²  where V = dx³ (voxel volume)
# F = -∇E  (negative gradient of energy density)


def compute_energy_density(rms_am):
    """Energy density at each x: E = ρ · dx³ · (f · A_rms)²."""
    return rho_qgam * dx_am**3 * (f_rHz * rms_am) ** 2


def compute_force_field(rms_am):
    """Force field: F = -∇E via central differences."""
    energy = compute_energy_density(rms_am)
    return -np.gradient(energy, dx_am)


# ================================================================
# Animation and Plotting
# ================================================================

ANIMATION_FRAMES = 100
ANIMATION_INTERVAL = 50  # ms (20 FPS)
START_PAUSED = False

_MODULE_DIR = Path(__file__).parent
PLOT_DIR = _MODULE_DIR / "_plots"


def plot_sandbox():
    """Animate the 1D base wave sandbox with 3 panels.

    Panel 1: Displacement ψ(x,t) + RMS envelope
    Panel 2: Energy density E(x)
    Panel 3: Force field F(x)

    Controls:
        SPACE: Pause/Resume
        LEFT/RIGHT: Step frame when paused
    """
    state = {"paused": START_PAUSED, "frame": 0}

    # Initial computation
    rms = compute_base_rms(x_am)
    energy = compute_energy_density(rms) * ENERGY_TO_J
    force = compute_force_field(rms) * FORCE_TO_N

    # Setup figure
    plt.style.use("dark_background")
    fig = plt.figure(figsize=(16, 9), facecolor=colormap.DARK_GRAY[1])

    gs = fig.add_gridspec(
        3,
        1,
        height_ratios=[3, 1, 1],
        hspace=0.30,
        top=0.91,
        bottom=0.06,
        left=0.06,
        right=0.98,
    )

    ax1 = fig.add_subplot(gs[0, 0])  # displacement + RMS
    ax2 = fig.add_subplot(gs[1, 0])  # energy
    ax3 = fig.add_subplot(gs[2, 0])  # force

    mode_name = BASE_WAVE_NAMES.get(BASE_WAVE_MODE, "Unknown")
    fig.suptitle(
        f"OPENWAVE — Phase 1b: Base Wave  [{mode_name}]",
        fontsize=16,
        family="Monospace",
        y=0.98,
    )

    # Mode info text
    info_parts = [f"Mode {BASE_WAVE_MODE}: {mode_name}"]
    if BASE_WAVE_MODE == "C":
        info_parts.append(f"N={STOCHASTIC_N_SOURCES}, seed={STOCHASTIC_SEED}")
    elif BASE_WAVE_MODE == "D":
        info_parts.append(
            f"spatial={np.degrees(DUAL_SPATIAL_OFFSET):.0f}°, "
            f"temporal={np.degrees(DUAL_TEMPORAL_OFFSET):.0f}°"
        )
    elif BASE_WAVE_MODE == "E":
        info_parts.append(f"init={LAPLACIAN_INIT}, warmup={LAPLACIAN_WARMUP_PERIODS}T")
    info_str = "  |  ".join(info_parts)
    fig.text(
        0.50,
        0.935,
        info_str,
        fontsize=10,
        family="Monospace",
        ha="center",
        va="center",
        color="#AAAAAA",
    )

    # --- Panel 1: Displacement + RMS ---
    (line_psi,) = ax1.plot(
        [],
        [],
        color=colormap.viridis_palette[2][1],
        linewidth=1.5,
        alpha=0.8,
        label="ψ(x,t)",
    )
    (line_rms_pos,) = ax1.plot(
        x_am,
        rms,
        color=colormap.viridis_palette[4][1],
        linewidth=2,
        label="RMS envelope",
    )
    (line_rms_neg,) = ax1.plot(
        x_am,
        -rms,
        color=colormap.viridis_palette[4][1],
        linewidth=2,
        alpha=0.5,
    )
    ax1.axhline(y=0, color="w", linestyle="--", alpha=0.2)
    ax1.set_xlim(x_am.min(), x_am.max())
    psi_max = max(np.max(rms) * np.sqrt(2.0) * 1.3, A0_am * 0.1)
    ax1.set_ylim(-psi_max, psi_max)
    ax1.set_ylabel("Displacement (am)", family="Monospace")
    ax1.legend(loc="upper right", fontsize=9)
    ax1.grid(True, alpha=0.2)

    time_text = ax1.text(
        0.02,
        0.95,
        "",
        transform=ax1.transAxes,
        fontsize=10,
        family="Monospace",
        verticalalignment="top",
    )

    # --- Panel 2: Energy Density ---
    (line_energy,) = ax2.plot(
        x_am,
        energy,
        color=colormap.viridis_palette[3][1],
        linewidth=1.5,
        label="E(x) = ρV(fA)² [J]",
    )
    fill_energy = [ax2.fill_between(x_am, energy, color=colormap.viridis_palette[3][1], alpha=0.5)]
    ax2.set_xlim(x_am.min(), x_am.max())
    ax2.set_ylabel("Energy (J)", family="Monospace", fontsize=9)
    ax2.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    ax2.legend(loc="upper right", fontsize=9)
    ax2.grid(True, alpha=0.2)

    e_max = np.max(energy) if np.max(energy) > 0 else 1e-10
    ax2.set_ylim(0, e_max * 1.2)

    # Energy uniformity annotation
    e_min_val = np.min(energy)
    e_max_val = np.max(energy)
    if e_max_val > 0:
        flatness = e_min_val / e_max_val
        flat_color = "#88FF00" if flatness > 0.95 else "#FF8800" if flatness > 0.5 else "#FF4444"
        flat_label = "FLAT" if flatness > 0.95 else f"ratio min/max = {flatness:.3f}"
    else:
        flat_color = "#666666"
        flat_label = "zero energy"
    energy_info = ax2.text(
        0.5,
        0.90,
        f"Energy: {flat_label}",
        transform=ax2.transAxes,
        fontsize=10,
        family="Monospace",
        fontweight="bold",
        ha="center",
        va="top",
        color=flat_color,
    )

    # --- Panel 3: Force Field ---
    (line_force,) = ax3.plot(x_am, force, color="w", linewidth=1, alpha=0.8)
    fill_force_r = [
        ax3.fill_between(
            x_am,
            force,
            where=(force > 0),
            color="#FF4444",
            alpha=0.3,
            label="→ Right (+)",
        )
    ]
    fill_force_l = [
        ax3.fill_between(
            x_am,
            force,
            where=(force < 0),
            color="#4488FF",
            alpha=0.3,
            label="← Left (−)",
        )
    ]
    ax3.axhline(y=0, color="w", linestyle="--", alpha=0.3)
    ax3.set_xlim(x_am.min(), x_am.max())
    ax3.set_xlabel("X (am)", family="Monospace")
    ax3.set_ylabel("Force (N)", family="Monospace", fontsize=9)
    ax3.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    ax3.legend(loc="upper right", fontsize=9)
    ax3.grid(True, alpha=0.2)

    f_max = np.max(np.abs(force)) if np.max(np.abs(force)) > 0 else 1e-10
    ax3.set_ylim(-f_max * 1.3, f_max * 1.3)

    # Force magnitude annotation
    f_peak = np.max(np.abs(force))
    force_info = ax3.text(
        0.5,
        0.90,
        f"|F|_max = {f_peak:.3e} N",
        transform=ax3.transAxes,
        fontsize=10,
        family="Monospace",
        ha="center",
        va="top",
        color="#AAAAAA",
    )

    # --- Animation ---

    def update_plot(frame):
        t_rs = (frame / ANIMATION_FRAMES) * period_rs

        if BASE_WAVE_MODE == "E":
            # Laplacian: step simulation forward
            _lap_step(_lap["steps_per_frame"])
            psi = _lap["psi"].copy()
            t_rs = _lap["t"]
        else:
            psi = compute_base_displacement(x_am, t_rs)

        line_psi.set_data(x_am, psi)

        pause_str = " [PAUSED]" if state["paused"] else ""
        time_text.set_text(f"t = {t_rs:.4f} rs  (frame {frame}/{ANIMATION_FRAMES}){pause_str}")

    def animate(frame):
        if not state["paused"]:
            state["frame"] = frame
        update_plot(state["frame"])

    def on_key(event):
        if event.key == " ":
            state["paused"] = not state["paused"]
            update_plot(state["frame"])
            fig.canvas.draw_idle()
        elif event.key == "right" and state["paused"]:
            state["frame"] = (state["frame"] + 1) % ANIMATION_FRAMES
            update_plot(state["frame"])
            fig.canvas.draw_idle()
        elif event.key == "left" and state["paused"]:
            state["frame"] = (state["frame"] - 1) % ANIMATION_FRAMES
            update_plot(state["frame"])
            fig.canvas.draw_idle()

    fig.canvas.mpl_connect("key_press_event", on_key)

    anim = FuncAnimation(
        fig,
        animate,
        frames=ANIMATION_FRAMES,
        interval=ANIMATION_INTERVAL,
        blit=False,
    )

    return anim


if __name__ == "__main__":
    anim = plot_sandbox()
    plt.show()
