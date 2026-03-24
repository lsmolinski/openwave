"""
1D WAVE ENGINE SANDBOX v3 — Phase 1b: Base Wave Research

Base wave candidate models for testing the fundamental energy wave field:
- Model A: Uniform oscillation (null baseline)
- Model B: Standing wave (counter-propagating)
- Model C: Stochastic multi-phase (statistical isotropy)
- Model D: Dual-phase standing wave (quadrature complementary)
- Model E: Laplacian propagation (time-stepped, reflecting BC)

WC interaction: Node-locking charge hypothesis (Phase 1b Step 2).
Charge is spatial (determined by node position), not intrinsic to the WC.

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
# Wave Center Configuration (Node-Locking Hypothesis)
# ================================================================
# Hypothesis: charge is NOT intrinsic to the WC — it emerges from
# the WC's position in the standing wave node lattice.
# - WCs move to energy minima (nodes) via F = -∇E
# - Adjacent nodes are λ/2 apart
# - Even nodes vs odd nodes have opposite displacement gradient
# - This alternating character IS the emergent charge
# - WC phase (source_offset) is set to 0 for ALL WCs — no imposed charge

WC_ENABLED = True  # False = base wave only (Step 1 mode)
WC_SNAP_TO_NODES = True  # snap WC positions to nearest standing wave node

# Weight function for WC standing/traveling wave transition (from v2, eq #5)
TRANSITION_LAM = 1.25
WEIGHT_POWER = 8


class WaveCenter:
    """A wave center — charge-neutral object that acquires charge from position."""

    def __init__(self, x_am: float, amplitude: float = A0_am):
        self.x_am = x_am
        self.amplitude = amplitude
        # Phase is ZERO — charge comes from node position, not from this parameter
        self.phase = 0.0


def _snap_to_grid(pos):
    """Snap position to nearest grid point."""
    idx = np.argmin(np.abs(x_am - pos))
    return float(x_am[idx])


def _snap_to_node(pos):
    """Snap position to nearest standing wave energy node.

    Nodes of cos²(kx) are at kx = π/2 + nπ, i.e., x = (2n+1)·λ/4.
    """
    # Node positions: x_n = (2n+1) · λ/4
    n = np.round((pos / (lam_am / 4) - 1) / 2)
    x_node = (2 * n + 1) * lam_am / 4
    return _snap_to_grid(x_node)


def _get_node_index(pos):
    """Return the node index n for a WC at position pos.

    Node n is at x = (2n+1)·λ/4. Even n = one 'charge', odd n = opposite.
    """
    return int(np.round((pos / (lam_am / 4) - 1) / 2))


# Default: 2 WCs separated by 3 nodes (1.5λ) — odd separation = "opposite charge"
wc_sep_nodes = 3  # separation in node-count (λ/2 units)
wave_centers = [
    WaveCenter(x_am=_snap_to_node(-wc_sep_nodes * lam_am / 4)),
    WaveCenter(x_am=_snap_to_node(+wc_sep_nodes * lam_am / 4)),
]


def _compute_weight(r_am):
    """In-wave weight: 1 near WC (standing), 0 far away (traveling)."""
    return 1.0 / (1.0 + (r_am / (TRANSITION_LAM * lam_am)) ** WEIGHT_POWER)


# ================================================================
# Combined Phasor: Base Wave + WC Waves
# ================================================================
# Both base wave and WC waves share the same ω, so phasor superposition
# is exact. Base wave contributes P_base, Q_base. Each WC contributes
# P_n, Q_n via the weighted partial standing wave (eq #5 from v2).
#
# Total: P = P_base + Σ P_n,  Q = Q_base + Σ Q_n
# RMS = √(P² + Q²) / √2


def compute_combined_rms(x):
    """Compute phasor RMS from base standing wave + all WC waves."""
    # Base wave phasor: ψ_base = A₀·cos(kx)·cos(ωt) → P = A₀·cos(kx), Q = 0
    P = A0_am * np.cos(k * x)
    Q = np.zeros_like(x)

    if WC_ENABLED:
        for wc in wave_centers:
            r = np.abs(x - wc.x_am)
            kr = k * r
            w = _compute_weight(r)

            # Eq #5: weighted partial standing wave phasor coefficients
            with np.errstate(divide="ignore", invalid="ignore"):
                C_n = np.where(kr < 1e-10, 2.0 * wc.amplitude,
                               wc.amplitude * (w + 1.0) * np.sin(kr) / kr)
                S_n = np.where(kr < 1e-10, 0.0,
                               wc.amplitude * (w - 1.0) * np.cos(kr) / kr)

            # Rotate by phase (always 0 in node-locking, but kept general)
            cos_phi = np.cos(wc.phase)
            sin_phi = np.sin(wc.phase)
            P += C_n * cos_phi + S_n * sin_phi
            Q += -C_n * sin_phi + S_n * cos_phi

    peak = np.sqrt(P**2 + Q**2)
    return peak / np.sqrt(2.0)


def compute_combined_displacement(x, t):
    """Compute total displacement: base standing wave + WC waves."""
    psi = A0_am * np.cos(k * x) * np.cos(omega * t)

    if WC_ENABLED:
        for wc in wave_centers:
            r = np.abs(x - wc.x_am)
            kr = k * r
            w = _compute_weight(r)
            wt_phi = omega * t + wc.phase

            with np.errstate(divide="ignore", invalid="ignore"):
                wave = np.where(kr < 1e-10, 2.0 * np.cos(wt_phi),
                                (w * np.sin(kr + wt_phi) + np.sin(kr - wt_phi)) / kr)

            psi += wc.amplitude * wave

    return psi


# ================================================================
# Base Wave Mode Configuration
# ================================================================
# uniform    = Uniform Oscillation:       ψ = A₀·cos(ωt), flat energy
# standing   = Standing Wave:             ψ = A₀·cos(kx)·cos(ωt), nodes at λ/2
# stochastic = Stochastic (N-source):     N broadband random-phase waves, ~flat energy
# quadrature = Dual-Phase Standing Wave:  two offset standing waves, flat energy (at 90°)
# laplacian  = Laplacian Propagation:     time-stepped wave equation, reflecting BC

BASE_WAVE_MODE = "standing"

BASE_WAVE_NAMES = {
    "uniform": "Uniform Oscillation",
    "standing": "Standing Wave",
    "stochastic": "Stochastic (N-source)",
    "quadrature": "Dual-Phase Standing Wave",
    "laplacian": "Laplacian Propagation",
}

# --- Model C parameters ---
STOCHASTIC_N_SOURCES = 100
STOCHASTIC_SEED = 42
# Bandwidth: spread of k values around base k.
# Monochromatic (σ=0) collapses to a single phase-shifted standing wave
# because Σ cos(kx+φᵢ) = |S|·cos(kx+arg(S)) for any random φᵢ at same k.
# Broadband (σ>0) breaks this degeneracy and produces quasi-uniform energy.
STOCHASTIC_K_SPREAD = 1  # fractional spread: k_i ∈ k·(1 ± σ)

# --- Model D parameters ---
# Spatial offset between channels: π/2 → flat energy, π → nodes
DUAL_SPATIAL_OFFSET = np.pi / 2
# Temporal offset between channels: π/2 → traveling wave sum, 0 → standing wave sum
DUAL_TEMPORAL_OFFSET = np.pi / 2

# --- Model E parameters ---
LAPLACIAN_WARMUP_PERIODS = 20  # warmup in wave periods before animation
LAPLACIAN_INIT = "standing"  # "gaussian" pulse or "standing" analytical


# ================================================================
# Stochastic Model Precomputation (Model C) — Broadband
# ================================================================
# Each source has its own k_i and ω_i = c·k_i (broadband).
# This breaks the monochromatic degeneracy where Σ cos(kx+φ) collapses
# to a single standing wave. With different k_i, spatial patterns don't
# combine into a single sinusoid — producing quasi-uniform RMS.

_rng_c = np.random.default_rng(STOCHASTIC_SEED)
_stochastic_phases = _rng_c.uniform(0, 2 * np.pi, STOCHASTIC_N_SOURCES)
_stochastic_k = k * (1.0 + STOCHASTIC_K_SPREAD * _rng_c.uniform(-1, 1, STOCHASTIC_N_SOURCES))
_stochastic_omega = c_am_rs * _stochastic_k  # ω_i = c · k_i

# Pre-compute per-source spatial patterns: cos(k_i·x + φ_i)  shape (N, N_POINTS)
_stochastic_spatial = np.cos(_stochastic_k[:, None] * x_am[None, :] + _stochastic_phases[:, None])

# Amplitude scale: A₀·√2/√N so that RMS ≈ A₀/√2 (matches Model A energy)
_stochastic_amp = A0_am * np.sqrt(2.0) / np.sqrt(STOCHASTIC_N_SOURCES)

# Pre-compute RMS: different ω_i → cross terms average to zero →
# RMS²(x) = (amp²/2) · Σ cos²(k_i·x + φ_i), RMS = amp/√2 · √(Σ cos²)
_stochastic_rms = (_stochastic_amp / np.sqrt(2.0)) * np.sqrt(
    np.sum(_stochastic_spatial**2, axis=0)
)


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
    "steps_per_period": 1,
    # Live RMS tracking: peak |ψ| over current period
    "peak_tracker": np.zeros(N_POINTS),
    "peak_steps": 0,
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
    _lap["steps_per_period"] = steps_per_period
    _lap["steps_per_frame"] = max(1, steps_per_period // 50)
    _lap["peak_tracker"][:] = 0.0
    _lap["peak_steps"] = 0
    _lap["initialized"] = True


def _lap_step(n=1):
    """Advance Laplacian simulation by n timesteps using leap-frog."""
    if not _lap["initialized"]:
        _lap_init()

    psi = _lap["psi"]
    psi_old = _lap["psi_old"]
    dt = _lap["dt"]
    c2dt2 = (c_am_rs * dt) ** 2

    peak_tracker = _lap["peak_tracker"]
    spp = _lap["steps_per_period"]

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

        # Live RMS: track peak |ψ| and update RMS every full period
        np.maximum(peak_tracker, np.abs(psi), out=peak_tracker)
        _lap["peak_steps"] += 1
        if _lap["peak_steps"] >= spp:
            _lap["rms"] = peak_tracker / np.sqrt(2.0)
            peak_tracker[:] = 0.0
            _lap["peak_steps"] = 0


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
    if BASE_WAVE_MODE == "uniform":
        # Uniform: every point oscillates together, no spatial structure
        return A0_am * np.cos(omega * t) * np.ones_like(x)

    elif BASE_WAVE_MODE == "standing":
        # Standing wave: counter-propagating left + right waves
        return A0_am * np.cos(k * x) * np.cos(omega * t)

    elif BASE_WAVE_MODE == "stochastic":
        # Broadband stochastic: each source has its own k_i and ω_i
        # ψ = amp · Σ cos(k_i·x + φ_i) · cos(ω_i·t)
        temporal = np.cos(_stochastic_omega * t)  # shape (N,)
        return _stochastic_amp * np.sum(_stochastic_spatial * temporal[:, None], axis=0)

    elif BASE_WAVE_MODE == "quadrature":
        # Dual-phase: two standing wave channels with spatial + temporal offsets
        # Channel 1: A₀ · cos(kx) · cos(ωt)
        # Channel 2: A₀ · cos(kx + δs) · cos(ωt + δt)
        # Displayed displacement is the field sum (what would be measured)
        ch1 = np.cos(k * x) * np.cos(omega * t)
        ch2 = np.cos(k * x + DUAL_SPATIAL_OFFSET) * np.cos(omega * t + DUAL_TEMPORAL_OFFSET)
        return A0_am * (ch1 + ch2)

    elif BASE_WAVE_MODE == "laplacian":
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
    if BASE_WAVE_MODE == "uniform":
        # Uniform: peak = A₀ everywhere, RMS = A₀/√2
        return (A0_am / np.sqrt(2.0)) * np.ones_like(x)

    elif BASE_WAVE_MODE == "standing":
        # Standing: peak at x = A₀|cos(kx)|, RMS = peak/√2
        return (A0_am / np.sqrt(2.0)) * np.abs(np.cos(k * x))

    elif BASE_WAVE_MODE == "stochastic":
        # Broadband: different ω_i → cross terms average to zero →
        # RMS²(x) = (amp²/2) · Σ cos²(k_i·x + φ_i) ≈ A₀²/2 (quasi-uniform)
        return _stochastic_rms

    elif BASE_WAVE_MODE == "quadrature":
        # Per-channel RMS², summed (independent energy channels):
        # RMS₁² = (A₀²/2)·cos²(kx),  RMS₂² = (A₀²/2)·cos²(kx + δs)
        # RMS_combined = √(RMS₁² + RMS₂²)
        # For δs = π/2: cos² + sin² = 1 → RMS = A₀/√2 (FLAT)
        # For δs = π:   2cos² → A₀|cos(kx)| (NODES, like B but √2 larger)
        rms_sq = (A0_am**2 / 2.0) * (np.cos(k * x) ** 2 + np.cos(k * x + DUAL_SPATIAL_OFFSET) ** 2)
        return np.sqrt(rms_sq)

    elif BASE_WAVE_MODE == "laplacian":
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
        Slider: WC separation in node-steps (λ/2 units)
    """
    from matplotlib.widgets import Slider
    from matplotlib.transforms import blended_transform_factory

    state = {"paused": START_PAUSED, "frame": 0}

    def _compute_rms():
        if WC_ENABLED and wave_centers:
            return compute_combined_rms(x_am)
        return compute_base_rms(x_am)

    def _compute_disp(t):
        if WC_ENABLED and wave_centers:
            return compute_combined_displacement(x_am, t)
        return compute_base_displacement(x_am, t)

    # Initial computation
    rms = _compute_rms()
    energy = compute_energy_density(rms) * ENERGY_TO_J
    force = compute_force_field(rms) * FORCE_TO_N

    # Setup figure
    plt.style.use("dark_background")
    fig = plt.figure(figsize=(16, 9), facecolor=colormap.DARK_GRAY[1])

    gs = fig.add_gridspec(
        3, 1, height_ratios=[3, 1, 1], hspace=0.30,
        top=0.91, bottom=0.10, left=0.06, right=0.98,
    )

    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0])
    ax3 = fig.add_subplot(gs[2, 0])
    ax_slider = fig.add_axes([0.15, 0.02, 0.70, 0.025])

    mode_name = BASE_WAVE_NAMES.get(BASE_WAVE_MODE, "Unknown")
    wc_label = "Node-Locking" if WC_ENABLED else "Base Only"
    fig.suptitle(
        f"OPENWAVE — Phase 1b  [{mode_name} + {wc_label}]",
        fontsize=16, family="Monospace", y=0.98,
    )

    # Mode info
    info_parts = [f"{mode_name}"]
    if WC_ENABLED:
        n1 = _get_node_index(wave_centers[0].x_am)
        n2 = _get_node_index(wave_centers[1].x_am)
        sep_nodes = abs(n2 - n1)
        parity = "ODD (opposite)" if sep_nodes % 2 == 1 else "EVEN (same)"
        info_parts.append(f"WC sep: {sep_nodes} nodes ({sep_nodes/2:.1f}λ) = {parity}")
    fig.text(
        0.50, 0.935, "  |  ".join(info_parts),
        fontsize=10, family="Monospace", ha="center", va="center", color="#AAAAAA",
    )
    # Store reference for updates
    info_text = fig.texts[-1]

    # --- Panel 1: Displacement + RMS ---
    (line_psi,) = ax1.plot(
        [], [], color=colormap.viridis_palette[2][1],
        linewidth=1.5, alpha=0.8, label="ψ(x,t)",
    )
    (line_rms_pos,) = ax1.plot(
        x_am, rms, color=colormap.viridis_palette[4][1],
        linewidth=2, label="RMS envelope",
    )
    (line_rms_neg,) = ax1.plot(
        x_am, -rms, color=colormap.viridis_palette[4][1],
        linewidth=2, alpha=0.5,
    )

    # Node markers on panel 1 (standing wave energy zeros)
    node_positions = []
    if BASE_WAVE_MODE == "standing":
        n_max = int(domain_half / (lam_am / 2)) + 1
        for n in range(-n_max, n_max + 1):
            x_node = (2 * n + 1) * lam_am / 4
            if x_am.min() < x_node < x_am.max():
                node_positions.append(x_node)
                node_type = "even" if n % 2 == 0 else "odd"
                color = "#FF444440" if node_type == "even" else "#4488FF40"
                ax1.axvline(x=x_node, color=color, linestyle=":", alpha=0.3, linewidth=0.8)

    # WC marker lines
    wc_vlines = []
    if WC_ENABLED:
        for i, wc in enumerate(wave_centers):
            ni = _get_node_index(wc.x_am)
            node_type = "even" if ni % 2 == 0 else "odd"
            color = "#FF4444" if node_type == "even" else "#4488FF"
            label = f"WC{i+1} (node {ni}, {node_type})"
            vl1 = ax1.axvline(x=wc.x_am, color=color, linestyle="--", alpha=0.6, label=label)
            vl2 = ax2.axvline(x=wc.x_am, color=color, linestyle="--", alpha=0.4)
            vl3 = ax3.axvline(x=wc.x_am, color=color, linestyle="--", alpha=0.4)
            wc_vlines.append((vl1, vl2, vl3))

    ax1.axhline(y=0, color="w", linestyle="--", alpha=0.2)
    ax1.set_xlim(x_am.min(), x_am.max())
    psi_max = max(np.max(rms) * np.sqrt(2.0) * 1.3, A0_am * 0.1)
    ax1.set_ylim(-psi_max, psi_max)
    ax1.set_ylabel("Displacement (am)", family="Monospace")
    ax1.legend(loc="upper right", fontsize=8)
    ax1.grid(True, alpha=0.2)

    time_text = ax1.text(
        0.02, 0.95, "", transform=ax1.transAxes,
        fontsize=10, family="Monospace", verticalalignment="top",
    )

    # --- Panel 2: Energy Density ---
    (line_energy,) = ax2.plot(
        x_am, energy, color=colormap.viridis_palette[3][1],
        linewidth=1.5, label="E(x) = ρV(fA)² [J]",
    )
    fill_energy = [ax2.fill_between(x_am, energy, color=colormap.viridis_palette[3][1], alpha=0.5)]
    ax2.set_xlim(x_am.min(), x_am.max())
    ax2.set_ylabel("Energy (J)", family="Monospace", fontsize=9)
    ax2.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    ax2.legend(loc="upper right", fontsize=9)
    ax2.grid(True, alpha=0.2)
    e_max = np.max(energy) if np.max(energy) > 0 else 1e-10
    ax2.set_ylim(0, e_max * 1.2)

    energy_info = ax2.text(
        0.5, 0.90, "", transform=ax2.transAxes, fontsize=10,
        family="Monospace", fontweight="bold", ha="center", va="top", color="#AAAAAA",
    )

    # --- Panel 3: Force Field ---
    (line_force,) = ax3.plot(x_am, force, color="w", linewidth=1, alpha=0.8)
    fill_force_r = [ax3.fill_between(
        x_am, force, where=(force > 0), color="#FF4444", alpha=0.3, label="→ Right (+)",
    )]
    fill_force_l = [ax3.fill_between(
        x_am, force, where=(force < 0), color="#4488FF", alpha=0.3, label="← Left (−)",
    )]
    ax3.axhline(y=0, color="w", linestyle="--", alpha=0.3)
    ax3.set_xlim(x_am.min(), x_am.max())
    ax3.set_xlabel("X (am)", family="Monospace")
    ax3.set_ylabel("Force (N)", family="Monospace", fontsize=9)
    ax3.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    ax3.legend(loc="upper right", fontsize=9)
    ax3.grid(True, alpha=0.2)
    f_max = np.max(np.abs(force)) if np.max(np.abs(force)) > 0 else 1e-10
    ax3.set_ylim(-f_max * 1.3, f_max * 1.3)

    # Force annotations at WC positions
    wc_force_texts = []
    if WC_ENABLED:
        for wc in wave_centers:
            txt = ax3.text(
                wc.x_am, 0.97, "",
                transform=blended_transform_factory(ax3.transData, ax3.transAxes),
                fontsize=9, family="Monospace", fontweight="bold",
                ha="center", va="top",
                bbox=dict(facecolor=colormap.DARK_GRAY[1], edgecolor="none", pad=2),
            )
            wc_force_texts.append(txt)

    # Coulomb / node-locking result text
    result_text = ax3.text(
        0.5, 0.02, "", transform=ax3.transAxes, fontsize=9, family="Monospace",
        ha="center", va="bottom", color="#FFCC00",
        bbox=dict(facecolor="#222222", edgecolor="#FFCC00", boxstyle="round,pad=0.3", alpha=0.7),
    )

    # --- Separation Slider (node steps) ---
    slider_sep = Slider(
        ax_slider, "WC Separation (nodes, λ/2 each)", 1, 30,
        valinit=wc_sep_nodes, valstep=1,
        color=colormap.viridis_palette[3][1],
    )

    # --- Refresh function ---

    def refresh_static():
        """Recompute and redraw all static elements."""
        new_rms = _compute_rms()
        new_energy = compute_energy_density(new_rms) * ENERGY_TO_J
        new_force = compute_force_field(new_rms) * FORCE_TO_N

        line_rms_pos.set_ydata(new_rms)
        line_rms_neg.set_ydata(-new_rms)
        pm = max(np.max(new_rms) * np.sqrt(2.0) * 1.3, A0_am * 0.1)
        ax1.set_ylim(-pm, pm)

        line_energy.set_ydata(new_energy)
        fill_energy[0].remove()
        fill_energy[0] = ax2.fill_between(
            x_am, new_energy, color=colormap.viridis_palette[3][1], alpha=0.5,
        )
        em = np.max(new_energy) if np.max(new_energy) > 0 else 1e-10
        ax2.set_ylim(0, em * 1.2)

        line_force.set_ydata(new_force)
        fill_force_r[0].remove()
        fill_force_l[0].remove()
        fill_force_r[0] = ax3.fill_between(
            x_am, new_force, where=(new_force > 0), color="#FF4444", alpha=0.3,
        )
        fill_force_l[0] = ax3.fill_between(
            x_am, new_force, where=(new_force < 0), color="#4488FF", alpha=0.3,
        )
        fm = np.max(np.abs(new_force)) if np.max(np.abs(new_force)) > 0 else 1e-10
        ax3.set_ylim(-fm * 1.3, fm * 1.3)

        # Update WC markers, force annotations, and node-locking analysis
        if WC_ENABLED and len(wave_centers) >= 2:
            for i, (vl1, vl2, vl3) in enumerate(wc_vlines):
                wc = wave_centers[i]
                ni = _get_node_index(wc.x_am)
                node_type = "even" if ni % 2 == 0 else "odd"
                color = "#FF4444" if node_type == "even" else "#4488FF"
                for vl in (vl1, vl2, vl3):
                    vl.set_xdata([wc.x_am, wc.x_am])
                    vl.set_color(color)
                vl1.set_label(f"WC{i+1} (node {ni}, {node_type})")

            ax1.legend(loc="upper right", fontsize=8)

            # Force at each WC
            idx_left = np.argmin(np.abs(x_am - wave_centers[0].x_am))
            idx_right = np.argmin(np.abs(x_am - wave_centers[1].x_am))
            f_left = new_force[idx_left]
            f_right = new_force[idx_right]

            is_attract = f_left > 0 and f_right < 0
            is_repel = f_left < 0 and f_right > 0

            # Node analysis
            n1 = _get_node_index(wave_centers[0].x_am)
            n2 = _get_node_index(wave_centers[1].x_am)
            sep_nodes = abs(n2 - n1)
            parity = "ODD" if sep_nodes % 2 == 1 else "EVEN"
            expected = "attract" if sep_nodes % 2 == 1 else "repel"
            actual = "attract" if is_attract else "repel" if is_repel else "unclear"
            match = actual == expected

            # Update info text
            info_parts_new = [f"{mode_name}"]
            info_parts_new.append(
                f"WC sep: {sep_nodes} nodes ({sep_nodes/2:.1f}λ) = {parity}"
            )
            info_text.set_text("  |  ".join(info_parts_new))

            # Force annotations
            for i, txt in enumerate(wc_force_texts):
                wc = wave_centers[i]
                txt.set_x(wc.x_am)
                f_val = f_left if i == 0 else f_right
                if is_attract:
                    label = f">>> Attract\n{f_val:.2e} N" if i == 0 else f"Attract <<<\n{f_val:.2e} N"
                    txt.set_color("#88FF00")
                elif is_repel:
                    label = f"<<< Repel\n{f_val:.2e} N" if i == 0 else f"Repel >>>\n{f_val:.2e} N"
                    txt.set_color("#FF8800")
                else:
                    label = f"{f_val:.2e} N"
                    txt.set_color("#AAAAAA")
                txt.set_text(label)

            # Result text: does node parity predict force direction?
            match_color = "#88FF00" if match else "#FF4444"
            match_label = "MATCH" if match else "MISMATCH"
            result_text.set_text(
                f"Node-Lock: {parity} sep → expect {expected} → got {actual} → {match_label}"
            )
            result_text.set_color(match_color)
            result_text.get_bbox_patch().set_edgecolor(match_color)
            result_text.set_visible(True)
        else:
            result_text.set_visible(False)
            energy_info.set_text("")

        fig.canvas.draw_idle()

    def on_slider_change(val):
        """Update WC positions from node-step slider."""
        sep_nodes = int(val)
        wave_centers[0].x_am = _snap_to_node(-sep_nodes * lam_am / 4)
        wave_centers[1].x_am = _snap_to_node(+sep_nodes * lam_am / 4)
        refresh_static()

    slider_sep.on_changed(on_slider_change)

    # --- Animation ---

    def update_plot(frame):
        t_rs = (frame / ANIMATION_FRAMES) * period_rs
        psi = _compute_disp(t_rs)
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
        fig, animate, frames=ANIMATION_FRAMES,
        interval=ANIMATION_INTERVAL, blit=False,
    )

    # Initial static draw
    refresh_static()

    return anim, slider_sep


if __name__ == "__main__":
    anim, _slider = plot_sandbox()
    plt.show()
