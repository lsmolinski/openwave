"""
1D WAVE ENGINE SANDBOX v3 — Phase 1b: Base Wave + WC Disturbance

Base wave models + wave center disturbance for testing force emergence:
- Base wave modes: uniform, standing, quadrature (+ stochastic, laplacian)
- WC disturbance: weighted partial standing wave (eq #5 from v2)
- Combined phasor: base wave + WC waves → total RMS → energy → force

3 panels:
  1. Displacement ψ(x,t) + Phasor RMS envelope (overlay)
  2. Energy density E(x) = ρ·V·(f·A)²
  3. Force field F(x) = -∇E

Controls:
  SPACE: Pause/Resume animation
  LEFT/RIGHT arrows: Step frame when paused
  Slider: WC separation (λ)
  Click: Phase toggle (opposite/same charge), WC on/off
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
# 1 qg·am²/rs² = 1e-33 kg · (1e-18 m)² / (1e-27 s)² = 1e-15 J
# 1 qg·am/rs²  = 1e-33 kg · (1e-18 m) / (1e-27 s)²  = 1e+3  N
ENERGY_TO_J = constants.QUECTOGRAM * constants.ATTOMETER**2 / constants.RONTOSECOND**2
FORCE_TO_N = constants.QUECTOGRAM * constants.ATTOMETER / constants.RONTOSECOND**2


# ================================================================
# Spatial Domain
# ================================================================

domain_half = 10 * lam_am  # 10λ each side
x_am = np.linspace(-domain_half, +domain_half, 4001)  # odd for exact center
dx_am = x_am[1] - x_am[0]
N_POINTS = len(x_am)


def snap_to_grid(pos):
    """Snap a position to the nearest grid point to avoid asymmetric gradient artifacts."""
    idx = np.argmin(np.abs(x_am - pos))
    return float(x_am[idx])


# ================================================================
# Wave Center Configuration
# ================================================================
# WC wave equation: #5 Weighted Partial Standing Wave (from v2)
# ψ_wc = A · [w(r)·sin(kr + ωt + φ) + sin(kr - ωt - φ)] / kr
# w(r) = 1/(1 + (r/transition·λ)^8) — sharp Lorentzian rolloff
# Standing waves near WC center (w ≈ 1), traveling waves far out (w → 0)
# Phase φ (source_offset) determines charge: 0 = positron, π = electron

TRANSITION_LAM = 1.25  # standing wave region in wavelengths
WEIGHT_POWER = 8  # sharpness of standing → traveling transition
SEPARATION_STEP = 0.5  # step size for separation slider (in wavelength)
DISP_SCALE = 0.5  # y-axis scale for displacement


class WaveCenter:
    """A single wave center with position, phase offset, and amplitude."""

    def __init__(self, x_am: float, phase: float = 0.0, amplitude: float = A0_am):
        self.x_am = x_am  # position on x-axis (am)
        self.phase = phase  # source_offset: 0 = positron, π = electron
        self.amplitude = amplitude  # base amplitude (am)


# Default: 2 opposite-charge WCs separated by 4λ
wc_separation = 4 * lam_am
wave_centers = [
    WaveCenter(x_am=snap_to_grid(-wc_separation / 2), phase=0.0),  # positron
    WaveCenter(x_am=snap_to_grid(+wc_separation / 2), phase=np.pi),  # electron
]


def _compute_weight(r_am):
    """In-wave weight: 1 near WC (standing), 0 far away (traveling)."""
    return 1.0 / (1.0 + (r_am / (TRANSITION_LAM * lam_am)) ** WEIGHT_POWER)


# ================================================================
# Base Wave Mode Configuration
# ================================================================
# uniform    = Uniform Oscillation:       ψ = A₀·cos(ωt), flat energy
# standing   = Standing Wave:             ψ = A₀·cos(kx)·cos(ωt), nodes at λ/2
# stochastic = Stochastic (N-source):     N broadband random-phase waves, ~flat energy
# quadrature = Quadrature Wave:  two offset standing waves, flat energy (at 90°)
# laplacian  = Laplacian Propagation:     time-stepped wave equation, reflecting BC

BASE_WAVE_MODE = "quadrature"

BASE_WAVE_NAMES = {
    "uniform": "Uniform Oscillation",
    "standing": "Standing Wave",
    "quadrature": "Quadrature Wave",
    "stochastic": "Stochastic (N-source)",
    "laplacian": "Laplacian Propagation",
}

# --- Stochastic parameters ---
STOCHASTIC_N_SOURCES = 100
STOCHASTIC_SEED = 42
# Bandwidth: spread of k values around base k.
# Monochromatic (σ=0) collapses to a single phase-shifted standing wave
# because Σ cos(kx+φᵢ) = |S|·cos(kx+arg(S)) for any random φᵢ at same k.
# Broadband (σ>0) breaks this degeneracy and produces quasi-uniform energy.
STOCHASTIC_K_SPREAD = 1  # fractional spread: k_i ∈ k·(1 ± σ)

# --- Quadrature parameters ---
# Spatial offset between channels: π/2 → flat energy, π → nodes
DUAL_SPATIAL_OFFSET = np.pi / 2
# Temporal offset between channels: π/2 → traveling wave sum, 0 → standing wave sum
DUAL_TEMPORAL_OFFSET = np.pi / 2

# --- Laplacian parameters ---
LAPLACIAN_WARMUP_PERIODS = 20  # warmup in wave periods before animation
LAPLACIAN_INIT = "standing"  # "gaussian" pulse or "standing" analytical


# ================================================================
# Stochastic Model Precomputation — Broadband
# ================================================================
# Each source has its own k_i and ω_i = c·k_i (broadband).
# This breaks the monochromatic degeneracy where Σ cos(kx+φ) collapses
# to a single standing wave. With different k_i, spatial patterns don't
# combine into a single sinusoid — producing quasi-uniform RMS.

_rng_c = np.random.default_rng(STOCHASTIC_SEED)
_stochastic_phases = _rng_c.uniform(0, 2 * np.pi, STOCHASTIC_N_SOURCES)
_stochastic_k = k * (1.0 + STOCHASTIC_K_SPREAD * _rng_c.uniform(-1, 1, STOCHASTIC_N_SOURCES))
_stochastic_omega = c_am_rs * _stochastic_k
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
# Laplacian Model State
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
    spp = max(1, int(period_rs / dt))
    _lap["steps_per_period"] = spp
    _lap["steps_per_frame"] = max(1, spp // 50)
    _lap["peak_tracker"][:] = 0.0
    _lap["peak_steps"] = 0
    _lap["initialized"] = True


def _lap_step(n=1):
    """Advance Laplacian simulation by n timesteps using leap-frog."""
    if not _lap["initialized"]:
        _lap_init()
    psi, psi_old, dt = _lap["psi"], _lap["psi_old"], _lap["dt"]
    c2dt2 = (c_am_rs * dt) ** 2
    pt, spp = _lap["peak_tracker"], _lap["steps_per_period"]
    for _ in range(n):
        # 3-point Laplacian on interior
        lap = np.empty_like(psi)
        lap[0] = lap[-1] = 0.0
        lap[1:-1] = (psi[2:] - 2.0 * psi[1:-1] + psi[:-2]) / dx_am**2
        psi_new = 2.0 * psi - psi_old + c2dt2 * lap
        psi_new[0] = psi_new[-1] = 0.0
        psi_old[:] = psi
        psi[:] = psi_new
        _lap["t"] += dt
        # Live RMS: track peak |ψ| and update RMS every full period
        np.maximum(pt, np.abs(psi), out=pt)
        _lap["peak_steps"] += 1
        if _lap["peak_steps"] >= spp:
            _lap["rms"] = pt / np.sqrt(2.0)
            pt[:] = 0.0
            _lap["peak_steps"] = 0


def _lap_warmup():
    """Warm up Laplacian simulation to reach steady-state RMS."""
    if not _lap["initialized"]:
        _lap_init()
    dt = _lap["dt"]
    # Measure RMS: track peak |ψ| at each point over one full period
    spp = max(1, int(period_rs / dt))
    _lap_step(LAPLACIAN_WARMUP_PERIODS * spp)
    peak = np.zeros(N_POINTS)
    for _ in range(spp):
        _lap_step(1)
        np.maximum(peak, np.abs(_lap["psi"]), out=peak)
    _lap["rms"] = peak / np.sqrt(2.0)


# ================================================================
# Base Wave Phasor Coefficients
# ================================================================
# Returns (P_base, Q_base) — the cos(ωt) and sin(ωt) coefficients
# of the base wave at each grid point. Used for combined phasor with WCs.


def _base_phasor(x):
    """Base wave phasor coefficients (P, Q) for current mode."""
    if BASE_WAVE_MODE == "uniform":
        return A0_am * np.ones_like(x), np.zeros_like(x)

    elif BASE_WAVE_MODE == "standing":
        return A0_am * np.cos(k * x), np.zeros_like(x)

    elif BASE_WAVE_MODE == "quadrature":
        # Channel 1: A₀·cos(kx)·cos(ωt) → P1 = A₀·cos(kx), Q1 = 0
        # Channel 2: A₀·cos(kx+δs)·cos(ωt+δt)
        #   = A₀·cos(kx+δs)·[cos(δt)·cos(ωt) - sin(δt)·sin(ωt)]
        #   → P2 = A₀·cos(kx+δs)·cos(δt), Q2 = -A₀·cos(kx+δs)·sin(δt)
        P = A0_am * np.cos(k * x) + A0_am * np.cos(k * x + DUAL_SPATIAL_OFFSET) * np.cos(
            DUAL_TEMPORAL_OFFSET
        )
        Q = -A0_am * np.cos(k * x + DUAL_SPATIAL_OFFSET) * np.sin(DUAL_TEMPORAL_OFFSET)
        return P, Q

    # For stochastic/laplacian, no clean phasor — fall back to base-only RMS
    return np.zeros_like(x), np.zeros_like(x)


# ================================================================
# Combined Computation: Base Wave + WC Waves
# ================================================================
# Both base wave and WC waves share the same ω, so phasor superposition
# is exact. Total phasor: P = P_base + Σ P_n,  Q = Q_base + Σ Q_n
# RMS = √(P² + Q²) / √2
#
# Per WC n: compute C_n, S_n from the weighted partial standing wave,
# then rotate by source_offset φ to shared cos(ωt)/sin(ωt) basis:
#   P += C_n·cos(φ) + S_n·sin(φ)
#   Q += -C_n·sin(φ) + S_n·cos(φ)


def compute_combined_rms(x, active_wcs=None):
    """Compute phasor RMS: base wave + WC waves (eq #5) via phasor superposition.

    Returns RMS = √(P² + Q²) / √2, where P and Q are the accumulated
    cos(ωt) and sin(ωt) coefficients from the base wave + all active WCs.

    For modes with clean phasors (uniform, standing, quadrature): exact analytical.
    For other modes: base RMS only (WC phasors not addable to non-monochromatic base).
    """
    if active_wcs is None:
        active_wcs = wave_centers

    # Check if mode supports phasor superposition with WCs
    if BASE_WAVE_MODE in ("uniform", "standing", "quadrature"):
        P, Q = _base_phasor(x)

        for wc in active_wcs:
            r = np.abs(x - wc.x_am)
            kr = k * r
            w = _compute_weight(r)

            # Phasor coefficients for eq #5 weighted partial standing wave:
            # C_n = (w+1)·sin(kr)/kr  (cos ωt coefficient)
            # S_n = (w-1)·cos(kr)/kr  (sin ωt coefficient)
            with np.errstate(divide="ignore", invalid="ignore"):
                C_n = np.where(
                    kr < 1e-10, 2.0 * wc.amplitude, wc.amplitude * (w + 1.0) * np.sin(kr) / kr
                )
                S_n = np.where(kr < 1e-10, 0.0, wc.amplitude * (w - 1.0) * np.cos(kr) / kr)

            # Rotate by source_offset φ to shared cos(ωt)/sin(ωt) basis
            cos_phi = np.cos(wc.phase)
            sin_phi = np.sin(wc.phase)
            P += C_n * cos_phi + S_n * sin_phi
            Q += -C_n * sin_phi + S_n * cos_phi

        return np.sqrt(P**2 + Q**2) / np.sqrt(2.0)

    # Fallback: base-only RMS for non-phasor modes
    return compute_base_rms(x)


def compute_combined_displacement(x, t, active_wcs=None):
    """Compute total displacement: base wave + WC waves (eq #5).

    ψ_total(x,t) = ψ_base(x,t) + Σ ψ_wc_n(x,t)
    where each WC contributes a weighted partial standing wave centered at its position.
    """
    if active_wcs is None:
        active_wcs = wave_centers

    psi = compute_base_displacement(x, t)

    for wc in active_wcs:
        r = np.abs(x - wc.x_am)
        kr = k * r
        w = _compute_weight(r)
        wt_phi = omega * t + wc.phase

        with np.errstate(divide="ignore", invalid="ignore"):
            wave = np.where(
                kr < 1e-10,
                2.0 * np.cos(wt_phi),
                (w * np.sin(kr + wt_phi) + np.sin(kr - wt_phi)) / kr,
            )

        psi += wc.amplitude * wave

    return psi


def compute_base_displacement(x, t):
    """Compute base wave displacement ψ(x,t) for current mode.

    Models A–D are analytical; Model E returns current Laplacian field state.
    This is the base wave ONLY — no WC contribution.
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
        temporal = np.cos(_stochastic_omega * t)
        return _stochastic_amp * np.sum(_stochastic_spatial * temporal[:, None], axis=0)
    elif BASE_WAVE_MODE == "quadrature":
        # Quadrature: two standing wave channels with spatial + temporal offsets
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
    """Compute base wave RMS envelope for current mode (no WCs).

    For Quadrature: per-channel energy sum gives combined RMS.
    Energy is independent per channel — RMS_combined = √(RMS₁² + RMS₂²).
    For Laplacian: measured empirically from simulation.
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
#
# Computing F from ∇E directly (not expanded A·∇A) so that future
# additions of variable ρ(x), f(x), λ(x) are automatically captured
# without changing the force computation logic.


def compute_energy_density(rms_am):
    """Energy density at each x: E = ρ · dx³ · (f · A_rms)²."""
    return rho_qgam * dx_am**3 * (f_rHz * rms_am) ** 2


def compute_force_field(rms_am):
    """Force field: F = -∇E.

    Computes energy per grid point, then takes the negative gradient.
    Uses central differences via np.gradient.
    """
    return -np.gradient(compute_energy_density(rms_am), dx_am)


# ================================================================
# Animation and Plotting
# ================================================================

ANIMATION_FRAMES = 100
ANIMATION_INTERVAL = 50  # ms (20 FPS)
START_PAUSED = False

_MODULE_DIR = Path(__file__).parent
PLOT_DIR = _MODULE_DIR / "_plots"


def plot_sandbox():
    """Animate the 1D wave sandbox with base wave + WC disturbance.

    Panel 1: Displacement ψ(x,t) + RMS envelope
    Panel 2: Energy density E(x)
    Panel 3: Force field F(x) with Coulomb comparison

    Controls:
        SPACE: Pause/Resume
        LEFT/RIGHT: Step frame when paused
        Slider: WC separation (λ)
        Click: Phase toggle, WC toggles
    """
    from matplotlib.widgets import Slider
    from matplotlib.transforms import blended_transform_factory

    state = {"paused": START_PAUSED, "frame": 0}

    # Track active WCs and base wave
    wc_active = [True] * len(wave_centers)
    base_active = [True]  # mutable list so inner functions can modify

    def get_active_centers():
        return [wc for wc, active in zip(wave_centers, wc_active) if active]

    def recompute_static():
        active = get_active_centers()
        if base_active[0] and active:
            rms = compute_combined_rms(x_am, active)
        elif base_active[0]:
            rms = compute_base_rms(x_am)
        elif active:
            # WCs only, no base wave — use v2-style phasor (zero base phasor)
            P = np.zeros_like(x_am)
            Q = np.zeros_like(x_am)
            for wc in active:
                r = np.abs(x_am - wc.x_am)
                kr = k * r
                w = _compute_weight(r)
                with np.errstate(divide="ignore", invalid="ignore"):
                    C_n = np.where(
                        kr < 1e-10, 2.0 * wc.amplitude, wc.amplitude * (w + 1.0) * np.sin(kr) / kr
                    )
                    S_n = np.where(kr < 1e-10, 0.0, wc.amplitude * (w - 1.0) * np.cos(kr) / kr)
                cos_phi = np.cos(wc.phase)
                sin_phi = np.sin(wc.phase)
                P += C_n * cos_phi + S_n * sin_phi
                Q += -C_n * sin_phi + S_n * cos_phi
            rms = np.sqrt(P**2 + Q**2) / np.sqrt(2.0)
        else:
            rms = np.zeros_like(x_am)
        energy = compute_energy_density(rms)
        force = compute_force_field(rms)
        return rms, energy, force

    # Initial computation
    rms, energy, force = recompute_static()
    energy_j = energy * ENERGY_TO_J
    force_n = force * FORCE_TO_N

    # Setup figure
    plt.style.use("dark_background")
    fig = plt.figure(figsize=(16, 9), facecolor=colormap.DARK_GRAY[1])
    gs = fig.add_gridspec(
        3,
        1,
        height_ratios=[3, 1, 1],
        hspace=0.30,
        top=0.91,
        bottom=0.10,
        left=0.06,
        right=0.98,
    )
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0])
    ax3 = fig.add_subplot(gs[2, 0])
    ax_slider = fig.add_axes([0.15, 0.02, 0.70, 0.025])

    mode_name = BASE_WAVE_NAMES.get(BASE_WAVE_MODE, "Unknown")
    fig.suptitle(
        f"OPENWAVE — Phase 1b: Base Wave + WC  [{mode_name}]",
        fontsize=16,
        family="Monospace",
        y=0.98,
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
        label="Phasor RMS",
    )
    (line_rms_neg,) = ax1.plot(
        x_am,
        -rms,
        color=colormap.viridis_palette[4][1],
        linewidth=2,
        alpha=0.5,
    )

    # WC marker lines
    wc_vlines = []
    for i, wc in enumerate(wave_centers):
        color = "#FF4444" if wc.phase == 0.0 else "#4488FF"
        charge = "+" if wc.phase == 0.0 else "−"
        vl1 = ax1.axvline(
            x=wc.x_am, color=color, linestyle="--", alpha=0.6, label=f"WC{i+1} ({charge})"
        )
        vl2 = ax2.axvline(x=wc.x_am, color=color, linestyle="--", alpha=0.4)
        vl3 = ax3.axvline(x=wc.x_am, color=color, linestyle="--", alpha=0.4)
        wc_vlines.append((vl1, vl2, vl3))

    ax1.axhline(y=0, color="w", linestyle="--", alpha=0.2)
    ax1.set_xlim(x_am.min(), x_am.max())
    psi_max = max(np.max(rms) * np.sqrt(2.0) * 1.3, A0_am * 0.1) * DISP_SCALE
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
        energy_j,
        color=colormap.viridis_palette[3][1],
        linewidth=1.5,
        label="E(x) = ρV(fA)² [J]",
    )
    fill_energy = [
        ax2.fill_between(x_am, energy_j, color=colormap.viridis_palette[3][1], alpha=0.5)
    ]
    ax2.set_xlim(x_am.min(), x_am.max())
    ax2.set_ylabel("Energy (J)", family="Monospace", fontsize=9)
    ax2.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    ax2.legend(loc="upper right", fontsize=9)
    ax2.grid(True, alpha=0.2)
    e_max = np.max(energy_j) if np.max(energy_j) > 0 else 1e-10
    ax2.set_ylim(0, e_max * 1.2)

    # Total energy annotation (below legend) — for energy conservation inspection
    total_e = np.sum(energy_j)
    total_energy_text = ax2.text(
        0.99,
        1.15,
        f"ΣE = {total_e:.3e} J",
        transform=ax2.transAxes,
        fontsize=9,
        family="Monospace",
        ha="right",
        va="top",
        color="#FFCC00",
    )

    # --- Panel 3: Force Field ---
    (line_force,) = ax3.plot(x_am, force_n, color="w", linewidth=1, alpha=0.8)
    fill_force_r = [
        ax3.fill_between(
            x_am, force_n, where=(force_n > 0), color="#FF4444", alpha=0.3, label="→ Right (+)"
        )
    ]
    fill_force_l = [
        ax3.fill_between(
            x_am, force_n, where=(force_n < 0), color="#4488FF", alpha=0.3, label="← Left (−)"
        )
    ]
    ax3.axhline(y=0, color="w", linestyle="--", alpha=0.3)
    ax3.set_xlim(x_am.min(), x_am.max())
    ax3.set_xlabel("X (am)", family="Monospace")
    ax3.set_ylabel("Force (N)", family="Monospace", fontsize=9)
    ax3.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    ax3.legend(loc="upper right", fontsize=9)
    ax3.grid(True, alpha=0.2)
    f_max = np.max(np.abs(force_n)) if np.max(np.abs(force_n)) > 0 else 1e-10
    ax3.set_ylim(-f_max * 1.3, f_max * 1.3)

    # --- Coulomb reference text (bottom of force panel) ---
    # Shows: Coulomb force, wave force, ratio, separation
    # Color: yellow = direction match, red = WRONG DIRECTION
    coulomb_text = ax3.text(
        0.5,
        0.02,
        "",
        transform=ax3.transAxes,
        fontsize=9,
        family="Monospace",
        ha="center",
        va="bottom",
        color="#FFCC00",
        visible=False,
        bbox=dict(facecolor="#222222", edgecolor="#FFCC00", boxstyle="round,pad=0.3", alpha=0.7),
    )

    # --- Force annotations at each WC position ---
    # Blended transform: x in data coords (am), y in axes fraction (0–1)
    wc_force_texts = []
    for wc in wave_centers:
        txt = ax3.text(
            wc.x_am,
            0.97,
            "",
            transform=blended_transform_factory(ax3.transData, ax3.transAxes),
            fontsize=9,
            family="Monospace",
            fontweight="bold",
            ha="center",
            va="top",
            bbox=dict(facecolor=colormap.DARK_GRAY[1], edgecolor="none", pad=2),
        )
        wc_force_texts.append(txt)

    # --- Separation Slider ---
    sep_init = wc_separation / lam_am
    slider_sep = Slider(
        ax_slider,
        "WC Separation (λ)",
        0.0,
        15.0,
        valinit=sep_init,
        valstep=SEPARATION_STEP,
        color=colormap.viridis_palette[3][1],
    )

    # --- BASE-WAVE Toggle Button ---
    base_txt = fig.text(
        0.10,
        0.93,
        "BASE",
        fontsize=11,
        family="Monospace",
        fontweight="bold",
        ha="center",
        va="center",
        color="yellow",
        picker=True,
        bbox=dict(facecolor="none", edgecolor="yellow", boxstyle="round,pad=0.3", alpha=0.5),
    )

    def update_base_text_style():
        if base_active[0]:
            base_txt.set_color("yellow")
            base_txt.get_bbox_patch().set_edgecolor("yellow")
            base_txt.get_bbox_patch().set_alpha(0.5)
        else:
            base_txt.set_color("#666666")
            base_txt.get_bbox_patch().set_edgecolor("#666666")
            base_txt.get_bbox_patch().set_alpha(0.2)

    # --- WC Toggle Buttons ---
    wc_labels = []
    wc_texts = []
    for i, wc in enumerate(wave_centers):
        charge = "+" if wc.phase == 0.0 else "−"
        label = f"WC{i+1} ({charge})"
        wc_labels.append(label)
        txt = fig.text(
            0.15 + i * 0.06,
            0.93,
            label,
            fontsize=11,
            family="Monospace",
            fontweight="bold",
            ha="center",
            va="center",
            color="yellow",
            picker=True,
            bbox=dict(facecolor="none", edgecolor="yellow", boxstyle="round,pad=0.3", alpha=0.5),
        )
        wc_texts.append(txt)

    def update_wc_text_style():
        for i, txt in enumerate(wc_texts):
            if wc_active[i]:
                txt.set_color("yellow")
                txt.get_bbox_patch().set_edgecolor("yellow")
                txt.get_bbox_patch().set_alpha(0.5)
            else:
                txt.set_color("#666666")
                txt.get_bbox_patch().set_edgecolor("#666666")
                txt.get_bbox_patch().set_alpha(0.2)

    # --- Phase Offset Toggle ---
    phase_state = {"opposite": True}
    phase_txt = fig.text(
        0.86,
        0.93,
        "Phase Δ: 180° (opposite)",
        fontsize=11,
        family="Monospace",
        fontweight="bold",
        ha="center",
        va="center",
        color="#88FF00",
        picker=True,
        bbox=dict(facecolor="none", edgecolor="#88FF00", boxstyle="round,pad=0.3", alpha=0.5),
    )

    def update_phase_text():
        if phase_state["opposite"]:
            phase_txt.set_text("Phase Δ: 180° (opposite)")
            phase_txt.set_color("#88FF00")
            phase_txt.get_bbox_patch().set_edgecolor("#88FF00")
        else:
            phase_txt.set_text("Phase Δ: 0° (same)")
            phase_txt.set_color("#FF8800")
            phase_txt.get_bbox_patch().set_edgecolor("#FF8800")

    def apply_phase_offset():
        if phase_state["opposite"]:
            if len(wave_centers) >= 2:
                wave_centers[0].phase = 0.0
                wave_centers[1].phase = np.pi
        else:
            for wc in wave_centers:
                wc.phase = 0.0
        for i, (wc, txt) in enumerate(zip(wave_centers, wc_texts)):
            charge = "+" if wc.phase == 0.0 else "−"
            wc_labels[i] = f"WC{i+1} ({charge})"
            txt.set_text(wc_labels[i])
            color = "#FF4444" if wc.phase == 0.0 else "#4488FF"
            for vl in wc_vlines[i]:
                vl.set_color(color)
            wc_vlines[i][0].set_label(wc_labels[i])
        ax1.legend(loc="upper right", fontsize=9)

    def on_pick(event):
        # Base wave toggle
        if event.artist == base_txt:
            base_active[0] = not base_active[0]
            update_base_text_style()
            refresh_static()
            return
        # Phase toggle
        if event.artist == phase_txt:
            phase_state["opposite"] = not phase_state["opposite"]
            apply_phase_offset()
            update_phase_text()
            update_wc_text_style()
            refresh_static()
            return
        # WC toggles
        for i, txt in enumerate(wc_texts):
            if event.artist == txt:
                wc_active[i] = not wc_active[i]
                update_wc_text_style()
                refresh_static()
                break

    fig.canvas.mpl_connect("pick_event", on_pick)

    # --- Refresh function ---
    # Recomputes all static elements when WC positions, phases, or active state change:
    # RMS envelope, energy density, force field, WC markers, force annotations,
    # and Coulomb reference comparison.

    def refresh_static():
        """Recompute and redraw all static elements (RMS, energy, force, markers)."""
        nonlocal rms, energy_j, force_n

        rms, energy, force = recompute_static()
        energy_j = energy * ENERGY_TO_J
        force_n = force * FORCE_TO_N

        line_rms_pos.set_ydata(rms)
        line_rms_neg.set_ydata(-rms)
        pm = max(np.max(rms) * np.sqrt(2.0) * 1.3, A0_am * 0.1) * DISP_SCALE
        ax1.set_ylim(-pm, pm)

        line_energy.set_ydata(energy_j)
        fill_energy[0].remove()
        fill_energy[0] = ax2.fill_between(
            x_am, energy_j, color=colormap.viridis_palette[3][1], alpha=0.5
        )
        em = np.max(energy_j) if np.max(energy_j) > 0 else 1e-10
        ax2.set_ylim(0, em * 1.2)

        # Update total energy annotation
        total_energy_text.set_text(f"ΣE = {np.sum(energy_j):.3e} J")

        line_force.set_ydata(force_n)
        fill_force_r[0].remove()
        fill_force_l[0].remove()
        fill_force_r[0] = ax3.fill_between(
            x_am, force_n, where=(force_n > 0), color="#FF4444", alpha=0.3
        )
        fill_force_l[0] = ax3.fill_between(
            x_am, force_n, where=(force_n < 0), color="#4488FF", alpha=0.3
        )
        fm = np.max(np.abs(force_n)) if np.max(np.abs(force_n)) > 0 else 1e-10
        ax3.set_ylim(-fm * 1.3, fm * 1.3)

        # Pre-pass: collect force values at each active WC position
        active_indices = [i for i in range(len(wave_centers)) if wc_active[i]]
        wc_f_vals = []
        for i in range(len(wave_centers)):
            if wc_active[i]:
                idx = int(np.argmin(np.abs(x_am - wave_centers[i].x_am)))
                wc_f_vals.append(force_n[idx])
            else:
                wc_f_vals.append(None)

        # Detect attraction / repulsion when exactly 2 WCs are both active.
        # Sort by position to identify left vs right WC regardless of list order.
        both_active = len(active_indices) == 2
        is_attraction = is_repulsion = False

        if both_active:
            i_left, i_right = sorted(active_indices, key=lambda i: wave_centers[i].x_am)
            f_left = wc_f_vals[i_left]
            f_right = wc_f_vals[i_right]
            is_attraction = f_left > 0 and f_right < 0  # forces point toward each other
            is_repulsion = f_left < 0 and f_right > 0  # forces point away from each other

            # Coulomb reference: F = ke · q1·q2 / r²
            ke = constants.COULOMB_CONSTANT  # N·m²/C²
            qe = constants.ELEMENTARY_CHARGE  # C
            sep_m = (
                abs(wave_centers[i_right].x_am - wave_centers[i_left].x_am) * constants.ATTOMETER
            )
            if sep_m > 0:
                opposite_phase = (
                    abs(wave_centers[i_left].phase - wave_centers[i_right].phase) > 1.0
                )
                charge_product = -(qe**2) if opposite_phase else qe**2
                coulomb_f = ke * charge_product / sep_m**2
                sep_lam = sep_m / constants.EWAVE_LENGTH
                # Check if wave force direction matches Coulomb prediction
                # Coulomb negative = attraction (forces point inward), positive = repulsion (outward)
                # For left WC: attraction means force points right (+), repulsion means left (-)
                wave_f = (
                    f_left if f_left is not None else 0.0
                )  # left WC force (+ = right, - = left)
                coulomb_dir = -coulomb_f  # flip: coulomb negative(attract) → left WC goes right(+)
                signs_match = (wave_f * coulomb_dir) > 0

                if signs_match:
                    ratio = abs(wave_f) / abs(coulomb_f) if coulomb_f != 0 else 0.0
                    coulomb_text.set_text(
                        f"Coulomb: {coulomb_f:.3e} N  |  Wave: {wave_f:.3e} N  |  "
                        f"Ratio: {ratio:.3e}  |  Sep: {sep_lam:.2f}λ"
                    )
                    coulomb_text.set_color("#FFCC00")
                    coulomb_text.get_bbox_patch().set_edgecolor("#FFCC00")
                else:
                    coulomb_text.set_text(
                        f"WRONG DIRECTION  |  Coulomb: {coulomb_f:.3e} N  |  "
                        f"Wave: {wave_f:.3e} N  |  Sep: {sep_lam:.2f}λ"
                    )
                    coulomb_text.set_color("#FF4444")
                    coulomb_text.get_bbox_patch().set_edgecolor("#FF4444")
                coulomb_text.set_visible(True)
            else:
                coulomb_text.set_visible(False)
        else:
            coulomb_text.set_visible(False)

        # Update WC markers and force annotations
        for i, (vl1, vl2, vl3) in enumerate(wc_vlines):
            x_pos = wave_centers[i].x_am
            visible = wc_active[i]
            for vl in (vl1, vl2, vl3):
                vl.set_xdata([x_pos, x_pos])
                vl.set_visible(visible)

            ftxt = wc_force_texts[i]
            ftxt.set_x(x_pos)
            if not visible:
                ftxt.set_text("")
                continue

            f_val = wc_f_vals[i]
            if both_active and is_attraction:
                label = (
                    f">>> Attract\n{f_val:.3e} N" if i == i_left else f"Attract <<<\n{f_val:.3e} N"
                )
                color = "#88FF00"
            elif both_active and is_repulsion:
                label = f"<<< Repel\n{f_val:.3e} N" if i == i_left else f"Repel >>>\n{f_val:.3e} N"
                color = "#FF8800"
            elif f_val > 0:
                label = f"→ {f_val:.3e} N"
                color = "#FF4444"
            elif f_val < 0:
                label = f"← {f_val:.3e} N"
                color = "#4488FF"
            else:
                label = f"0 N"
                color = "#aaaaaa"
            ftxt.set_text(label)
            ftxt.set_color(color)

        fig.canvas.draw_idle()

    def on_slider_change(val):
        new_sep = val * lam_am
        if len(wave_centers) >= 2:
            wave_centers[0].x_am = snap_to_grid(-new_sep / 2)
            wave_centers[1].x_am = snap_to_grid(+new_sep / 2)
        refresh_static()

    slider_sep.on_changed(on_slider_change)

    # --- Animation ---

    def update_plot(frame):
        t_rs = (frame / ANIMATION_FRAMES) * period_rs
        active = get_active_centers()
        if base_active[0] and active:
            psi = compute_combined_displacement(x_am, t_rs, active)
        elif base_active[0]:
            psi = compute_base_displacement(x_am, t_rs)
        elif active:
            # WCs only, no base wave
            psi = np.zeros_like(x_am)
            for wc in active:
                r = np.abs(x_am - wc.x_am)
                kr = k * r
                w = _compute_weight(r)
                wt_phi = omega * t_rs + wc.phase
                with np.errstate(divide="ignore", invalid="ignore"):
                    wave = np.where(
                        kr < 1e-10,
                        2.0 * np.cos(wt_phi),
                        (w * np.sin(kr + wt_phi) + np.sin(kr - wt_phi)) / kr,
                    )
                psi += wc.amplitude * wave
        else:
            psi = np.zeros_like(x_am)
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

    refresh_static()

    return anim, slider_sep, wc_texts, phase_txt


if __name__ == "__main__":
    anim, _slider, _wc_texts, _phase_txt = plot_sandbox()
    plt.show()
