"""
1D WAVE ENGINE SANDBOX v2 — Force Unification Research Tool

Lightweight 1D wave engine for rapid prototyping and validation of:
- Wave equations (weighted partial standing wave as primary)
- Phasor superposition (exact analytical amplitude)
- Energy density and force field computation

3 panels:
  1. Displacement ψ(x,t) + Phasor RMS envelope (overlay)
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

# SI unit conversion factors (internal → SI)
# 1 qg·am²/rs² = 1e-33 kg · (1e-18 m)² / (1e-27 s)² = 1e-15 J
# 1 qg·am/rs²  = 1e-33 kg · (1e-18 m) / (1e-27 s)²  = 1e+3  N
ENERGY_TO_J = constants.QUECTOGRAM * constants.ATTOMETER**2 / constants.RONTOSECOND**2
FORCE_TO_N = constants.QUECTOGRAM * constants.ATTOMETER / constants.RONTOSECOND**2


# ================================================================
# Wave Center Configuration
# ================================================================


class WaveCenter:
    """A single wave center with position, phase offset, and amplitude."""

    def __init__(self, x_am: float, phase: float = 0.0, amplitude: float = A0_am):
        self.x_am = x_am  # position on x-axis (am)
        self.phase = phase  # source_offset: 0 = positron, π = electron
        self.amplitude = amplitude  # base amplitude (am)


# Default: 2 opposite-charge WCs separated by 4λ
wc_separation = 4 * lam_am
wave_centers = [
    WaveCenter(x_am=-wc_separation / 2, phase=0.0),  # positron
    WaveCenter(x_am=+wc_separation / 2, phase=np.pi),  # electron
]


# ================================================================
# Spatial Domain
# ================================================================

# Domain must be wide enough for max slider separation (10λ) + padding
max_sep_am = 1 * lam_am  # max slider value in am
domain_half = max_sep_am / 2 + 5 * lam_am  # WC at edge + 5λ padding each side
x_am = np.linspace(-domain_half, +domain_half, 4000)
dx_am = x_am[1] - x_am[0]  # grid spacing


def snap_to_grid(pos: float) -> float:
    """Snap a position to the nearest grid point to avoid asymmetric gradient artifacts."""
    idx = np.argmin(np.abs(x_am - pos))
    return float(x_am[idx])


# Snap initial WC positions to grid
for wc in wave_centers:
    wc.x_am = snap_to_grid(wc.x_am)


# ================================================================
# Wave Equation: Weighted Partial Standing Wave
# ================================================================
# ψ = A · [w(r)·sin(kr + ωt + φ) + sin(kr - ωt - φ)] / kr
#
# w(r) = 1 / (1 + (r / (transition·λ))^power)
#
# Standing limit (w=1): 2·sin(kr)·cos(ωt+φ) / kr
# Traveling limit (w=0): sin(kr - ωt - φ) / kr

TRANSITION_LAM = 1.25  # standing wave region in wavelengths
WEIGHT_POWER = 8  # sharpness of standing → traveling transition


def compute_weight(r_am: np.ndarray) -> np.ndarray:
    """In-wave weight: 1 near center (standing), 0 far away (traveling)."""
    return 1.0 / (1.0 + (r_am / (TRANSITION_LAM * lam_am)) ** WEIGHT_POWER)


def compute_displacement(x_am: np.ndarray, t_rs: float) -> np.ndarray:
    """Compute total displacement from all wave centers at time t.

    Returns superposition of all WCs using weighted partial standing wave.
    """
    psi_total = np.zeros_like(x_am)

    for wc in wave_centers:
        r_am = np.abs(x_am - wc.x_am)
        kr = k * r_am
        wt_phi = omega * t_rs + wc.phase
        w = compute_weight(r_am)

        # Combined partially standing wave
        # Center safe: kr=0 → limit is 2·cos(ωt+φ)
        with np.errstate(divide="ignore", invalid="ignore"):
            wave = np.where(
                kr < 1e-10,
                2.0 * np.cos(wt_phi),  # center limit
                (w * np.sin(kr + wt_phi) + np.sin(kr - wt_phi)) / kr,
            )

        psi_total += wc.amplitude * wave

    return psi_total


# ================================================================
# Phasor Superposition (Analytical Amplitude)
# ================================================================
# ψ_total(t) = P·cos(ωt) + Q·sin(ωt)
# Peak = √(P² + Q²),  RMS = Peak / √2
#
# Per WC: ψ = A·[(w+1)·sin(kr)·cos(ωt+φ) + (w-1)·cos(kr)·sin(ωt+φ)] / kr
# Rotate by φ to shared cos(ωt)/sin(ωt) basis.


def compute_phasor_rms(x_am: np.ndarray) -> np.ndarray:
    """Compute exact phasor RMS amplitude at each x position.

    Returns RMS = √(P² + Q²) / √2, where P and Q are the
    accumulated cos(ωt) and sin(ωt) coefficients from all WCs.
    """
    P = np.zeros_like(x_am)
    Q = np.zeros_like(x_am)

    for wc in wave_centers:
        r_am = np.abs(x_am - wc.x_am)
        kr = k * r_am
        w = compute_weight(r_am)

        # Phasor coefficients for this WC
        with np.errstate(divide="ignore", invalid="ignore"):
            C_n = np.where(
                kr < 1e-10,
                2.0 * wc.amplitude,  # center limit: (w+1)·sin(kr)/kr → 2
                wc.amplitude * (w + 1.0) * np.sin(kr) / kr,
            )
            S_n = np.where(
                kr < 1e-10,
                0.0,  # center limit: (w-1)·cos(kr)/kr → 0
                wc.amplitude * (w - 1.0) * np.cos(kr) / kr,
            )

        # Rotate by source_offset to shared cos(ωt)/sin(ωt) basis
        cos_phi = np.cos(wc.phase)
        sin_phi = np.sin(wc.phase)
        P += C_n * cos_phi + S_n * sin_phi
        Q += -C_n * sin_phi + S_n * cos_phi

    # RMS = peak / √2
    peak = np.sqrt(P**2 + Q**2)
    return peak / np.sqrt(2.0)


# ================================================================
# Energy Density and Force
# ================================================================
# E = ρ · V · (f · A)²  where V = dx³ (voxel volume)
# F = -∇E = -2 · ρ · V · f² · A · ∇A


def compute_energy_density(rms_am: np.ndarray) -> np.ndarray:
    """Energy density at each x: E = ρ · dx³ · (f · A)²."""
    return rho_qgam * dx_am**3 * (f_rHz * rms_am) ** 2


def compute_force_field(rms_am: np.ndarray) -> np.ndarray:
    """Force field: F = -2 · ρ · V · f² · A · ∇A.

    Uses central differences for gradient.
    """
    grad_A = np.gradient(rms_am, dx_am)
    force = -2.0 * rho_qgam * dx_am**3 * f_rHz**2 * rms_am * grad_A
    return force


# ================================================================
# Animation and Plotting
# ================================================================

# Animation parameters
ANIMATION_FRAMES = 100
ANIMATION_INTERVAL = 50  # ms (20 FPS)
START_PAUSED = False

_MODULE_DIR = Path(__file__).parent
PLOT_DIR = _MODULE_DIR / "_plots"


def plot_sandbox():
    """Animate the 1D wave sandbox with 3 panels + interactive controls.

    Panel 1: Displacement ψ(x,t) + Phasor RMS envelope
    Panel 2: Energy density E(x)
    Panel 3: Force field F(x)

    Controls:
        SPACE: Pause/Resume
        LEFT/RIGHT: Step frame when paused
        Slider: WC separation distance
        Checkboxes: Toggle individual WCs on/off
    """
    from matplotlib.widgets import Slider

    state = {"paused": START_PAUSED, "frame": 0}

    # Track active WCs (all on by default)
    wc_active = [True] * len(wave_centers)
    # Store original positions for separation slider
    wc_original_positions = [wc.x_am for wc in wave_centers]

    def get_active_centers():
        """Return list of currently active wave centers."""
        return [wc for wc, active in zip(wave_centers, wc_active) if active]

    def recompute_static():
        """Recompute phasor RMS, energy, force from active WCs."""
        active = get_active_centers()
        if not active:
            return np.zeros_like(x_am), np.zeros_like(x_am), np.zeros_like(x_am)
        # Temporarily swap wave_centers for computation
        orig = wave_centers.copy()
        wave_centers.clear()
        wave_centers.extend(active)
        rms = compute_phasor_rms(x_am)
        energy = compute_energy_density(rms)
        force = compute_force_field(rms)
        wave_centers.clear()
        wave_centers.extend(orig)
        return rms, energy, force

    # Initial computation (internal units → convert to SI for display)
    rms, energy, force = recompute_static()
    energy = energy * ENERGY_TO_J
    force = force * FORCE_TO_N

    # Setup figure with space for widgets
    plt.style.use("dark_background")
    fig = plt.figure(figsize=(16, 9), facecolor=colormap.DARK_GRAY[1])

    # GridSpec: 3 plot rows (slider placed manually below)
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

    ax1 = fig.add_subplot(gs[0, 0])  # displacement + RMS
    ax2 = fig.add_subplot(gs[1, 0])  # energy
    ax3 = fig.add_subplot(gs[2, 0])  # force
    # Slider: manually positioned, narrower than plots
    ax_slider = fig.add_axes([0.15, 0.02, 0.70, 0.025])  # [left, bottom, width, height]

    fig.suptitle(
        "OPENWAVE Analytics — 1D Wave Force Research",
        fontsize=18,
        family="Monospace",
        y=0.98,
    )

    # --- Panel 1: Displacement + Phasor RMS ---
    (line_psi,) = ax1.plot(
        [], [], color=colormap.viridis_palette[2][1], linewidth=1.5, alpha=0.8, label="ψ(x,t)"
    )
    (line_rms_pos,) = ax1.plot(
        x_am, rms, color=colormap.viridis_palette[4][1], linewidth=2, label="Phasor RMS"
    )
    (line_rms_neg,) = ax1.plot(
        x_am, -rms, color=colormap.viridis_palette[4][1], linewidth=2, alpha=0.5
    )

    # WC marker lines (stored for update)
    wc_vlines = []
    for i, wc in enumerate(wave_centers):
        color = "#FF4444" if wc.phase == 0.0 else "#4488FF"
        charge = "+" if wc.phase == 0.0 else "−"
        vl = ax1.axvline(
            x=wc.x_am, color=color, linestyle="--", alpha=0.6, label=f"WC{i+1} ({charge})"
        )
        vl2 = ax2.axvline(x=wc.x_am, color=color, linestyle="--", alpha=0.4)
        vl3 = ax3.axvline(x=wc.x_am, color=color, linestyle="--", alpha=0.4)
        wc_vlines.append((vl, vl2, vl3))

    ax1.axhline(y=0, color="w", linestyle="--", alpha=0.2)
    ax1.set_xlim(x_am.min(), x_am.max())
    psi_max = A0_am * 2
    ax1.set_ylim(-psi_max, psi_max)
    ax1.set_ylabel("Displacement (am)", family="Monospace")
    ax1.legend(loc="upper right", fontsize=9)
    ax1.grid(True, alpha=0.2)

    # --- Panel 2: Energy Density ---
    (line_energy,) = ax2.plot(
        x_am,
        energy,
        color=colormap.viridis_palette[3][1],
        linewidth=1.5,
        label="E(x) = ρV(fA)² [J]",
    )
    # Mutable list so fill can be replaced in refresh_static
    fill_energy = [ax2.fill_between(x_am, energy, color=colormap.viridis_palette[3][1], alpha=0.5)]
    ax2.set_xlim(x_am.min(), x_am.max())
    ax2.set_ylabel("Energy (J)", family="Monospace", fontsize=9)
    ax2.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    ax2.legend(loc="upper right", fontsize=9)
    ax2.grid(True, alpha=0.2)

    # --- Panel 3: Force Field ---
    (line_force,) = ax3.plot(x_am, force, color="w", linewidth=1, alpha=0.8)
    # Mutable lists so fills can be replaced in refresh_static
    fill_force_r = [
        ax3.fill_between(
            x_am,
            force,
            where=(force > 0),
            color="#FF4444",
            alpha=0.3,
            label="→ Right Force Vector (+)",
        )
    ]
    fill_force_l = [
        ax3.fill_between(
            x_am,
            force,
            where=(force < 0),
            color="#4488FF",
            alpha=0.3,
            label="← Left Force Vector (−)",
        )
    ]
    ax3.axhline(y=0, color="w", linestyle="--", alpha=0.3)
    ax3.set_xlim(x_am.min(), x_am.max())
    ax3.set_xlabel("X (am)", family="Monospace")
    ax3.set_ylabel("Force (N)", family="Monospace", fontsize=9)
    ax3.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    ax3.legend(loc="upper right", fontsize=9)
    ax3.grid(True, alpha=0.2)

    # --- Coulomb reference text (bottom of force panel) ---
    coulomb_text = ax3.text(
        0.5, 0.02, "",
        transform=ax3.transAxes,
        fontsize=9, family="Monospace",
        ha="center", va="bottom",
        color="#FFCC00",
        visible=False,
        bbox=dict(facecolor="#222222", edgecolor="#FFCC00", boxstyle="round,pad=0.3", alpha=0.7),
    )

    # --- Force annotations at each WC position ---
    # Blended transform: x in data coords (am), y in axes fraction (0–1)
    from matplotlib.transforms import blended_transform_factory

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
    sep_min, sep_max = 0.0, 10.0  # in wavelengths
    sep_init = wc_separation / lam_am
    slider_sep = Slider(
        ax_slider,
        "WC Separation (λ)",
        sep_min,
        sep_max,
        valinit=sep_init,
        valstep=0.25,
        color=colormap.viridis_palette[3][1],
    )

    # --- WC Toggle Buttons (click text to toggle on/off, top bar) ---
    wc_labels = []
    wc_texts = []
    n_wcs = len(wave_centers)
    for i, wc in enumerate(wave_centers):
        charge = "+" if wc.phase == 0.0 else "−"
        label = f"WC{i+1} ({charge})"
        wc_labels.append(label)
        # Lay out horizontally in the title area, left-aligned
        x_pos = 0.10 + i * 0.10  # left-to-right spacing
        txt = fig.text(
            x_pos,
            0.95,
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
        """Update text colors based on active state."""
        for i, txt in enumerate(wc_texts):
            if wc_active[i]:
                txt.set_color("yellow")
                txt.get_bbox_patch().set_edgecolor("yellow")
                txt.get_bbox_patch().set_alpha(0.5)
            else:
                txt.set_color("#666666")
                txt.get_bbox_patch().set_edgecolor("#666666")
                txt.get_bbox_patch().set_alpha(0.2)

    # --- Phase Offset Toggle (top bar, after WC labels) ---
    phase_state = {"opposite": True}  # True = 0 vs π (opposite), False = same phase
    x_phase = 0.98 - 0.12  # right-to-left spacing
    phase_label_text = "Phase Δ: 180° (opposite charges)"
    phase_txt = fig.text(
        x_phase,
        0.95,
        phase_label_text,
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
        """Update phase toggle label."""
        if phase_state["opposite"]:
            phase_txt.set_text("Phase Δ: 180° (opposite charges)")
            phase_txt.set_color("#88FF00")
            phase_txt.get_bbox_patch().set_edgecolor("#88FF00")
        else:
            phase_txt.set_text("Phase Δ: 0° (same charges)")
            phase_txt.set_color("#FF8800")
            phase_txt.get_bbox_patch().set_edgecolor("#FF8800")

    def apply_phase_offset():
        """Set WC phases based on toggle state."""
        if phase_state["opposite"]:
            # Opposite: first WC = 0, second WC = π
            if len(wave_centers) >= 2:
                wave_centers[0].phase = 0.0
                wave_centers[1].phase = np.pi
        else:
            # Same: all WCs share phase 0
            for wc in wave_centers:
                wc.phase = 0.0
        # Update WC toggle labels and ax1 vline colors/legend to reflect charge
        for i, (wc, txt) in enumerate(zip(wave_centers, wc_texts)):
            charge = "+" if wc.phase == 0.0 else "−"
            wc_labels[i] = f"WC{i+1} ({charge})"
            txt.set_text(wc_labels[i])
            color = "#FF4444" if wc.phase == 0.0 else "#4488FF"
            vl1, vl2, vl3 = wc_vlines[i]
            for vl in (vl1, vl2, vl3):
                vl.set_color(color)
            vl1.set_label(wc_labels[i])
        ax1.legend(loc="upper right", fontsize=9)

    def on_pick(event):
        """Handle clicks on WC toggles and phase toggle."""
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

    # Time display
    time_text = ax1.text(
        0.02,
        0.95,
        "",
        transform=ax1.transAxes,
        fontsize=10,
        family="Monospace",
        verticalalignment="top",
    )

    def refresh_static():
        """Recompute and redraw all static elements (RMS, energy, force, markers)."""
        nonlocal rms, energy, force, psi_max

        rms, energy, force = recompute_static()
        energy = energy * ENERGY_TO_J
        force = force * FORCE_TO_N

        # Update RMS lines
        line_rms_pos.set_ydata(rms)
        line_rms_neg.set_ydata(-rms)

        # Update y-limits for panel 1
        psi_max = A0_am * 2
        ax1.set_ylim(-psi_max, psi_max)

        # Update energy line + fill
        line_energy.set_ydata(energy)
        fill_energy[0].remove()
        fill_energy[0] = ax2.fill_between(
            x_am, energy, color=colormap.viridis_palette[3][1], alpha=0.5
        )
        e_max = np.max(energy) if np.max(energy) > 0 else 1e-10
        ax2.set_ylim(0, e_max * 1.2)

        # Update force line + fills
        line_force.set_ydata(force)
        fill_force_r[0].remove()
        fill_force_l[0].remove()
        fill_force_r[0] = ax3.fill_between(
            x_am,
            force,
            where=(force > 0),
            color="#FF4444",
            alpha=0.3,
        )
        fill_force_l[0] = ax3.fill_between(
            x_am,
            force,
            where=(force < 0),
            color="#4488FF",
            alpha=0.3,
        )
        f_max = np.max(np.abs(force)) if np.max(np.abs(force)) > 0 else 1e-10
        ax3.set_ylim(-f_max * 1.3, f_max * 1.3)

        # Pre-pass: collect force values at each active WC position
        wc_f_vals = []
        for i in range(len(wave_centers)):
            if wc_active[i]:
                idx = int(np.argmin(np.abs(x_am - wave_centers[i].x_am)))
                wc_f_vals.append(force[idx])
            else:
                wc_f_vals.append(None)

        # Detect attraction / repulsion when exactly 2 WCs are both active.
        # Sort by position to identify left vs right WC regardless of list order.
        active_indices = [i for i in range(len(wave_centers)) if wc_active[i]]
        both_active = len(active_indices) == 2
        if both_active:
            i_left, i_right = sorted(active_indices, key=lambda i: wave_centers[i].x_am)
            f_left = wc_f_vals[i_left]
            f_right = wc_f_vals[i_right]
            is_attraction = f_left > 0 and f_right < 0  # forces point toward each other
            is_repulsion = f_left < 0 and f_right > 0  # forces point away from each other

            # Coulomb reference: F = ke · q1·q2 / r²
            ke = constants.COULOMB_CONSTANT  # N·m²/C²
            qe = constants.ELEMENTARY_CHARGE  # C
            sep_m = abs(wave_centers[i_right].x_am - wave_centers[i_left].x_am) * constants.ATTOMETER
            if sep_m > 0:
                opposite_phase = abs(wave_centers[i_left].phase - wave_centers[i_right].phase) > 1.0
                charge_product = -qe**2 if opposite_phase else qe**2
                coulomb_f = ke * charge_product / sep_m**2  # Newtons (signed)
                sep_lam = sep_m / constants.EWAVE_LENGTH
                wave_f_signed = f_left if f_left is not None else 0.0  # left WC force (+ = right, - = left)
                # Check if wave force direction matches Coulomb prediction
                # Coulomb negative = attraction (forces point inward), positive = repulsion (outward)
                # For left WC: attraction means force points right (+), repulsion means left (-)
                coulomb_direction_for_left = -coulomb_f  # flip: coulomb negative(attract) → left WC goes right(+)
                signs_match = (wave_f_signed * coulomb_direction_for_left) > 0

                if signs_match:
                    ratio = abs(wave_f_signed) / abs(coulomb_f) if coulomb_f != 0 else 0.0
                    coulomb_text.set_text(
                        f"Coulomb: {coulomb_f:.3e} N  |  Wave: {wave_f_signed:.3e} N  |  "
                        f"Ratio: {ratio:.3e}  |  Sep: {sep_lam:.2f}λ"
                    )
                    coulomb_text.set_color("#FFCC00")
                    coulomb_text.get_bbox_patch().set_edgecolor("#FFCC00")
                else:
                    coulomb_text.set_text(
                        f"WRONG DIRECTION  |  Coulomb: {coulomb_f:.3e} N  |  "
                        f"Wave: {wave_f_signed:.3e} N  |  Sep: {sep_lam:.2f}λ"
                    )
                    coulomb_text.set_color("#FF4444")
                    coulomb_text.get_bbox_patch().set_edgecolor("#FF4444")
                coulomb_text.set_visible(True)
            else:
                coulomb_text.set_visible(False)
        else:
            is_attraction = is_repulsion = False
            coulomb_text.set_visible(False)

        # Update WC marker positions, visibility, and force annotations
        for i, (vl1, vl2, vl3) in enumerate(wc_vlines):
            x_pos = wave_centers[i].x_am
            visible = wc_active[i]
            for vl in (vl1, vl2, vl3):
                vl.set_xdata([x_pos, x_pos])
                vl.set_visible(visible)

            # Force annotation at this WC position
            ftxt = wc_force_texts[i]
            ftxt.set_x(x_pos)
            if not visible:
                ftxt.set_text("")
                continue

            f_val = wc_f_vals[i]
            if both_active and is_attraction:
                label = (
                    f">>> Attraction\n{f_val:.3e} N"
                    if i == i_left
                    else f"Attraction <<<\n{f_val:.3e} N"
                )
                color = "#88FF00"  # green — opposite charges attract
            elif both_active and is_repulsion:
                label = (
                    f"<<< Repulsion\n{f_val:.3e} N"
                    if i == i_left
                    else f"Repulsion >>>\n{f_val:.3e} N"
                )
                color = "#FF8800"  # orange — same charges repel
            elif f_val > 0:
                label = f"→ Right\n{f_val:.3e} N"
                color = "#FF4444"
            elif f_val < 0:
                label = f"← Left\n{f_val:.3e} N"
                color = "#4488FF"
            else:
                label = f"Neutral\n{f_val:.3e} N"
                color = "#aaaaaa"
            ftxt.set_text(label)
            ftxt.set_color(color)

        fig.canvas.draw_idle()

    def on_slider_change(val):
        """Update WC positions from separation slider."""
        new_sep = val * lam_am
        # Reposition WCs symmetrically, snapped to grid
        if len(wave_centers) >= 2:
            wave_centers[0].x_am = snap_to_grid(-new_sep / 2)
            wave_centers[1].x_am = snap_to_grid(+new_sep / 2)
        refresh_static()

    slider_sep.on_changed(on_slider_change)

    def update_plot(frame):
        """Update displacement for given frame."""
        t_rs = (frame / ANIMATION_FRAMES) * period_rs

        # Compute displacement from active WCs only
        active = get_active_centers()
        if active:
            orig = wave_centers.copy()
            wave_centers.clear()
            wave_centers.extend(active)
            psi = compute_displacement(x_am, t_rs)
            wave_centers.clear()
            wave_centers.extend(orig)
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

    # blit=False required for widget redraws
    anim = FuncAnimation(
        fig,
        animate,
        frames=ANIMATION_FRAMES,
        interval=ANIMATION_INTERVAL,
        blit=False,
    )

    # Initial static draw
    refresh_static()

    return anim, slider_sep, wc_texts, phase_txt  # keep references to prevent GC


if __name__ == "__main__":
    anim, _slider, _wc_texts, _phase_txt = plot_sandbox()
    plt.show()
