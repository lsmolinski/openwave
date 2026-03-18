"""
Phase 1 Diagnostic: Force vs Separation Parameter Sweep

Sweeps WC separation from 2λ to 10λ and plots:
- Force at left WC position (wave-computed) for both force methods
- Coulomb reference (1/r²)
- Direction match/mismatch markers

Compares two force computation methods:
  1. A·∇A  (current):  F = -2ρVf² · A · ∇A
  2. ∇E    (alternative): F = -∇(ρVf²A²)

Usage:
    python sweep_force_vs_separation.py
"""

import numpy as np
import matplotlib.pyplot as plt

import wave_engine_1D_v2 as we
from openwave.common import constants, colormap

# ================================================================
# Sweep Configuration
# ================================================================

SEP_MIN_LAM = 2.0  # minimum separation (wavelengths)
SEP_MAX_LAM = 10.0  # maximum separation (wavelengths)
SEP_STEP_LAM = 0.1  # step size (wavelengths)

# Domain sizing: enough room for max separation + padding
DOMAIN_HALF_LAM = SEP_MAX_LAM / 2 + 5  # half-domain in wavelengths
DOMAIN_NPOINTS = 16001  # odd for exact center

lam = we.lam_am
k = we.k
A0 = we.A0_am
f_rHz = we.f_rHz
rho = we.rho_qgam


# ================================================================
# Force Methods
# ================================================================


def force_A_grad_A(rms, dx):
    """Current method: F = -2ρVf² · A · ∇A."""
    grad_A = np.gradient(rms, dx)
    return -2.0 * rho * dx**3 * f_rHz**2 * rms * grad_A


def force_grad_E(rms, dx):
    """Alternative: F = -∇E where E = ρVf²A²."""
    energy = rho * dx**3 * (f_rHz * rms) ** 2
    return -np.gradient(energy, dx)


# ================================================================
# Sweep
# ================================================================


def run_sweep(phase1=0.0, phase2=np.pi):
    """Sweep separation and compute force at left WC for both methods.

    Returns dict with arrays: sep_lam, f_AgradA, f_gradE, f_coulomb.
    """
    x = np.linspace(-DOMAIN_HALF_LAM * lam, DOMAIN_HALF_LAM * lam, DOMAIN_NPOINTS)
    dx = x[1] - x[0]

    seps = np.arange(SEP_MIN_LAM, SEP_MAX_LAM + SEP_STEP_LAM / 2, SEP_STEP_LAM)
    f_AgradA = np.zeros_like(seps)
    f_gradE = np.zeros_like(seps)
    f_coulomb = np.zeros_like(seps)

    opposite_phase = abs(phase1 - phase2) > 1.0
    ke = constants.COULOMB_CONSTANT
    qe = constants.ELEMENTARY_CHARGE
    charge_product = -(qe**2) if opposite_phase else qe**2

    for i, sep_lam in enumerate(seps):
        sep = sep_lam * lam
        x_left = -sep / 2
        x_right = +sep / 2

        # Snap to grid
        idx_left = np.argmin(np.abs(x - x_left))
        idx_right = np.argmin(np.abs(x - x_right))

        we.wave_centers.clear()
        we.wave_centers.extend([
            we.WaveCenter(x_am=x[idx_left], phase=phase1),
            we.WaveCenter(x_am=x[idx_right], phase=phase2),
        ])

        rms = we.compute_phasor_rms(x)

        # Force at left WC position (positive = rightward)
        force1 = force_A_grad_A(rms, dx)
        force2 = force_grad_E(rms, dx)

        f_AgradA[i] = force1[idx_left] * we.FORCE_TO_N
        f_gradE[i] = force2[idx_left] * we.FORCE_TO_N

        # Coulomb reference (force on left WC from right WC)
        sep_m = sep * constants.ATTOMETER
        coulomb_f = ke * charge_product / sep_m**2
        # Coulomb convention: negative = attractive. For left WC:
        # attraction → force points right (+), repulsion → force points left (-)
        f_coulomb[i] = -coulomb_f  # flip to match left-WC sign convention

    return {
        "sep_lam": seps,
        "f_AgradA": f_AgradA,
        "f_gradE": f_gradE,
        "f_coulomb": f_coulomb,
        "opposite_phase": opposite_phase,
    }


# ================================================================
# Plotting
# ================================================================


def plot_sweep(results):
    """Plot force vs separation for both methods + Coulomb reference."""
    sep = results["sep_lam"]
    f1 = results["f_AgradA"]
    f2 = results["f_gradE"]
    fc = results["f_coulomb"]
    opp = results["opposite_phase"]

    charge_label = "opposite charge (attract)" if opp else "same charge (repel)"

    plt.style.use("dark_background")
    fig, axes = plt.subplots(3, 1, figsize=(14, 10), facecolor=colormap.DARK_GRAY[1])
    fig.suptitle(
        f"Force vs Separation — {charge_label}",
        fontsize=16,
        family="Monospace",
        y=0.98,
    )

    # --- Panel 1: Raw force values ---
    ax1 = axes[0]
    ax1.plot(sep, f1, "o-", color="#FF6644", markersize=2, linewidth=1, label="A·∇A (current)")
    ax1.plot(sep, f2, "s-", color="#44AAFF", markersize=2, linewidth=1, label="∇E (alternative)")
    ax1.plot(sep, fc, "--", color="#FFCC00", linewidth=2, label="Coulomb reference")
    ax1.axhline(0, color="w", alpha=0.3, linewidth=0.5)
    ax1.set_ylabel("Force on left WC (N)", family="Monospace", fontsize=10)
    ax1.set_title("Raw force values", fontsize=11, family="Monospace")
    ax1.legend(fontsize=9)
    ax1.grid(True, alpha=0.2)
    ax1.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))

    # --- Panel 2: Force sign (direction match) ---
    ax2 = axes[1]
    # Coulomb sign determines expected direction
    sign_match_1 = np.sign(f1) == np.sign(fc)
    sign_match_2 = np.sign(f2) == np.sign(fc)

    ax2.fill_between(
        sep, 0, 1, where=sign_match_1,
        color="#44FF44", alpha=0.3, transform=ax2.get_xaxis_transform(), label="A·∇A correct"
    )
    ax2.fill_between(
        sep, 0, 1, where=~sign_match_1,
        color="#FF4444", alpha=0.3, transform=ax2.get_xaxis_transform(), label="A·∇A WRONG"
    )

    # Count matches
    n1 = np.sum(sign_match_1)
    n2 = np.sum(sign_match_2)
    ntot = len(sep)
    ax2.set_title(
        f"Direction match: A·∇A = {n1}/{ntot} ({100*n1/ntot:.0f}%)  |  "
        f"∇E = {n2}/{ntot} ({100*n2/ntot:.0f}%)",
        fontsize=11,
        family="Monospace",
    )
    ax2.plot(sep, np.where(sign_match_1, 1, -1), "o", color="#FF6644", markersize=3, label="A·∇A")
    ax2.plot(
        sep, np.where(sign_match_2, 0.8, -0.8), "s", color="#44AAFF", markersize=3, label="∇E"
    )
    ax2.axhline(0, color="w", alpha=0.3)
    ax2.set_ylabel("Direction", family="Monospace", fontsize=10)
    ax2.set_yticks([-1, 0, 1])
    ax2.set_yticklabels(["WRONG", "", "CORRECT"])
    ax2.legend(fontsize=9, loc="lower right")
    ax2.grid(True, alpha=0.2)

    # --- Panel 3: Log magnitude comparison ---
    ax3 = axes[2]
    ax3.semilogy(sep, np.abs(f1), "o-", color="#FF6644", markersize=2, linewidth=1, label="|A·∇A|")
    ax3.semilogy(sep, np.abs(f2), "s-", color="#44AAFF", markersize=2, linewidth=1, label="|∇E|")
    ax3.semilogy(sep, np.abs(fc), "--", color="#FFCC00", linewidth=2, label="|Coulomb|")
    ax3.set_xlabel("Separation (λ)", family="Monospace", fontsize=10)
    ax3.set_ylabel("|Force| (N)", family="Monospace", fontsize=10)
    ax3.set_title("Log magnitude — should follow Coulomb 1/r² slope", fontsize=11, family="Monospace")
    ax3.legend(fontsize=9)
    ax3.grid(True, alpha=0.2)

    plt.tight_layout()
    return fig


# ================================================================
# Main
# ================================================================

if __name__ == "__main__":
    print("Running force vs separation sweep...")

    # Opposite charge (should attract: positive force on left WC)
    print("  Opposite charge sweep...")
    results_opp = run_sweep(phase1=0.0, phase2=np.pi)

    # Same charge (should repel: negative force on left WC)
    print("  Same charge sweep...")
    results_same = run_sweep(phase1=0.0, phase2=0.0)

    fig1 = plot_sweep(results_opp)
    fig2 = plot_sweep(results_same)

    print("Done. Close plots to exit.")
    plt.show()
