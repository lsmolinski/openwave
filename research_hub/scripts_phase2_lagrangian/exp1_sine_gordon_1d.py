"""
Experiment 1: Sine-Gordon 1D Solitons

Numerical Experiment 1 in the Lagrangian Framework Sub-Project.

Hypothesis:
The Sine-Gordon equation produces stable kink solitons that exhibit:
  - Topological stability (kinks cannot dissipate)
  - Pair creation/annihilation (kink + anti-kink collisions)
  - Lorentz contraction (relativistic kinematics from wave equation)

Equation:
    ∂²φ/∂t² - c²·∂²φ/∂x² + (m²c²/ℏ²) · sin(φ) = 0

This comes from the Lagrangian:
    L = ½(∂φ/∂t)² - ½c²(∂φ/∂x)² - (m²c⁴/ℏ²)·(1 - cos(φ))

Spec: ../2a_lagrangian_eval.md  (Experiment 1 section)
Results: ../2b_lagrangian_experiments.md  (Experiment 1 section)
"""

import numpy as np
import matplotlib.pyplot as plt  # noqa: F401  used when plots are added


# ---------------------------------------------------------------------------
# PARAMETERS (fill in based on Setup section in 2b_lagrangian_experiments.md)
# ---------------------------------------------------------------------------

# Grid
N_X = 1024            # number of spatial points
DX = 0.1              # spatial step
DOMAIN = N_X * DX     # total spatial extent

# Physics
C = 1.0               # wave speed (natural units)
M = 1.0               # mass parameter
HBAR = 1.0            # ℏ (natural units)

# Time
DT = 0.05             # time step (must satisfy CFL: c·dt/dx < 1)
N_STEPS = 2000        # total simulation steps

# Initial conditions
KINK_POSITION = -10.0     # initial center of kink
KINK_VELOCITY = 0.0       # initial kink velocity (units of c)
KINK_WIDTH = HBAR / (M * C)   # natural width


# ---------------------------------------------------------------------------
# INITIAL CONDITION: kink ansatz
#   φ(x, 0) = 4 · arctan( exp( γ·(x - x₀) / L ) )
#   γ = 1/√(1-v²/c²) for moving kink (Lorentz boost)
# ---------------------------------------------------------------------------

def kink_profile(x, x0, v, L, c):
    """Static or moving Sine-Gordon kink."""
    gamma = 1.0 / np.sqrt(1.0 - (v / c) ** 2)
    return 4.0 * np.arctan(np.exp(gamma * (x - x0) / L))


# ---------------------------------------------------------------------------
# TIME EVOLUTION: leapfrog
#   φ_{n+1} = 2·φ_n - φ_{n-1} + dt²·[c²·∂²φ/∂x² - (m²c²/ℏ²)·sin(φ)]
# ---------------------------------------------------------------------------

def laplacian_1d(phi, dx):
    """Second spatial derivative via 3-point stencil."""
    lap = np.zeros_like(phi)
    lap[1:-1] = (phi[2:] - 2.0 * phi[1:-1] + phi[:-2]) / dx ** 2
    return lap


def step_sine_gordon(phi, phi_prev, dt, dx, c, m, hbar):
    """One leapfrog step."""
    lap = laplacian_1d(phi, dx)
    mass_term = (m * c / hbar) ** 2 * np.sin(phi)
    accel = c ** 2 * lap - mass_term * c ** 2
    phi_next = 2.0 * phi - phi_prev + dt ** 2 * accel
    return phi_next


# ---------------------------------------------------------------------------
# ENERGY DIAGNOSTICS
# ---------------------------------------------------------------------------

def kink_energy(phi, phi_dot, dx, c, m, hbar):
    """Total Sine-Gordon energy."""
    kinetic = 0.5 * phi_dot ** 2
    gradient = 0.5 * c ** 2 * np.gradient(phi, dx) ** 2
    potential = (m * c ** 2 / hbar) ** 2 * (1.0 - np.cos(phi))
    return np.sum(kinetic + gradient + potential) * dx


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------

def main():
    """Run the experiment and produce diagnostics."""
    x = np.linspace(-DOMAIN / 2, DOMAIN / 2, N_X)
    phi = kink_profile(x, KINK_POSITION, KINK_VELOCITY, KINK_WIDTH, C)
    phi_prev = kink_profile(
        x, KINK_POSITION - KINK_VELOCITY * DT, KINK_VELOCITY, KINK_WIDTH, C
    )

    # TODO: implement evolution loop using step_sine_gordon(phi, phi_prev, ...)
    # TODO: energy tracking via kink_energy()
    # TODO: kink position tracking
    # TODO: tests — kink stability, pair annihilation, Lorentz contraction
    # TODO: plots — phi(x,t) snapshots, energy vs time, kink position vs time

    print("Experiment 1 skeleton — implementation pending")
    print(f"Initial phi range: [{phi.min():.3f}, {phi.max():.3f}]")
    print(f"Phi_prev range:    [{phi_prev.min():.3f}, {phi_prev.max():.3f}]")
    print(f"Grid: {N_X} points, dx={DX}, domain=[{-DOMAIN/2:.1f}, {DOMAIN/2:.1f}]")
    print(f"Initial kink at x={KINK_POSITION}, v={KINK_VELOCITY}c, L={KINK_WIDTH:.3f}")
    print(f"Predicted kink rest energy: E = 8·m·c² = {8 * M * C ** 2:.3f}")


if __name__ == "__main__":
    main()
