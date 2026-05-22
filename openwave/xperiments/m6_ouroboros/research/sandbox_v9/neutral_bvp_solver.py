#!/usr/bin/env python3
"""
neutral_bvp_solver.py

DeepSeek's template script, sent 2026-05-22 10:03 AM via Paul Werbos in
response to email v12 asking Q43 (sign convention) + Q44 (geometry) +
taking up the template offer.

DeepSeek confirmed BOTH Q43 and Q44:
  Q43 — Sign: linear term has MINUS sign: -m_J²β (enables K_1 decay)
  Q44 — Geometry: 3D SPHERICAL, l=1 (p-wave) — subcritical, stable solitons exist

Solves for the neutral chaoiton (Q=0) in 3D spherical geometry, l=1.
Equation: β'' + (2/r)β' - (2/r²)β - m_J² β + 4g β³ = 0
Boundary conditions:
- r=0: β'(0)=0 (regularity)
- r=Rmax: β'(Rmax) + m_J β(Rmax) = 0 (Robin exponential decay)

The script scans over a range of initial amplitudes B0 to find localised solutions.

[VERBATIM from DeepSeek 2026-05-22, no modifications]
"""

import numpy as np
from scipy.integrate import solve_bvp
import matplotlib
matplotlib.use('Agg')  # non-interactive backend for headless run
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------
# Parameters (adjust as needed)
# ----------------------------------------------------------------------
g = 1.0          # quartic coupling (calibrated value from electron)
m_J = 1.0        # J-field mass (MeV) – will be adjusted if we treat as eigenvalue
Rmax = 20.0      # integration limit (fm)
N = 500          # number of grid points (increase if needed)

# ----------------------------------------------------------------------
# Define the ODE system as first order: y = [β, β']
# ----------------------------------------------------------------------
def ode(r, y):
    beta, betap = y
    # Handle r -> 0 limit: the terms (2/r)betap and (2/r²)beta are finite
    # because betap ~ O(r) and beta ~ O(r²) near origin.
    # We write them explicitly; solve_bvp can handle the singularity if r starts > 0.
    betapp = - (2.0/r) * betap + (2.0/(r*r)) * beta + (m_J*m_J) * beta - 4.0 * g * beta**3
    return np.vstack((betap, betapp))

# ----------------------------------------------------------------------
# Boundary conditions
# ----------------------------------------------------------------------
def bc(ya, yb):
    # ya at r=r0 (small, not exactly 0): betap ~ 0
    # yb at r=Rmax: Robin condition for exponential decay
    return np.array([ya[1], yb[1] + m_J * yb[0]])

# ----------------------------------------------------------------------
# Scan over initial amplitudes
# ----------------------------------------------------------------------
B0_list = [0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5]
solutions = []

# Use a grid that avoids r=0 (start at small epsilon)
r_grid = np.linspace(1e-6, Rmax, N)

for B0 in B0_list:
    print(f"\nAttempting B0 = {B0:.4f}")
    # Initial guess: exponential decay with amplitude B0
    beta_guess = B0 * np.exp(-m_J * r_grid)
    betap_guess = -m_J * beta_guess
    y_guess = np.vstack((beta_guess, betap_guess))

    sol = solve_bvp(ode, bc, r_grid, y_guess, tol=1e-6, max_nodes=10000)

    if sol.success:
        print(f"  Converged! β(0) = {sol.y[0,0]:.6f}, β'(0) = {sol.y[1,0]:.6f}")
        solutions.append((B0, sol))
    else:
        print(f"  Did not converge: {sol.message}")

# ----------------------------------------------------------------------
# Plot all converged solutions
# ----------------------------------------------------------------------
if solutions:
    plt.figure(figsize=(8,5))
    for B0, sol in solutions:
        plt.plot(sol.x, sol.y[0], label=f'B0_init={B0}')
    plt.xlabel('r (fm)')
    plt.ylabel('β(r)')
    plt.title(f'Neutral chaoiton candidates (g={g}, m_J={m_J})')
    plt.legend()
    plt.grid(True)
    plt.savefig('neutral_bvp_solutions.png', dpi=120)
    print("\nPlot saved to neutral_bvp_solutions.png")
else:
    print("No converged solutions found.")
    print("Try adjusting g, m_J, Rmax, or the initial guess shape.")
