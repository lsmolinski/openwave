"""Patched version of DeepSeek's reference script.

The original (`deepseek_reference.py`) does not run as-is:
  1. r = np.linspace(0, Rmax, N) starts at r=0, but ODE uses J/r, V/r terms
     → divide-by-zero at the first grid point.
  2. BC returns 10 residuals for a 9-state ODE with no free params
     → solve_bvp expects 9.

Minimal patches to make it runnable, preserving DeepSeek's intent:
  - Start grid at r = R_MIN = 0.05 (matches our v6).
  - Add lambda_lm as a free parameter `p[0]` so solve_bvp can balance the
    integral constraint (gives 9 states + 1 free param = 10 BCs = matches).
  - omega kept fixed at 1.0 per DeepSeek's parameter list.

This isolates the ODE+H+Q_CS structure (DeepSeek's actual contribution) from
the runnability bugs, allowing a clean H/Q comparison against our v6.
"""

import numpy as np
from scipy.integrate import solve_bvp

# Physical parameters
omega = 1.0
g = 1.0625

# Radial grid — patched to start at R_MIN > 0
R_MIN = 0.05
Rmax = 20.0
N = 500
r = np.linspace(R_MIN, Rmax, N)

# Initial guess (DeepSeek's recipe: exp decay, J faster to break A∝J degeneracy)
V0 = 0.1 * np.exp(-r)
A0 = -0.1 * np.exp(-r)
Q0 = 0.1 * np.exp(-r)
J0 = -0.1 * np.exp(-1.5 * r)
Vp0 = -V0
Ap0 = -A0
Qp0 = -Q0
Jp0 = 1.5 * (-J0)

# State: y = [V, V', A, A', Q, Q', J, J', I]
y0 = np.vstack([V0, Vp0, A0, Ap0, Q0, Qp0, J0, Jp0, np.zeros_like(r)])
p0 = np.array([1.0])  # lambda_lm initial guess


def ode(r, y, p):
    """DeepSeek's ODE, verbatim except lambda_lm now a free parameter."""
    lambda_lm = p[0]
    V, Vp, A, Ap, Q, Qp, J, Jp, I = y
    dV = Vp
    dVp = -(omega**2)*V + Q - lambda_lm*(Jp + J/r) - g*((V**2+Q**2)*V + (V*A - Q*J)*A)
    dA = Ap
    dAp = -(omega**2)*A - J + lambda_lm*(-(Vp + V/r)) - g*((A**2+J**2)*A + (V*A - Q*J)*V)
    dQ = Qp
    dQp = -(omega**2)*Q + V - lambda_lm*(Jp + J/r) - g*((V**2+Q**2)*Q - (V*A - Q*J)*J)
    dJ = Jp
    dJp = -(omega**2)*J - A + lambda_lm*(-(Vp + V/r)) - g*((A**2+J**2)*J - (V*A - Q*J)*Q)
    f_Q = r * (A * Jp - J * Ap)
    dI = f_Q
    return np.array([dV, dVp, dA, dAp, dQ, dQp, dJ, dJp, dI])


def bc(ya, yb, p):
    """DeepSeek's BCs verbatim. 10 residuals = 9 states + 1 free param."""
    k = omega
    return np.array([
        ya[1], ya[3], ya[5], ya[7],
        ya[8],
        yb[0]*k + yb[1],
        yb[2]*k + yb[3],
        yb[4]*k + yb[5],
        yb[6]*k + yb[7],
        yb[8] - 1.0,
    ])


print(f"Running DeepSeek reference (patched) — R_MIN={R_MIN}, Rmax={Rmax}, N={N}")
print(f"  omega = {omega}, g = {g}, lambda_lm = free (init {p0[0]})")

sol = solve_bvp(ode, bc, r, y0, p=p0, tol=1e-6, max_nodes=5000, verbose=2)
print(f"\nSuccess: {sol.success}, status: {sol.status}, message: {sol.message}")

if sol.success or sol.status in (0, 1):
    r_plot = sol.x
    V, Vp, A, Ap, Q, Qp, J, Jp, I = sol.y
    H = np.trapezoid(0.5*(Vp**2 + Ap**2 + Qp**2 + Jp**2) +
                     0.5*omega**2*(V**2 + A**2 + Q**2 + J**2) +
                     (-V*Q + A*J) +
                     (g/4)*((V**2+Q**2)**2 + (A**2+J**2)**2 + 2*(V*A - Q*J)**2),
                     x=r_plot)
    Q_CS = I[-1]
    Q_CS_grid = float(np.trapezoid(r_plot * (A * Jp - J * Ap), x=r_plot))
    print(f"\n  Final lambda_lm = {sol.p[0]:.6f}")
    print(f"  Q_CS (from I-state)  = {Q_CS:.6f}")
    print(f"  Q_CS (from grid)     = {Q_CS_grid:.6f}  (should match I-state)")
    print(f"  H                    = {H:.6f}")
    print(f"  H/Q_CS               = {H/Q_CS:.6f}")
    print(f"  Target H/Q           = 1.6969")
    print(f"  Ratio                = {(H/Q_CS)/1.6969:.4f}")
    print(f"  Final grid size      = {len(sol.x)}")
    print(f"  Peak V/A/Q/J         = {abs(V).max():.4f} / {abs(A).max():.4f} / "
          f"{abs(Q).max():.4f} / {abs(J).max():.4f}")
    print(f"  V(R_MIN), A, Q, J    = {V[0]:+.4f}, {A[0]:+.4f}, {Q[0]:+.4f}, {J[0]:+.4f}")
else:
    print(f"Solver failed: {sol.message}")
