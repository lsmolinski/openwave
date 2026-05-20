"""DeepSeek's reference Python script — received 2026-05-20 ~5:30 PM via Paul.
Saved AS-IS for diff against our v6. Modifications documented in companion
files (deepseek_reference_patched.py if patches are needed to run).
"""

import numpy as np
from scipy.integrate import solve_bvp

# Physical parameters (calibrated for electron)
omega = 1.0          # frequency
g = 1.0625           # quartic coupling
lambda_lm = 1.0      # Lagrange multiplier (initial guess)

# Radial grid
Rmax = 20.0
N = 500
r = np.linspace(0, Rmax, N)

# Initial guess (exponential decay, J decays faster to break proportionality)
V0 = 0.1 * np.exp(-r)
A0 = -0.1 * np.exp(-r)
Q0 = 0.1 * np.exp(-r)
J0 = -0.1 * np.exp(-1.5 * r)

# Derivatives for initial guess (approximate)
Vp0 = -V0
Ap0 = -A0
Qp0 = -Q0
Jp0 = 1.5 * (-J0)   # because d/dr exp(-1.5r) = -1.5 exp(-1.5r)

# Initial state vector: y = [V, V', A, A', Q, Q', J, J', I]
y0 = np.vstack([V0, Vp0, A0, Ap0, Q0, Qp0, J0, Jp0, np.zeros_like(r)])


# ODE system with auxiliary integral I(r) = ∫_0^r f(r') dr'
def ode(r, y):
    V, Vp, A, Ap, Q, Qp, J, Jp, I = y
    # Note: simplified quartic terms; full form may need sign adjustments.
    # The key is the integral constraint.
    dV = Vp
    dVp = -(omega**2)*V + Q - lambda_lm*(Jp + J/r) - g*( (V**2+Q**2)*V + (V*A - Q*J)*A )
    dA = Ap
    dAp = -(omega**2)*A - J + lambda_lm*( - (Vp + V/r) ) - g*( (A**2+J**2)*A + (V*A - Q*J)*V )
    dQ = Qp
    dQp = -(omega**2)*Q + V - lambda_lm*(Jp + J/r) - g*( (V**2+Q**2)*Q - (V*A - Q*J)*J )
    dJ = Jp
    dJp = -(omega**2)*J - A + lambda_lm*( -(Vp + V/r) ) - g*( (A**2+J**2)*J - (V*A - Q*J)*Q )
    # Integrand for Q_CS (mutual linking density)
    f_Q = r * (A * Jp - J * Ap)
    dI = f_Q
    return np.array([dV, dVp, dA, dAp, dQ, dQp, dJ, dJp, dI])


# Boundary conditions: regularity at r=0, Robin decay at r=Rmax, integral constraint
def bc(ya, yb):
    k = omega   # decay constant approx
    return np.array([
        ya[1], ya[3], ya[5], ya[7],    # V'(0)=0, A'(0)=0, Q'(0)=0, J'(0)=0
        ya[8],                         # I(0)=0
        yb[0]*k + yb[1],               # V'(Rmax) + k V(Rmax) = 0
        yb[2]*k + yb[3],               # A'(Rmax) + k A(Rmax) = 0
        yb[4]*k + yb[5],               # Q'(Rmax) + k Q(Rmax) = 0
        yb[6]*k + yb[7],               # J'(Rmax) + k J(Rmax) = 0
        yb[8] - 1.0                    # I(Rmax) = 1   (integral constraint)
    ])


# Solve BVP
sol = solve_bvp(ode, bc, r, y0, tol=1e-6, max_nodes=5000, verbose=2)
print(f"Success: {sol.success}, message: {sol.message}")

if sol.success:
    r_plot = sol.x
    V, Vp, A, Ap, Q, Qp, J, Jp, I = sol.y
    # Compute energy H (Hamiltonian integral)
    H = np.trapz(0.5*(Vp**2 + Ap**2 + Qp**2 + Jp**2) +
                 0.5*omega**2*(V**2 + A**2 + Q**2 + J**2) +
                 (-V*Q + A*J) +
                 (g/4)*((V**2+Q**2)**2 + (A**2+J**2)**2 + 2*(V*A - Q*J)**2),
                 x=r_plot)
    Q_CS = I[-1]
    print(f"H = {H:.6f}, Q_CS = {Q_CS:.6f}, H/Q = {H/Q_CS:.6f}")
    print(f"Target H/Q = 1.6969, ratio = {(H/Q_CS)/1.6969}")
else:
    print("Solver failed; adjust parameters.")
