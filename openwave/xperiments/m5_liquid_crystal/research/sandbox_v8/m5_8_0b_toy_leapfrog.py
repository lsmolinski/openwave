"""
M5.8.0b-3/0b-4 — 1+1D field leapfrog for Duda's time-crystal toy model.

The dynamical pre-check (the field-stepper dry-run before the 4×4 promotion). Uses
the sympy-derived non-canonical EOM from m5_8_0b_eom_derivation (single source of
truth): build H (velocity-Hessian / mass matrix) and b, then per grid point solve
H @ [phi_tt, psi_tt] = b.

Two known subtleties from the derivation, both handled here:
  * det(H) ∝ phi_x^2  -> psi is inertia-less in the vacuum (the 1D shadow of "the
    clock lives on the hedgehog"). Regularize with a small standard psi-kinetic
    eps*(psi_t^2 - psi_x^2); eps -> 0 recovers the pure toy.
  * H[1,1] = -2 alpha phi_x^2 < 0 -> the mass matrix is INDEFINITE on the kink core
    (the negative-energy mechanism). det(H) can change sign -> a singular ring; we
    floor |det| to keep the inversion finite and watch for the pathology honestly.

Tests (gate M5.8.0b):
  (A) ROBUST residual + energy check (no long-time stability needed): seed the
      static-kink + psi=omega*t ansatz; the EOM residual ||accel|| should be small
      at omega* and the energy E(omega*) < E(0).
  (B) DYNAMICAL run: evolve from the ansatz with RK4; report whether the clock holds
      at omega ≈ 1.2898 and energy is conserved, or whether it hits the indefinite-H
      pathology (an honest, documentable M5.8.2 finding either way).

Run:  python m5_8_0b_toy_leapfrog.py
"""

import numpy as np
import sympy as sp

# ----------------------------------------------------------------------------
# 1. Symbolic H, b (regularized) -> numpy callables. Reuses the derivation logic.
# ----------------------------------------------------------------------------
phi, phit, phix, phitt, phitx, phixx = sp.symbols(
    "phi phit phix phitt phitx phixx", real=True
)
psi, psit, psix, psitt, psitx, psixx = sp.symbols(
    "psi psit psix psitt psitx psixx", real=True
)
A, B, E = sp.symbols("A B E", real=True)  # alpha, beta, eps

R = phit * psix - phix * psit
L = (
    phit**2 - phix**2 - (1 - phi**2) ** 2
    - A * R**2 + (B / 3) * R**4
    + E * (psit**2 - psix**2)               # <-- regularizing psi-kinetic
)

def Dt(e):
    return (
        sp.diff(e, phi) * phit + sp.diff(e, phit) * phitt + sp.diff(e, phix) * phitx
        + sp.diff(e, psi) * psit + sp.diff(e, psit) * psitt + sp.diff(e, psix) * psitx
    )

def Dx(e):
    return (
        sp.diff(e, phi) * phix + sp.diff(e, phit) * phitx + sp.diff(e, phix) * phixx
        + sp.diff(e, psi) * psix + sp.diff(e, psit) * psitx + sp.diff(e, psix) * psixx
    )

EL_phi = sp.expand(Dt(sp.diff(L, phit)) + Dx(sp.diff(L, phix)) - sp.diff(L, phi))
EL_psi = sp.expand(Dt(sp.diff(L, psit)) + Dx(sp.diff(L, psix)) - sp.diff(L, psi))

H00 = sp.diff(EL_phi, phitt); H01 = sp.diff(EL_phi, psitt)
H10 = sp.diff(EL_psi, phitt); H11 = sp.diff(EL_psi, psitt)
b0 = -sp.expand(EL_phi - H00 * phitt - H01 * psitt)
b1 = -sp.expand(EL_psi - H10 * phitt - H11 * psitt)

# Hamiltonian density (Legendre transform, incl. eps term)
Hdens = phit * sp.diff(L, phit) + psit * sp.diff(L, psit) - L

ARGS = (phi, phit, phix, phixx, phitx, psit, psix, psixx, psitx, A, B, E)
f_H00 = sp.lambdify(ARGS, H00, "numpy")
f_H01 = sp.lambdify(ARGS, H01, "numpy")
f_H11 = sp.lambdify(ARGS, H11, "numpy")
f_b0 = sp.lambdify(ARGS, b0, "numpy")
f_b1 = sp.lambdify(ARGS, b1, "numpy")
f_Hd = sp.lambdify(ARGS, Hdens, "numpy")

# ----------------------------------------------------------------------------
# 2. Grid + finite differences (non-periodic; np.gradient one-sided at edges).
# ----------------------------------------------------------------------------
Lbox, N = 12.0, 400
xg = np.linspace(-Lbox, Lbox, N)
dx = xg[1] - xg[0]
ALPHA, BETA, EPS = 1.0, 1.0, 0.0   # eps=0 = pure toy (psi auto-localized: R=0 in vacuum)
DETFLOOR = 1e-3

def ddx(u):
    return np.gradient(u, dx, edge_order=2)

def d2dx2(u):
    return np.gradient(np.gradient(u, dx, edge_order=2), dx, edge_order=2)

def accel(phiv, psiv, vphi, vpsi):
    """Return (a_phi, a_psi) by solving H @ a = b per grid point."""
    px, pxx = ddx(phiv), d2dx2(phiv)
    sx, sxx = ddx(psiv), d2dx2(psiv)
    ptx, stx = ddx(vphi), ddx(vpsi)
    args = (phiv, vphi, px, pxx, ptx, vpsi, sx, sxx, stx, ALPHA, BETA, EPS)
    H00v, H01v, H11v = f_H00(*args), f_H01(*args), f_H11(*args)
    b0v, b1v = f_b0(*args), f_b1(*args)
    # broadcast scalars to arrays
    H00v, H01v, H11v = (np.broadcast_to(v, phiv.shape).astype(float) for v in (H00v, H01v, H11v))
    b0v, b1v = (np.broadcast_to(v, phiv.shape).astype(float) for v in (b0v, b1v))
    det = H00v * H11v - H01v**2
    det = np.where(np.abs(det) < DETFLOOR, np.sign(det + 1e-30) * DETFLOOR, det)
    a_phi = (H11v * b0v - H01v * b1v) / det
    a_psi = (-H01v * b0v + H00v * b1v) / det
    return a_phi, a_psi

def energy(phiv, psiv, vphi, vpsi):
    px = ddx(phiv); sx = ddx(psiv)
    args = (phiv, vphi, px, 0 * px, 0 * px, vpsi, sx, 0 * px, 0 * px, ALPHA, BETA, EPS)
    trapz = getattr(np, "trapezoid", getattr(np, "trapz", None))
    return trapz(f_Hd(*args), xg)

# ----------------------------------------------------------------------------
# 3. Seed: static kink + psi = omega*t  (so psi_x = 0, psi_t = omega initially).
# ----------------------------------------------------------------------------
W = np.sqrt(96 / 61)            # standard-tanh kink width (5a §10a Eq.5)
OMEGA_STAR = 1.2898             # the time-crystal clock (5a §10a)

def seed(omega):
    phiv = np.tanh(xg / W)
    psiv = np.zeros_like(xg)     # psi grows in time; only derivatives matter
    vphi = np.zeros_like(xg)
    vpsi = omega * np.ones_like(xg)
    return phiv, psiv, vphi, vpsi

# ============================================================================
# TEST A — robust residual + energy check (no time-stepping)
# ============================================================================
print("=" * 70)
print("TEST A — static-ansatz residual + energy (robust; gate 0b-(i),(iii))")
print("=" * 70)
core = np.abs(xg) < 2.0
for om in (OMEGA_STAR, 0.0):
    p, s, vp, vs = seed(om)
    a_phi, a_psi = accel(p, s, vp, vs)
    res = np.sqrt(np.mean(a_phi[core] ** 2 + a_psi[core] ** 2))
    Etot = energy(p, s, vp, vs)
    print(f"  omega={om:6.4f}:  core EOM residual ||accel|| = {res:8.4e}   E = {Etot:8.5f}")
E_star = energy(*seed(OMEGA_STAR)); E_zero = energy(*seed(0.0))
print(f"\n  E(omega*) < E(0)?  {E_star:.5f} < {E_zero:.5f}  -> "
      f"{'PASS ✓ (time crystal: oscillating wins)' if E_star < E_zero else 'FAIL ✗'}")

# ============================================================================
# TEST B — dynamical RK4 leapfrog from the ansatz
# ============================================================================
print("\n" + "=" * 70)
print("TEST B — dynamical RK4 evolution (gate 0b-(i) energy, 0b-(ii) clock holds)")
print("=" * 70)

def rhs(state):
    p, s, vp, vs = state
    a_phi, a_psi = accel(p, s, vp, vs)
    return np.array([vp, vs, a_phi, a_psi])

def rk4_step(state, dt):
    k1 = rhs(state)
    k2 = rhs(state + 0.5 * dt * k1)
    k3 = rhs(state + 0.5 * dt * k2)
    k4 = rhs(state + dt * k1 * 0 + dt * k3)  # standard RK4
    return state + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)

dt = 0.2 * dx          # CFL-ish
nsteps = 600
state = np.array(seed(OMEGA_STAR))
E0 = energy(*state)
ic = np.argmin(np.abs(xg))     # core index
psi_core_hist = [state[1][ic]]
t = 0.0
blew_up = False
for n in range(nsteps):
    state = rk4_step(state, dt)
    t += dt
    psi_core_hist.append(state[1][ic])
    if not np.all(np.isfinite(state)) or np.max(np.abs(state[0])) > 5:
        print(f"  ⚠ blow-up / pathology at step {n} (t={t:.3f}) — indefinite-H singular ring")
        blew_up = True
        break

if not blew_up:
    Ef = energy(*state)
    drift = abs(Ef - E0) / abs(E0)
    # measure omega from psi at the core (slope of psi(t))
    psi_core_hist = np.array(psi_core_hist)
    tg = np.arange(len(psi_core_hist)) * dt
    omega_meas = np.polyfit(tg, psi_core_hist, 1)[0]
    print(f"  ran {nsteps} steps to t={t:.3f}")
    print(f"  energy drift |dE/E|        = {drift:.3e}   "
          f"({'PASS ✓' if drift < 0.05 else 'high — see notes'})")
    print(f"  measured clock omega(core) = {omega_meas:.4f}   (target {OMEGA_STAR})   "
          f"({'PASS ✓' if abs(omega_meas - OMEGA_STAR) / OMEGA_STAR < 0.1 else 'off'})")
    print(f"  max|phi|                   = {np.max(np.abs(state[0])):.4f}  (kink intact if ≈1)")
