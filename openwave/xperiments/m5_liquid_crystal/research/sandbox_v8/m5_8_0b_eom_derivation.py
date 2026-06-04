"""
M5.8.0b — step 0b-1/0b-2: symbolic derivation of the coupled (phi, psi) EOM for
Duda's 1+1D time-crystal toy model, and the Legendre-inverted accelerations the
leapfrog kernel needs.

Toy Lagrangian density (arXiv:2501.04036 Eq.1; 5a §10a), 1+1D, signature (+,-):

    L = phi_t^2 - phi_x^2 - (1 - phi^2)^2 - alpha * R^2 + (beta/3) * R^4
    R = phi_t*psi_x - phi_x*psi_t        (the curvature coupling; psi has NO
                                          standalone kinetic/potential — it enters
                                          ONLY through R)

Why this script exists (math-first, before any kernel):
  * psi couples to phi only through R, which mixes phi_t and psi_t, so the momenta
    are NON-CANONICAL -> the (phi_tt, psi_tt) accelerations are coupled by a
    field-dependent "mass matrix" H. A vanilla wave leapfrog is wrong; the stepper
    must invert H each step (the Legendre inversion).
  * SANITY CHECK (gate 0b-1): the conserved Hamiltonian density, restricted to the
    static-kink + psi=omega*t ansatz (phi_t=0, psi_x=0, psi_t=omega), MUST collapse
    to the M5.8.0a energy
        H = phi_x^2 (1 - alpha*omega^2) + (1 - phi^2)^2 + beta * omega^4 * phi_x^4
    whose minimum is omega* = 1.2897, E* = 2.0252 (5a §10a). If it does, the EOM is
    right and we can code the kernel (0b-3/0b-4).

Run:  python m5_8_0b_eom_derivation.py
"""

import sympy as sp

# ----------------------------------------------------------------------------
# 1. Jet variables (treat each derivative as an independent symbol so we can take
#    partials, then apply total-derivative operators by hand).
# ----------------------------------------------------------------------------
phi, phit, phix, phitt, phitx, phixx = sp.symbols(
    "phi phi_t phi_x phi_tt phi_tx phi_xx", real=True
)
psi, psit, psix, psitt, psitx, psixx = sp.symbols(
    "psi psi_t psi_x psi_tt psi_tx psi_xx", real=True
)
alpha, beta, omega, fx = sp.symbols("alpha beta omega phi_x_kink", positive=True)

R = phit * psix - phix * psit
L = phit**2 - phix**2 - (1 - phi**2) ** 2 - alpha * R**2 + (beta / 3) * R**4

# ----------------------------------------------------------------------------
# 2. Momenta (first derivatives of L w.r.t. the field velocities/gradients).
# ----------------------------------------------------------------------------
pL = {
    "phi_t": sp.diff(L, phit),
    "phi_x": sp.diff(L, phix),
    "psi_t": sp.diff(L, psit),
    "psi_x": sp.diff(L, psix),
}

# ----------------------------------------------------------------------------
# 3. Total-derivative operators Dt, Dx acting on jet expressions.
# ----------------------------------------------------------------------------
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

# Euler-Lagrange:  Dt(dL/dq_t) + Dx(dL/dq_x) - dL/dq = 0
EL_phi = sp.expand(Dt(pL["phi_t"]) + Dx(pL["phi_x"]) - sp.diff(L, phi))
EL_psi = sp.expand(Dt(pL["psi_t"]) + Dx(pL["psi_x"]) - sp.diff(L, psi))

# ----------------------------------------------------------------------------
# 4. The non-canonical mass matrix H (coefficients of the accelerations) and the
#    force vector b, so that  H @ [phi_tt, psi_tt] = b.
# ----------------------------------------------------------------------------
accel = [phitt, psitt]
H = sp.Matrix(
    [[sp.diff(EL_phi, a) for a in accel], [sp.diff(EL_psi, a) for a in accel]]
)
# everything that is NOT the H@accel part (move to RHS, flip sign):
rest_phi = sp.expand(EL_phi - (H[0, 0] * phitt + H[0, 1] * psitt))
rest_psi = sp.expand(EL_psi - (H[1, 0] * phitt + H[1, 1] * psitt))
b = sp.Matrix([-rest_phi, -rest_psi])

print("=" * 78)
print("EULER-LAGRANGE EQUATIONS  (general dynamical, = 0)")
print("=" * 78)
print("EL_phi =", EL_phi)
print()
print("EL_psi =", EL_psi)

print("\n" + "=" * 78)
print("MASS MATRIX H  (coeff of [phi_tt, psi_tt]) — the Legendre / non-canonical part")
print("=" * 78)
H = sp.simplify(H)
sp.pprint(H)
print("\noff-diagonal H[0,1] =", sp.simplify(H[0, 1]),
      "   <-- non-zero => phi_tt and psi_tt are COUPLED (the gotcha)")
print("det(H) =", sp.simplify(H.det()))

print("\n" + "=" * 78)
print("SOLVED ACCELERATIONS  [phi_tt, psi_tt] = H^{-1} b  (the kernel update)")
print("=" * 78)
acc = sp.simplify(H.inv() * b)
print("phi_tt =")
sp.pprint(sp.simplify(acc[0]))
print("\npsi_tt =")
sp.pprint(sp.simplify(acc[1]))

# ----------------------------------------------------------------------------
# 5. GATE 0b-1 SANITY CHECK: Hamiltonian density reduces to the 0a energy.
# ----------------------------------------------------------------------------
H_density = phit * pL["phi_t"] + psit * pL["psi_t"] - L  # Legendre transform
ansatz = {phit: 0, psix: 0, psit: omega, phix: fx}        # static kink + psi=omega*t
H_reduced = sp.simplify(H_density.subs(ansatz))
target = fx**2 * (1 - alpha * omega**2) + (1 - phi**2) ** 2 + beta * omega**4 * fx**4

print("\n" + "=" * 78)
print("SANITY CHECK (gate 0b-1): Hamiltonian under static-kink + psi=omega*t ansatz")
print("=" * 78)
print("H_reduced =", H_reduced)
print("target    =", sp.expand(target), "   (5a §10a Eq.3)")
ok = sp.simplify(H_reduced - target) == 0
print("\nMATCH 0a energy form?  ", "PASS ✓" if ok else "FAIL ✗")

# psi-EL must be a pure conservation law (dL/dpsi = 0 -> no psi source term):
psi_has_source = sp.diff(L, psi) != 0
print("psi-EL is a conservation law (dL/dpsi == 0)? ",
      "PASS ✓" if not psi_has_source else "FAIL ✗")

# Numeric minimum of E(omega) for the standard tanh kink (closed-form anchor):
#   integral phi_x^2 dx and phi_x^4 dx for phi=tanh(x/w), w chosen by 0a.
# Here just verify the *frequency-selection* formula omega^2 = a/(2 b I4/I2):
I2, I4 = sp.symbols("I2 I4", positive=True)  # ∫phi_x^2, ∫phi_x^4
E_of_omega = I2 * (1 - alpha * omega**2) + sp.Symbol("Vint", positive=True) + beta * omega**4 * I4
omega_star_sq = sp.solve(sp.diff(E_of_omega, omega), omega**2)
print("\nFrequency-selection (dE/domega = 0):  omega*^2 =", omega_star_sq, " = alpha*I2/(2*beta*I4)")
print("(5a §10a Eq.4 — the nonzero clock frequency = the time crystal)")
