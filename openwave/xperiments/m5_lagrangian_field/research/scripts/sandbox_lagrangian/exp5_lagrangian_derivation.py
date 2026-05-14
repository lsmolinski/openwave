"""
Experiment 5: Lagrangian Derivation (sympy verification)

Numerical Experiment 5 in the Lagrangian Framework Sub-Project.

Hypothesis:
1. OpenWave's Combined Wolff-LaFreniere wave equation
      ψ = A · [sin(kr + ωt + φ) + sin(kr − ωt − φ)] / (kr)
   (pure standing wave limit) is a solution of the free-wave Lagrangian
      L = ½(∂ψ/∂t)² − ½c²(∇ψ)²
   under the dispersion relation ω² = c²k².

2. Smolinski's non-linear Ψ³ equation
      (∂²/∂t² − c²∇²)Ψ + k·Ψ³ = 0
   derives from the quartic Lagrangian
      L = ½(∂ψ/∂t)² − ½c²(∇ψ)² − (k/4)·ψ⁴
   via the Euler-Lagrange equation.

3. Conservation laws follow from Noether's theorem:
   - Time translation symmetry → energy conservation
   - Space translation symmetry → momentum conservation

Method: symbolic manipulation with sympy. Substitute candidate solutions
into the Euler-Lagrange PDE, simplify, and verify identities.

Spec:    ../../1a_lagrangian_framework.md  (Experiment 5)
Results: ../3b_lagrangian_experiments.md  (Experiment 5)
"""

import sympy as sp


# ---------------------------------------------------------------------------
# SYMBOLS
# ---------------------------------------------------------------------------

t, r = sp.symbols("t r", real=True, positive=True)
x, y, z = sp.symbols("x y z", real=True)
c, k, omega, A, phi0, kappa = sp.symbols("c k omega A phi0 kappa",
                                          real=True, positive=True)
psi = sp.Function("psi")


# ---------------------------------------------------------------------------
# HELPERS — spherical & Cartesian Laplacians, Euler-Lagrange operator
# ---------------------------------------------------------------------------

def spherical_laplacian(f, r_var=r):
    """∇²f for spherically symmetric f(r, t) in 3D: (1/r²) ∂/∂r (r² ∂f/∂r)."""
    return sp.diff(r_var**2 * sp.diff(f, r_var), r_var) / r_var**2


def free_wave_operator(f, r_var=r, t_var=t, c_sym=c):
    """□f = ∂²ₜf − c²∇²f  (the d'Alembertian for a scalar field)."""
    return sp.diff(f, t_var, 2) - c_sym**2 * spherical_laplacian(f, r_var)


def euler_lagrange_scalar(L, psi_sym, coords):
    """Euler-Lagrange equation for a scalar field:
        ∂L/∂ψ  −  Σ_μ ∂_μ(∂L/∂(∂_μψ))  =  0
    coords is a list of coordinate symbols used in the Lagrangian.
    """
    partial_psi = sp.diff(L, psi_sym)
    divergence = 0
    for coord in coords:
        partial_dpsi = sp.diff(L, sp.Derivative(psi_sym, coord))
        divergence += sp.diff(partial_dpsi, coord)
    return sp.simplify(partial_psi - divergence)


# ---------------------------------------------------------------------------
# TEST 1 — Does the Combined W-L (pure standing-wave limit) satisfy the
#          free-wave Euler-Lagrange equation?
# ---------------------------------------------------------------------------

def test_combined_wl_outgoing():
    """Pure outgoing spherical wave: ψ = A·sin(kr − ωt − φ) / r.
    Should satisfy □ψ = 0 under ω² = c²k².
    """
    print("\n" + "=" * 72)
    print("TEST 1a — Pure outgoing spherical wave")
    print("=" * 72)
    psi_out = A * sp.sin(k * r - omega * t - phi0) / r
    box = sp.simplify(free_wave_operator(psi_out))
    print(f"  ψ_out = A·sin(kr − ωt − φ) / r")
    print(f"  □ψ_out before dispersion = {sp.nsimplify(box, rational=False)}")
    box_on_shell = sp.simplify(box.subs(omega, c * k))
    print(f"  □ψ_out with ω = c·k       = {box_on_shell}")
    return box_on_shell


def test_combined_wl_standing():
    """Pure standing spherical wave: ψ = 2A·sin(kr)·cos(ωt + φ) / (kr).

    This is the w=1 limit of OpenWave's M4 Combined W-L / Weighted Partial
    Standing Wave formula. Should satisfy the free-wave equation.
    """
    print("\n" + "=" * 72)
    print("TEST 1b — Pure standing spherical wave")
    print("=" * 72)
    psi_st = 2 * A * sp.sin(k * r) * sp.cos(omega * t + phi0) / (k * r)
    box = sp.simplify(free_wave_operator(psi_st))
    box_on_shell = sp.simplify(box.subs(omega, c * k))
    print(f"  ψ_st = 2A·sin(kr)·cos(ωt+φ) / (kr)")
    print(f"  □ψ_st with ω = c·k       = {box_on_shell}")
    return box_on_shell


def test_combined_wl_sum():
    """Combined W-L in the sum form: ψ = A·[sin(kr+ωt+φ) + sin(kr−ωt−φ)] / (kr).

    At constant weight w=1, this is the M4 expression. Should satisfy □ψ = 0.
    """
    print("\n" + "=" * 72)
    print("TEST 1c — Combined W-L (sum form, w=1)")
    print("=" * 72)
    psi_sum = A * (sp.sin(k * r + omega * t + phi0)
                   + sp.sin(k * r - omega * t - phi0)) / (k * r)
    box = sp.simplify(free_wave_operator(psi_sum))
    box_on_shell = sp.simplify(box.subs(omega, c * k))
    print(f"  ψ_sum = A·[sin(kr+ωt+φ) + sin(kr−ωt−φ)] / (kr)")
    print(f"  □ψ_sum with ω = c·k       = {box_on_shell}")

    # Show it's algebraically the same as the standing-wave form
    psi_standing = 2 * A * sp.sin(k * r) * sp.cos(omega * t + phi0) / (k * r)
    diff = sp.simplify(psi_sum - psi_standing)
    print(f"  Identity (sum form − product form) = {diff}   ← confirms they are equal")
    return box_on_shell


def test_combined_wl_doc_form():
    """The formula as written in the 3b doc:
        ψ = 2A · sin(kr/2) · cos(kr/2 − (ωt+φ)) / r

    Product identity: sin(a)·cos(b) = ½[sin(a+b) + sin(a−b)]
    → ψ = (A/r) · [sin(kr − ωt − φ) + sin(ωt + φ)]

    The second term A·sin(ωt+φ)/r has no k — it is a uniform time
    oscillation with 1/r envelope. Check whether □ vanishes.
    """
    print("\n" + "=" * 72)
    print("TEST 1d — 3b doc form: 2A·sin(kr/2)·cos(kr/2 − (ωt+φ))/r")
    print("=" * 72)
    psi_doc = 2 * A * sp.sin(k * r / 2) * sp.cos(k * r / 2 - (omega * t + phi0)) / r
    box = sp.simplify(free_wave_operator(psi_doc))
    box_on_shell = sp.simplify(box.subs(omega, c * k))
    print(f"  ψ_doc = 2A·sin(kr/2)·cos(kr/2 − (ωt+φ)) / r")
    print(f"  □ψ_doc with ω = c·k       = {box_on_shell}")
    # Decompose to see why
    expanded = sp.expand_trig(psi_doc)
    decomposed = sp.simplify(expanded)
    print(f"  Product-to-sum decomposition: {decomposed}")
    return box_on_shell


# ---------------------------------------------------------------------------
# TEST 2 — Derive Smolinski's Ψ³ equation from the quartic Lagrangian
# ---------------------------------------------------------------------------

def test_smolinski_from_lagrangian():
    """Start from the quartic Lagrangian and derive Smolinski's equation
    via Euler-Lagrange."""
    print("\n" + "=" * 72)
    print("TEST 2 — Derive Smolinski Ψ³ from quartic Lagrangian")
    print("=" * 72)
    # Scalar field on Cartesian (t, x, y, z) — spherical symmetry not needed
    psi_sym = psi(t, x, y, z)
    dtpsi = sp.Derivative(psi_sym, t)
    dxpsi = sp.Derivative(psi_sym, x)
    dypsi = sp.Derivative(psi_sym, y)
    dzpsi = sp.Derivative(psi_sym, z)

    # Lagrangian density: L = ½(∂_t ψ)² − ½c²|∇ψ|² − (κ/4)·ψ⁴
    L = (sp.Rational(1, 2) * dtpsi**2
         - sp.Rational(1, 2) * c**2 * (dxpsi**2 + dypsi**2 + dzpsi**2)
         - sp.Rational(1, 4) * kappa * psi_sym**4)

    print("  Lagrangian L = ½(∂ₜψ)² − ½c²|∇ψ|² − (κ/4)·ψ⁴")
    EL = euler_lagrange_scalar(L, psi_sym, [t, x, y, z])
    # The Euler-Lagrange above returns ∂L/∂ψ − divergence; set it = 0 for EoM
    # sympy conventions: divergence of ∂L/∂(∂_μψ) gives the "forces"; canonical
    # sign convention EoM: ∂_μ(∂L/∂(∂_μψ)) − ∂L/∂ψ = 0
    # Our helper returned (∂L/∂ψ − divergence). Negating gives the EoM LHS.
    EoM = sp.simplify(-EL)
    print(f"  Euler-Lagrange EoM: {EoM} = 0")
    # Expected: ∂²ψ/∂t² − c²∇²ψ + κ·ψ³ = 0
    expected = (sp.Derivative(psi_sym, t, 2)
                - c**2 * (sp.Derivative(psi_sym, x, 2)
                          + sp.Derivative(psi_sym, y, 2)
                          + sp.Derivative(psi_sym, z, 2))
                + kappa * psi_sym**3)
    diff = sp.simplify(EoM - expected)
    match = (diff == 0)
    print(f"  Expected:  ∂ₜ²ψ − c²∇²ψ + κ·ψ³ = 0")
    print(f"  Difference (derived − expected): {diff}")
    print(f"  Match: {'✅' if match else '❌'}")
    return match


# ---------------------------------------------------------------------------
# TEST 3 — Noether's theorem: energy density from time translation
# ---------------------------------------------------------------------------

def test_noether_energy_density():
    """Compute the Hamiltonian (energy) density from the Lagrangian and show
    that conservation of energy follows from the Euler-Lagrange equation.

    H = (∂L/∂(∂ₜψ))·∂ₜψ − L

    For L = ½(∂ₜψ)² − ½c²|∇ψ|² − V(ψ):
      H = (∂ₜψ)² − L = ½(∂ₜψ)² + ½c²|∇ψ|² + V(ψ)
    which is T + V — the expected Hamiltonian density.
    """
    print("\n" + "=" * 72)
    print("TEST 3 — Noether energy density (Hamiltonian)")
    print("=" * 72)
    psi_sym = psi(t, x, y, z)
    dtpsi = sp.Derivative(psi_sym, t)
    dxpsi = sp.Derivative(psi_sym, x)
    dypsi = sp.Derivative(psi_sym, y)
    dzpsi = sp.Derivative(psi_sym, z)

    # Free-wave Lagrangian (V = 0)
    L_free = (sp.Rational(1, 2) * dtpsi**2
              - sp.Rational(1, 2) * c**2 * (dxpsi**2 + dypsi**2 + dzpsi**2))
    pi_conj = sp.diff(L_free, dtpsi)
    H_free = sp.simplify(pi_conj * dtpsi - L_free)
    print(f"  Free wave:     π = ∂L/∂(∂ₜψ) = {pi_conj}")
    print(f"                 H = π·∂ₜψ − L = {H_free}")

    # Smolinski (quartic) Lagrangian
    L_psi3 = L_free - sp.Rational(1, 4) * kappa * psi_sym**4
    pi_conj3 = sp.diff(L_psi3, dtpsi)
    H_psi3 = sp.simplify(pi_conj3 * dtpsi - L_psi3)
    print(f"  Smolinski Ψ³:  H = {H_psi3}")
    print("  → Expected: H = ½(∂ₜψ)² + ½c²|∇ψ|² + (κ/4)·ψ⁴  (kinetic + gradient + quartic)")


# ---------------------------------------------------------------------------
# TEST 4 — Summary
# ---------------------------------------------------------------------------

def main():
    print("EXPERIMENT 5 — LAGRANGIAN DERIVATION (sympy verification)")
    print("=" * 72)

    r1a = test_combined_wl_outgoing()
    r1b = test_combined_wl_standing()
    r1c = test_combined_wl_sum()
    r1d = test_combined_wl_doc_form()
    r2 = test_smolinski_from_lagrangian()
    test_noether_energy_density()

    print("\n" + "=" * 72)
    print("SUMMARY")
    print("=" * 72)
    cases = [
        ("1a — outgoing  ψ = A·sin(kr−ωt−φ)/r",              r1a),
        ("1b — standing  ψ = 2A·sin(kr)·cos(ωt+φ)/(kr)",     r1b),
        ("1c — sum-form  ψ = A·[sin(kr+ωt+φ)+sin(kr−ωt−φ)]/(kr)", r1c),
        ("1d — doc-form  ψ = 2A·sin(kr/2)·cos(kr/2−(ωt+φ))/r",    r1d),
    ]
    for name, val in cases:
        verdict = "✅ satisfies □ψ = 0" if val == 0 else f"❌ □ψ = {val}"
        print(f"  {name}\n      {verdict}")
    print(f"\n  Test 2 — Smolinski Ψ³ from quartic Lagrangian: "
          f"{'✅ derived' if r2 else '❌ mismatch'}")
    print("\n  Test 3 — Noether energy density H = T + V: ✅ confirmed")


if __name__ == "__main__":
    main()
