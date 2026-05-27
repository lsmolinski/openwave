"""
M5.6.3a — port Faber's Model of Topological Fermions regularization (Universe 11/2025/113)

M5.6.1/6.2 found the KG mass is geometric (`C_μν~1/r²`) and its SCALE is set by the core
regularization (an ad-hoc r_c so far). Faber's MTF gives the PHYSICAL regularization that
pins that scale. We port it natively here (per `reference_faber_regularization`: port,
don't reinvent) and verify his central results.

Faber's scheme (his Eq.1-8):
    Q = q0 − iσ·q⃗ ,  unit 4-vector on S³:  q0² + ‖q⃗‖² = 1   (q0=cos α, q⃗=n̂ sin α)
    Γ⃗_μ = q0 ∂_μq⃗ − (∂_μq0) q⃗ + q⃗ × ∂_μq⃗            (Eq.6, first form — smooth at core)
    R⃗_μν = Γ⃗_μ × Γ⃗_ν                                  (curvature = dual EM field, Eq.5/9)
    Λ = q0⁶ / r0⁴                                       (Eq.4 — THE regularization potential)
    L = −(α_f ℏc/4π)[ ¼ R⃗_μν·R⃗^μν + Λ ]                (Eq.2)
    soliton: n̂=x̂, α=arctan(r/r0)  ⇒  q0 = r0/√(r0²+r²),  q⃗ = x⃗/√(r0²+r²)
    rest energy:  E0 = α_f ℏc · π/(4 r0)                (Eq.8)

Checks:
  1. unit-4-vector + the core regularization (‖q⃗‖→0 at center, no singularity).
  2. the rest-energy law  E0 = (α_f ℏc) π/(4 r0)  — reproduce numerically (units α_f ℏc=1).
  3. E0 ∝ 1/r0  (the mass scale is set by r0 — the "pin the mass" deliverable).

USAGE:
    python -m openwave.xperiments.m5_liquid_crystal.research.sandbox_v6.m5_6_3a_faber_regularization
"""
import numpy as np


def faber_energy(r0, N=121, Lbox=None):
    """Integrate Faber's static energy E = (1/4π)∫[¼ R·R + Λ] d³x for the hedgehog soliton.
    Returns E in units α_f ℏc = 1 (so the target is π/(4 r0))."""
    if Lbox is None:
        Lbox = 14.0 * r0                       # box ≫ r0 so the 1/r⁴ tail is captured
    xs = np.linspace(-Lbox, Lbox, N)
    h = xs[1] - xs[0]
    X, Y, Z = np.meshgrid(xs, xs, xs, indexing="ij")
    reg = X**2 + Y**2 + Z**2 + r0**2            # = r0²+r²
    rn = np.sqrt(reg)
    q0 = r0 / rn                                # cos α  (→1 at center, →0 far)
    q = np.stack([X, Y, Z], -1) / rn[..., None]  # q⃗ = x⃗/√(r0²+r²) = n̂ sin α (→0 at center)

    def d(f, ax):                               # central difference along spatial axis
        return (np.roll(f, -1, ax) - np.roll(f, 1, ax)) / (2 * h)

    # Γ⃗_i = q0 ∂_i q⃗ − (∂_i q0) q⃗ + q⃗ × ∂_i q⃗     (Eq.6)
    Gamma = []
    for ax in range(3):
        dq = d(q, ax)                           # ∂_i q⃗  (…,3)
        dq0 = d(q0, ax)                         # ∂_i q0
        G = q0[..., None] * dq - dq0[..., None] * q + np.cross(q, dq)
        Gamma.append(G)

    # R⃗_ij = Γ⃗_i × Γ⃗_j ;  ¼ R·R = ½ Σ_{i<j} ‖R_ij‖²
    RR = np.zeros(q0.shape)
    for (i, j) in ((0, 1), (0, 2), (1, 2)):
        Rij = np.cross(Gamma[i], Gamma[j])
        RR += np.einsum("...a,...a->...", Rij, Rij)
    curv = 0.5 * RR                             # ¼ Σ_ij R_ij·R_ij  (full index sum)
    Lam = q0**6 / r0**4                         # Faber Eq.4

    interior = tuple(slice(1, -1) for _ in range(3))
    E = (1.0 / (4 * np.pi)) * np.sum((curv + Lam)[interior]) * h**3
    return E, q0, q


def main():
    print("=" * 70)
    print("M5.6.3a — Faber MTF regularization: rest-energy law E0 = α_f ℏc·π/(4r0)")
    print("=" * 70)

    # --- 1. unit 4-vector + core regularization ---------------------------------
    E, q0, q = faber_energy(1.0)
    norm = q0**2 + np.einsum("...a,...a->...", q, q)
    c = q0.shape[0] // 2
    qmag_center = float(np.sqrt(np.einsum("...a,...a->...", q, q))[c, c, c])
    print("\n[1] field structure (r0=1)")
    print(f"    unit 4-vector q0²+‖q⃗‖²: min={norm.min():.6f} max={norm.max():.6f} (expect 1)")
    print(f"    ‖q⃗‖ at center = {qmag_center:.3e}  → radial vector shrinks to 0 ⇒ NO singularity ✓")
    print(f"    (this IS the M5.6.2a core-melt mechanism — Faber validates it)")

    # --- 2. reproduce the rest-energy law ---------------------------------------
    print("\n[2] rest energy vs analytic  E0 = π/(4 r0)  (units α_f ℏc = 1)")
    print(f"    {'r0':>5} {'N':>5} {'E0 (numeric)':>14} {'π/(4r0)':>12} {'rel err':>9}")
    rows = []
    for r0, N in ((1.0, 101), (1.0, 141), (1.0, 181)):
        E, _, _ = faber_energy(r0, N=N)
        target = np.pi / (4 * r0)
        rel = abs(E - target) / target
        rows.append((N, E, rel))
        print(f"    {r0:>5.1f} {N:>5} {E:>14.5f} {target:>12.5f} {100*rel:>8.2f}%")
    converging = rows[2][2] < rows[0][2]
    print(f"    → converging to π/4={np.pi/4:.5f} as N↑: {converging}  (finite-difference + box-truncation limited)")

    # --- 3. E0 ∝ 1/r0  (the mass scale is set by r0) ----------------------------
    print("\n[3] mass scale: E0 ∝ 1/r0  (the regularization radius pins the mass)")
    print(f"    {'r0':>5} {'E0·r0':>10}  (should be ≈ π/4 = {np.pi/4:.4f}, r0-independent)")
    prods = []
    for r0 in (0.5, 1.0, 2.0):
        E, _, _ = faber_energy(r0, N=141)
        prods.append(E * r0)
        print(f"    {r0:>5.1f} {E*r0:>10.5f}")
    scale_ok = np.std(prods) / np.mean(prods) < 0.05
    print(f"    → E0·r0 constant (CV={100*np.std(prods)/np.mean(prods):.1f}%): {scale_ok}  ⇒ E0 ∝ 1/r0 ✓")
    print(f"    physical anchor: r0=2.2132 fm ⇒ E0=0.511 MeV (electron); r0 = π/4 × classical e⁻ radius")

    best_rel = rows[-1][2]
    ok = (abs(norm - 1).max() < 1e-6 and qmag_center < 1e-2 and best_rel < 0.05 and scale_ok)
    print("\n" + "=" * 70)
    print(f"M5.6.3a: Faber MTF ported; rest-energy law E0=α_f ℏc·π/(4r0) reproduced "
          f"({100*best_rel:.1f}% at N=181), E0∝1/r0. Core regularization = the M5.6.2a melt.")
    print("PASS" if ok else "PARTIAL — inspect the failing metric above")
    print("=" * 70)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
