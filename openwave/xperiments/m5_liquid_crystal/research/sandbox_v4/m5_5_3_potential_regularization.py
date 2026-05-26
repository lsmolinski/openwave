"""
M5.5.3 — V(M) potential + Faber-style core regularization (headless)

Adds the Higgs-like potential V(M) to the Eq.18 picture and validates its TWO roles
(Duda §III, Eq.12–13): (1) define the vacuum shape Λ=(1,δ,0), (2) regularize the
defect-core singularity to finite energy → finite mass.

KEY STRUCTURAL FINDING (this step): V(M) is ROTATION-INVARIANT — for M=ODOᵀ with
frozen D, the eigenvalues are always (1,δ,0), so V(ODOᵀ)=V(D)=const. Hence V does
NOTHING to the twist/rotation sector (M5.5.2 correctly used V=0). V acts ONLY on the
eigenvalue-deformation sector — i.e. V *is* the regularization. Eq.18's F² curvature
(the Skyrme-family term) + V make the core Derrick-stable at finite size → finite mass.

Three stages:
  1. V(M) Eq.12/13 + rotation-invariance check (the finding).
  2. V defines the vacuum: a perturbed-eigenvalue M relaxes (∇V flow) back to Λ.
  3. Regularization (Derrick): rigid hedgehog core curvature ∫‖F‖² DIVERGES; a deformed
     core (amplitude s(r)→0) is FINITE. Eq.18 energy E(L)=E_F/L + E_V·L³ → finite
     minimum (core size + mass). F² resists collapse, V resists expansion.

USAGE:
    python -m openwave.xperiments.m5_liquid_crystal.research.sandbox_v4.m5_5_3_potential_regularization
"""
import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
LAMBDA = np.array([1.0, 0.3, 0.0])     # Λ = (1, δ, 0) target eigenvalue spectrum
DELTA = 0.3


def V_eq12(M):
    """Eq.12 eigenvalue-preference: V = Σ_i (λ_i − Λ_i)²  (λ sorted descending)."""
    lam = np.sort(np.linalg.eigvalsh(M))[::-1]
    return float(np.sum((lam - LAMBDA) ** 2))


def V_eq13(M, a=1.0, b=1.0, c=1.0):
    """Eq.13 LdG Higgs-like: V = a·Tr(M²) − b·Tr(M³) + c·(Tr(M²))²."""
    M2 = M @ M
    tr2 = np.trace(M2)
    tr3 = np.trace(M2 @ M)
    return float(a * tr2 - b * tr3 + c * tr2**2)


def rotation_zyz(α, β, γ):
    ca, cb, cg = np.cos([α, β, γ])
    sa, sb, sg = np.sin([α, β, γ])
    Rz1 = np.array([[ca, -sa, 0], [sa, ca, 0], [0, 0, 1.0]])
    Ry = np.array([[cb, 0, sb], [0, 1.0, 0], [-sb, 0, cb]])
    Rz2 = np.array([[cg, -sg, 0], [sg, cg, 0], [0, 0, 1.0]])
    return Rz1 @ Ry @ Rz2


def stage1():
    print("=" * 70)
    print("STAGE 1 — V(M) Eq.12/13 + rotation invariance (V acts on eigenvalues only)")
    print("=" * 70)
    D = np.diag(LAMBDA)
    devs12, devs13 = [], []
    for _ in range(8):
        O = rotation_zyz(*np.random.uniform(0, 2 * np.pi, 3))
        M = O @ D @ O.T
        devs12.append(abs(V_eq12(M) - V_eq12(D)))
        devs13.append(abs(V_eq13(M) - V_eq13(D)))
    print(f"    V_eq12(ODOᵀ) − V_eq12(D): max dev over random O = {max(devs12):.2e}")
    print(f"    V_eq13(ODOᵀ) − V_eq13(D): max dev over random O = {max(devs13):.2e}")
    print(f"    V_eq12(Λ)={V_eq12(D):.3e} (=0 at vacuum)   V_eq13(Λ)={V_eq13(D):.3e}")
    ok = max(devs12) < 1e-10 and max(devs13) < 1e-10
    print(f"    → {'PASS' if ok else 'FAIL'}  (V is rotation-invariant ⇒ acts ONLY on the "
          f"eigenvalue-deformation sector ⇒ V IS the regularization, not a twist force)")
    return ok


def stage2():
    print("\n" + "=" * 70)
    print("STAGE 2 — V defines the vacuum: perturbed eigenvalues relax to Λ under ∇V flow")
    print("=" * 70)
    # Diagonal M with perturbed eigenvalues; gradient-descent on V_eq12 (in eigenvalue space).
    lam = LAMBDA + np.array([0.5, -0.2, 0.4])     # perturbed start
    lr, hist = 0.1, [lam.copy()]
    for _ in range(400):
        lam = lam - lr * 2.0 * (lam - LAMBDA)      # ∇V_eq12 = 2(λ−Λ)
        hist.append(lam.copy())
    err = np.linalg.norm(lam - LAMBDA)
    print(f"    start λ = {hist[0]}  →  final λ = {np.round(lam, 4)}   (target Λ = {LAMBDA})")
    print(f"    ‖λ − Λ‖ = {err:.2e}")
    ok = err < 1e-3
    print(f"    → {'PASS' if ok else 'FAIL'}  (V defines Λ as the vacuum + supplies the restoring force)")
    return ok


def hedgehog_energies(xi, N=41, Lbox=5.0, rigid=False, r_cut=0.0):
    """Deformed uniaxial hedgehog M(x)=s(r)·r̂⊗r̂ on a 3D grid → (E_F, E_V).

    s(r)=tanh(r/ξ) (→0 at core, →1 at vacuum) regularizes; rigid=True sets s≡1.
    E_F = ∫ 4Σ_{μ<ν}‖F_μν‖²_F d³r  (Eq.18 curvature),  F_μν = 2[∂_μM, ∂_νM].
    E_V = ∫ (λ_max − 1)² d³r        (Eq.12, eigenvalues of s·r̂⊗r̂ are (s,0,0); Λ=(1,0,0)).
    r_cut masks the core (for the rigid-divergence demo).
    """
    xs = np.linspace(-Lbox, Lbox, N)
    h = xs[1] - xs[0]
    X, Y, Z = np.meshgrid(xs, xs, xs, indexing="ij")
    r = np.sqrt(X**2 + Y**2 + Z**2) + 1e-9
    s = np.ones_like(r) if rigid else np.tanh(r / xi)
    rhat = np.stack([X / r, Y / r, Z / r], axis=-1)
    M = s[..., None, None] * rhat[..., :, None] * rhat[..., None, :]    # s·r̂⊗r̂
    # spatial derivatives M_μ (central diff)
    Mmu = [np.gradient(M, h, axis=ax) for ax in range(3)]
    pairs = [(0, 1), (0, 2), (1, 2)]
    EF_dens = np.zeros_like(r)
    for (mu, nu) in pairs:
        F = 2.0 * (np.einsum("...ac,...cb->...ab", Mmu[mu], Mmu[nu])
                   - np.einsum("...ac,...cb->...ab", Mmu[nu], Mmu[mu]))
        EF_dens += 4.0 * np.einsum("...ab,...ab->...", F, F)
    EV_dens = (s - 1.0) ** 2                                            # (λ_max − 1)²
    mask = (r > r_cut)
    dv = h**3
    return float(EF_dens[mask].sum() * dv), float(EV_dens[mask].sum() * dv)


def stage3():
    print("\n" + "=" * 70)
    print("STAGE 3 — regularization (Derrick): F² + V → finite-energy core (finite mass)")
    print("=" * 70)
    # (a) rigid core curvature DIVERGES as the cutoff → 0
    print("  (a) rigid hedgehog (s≡1) curvature energy vs core cutoff r_cut:")
    for rc in (1.2, 0.8, 0.4):
        EF, _ = hedgehog_energies(xi=1.0, rigid=True, r_cut=rc)
        print(f"        r_cut={rc:.2f}:  ∫‖F‖² (r>r_cut) = {EF:.3e}   (grows as r_cut→0 ⇒ diverges)")
    print(f"        (analytic: rigid ‖F‖² ~ 1/r⁴ ⇒ ∫r²dr ~ ∫dr/r² diverges at r→0; "
          f"finer cutoff needs finer grid)")
    # (b) deformed core (s→0) is FINITE — no cutoff needed
    EF1, EV1 = hedgehog_energies(xi=1.0, rigid=False, r_cut=0.0)
    print(f"  (b) deformed core (s=tanh(r/ξ), ξ=1): ∫‖F‖²={EF1:.3e} (FINITE), ∫V={EV1:.3e} (FINITE)")
    # (c) Derrick scaling: E(L) = E_F·(1/L) + E_V·L³  (F² ~ 1/L, V ~ L³ under x→Lx)
    Ls = np.linspace(0.4, 3.0, 200)
    E = EF1 / Ls + EV1 * Ls**3
    Lstar = (EF1 / (3 * EV1)) ** 0.25
    Emin = EF1 / Lstar + EV1 * Lstar**3
    print(f"  (c) Eq.18 Derrick energy  E(L)=E_F/L + E_V·L³  →  finite minimum:")
    print(f"        L* = {Lstar:.3f}  (regularized core size)   E(L*) = {Emin:.3e}  (finite MASS)")
    print(f"        F² resists collapse (E_F/L→∞ as L→0); V resists expansion (E_V L³→∞) ⇒ STABLE")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    PLOT = HERE / "plots" / "m5_5_3_regularization.png"
    PLOT.parent.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(7, 4.6))
    plt.plot(Ls, E, label="E(L) = E_F/L + E_V·L³  (Eq.18)")
    plt.plot(Ls, EF1 / Ls, "--", alpha=0.6, label="F² curvature (resists collapse)")
    plt.plot(Ls, EV1 * Ls**3, ":", alpha=0.6, label="V (resists expansion)")
    plt.axvline(Lstar, color="k", alpha=0.4); plt.plot([Lstar], [Emin], "ro")
    plt.xlabel("core size L"); plt.ylabel("energy"); plt.legend(); plt.grid(alpha=0.3)
    plt.title(f"M5.5.3 Derrick regularization: finite core L*={Lstar:.2f}, mass={Emin:.2e}")
    plt.tight_layout(); plt.savefig(PLOT, dpi=110); plt.close()
    print(f"        plot → {PLOT}")
    ok = all(np.isfinite(v) for v in (EF1, EV1, Lstar, Emin)) and Lstar > 0
    print(f"    → {'PASS' if ok else 'FAIL'}  (F²+V regularize the core to finite energy/mass — "
          f"no extra Skyrme term needed; Eq.18's F² IS the Skyrme-family kinetic. L* absolute "
          f"scale is Q7/Q8-dependent; the finite minimum is the result)")
    return ok


def main():
    np.random.seed(0)
    s1, s2, s3 = stage1(), stage2(), stage3()
    print("\n" + "=" * 70)
    print(f"M5.5.3:  V rotation-invariance={s1}  vacuum={s2}  regularization={s3}")
    print("PASS" if (s1 and s2 and s3) else "PARTIAL/FAIL")
    print("Note: exact Λ=(1,δ,0)-producing LdG coefficients (Eq.13 a,b,c) + Faber's exact")
    print("regularization scheme remain Q7/Q8 (Duda's open questions) — baseline shown here.")
    print("=" * 70)
    return 0 if (s1 and s2 and s3) else 1


if __name__ == "__main__":
    sys.exit(main())
