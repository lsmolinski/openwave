"""
M5.6.5c — turning V on: amplitude confinement vs the biaxiality it erodes

The user's M5.5.4 observation: with V off the Evolve-PDE hedgehog *sloshes* and its
energy DILUTES over a growing radius (bounded, energy-conserving, but not localized —
nothing supplies a restoring force against amplitude spread). This step turns the
production V_M (Eq.13 LdG, already in engine2_pde.py, off by default) ON and shows it
confines the energy — plus an important structural finding about WHICH V can do it.

THE FINDING (Stage 1).  The production potential is

    V(M) = a·Tr(M²) − b·Tr(M³) + c·(Tr(M²))²            (Eq.13, rotation-invariant)
    ∂V/∂λ_i = λ_i·(2a − 3b·λ_i + 4c·s₂)  ,  s₂ = Tr(M²) = Σλ²

At a critical point each λ_i is either 0 OR the SINGLE root λ* = (2a+4c·s₂)/(3b) (one
linear equation, shared s₂). So every nonzero eigenvalue equals λ*: the anisotropic
critical points are uniaxial — (λ*,λ*,0), (λ*,0,0). **The canonical 3-term LdG cannot
have a biaxial (1,δ,0) minimum with three DISTINCT eigenvalues.** Consequence: turning
on the full Eq.13 (b≠0) confines the amplitude but PULLS δ toward a uniaxial value —
it erodes the very biaxiality M5.6.2 needs (the C_μν≠0 mass source). A genuinely
biaxial-stable vacuum needs an extra invariant in V (a biaxiality term) — a Q7 refinement
to flag to Duda.

THE CLEAN CONFINEMENT (Stage 2).  Set b=0:  V = a·s₂ + c·s₂²  depends ONLY on s₂.
Minimum at s₂* = −a/(2c). Choose s₂* = Tr(diag(1,δ,0)²) = 1+δ². This pins the amplitude
(confines) and is EXACTLY FLAT in the biaxiality direction (V is constant on the s₂
sphere) — confine without uniaxializing. This is the coefficient set we wire to production.

THE DEMONSTRATION (Stage 3).  Evolve the full biaxial hedgehog M with the production
leapfrog  M_tt = c²·div(G) − dV_M(M),  G_α = 8Σ_ν[[M_α,M_ν],M_ν]  (mirrors evolve_M +
compute_curvature_flux). V-off vs V-on (b=0 well). Metric = energy-weighted RMS radius
R_rms(t): V-off grows (dilution, the user's observation), V-on stays bounded (contained).

USAGE:
    python -m openwave.xperiments.m5_liquid_crystal.research.sandbox_v5.m5_6_5c_potential_confinement
"""
import numpy as np

DELTA = 0.3
S2_STAR = 1.0 + DELTA**2                # Tr(diag(1,δ,0)²) — the target amplitude
N = 40
L = 5.0
RC = 0.9                                # core melt (point defect)
RHOC = 0.9                              # disclination melt (z-axis)
C_WAVE = 0.3                            # matches the GUI c (5a §… bounded-not-bug run)
C_COEF = 2.0                            # well stiffness (a = −2·c·s₂* keeps the min at s₂*)
DT = 0.006
N_STEPS = 900


# ----- matrix helpers (numpy mirror of the @ti.func ops) ----------------------
def commf(A, B):
    return np.einsum("...ac,...cb->...ab", A, B) - np.einsum("...ac,...cb->...ab", B, A)


def frob2(A):
    return np.einsum("...ab,...ab->...", A, A)


def tr(A):
    return np.einsum("...aa->...", A)


def central(f, axis, h):
    out = np.zeros_like(f)
    sp, sm, so = [slice(None)] * f.ndim, [slice(None)] * f.ndim, [slice(None)] * f.ndim
    sp[axis], sm[axis], so[axis] = slice(2, None), slice(0, -2), slice(1, -1)
    out[tuple(so)] = (f[tuple(sp)] - f[tuple(sm)]) / (2 * h)
    return out


def V_M(M, a, b, c):
    M2 = np.einsum("...ac,...cb->...ab", M, M)
    tr2 = tr(M2)
    tr3 = tr(np.einsum("...ac,...cb->...ab", M2, M))
    return a * tr2 - b * tr3 + c * tr2 * tr2


def dV_M(M, a, b, c):
    M2 = np.einsum("...ac,...cb->...ab", M, M)
    tr2 = tr(M2)[..., None, None]
    return 2.0 * a * M - 3.0 * b * M2 + 4.0 * c * tr2 * M


# ============================================================================
def stage1():
    print("=" * 74)
    print("STAGE 1 — Eq.13 critical points are UNIAXIAL (no biaxial minimum)")
    print("=" * 74)
    print("  ∂V/∂λ_i = λ_i·(2a − 3b·λ_i + 4c·s₂) = 0 ⇒ λ_i ∈ {0, λ*} (one shared root λ*)")
    print("  ⇒ every nonzero eigenvalue equals λ* ⇒ anisotropic minima are uniaxial.\n")
    rng = np.random.default_rng(0)
    max_spread = 0.0
    n_biaxial = 0
    for (a, b, c) in [(1.0, 1.0, 1.0), (-2.18, 1.0, 1.0), (1.0, 2.0, 0.5), (-1.0, 0.5, 1.0)]:
        endpoints = []
        for _ in range(120):
            lam = rng.uniform(-0.5, 1.6, 3)
            for _ in range(3000):
                s2 = np.sum(lam**2)
                g = lam * (2 * a - 3 * b * lam + 4 * c * s2)
                lam = np.clip(lam - 0.01 * g, -3, 3)
            endpoints.append(np.sort(lam)[::-1])
        endpoints = np.array(endpoints)
        # spread of the NONZERO eigenvalues at each endpoint (biaxial ⇒ two distinct nonzeros)
        spreads = []
        for lam in endpoints:
            nz = lam[np.abs(lam) > 1e-3]
            spreads.append(np.ptp(nz) if len(nz) >= 2 else 0.0)
        spreads = np.array(spreads)
        # a biaxial minimum would have two DISTINCT nonzero eigenvalues differing by ≫ tol
        biax = int(np.sum(spreads > 1e-2))
        n_biaxial += biax
        max_spread = max(max_spread, spreads.max())
        lam_ex = endpoints[np.argmin([V_M(np.diag(l), a, b, c) for l in endpoints])]
        print(f"  (a,b,c)=({a:+.2f},{b:.2f},{c:.2f}): min eig-set ≈ {np.round(lam_ex, 3)}  "
              f"max nonzero-spread={spreads.max():.1e}  biaxial-minima={biax}/120")
    ok = (n_biaxial == 0) and (max_spread < 1e-2)
    print(f"\n  → {'PASS' if ok else 'FAIL'}  no biaxial minimum found for any (a,b,c): the "
          f"nonzero eigenvalues always collapse to a single λ*.")
    print("    ⇒ full Eq.13 (b≠0) confines amplitude but erodes δ → uniaxial. Biaxiality-")
    print("      stable V needs an extra invariant (Q7 flag). Use b=0 amplitude well below.")
    return ok


def stage2():
    print("\n" + "=" * 74)
    print("STAGE 2 — b=0 amplitude well: confine Tr(M²)→1+δ², leave δ exactly flat")
    print("=" * 74)
    c_coef = C_COEF
    a_coef = -2.0 * c_coef * S2_STAR          # min of a·x+c·x² at x=−a/2c = S2_STAR
    s2 = np.linspace(0.0, 2.5, 400)
    V = a_coef * s2 + c_coef * s2**2
    s2_min = s2[np.argmin(V)]
    curv = 2 * c_coef                          # d²V/d(s₂)² = 2c  (restoring stiffness)
    print(f"  V(s₂) = {a_coef:.3f}·s₂ + {c_coef:.3f}·s₂² , target s₂* = 1+δ² = {S2_STAR:.3f}")
    print(f"    numeric min at s₂ = {s2_min:.3f}  (curvature d²V/ds₂² = {curv:.2f} > 0 ⇒ restoring)")
    # biaxiality-flatness check: V identical for any eigenvalue set with the same s₂
    M_biax = np.diag([1.0, DELTA, 0.0])
    M_uni = np.diag([np.sqrt(S2_STAR / 2), np.sqrt(S2_STAR / 2), 0.0])   # same s₂, uniaxial
    vb, vu = V_M(M_biax, a_coef, 0.0, c_coef), V_M(M_uni, a_coef, 0.0, c_coef)
    print(f"    V(biaxial 1,δ,0) = {vb:.4f}  vs  V(uniaxial same s₂) = {vu:.4f}  "
          f"|Δ|={abs(vb - vu):.1e} ⇒ δ flat ✓")
    ok = abs(s2_min - S2_STAR) < 0.05 and abs(vb - vu) < 1e-9
    print(f"  → {'PASS' if ok else 'FAIL'}  (a,b,c)=({a_coef:.3f}, 0, {c_coef:.3f}) "
          f"confines amplitude without touching biaxiality.")
    return ok, (a_coef, 0.0, c_coef)


def build_biaxial_M():
    """Production-equivalent biaxial hedgehog M = O·D(s(r))·Oᵀ with radial + disclination melt."""
    xs = np.linspace(-L, L, N)
    h = xs[1] - xs[0]
    X, Y, Z = np.meshgrid(xs, xs, xs, indexing="ij")
    r = np.sqrt(X**2 + Y**2 + Z**2) + 1e-9
    rho = np.sqrt(X**2 + Y**2) + 1e-9
    rhat = np.stack([X / r, Y / r, Z / r], -1)
    azim = np.stack([-Y / rho, X / rho, np.zeros_like(r)], -1)
    sdisc = np.clip(rho / RHOC, 0, 1)
    shrink = (sdisc**2 * (3 - 2 * sdisc))[..., None]
    ephi = azim * shrink
    ephi = ephi - np.einsum("...a,...a->...", ephi, rhat)[..., None] * rhat
    etheta = np.cross(ephi, rhat)
    srad = (r / np.sqrt(r**2 + RC**2))[..., None, None]
    diso = (1.0 + DELTA) / 3.0
    d0 = diso + srad[..., 0, 0] * (1.0 - diso)
    d1 = diso + srad[..., 0, 0] * (DELTA - diso)
    d2 = diso + srad[..., 0, 0] * (0.0 - diso)
    oprod = lambda v: v[..., :, None] * v[..., None, :]
    M = (d0[..., None, None] * oprod(rhat) + d1[..., None, None] * oprod(etheta)
         + d2[..., None, None] * oprod(ephi))
    interior = np.zeros(r.shape, bool)
    interior[2:-2, 2:-2, 2:-2] = True
    active = (r > 2 * RC) & (rho > RHOC) & interior     # off point-core + disclination
    return dict(M=M, h=h, r=r, active=active)


def force(M, h, a, b, c):
    Mx, My, Mz = central(M, 0, h), central(M, 1, h), central(M, 2, h)
    cc = lambda A, B: commf(commf(A, B), B)
    Gx = 8.0 * (cc(Mx, My) + cc(Mx, Mz))
    Gy = 8.0 * (cc(My, Mx) + cc(My, Mz))
    Gz = 8.0 * (cc(Mz, Mx) + cc(Mz, My))
    divG = central(Gx, 0, h) + central(Gy, 1, h) + central(Gz, 2, h)
    return C_WAVE**2 * divG - dV_M(M, a, b, c)


def energy_density(M, M_prev, h, a, b, c):
    Mdot = (M - M_prev) / DT
    kinetic = 0.5 * frob2(Mdot)
    Mx, My, Mz = central(M, 0, h), central(M, 1, h), central(M, 2, h)
    curv = 4.0 * (frob2(commf(Mx, My)) + frob2(commf(Mx, Mz)) + frob2(commf(My, Mz)))
    return kinetic + C_WAVE**2 * curv + V_M(M, a, b, c)


def evolve(M0, bg, a, b, c, n_steps):
    """Returns (amp_dev[t], R_rms[t], finite). amp_dev = ⟨|Tr(M²)−s₂*|⟩ over the active
    region — the DIRECT confinement metric (V pins the amplitude; V-off lets it wander)."""
    h, r, active = bg["h"], bg["r"], bg["active"]
    M = M0.copy()
    M_prev = M0.copy()
    mask = active[..., None, None]
    amp_dev, R_rms = [], []
    for step in range(n_steps):
        f = force(M, h, a, b, c)
        M_new = np.where(mask, 2 * M - M_prev + DT**2 * f, M0)   # Dirichlet: boundary fixed
        M_prev, M = M, M_new
        if step % 30 == 0:
            s2 = tr(np.einsum("...ac,...cb->...ab", M, M))
            amp_dev.append(float(np.mean(np.abs(s2 - S2_STAR)[active])))
            H = energy_density(M, M_prev, h, a, b, c)
            Hf = np.clip(H - H.min(), 0, None) * active           # localize-able energy
            R_rms.append(np.sqrt((Hf * r**2).sum() / (Hf.sum() + 1e-30)))
    return np.array(amp_dev), np.array(R_rms), np.isfinite(M).all()


def stage3(coeffs):
    a, b, c = coeffs
    print("\n" + "=" * 74)
    print("STAGE 3 — dynamical containment: full-M leapfrog, V-off vs V-on (b=0 well)")
    print("=" * 74)
    bg = build_biaxial_M()
    print(f"  grid {N}³  L={L}  c={C_WAVE}  dt={DT}  steps={N_STEPS}  active={100*bg['active'].mean():.0f}%")
    Aoff, Roff, fin_off = evolve(bg["M"], bg, 0.0, 0.0, 0.0, N_STEPS)   # V off
    Aon, Ron, fin_on = evolve(bg["M"], bg, a, b, c, N_STEPS)            # V on (amplitude well)
    print(f"  amplitude deviation ⟨|Tr(M²)−s₂*|⟩  (the confinement metric, start → max → end):")
    print(f"    V OFF : {Aoff[0]:.4f} → {Aoff.max():.4f} → {Aoff[-1]:.4f}   finite={fin_off}")
    print(f"    V ON  : {Aon[0]:.4f} → {Aon.max():.4f} → {Aon[-1]:.4f}   finite={fin_on}")
    print(f"  energy-weighted RMS radius R_rms (start → end):")
    print(f"    V OFF : {Roff[0]:.3f} → {Roff[-1]:.3f}    V ON : {Ron[0]:.3f} → {Ron[-1]:.3f}")
    # V-on pins the amplitude: its max deviation stays well below V-off's
    contained = fin_on and (Aon.max() < 0.6 * Aoff.max()) and (Aon[-1] < Aoff[-1])
    print(f"  → {'PASS' if contained else 'PARTIAL'}  V-on pins Tr(M²) near s₂* "
          f"(max dev {Aon.max():.3f} vs {Aoff.max():.3f}); V-off amplitude wanders = the dilution.")
    return contained


def main():
    s1 = stage1()
    s2, coeffs = stage2()
    s3 = stage3(coeffs)
    print("\n" + "=" * 74)
    print(f"M5.6.5c:  Eq.13-uniaxial-only finding={s1}   b=0 amplitude well={s2}   containment={s3}")
    print("Production coeffs (b=0 amplitude well, confines without uniaxializing):")
    print(f"    ldg_a = {coeffs[0]:.3f}   ldg_b = 0.0   ldg_c = {coeffs[2]:.3f}   (s₂*=1+δ²={S2_STAR:.3f})")
    print("Q7 flag for Duda: a fully biaxial-STABLE vacuum needs an extra invariant in V")
    print("(the 3-term Eq.13 minima are uniaxial). b=0 confines + leaves δ flat as the interim.")
    print("PASS" if (s1 and s2 and s3) else "PARTIAL — inspect stages above")
    print("=" * 74)
    return 0 if (s1 and s2 and s3) else 1


if __name__ == "__main__":
    raise SystemExit(main())
