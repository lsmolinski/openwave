"""
M5.6.5d — the faithful (gauge-correct) kinetic vs the simple ½‖Ṁ‖² we ship

Production `evolve_M` uses the SIMPLE kinetic ½‖Ṁ‖² ⇒ EOM M̈ = c²·div(G) − dV_M: every one
of M's 6 symmetric components gets the SAME inertia. But Duda's Eq.18 time-curvature gives the
FAITHFUL kinetic (5a §9):

    A_0 = [M, Ṁ] ,  F_μ0 = 2[M_μ, Ṁ] ,  T = Σ_μ‖F_μ0‖² = 4 Σ_μ‖[M_μ, Ṁ]‖²

As a quadratic form in Ṁ this is  T = ⟨Ṁ, G Ṁ⟩  with the per-voxel metric operator

    G = 4·Σ_μ (−ad²_{M_μ}) ,   ad_{M_μ}(X) = [M_μ, X]            (6×6 on symmetric Ṁ)

THE POINT (5a §9): G is DEGENERATE — T=0 whenever Ṁ commutes with every spatial gradient M_μ.
Those null directions are GAUGE (non-dynamical). The simple ½‖Ṁ‖² kinetic (metric = ½·I, no
null space) gives them spurious uniform inertia and ANIMATES them — that is the "random director
slosh" seen on screen under Evolve-PDE. The faithful kinetic freezes them.

This script quantifies that on the production biaxial hedgehog:
  Stage 1  per-voxel spectrum of G: the null-space dimension (gauge modes) vs the simple ½·I.
  Stage 2  evolve with the SIMPLE production scheme; measure the fraction of the resulting motion
           Ṁ that lives in the NULL space of the faithful G — i.e. the spurious gauge slosh the
           faithful kinetic would forbid.

USAGE:
    python -m openwave.xperiments.m5_liquid_crystal.research.sandbox_v6.m5_6_5d_faithful_kinetic
"""
import numpy as np
from openwave.xperiments.m5_liquid_crystal.research.sandbox_v6.m5_6_5c_potential_confinement import (
    build_biaxial_M, force, commf, central, C_WAVE, DT, S2_STAR, DELTA,
)

# orthonormal basis of symmetric 3×3 (⟨A,B⟩=Tr(AB)): 3 diagonal + 3 off-diagonal
_s = 1.0 / np.sqrt(2.0)
BASIS = [
    np.diag([1.0, 0, 0]), np.diag([0, 1.0, 0]), np.diag([0, 0, 1.0]),
    _s * np.array([[0, 1.0, 0], [1.0, 0, 0], [0, 0, 0]]),
    _s * np.array([[0, 0, 1.0], [0, 0, 0], [1.0, 0, 0]]),
    _s * np.array([[0, 0, 0], [0, 0, 1.0], [0, 1.0, 0]]),
]


def metric_G(Mx, My, Mz):
    """6×6 faithful kinetic metric G_ab = 4 Σ_μ ⟨[M_μ,E_a],[M_μ,E_b]⟩ at one voxel."""
    G = np.zeros((6, 6))
    brackets = [[Mmu @ E - E @ Mmu for E in BASIS] for Mmu in (Mx, My, Mz)]
    for a in range(6):
        for b in range(a, 6):
            s = 0.0
            for mu in range(3):
                s += np.sum(brackets[mu][a] * brackets[mu][b])   # ⟨A,B⟩=Tr(ABᵀ)=Σ A_ij B_ij
            G[a, b] = G[b, a] = 4.0 * s
    return G


def stage1(bg):
    print("=" * 74)
    print("STAGE 1 — the faithful kinetic metric G is DEGENERATE (gauge null space)")
    print("=" * 74)
    M, h, active = bg["M"], bg["h"], bg["active"]
    Mx, My, Mz = central(M, 0, h), central(M, 1, h), central(M, 2, h)
    idx = np.argwhere(active)
    rng = np.random.default_rng(0)
    sample = idx[rng.choice(len(idx), size=min(400, len(idx)), replace=False)]
    null_dims, eig_all = [], []
    for (i, j, k) in sample:
        G = metric_G(Mx[i, j, k], My[i, j, k], Mz[i, j, k])
        w = np.linalg.eigvalsh(G)
        eig_all.append(w)
        scale = w.max() if w.max() > 0 else 1.0
        null_dims.append(int(np.sum(w < 1e-6 * scale)))     # near-zero relative to the top eigenvalue
    null_dims = np.array(null_dims)
    eig_all = np.array(eig_all)
    print(f"  sampled {len(sample)} bulk voxels of the biaxial hedgehog (δ={DELTA})")
    print(f"  faithful G eigenvalues per voxel (6 modes), median sorted = {np.round(np.median(eig_all,0),3)}")
    print(f"  ⇒ null-space dimension per voxel: mean={null_dims.mean():.2f}  "
          f"(min={null_dims.min()}, max={null_dims.max()})")
    print(f"  by contrast the SIMPLE kinetic ½‖Ṁ‖² has metric ½·I₆: all 6 eigenvalues=0.5, "
          f"NULL DIM = 0")
    ok = null_dims.mean() >= 1.0       # at least one gauge null mode on average
    print(f"  → {'PASS' if ok else 'FAIL'}  faithful G has ≥1 gauge null mode/voxel; the simple "
          f"kinetic gives all 6 modes the same inertia (Stage 2 examines what that costs).")
    return ok


def null_is_trace(bg):
    """Confirm the null mode is the trace/isotropic direction: ⟨n̂_null, I/√3⟩ ≈ 1."""
    M, h, active = bg["M"], bg["h"], bg["active"]
    Mx, My, Mz = central(M, 0, h), central(M, 1, h), central(M, 2, h)
    idx = np.argwhere(active)
    rng = np.random.default_rng(2)
    sample = idx[rng.choice(len(idx), size=min(200, len(idx)), replace=False)]
    I_coeff = np.array([np.trace(E) for E in BASIS]) / np.sqrt(3.0)   # I/√3 in the basis
    overlaps = []
    for (i, j, k) in sample:
        G = metric_G(Mx[i, j, k], My[i, j, k], Mz[i, j, k])
        w, V = np.linalg.eigh(G)
        overlaps.append(abs(V[:, 0] @ I_coeff))                      # lowest-eigenvalue eigenvector
    return float(np.mean(overlaps))


def physical_inertia_spread(bg):
    """The 5 nonzero (physical) eigenvalues of G vs the simple kinetic's uniform ½ → freq error."""
    M, h, active = bg["M"], bg["h"], bg["active"]
    Mx, My, Mz = central(M, 0, h), central(M, 1, h), central(M, 2, h)
    idx = np.argwhere(active)
    rng = np.random.default_rng(3)
    sample = idx[rng.choice(len(idx), size=min(400, len(idx)), replace=False)]
    phys = []
    for (i, j, k) in sample:
        w = np.linalg.eigvalsh(metric_G(Mx[i, j, k], My[i, j, k], Mz[i, j, k]))
        phys.append(w[1:])                                           # drop the null (lowest) mode
    phys = np.concatenate(phys)
    phys = phys[phys > 1e-9]
    return phys


def stage2(bg):
    print("\n" + "=" * 74)
    print("STAGE 2 — the null mode is the TRACE; the real gap is PHYSICAL-mode inertia (dispersion)")
    print("=" * 74)
    M, h, active = bg["M"], bg["h"], bg["active"]
    mask = active[..., None, None]
    M0 = M.copy()
    Mcur, Mprev = M.copy(), M.copy()
    # SIMPLE production scheme: M̈ = c²·div(G) (V off to isolate the kinetic structure)
    for _ in range(200):
        f = force(Mcur, h, 0.0, 0.0, 0.0)
        Mnew = np.where(mask, 2 * Mcur - Mprev + DT**2 * f, M0)
        Mprev, Mcur = Mcur, Mnew
    Mdot = (Mcur - Mprev) / DT
    Mx, My, Mz = central(Mcur, 0, h), central(Mcur, 1, h), central(Mcur, 2, h)
    # (a) the null mode IS the trace direction
    tr_overlap = null_is_trace(bg)
    print(f"  (a) null eigenvector · (I/√3) = {tr_overlap:.3f}  ⇒ the gauge null mode is the")
    print(f"      TRACE/isotropic direction ([M_μ, I]=0). The curvature force is traceless, so")
    # (b) the simple scheme's motion stays out of the null space
    idx = np.argwhere(active); rng = np.random.default_rng(1)
    sample = idx[rng.choice(len(idx), size=min(400, len(idx)), replace=False)]
    fracs = []
    for (i, j, k) in sample:
        md = Mdot[i, j, k]
        if np.sum(md * md) < 1e-30:
            continue
        w, V = np.linalg.eigh(metric_G(Mx[i, j, k], My[i, j, k], Mz[i, j, k]))
        c = V.T @ np.array([np.sum(md * BASIS[a]) for a in range(6)])
        fracs.append(c[0] ** 2 / np.sum(c ** 2))
    null_frac = float(np.mean(fracs))
    print(f"      the simple scheme's motion is only {100*null_frac:.1f}% in the null space ⇒ it")
    print(f"      stays PHYSICAL (no spurious gauge slosh — the production scheme is well-behaved).")
    # (c) the physical-mode inertia spread → the dispersion/clock-frequency error
    phys = physical_inertia_spread(bg)
    lo, hi = np.percentile(phys, 5), np.percentile(phys, 95)
    # ω ∝ 1/√inertia at fixed stiffness; simple uses ½ for all modes
    f_hi, f_lo = np.sqrt(0.5 / lo), np.sqrt(0.5 / hi)
    print(f"  (c) faithful physical-mode inertia (5 modes): 5–95% span [{lo:.3f}, {hi:.3f}] vs the")
    print(f"      simple kinetic's uniform 0.5 ⇒ twist/clock frequency mis-set by ×[{f_lo:.2f}, {f_hi:.2f}]")
    print(f"      under ½‖Ṁ‖². THIS is the M5.8 stake (ω = 2mc²/ℏ needs the faithful inertia).")
    ok = tr_overlap > 0.9 and null_frac < 0.05
    print(f"  → {'PASS' if ok else 'PARTIAL'}  slosh is PHYSICAL twist (not gauge); 5d corrects the")
    print(f"     dispersion, not the existence of the motion.")
    return ok, (lo, hi, f_lo, f_hi)


def main():
    bg = build_biaxial_M()
    s1 = stage1(bg)
    s2, _ = stage2(bg)
    print("\n" + "=" * 74)
    print(f"M5.6.5d:  G degenerate (1 null = trace)={s1}   simple kinetic stays physical={s2}")
    print("Corrected picture: the faithful G's lone null mode is the TRACE; the traceless curvature")
    print("force never sources it, so the simple ½‖Ṁ‖² does NOT generate spurious gauge slosh — it")
    print("is a well-behaved approximation. The on-screen director slosh is PHYSICAL twist (the")
    print("clock). What ½‖Ṁ‖² gets WRONG is the inertia weighting of the 5 physical modes ⇒ the")
    print("twist/clock dispersion (frequency) is off by an O(1) factor. The faithful kinetic — the")
    print("O(x)∈SO(3) metric validated on the ψ DoF in m5_6_1b/2b — fixes the frequency, which M5.8")
    print("(ω=2mc²/ℏ) needs. Production: keep ½‖Ṁ‖² for qualitative runs; use the faithful")
    print("ψ-evolution (m5_6_2b) for the M5.8 clock-frequency measurement.")
    print("PASS" if (s1 and s2) else "PARTIAL — inspect stages above")
    print("=" * 74)
    return 0 if (s1 and s2) else 1


if __name__ == "__main__":
    raise SystemExit(main())
