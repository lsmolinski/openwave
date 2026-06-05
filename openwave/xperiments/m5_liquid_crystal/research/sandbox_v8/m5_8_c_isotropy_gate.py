"""
M5.8 carry-over — SR-emergence gate: is the propagation speed isotropic? (c-isotropy)

PROVENANCE (2026-06-05): an external "fabric of spacetime" comparison table (LC vs elastic
continuum, Close-school) critiques liquid-crystal media with "no fixed c unless crosslinked;
anisotropic otherwise". The one-constant Frank approximation answers this for the QUADRATIC
gradient energy — but production runs Duda's Eq.18 F² action, whose curvature force
G_α = 8 Σ_ν [[M_α,M_ν],M_ν] is CUBIC in gradients (engine2_pde.compute_curvature_flux), so
isotropy of the emergent wave cone is NOT guaranteed by construction. Nobody had measured it.

Three parts, each answerable without the full production engine (V(M) is rotation-invariant
⇒ doesn't touch the twist sector; the M5.8.1 frozen g-axis decouples ⇒ 3×3 block suffices):

  A  UNIFORM VACUUM NULL: around a uniform biaxial vacuum M₀ (∂M₀=0) the F² force is cubic
     in the perturbation gradient — there is NO linear-order propagation at all. Sharper: a
     single-matrix-direction perturbation δM = ε·f(x)·E has [f_αE, f_νE] ≡ 0 ⇒ EXACTLY zero
     force (machine-exact, all orders); a two-direction perturbation scales as ε³ (log-log
     slope 3). ⇒ "what is c in the vacuum?" is the WRONG question for Eq.18 — the uniform
     vacuum is transparent; propagation is background-gradient-mediated (defect dressings).
     Hooks: 9a (what carries EM packets in empty space?) + Q9 (deeper substrate beneath M).

  B  STATIC CONE TENSOR (far-field hedgehog): the verified twist EOM (m5_6_2b, 5a §9)
     2K ψ_tt = Σ_μ ∂_μ J_μ has principal symbol  c²(k̂) = 4·(1 − k̂·Gram·k̂ / T),
     Gram_μν = P_μ•P_ν (Frobenius), T = Tr Gram, K = 4T. Isotropy ⟺ Gram ∝ I. For the
     biaxial hedgehog the far-field frame depends on ANGLES only ⇒ P_radial ≈ 0 ⇒ Gram has
     a soft radial eigenvalue ⇒ predicted anisotropic cone: c(r̂)/c(tangential) ≈ √2, with
     the soft Gram eigenvector aligned to r̂. Measured per voxel on a far-field shell, plus
     r-binning (scale-free cone?) — pure background algebra, no evolution.

  C  DYNAMIC CROSS-CHECK: with use_C=False the m5_6_2b EOM is exactly the homogeneous
     linear operator (force is affine in ψ; dropping C removes the source term only), so a
     seeded blob propagates with NO background C-drive contamination. Seed a Gaussian twist
     blob at (R,0,0); along 5 rays (radial +x̂; tangential ±ŷ, ±ẑ) place TWO probes (d₁, d₂)
     and time the PEAK transit: c = (d₂−d₁)/(t_peak₂−t_peak₁). (v1 used single-probe
     threshold crossing — systematically ~6× fast: a fixed small threshold detects the
     Gaussian LEADING EDGE at t≈(d−σ√(2ln(amp/2thr)))/c, not the front; peak-to-peak
     differencing cancels the flank systematics.) Compare against B's prediction averaged
     over the [d₁,d₂] segment. Repeat at two resolutions (N=48, N=64) — grid convergence
     separates LATTICE anisotropy from PHYSICAL cone structure (the M5.7 lesson).

GATE SEMANTICS: PASS = the measurement is trustworthy (A null exact + cubic slope; B soft
axis radial; C matches B within tolerance at both resolutions). The isotropy VERDICT itself
(isotropic vs anisotropic cone, ratio) is the reported RESULT, not the pass criterion —
an anisotropic cone is a finding to take to Duda (faithful-kinetic 5d re-check + whether
the tilt/EM cone shares it), not a script failure. Ratios are the deliverable; absolute
speed maps to production via c_amrs².

USAGE:
    python -m openwave.xperiments.m5_liquid_crystal.research.sandbox_v8.m5_8_c_isotropy_gate
"""
import numpy as np

from openwave.xperiments.m5_liquid_crystal.research.sandbox_v6.m5_6_2a_biaxial_hedgehog import (
    DELTA, D, L, RC, RHOC, central, commf, frob2, matmul, normalize,
)
from openwave.xperiments.m5_liquid_crystal.research.sandbox_v6.m5_6_2b_biaxial_evolution import (
    Gx, force_and_U, frob_inner,
)

DT = 0.015
N_STEPS = 160                    # t_max=2.4 — covers the slowest predicted d₂ peak + margin
AMP = 0.05                       # seed amplitude (EOM is linear with use_C=False — shape only)
SIGMA = 0.45                     # seed width (≥1.7h at N=48, ≥2.4h at N=64)
R_SEED = 3.0                     # blob center on +x axis (bulk: r>2RC, ρ>2ρc)
D1, D2 = 1.0, 1.6                # two-probe transit distances; d₂ peak lands BEFORE the
#                                  core-mask scatter (~t=1.7) + boundary reflections (~t=2.0)
TOL_C = 0.25                     # measured-vs-predicted speed tolerance (dispersion slack)


def build_bg(n):
    """m5_6_2a build_frame + m5_6_2b build_bg, parameterized by resolution n (same L/RC/ρc/δ)."""
    xs = np.linspace(-L, L, n)
    h = xs[1] - xs[0]
    X, Y, Z = np.meshgrid(xs, xs, xs, indexing="ij")
    r = np.sqrt(X**2 + Y**2 + Z**2)
    rho = np.sqrt(X**2 + Y**2)
    rhat = normalize(np.stack([X, Y, Z], -1) / np.sqrt(X**2 + Y**2 + Z**2 + RC**2)[..., None])
    azim = np.stack([-Y, X, np.zeros_like(X)], -1) / np.sqrt(rho**2 + 1e-12)[..., None]
    s = np.clip(rho / RHOC, 0.0, 1.0)
    shrink = (s * s * (3.0 - 2.0 * s))[..., None]
    ePhi = azim * shrink
    ePhi = ePhi - np.einsum("...i,...i->...", ePhi, rhat)[..., None] * rhat
    eTheta = np.cross(ePhi, rhat)
    O = np.stack([rhat, eTheta, ePhi], axis=-1)
    OT = np.swapaxes(O, -1, -2)
    Mbg = matmul(matmul(O, np.broadcast_to(D, O.shape)), OT)
    Mpsi = matmul(matmul(O, np.broadcast_to(Gx @ D - D @ Gx, O.shape)), OT)
    Mmu = [central(Mbg, ax, h) for ax in range(3)]
    P = [commf(Mmu[ax], Mpsi) for ax in range(3)]
    K = 4.0 * sum(frob2(P[ax]) for ax in range(3))
    pairs = [(0, 1), (0, 2), (1, 2)]
    C = {p: commf(Mmu[p[0]], Mmu[p[1]]) for p in pairs}
    interior = np.zeros(r.shape, bool)
    interior[2:-2, 2:-2, 2:-2] = True
    geom = (r > 2 * RC) & (rho > RHOC) & interior
    return dict(n=n, xs=xs, h=h, K=K, P=P, C=C, pairs=pairs, r=r, rho=rho, rhat=rhat, geom=geom)


def gram_tensor(bg):
    """Gram_μν = P_μ•P_ν per voxel ((...,3,3) sym) — the twist-cone tensor; c²(k̂)=4(1−k̂Gk̂/T)."""
    P = bg["P"]
    gram = np.zeros(bg["r"].shape + (3, 3))
    for mu in range(3):
        for nu in range(mu, 3):
            g = frob_inner(P[mu], P[nu])
            gram[..., mu, nu] = g
            gram[..., nu, mu] = g
    return gram


def c2_along(gram, T, khat):
    """Principal-symbol speed² along unit k̂ (per voxel)."""
    q = np.einsum("...ab,a,b->...", gram, khat, khat)
    return 4.0 * (1.0 - q / np.maximum(T, 1e-30))


# ----------------------------------------------------------------------------------
# Part A — uniform-vacuum null: the EXACT production force (numpy mirror of
# engine2_pde.compute_curvature_flux + the div in evolve_M), V off, 3×3 block.
# ----------------------------------------------------------------------------------
def production_force(M, h):
    Mmu = [central(M, ax, h) for ax in range(3)]
    G = []
    for al in range(3):
        g = np.zeros_like(M)
        for nu in range(3):
            if nu != al:
                g += commf(commf(Mmu[al], Mmu[nu]), Mmu[nu])
        G.append(8.0 * g)
    return sum(central(G[ax], ax, h) for ax in range(3))


def part_A():
    n = 32
    xs = np.linspace(-L, L, n)
    h = xs[1] - xs[0]
    X, Y, Z = np.meshgrid(xs, xs, xs, indexing="ij")
    f1 = np.exp(-((X - 0.8) ** 2 + Y**2 + Z**2) / 2.0)
    f2 = np.exp(-(X**2 + (Y - 0.8) ** 2 + Z**2) / 2.0)        # overlapping blobs
    E1 = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]], float)   # symmetric, [E1,E2]≠0
    E2 = np.array([[0, 0, 1], [0, 0, 0], [1, 0, 0]], float)
    assert frob2(commf(E1, E2)) > 0
    M0 = np.broadcast_to(D, X.shape + (3, 3))                  # uniform biaxial vacuum

    print("\n[A] uniform-vacuum null (production F² force, V off)")
    eps = 1e-2
    F_single = production_force(M0 + eps * f1[..., None, None] * E1, h)
    fmax_single = np.sqrt(frob2(F_single)).max()
    print(f"    single-direction δM=ε·f·E:        max‖force‖ = {fmax_single:.3e}  (expect EXACT 0)")
    fmaxes = []
    eps_list = [1e-3, 1e-2, 1e-1]
    for ep in eps_list:
        dM = ep * (f1[..., None, None] * E1 + f2[..., None, None] * E2)
        fmaxes.append(np.sqrt(frob2(production_force(M0 + dM, h))).max())
    slope = np.polyfit(np.log(eps_list), np.log(fmaxes), 1)[0]
    for ep, fm in zip(eps_list, fmaxes):
        print(f"    two-direction ε={ep:.0e}:             max‖force‖ = {fm:.3e}")
    print(f"    log-log slope = {slope:.3f}  (expect 3 — force is CUBIC in the perturbation)")
    ratio = fmax_single / max(fmaxes[1], 1e-300)
    ok = ratio < 1e-10 and 2.8 < slope < 3.2
    print(f"    ⇒ NO linear-order propagation in the uniform vacuum — 'vacuum c' is undefined")
    print(f"      at linear order; waves ride on background gradients (defect dressings).  "
          f"[{'OK' if ok else 'FAIL'}]")
    return ok


# ----------------------------------------------------------------------------------
# Part B — static cone tensor on the far-field hedgehog shell.
# ----------------------------------------------------------------------------------
def part_B(bg, label=""):
    gram = gram_tensor(bg)
    T = np.einsum("...aa->...", gram)
    shell = bg["geom"] & (bg["r"] > 2.5) & (bg["r"] < 4.5) & (bg["rho"] > 2 * RHOC)
    lam, vec = np.linalg.eigh(gram[shell])                     # ascending eigenvalues
    Ts = lam.sum(1)
    c2_fast = 4.0 * (1.0 - lam[:, 0] / Ts)                     # k̂ along SOFT Gram axis
    c2_slow = 4.0 * (1.0 - lam[:, 2] / Ts)                     # k̂ along STIFF Gram axis
    cone = np.sqrt(c2_fast / np.maximum(c2_slow, 1e-30))
    v_soft = vec[..., 0]                                       # eigenvector of λ_min
    align = np.abs(np.einsum("...i,...i->...", v_soft, bg["rhat"][shell]))

    print(f"\n[B{label}] static twist-cone tensor, far-field shell r∈(2.5,4.5)  "
          f"({shell.sum()} voxels, N={bg['n']})")
    print(f"    Gram eigen-fractions λ/T (mean): {np.round((lam / Ts[:, None]).mean(0), 3)}"
          f"   (soft≈0 ⇒ P_radial≈0)")
    print(f"    soft-axis · r̂ alignment: mean={align.mean():.4f}  (expect →1: the fast direction"
          f" is RADIAL)")
    print(f"    cone ratio c_max/c_min: mean={cone.mean():.3f} ± {cone.std():.3f}"
          f"   (isotropic ⟺ 1; pure-tangential-Gram prediction ≈ √2 = 1.414)")
    rs = bg["r"][shell]
    bins = np.linspace(2.5, 4.5, 5)
    means = [cone[(rs >= a) & (rs < b)].mean() for a, b in zip(bins[:-1], bins[1:])]
    print(f"    r-binned cone ratio {np.round(bins[:-1], 1)}→: {np.round(means, 3)}"
          f"   (flat ⇒ scale-free cone)")
    ok = align.mean() > 0.9
    return ok, gram, T, float(cone.mean())


# ----------------------------------------------------------------------------------
# Part C — dynamic two-probe peak-transit timing vs the Part-B prediction (use_C=False
# ⇒ exactly the homogeneous linear twist operator: no C-drive contamination).
# c = (d₂−d₁)/(t_peak₂−t_peak₁) — flank/threshold systematics cancel between probes.
# ----------------------------------------------------------------------------------
def peak_time(series, t_lo=0.0):
    """Sub-step peak time of |ψ|(t): descend the initial seed-tail decay to its first
    local minimum, THEN take the argmax (parabolic-refined). Without the descent, a probe
    whose arriving packet is weaker than its initial Gaussian-tail value (the φ̂ probes)
    returns t≈0 — the tail's decay, not the packet.

    Returns (t_peak, at_edge): at_edge=True flags an argmax pinned to the window edge."""
    s = np.abs(np.asarray(series))
    i0 = max(1, int(round(t_lo / DT)))
    while i0 < len(s) - 1 and s[i0 + 1] <= s[i0]:              # descend the initial decay
        i0 += 1
    i = i0 + int(np.argmax(s[i0:]))
    if i >= len(s) - 1 or i <= 0:
        return i * DT, True
    a, b, c = s[i - 1], s[i], s[i + 1]
    denom = a - 2 * b + c
    off = 0.5 * (a - c) / denom if abs(denom) > 1e-30 else 0.0
    return (i + off) * DT, False


def part_C(bg, gram, T):
    n, xs, h = bg["n"], bg["xs"], bg["h"]
    idx = lambda v: int(round((v + L) / h))
    ic = (idx(R_SEED), idx(0.0), idx(0.0))
    c0 = np.array([xs[ic[0]], xs[ic[1]], xs[ic[2]]])           # snapped blob center
    dirs = {"+x̂ (radial)": np.array([1.0, 0, 0]), "+ŷ (tang φ)": np.array([0, 1.0, 0]),
            "−ŷ (tang φ)": np.array([0, -1.0, 0]), "+ẑ (tang θ)": np.array([0, 0, 1.0]),
            "−ẑ (tang θ)": np.array([0, 0, -1.0])}
    rays = {}
    for name, d in dirs.items():
        pis, dists = [], []
        for dp in (D1, D2):
            p = c0 + dp * d
            pi = (idx(p[0]), idx(p[1]), idx(p[2]))
            pis.append(pi)
            dists.append(np.linalg.norm(np.array([xs[pi[0]], xs[pi[1]], xs[pi[2]]]) - c0))
        # predicted speed: ⟨c(k̂=d̂)⟩ over the measured [d₁,d₂] segment (principal symbol, B)
        cs = []
        for t in np.linspace(dists[0], dists[1], 7):
            q = c0 + t * d
            qi = (idx(q[0]), idx(q[1]), idx(q[2]))
            cs.append(np.sqrt(max(c2_along(gram[qi], T[qi], d), 0.0)))
        rays[name] = dict(pis=pis, dists=dists, c_pred=float(np.mean(cs)),
                          series=[[], []])

    X, Y, Z = np.meshgrid(xs, xs, xs, indexing="ij")
    seed = AMP * np.exp(-((X - c0[0]) ** 2 + (Y - c0[1]) ** 2 + (Z - c0[2]) ** 2)
                        / (2 * SIGMA**2))
    K = bg["K"]
    active = (K > 1e-6 * K.max()) & bg["geom"]
    inv2K = np.where(active, 1.0 / (2.0 * K + 1e-30), 0.0)
    psi = (seed * active).copy()
    psi_prev = psi.copy()                                       # zero initial velocity
    Hs = []
    for step in range(N_STEPS):
        force, _ = force_and_U(psi, bg, use_C=False)
        psi_new = (2 * psi - psi_prev + DT**2 * force * inv2K) * active
        psi_prev, psi = psi, psi_new
        for ray in rays.values():
            ray["series"][0].append(psi[ray["pis"][0]])
            ray["series"][1].append(psi[ray["pis"][1]])
        if step % 40 == 0:
            _, U = force_and_U(psi, bg, use_C=False)
            pt = (psi - psi_prev) / DT
            Hs.append(float(np.sum((K * pt**2 + U)[active])))
    drift = (max(Hs) - min(Hs)) / abs(Hs[0])

    print(f"\n[C] two-probe peak-transit timing, N={n}  (blob at r={R_SEED}, probes at "
          f"d₁={D1}, d₂={D2}, use_C=False)")
    print(f"    {'direction':<14} {'t_pk1':>6} {'t_pk2':>6} {'c_meas':>7} {'c_pred(B)':>9} "
          f"{'ratio':>6}")
    # PASS uses radial + θ̂(±ẑ) only — straight-chord timing is valid where rays don't
    # bend. φ̂(±ŷ) is the SLOW axis of an anisotropic cone (fast axis radial): rays
    # launched along φ̂ refract outward, so its chord numbers are reported but
    # informational (the static Part-B tensor is the primary φ̂ measurement).
    ok, meas, preds = True, {}, {}
    for name, ray in rays.items():
        t1, e1 = peak_time(ray["series"][0])
        t2, e2 = peak_time(ray["series"][1], t_lo=t1)           # d₂ peak comes after d₁'s
        edge = e1 or e2
        cm = (ray["dists"][1] - ray["dists"][0]) / (t2 - t1) if t2 > t1 else 0.0
        rat = cm / ray["c_pred"] if ray["c_pred"] > 0 else 0.0
        meas[name], preds[name] = cm, ray["c_pred"]
        gating = "ŷ" not in name
        if gating:
            ok = ok and not edge and abs(rat - 1.0) < TOL_C
        flag = ("  ⚠️ edge" if edge else "") + (
            "  ⚠️" if abs(rat - 1.0) >= TOL_C and gating else
            ("  (info — refracting axis)" if not gating else ""))
        print(f"    {name:<14} {t1:>6.3f} {t2:>6.3f} {cm:>7.3f} {ray['c_pred']:>9.3f} "
              f"{rat:>6.2f}{flag}")
    rad = meas["+x̂ (radial)"]
    th = np.mean([meas[k] for k in meas if "ẑ" in k])
    th_p = np.mean([preds[k] for k in preds if "ẑ" in k])
    aniso = rad / th if th > 0 else 0.0
    aniso_p = preds["+x̂ (radial)"] / th_p
    print(f"    radial/θ̂ speed ratio: measured = {aniso:.3f}  vs symbol-predicted = "
          f"{aniso_p:.3f}   energy drift = {100 * drift:.2f}%")
    return ok, aniso


def main():
    print("=" * 78)
    print("M5.8 — c-ISOTROPY GATE (SR emergence): vacuum null + twist-cone tensor + fronts")
    print(f"  δ={DELTA}  r_c={RC}  ρ_c={RHOC}  dt={DT}  (Eq.18 F² dynamics, V off, 3×3 block)")
    print("=" * 78)

    okA = part_A()

    bg48 = build_bg(48)
    okB, gram48, T48, cone48 = part_B(bg48, label=", N=48")
    okC48, aniso48 = part_C(bg48, gram48, T48)

    bg64 = build_bg(64)
    okB64, gram64, T64, cone64 = part_B(bg64, label=", N=64")
    okC64, aniso64 = part_C(bg64, gram64, T64)

    conv = abs(aniso64 - aniso48) / max(aniso48, 1e-30)
    print(f"\n[D] grid convergence: measured radial/θ̂ ratio  N48={aniso48:.3f}  "
          f"N64={aniso64:.3f}  (Δ={100 * conv:.1f}%, want <10% — lattice vs physical separation)")

    verdict_iso = abs(cone64 - 1.0) < 0.05
    ok = okA and okB and okB64 and okC48 and okC64 and conv < 0.10
    print("\n" + "=" * 78)
    print("RESULT — the SR question restated for Eq.18:")
    print("  (1) the uniform biaxial vacuum has NO linear c (transparent; cubic-only force);")
    print(f"  (2) the twist cone on the defect background is "
          f"{'ISOTROPIC' if verdict_iso else 'ANISOTROPIC'} — static cone ratio "
          f"{cone64:.3f} (radial fast / φ̂ slow), dynamic radial/θ̂ {aniso64:.3f};")
    print("  (3) interpretation (Duda flag if anisotropic): the cone is tied to the local")
    print("      defect frame (Gram(P) soft axis = r̂) — check whether the tilt/EM sector")
    print("      shares it (M5.6.4) + whether the faithful kinetic (5d) reshapes it.")
    print("MEASUREMENT GATE: " + ("PASS" if ok else "PARTIAL — inspect the failing metric above"))
    print("=" * 78)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
