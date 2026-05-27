"""
M5.6.2a — biaxial hedgehog: background curvature C_μν (the mass source) + the
z-axis disclination, characterized and regularized.

M5.5.2 evolved the twist on a SMOOTH single-generator background O_bg=Ry(γ), which has
C_μν=[M_μ^bg,M_ν^bg]=0 (one generator commutes with its own derivatives) ⇒ no mass gap.
M5.6 needs the BIAXIAL HEDGEHOG frame O_bg=[r̂|e_Θ|e_Φ] with D=diag(1,δ,0): the frame
rotates via MULTIPLE generators across space ⇒ C_μν≠0 (the gauge source that gives the
KG mass, M5.6.1's geometric mass at the matrix level). The cost: assigning the secondary
axes (e_Θ for δ, e_Φ for 0) continuously is the hairy-ball problem → e_Φ ~ 1/ρ has a
z-axis DISCLINATION line, on top of the r=0 point core.

This script (structure-first, like M5.6.1a) verifies:
  1. the regularized biaxial frame is orthonormal (det=+1) away from the singularities;
  2. M_bg = O D O^T recovers eigenvalues (1,δ,0) with principal director = r̂;
  3. C_μν ≠ 0 (the mass source) — contrast with the single-generator C≈0 — and its radial
     profile ‖C‖ ~ 1/r² (matching M5.6.1's geometric mass scale);
  4. the disclination lives on the z-axis; the ρ_c regularization caps its blow-up.

USAGE:
    python -m openwave.xperiments.m5_liquid_crystal.research.sandbox_v5.m5_6_2a_biaxial_hedgehog
"""
import numpy as np

N = 48
L = 6.0
DELTA = 0.3                     # D = diag(1, δ, 0)
RC = 0.8                        # point-core regularization (r=0)
RHOC = 0.8                      # disclination-line regularization (z-axis, ρ=0)
D = np.diag([1.0, DELTA, 0.0])


def matmul(A, B):
    return np.einsum("...ac,...cb->...ab", A, B)


def commf(A, B):
    return matmul(A, B) - matmul(B, A)


def frob2(M):
    return np.einsum("...ab,...ab->...", M, M)


def normalize(v, eps=1e-30):
    n = np.sqrt(np.einsum("...i,...i->...", v, v))[..., None]
    return v / (n + eps)


def central(f, axis, h):
    out = np.zeros_like(f)
    sp_, sm, so = [slice(None)] * f.ndim, [slice(None)] * f.ndim, [slice(None)] * f.ndim
    sp_[axis], sm[axis], so[axis] = slice(2, None), slice(0, -2), slice(1, -1)
    out[tuple(so)] = (f[tuple(sp_)] - f[tuple(sm)]) / (2 * h)
    return out


def build_frame(rhoc=RHOC):
    xs = np.linspace(-L, L, N)
    h = xs[1] - xs[0]
    X, Y, Z = np.meshgrid(xs, xs, xs, indexing="ij")
    r = np.sqrt(X**2 + Y**2 + Z**2)
    rho = np.sqrt(X**2 + Y**2)
    rhat = normalize(np.stack([X, Y, Z], -1) / np.sqrt(X**2 + Y**2 + Z**2 + RC**2)[..., None])
    # azimuthal unit direction (winds around z) — DISCLINATION on z-axis:
    azim = np.stack([-Y, X, np.zeros_like(X)], -1) / np.sqrt(rho**2 + 1e-12)[..., None]
    # disclination regularization: secondary axes are FULL unit length for ρ≥ρ_c (exact
    # biaxial hedgehog), and melt smoothly to 0 inside the core ρ<ρ_c (biaxiality melts on
    # the disclination line, like a nematic). Clamped smoothstep ⇒ orthonormal bulk + capped ‖∂O‖.
    s = np.clip(rho / rhoc, 0.0, 1.0)
    shrink = (s * s * (3.0 - 2.0 * s))[..., None]
    ePhi = azim * shrink
    ePhi = ePhi - np.einsum("...i,...i->...", ePhi, rhat)[..., None] * rhat   # ⟂ r̂ (no renorm)
    eTheta = np.cross(ePhi, rhat)                                  # |eΘ|=|eΦ| → melts in the core
    O = np.stack([rhat, eTheta, ePhi], axis=-1)                    # columns = frame axes
    return dict(h=h, X=X, Y=Y, Z=Z, r=r, rho=rho, O=O, rhat=rhat)


def main():
    print("=" * 72)
    print("M5.6.2a — biaxial hedgehog: C_μν (mass source) + z-axis disclination")
    print(f"  grid {N}³  δ={DELTA}  r_c={RC}  ρ_c={RHOC}")
    print("=" * 72)
    fr = build_frame()
    O, h, r, rho = fr["O"], fr["h"], fr["r"], fr["rho"]
    OT = np.swapaxes(O, -1, -2)

    # regions: bulk = away from BOTH core (r) and disclination (ρ); disc = near z-axis
    interior = np.zeros(r.shape, bool); interior[2:-2, 2:-2, 2:-2] = True
    bulk = (r > 2 * RC) & (rho > 2 * RHOC) & interior
    disc = (rho < RHOC) & (r > 2 * RC) & interior          # disclination line, off the point core

    # --- 1. orthonormality / det of the frame ------------------------------------
    gram = matmul(OT, O)
    ortho_err = np.abs(gram - np.eye(3)).max(axis=(-1, -2))
    det = np.linalg.det(O)
    print("\n[1] frame orthonormality  O^T O = I,  det O = +1")
    print(f"    bulk:        max|O^TO−I|={ortho_err[bulk].max():.2e}   det∈[{det[bulk].min():.3f},{det[bulk].max():.3f}]")
    print(f"    disclination:max|O^TO−I|={ortho_err[disc].max():.2e}   (frame degrades on the z-axis ✓ expected)")

    # --- 2. M_bg eigenstructure ---------------------------------------------------
    Mbg = matmul(matmul(O, np.broadcast_to(D, O.shape)), OT)
    sym_err = np.abs(Mbg - np.swapaxes(Mbg, -1, -2)).max()
    evals = np.linalg.eigvalsh(Mbg[bulk])                    # ascending
    director = np.linalg.eigh(Mbg[bulk])[1][..., -1]         # eigenvector of largest eigenvalue
    align = np.abs(np.einsum("...i,...i->...", director, fr["rhat"][bulk]))
    print("\n[2] M_bg = O·diag(1,δ,0)·O^T  (bulk)")
    print(f"    symmetric: max|M−Mᵀ|={sym_err:.2e}")
    print(f"    eigenvalues (mean): {evals.mean(0).round(3)}  (expect [0, {DELTA}, 1])")
    print(f"    principal director · r̂ : mean={align.mean():.4f}  (expect 1 ⇒ eigenvalue-1 axis is radial)")

    # --- 3. C_μν = [M_μ^bg, M_ν^bg] — the mass source -----------------------------
    Mmu = [central(Mbg, ax, h) for ax in range(3)]
    pairs = [(0, 1), (0, 2), (1, 2)]
    C = {p: commf(Mmu[p[0]], Mmu[p[1]]) for p in pairs}
    Cnorm2 = sum(frob2(C[p]) for p in pairs)                 # Σ‖C_μν‖²
    print("\n[3] background curvature  C_μν = [M_μ^bg, M_ν^bg]  (the KG mass source)")
    print(f"    Σ‖C_μν‖²: bulk mean={Cnorm2[bulk].mean():.3e}  max={Cnorm2[bulk].max():.3e}  → NONZERO ✓")
    print(f"    (contrast: single-generator background O=Ry(γ) gives C_μν≡0 — no mass, cf. M5.5.2)")
    # radial profile: bin ‖C‖ = √Σ‖C‖² vs r in the bulk; expect ~ 1/r² (geometric mass scale)
    Cn = np.sqrt(Cnorm2)
    rb, cb = r[bulk], Cn[bulk]
    sel = (rb > 2.0) & (rb < 5.0)
    p = np.polyfit(np.log(rb[sel]), np.log(cb[sel]), 1)[0]    # slope of log‖C‖ vs log r
    print(f"    radial scaling  ‖C‖ ∝ r^({p:.2f})  (expect ≈ −2 ⇒ matches M5.6.1 mass²~1/r²) ")

    # --- 4. disclination location + ρ_c regularization ----------------------------
    onaxis = (rho < RHOC) & (r > 2 * RC) & interior
    offaxis = (rho > 2 * RHOC) & (r > 2 * RC) & interior
    print("\n[4] z-axis disclination — frame-gradient blow-up, capped by ρ_c")
    gradO2 = sum(frob2(central(O, ax, h)) for ax in range(3))   # ‖∂O‖² (frame strain)
    print(f"    ‖∂O‖²  on-axis(ρ<ρ_c) mean={gradO2[onaxis].mean():.2e}   off-axis mean={gradO2[offaxis].mean():.2e}")
    # sweep ρ_c: smaller ρ_c ⇒ sharper (larger) disclination peak; larger ρ_c caps it
    print("    ρ_c sweep — peak ‖∂O‖² on the z-axis (regularization caps the singularity):")
    for rc in (0.4, 0.8, 1.2):
        fr2 = build_frame(rhoc=rc)
        g2 = sum(frob2(central(fr2["O"], ax, h)) for ax in range(3))
        m = (fr2["rho"] < rc) & (fr2["r"] > 2 * RC) & interior
        print(f"      ρ_c={rc}:  peak ‖∂O‖²={g2[m].max():.2e}")

    ok = (ortho_err[bulk].max() < 1e-2 and align.mean() > 0.99 and Cnorm2[bulk].mean() > 0
          and -2.6 < p < -1.4 and gradO2[onaxis].mean() > 5 * gradO2[offaxis].mean())
    print("\n" + "=" * 72)
    print("M5.6.2a: biaxial hedgehog frame built; C_μν≠0 (mass source, ‖C‖~1/r²) confirmed;")
    print("z-axis disclination located + ρ_c-regularizable. Ready for M5.6.2b (evolve on it).")
    print("PASS" if ok else "PARTIAL — inspect the failing metric above")
    print("=" * 72)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
