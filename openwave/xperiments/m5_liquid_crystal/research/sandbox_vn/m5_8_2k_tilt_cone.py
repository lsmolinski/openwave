"""
M5.8.2k — N-4: THE EM-TILT CONE CHECK — do the tilt/EM generators ride the
SAME wave cone as the twist?

THE QUESTION (the SR question, upgraded from the backlog 2026-06-07): the
c-isotropy gate (m5_8_c_isotropy_gate.py) measured the TWIST channel (Gx,
the δ–0 block = the QM/KG axis of the M5.6 eigenvalue map): anisotropic
cone, soft Gram axis radial, ratio ≈ √2, scale-free, tied to the defect
frame. The EM sector lives on the TILT generators (the 1-eigenvalue axis →
Maxwell, M5.6.4):
    Gy  rotation about e_Θ — mixes r̂ ↔ e_Φ (the 1–0 block)
    Gz  rotation about e_Φ — mixes r̂ ↔ e_Θ (the 1–δ block)
SAME cone ⇒ local SR survives as cone-relative kinematics (one cone for
matter + light on the defect background); DIFFERENT ⇒ falsifier-candidate
flag — honest + report-worthy either way.

METHOD (the gate's Part B verbatim, per generator): channel field
ψ(x)·M_G(x) with M_G = O(G·D − D·G)Oᵀ; P_μ = [∂_μM_bg, M_G];
Gram_μν = P_μ•P_ν; principal-symbol cone c²(k̂) = 4(1 − k̂·Gram·k̂/T).
The K = 4T identity holds per channel BY CONSTRUCTION (same algebra), so
the symbol formula carries; the geometric prediction P_radial ≈ 0 (far-
field frame depends on angles only) is also channel-independent — the
DISCRIMINATING measurables are the tangential eigen-structure + cone ratio.
Two resolutions (N = 48, 64) per generator — the M5.7 lattice-vs-physics
lesson. CPU-only (numpy); deliverable = the 3-row SAME/DIFFERENT table.

PRE-REGISTERED READING: SAME = soft axis radial (align → 1) at both N for
all three generators AND cone ratios agreeing within a few % AND r-binned
flat. Any generator breaking any leg ⇒ DIFFERENT on that leg — report the
table, no auto-verdict beyond it.

RESULTS (2026-06-07 — N-4 COMPLETE: SHARED CEILING, DIFFERENT CONE SHAPES):
  | generator | eigen-fractions (N=64) | stiff·r̂ | c(r̂)/2 | cone ratio |
  | Gx twist  | (0.000, 0.329, 0.671) | 0.0008 | 1.0000 | 1.743 (N-stable) |
  | Gy tilt   | (0.000, 0.089, 0.911) | 0.0008 | 1.0000 | ~10–11 (near-rank-1) |
  | Gz tilt   | (0.000, 0.000, 1.000) | 0.0008 | 1.0000 | unbounded BY CONSTRUCTION |
  · THE SAME LEG (the SR-relevant statement): all three channels share the
    MAXIMAL cone — the fast direction is RADIAL with the identical ceiling
    c(r̂)/2 = 1.0000 at both N, stiff axis exactly tangential (→0 with N).
    Matter (twist/KG) and light (tilt/EM) signals have the SAME maximal
    speed on the defect background: local SR survives as cone-relative
    kinematics in the radial direction.
  · THE DIFFERENT LEG: the tangential cone SHAPES are channel-dependent —
    twist couples to both tangential gradients (1/3–2/3 split, mild 1.74
    cone); Gy is near-rank-1 (one tangential direction dominates 10:1);
    Gz is rank-1 to ≤1e-3 at BOTH N (consistent with EXACT): its Gram has
    P_θ ∥ P_φ ⇒ the channel CANNOT propagate along one tangential
    direction at all. For rank-1 channels the cone RATIO diverges by
    construction — the converged statement is the eigen-fraction, not the
    ratio (the N-drift of Gy's 10.4→11.5 ratio is the same small-λ_mid
    sensitivity; its fraction 0.088→0.089 is converged).
  · Consistency: the Gx row reproduces the original c-isotropy gate
    (ratio 1.742/1.743, soft axis radial 1.0000) — the rig is validated
    against the prior measurement.
  · Honest scope: static principal-symbol analysis on the far-field
    biaxial-hedgehog background (the gate's Part B); no dynamic transit
    cross-check for the tilt channels (the gate's Part C pattern would
    need the tilt EOM ported — a refinement item, not run here).
  ⇒ REPORT FRAMING: "SAME ceiling / DIFFERENT tangential optics" — light
    and matter share the maximal signal speed on the defect frame; the EM
    sector's tangential propagation is strongly anisotropic, one tilt
    polarization being non-propagating along one tangential direction.

USAGE:  python m5_8_2k_tilt_cone.py
"""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[4]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from openwave.xperiments.m5_liquid_crystal.research.sandbox_v6.m5_6_2a_biaxial_hedgehog import (  # noqa: E402
    D, L, RC, RHOC, central, commf, frob2, matmul, normalize,
)
from openwave.xperiments.m5_liquid_crystal.research.sandbox_v6.m5_6_2b_biaxial_evolution import (  # noqa: E402
    Gx,
)
from openwave.xperiments.m5_liquid_crystal.research.sandbox_v8.m5_8_c_isotropy_gate import (  # noqa: E402
    gram_tensor, part_B,
)

GY = np.array([[0, 0, 1], [0, 0, 0], [-1, 0, 0]], float)    # tilt: r̂↔e_Φ (1–0)
GZ = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 0]], float)    # tilt: r̂↔e_Θ (1–δ)
GENS = (("Gx twist (δ–0, the QM/KG axis)", Gx),
        ("Gy tilt  (1–0, the EM axis)", GY),
        ("Gz tilt  (1–δ, the EM axis)", GZ))


def build_bg_G(n, G):
    """The gate's build_bg with the channel generator parameterized."""
    xs = np.linspace(-L, L, n)
    h = xs[1] - xs[0]
    X, Y, Z = np.meshgrid(xs, xs, xs, indexing="ij")
    r = np.sqrt(X**2 + Y**2 + Z**2)
    rho = np.sqrt(X**2 + Y**2)
    rhat = normalize(np.stack([X, Y, Z], -1)
                     / np.sqrt(X**2 + Y**2 + Z**2 + RC**2)[..., None])
    azim = np.stack([-Y, X, np.zeros_like(X)], -1) / np.sqrt(rho**2 + 1e-12)[..., None]
    s = np.clip(rho / RHOC, 0.0, 1.0)
    shrink = (s * s * (3.0 - 2.0 * s))[..., None]
    ePhi = azim * shrink
    ePhi = ePhi - np.einsum("...i,...i->...", ePhi, rhat)[..., None] * rhat
    eTheta = np.cross(ePhi, rhat)
    O = np.stack([rhat, eTheta, ePhi], axis=-1)
    OT = np.swapaxes(O, -1, -2)
    Mbg = matmul(matmul(O, np.broadcast_to(D, O.shape)), OT)
    Mch = matmul(matmul(O, np.broadcast_to(G @ D - D @ G, O.shape)), OT)
    Mmu = [central(Mbg, ax, h) for ax in range(3)]
    P = [commf(Mmu[ax], Mch) for ax in range(3)]
    K = 4.0 * sum(frob2(P[ax]) for ax in range(3))
    interior = np.zeros(r.shape, bool)
    interior[2:-2, 2:-2, 2:-2] = True
    geom = (r > 2 * RC) & (rho > RHOC) & interior
    return dict(n=n, xs=xs, h=h, K=K, P=P, r=r, rho=rho, rhat=rhat, geom=geom)


def main():
    print("=" * 78)
    print("M5.8.2k — the EM-tilt cone check (the gate's Part B per generator,"
          " N = 48 & 64)")
    print("=" * 78)
    ch_norm = np.linalg.norm(Gx @ D - D @ Gx)
    print("  channel-generator norms ‖G·D − D·G‖:"
          + ",  ".join(f" {name.split()[0]} {np.linalg.norm(G @ D - D @ G):.3f}"
                       for name, G in GENS)
          + f"   (all ≠ 0 — every channel exists; twist ref {ch_norm:.3f})")
    rows = []
    for name, G in GENS:
        for n in (48, 64):
            print(f"\n── {name}, N = {n}")
            bg = build_bg_G(n, G)
            ok, gram, T, cone = part_B(bg, label=f" {name.split()[0]} N{n}")
            # tangential eigen-structure (the discriminator beyond the ratio)
            shell = bg["geom"] & (bg["r"] > 2.5) & (bg["r"] < 4.5) & (bg["rho"] > 2 * RHOC)
            lam, vec = np.linalg.eigh(gram[shell])
            fr = (lam / lam.sum(1)[:, None]).mean(0)
            al = np.abs(np.einsum("...i,...i->...", vec[..., 0],
                                  bg["rhat"][shell])).mean()
            # degenerate-safe invariants: the STIFF axis must be tangential
            # (|v_stiff·r̂| → 0) and the RADIAL direction must hit the shared
            # ceiling c(r̂) → 2 — both well-defined even for rank-1 Gram
            # (where the soft-axis·r̂ readout is an arbitrary in-plane vector)
            stiff = np.abs(np.einsum("...i,...i->...", vec[..., 2],
                                     bg["rhat"][shell])).mean()
            rh = bg["rhat"][shell]
            q = np.einsum("...ab,...a,...b->...", gram[shell], rh, rh)
            ceil = np.sqrt(np.maximum(4.0 * (1.0 - q / lam.sum(1)), 0.0)).mean() / 2.0
            rows.append((name.split()[0], n, fr, al, stiff, ceil, cone, ok))
    print("\n" + "=" * 78)
    print("[N-4 TABLE] cross-generator cone comparison")
    print("=" * 78)
    print("| generator | N | Gram eigen-fractions (soft, mid, stiff) |"
          " soft·r̂ | stiff·r̂ | c(r̂)/2 ceiling | cone ratio |")
    print("| --- | --- | --- | --- | --- | --- | --- |")
    for g, n, fr, al, stiff, ceil, cone, ok in rows:
        deg = " (degenerate soft plane — soft·r̂ unphysical)" if fr[1] < 0.01 else ""
        print(f"| {g} | {n} | ({fr[0]:.3f}, {fr[1]:.3f}, {fr[2]:.3f}) |"
              f" {al:.4f}{deg} | {stiff:.4f} | {ceil:.4f} | {cone:.3f} |")
    print("\n  degenerate-safe reading: stiff·r̂ → 0 (stiff axis tangential)"
          " + c(r̂)/2 → 1 (the SHARED")
    print("  radial ceiling) are the channel-independent legs; the eigen-"
          "fractions + cone ratio are")
    print("  the channel-DEPENDENT shape. A rank-1 channel (mid fraction ≈ 0)"
          " has an unbounded ratio")
    print("  BY CONSTRUCTION — the fraction, not the ratio, is its converged"
          " statement.")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
