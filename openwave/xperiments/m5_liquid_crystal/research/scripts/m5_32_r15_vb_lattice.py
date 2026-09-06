"""M5.32 R15-V-b: the author's floor witness V1 on the lattice.

Equations first (the conventions of m5_32_r15_vh_symbolic.py).  The witness field:
    TWIST INSIDE   M(x) = L(n(x), chi(r)) R_12(k z) D R_12(k z)^T L^T
    TWIST AFTER    M(x) = R_12(k z) L(n(x), chi(r)) D L^T R_12(k z)^T
    L(n, chi) = 1 + sinh(chi) K(n) + (cosh(chi) - 1) K2(n)   (the boost along the radial unit vector n,
                 the m5_21_8_b_lattice.dressed form),  chi(r) = chi_0 exp(-r^2 / (2 w^2)), chi_0 = 0.5, w = 2
    D = diag(g, 1, delta, delta) (degenerate) or diag(g, 1, delta, 0) (certified);  k in {0.5, 1, 2}
Reads (static, E-orientation, h^3-weighted over the whole box, the sym stencil with one-sided edges):
    E_eta = E_u = 4 sum <F_ij, F_ij>_eta   (the certified -4 I1 read by m5_21_3_a_4d.e_parts)
    E_h   = 4 x I1_h                        (the author's -4 I1^h read by m5_32_terms_ext, C15.i1h_static)
    DeltaE(k) = E[witness](k) - E[dressed, k = 0]   (the pure twist alone has E = 0 exactly: one gradient direction)
The author's 64^3 numbers (V1, 09-05): DeltaE_eta = -114 / -388 / -1091 and DeltaE_h = +183 / +594 / +1549
for k = 0.5 / 1 / 2 (twist inside); twist after: +160 / +636 / +2473, "both forms".
Claims: V1 DeltaE_eta < 0 and growing with k (twist inside); V2 DeltaE_h > 0; V3 twist after positive in both
forms (the jet rung: equal only at leading order in chi); V5 the same signs on the relaxed L_cert hedgehog
(R10 n32 L48) dressed and twisted the same way.
Grids: n64 L48 (h 0.75, the author's), n32 L48 (h 1.5) and n48 L72 (h 1.5, the box-size control).

usage: python3 m5_32_r15_vb_lattice.py
"""
import sys
ARGS = list(sys.argv[1:])
import os, json, time
import numpy as np
import m5_32_r15_common as C15

INS4, C13, B8 = C15.INS4, C15.C13, C15.B8
log = C15.log
G1, G2, G3 = B8.G1, B8.G2, B8.G3


def boost_field(cfg, chi0=0.5, w=2.0):
    n, h = cfg["n"], cfg["h"]
    X, Y, Z = INS4.coords(n, h)
    R = np.sqrt(X * X + Y * Y + Z * Z)
    assert R.min() > 0, "a grid point at the origin: the radial unit vector is undefined there"
    chi = chi0 * np.exp(-R * R / (2 * w * w))
    nx, ny, nz = X / R, Y / R, Z / R
    K = np.zeros(X.shape + (4, 4))
    K[..., 0, 1], K[..., 0, 2], K[..., 0, 3] = nx, ny, nz
    K[..., 1, 0], K[..., 2, 0], K[..., 3, 0] = nx, ny, nz
    K2 = np.zeros_like(K)
    K2[..., 0, 0] = 1.0
    for i, a in enumerate((nx, ny, nz)):
        for j, b in enumerate((nx, ny, nz)):
            K2[..., 1 + i, 1 + j] = a * b
    Lb = np.eye(4)[None, None, None] + np.sinh(chi)[..., None, None] * K + (np.cosh(chi) - 1.0)[..., None, None] * K2
    return Lb, Z


def twist_field(Z, k):
    return B8.rot_field(G3, k * Z)              # G3 = the (1,2) rotation generator


def conj(Q, M):
    return np.einsum("...ab,...bc,...dc->...ad", Q, M, Q)


def reads(M, cfg):
    e_u = float(INS4.e_parts(M, cfg)[0])
    e_h = float(C15.i1h_static(M, cfg))
    return e_u, e_h


def run_grid(n, L, base, label, cfg, Mbase):
    Lb, Z = boost_field(cfg)
    Md = INS4.sym4(conj(Lb, Mbase))
    ed = reads(Md, cfg)
    e0 = reads(Mbase, cfg)
    rec = {"n": n, "L": L, "h": cfg["h"], "base": label, "E_base_[eta,h]": e0, "E_dressed_k0_[eta,h]": ed, "k": {}}
    log(f"{label} n{n} L{L:g}: base E_eta {e0[0]:.4f} E_h {e0[1]:.4f}; dressed (k 0) E_eta {ed[0]:.4f} E_h {ed[1]:.4f}")
    for k in (0.5, 1.0, 2.0):
        Rt = twist_field(Z, k)
        Mt = INS4.sym4(conj(Rt, Mbase))
        et = reads(Mt, cfg)
        Min = INS4.sym4(conj(Lb, conj(Rt, Mbase)))
        Maf = INS4.sym4(conj(Rt, conj(Lb, Mbase)))
        ein, eaf = reads(Min, cfg), reads(Maf, cfg)
        rec["k"][f"{k:g}"] = {"E_twist_alone_[eta,h]": et,
                              "inside_[eta,h]": ein, "after_[eta,h]": eaf,
                              "DeltaE_inside_[eta,h]": [ein[0] - ed[0], ein[1] - ed[1]],
                              "DeltaE_after_[eta,h]": [eaf[0] - ed[0], eaf[1] - ed[1]]}
        log(f"   k {k:g}: twist alone {et[0]:.3e}/{et[1]:.3e}; inside dE_eta {ein[0] - ed[0]:+.3f} dE_h {ein[1] - ed[1]:+.3f}; "
            f"after dE_eta {eaf[0] - ed[0]:+.3f} dE_h {eaf[1] - ed[1]:+.3f}")
    return rec


def main():
    out = {"rung": "R15-V-b", "author_numbers": {"inside_eta": [-114, -388, -1091], "inside_h": [183, 594, 1549], "after_both": [160, 636, 2473], "k": [0.5, 1, 2]},
           "profile": "chi(r) = 0.5 exp(-r^2 / 8), radial boost; twist (1,2) psi = k z", "runs": []}
    for n, L in ((32, 48.0), (64, 48.0), (48, 72.0)):
        for vac in ("degenerate", "certified"):
            cfg = C15.cfg_dd(n, L) if vac == "degenerate" else C15.cfg_cert(n, L)
            d = INS4.vac4(cfg)
            Mbase = np.broadcast_to(d, (n, n, n, 4, 4)).copy()
            out["runs"].append(run_grid(n, L, "vacuum", f"{vac} vacuum", cfg, Mbase))
    # V5: the relaxed L_cert hedgehog (R10, n32 L48)
    Mhh, cfg, src = C13.seed_hedgehog(32, 48)
    out["runs"].append(run_grid(32, 48.0, "hedgehog", "relaxed L_cert hedgehog (R10 n32 L48)", cfg, Mhh))
    json.dump(out, open(os.path.join(C15.DATA, "m5_32_r15_vb_lattice.json"), "w"), indent=1)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
    ks = [0.5, 1.0, 2.0]
    for rec in out["runs"]:
        lab = f"{rec['base']} n{rec['n']} L{rec['L']:g}"
        ax[0].plot(ks, [rec["k"][f"{k:g}"]["DeltaE_inside_[eta,h]"][0] for k in ks], marker="o", label=lab + " eta")
        ax[0].plot(ks, [rec["k"][f"{k:g}"]["DeltaE_inside_[eta,h]"][1] for k in ks], marker="s", ls="--", label=lab + " h")
        ax[1].plot(ks, [rec["k"][f"{k:g}"]["DeltaE_after_[eta,h]"][0] for k in ks], marker="o", label=lab + " eta")
        ax[1].plot(ks, [rec["k"][f"{k:g}"]["DeltaE_after_[eta,h]"][1] for k in ks], marker="s", ls="--", label=lab + " h")
    ax[0].plot(ks, out["author_numbers"]["inside_eta"], "k:", marker="x", label="author eta (64^3)")
    ax[0].plot(ks, out["author_numbers"]["inside_h"], "k-.", marker="+", label="author h (64^3)")
    ax[1].plot(ks, out["author_numbers"]["after_both"], "k:", marker="x", label="author both (64^3)")
    ax[0].set_title("twist INSIDE the boost dressing: DeltaE(k)"); ax[1].set_title("twist AFTER the dressing: DeltaE(k)")
    for a in ax:
        a.axhline(0, color="gray", lw=0.8); a.set_xlabel("k"); a.grid(alpha=0.3); a.legend(fontsize=6)
    fig.suptitle("M5.32 R15-V-b: the floor witness on the lattice (E-orientation, E = 4 I1 / 4 I1^h)")
    fig.tight_layout()
    fig.savefig(os.path.join(C15.PLOTS, "m5_32_r15_vb_lattice.png"), dpi=110)
    log("wrote data/m5_32_r15_vb_lattice.json + plots/m5_32_r15_vb_lattice.png")


if __name__ == "__main__":
    main()
