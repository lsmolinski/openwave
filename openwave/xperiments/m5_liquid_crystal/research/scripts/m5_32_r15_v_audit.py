"""M5.32 R15-V adversarial audit: the floor witness (radial boost x (1,2) twist),
jets (own sympy) and lattice (own numpy fields on the certified stack).

Independent of the producer's scripts: the jets are built from the definitions
(L_a(chi), R_ij(psi), M = L R D R^T L^T or R L D L^T R^T, F_st = b k
[M_chi, M_psi]_eta, U_G = 4 <F_st, F_st>_G, h = eta + 2 (eta u)(eta u)^T with
u = Q e0), the lattice witness from coords(n, h), the radial boost
1 + sinh(chi) K(n) + (cosh(chi) - 1) K2(n), chi(r) = a exp(-r^2 / (2 w^2)),
the twist rot_field(G3, k Z).  Energies: E_eta = e_parts(M, cfg)[0] (certified
instrument) and E_h by OWN pairwise numpy (u analytic = Q e0 on the witness,
own eigen-solve on the hedgehog), cross-checked against
4 x term_lagrangian(REGISTRY_EXT["I1_h"]).

Modes:  jets | lattice | hedgehog | all      (default all)
Out:    data/m5_32_r15_v_audit.json (merged per mode)
"""
from __future__ import annotations

import sys
ARGV = list(sys.argv)                      # captured before any import wipes it
import os
os.environ.setdefault("OMP_NUM_THREADS", "2")
import importlib.util
import json
import time

import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.dirname(HERE)
DATA = os.path.join(RES, "data")
OUT = os.path.join(DATA, "m5_32_r15_v_audit.json")
T0 = time.time()


def log(m):
    print(f"[{time.time() - T0:8.1f}s] {m}", flush=True)


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


INS4 = _load("m5_21_3_a_4d", "m5_21_3_a_4d.py")
B8 = _load("m5_21_8_b_lattice", "m5_21_8_b_lattice.py")
LAG = _load("m5_32_lagrangian", "m5_32_lagrangian.py")
EXT = _load("m5_32_terms_ext", "m5_32_terms_ext.py")
ETA = np.diag([-1.0, 1.0, 1.0, 1.0])
G, DELTA = 8.0, 0.3


def save(section, payload):
    rec = json.load(open(OUT)) if os.path.exists(OUT) else {}
    rec[section] = payload
    rec["_meta"] = {"script": os.path.basename(__file__), "argv": ARGV,
                    "written": time.strftime("%Y-%m-%d %H:%M:%S")}
    json.dump(rec, open(OUT, "w"), indent=1)


# ====================================================================
# PART A: jets, own sympy
# ====================================================================
gS, dS, d4S, chiS, psiS, bS, kS = sp.symbols("g delta delta4 chi psi b k", real=True)
XI = sp.diag(-1, 1, 1, 1)


def boost_sym(a, chi):
    L = sp.eye(4)
    L[0, 0] = L[a, a] = sp.cosh(chi)
    L[0, a] = L[a, 0] = sp.sinh(chi)
    return L


def rot_sym(i, j, psi):
    R = sp.eye(4)
    R[i, i] = R[j, j] = sp.cos(psi)
    R[i, j] = -sp.sin(psi)
    R[j, i] = sp.sin(psi)
    return R


def witness_sym(a, plane, D, order, chi=chiS, psi=psiS):
    L, R = boost_sym(a, chi), rot_sym(plane[0], plane[1], psi)
    Q = L * R if order == "inside" else R * L
    return Q * D * Q.T, Q


def jets_densities(a, plane, D, order, s_axis=1, t_axis=3):
    """U_eta, U_h (symbolic in chi, psi, b, k) for chi-gradient b on axis s and
    psi-gradient k on axis t; the only nonzero F is F_st."""
    M, Q = witness_sym(a, plane, D, order)
    Mchi, Mpsi = sp.diff(M, chiS), sp.diff(M, psiS)
    F = bS * kS * (Mchi * XI * Mpsi - Mpsi * XI * Mchi)
    u = Q[:, 0]                              # Q e0, the timelike eigenvector of M xi
    hu = XI * u
    h = XI + 2 * hu * hu.T
    w = XI[s_axis, s_axis] * XI[t_axis, t_axis]
    U_eta = 4 * w * sp.trace(XI * F * XI * F.T)
    U_h = 4 * w * sp.trace(h * F * h * F.T)
    return U_eta, U_h, M, u, F


def simp(e):
    return sp.simplify(sp.expand(sp.expand_trig(e).rewrite(sp.exp)))


def num(e, **sub):
    d = {gS: 8, dS: sp.Rational(3, 10), d4S: sp.Rational(3, 10), bS: 1, kS: 1,
         chiS: 0, psiS: 0}
    names = {"g": gS, "delta": dS, "delta4": d4S, "b": bS, "k": kS, "chi": chiS, "psi": psiS}
    d.update({names[k]: v for k, v in sub.items()})
    return float(e.subs(d))


def part_jets():
    out = {}
    D = sp.diag(gS, 1, dS, d4S)
    dvac = {1: 1, 2: dS, 3: 0}
    claimed = lambda a: -8 * (dS - 1) ** 2 * (gS + dvac[a]) ** 2

    # ---- V1 / V2 / V4: twist inside, (1,2) plane, boost axis 1, 2, 3
    log("jets: V1/V2 twist inside")
    inside, c_sym = {}, {}
    for a in (1, 2, 3):
        Ue, Uh, M, u, F = jets_densities(a, (1, 2), D, "inside")
        # eigenvector check: (M xi) u = -g u
        ev = simp((M * XI * u + gS * u).norm() ** 2)
        Ue0 = simp(Ue.subs(psiS, 0))
        Uh0 = simp(Uh.subs(psiS, 0))
        dchi = simp(sp.diff(Ue, chiS))          # chi-derivative at GENERAL psi
        dchi_h = simp(sp.diff(Uh, chiS))
        Ue_psi = simp(Ue.subs(chiS, 0))          # general psi, closed form
        Uh_psi = simp(Uh.subs(chiS, 0))
        c_e = simp(Ue0 / (bS * kS) ** 2)
        c_h = simp(Uh0 / (bS * kS) ** 2)
        c_sym[a] = c_e
        inside[f"a{a}"] = {
            "u_is_eigvec_resid": str(ev),
            "c_eta_psi0_symbolic": str(c_e),
            "c_h_psi0_symbolic": str(c_h),
            "claimed_V1": str(sp.expand(claimed(a))),
            "c_eta_minus_claimed": str(simp(c_e - claimed(a))),
            "c_h_plus_c_eta": str(simp(c_h + c_e)),
            "dU_eta_dchi_general_psi": str(dchi),
            "dU_h_dchi_general_psi": str(dchi_h),
            "U_eta_general_psi_chi0": str(Ue_psi),
            "U_h_general_psi_chi0": str(Uh_psi),
            "c_eta_num_g8_d0.3_d4_0.3": num(c_e),
            "c_h_num_g8_d0.3_d4_0.3": num(c_h),
            "c_eta_num_d4_0": num(c_e, delta4=0),
            "c_eta_num_chi1.3": num(Ue, chi=1.3),
            "U_eta_psi_pi4": num(Ue, psi=sp.pi / 4),
            "U_h_psi_pi4": num(Uh, psi=sp.pi / 4),
            "U_eta_psi_pi4_chi0.5": num(Ue, psi=sp.pi / 4, chi=0.5),
        }
        log(f"  a={a}: c_eta={inside[f'a{a}']['c_eta_num_g8_d0.3_d4_0.3']:.4f} "
            f"c_h={inside[f'a{a}']['c_h_num_g8_d0.3_d4_0.3']:.4f} "
            f"dchi={dchi} pi/4: {inside[f'a{a}']['U_eta_psi_pi4']:.4f}")
    out["V1_V2_inside"] = inside
    avg = 0.5 * (inside["a1"]["c_eta_num_g8_d0.3_d4_0.3"] + inside["a2"]["c_eta_num_g8_d0.3_d4_0.3"])
    out["V1_pi4_average_check"] = {"avg_a1_a2": avg,
                                   "U_eta_pi4_a1": inside["a1"]["U_eta_psi_pi4"],
                                   "U_eta_pi4_a2": inside["a2"]["U_eta_psi_pi4"]}
    # large-g leading form (V4)
    c1 = c_sym[1]
    lead = sp.limit(c1 / gS ** 2, gS, sp.oo)
    out["V4_large_g"] = {"lim_c_eta_over_g2": str(sp.simplify(lead)),
                         "author_over_g2": str(-2 * (dS - 1) ** 2),
                         "ratio": str(sp.simplify(lead / (-2 * (dS - 1) ** 2))),
                         "exact_ratio_g8_a1": float(sp.N(c1.subs({gS: 8, dS: sp.Rational(3, 10)})
                                                        / (-2 * 8 ** 2 * (sp.Rational(3, 10) - 1) ** 2))),
                         "exact_ratio_g8_a2": float(sp.N(c_sym[2].subs({gS: 8, dS: sp.Rational(3, 10)})
                                                        / (-2 * 8 ** 2 * (sp.Rational(3, 10) - 1) ** 2)))}

    # ---- V3: twist after
    log("jets: V3 twist after")
    after = {}
    for a in (1, 2, 3):
        Ue, Uh, M, u, F = jets_densities(a, (1, 2), D, "after")
        ev = simp((M * XI * u + gS * u).norm() ** 2)
        Ue0 = simp(Ue.subs({psiS: 0, chiS: 0}) / (bS * kS) ** 2)
        Uh0 = simp(Uh.subs({psiS: 0, chiS: 0}) / (bS * kS) ** 2)
        Ue_chi = simp(Ue.subs(psiS, 0) / (bS * kS) ** 2)
        Uh_chi = simp(Uh.subs(psiS, 0) / (bS * kS) ** 2)
        after[f"a{a}"] = {
            "u_is_eigvec_resid": str(ev),
            "c_eta_chi0_psi0": str(Ue0), "c_h_chi0_psi0": str(Uh0),
            "c_eta_chi0_num": num(Ue0), "c_h_chi0_num": num(Uh0),
            "U_eta_chi_general_psi0": str(Ue_chi),
            "U_h_chi_general_psi0": str(Uh_chi),
            "U_eta_chi0.5": num(Ue, chi=0.5), "U_h_chi0.5": num(Uh, chi=0.5),
            "U_eta_chi0.5_d4_0": num(Ue, chi=0.5, delta4=0),
            "U_h_chi0.5_d4_0": num(Uh, chi=0.5, delta4=0),
            "U_eta_chi0.25": num(Ue, chi=0.25), "U_h_chi0.25": num(Uh, chi=0.25),
            "U_eta_chi0.5_psi_pi4": num(Ue, chi=0.5, psi=sp.pi / 4),
            "U_h_chi0.5_psi_pi4": num(Uh, chi=0.5, psi=sp.pi / 4),
        }
        log(f"  a={a}: chi0 c_eta={after[f'a{a}']['c_eta_chi0_num']:.4f} c_h={after[f'a{a}']['c_h_chi0_num']:.4f}; "
            f"chi0.5: U_eta={after[f'a{a}']['U_eta_chi0.5']:.3f} U_h={after[f'a{a}']['U_h_chi0.5']:.3f}")
    out["V3_after"] = after

    # ---- mutations
    log("jets: mutations")
    mut = {}
    # (m1) general vacuum D = diag(g, d1, d2, d3): the closed form
    d1S, d2S, d3S = sp.symbols("d1 d2 d3", real=True)
    Dg = sp.diag(gS, d1S, d2S, d3S)
    gen = {}
    for a in (1, 2, 3):
        Ue, Uh, M, u, F = jets_densities(a, (1, 2), Dg, "inside")
        ce = simp(Ue.subs(psiS, 0) / (bS * kS) ** 2)
        ch = simp(Uh.subs(psiS, 0) / (bS * kS) ** 2)
        gen[f"a{a}"] = {"c_eta": str(sp.factor(ce)), "c_h": str(sp.factor(ch))}
    mut["general_vacuum_inside_plane12"] = gen
    # (m2) twist plane (1,3), boost axes 1, 2, 3, twist inside
    p13 = {}
    for a in (1, 2, 3):
        Ue, Uh, M, u, F = jets_densities(a, (1, 3), Dg, "inside")
        ce = simp(Ue.subs(psiS, 0) / (bS * kS) ** 2)
        ch = simp(Uh.subs(psiS, 0) / (bS * kS) ** 2)
        dchi = simp(sp.diff(Ue, chiS))
        p13[f"a{a}"] = {"c_eta": str(sp.factor(ce)), "c_h": str(sp.factor(ch)),
                        "dU_eta_dchi": str(dchi),
                        "c_eta_num_g8_d0.3_d4_0.3": float(ce.subs({gS: 8, d1S: 1, d2S: sp.Rational(3, 10), d3S: sp.Rational(3, 10)})),
                        "c_eta_num_g8_d0.3_d4_0": float(ce.subs({gS: 8, d1S: 1, d2S: sp.Rational(3, 10), d3S: 0}))}
    mut["twist_plane_13_inside"] = p13
    # (m3) the gradient-axis labels: F_st = b k [M_chi, M_psi] does not depend on s, t
    #      beyond the weight eta^s eta^t; s = 0 (time gradient) flips the sign
    Ue_s1, _, _, _, _ = jets_densities(1, (1, 2), D, "inside", s_axis=1, t_axis=3)
    Ue_s2, _, _, _, _ = jets_densities(1, (1, 2), D, "inside", s_axis=2, t_axis=3)
    Ue_s0, _, _, _, _ = jets_densities(1, (1, 2), D, "inside", s_axis=0, t_axis=3)
    mut["gradient_axis"] = {"s1_minus_s2": str(simp(Ue_s1 - Ue_s2)),
                            "s0_plus_s1": str(simp(Ue_s0 + Ue_s1)),
                            "note": "F_st = b k [M_chi, M_psi]_eta carries no s, t beyond the weight eta^s eta^t; "
                                    "any spatial s != t gives the same density, s = 0 flips the sign"}
    # (m4) opposite vacuum sign (D00 = -g, the s = +1 convention)
    Ue, Uh, M, u, F = jets_densities(1, (1, 2), sp.diag(-gS, 1, dS, d4S), "inside")
    mut["D00_minus_g"] = {"c_eta_psi0": str(sp.factor(simp(Ue.subs(psiS, 0) / (bS * kS) ** 2)))}
    # (m5) boost axis 3 (commutes with R_12) at general psi and chi
    Ue, Uh, M, u, F = jets_densities(3, (1, 2), Dg, "inside")
    mut["boost3_general"] = {"U_eta": str(simp(Ue)), "U_h": str(simp(Uh))}
    out["mutations"] = mut

    # ---- cross-check of the registry instrument on the jets (numpy)
    log("jets: registry cross-check at one point")
    p = LAG.default_params(s=-1, g=8, delta=0.3)
    reg = {}
    for order in ("inside", "after"):
        Ue, Uh, M, u, F = jets_densities(1, (1, 2), D, order)
        sub = {gS: 8, dS: 0.3, d4S: 0.3, bS: 1, kS: 1, chiS: 0.5, psiS: 0.7}
        Mn = np.array(M.subs(sub).evalf(), dtype=float)[None]
        Mchi = np.array(sp.diff(M, chiS).subs(sub).evalf(), dtype=float)
        Mpsi = np.array(sp.diff(M, psiS).subs(sub).evalf(), dtype=float)
        A = np.zeros((4, 1, 4, 4))
        A[1, 0], A[3, 0] = Mchi, Mpsi
        de = float(4 * LAG.REGISTRY["I1"].density(A, Mn, p)[0])
        dh = float(4 * EXT.REGISTRY_EXT["I1_h"].density(A, Mn, p)[0])
        reg[order] = {"own_U_eta": float(Ue.subs(sub)), "registry_4xI1": de,
                      "own_U_h": float(Uh.subs(sub)), "registry_4xI1_h": dh}
        log(f"  {order}: own {reg[order]['own_U_eta']:.4f}/{reg[order]['own_U_h']:.4f} "
            f"registry {de:.4f}/{dh:.4f}")
    out["registry_crosscheck_chi0.5_psi0.7"] = reg
    save("jets", out)
    return out


# ====================================================================
# PART B: lattice, own fields
# ====================================================================
def radial_boost(X, Y, Z, chi, smooth=False, amp=0.5, w2=8.0):
    """1 + sinh(chi) K(n) + (cosh(chi) - 1) K2(n); chi(r) = amp exp(-r^2 / w2);
    smooth: chi(r) = amp (r / 2) exp((4 - r^2) / w2) (vanishes linearly at r = 0,
    peak amp at r = 2 for w2 = 8)."""
    r = np.sqrt(X * X + Y * Y + Z * Z)
    if chi is None:
        chi = amp * (r / 2.0) * np.exp((4.0 - r * r) / w2) if smooth else amp * np.exp(-r * r / w2)
    nx, ny, nz = X / r, Y / r, Z / r
    K = np.zeros(X.shape + (4, 4))
    for i, c in enumerate((nx, ny, nz)):
        K[..., 0, 1 + i] = K[..., 1 + i, 0] = c
    K2 = np.zeros_like(K)
    K2[..., 0, 0] = 1.0
    for i, c in enumerate((nx, ny, nz)):
        for j, d in enumerate((nx, ny, nz)):
            K2[..., 1 + i, 1 + j] = c * d
    return (np.eye(4)[None, None, None] + np.sinh(chi)[..., None, None] * K
            + (np.cosh(chi) - 1.0)[..., None, None] * K2), chi


def witness(cfg, D, k, order, amp=0.5, w2=8.0, smooth=False, Mcore=None):
    """M = L R D R^T L^T (inside) or R L D L^T R^T (after); D a 4-vector (uniform
    vacuum) or Mcore a per-cell field (the hedgehog).  Returns M and u = Q e0."""
    X, Y, Z = INS4.coords(cfg["n"], cfg["h"])
    L, chi = radial_boost(X, Y, Z, None, smooth=smooth, amp=amp, w2=w2)
    R = B8.rot_field(B8.G3, k * Z)
    Q = L @ R if order == "inside" else R @ L
    if Mcore is None:
        core = np.diag(D)[None, None, None]
        M = np.einsum("...ab,...bc,...dc->...ad", Q, np.broadcast_to(core, Q.shape), Q)
    else:
        M = np.einsum("...ab,...bc,...dc->...ad", Q, Mcore, Q)
    return INS4.sym4(M), Q[..., :, 0]


def own_timelike_u(M):
    """own eigen-solve: the eigenvector of M eta with negative eta-norm, unit."""
    lam, V = np.linalg.eig(M @ ETA)
    V = V.real
    n2 = np.einsum("...ak,a,...ak->...k", V, np.diag(ETA), V)
    k0 = np.argmin(n2, axis=-1)
    u = np.take_along_axis(V, k0[..., None, None], axis=-1)[..., 0]
    return u / np.sqrt(-np.take_along_axis(n2, k0[..., None], axis=-1))


def E_h_own(M, cfg, u):
    """4 h^3 sum_branches wt sum_cells sum_{i<j} tr(h F_ij h F_ij^T), h = eta + 2 (eta u)(eta u)^T."""
    hu = u @ ETA
    hc = ETA[None, None, None] + 2.0 * hu[..., :, None] * hu[..., None, :]
    tot = 0.0
    for br, wt in INS4.branches(cfg["stencil"]):
        A = [INS4.d1(M, ax, cfg["h"], br) for ax in range(3)]
        for i in range(3):
            for j in range(i + 1, 3):
                F = INS4.comm_eta(A[i], A[j])
                tot += wt * 4.0 * np.sum(np.einsum("...ab,...bc,...cd,...ad->...", hc, F, hc, F, optimize=True))
    return cfg["h"] ** 3 * tot


def E_h_reg(M, cfg):
    p = LAG.default_params(s=-1, g=8, delta=0.3)
    return 4.0 * LAG.term_lagrangian(EXT.REGISTRY_EXT["I1_h"], M, cfg, p)


def energies(M, cfg, u=None, with_reg=True):
    eu, ev = INS4.e_parts(M, cfg)
    if u is None:
        u = own_timelike_u(M)
    eh = E_h_own(M, cfg, u)
    rec = {"E_eta": float(eu), "V4": float(ev), "E_h_own": float(eh)}
    if with_reg:
        rec["E_h_reg"] = float(E_h_reg(M, cfg))
    return rec


def part_lattice():
    out = {}
    Ddeg = np.array([G, 1.0, DELTA, DELTA])
    Dcert = np.array([G, 1.0, DELTA, 0.0])

    def run(n, L, D, order, k, amp=0.5, w2=8.0, smooth=False, with_reg=True):
        cfg = INS4.base_cfg(s=-1, g=8, n=n, L=L, delta=0.3)
        M, u = witness(cfg, D, k, order, amp, w2, smooth)
        e = energies(M, cfg, u, with_reg)
        if k == 0.0 and n <= 64:
            e["E_h_own_eig_u"] = float(E_h_own(M, cfg, own_timelike_u(M)))
        e.update({"n": n, "L": L, "h": cfg["h"], "order": order, "k": k, "amp": amp, "w2": w2,
                  "smooth": smooth, "D": list(map(float, D))})
        log(f"  n{n} L{L:g} {order} k={k:+g} amp={amp} w2={w2} smooth={smooth} D3={D[3]:g}: "
            f"E_eta={e['E_eta']:.2f} V4={e['V4']:.2e} E_h_own={e['E_h_own']:.2f}"
            + (f" E_h_reg={e['E_h_reg']:.2f}" if with_reg else ""))
        return e

    # ---- B0: pure twist on the uniform vacuum
    log("lattice: B0 pure twist on the vacuum (no boost)")
    cfg = INS4.base_cfg(s=-1, g=8, n=32, L=48.0, delta=0.3)
    X, Y, Z = INS4.coords(32, cfg["h"])
    R = B8.rot_field(B8.G3, 1.0 * Z)
    M = INS4.sym4(np.einsum("...ab,bc,...dc->...ad", R, np.diag(Ddeg), R))
    e0 = energies(M, cfg, R[..., :, 0])
    out["B0_pure_twist_vacuum_k1"] = e0
    log(f"  pure twist: {e0}")

    # ---- B1: n64 L48 degenerate, inside + after, and the certified vacuum
    log("lattice: B1 n64 L48")
    b1 = {"baseline_deg": run(64, 48.0, Ddeg, "inside", 0.0),
          "baseline_cert": run(64, 48.0, Dcert, "inside", 0.0)}
    for D, tag in ((Ddeg, "deg"), (Dcert, "cert")):
        for order in ("inside", "after"):
            for k in (0.5, 1.0, 2.0):
                e = run(64, 48.0, D, order, k)
                base = b1[f"baseline_{tag}"]
                e["dE_eta"] = e["E_eta"] - base["E_eta"]
                e["dE_h_own"] = e["E_h_own"] - base["E_h_own"]
                e["dE_h_reg"] = e["E_h_reg"] - base["E_h_reg"]
                b1[f"{tag}_{order}_k{k:g}"] = e
                log(f"    dE_eta={e['dE_eta']:+.2f} dE_h_own={e['dE_h_own']:+.2f} dE_h_reg={e['dE_h_reg']:+.2f}")
                save("lattice", {"B1": b1})
    out["B1"] = b1

    # ---- B2: spacing 1.5, n32 L48 and n48 L72 (V6)
    log("lattice: B2 spacing 1.5")
    b2 = {}
    for n, L in ((32, 48.0), (48, 72.0)):
        base = run(n, L, Ddeg, "inside", 0.0)
        b2[f"n{n}_baseline"] = base
        for k in (0.5, 1.0, 2.0):
            e = run(n, L, Ddeg, "inside", k)
            e["dE_eta"] = e["E_eta"] - base["E_eta"]
            e["dE_h_own"] = e["E_h_own"] - base["E_h_own"]
            e["dE_h_reg"] = e["E_h_reg"] - base["E_h_reg"]
            b2[f"n{n}_inside_k{k:g}"] = e
            log(f"    dE_eta={e['dE_eta']:+.2f} dE_h_own={e['dE_h_own']:+.2f}")
    out["B2"] = b2
    save("lattice", out)

    # ---- B3: mutations on n64 (amplitude, width, smooth core, sign of k, small k)
    log("lattice: B3 mutations")
    b3 = {}
    base = b1["baseline_deg"]
    for amp in (0.25, 1.0):
        bb = run(64, 48.0, Ddeg, "inside", 0.0, amp=amp, with_reg=False)
        e = run(64, 48.0, Ddeg, "inside", 1.0, amp=amp, with_reg=False)
        b3[f"amp{amp:g}_k1"] = {"dE_eta": e["E_eta"] - bb["E_eta"], "dE_h_own": e["E_h_own"] - bb["E_h_own"],
                                "E_eta_base": bb["E_eta"], "E_h_base": bb["E_h_own"]}
        log(f"    amp {amp}: dE_eta={b3[f'amp{amp:g}_k1']['dE_eta']:+.2f} dE_h={b3[f'amp{amp:g}_k1']['dE_h_own']:+.2f}")
    bb = run(64, 48.0, Ddeg, "inside", 0.0, w2=32.0, with_reg=False)
    e = run(64, 48.0, Ddeg, "inside", 1.0, w2=32.0, with_reg=False)
    b3["width4_k1"] = {"dE_eta": e["E_eta"] - bb["E_eta"], "dE_h_own": e["E_h_own"] - bb["E_h_own"],
                       "E_eta_base": bb["E_eta"], "E_h_base": bb["E_h_own"]}
    log(f"    width 4: {b3['width4_k1']}")
    bb = run(64, 48.0, Ddeg, "inside", 0.0, smooth=True, with_reg=False)
    for order in ("inside", "after"):
        e = run(64, 48.0, Ddeg, order, 1.0, smooth=True, with_reg=False)
        b3[f"smooth_{order}_k1"] = {"dE_eta": e["E_eta"] - bb["E_eta"], "dE_h_own": e["E_h_own"] - bb["E_h_own"],
                                    "E_eta_base": bb["E_eta"], "E_h_base": bb["E_h_own"]}
        log(f"    smooth {order}: {b3[f'smooth_{order}_k1']}")
    e = run(64, 48.0, Ddeg, "inside", -1.0, with_reg=False)
    b3["kminus1_inside"] = {"dE_eta": e["E_eta"] - base["E_eta"], "dE_h_own": e["E_h_own"] - base["E_h_own"]}
    for k in (0.125, 0.25):
        e = run(64, 48.0, Ddeg, "inside", k, with_reg=False)
        b3[f"smallk_{k:g}_inside"] = {"dE_eta": e["E_eta"] - base["E_eta"], "dE_h_own": e["E_h_own"] - base["E_h_own"],
                                      "dE_eta_over_k2": (e["E_eta"] - base["E_eta"]) / k ** 2,
                                      "dE_h_over_k2": (e["E_h_own"] - base["E_h_own"]) / k ** 2}
        log(f"    k={k}: {b3[f'smallk_{k:g}_inside']}")
    # finer grid, h = 0.5 (n96 L48), k = 0.5 and 1
    bb = run(96, 48.0, Ddeg, "inside", 0.0, with_reg=False)
    b3["n96_baseline"] = bb
    for k in (0.5, 1.0):
        e = run(96, 48.0, Ddeg, "inside", k, with_reg=False)
        b3[f"n96_inside_k{k:g}"] = {"dE_eta": e["E_eta"] - bb["E_eta"], "dE_h_own": e["E_h_own"] - bb["E_h_own"]}
        log(f"    n96 k={k}: {b3[f'n96_inside_k{k:g}']}")
    out["B3"] = b3
    save("lattice", out)
    return out


# ====================================================================
# PART C: the relaxed hedgehog (V7)
# ====================================================================
def part_hedgehog():
    R13W = _load("m5_32_r13w_common", "m5_32_r13w_common.py")
    Mhh, cfg, src = R13W.seed_hedgehog(32, 48)
    log(f"hedgehog seed: {src}, n={cfg['n']} h={cfg['h']}")
    out = {"source": src, "n": cfg["n"], "L": cfg["L"]}
    X, Y, Z = INS4.coords(cfg["n"], cfg["h"])
    Lb, chi = radial_boost(X, Y, Z, None)
    e_hh = energies(Mhh, cfg)
    Md = INS4.sym4(np.einsum("...ab,...bc,...dc->...ad", Lb, Mhh, Lb))
    e_d = energies(Md, cfg)
    out["E_hedgehog"], out["E_dressed"] = e_hh, e_d
    log(f"  hedgehog: {e_hh}\n  dressed: {e_d}")
    rows = {}
    for k in (0.5, 1.0, 2.0):
        R = B8.rot_field(B8.G3, k * Z)
        Mt = INS4.sym4(np.einsum("...ab,...bc,...dc->...ad", R, Mhh, R))          # twist alone
        Mi = INS4.sym4(np.einsum("...ab,...bc,...dc->...ad", Lb @ R, Mhh, Lb @ R))  # twist inside
        Ma = INS4.sym4(np.einsum("...ab,...bc,...dc->...ad", R @ Lb, Mhh, R @ Lb))  # twist after
        et, ei, ea = energies(Mt, cfg), energies(Mi, cfg), energies(Ma, cfg)
        r = {"E_twist_alone_eta": et["E_eta"] - e_hh["E_eta"],
             "E_twist_alone_h": et["E_h_own"] - e_hh["E_h_own"],
             "dE_eta_inside": ei["E_eta"] - e_d["E_eta"],
             "dE_h_inside": ei["E_h_own"] - e_d["E_h_own"],
             "dE_h_reg_inside": ei["E_h_reg"] - e_d["E_h_reg"],
             "dE_eta_after": ea["E_eta"] - e_d["E_eta"],
             "dE_h_after": ea["E_h_own"] - e_d["E_h_own"],
             "V4_twist": et["V4"], "V4_inside": ei["V4"]}
        r["cross_eta_inside"] = r["dE_eta_inside"] - r["E_twist_alone_eta"]
        r["cross_h_inside"] = r["dE_h_inside"] - r["E_twist_alone_h"]
        r["cross_eta_after"] = r["dE_eta_after"] - r["E_twist_alone_eta"]
        r["cross_h_after"] = r["dE_h_after"] - r["E_twist_alone_h"]
        rows[f"k{k:g}"] = r
        log(f"  k={k}: " + " ".join(f"{a}={b:+.1f}" for a, b in r.items()))
    out["rows"] = rows
    # mutation: the hedgehog dressed with the smooth-core boost
    Ls, _ = radial_boost(X, Y, Z, None, smooth=True)
    Mds = INS4.sym4(np.einsum("...ab,...bc,...dc->...ad", Ls, Mhh, Ls))
    e_ds = energies(Mds, cfg, with_reg=False)
    R = B8.rot_field(B8.G3, 1.0 * Z)
    Mis = INS4.sym4(np.einsum("...ab,...bc,...dc->...ad", Ls @ R, Mhh, Ls @ R))
    e_is = energies(Mis, cfg, with_reg=False)
    out["smooth_core_k1"] = {"dE_eta_inside": e_is["E_eta"] - e_ds["E_eta"],
                             "dE_h_inside": e_is["E_h_own"] - e_ds["E_h_own"],
                             "cross_eta": e_is["E_eta"] - e_ds["E_eta"] - rows["k1"]["E_twist_alone_eta"],
                             "cross_h": e_is["E_h_own"] - e_ds["E_h_own"] - rows["k1"]["E_twist_alone_h"]}
    log(f"  smooth core k1: {out['smooth_core_k1']}")
    save("hedgehog", out)
    return out


if __name__ == "__main__":
    mode = ARGV[1] if len(ARGV) > 1 else "all"
    if mode in ("jets", "all"):
        part_jets()
    if mode in ("lattice", "all"):
        part_lattice()
    if mode in ("hedgehog", "all"):
        part_hedgehog()
    log("done")
