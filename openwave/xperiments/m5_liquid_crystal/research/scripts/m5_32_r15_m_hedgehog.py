"""M5.32 R15-M: admissibility of the author's degenerate-vacuum object L_P on the lattice.

Equations first.  The field M(x) (4x4 symmetric), N = M eta, eta = diag(-1, 1, 1, 1).  The
degenerate vacuum d = diag(g, 1, delta, delta), g = 8, delta = 0.3, N-spectrum (-g, 1, delta, delta).

    L_P (our reading of the author's 09-05 object, E-orientation, static):
        E_stat = E_u + V4^dd + mu * SPLIT + c_P * K_P^23
        E_u    = 4 h^3 sum_cells sum_{i<j} <F_ij, F_ij>_eta,  F_ij = [A_i, A_j]_eta  (certified -4 I1)
        V4^dd  = W1 h^3 sum_cells sum_{p=1..4} (tr N^p - C_p)^2,  C_p = (-g)^p + 1 + 2 delta^p
        SPLIT  = h^3 sum_cells (lambda_2 - lambda_3)^2 = h^3 sum_cells (s^2 - 4 p)
        K_P^23 = (1/2) h^3 sum_cells sum_i tr(Om_i^T eta Om_i eta),  Om_i = P23 A_i eta P23,
                 P23 = I - P_g - P_1 (the exact Lagrange projector onto the (2,3) eigenplane)
    (definitions, projector formulas and the exact gradients: m5_32_r15_common.py, selftest 19/19)

M-a (the Hessian at the vacuum).  H_ij = d^2 (V4^dd + mu SPLIT) / dM_i dM_j on the 10-dim
symmetric-matrix space at M = d, per unit volume, by central differences of the ANALYTIC
gradients (eps 1e-5).  Prediction stated before the run: at mu = 0 the null space is
7-dimensional (5 Lorentz-orbit directions: 3 boosts + the two rotations mixing axis 1 with
the (2,3) plane; the (2,3) rotation is a stabilizer; + the 2-dim traceless (2,3) block, which
V4 sees only at quartic order since it pins the spectrum as a set); for mu > 0 it is
5-dimensional (the orbit alone).  A different count is a finding.

M-b (the relaxed hedgehog).  Seed: the uniaxial radial hedgehog B8.dressed(cfg, 0) on the
degenerate vacuum (the orientation texture Q_h d Q_h^T, no boost).  FIRE descent of E_stat with
the vacuum shell pinned (depth 1.6), dt0 0.01, dt_max 0.1, max_iter as given, on (n, L) in
{(32, 48), (48, 72)} x mu in {0, 1e-2} x c_P in {0, 1}.  Reads on the end field: the energy
ladder, the sorted N-spectrum along the +x axis through the center, the split
(lambda_3 - lambda_2) along that axis, the static energy density in radial shells, the
exterior spectrum drift (shell 0.35 L < r < 0.45 L, max |sorted spectrum - (-g, delta, delta, 1)|),
the (2,3)-clock inertias of I1 and K_P^23 with a0 = a0_local (J M - M J about the leading
spatial eigenvector).

Pre-registered verdict (fixed before any run): ADMISSIBLE if the descent is finite AND the
energy inside r < L/4 holds >= 0.8 of the total static energy AND the exterior drift < 0.05;
NOT_LOCALIZED if the descent is finite but either localization test fails (energy migrating to
the pin shell, or the exterior leaving the degenerate vacuum); RUNAWAY if non-finite.  The
certified L_cert hedgehog (R10 n32 L48 seed, the same reads) is the reference row.

usage:  python3 m5_32_r15_m_hedgehog.py hess
        python3 m5_32_r15_m_hedgehog.py relax <n> <L> <mu> <cP> <maxit>
        python3 m5_32_r15_m_hedgehog.py ref
        python3 m5_32_r15_m_hedgehog.py collect
"""
import sys
ARGS = list(sys.argv[1:])            # captured before any import (the R7 argv-wipe lesson)
import os, json, time
import numpy as np
import m5_32_r15_common as C15

INS4, C13, B8 = C15.INS4, C15.C13, C15.B8
G, DELTA = C15.G, C15.DELTA
OUT = os.path.join(C15.CK, "m_hedgehog")
os.makedirs(OUT, exist_ok=True)
log = C15.log


# ------------------------------------------------ densities
def static_density(M, cfg):
    """per-cell h^3-weighted static E-density of L_P: E_u + V4^dd + mu SPLIT + c_P K_P^23."""
    h3 = cfg["h"] ** 3
    du = np.zeros(M.shape[:3])
    for br, (A, wt) in INS4.a_fields(M, cfg).items():
        for i in range(3):
            for j in range(i + 1, 3):
                F = INS4.comm_eta(A[i], A[j])
                du += wt * 4.0 * INS4.inner_eta(F, F)
    Me = M @ C15.ETA
    P, t = Me, []
    for p in range(4):
        if p:
            P = P @ Me
        t.append(np.einsum("...kk->...", P))
    cp = INS4.c4_of(cfg)
    dv = C15.W1 * sum((t[p] - cp[p]) ** 2 for p in range(4)) if hasattr(C15, "W1") else None
    if dv is None:
        dv = INS4.W1 * sum((t[p] - cp[p]) ** 2 for p in range(4))
    ds = C15.split_cells(M, need_grad=False)[0]
    dk = C15.kp23_static_density(M, cfg)
    return {"E_u": h3 * du, "V4": h3 * dv, "split": cfg["mu"] * h3 * ds, "KP": cfg["cP"] * dk}


def reads(M, cfg, tag=""):
    n, h, L = cfg["n"], cfg["h"], cfg["L"]
    X, Y, Z = INS4.coords(n, h)
    r = np.sqrt(X * X + Y * Y + Z * Z)
    dens = static_density(M, cfg)
    tot = sum(dens.values())
    T = float(np.sum(tot))
    fr = {f"frac_r_lt_L/{k}": float(np.sum(tot[r < L / k]) / T) if T != 0 else None for k in (4, 8)}
    edges = np.linspace(0.0, L / 2, 13)
    shells = [float(np.sum(tot[(r >= a) & (r < b)])) for a, b in zip(edges[:-1], edges[1:])]
    lam = C15.sorted_spectrum(M)
    c = n // 2
    axis_r = X[c:, c, c]
    axis_lam = lam[c:, c, c, :]
    split_axis = axis_lam[:, 2] - axis_lam[:, 1]
    ext = (r > 0.35 * L) & (r < 0.45 * L)
    drift = float(np.max(np.abs(lam[ext] - np.array([-G, DELTA, DELTA, 1.0])[None, :])))
    a0 = C13.a0_local(M)
    kin_i1 = float(INS4.kin_of(M, a0, cfg))
    kin_kp = float(C15.kp23_energy_grad(M, cfg, a0, need_grad=False)[2])
    pp = C15.lp_parts(M, cfg)
    pp["I1h_static"] = C15.i1h_static(M, cfg)
    return {"tag": tag, "parts": pp, "energy_fractions": fr, "shell_energy": shells, "shell_edges": edges.tolist(),
            "axis_r": axis_r.tolist(), "axis_spectrum": axis_lam.tolist(), "axis_split": split_axis.tolist(),
            "exterior_drift": drift, "kin_I1_a0local": kin_i1, "kin_KP23_raw_a0local": kin_kp,
            "max_split": float(np.max(lam[..., 2] - lam[..., 1])),
            "density_max": {k: float(np.max(v)) for k, v in dens.items()}}


def verdict(rd, stop):
    if stop == "non-finite":
        return "RUNAWAY"
    ok = rd["energy_fractions"]["frac_r_lt_L/4"] is not None and rd["energy_fractions"]["frac_r_lt_L/4"] >= 0.8 \
        and rd["exterior_drift"] < 0.05
    return "ADMISSIBLE" if ok else "NOT_LOCALIZED"


# ------------------------------------------------ M-a
def sym_basis():
    B = []
    for a in range(4):
        for b in range(a, 4):
            X = np.zeros((4, 4)); X[a, b] = X[b, a] = 1.0
            B.append(X / np.sqrt(np.sum(X * X)))
    return B


def labeled_dirs(d):
    """the tangent directions of the Lorentz orbit dM = K d + d K^T (K eta-antisymmetric) and the
    (2,3) traceless block, at the vacuum d, as unit symmetric matrices (zero maps dropped)."""
    out = {}
    def gen(i, j, boost):
        K = np.zeros((4, 4))
        if boost:
            K[0, j] = K[j, 0] = 1.0
        else:
            K[i, j], K[j, i] = 1.0, -1.0
        return K
    for j in (1, 2, 3):
        out[f"boost_{j}"] = gen(0, j, True)
    for (i, j) in ((1, 2), (1, 3), (2, 3)):
        out[f"rot_{i}{j}"] = gen(i, j, False)
    dirs = {}
    for k, K in out.items():
        dM = K @ d + d @ K.T
        nrm = np.sqrt(np.sum(dM * dM))
        if nrm > 1e-12:
            dirs[k] = dM / nrm
    X = np.zeros((4, 4)); X[2, 2], X[3, 3] = 1.0, -1.0
    dirs["split_diag_23"] = X / np.sqrt(2.0)
    X = np.zeros((4, 4)); X[2, 3] = X[3, 2] = 1.0
    dirs["split_offdiag_23"] = X / np.sqrt(2.0)
    return dirs


def hess_mode():
    n, L = 4, 6.0
    res = {}
    for mu in (0.0, 1e-3, 1e-2, 1e-1):
        cfg = C15.cfg_dd(n, L, mu=mu, cP=0.0)
        d = INS4.vac4(cfg)
        M0 = np.broadcast_to(d, (n, n, n, 4, 4)).copy()
        B = sym_basis()
        vol = n ** 3 * cfg["h"] ** 3
        def grad_pot(M):
            gu = INS4.grad(M, cfg)                  # E_u + V4 gradient (E_u part = 0 on a uniform field)
            gs = C15.split_energy_grad(M, cfg)[1]
            return gu + gs
        eps = 1e-5
        H = np.zeros((10, 10))
        for j, Xj in enumerate(B):
            gp = grad_pot(M0 + eps * Xj[None, None, None])
            gm = grad_pot(M0 - eps * Xj[None, None, None])
            col = (gp - gm) / (2 * eps)
            for i, Xi in enumerate(B):
                H[i, j] = float(np.sum(col * Xi[None, None, None])) / vol
        H = 0.5 * (H + H.T)
        w, V = np.linalg.eigh(H)
        thr = 1e-6                                   # absolute (per unit volume): the FD floor is 1e-8, the smallest true stiffness 2e-3
        null = int(np.sum(np.abs(w) < thr))
        dirs = labeled_dirs(d)
        Vn = V[:, np.abs(w) < thr]
        overlaps = {}
        for k, X in dirs.items():
            c = np.array([np.sum(X * Bi) for Bi in B])
            overlaps[k] = float(np.sqrt(np.sum((Vn.T @ c) ** 2))) if Vn.shape[1] else 0.0
        res[f"mu_{mu:g}"] = {"eigenvalues": w.tolist(), "null_count": null, "threshold": float(thr),
                             "null_overlap_of_labeled_dirs": overlaps}
        log(f"M-a mu {mu:g}: null count {null}, eigenvalues {np.array2string(w, precision=3)}")
        log(f"      overlaps of the labeled directions with the null space: " + ", ".join(f"{k} {v:.3f}" for k, v in overlaps.items()))
    json.dump(res, open(os.path.join(OUT, "hess.json"), "w"), indent=1)
    return res


# ------------------------------------------------ M-b
def relax_mode(n, L, mu, cP, maxit):
    tag = f"relax_n{n}_L{L:g}_mu{mu:g}_cP{cP:g}"
    cfg = C15.cfg_dd(n, L, mu=mu, cP=cP)
    M0 = C15.seed_uniaxial(cfg)
    free = ~INS4.pin_shell(n, cfg["h"])
    r0 = reads(M0, cfg, "seed")
    log(f"{tag} seed: {json.dumps(r0['parts'])} frac {r0['energy_fractions']} drift {r0['exterior_drift']:.3e}")
    ck = os.path.join(OUT, tag + ".npy")
    M, info = C15.fire_lp(M0, cfg, free, maxit, log_every=100, tag=tag, ck_path=ck, ck_every=500)
    np.save(ck, M)
    r1 = reads(M, cfg, "end")
    v = verdict(r1, info["stop"])
    rec = {"n": n, "L": L, "h": cfg["h"], "mu": mu, "cP": cP, "maxit": maxit, "stop": info["stop"], "iters": info["iters"],
           "wall_s": info["wall_s"], "seed": r0, "end": r1, "verdict": v, "trace": info["trace"], "field": ck}
    json.dump(rec, open(os.path.join(OUT, tag + ".json"), "w"), indent=1)
    log(f"{tag} END {v}: stop {info['stop']} it {info['iters']} E_stat {r1['parts']['E_stat']:.6f} "
        f"frac(L/4) {r1['energy_fractions']['frac_r_lt_L/4']:.3f} drift {r1['exterior_drift']:.3e} "
        f"max split {r1['max_split']:.4f} kin_I1 {r1['kin_I1_a0local']:.4f} kin_KP {r1['kin_KP23_raw_a0local']:.4f}")
    return rec


def ref_mode():
    M, cfg, src = C13.seed_hedgehog(32, 48)
    cfg = dict(cfg); cfg["mu"], cfg["cP"] = 0.0, 0.0
    # the certified reads: on the certified vacuum the split is the certified (delta, 0) gap, not a defect
    r = reads(M, cfg, "L_cert n32 L48 hedgehog (R10)")
    rec = {"source": src, "reads": r}
    json.dump(rec, open(os.path.join(OUT, "ref_cert.json"), "w"), indent=1)
    log(f"ref: frac {r['energy_fractions']} drift(vs degenerate targets, expected ~delta) {r['exterior_drift']:.3f} E_stat {r['parts']['E_stat']:.4f}")
    return rec


def collect_mode():
    import glob
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    recs = []
    for f in sorted(glob.glob(os.path.join(OUT, "relax_*.json"))):
        rec = json.load(open(f))
        recs.append(rec)
    hess = json.load(open(os.path.join(OUT, "hess.json"))) if os.path.exists(os.path.join(OUT, "hess.json")) else None
    ref = json.load(open(os.path.join(OUT, "ref_cert.json"))) if os.path.exists(os.path.join(OUT, "ref_cert.json")) else None
    def outer_slope(e, L):
        ed = np.array(e["shell_edges"]); sh = np.array(e["shell_energy"])
        rc = 0.5 * (ed[1:] + ed[:-1])
        sel = (rc > L / 4) & (sh > 0)
        return float(np.polyfit(np.log(rc[sel]), np.log(sh[sel]), 1)[0]) if np.sum(sel) >= 3 else None

    def calibrated(e, L, stop):
        """the rule calibrated on the certified reference BEFORE any M-b end state existed (progress log
        2026-09-05 20:36 UTC): LOCALIZED iff the log-log slope of the shell energy over the outer shells
        (r > L/4) is <= 0 and the exterior drift < 0.05."""
        if stop == "non-finite":
            return "RUNAWAY", None
        sl = outer_slope(e, L)
        ok = sl is not None and sl <= 0.0 and e["exterior_drift"] < 0.05
        return ("ADMISSIBLE" if ok else "NOT_LOCALIZED"), sl

    rows = []
    for rec in recs:
        e = rec["end"]
        vcal, slope = calibrated(e, rec["L"], rec["stop"])
        rows.append({"verdict_calibrated": vcal, "outer_shell_slope": slope, "verdict_raw_0.8_rule": rec["verdict"],"n": rec["n"], "L": rec["L"], "mu": rec["mu"], "cP": rec["cP"], "stop": rec["stop"], "iters": rec["iters"],
                     "E_stat_seed": rec["seed"]["parts"]["E_stat"], "E_stat_end": e["parts"]["E_stat"],
                     "E_u": e["parts"]["E_u"], "V4": e["parts"]["V4"], "split": e["parts"]["split"], "KP": e["parts"]["KP"],
                     "I1h_static": e["parts"]["I1h_static"], "frac_L4": e["energy_fractions"]["frac_r_lt_L/4"],
                     "frac_L8": e["energy_fractions"]["frac_r_lt_L/8"], "drift": e["exterior_drift"], "max_split": e["max_split"],
                     "kin_I1": e["kin_I1_a0local"], "kin_KP": e["kin_KP23_raw_a0local"], "wall_s": rec["wall_s"]})
    if ref is not None:
        ref["outer_shell_slope"] = outer_slope(ref["reads"], 48.0)
    out = {"rung": "R15-M", "hessian": hess, "reference_cert": ref, "runs": rows,
           "verdict_rule_raw": "pre-registered in the docstring: finite, frac(r < L/4) >= 0.8, drift < 0.05 (calls the certified reference itself NOT_LOCALIZED: frac 0.62)",
           "verdict_rule_calibrated": "LOCALIZED iff the log-log slope of the shell energy over the outer shells (r > L/4) <= 0 and drift < 0.05 (fixed 2026-09-05 20:36 UTC before any M-b end state)",
           "fields": [os.path.relpath(r["field"], C15.RES) for r in recs]}
    json.dump(out, open(os.path.join(C15.DATA, "m5_32_r15_m_hedgehog.json"), "w"), indent=1)
    if recs:
        fig, axes = plt.subplots(2, 2, figsize=(12, 9))
        for rec in recs:
            e = rec["end"]
            lab = f"n{rec['n']} L{rec['L']:g} mu {rec['mu']:g} cP {rec['cP']:g}"
            ar = np.array(e["axis_r"]); al = np.array(e["axis_spectrum"])
            axes[0, 0].plot(ar, al[:, 1], label=lab); axes[0, 0].plot(ar, al[:, 2], ls="--", color=axes[0, 0].lines[-1].get_color())
            axes[0, 1].plot(ar, e["axis_split"], label=lab)
            ed = np.array(e["shell_edges"]); sh = np.array(e["shell_energy"])
            axes[1, 0].semilogy(0.5 * (ed[1:] + ed[:-1]), np.maximum(sh, 1e-12), marker="o", label=lab)
            tr = rec["trace"]
            axes[1, 1].plot([t["it"] for t in tr], [t["E_stat"] for t in tr], label=lab)
        axes[0, 0].set_title("pair eigenvalues lambda_2 (solid), lambda_3 (dashed) along +x"); axes[0, 0].set_xlabel("r")
        axes[0, 1].set_title("split lambda_3 - lambda_2 along +x"); axes[0, 1].set_xlabel("r")
        axes[1, 0].set_title("static energy per radial shell"); axes[1, 0].set_xlabel("r")
        axes[1, 1].set_title("E_stat descent"); axes[1, 1].set_xlabel("iteration")
        for ax in axes.flat:
            ax.grid(alpha=0.3)
        axes[0, 1].legend(fontsize=7)
        fig.suptitle("M5.32 R15-M: the uniaxial hedgehog relaxed on the degenerate vacuum under L_P")
        fig.tight_layout()
        fig.savefig(os.path.join(C15.PLOTS, "m5_32_r15_m_hedgehog.png"), dpi=110)
    log(f"collected {len(rows)} runs -> data/m5_32_r15_m_hedgehog.json")
    for r in rows:
        log(f"  n{r['n']} L{r['L']:g} mu {r['mu']:g} cP {r['cP']:g}: {r['verdict_calibrated']} (raw rule {r['verdict_raw_0.8_rule']}) E {r['E_stat_end']:.4f} E_u {r['E_u']:.4f} KP {r['KP']:.4f} frac {r['frac_L4']:.3f} slope {r['outer_shell_slope']:.2f} drift {r['drift']:.3e} split {r['max_split']:.4f} kin_I1 {r['kin_I1']:.2e} kin_KP {r['kin_KP']:.2e} stop {r['stop']}")
    if ref is not None:
        log(f"  reference (certified L_cert hedgehog n32 L48): frac {ref['reads']['energy_fractions']['frac_r_lt_L/4']:.3f} outer slope {ref['outer_shell_slope']:.2f}")
    return out


if __name__ == "__main__":
    mode = ARGS[0]
    if mode == "hess":
        hess_mode()
    elif mode == "relax":
        relax_mode(int(ARGS[1]), float(ARGS[2]), float(ARGS[3]), float(ARGS[4]), int(ARGS[5]))
    elif mode == "ref":
        ref_mode()
    elif mode == "collect":
        collect_mode()
