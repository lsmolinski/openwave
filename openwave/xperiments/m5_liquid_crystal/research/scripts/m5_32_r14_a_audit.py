"""M5.32 R14-A adversarial audit: the coexistence conjecture as a linear program.

Independent re-assembly of the R14-A LP from the measured rows JSON
(data/m5_32_r14_a_rows.json) with the auditor's own sign conventions, tolerances,
solver calls (HiGHS dual simplex AND interior point), cutting planes on the exact
UV quadratic forms, an own Farkas certificate verified in exact rational arithmetic,
row-group dropping, mutations, and own re-measurements (plane-wave Q forms of T1 and
K_P_h with an own slab construction; tail plateaus of K_P_h, T1, R_etaMeta on the
n48 L72 hedgehog with own stencil and own shell binning).

The producer's LP script and its two result JSONs were never opened.

Usage: python3 m5_32_r14_a_audit.py [all | lp | uvq | tail]
Writes data/m5_32_r14_a_audit.json.
"""

import sys
ARGV = list(sys.argv)                     # captured BEFORE any import (kt_form wipes it)

import importlib.util
import json
import os
import time
from fractions import Fraction

import numpy as np
from scipy.optimize import linprog

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.dirname(HERE)
DATA = os.path.join(RES, "data")
CK10 = os.path.join(RES, "checkpoints", "m5_32_r10")
CK13 = os.path.join(RES, "checkpoints", "m5_32_r13w")
OUT_JSON = os.path.join(DATA, "m5_32_r14_a_audit.json")
ROWS_JSON = os.path.join(DATA, "m5_32_r14_a_rows.json")

ETA = np.diag([-1.0, 1.0, 1.0, 1.0])
ETA_D = np.array([-1.0, 1.0, 1.0, 1.0])
T0 = time.time()


def log(m):
    print(f"[{time.time() - T0:7.1f}s] {m}", flush=True)


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


INS4 = _load("m5_21_3_a_4d", "m5_21_3_a_4d.py")            # certified stack (control only)
R13W = _load("m5_32_r13w_common", "m5_32_r13w_common.py")   # cfg_of
R14T = _load("m5_32_r14_terms", "m5_32_r14_terms.py")       # cross-check of own entrants only
R8 = _load("m5_32_r8_a_quartics", "m5_32_r8_a_quartics.py")  # quartic densities for own UV rows
sys.argv = ARGV

ROWS = json.load(open(ROWS_JSON))
NAMES = [b["name"] for b in ROWS["basis"]]
FIXED = {"I1": 4.0}
TWO_DERIV = ["K_T", "K_lambda", "R_etaMeta", "R_hcov", "K_P_h", "T1", "T2", "T3", "T4"]
QUARTICS = ["Q_I1sq", "C6a", "C6b"]
BPLUS = ["I2", "I3", "I4", "I5", "I6", "I1_h", "J1", "J2", "Pgrad"] + TWO_DERIV + QUARTICS
PARITY_ODD = ["E1", "E2", "E3"]
DIRS = ["1|0|0", "0|1|0", "0|0|1", "1|1|0", "1|0|1", "0|1|1", "1|1|1"]
RESULT = {}


# ============================================================ row assembly (own)
class Row:
    """sum_k a[k] x_k + aK * K <= rhs, with x over ALL basis names (I1 folded later)."""

    def __init__(self, label, group, a, rhs, aK=0.0):
        self.label, self.group, self.a, self.rhs, self.aK = label, group, dict(a), float(rhs), float(aK)


def vec_of(d):
    return {k: float(d.get(k, 0.0)) for k in NAMES}


def tail_plateau(field_key, lo, hi):
    """a_k = mean of the shell energy per unit r over the shells whose CENTERS lie in [lo, hi]."""
    T = ROWS[field_key]
    h = T["_shell_edges_h"]
    n_sh = len(T["I1"]["shell_per_r"])
    centers = (np.arange(n_sh) + 0.5) * h
    idx = [i for i in range(n_sh) if lo <= centers[i] <= hi]
    return {k: float(np.mean([T[k]["shell_per_r"][i] for i in idx])) for k in NAMES}, idx


def sheet_S(kind, w, mask_ratio=2.5, fine="nz256", coarse="nz128"):
    """S_k(w) at the fine resolution; entries whose coarse/fine ratio exceeds mask_ratio are 0."""
    Sf = ROWS["sheets"][f"{kind}|{fine}|w{w:g}"]["S_per_area"]
    Sc = ROWS["sheets"][f"{kind}|{coarse}|w{w:g}"]["S_per_area"]
    out, masked = {}, []
    for k in NAMES:
        f, c = float(Sf[k]), float(Sc[k])
        if mask_ratio is not None and abs(f) > 0 and abs(c) / abs(f) > mask_ratio:
            out[k] = 0.0
            masked.append(k)
        elif mask_ratio is not None and f == 0.0 and abs(c) > 0:
            out[k] = 0.0
            masked.append(k)
        else:
            out[k] = f
    return out, masked


def own_quartic_rows(n_rows=150, seed=20260905):
    if isinstance(seed, int) and seed == 2:
        seed = 777
    """UV quartic rows on random symmetric SPATIAL jets (A_0 = 0): the exact degree-4 forms."""
    rng = np.random.default_rng(seed)
    rows = []
    for i in range(n_rows):
        A = np.zeros((4, 4, 4))
        for mu in (1, 2, 3):
            X = rng.standard_normal((4, 4))
            A[mu] = 0.5 * (X + X.T)
        A4 = A[:, None, :, :]                         # (4, 1, 4, 4): one cell
        i1 = float(R8.d_I1(A4)[0])
        a = {"Q_I1sq": i1 * i1, "C6a": float(R8.d_C6a(A4)[0]), "C6b": float(R8.d_C6b(A4)[0])}
        rows.append(Row(f"uv4_{i}", "uv_quartic", {k: -v for k, v in a.items()}, 0.0))
    return rows


def assemble(opts):
    """opts: tail_tol (3.0), tail_fields ('both'|'n48'|'n32'), sheet_mask (2.5|None),
    sheet_tol (1e-3), pair_tol (1e-3), relaxed_eps ('tol'|0), coulomb (bool), channels (bool),
    uv (bool), sheets ('both'|'twist'|'zigzag'|None), pairs (bool), quartic_rows (bool),
    flip_row (label or None), kph_tail_sign (+1|-1), tail_window ('centers'|'overlap')."""
    rows = []
    # --- tail rows
    tails = {"n48": ("tail_n48_L72_it3000", 12.0, 21.0), "n32": ("tail_n32_L48_it12000", 10.0, 20.0)}
    for tag, (key, lo, hi) in tails.items():
        if opts.get("tail_fields", "both") not in ("both", tag):
            continue
        if opts.get("tail_window", "centers") == "overlap":
            a, idx = tail_plateau(key, lo - 0.75, hi + 0.75)
        else:
            a, idx = tail_plateau(key, lo, hi)
        a["K_P_h"] *= opts.get("kph_tail_sign", 1.0)
        tol = opts.get("tail_tol", 3.0) * abs(FIXED["I1"] * a["I1"])
        rows.append(Row(f"tail_{tag}_up", f"tail_{tag}", a, tol))
        rows.append(Row(f"tail_{tag}_dn", f"tail_{tag}", {k: -v for k, v in a.items()}, tol))
    # --- sheets
    if opts.get("sheets", "both"):
        kinds = ["twist", "zigzag"] if opts["sheets"] == "both" else [opts["sheets"]]
        for kind in kinds:
            for w in (1.5, 3.0, 6.0):
                S, masked = sheet_S(kind, w, opts.get("sheet_mask", 2.5))
                cols = NAMES if opts.get("sheet_tol_all") else opts["free"] + list(FIXED)
                mx = max(abs(S[k]) for k in cols)
                tol = opts.get("sheet_tol", 1e-3) * mx if w == 1.5 else 0.0
                rows.append(Row(f"{kind}_w{w:g}", f"sheet_{kind}", {k: -v for k, v in S.items()}, -tol))
    # --- pairs
    if opts.get("pairs", True):
        P = ROWS["pairs"]
        D = {k: float(P["ansatz_same_d16"][k]) - float(P["ansatz_same_d8"][k]) for k in NAMES}
        mx = max(abs(D[k]) for k in opts["free"] + list(FIXED))
        rows.append(Row("pair_ansatz_16_8", "pair_ansatz", {k: -v for k, v in D.items()},
                        -opts.get("pair_tol", 1e-3) * mx))
        D = {k: float(P["relaxed_lam0_same_d24"][k]) - float(P["relaxed_lam0_same_d10"][k]) for k in NAMES}
        mx = max(abs(D[k]) for k in opts["free"] + list(FIXED))
        eps = opts.get("pair_tol", 1e-3) * mx if opts.get("relaxed_eps", "tol") == "tol" else 0.0
        rows.append(Row("pair_relaxed_24_10", "pair_relaxed", {k: -v for k, v in D.items()}, -eps))
        if opts.get("pairs_lam1"):
            D = {k: float(P["relaxed_lam1_same_d24"][k]) - float(P["relaxed_lam1_same_d10"][k]) for k in NAMES}
            mx = max(abs(D[k]) for k in opts["free"] + list(FIXED))
            rows.append(Row("pair_relaxed_lam1_24_10", "pair_relaxed", {k: -v for k, v in D.items()}, -opts.get("pair_tol", 1e-3) * mx))
    # --- coulomb
    if opts.get("coulomb", True):
        C = ROWS["coulomb"]
        E = {d: {k: float(C[f"coulomb_same_d{d}"][k]) for k in NAMES} for d in (12, 18, 24)}
        dI = E[12]["I1"] - E[24]["I1"]
        rho = {(i, j): 1.0 / i - 1.0 / j for (i, j) in ((12, 18), (12, 24), (18, 24))}
        K_cert = FIXED["I1"] * dI / rho[(12, 24)]
        RESULT.setdefault("coulomb", {})["K_cert"] = K_cert
        D = {k: E[12][k] - E[24][k] for k in NAMES}
        rows.append(Row("coulomb_rep", "coulomb", {k: -v for k, v in D.items()}, -0.25 * FIXED["I1"] * dI))
        for (i, j), r in (rho.items() if opts.get("coulomb_form", True) else []):
            D = {k: E[i][k] - E[j][k] for k in NAMES}
            rows.append(Row(f"coulomb_form_{i}_{j}_up", "coulomb", D, 0.1 * K_cert * r, aK=-r))
            rows.append(Row(f"coulomb_form_{i}_{j}_dn", "coulomb", {k: -v for k, v in D.items()},
                            0.1 * K_cert * r, aK=r))
        if opts.get("coulomb_form", True):
            rows.append(Row("coulomb_Kmin", "coulomb", {k: 0.0 for k in NAMES}, -0.5 * K_cert, aK=-1.0))
    # --- positivity channels
    if opts.get("channels", True):
        for name, q in ROWS["channels"].items():
            rows.append(Row(f"pos_{name}", "channels", {k: -float(q[k]) for k in NAMES}, 0.0))
        for kind in ("twist", "zigzag"):
            q = ROWS["sheets"][f"{kind}|nz256|w1.5"]["q_per_area"]
            rows.append(Row(f"pos_clock_{kind}", "channels", {k: -float(q[k]) for k in NAMES}, 0.0))
        if opts.get("channels_all_sheets"):
            for key, sh in ROWS["sheets"].items():
                rows.append(Row(f"pos_clock_{key}", "channels", {k: -float(sh["q_per_area"][k]) for k in NAMES}, 0.0))
    # --- UV plane-wave rows
    if opts.get("uv", True):
        for i, u in enumerate(ROWS["uv"]):
            rows.append(Row(f"uv_{i}_m{u['m']}_e{u['eps']:g}", "uv",
                            {k: -float(u["avg"][k]) for k in NAMES}, 0.0))
    if opts.get("quartic_rows", True) and any(k in opts["free"] for k in QUARTICS):
        rows += own_quartic_rows(seed=opts.get("quartic_seed", 20260905))
    # --- seeded UV polarization rows (an outer approximation of the PSD cone, own random vectors)
    if opts.get("uvq_seed", 0):
        rng = np.random.default_rng(7)
        Q = uvq_matrices(opts["free"])
        for d in DIRS:
            for c in range(opts["uvq_seed"]):
                v = rng.standard_normal(10); v /= np.linalg.norm(v)
                a = {k: -float(v @ Q[d][k] @ v) if k in Q[d] else 0.0 for k in NAMES}
                rows.append(Row(f"uvseed_{d}_{c}", "uvq_seed", a, 0.0))
    # --- mutation: flip one row
    if opts.get("flip_row"):
        for r in rows:
            if r.label == opts["flip_row"]:
                r.a = {k: -v for k, v in r.a.items()}
                r.aK = -r.aK
                r.rhs = -r.rhs
    return rows


# ============================================================ LP core
def to_matrix(rows, free, use_K):
    cols = list(free) + (["K"] if use_K else [])
    A = np.zeros((len(rows), len(cols)))
    b = np.zeros(len(rows))
    for i, r in enumerate(rows):
        for j, k in enumerate(free):
            A[i, j] = r.a.get(k, 0.0)
        if use_K:
            A[i, -1] = r.aK
        b[i] = r.rhs - sum(FIXED[f] * r.a.get(f, 0.0) for f in FIXED)
    # row scaling by the max abs coefficient (positive factor: inequality preserved)
    sc = np.maximum(np.max(np.abs(A), axis=1), 1e-300)
    sc = np.where(np.max(np.abs(A), axis=1) > 0, sc, 1.0)
    return A / sc[:, None], b / sc, cols, sc


def solve_lp(A, b, cols, objective="feas", norm_bound=None, method="highs-ds", warm_c=None):
    """x free; objective 'feas' (c = 0), 'l1' (min sum |x_k| over the basis columns, K excluded)."""
    m, n = A.shape
    nb = n - (1 if cols[-1] == "K" else 0)
    if objective == "feas" and norm_bound is None:
        res = linprog(c=np.zeros(n), A_ub=A, b_ub=b, bounds=[(None, None)] * n, method=method)
        return res, (res.x if res.status == 0 else None)
    # split x = p - q, p,q >= 0 on the basis columns; K free
    Ap = np.hstack([A[:, :nb], -A[:, :nb], A[:, nb:]])
    c = np.concatenate([np.ones(nb), np.ones(nb), np.zeros(n - nb)]) if objective == "l1" \
        else np.zeros(2 * nb + n - nb)
    bounds = [(0, None)] * (2 * nb) + [(None, None)] * (n - nb)
    A_ub, b_ub = Ap, b
    if norm_bound is not None:
        A_ub = np.vstack([Ap, np.concatenate([np.ones(2 * nb), np.zeros(n - nb)])[None, :]])
        b_ub = np.concatenate([b, [norm_bound]])
    res = linprog(c=c, A_ub=A_ub, b_ub=b_ub, bounds=bounds, method=method)
    if res.status != 0:
        return res, None
    x = np.concatenate([res.x[:nb] - res.x[nb:2 * nb], res.x[2 * nb:]])
    return res, x


KPH_CLIP = False


def uvq_matrices(free):
    Q = {}
    for d in DIRS:
        Q[d] = {k: np.array(ROWS["uvq"][d][k]) for k in ROWS["uvq"][d]}
        if KPH_CLIP:
            S = 0.5 * (Q[d]["K_P_h"] + Q[d]["K_P_h"].T)
            w, V = np.linalg.eigh(S)
            Q[d]["K_P_h"] = (V * np.maximum(w, 0.0)) @ V.T
    return Q


PSD_TOL = 1e-7


def psd_check(x, cols, Q, tol_rel=None):
    tol_rel = PSD_TOL if tol_rel is None else tol_rel
    """returns list of (dir, eigenvalue, eigenvector) for negative eigen-directions of sum chat_k Q_k."""
    cuts = []
    worst = 0.0
    for d in DIRS:
        S = FIXED["I1"] * Q[d]["I1"]
        for j, k in enumerate(cols):
            if k in Q[d]:
                S = S + x[j] * Q[d][k]
        S = 0.5 * (S + S.T)
        w, V = np.linalg.eigh(S)
        scale = max(np.max(np.abs(w)), 1e-30)
        worst = min(worst, w[0] / scale)
        for i in range(len(w)):
            if w[i] < -tol_rel * scale:
                cuts.append((d, float(w[i]), V[:, i]))
    return cuts, worst


def cut_rows(cuts, Q, tag):
    rows = []
    for c, (d, w, v) in enumerate(cuts):
        a = {k: -float(v @ Q[d][k] @ v) if k in Q[d] else 0.0 for k in NAMES}
        rows.append(Row(f"uvcut_{tag}_{d}_{c}", "uvq_cut", a, 0.0))
    return rows


def run_lp(free, opts, objective="feas", norm_bound=None, methods=("highs-ds", "highs-ipm"),
           use_cuts=True, max_rounds=300, verbose=True):
    """full loop: assemble, solve, cutting planes until PSD in all seven directions or infeasible.
    Returns dict with status, x (dict), norm, rows (final list), cuts added, per-method agreement."""
    opts = dict(opts, free=free)
    rows = assemble(opts)
    Q = uvq_matrices(free)
    use_K = opts.get("coulomb", True)
    out = {"n_rows_initial": len(rows), "cuts": 0, "rounds": 0}
    x_final, status = None, None
    for rnd in range(max_rounds + 1):
        A, b, cols, _ = to_matrix(rows, free, use_K)
        res_by = {}
        for mth in methods:
            res, x = solve_lp(A, b, cols, objective, norm_bound, mth)
            res_by[mth] = (res.status, x, res.fun if res.status == 0 else None)
        st = {mth: r[0] for mth, r in res_by.items()}
        out["status_by_method"] = st
        m0 = methods[0]
        x = res_by[m0][1]
        if x is None:
            # try the other method before declaring
            for mth in methods[1:]:
                if res_by[mth][1] is not None:
                    x = res_by[mth][1]
                    break
        if x is None:
            status = "infeasible" if all(s == 2 for s in st.values()) else f"nonopt:{st}"
            break
        if not use_cuts:
            status = "feasible"
            x_final = x
            break
        cuts, worst = psd_check(x, cols, Q)
        out["worst_rel_eig"] = worst
        if not cuts:
            status = "feasible_psd"
            x_final = x
            break
        keep = cuts
        nrm = float(sum(abs(x[j]) for j, k in enumerate(cols) if k != "K"))
        out.setdefault("norm_trajectory", []).append(nrm)
        if objective == "l1" and nrm > 1e8:
            status = "diverging"
            break
        rows += cut_rows(keep, Q, f"r{rnd}")
        out["cuts"] += len(keep)
        out["rounds"] = rnd + 1
    out["status"] = status if status is not None else "cut_loop_stalled"
    out["n_rows_final"] = len(rows)
    if x_final is not None:
        A, b, cols, _ = to_matrix(rows, free, use_K)
        xd = {k: float(x_final[j]) for j, k in enumerate(cols)}
        out["x"] = xd
        out["norm_l1"] = float(sum(abs(xd[k]) for k in free))
        out["max_violation"] = float(np.max(A @ x_final - b))
        # per-method agreement of the objective
        if objective == "l1":
            vals = {}
            for mth in methods:
                res, xx = solve_lp(A, b, cols, objective, norm_bound, mth)
                vals[mth] = float(sum(abs(xx[j]) for j, k in enumerate(cols) if k != "K")) if xx is not None else None
            out["norm_by_method"] = vals
    out["_rows"] = rows
    out["_cols"] = list(free) + (["K"] if use_K else [])
    if verbose:
        log(f"  LP free={len(free)} obj={objective} bound={norm_bound}: {status}, rows {len(rows)}, "
            f"cuts {out['cuts']}, norm {out.get('norm_l1')}")
    return out


# ============================================================ certificate
def farkas(rows, free, use_K, B=1e4, method="highs-ds"):
    """dual LP: min b^T y  s.t. A^T y = 0, sum y = 1, y >= 0 on the SCALED rows; then exact
    rational verification on the same scaled rows: beta + B ||r||_inf < 0."""
    A, b, cols, sc = to_matrix(rows, free, use_K)
    m, n = A.shape
    A_eq = np.vstack([A.T, np.ones((1, m))])
    b_eq = np.concatenate([np.zeros(n), [1.0]])
    res = linprog(c=b, A_eq=A_eq, b_eq=b_eq, bounds=[(0, None)] * m, method=method)
    if res.status != 0:
        return {"status": f"dual_status_{res.status}"}
    y = res.x
    yF = [Fraction(float(v)) for v in y]
    AF = [[Fraction(float(A[i, j])) for j in range(n)] for i in range(m)]
    bF = [Fraction(float(v)) for v in b]
    r = [sum(AF[i][j] * yF[i] for i in range(m)) for j in range(n)]
    beta = sum(bF[i] * yF[i] for i in range(m))
    r_inf = max(abs(v) for v in r)
    margin = beta + Fraction(B) * r_inf
    supp = [(rows[i].label, float(y[i])) for i in np.argsort(-y) if y[i] > 1e-9]
    r_by_col = {cols[j]: float(r[j]) for j in range(n)}
    return {"status": "certificate" if margin < 0 else "no_certificate",
            "beta": float(beta), "r_inf": float(r_inf), "margin": float(margin), "B": B,
            "support": supp[:12], "n_support": len(supp), "r_by_col": r_by_col,
            "y_min": float(np.min(y)), "sum_y": float(np.sum(y))}


# ============================================================ stage LP: A1..A6, A8
def stage_lp():
    R = {}
    base = dict(tail_tol=3.0, tail_fields="both", sheet_mask=2.5, sheet_tol=1e-3, pair_tol=1e-3,
                relaxed_eps="tol", coulomb=True, channels=True, uv=True, sheets="both", pairs=True,
                quartic_rows=True)
    # ---------- A1: two-derivative class
    log("A1: two-derivative class")
    a1 = run_lp(TWO_DERIV, base, "feas")
    R["A1"] = {k: v for k, v in a1.items() if not k.startswith("_")}
    if a1["status"] == "infeasible":
        cert = farkas(a1["_rows"], TWO_DERIV, True)
        R["A1"]["certificate"] = cert
        log(f"  certificate: {cert['status']} beta={cert['beta']:.3e} r_inf={cert['r_inf']:.3e} "
            f"margin={cert['margin']:.3e} support={cert['support'][:6]}")
        # sanity: a feasible-LP dual must return no certificate
    # A1 without cuts (is the infeasibility already there before the PSD cuts?)
    a1nc = run_lp(TWO_DERIV, base, "feas", use_cuts=False)
    R["A1"]["without_uvq_cuts"] = a1nc["status"]
    if a1nc["status"] != "infeasible" and "x" in a1nc:
        cuts, worst = psd_check(np.array([a1nc["x"][k] for k in a1nc["_cols"]]), a1nc["_cols"],
                                uvq_matrices(TWO_DERIV))
        R["A1"]["without_uvq_cuts_point"] = a1nc["x"]
        R["A1"]["without_uvq_cuts_worst_rel_eig"] = worst
    # any-means search: objective variations and row-group dropping
    groups = ["tail_n48", "tail_n32", "sheet_twist", "sheet_zigzag", "pair_ansatz", "pair_relaxed",
              "coulomb", "channels", "uv", "uvq_cut"]
    drop = {}
    for g in groups:
        o = dict(base)
        if g == "tail_n48":
            o["tail_fields"] = "n32"
        elif g == "tail_n32":
            o["tail_fields"] = "n48"
        elif g == "sheet_twist":
            o["sheets"] = "zigzag"
        elif g == "sheet_zigzag":
            o["sheets"] = "twist"
        elif g == "pair_ansatz" or g == "pair_relaxed":
            o["pairs"] = False          # both pair rows (finer split below)
        elif g == "coulomb":
            o["coulomb"] = False
        elif g == "channels":
            o["channels"] = False
        elif g == "uv":
            o["uv"] = False
        r = run_lp(TWO_DERIV, o, "l1", use_cuts=(g != "uvq_cut"), verbose=False)
        drop[g] = {"status": r["status"], "norm": r.get("norm_l1"), "x": r.get("x")}
        log(f"  drop {g}: {r['status']} norm {r.get('norm_l1')}")
    R["A1"]["drop_one_group"] = drop
    # drop two groups: tails together, sheets together, tail+twist
    for combo, o in (("both_tails", dict(base, tail_tol=1e9)),
                     ("both_sheets", dict(base, sheets=None)),
                     ("no_tails_no_sheets", dict(base, tail_tol=1e9, sheets=None)),
                     ("no_uv_no_cuts", dict(base, uv=False))):
        r = run_lp(TWO_DERIV, o, "l1", use_cuts=(combo != "no_uv_no_cuts"), verbose=False)
        drop[combo] = {"status": r["status"], "norm": r.get("norm_l1"), "x": r.get("x")}
        log(f"  drop {combo}: {r['status']} norm {r.get('norm_l1')}")
    # subsets: the producer's stated certificate support (tails + twist + UV + cuts) on its own;
    # Coulomb repulsion only (no form rows); the feasible point without Coulomb
    sub = {}
    for name, o in (("producer_support_tails_twist_uv", dict(base, sheets="twist", pairs=False, coulomb=False, channels=False)),
                    ("tails_twist_uv_channels", dict(base, sheets="twist", pairs=False, coulomb=False)),
                    ("coulomb_rep_only_no_form", dict(base, coulomb_form=False)),
                    ("coulomb_only", dict(base, tail_tol=1e9, sheets=None, pairs=False, channels=False, uv=False)),
                    ("coulomb_plus_tails", dict(base, sheets=None, pairs=False, channels=False, uv=False)),
                    ("coulomb_plus_uv", dict(base, tail_tol=1e9, sheets=None, pairs=False, channels=False)),
                    ("coulomb_plus_tails_plus_uv", dict(base, sheets=None, pairs=False, channels=False)),
                    ("coulomb_plus_zigzag", dict(base, tail_tol=1e9, sheets="zigzag", pairs=False, channels=False, uv=False))):
        r = run_lp(TWO_DERIV, o, "l1", verbose=False)
        sub[name] = {"status": r["status"], "norm": r.get("norm_l1"), "x": r.get("x"), "max_violation": r.get("max_violation")}
        if r["status"] == "infeasible":
            c = farkas(r["_rows"], TWO_DERIV, o.get("coulomb", True))
            sub[name]["certificate"] = {k: c[k] for k in ("status", "beta", "r_inf", "margin", "support")}
        log(f"  subset {name}: {r['status']} norm {r.get('norm_l1')} x {({k: round(v, 5) for k, v in r['x'].items() if abs(v) > 1e-6} if 'x' in r else '')}")
    R["A1"]["subsets"] = sub
    R["A1"]["drop_coulomb_point"] = drop["coulomb"]
    log(f"  drop-coulomb point: {({k: round(v, 5) for k, v in drop['coulomb']['x'].items() if abs(v) > 1e-6} if drop['coulomb']['x'] else None)}")
    # tail tolerance ladder and window sensitivity
    lad = {}
    for tol in (3.0, 10.0, 30.0, 100.0, 1000.0):
        r = run_lp(TWO_DERIV, dict(base, tail_tol=tol), "l1", verbose=False)
        lad[str(tol)] = {"status": r["status"], "norm": r.get("norm_l1"), "x": r.get("x")}
        log(f"  tail_tol {tol}: {r['status']} norm {r.get('norm_l1')}")
    R["A1"]["tail_tol_ladder"] = lad
    r = run_lp(TWO_DERIV, dict(base, tail_window="overlap"), "feas", verbose=False)
    R["A1"]["tail_window_overlap"] = r["status"]
    # mutations
    mut = {}
    for name, o in (("flip_twist_w1.5", dict(base, flip_row="twist_w1.5")),
                    ("flip_zigzag_w1.5", dict(base, flip_row="zigzag_w1.5")),
                    ("flip_tail_n48_up", dict(base, flip_row="tail_n48_up")),
                    ("kph_tail_sign_flipped", dict(base, kph_tail_sign=-1.0)),
                    ("no_sheet_artifact_mask", dict(base, sheet_mask=None)),
                    ("sheet_tol_0", dict(base, sheet_tol=0.0)),
                    ("relaxed_eps_0", dict(base, relaxed_eps=0))):
        r = run_lp(TWO_DERIV, o, "l1", verbose=False)
        mut[name] = {"status": r["status"], "norm": r.get("norm_l1"), "x": r.get("x")}
        if r["status"] == "infeasible":
            c = farkas(r["_rows"], TWO_DERIV, True)
            mut[name]["certificate"] = {k: c[k] for k in ("status", "beta", "r_inf", "margin", "support")}
        log(f"  mutation {name}: {r['status']} norm {r.get('norm_l1')} "
            f"{'cert ' + str(mut[name]['certificate']['support'][:4]) if 'certificate' in mut[name] else ''}")
    R["A1"]["mutations"] = mut

    # ---------- A2: + quartics
    log("A2: two-derivative + quartics")
    a2 = run_lp(TWO_DERIV + QUARTICS, base, "feas")
    R["A2"] = {k: v for k, v in a2.items() if not k.startswith("_")}
    if a2["status"] == "infeasible":
        c = farkas(a2["_rows"], TWO_DERIV + QUARTICS, True)
        R["A2"]["certificate"] = {k: c[k] for k in ("status", "beta", "r_inf", "margin", "support", "r_by_col")}
        log(f"  certificate: {c['status']} beta={c['beta']:.3e} margin={c['margin']:.3e} support={c['support'][:6]}")
    a2b = run_lp(TWO_DERIV + QUARTICS, dict(base, quartic_rows=False), "feas", verbose=False)
    R["A2"]["without_own_quartic_rows"] = a2b["status"]

    # ---------- A3: full B+ ladder
    log("A3: full B+ norm ladder")
    A3 = {}
    seeded = dict(base, uvq_seed=150)
    for bound in (30.0, 100.0, 300.0, 1000.0, None):
        r = run_lp(BPLUS, seeded, "l1", norm_bound=bound)
        A3[str(bound)] = {"status": r["status"], "norm": r.get("norm_l1"), "x": r.get("x"),
                          "cuts": r["cuts"], "norm_by_method": r.get("norm_by_method"),
                          "worst_rel_eig": r.get("worst_rel_eig")}
    R["A3"] = A3
    # feasibility point without cuts at bound 300 (is the PSD requirement what pushes the norm up?)
    r = run_lp(BPLUS, base, "l1", norm_bound=None, use_cuts=False)
    R["A3"]["min_norm_without_psd"] = {"status": r["status"], "norm": r.get("norm_l1"), "x": r.get("x")}
    r = run_lp(BPLUS, dict(seeded, quartic_rows=False), "l1", norm_bound=None)
    R["A3"]["min_norm_without_own_quartic_rows"] = {"status": r["status"], "norm": r.get("norm_l1"), "x": r.get("x")}

    # ---------- A4: without the zigzag sheet
    log("A4: B+ without the zigzag rows")
    r = run_lp(BPLUS, dict(base, sheets="twist"), "l1")
    R["A4"] = {"status": r["status"], "norm": r.get("norm_l1"), "x": r.get("x"), "cuts": r["cuts"]}
    r2 = run_lp(BPLUS, dict(base, sheets="zigzag", uvq_seed=150), "l1")
    R["A4"]["without_twist_instead"] = {"status": r2["status"], "norm": r2.get("norm_l1"), "x": r2.get("x")}
    r4 = run_lp(TWO_DERIV, dict(base, sheets="twist"), "l1", verbose=False)
    R["A4"]["two_deriv_without_zigzag"] = {"status": r4["status"], "norm": r4.get("norm_l1"), "x": r4.get("x")}

    # ---------- A5: one-field tail (n48 only)
    log("A5: n48-only tail, two-derivative class")
    r = run_lp(TWO_DERIV, dict(base, tail_fields="n48"), "l1")
    R["A5"] = {"status": r["status"], "norm": r.get("norm_l1"), "x": r.get("x"), "cuts": r["cuts"],
               "worst_rel_eig": r.get("worst_rel_eig"), "norm_by_method": r.get("norm_by_method")}
    if "x" in r:
        # tail residuals of that vertex on both fields
        xs = r["x"]
        res = {}
        for tag, key, lo, hi in (("n48", "tail_n48_L72_it3000", 12.0, 21.0), ("n32", "tail_n32_L48_it12000", 10.0, 20.0)):
            a, idx = tail_plateau(key, lo, hi)
            s = FIXED["I1"] * a["I1"] + sum(xs[k] * a[k] for k in TWO_DERIV)
            res[tag] = {"sum": s, "tol": 3.0 * abs(FIXED["I1"] * a["I1"]),
                        "terms": {k: xs[k] * a[k] for k in TWO_DERIV if abs(xs[k] * a[k]) > 1e-4}}
        R["A5"]["tail_residuals"] = res
    r32 = run_lp(TWO_DERIV, dict(base, tail_fields="n32"), "l1")
    R["A5"]["n32_only"] = {"status": r32["status"], "norm": r32.get("norm_l1"), "x": r32.get("x")}

    # ---------- A6: Coulomb Delta/rho
    log("A6: Coulomb")
    C = ROWS["coulomb"]
    E = {d: {k: float(C[f"coulomb_same_d{d}"][k]) for k in NAMES} for d in (12, 18, 24)}
    dr = {}
    for (i, j) in ((12, 18), (12, 24), (18, 24)):
        rho = 1.0 / i - 1.0 / j
        dr[f"{i}_{j}"] = FIXED["I1"] * (E[i]["I1"] - E[j]["I1"]) / rho
    R["A6"] = {"I1_delta_over_rho": dr, "K_cert": RESULT["coulomb"]["K_cert"]}
    # exact-1/d family with {K_lambda, R_etaMeta, R_hcov, K_P_h} at K = 2.655: 2 independent
    # equations (the third pair is their sum), 4 unknowns: minimum-norm solution
    els = ["K_lambda", "R_etaMeta", "R_hcov", "K_P_h"]
    Kt = 2.655
    M = np.array([[E[12][k] - E[18][k] for k in els], [E[18][k] - E[24][k] for k in els]])
    rhs = np.array([Kt * (1 / 12 - 1 / 18) - FIXED["I1"] * (E[12]["I1"] - E[18]["I1"]),
                    Kt * (1 / 18 - 1 / 24) - FIXED["I1"] * (E[18]["I1"] - E[24]["I1"])])
    xmin, *_ = np.linalg.lstsq(M, rhs, rcond=None)
    chk = M @ xmin - rhs
    R["A6"]["exact_form_K2.655_minnorm"] = {"x": dict(zip(els, map(float, xmin))), "residual": float(np.max(np.abs(chk))),
                                            "K_ge_Kcert_half": bool(Kt >= 0.5 * RESULT["coulomb"]["K_cert"]),
                                            "Kcert_half": 0.5 * RESULT["coulomb"]["K_cert"]}
    # the K range reachable by exact form with these four elements under the LP's own K bound
    # (LP over the coulomb rows only, free = els): min and max K
    o = dict(base, tail_tol=1e9, sheets=None, pairs=False, channels=False, uv=False, quartic_rows=False)
    rows = assemble(dict(o, free=els))
    A, b, cols, _ = to_matrix(rows, els, True)
    kr = {}
    for sgn, nm in ((1.0, "min_K"), (-1.0, "max_K")):
        c = np.zeros(len(cols)); c[-1] = sgn
        res = linprog(c=c, A_ub=A, b_ub=b, bounds=[(None, None)] * len(cols), method="highs-ds")
        kr[nm] = float(res.x[-1]) if res.status == 0 else f"status_{res.status}"
    R["A6"]["K_range_coulomb_rows_only"] = kr

    # ---------- A8: null columns
    log("A8: parity-odd columns")
    rows = assemble(dict(base, free=BPLUS + PARITY_ODD))
    rel_rowmax, rel_I1 = {k: 0.0 for k in PARITY_ODD}, {k: 0.0 for k in PARITY_ODD}
    where = {}
    for r in rows:
        mx = max(abs(r.a.get(k, 0.0)) for k in BPLUS + list(FIXED))
        i1 = abs(FIXED["I1"] * r.a.get("I1", 0.0))
        for k in PARITY_ODD:
            v = abs(r.a.get(k, 0.0))
            if mx > 0 and v / mx > rel_rowmax[k]:
                rel_rowmax[k] = v / mx
                where[k + "_rowmax"] = r.label
            if i1 > 0 and v / i1 > rel_I1[k]:
                rel_I1[k] = v / i1
                where[k + "_I1"] = r.label
    R["A8"] = {"max_rel_to_row_max": rel_rowmax, "max_rel_to_4I1_column": rel_I1, "where": where,
               "raw_pairs": {lab: {k: float(ROWS["pairs"][lab][k]) for k in PARITY_ODD}
                             for lab in ROWS["pairs"]}}
    # does admitting E1..E3 as free columns change A1 / A3?
    r = run_lp(TWO_DERIV + PARITY_ODD, base, "l1", verbose=False)
    R["A8"]["two_deriv_plus_E"] = {"status": r["status"], "norm": r.get("norm_l1"), "x": r.get("x")}
    r = run_lp(BPLUS + PARITY_ODD, base, "l1", norm_bound=300.0, verbose=False)
    R["A8"]["bplus_plus_E_bound300"] = {"status": r["status"], "norm": r.get("norm_l1"), "x": r.get("x")}
    return R


# ============================================================ own entrant densities
def own_jets(M, h, br):
    """spatial jets (4, n, n, n, 4, 4) with A_0 = 0 on the fwd or bwd branch (own stencil)."""
    A = np.zeros((4,) + M.shape)
    for ax in range(3):
        d = np.zeros_like(M)
        sl_hi = [slice(None)] * 3; sl_lo = [slice(None)] * 3
        sl_hi[ax] = slice(1, None); sl_lo[ax] = slice(0, -1)
        diff = (M[tuple(sl_hi)] - M[tuple(sl_lo)]) / h
        if br == "fwd":
            d[tuple(sl_lo)] = diff
        else:
            d[tuple(sl_hi)] = diff
        A[1 + ax] = d
    return A


def own_eig(M):
    N = M @ ETA
    lam, V = np.linalg.eig(N)
    lam, V = lam.real, V.real
    nrm = np.einsum("...ak,a,...ak->...k", V, ETA_D, V)
    V = V / np.sqrt(np.abs(nrm))[..., None, :]
    return lam, V, np.sign(nrm)


def own_T1(A):
    return sum(np.einsum("...ab,bc,...cd,da->...", A[i], ETA, A[i], ETA) for i in (1, 2, 3))


def own_T3(A):
    div = sum(A[mu][..., mu, :] for mu in (1, 2, 3))          # div^b = sum_mu (A_mu)^{mu b}
    return np.einsum("...b,b,...b->...", div, ETA_D, div)


def own_KPh(A, M, roots=(-8.0, 1.0)):
    lam, V, sig = own_eig(M)
    f = (lam - roots[0]) * (lam - roots[1])
    tot = 0.0
    for i in (1, 2, 3):
        B = np.einsum("...ka,ab,...bc,cd,...dl->...kl", V.swapaxes(-1, -2), ETA, A[i], ETA, V)
        B = B * sig[..., :, None]
        W = f[..., :, None] * f[..., None, :] * B
        tot = tot + np.sum(W * W, axis=(-1, -2))
    return 0.5 * tot


def own_RG(A, G):
    """R_G = G_cd T^cd, T^{cd} = sum_{mu nu} [(A_mu)^{nu c}(A_nu)^{mu d} - (A_mu)^{mu c}(A_nu)^{nu d}], mu,nu spatial."""
    first = 0.0
    div = 0.0
    for mu in (1, 2, 3):
        div = div + A[mu][..., mu, :]
        for nu in (1, 2, 3):
            first = first + np.einsum("...c,...d->...cd", A[mu][..., nu, :], A[nu][..., mu, :])
    T = first - np.einsum("...c,...d->...cd", div, div)
    return np.einsum("...cd,...cd->...", G, T)


def own_I1(A):
    tot = 0.0
    for i in (1, 2, 3):
        for j in range(i + 1, 4):
            F = A[i] @ ETA @ A[j] - A[j] @ ETA @ A[i]
            tot = tot + np.einsum("...ab,...cd,ac,bd->...", F, F, ETA, ETA)
    return tot


def own_static_density(name, M, h):
    """h^3-weighted per-cell static density on the sym stencil (own fwd/bwd average)."""
    d = 0.0
    for br in ("fwd", "bwd"):
        A = own_jets(M, h, br)
        if name == "T1":
            v = own_T1(A)
        elif name == "T3":
            v = own_T3(A)
        elif name == "K_P_h":
            v = own_KPh(A, M)
        elif name == "R_etaMeta":
            v = own_RG(A, ETA @ M @ ETA)
        elif name == "I1":
            v = own_I1(A)
        else:
            raise ValueError(name)
        d = d + 0.5 * v
    return h ** 3 * d


# ============================================================ stage tail: A9
def stage_tail():
    R = {}
    fields = {"n48": (os.path.join(CK13, "seed_n48_L72_it3000.npy"), 48, 72.0, 12.0, 21.0,
                      "tail_n48_L72_it3000"),
              "n32": (os.path.join(CK10, "relax_g8_n32_L48_it12000.npy"), 32, 48.0, 10.0, 20.0,
                      "tail_n32_L48_it12000")}
    for tag, (path, n, L, lo, hi, key) in fields.items():
        log(f"A9 tail re-measure on {tag}")
        M = np.load(path)
        h = L / n
        cfg = R13W.cfg_of(n, L)
        x = (np.arange(n) - (n - 1) / 2.0) * h
        X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
        r = np.sqrt(X * X + Y * Y + Z * Z)
        out = {"n": n, "L": L, "h": h, "corner_spectrum": sorted(np.linalg.eigvals(M[0, 0, 0] @ ETA).real.tolist())}
        for name in ("I1", "K_P_h", "T1", "T3", "R_etaMeta"):
            dens = own_static_density(name, M, h)
            total = float(np.sum(dens))
            # own shell binning: width 1.5 (edges at multiples of 1.5), and width 1.0 as a second binning
            rec = {"total": total}
            for wsh in (1.5, 1.0):
                nb = int(np.ceil(r.max() / wsh)) + 1
                idx = np.minimum((r / wsh).astype(int), nb - 1)
                per = np.bincount(idx.ravel(), weights=dens.ravel(), minlength=nb) / wsh
                centers = (np.arange(nb) + 0.5) * wsh
                sel = (centers >= lo) & (centers <= hi)
                rec[f"plateau_w{wsh:g}"] = float(np.mean(per[sel]))
                rec[f"plateau_w{wsh:g}_std"] = float(np.std(per[sel]))
                if wsh == 1.5:
                    rec["shell_per_r_w1.5"] = per[:len(ROWS[key][name]["shell_per_r"])].tolist()
                    prod = np.array(ROWS[key][name]["shell_per_r"])
                    mine = per[:len(prod)]
                    rec["max_rel_diff_vs_rows_json"] = float(np.max(np.abs(mine - prod)) / max(np.max(np.abs(prod)), 1e-30))
            rec["rows_json_s"] = float(ROWS[key][name]["s"])
            rec["rows_json_plateau_centers"] = tail_plateau(key, lo, hi)[0][name]
            # control against the producer's registry (allowed): static_energy
            if name in ("K_P_h", "R_etaMeta", "K_lambda"):
                rec["registry_static_energy"] = float(R14T.static_energy(name, M, cfg))
            if name == "I1":
                rec["certified_E_u_over_4"] = float(INS4.e_parts(M, cfg)[0] / 4.0)
            out[name] = rec
            log(f"  {name}: total {total:.6g} plateau(w1.5) {rec['plateau_w1.5']:.6g} "
                f"plateau(w1.0) {rec['plateau_w1']:.6g} rows.json {rec['rows_json_plateau_centers']:.6g}")
        # mutation: the wrong-sign 2uu term in H (plain eta in place of the eta-orthonormal frame) and
        # a literal-root K_P_h: the tail plateau must change
        for mname, roots in (("literal_roots_(8,1)", (8.0, 1.0)),):
            d = 0.0
            for br in ("fwd", "bwd"):
                d = d + 0.5 * own_KPh(own_jets(M, h, br), M, roots)
            d = h ** 3 * d
            nb = int(np.ceil(r.max() / 1.5)) + 1
            idx = np.minimum((r / 1.5).astype(int), nb - 1)
            per = np.bincount(idx.ravel(), weights=d.ravel(), minlength=nb) / 1.5
            centers = (np.arange(nb) + 0.5) * 1.5
            out[f"K_P_h_mutant_{mname}_plateau"] = float(np.mean(per[(centers >= lo) & (centers <= hi)]))
        R[tag] = out
    return R


# ============================================================ stage uvq: A7
def slab_Q(name, dirvec, nz=64, Lz=24.0, m=2, eps=0.01, stencil="sym", periodic=False, M0=None):
    """own plane-wave construction: M = M0 + eps cos(k z') X_I along the direction dirvec on a
    4 x 4 x nz slab (transverse cells 4, h = Lz/nz); returns the 10x10 polarization matrix of the
    cycle-averaged static density divided by eps^2 (off-diagonal basis directions normalized 1/sqrt2)."""
    if M0 is None:
        M0 = np.diag([8.0, 1.0, 0.3, 0.0])
    h = Lz / nz
    # basis of symmetric directions, order 00 01 02 03 11 12 13 22 23 33
    basis, labels = [], []
    for a in range(4):
        for b in range(a, 4):
            E = np.zeros((4, 4))
            if a == b:
                E[a, a] = 1.0
            else:
                E[a, b] = E[b, a] = 1.0 / np.sqrt(2.0)
            basis.append(E); labels.append(f"{a}{b}")
    dv = np.array(dirvec, float)
    dv = dv / np.linalg.norm(dv)
    # the slab is 3D: the wave runs along dv in a periodic box so that k.dv is a lattice mode
    nx = nz if abs(dv[0]) > 0 else 4
    ny = nz if abs(dv[1]) > 0 else 4
    nzz = nz if abs(dv[2]) > 0 else 4
    xs = [(np.arange(nn) - (nn - 1) / 2.0) * h for nn in (nx, ny, nzz)]
    X, Y, Z = np.meshgrid(*xs, indexing="ij")
    k = 2 * np.pi * m / Lz
    phase = k * (dv[0] * X + dv[1] * Y + dv[2] * Z)
    c = np.cos(phase)

    def dens_of(Mf):
        d = 0.0
        brs = ("fwd", "bwd") if stencil == "sym" else (stencil,)
        for br in brs:
            A = own_jets(Mf, h, br)
            if periodic:
                for ax in range(3):
                    A[1 + ax] = (np.roll(Mf, -1, axis=ax) - Mf) / h if br == "fwd" else (Mf - np.roll(Mf, 1, axis=ax)) / h
            if name == "T1":
                v = own_T1(A)
            elif name == "K_P_h":
                v = own_KPh(A, Mf)
            elif name == "T3":
                v = own_T3(A)
            else:
                raise ValueError(name)
            d = d + v / len(brs)
        return d

    Q = np.zeros((10, 10))
    vals = {}
    # interior average (the fwd/bwd stencils lose one boundary layer; average over interior cells)
    inner = tuple(slice(1, -1) if nn > 4 else slice(None) for nn in (nx, ny, nzz))
    for I in range(10):
        Mf = M0[None, None, None] + eps * c[..., None, None] * basis[I][None, None, None]
        vals[I] = float(np.mean(dens_of(Mf)[inner])) / eps ** 2
        Q[I, I] = vals[I]
    for I in range(10):
        for J in range(I + 1, 10):
            Mf = M0[None, None, None] + eps * c[..., None, None] * (basis[I] + basis[J])[None, None, None]
            v = float(np.mean(dens_of(Mf)[inner])) / eps ** 2
            Q[I, J] = Q[J, I] = 0.5 * (v - vals[I] - vals[J])
    return Q, labels, k, h


def stage_uvq():
    R = {}
    Qp = {d: {k: np.array(ROWS["uvq"][d][k]) for k in ROWS["uvq"][d]} for d in DIRS}
    # (a) PSD / indefinite bookkeeping of the producer's matrices, all seven directions
    eig = {}
    for k in ("K_P_h", "T1", "K_T", "T3", "T2", "T4", "K_lambda", "R_etaMeta", "R_hcov", "I1"):
        eig[k] = {}
        for d in DIRS:
            S = 0.5 * (Qp[d][k] + Qp[d][k].T)
            w = np.linalg.eigvalsh(S)
            eig[k][d] = {"min": float(w[0]), "max": float(w[-1]), "asym": float(np.max(np.abs(Qp[d][k] - Qp[d][k].T)))}
    R["producer_eigs"] = eig
    R["K_P_h_psd_all_dirs"] = bool(all(eig["K_P_h"][d]["min"] >= -1e-9 * eig["K_P_h"][d]["max"] for d in DIRS))
    R["T1_indefinite_all_dirs"] = bool(all(eig["T1"][d]["min"] < 0 < eig["T1"][d]["max"] for d in DIRS))
    R["K_T_indefinite_dirs"] = [d for d in DIRS if eig["K_T"][d]["min"] < -1e-9 * max(eig["K_T"][d]["max"], 1e-30)]
    # where is T1 negative, and is any other class element positive there? (the UV lever of A1)
    lever = {}
    for d in DIRS:
        S = 0.5 * (Qp[d]["T1"] + Qp[d]["T1"].T)
        w, V = np.linalg.eigh(S)
        v = V[:, 0]
        lever[d] = {"T1_min_eig": float(w[0]), "vec": v.round(4).tolist(),
                    "others_on_vec": {k: float(v @ Qp[d][k] @ v) for k in TWO_DERIV}}
    R["T1_negative_direction_lever"] = lever
    # (b) own re-measure along z, T1 and K_P_h, amplitudes 0.01 and 0.003
    own = {}
    for name in ("T1", "K_P_h"):
        own[name] = {}
        for eps in (0.01, 0.003):
            Q, labels, k, h = slab_Q(name, (0, 0, 1), eps=eps)
            own[name][str(eps)] = Q.tolist()
        Qa, Qb = np.array(own[name]["0.01"]), np.array(own[name]["0.003"])
        own[name]["eps_convergence_max_abs_diff"] = float(np.max(np.abs(Qa - Qb)))
        own[name]["labels"] = labels
        own[name]["k"] = k
        own[name]["h"] = h
    R["own_z"] = own
    # normalization: my T1 [00] entry vs the producer's; then compare the STRUCTURE (ratio-free)
    Pz = {k: Qp["0|0|1"][k] for k in ("T1", "K_P_h")}
    myT1 = np.array(own["T1"]["0.003"])
    myKP = np.array(own["K_P_h"]["0.003"])
    scale = Pz["T1"][0, 0] / myT1[0, 0]
    R["norm_scale_producer_over_own_T1_00"] = float(scale)
    R["T1_z_max_abs_diff_after_scale"] = float(np.max(np.abs(myT1 * scale - Pz["T1"])))
    R["KPh_z_max_abs_diff_after_scale"] = float(np.max(np.abs(myKP * scale - Pz["K_P_h"])))
    R["KPh_z_max_rel_diff_after_scale"] = float(np.max(np.abs(myKP * scale - Pz["K_P_h"])) / np.max(np.abs(Pz["K_P_h"])))
    R["KPh_z_diag_ratio_producer_over_own_scaled"] = [float(Pz["K_P_h"][i, i] / (myKP[i, i] * scale)) if abs(myKP[i, i] * scale) > 1e-6 else None for i in range(10)]
    # analytic expectation at the vacuum for K_P_h (quadratic order): (1/2) f_a^2 f_b^2 |B|^2
    f2, f3 = (0.3 + 8) * (0.3 - 1), (0.0 + 8) * (0.0 - 1)
    R["KPh_analytic_over_T1_00"] = {"22": 0.5 * f2 ** 4, "33": 0.5 * f3 ** 4, "23": f2 ** 2 * f3 ** 2 / 2}
    R["KPh_producer_over_T1_00"] = {"22": float(Pz["K_P_h"][7, 7] / Pz["T1"][0, 0]),
                                    "33": float(Pz["K_P_h"][9, 9] / Pz["T1"][0, 0]),
                                    "23": float(Pz["K_P_h"][8, 8] / Pz["T1"][0, 0])}
    R["KPh_own_over_T1_00"] = {e: {"22": float(np.array(own["K_P_h"][e])[7, 7] / np.array(own["T1"][e])[0, 0]),
                                   "33": float(np.array(own["K_P_h"][e])[9, 9] / np.array(own["T1"][e])[0, 0]),
                                   "23": float(np.array(own["K_P_h"][e])[8, 8] / np.array(own["T1"][e])[0, 0])}
                               for e in ("0.01", "0.003")}
    # which (nz, Lz, m, stencil) reproduces the producer's T1 [00] = 0.14402 (normalization hunt)
    hunt = []
    target = float(Pz["T1"][0, 0])
    for nz, Lz in ((64, 24.0), (64, 12.0), (128, 24.0), (32, 24.0)):
        for m in (1, 2, 3, 4):
            for st in ("sym", "central"):
                h = Lz / nz
                k = 2 * np.pi * m / Lz
                if st == "sym":
                    val = (1 - np.cos(k * h)) / h ** 2
                else:
                    val = 0.5 * (np.sin(k * h) / h) ** 2
                hunt.append((abs(val - target), nz, Lz, m, st, val))
    hunt.sort()
    R["T1_00_normalization_hunt_top5"] = [{"nz": a[1], "Lz": a[2], "m": a[3], "stencil": a[4], "value": a[5]} for a in hunt[:5]]
    R["T1_00_producer"] = target
    R["T1_00_continuum_k2_over_2_for_k"] = float(np.sqrt(2 * target))
    # own PSD / indefinite along z at eps 0.003
    R["own_KPh_z_min_eig"] = float(np.linalg.eigvalsh(0.5 * (myKP + myKP.T))[0])
    R["own_KPh_z_max_eig"] = float(np.linalg.eigvalsh(0.5 * (myKP + myKP.T))[-1])
    R["own_T1_z_eigs"] = np.linalg.eigvalsh(0.5 * (myT1 + myT1.T)).round(6).tolist()
    # mutation: literal roots (8, 1): K_P_h's Q along z must acquire boost entries / lose PSD
    Qm, *_ = slab_Q("K_P_h", (0, 0, 1), eps=0.003)
    Ql = np.zeros((10, 10))
    # rebuild with literal roots through a tiny closure
    M0 = np.diag([8.0, 1.0, 0.3, 0.0])
    def kph_lit(A, Mf):
        return own_KPh(A, Mf, roots=(8.0, 1.0))
    # quick diagonal-only literal read at (0,3) and (0,0)
    lit = {}
    for I, lab in ((3, "03"), (0, "00"), (9, "33")):
        E = np.zeros((4, 4))
        if lab == "03":
            E[0, 3] = E[3, 0] = 1 / np.sqrt(2)
        elif lab == "00":
            E[0, 0] = 1.0
        else:
            E[3, 3] = 1.0
        h = 24.0 / 64
        z = (np.arange(64) - 31.5) * h
        c = np.cos(2 * np.pi * 2 / 24.0 * z)
        Mf = M0[None, None, None] + 0.003 * c[None, None, :, None, None] * E[None, None, None]
        Mf = np.broadcast_to(Mf, (4, 4, 64, 4, 4)).copy()
        d = 0.0
        for br in ("fwd", "bwd"):
            d = d + 0.5 * kph_lit(own_jets(Mf, h, br), Mf)
        lit[lab] = float(np.mean(d[:, :, 1:-1])) / 0.003 ** 2
    R["mutation_literal_roots_KPh_diag_z"] = lit
    return R


# ============================================================ stage extra: robustness
def stage_extra():
    R = {}
    base = dict(tail_tol=3.0, tail_fields="both", sheet_mask=2.5, sheet_tol=1e-3, pair_tol=1e-3,
                relaxed_eps="tol", coulomb=True, channels=True, uv=True, sheets="both", pairs=True,
                quartic_rows=True)
    # A3 robustness: more seeded rows, tighter PSD tolerance, second quartic seed, lam1 pairs added
    rob = {}
    for name, o, tol in (("seed150_tol1e-7", dict(base, uvq_seed=150), 1e-7),
                         ("seed600_tol1e-9", dict(base, uvq_seed=600), 1e-9),
                         ("seed600_tol1e-9_quartic_seed2", dict(base, uvq_seed=600, quartic_seed=2), 1e-9)):
        global PSD_TOL
        PSD_TOL = tol
        r = run_lp(BPLUS, o, "l1", norm_bound=None, verbose=False, max_rounds=400)
        x = r.get("x")
        rec = {"status": r["status"], "norm": r.get("norm_l1"), "cuts": r["cuts"], "x": x,
               "worst_rel_eig": r.get("worst_rel_eig"), "norm_by_method": r.get("norm_by_method")}
        if x is not None:
            # exact PSD margins of the point in all seven directions
            Q = uvq_matrices(BPLUS)
            cols = r["_cols"]
            xv = np.array([x[k] for k in cols])
            margins = {}
            for d in DIRS:
                S = FIXED["I1"] * Q[d]["I1"] + sum(xv[j] * Q[d][k] for j, k in enumerate(cols) if k in Q[d])
                w = np.linalg.eigvalsh(0.5 * (S + S.T))
                margins[d] = {"min_eig": float(w[0]), "max_eig": float(w[-1])}
            rec["psd_margins"] = margins
            rec["quartic_coeffs"] = {k: x[k] for k in QUARTICS}
            rec["max_violation"] = r.get("max_violation")
        rob[name] = rec
        log(f"  A3 robustness {name}: {r['status']} norm {r.get('norm_l1')} x {({k: round(v, 3) for k, v in x.items() if abs(v) > 1e-3} if x else None)}")
    PSD_TOL = 1e-7
    R["A3_robustness"] = rob
    # witness check: the producer's stated support rows (tails + twist + uv + full PSD) on the point
    # x = {K_T: 0.0654, T1: 0.001} found with the channels added
    o = dict(base, sheets="twist", pairs=False, coulomb=False, channels=False, free=TWO_DERIV)
    rows = assemble(o)
    A, b, cols, _ = to_matrix(rows, TWO_DERIV, False)
    xw = np.zeros(len(cols))
    xw[cols.index("K_T")] = 0.0654041
    xw[cols.index("T1")] = 0.0010
    Q = uvq_matrices(TWO_DERIV)
    cuts, worst = psd_check(xw, cols, Q, tol_rel=1e-12)
    R["producer_support_witness"] = {"x": {"K_T": 0.0654041, "T1": 0.001},
                                     "max_row_violation": float(np.max(A @ xw - b)),
                                     "worst_rel_eig_all_dirs": worst, "n_negative_directions": len(cuts),
                                     "n_rows": len(rows)}
    log(f"  witness on producer-support rows: max violation {R['producer_support_witness']['max_row_violation']:.3e}, worst rel eig {worst:.3e}")
    # the scan of the twist-row lever: min T1 forced by the twist row, and K_T needed per unit T1 on UV
    o = dict(base, free=TWO_DERIV, sheets="twist", pairs=False, coulomb=False, channels=False, tail_tol=1e9)
    rows = assemble(o)
    A, b, cols, _ = to_matrix(rows, TWO_DERIV, False)
    for tag, cvec in (("min_T1", {"T1": 1.0}), ("min_K_T", {"K_T": 1.0})):
        c = np.zeros(len(cols))
        for k, v in cvec.items():
            c[cols.index(k)] = v
        res = linprog(c=c, A_ub=A, b_ub=b, bounds=[(None, None)] * len(cols), method="highs-ds")
        R[f"twist_uv_{tag}"] = float(res.fun) if res.status == 0 else f"status_{res.status}"
    # uv_66 row content (the UV row in the audit certificates)
    u = ROWS["uv"][66]
    R["uv_66"] = {"m": u["m"], "eps": u["eps"], "avg_class": {k: float(u["avg"][k]) for k in TWO_DERIV + ["I1"]}}
    return R


# ============================================================ stage extra2
def stage_extra2():
    global KPH_CLIP, PSD_TOL
    R = {}
    base = dict(tail_tol=3.0, tail_fields="both", sheet_mask=2.5, sheet_tol=1e-3, pair_tol=1e-3,
                relaxed_eps="tol", coulomb=True, channels=True, uv=True, sheets="both", pairs=True,
                quartic_rows=True, uvq_seed=300)
    # stalled subsets, seeded
    sub = {}
    for name, o in (("producer_support_tails_twist_uv", dict(base, sheets="twist", pairs=False, coulomb=False, channels=False)),
                    ("coulomb_only", dict(base, tail_tol=1e9, sheets=None, pairs=False, channels=False, uv=False)),
                    ("coulomb_plus_tails", dict(base, sheets=None, pairs=False, channels=False, uv=False)),
                    ("coulomb_plus_uv", dict(base, tail_tol=1e9, sheets=None, pairs=False, channels=False)),
                    ("coulomb_plus_tails_plus_uv", dict(base, sheets=None, pairs=False, channels=False)),
                    ("coulomb_plus_tails_plus_channels", dict(base, sheets=None, pairs=False, uv=False)),
                    ("coulomb_plus_tails_plus_uv_plus_channels", dict(base, sheets=None, pairs=False)),
                    ("no_tails_no_sheets", dict(base, tail_tol=1e9, sheets=None))):
        r = run_lp(TWO_DERIV, o, "l1", verbose=False, max_rounds=200)
        traj = r.get("norm_trajectory", [])
        sub[name] = {"status": r["status"], "norm": r.get("norm_l1"), "x": r.get("x"), "cuts": r["cuts"],
                     "trajectory_head": traj[:3], "trajectory_tail": traj[-3:]}
        if r["status"] == "infeasible":
            c = farkas(r["_rows"], TWO_DERIV, o.get("coulomb", True))
            sub[name]["certificate"] = {k: c[k] for k in ("status", "beta", "r_inf", "margin", "support")}
        log(f"  seeded subset {name}: {r['status']} norm {r.get('norm_l1')} traj {traj[:2]}..{traj[-2:]} "
            f"x {({k: round(v, 4) for k, v in r['x'].items() if abs(v) > 1e-4} if 'x' in r else '')}")
    R["seeded_subsets"] = sub
    # A3 robustness under extra rows
    rob = {}
    for name, o in (("all_sheet_clock_channels", dict(base, channels_all_sheets=True)),
                    ("plus_lam1_pair", dict(base, pairs_lam1=True)),
                    ("sheet_tol_over_all_columns", dict(base, sheet_tol_all=True)),
                    ("all_three", dict(base, channels_all_sheets=True, pairs_lam1=True, sheet_tol_all=True)),
                    ("quartic_seed2_tol1e-7", dict(base, quartic_seed=2))):
        r = run_lp(BPLUS, o, "l1", verbose=False, max_rounds=300)
        rob[name] = {"status": r["status"], "norm": r.get("norm_l1"), "x": r.get("x"), "cuts": r["cuts"],
                     "worst_rel_eig": r.get("worst_rel_eig"), "trajectory_tail": r.get("norm_trajectory", [])[-3:]}
        log(f"  A3 extra rows {name}: {r['status']} norm {r.get('norm_l1')} x {({k: round(v, 3) for k, v in r['x'].items() if abs(v) > 1e-3} if 'x' in r else '')}")
    R["A3_extra_rows"] = rob
    # K_T forced by the boost channel? min K_T over the full B+ rows (no PSD cuts needed for a lower bound)
    o = dict(base, free=BPLUS)
    rows = assemble(o)
    A, b, cols, _ = to_matrix(rows, BPLUS, True)
    c = np.zeros(len(cols)); c[cols.index("K_T")] = 1.0
    res = linprog(c=c, A_ub=A, b_ub=b, bounds=[(None, None)] * len(cols), method="highs-ds")
    q = ROWS["channels"]["hedgehog_boost_z"]
    R["K_T_min_over_bplus_rows"] = float(res.fun) if res.status == 0 else f"status_{res.status}"
    R["K_T_min_two_deriv_rows"] = None
    rows2 = assemble(dict(o, free=TWO_DERIV))
    A2, b2, cols2, _ = to_matrix(rows2, TWO_DERIV, True)
    c2 = np.zeros(len(cols2)); c2[cols2.index("K_T")] = 1.0
    res2 = linprog(c=c2, A_ub=A2, b_ub=b2, bounds=[(None, None)] * len(cols2), method="highs-ds")
    R["K_T_min_two_deriv_rows"] = float(res2.fun) if res2.status == 0 else f"status_{res2.status}"
    R["hedgehog_boost_z_bound"] = {"4I1": 4 * float(q["I1"]), "K_T": float(q["K_T"]), "T1": float(q["T1"]), "T3": float(q["T3"]),
                                   "K_T_min_if_alone": -4 * float(q["I1"]) / float(q["K_T"])}
    log(f"  min K_T: B+ rows {R['K_T_min_over_bplus_rows']}, two-deriv rows {R['K_T_min_two_deriv_rows']}, boost_z alone {R['hedgehog_boost_z_bound']['K_T_min_if_alone']:.5f}")
    # K_P_h Q clipped to PSD: does A1 / A3 change?
    KPH_CLIP = True
    r1 = run_lp(TWO_DERIV, dict(base, uvq_seed=0), "feas", verbose=False)
    r3 = run_lp(BPLUS, base, "l1", verbose=False)
    KPH_CLIP = False
    R["kph_clipped"] = {"A1_status": r1["status"], "A3_status": r3["status"], "A3_norm": r3.get("norm_l1"), "A3_x": r3.get("x")}
    log(f"  K_P_h Q clipped to PSD: A1 {r1['status']}, A3 {r3['status']} norm {r3.get('norm_l1')}")
    # eps scaling of the spurious negativity of the own K_P_h Q along z
    sc = {}
    for eps in (0.003, 0.01, 0.03):
        Q, labels, k, h = slab_Q("K_P_h", (0, 0, 1), eps=eps)
        w = np.linalg.eigvalsh(0.5 * (Q + Q.T))
        sc[str(eps)] = {"min_eig": float(w[0]), "max_eig": float(w[-1]), "diag_22": float(Q[7, 7]), "diag_33": float(Q[9, 9]), "diag_23": float(Q[8, 8]),
                        "diag_12": float(Q[5, 5]), "diag_13": float(Q[6, 6])}
    R["own_KPh_z_eps_scaling"] = sc
    Pz = np.array(ROWS["uvq"]["0|0|1"]["K_P_h"])
    R["producer_KPh_z_min_eig"] = float(np.linalg.eigvalsh(0.5 * (Pz + Pz.T))[0])
    log(f"  own K_P_h z min eig by eps: {sc}; producer {R['producer_KPh_z_min_eig']:.4e}")
    return R


# ============================================================ stage extra3: subsets at PSD tol 1e-6
def stage_extra3():
    global PSD_TOL
    R = {}
    base = dict(tail_tol=3.0, tail_fields="both", sheet_mask=2.5, sheet_tol=1e-3, pair_tol=1e-3,
                relaxed_eps="tol", coulomb=True, channels=True, uv=True, sheets="both", pairs=True,
                quartic_rows=True, uvq_seed=300)
    PSD_TOL = 1e-6
    sub = {}
    for name, o in (("producer_support_tails_twist_uv", dict(base, sheets="twist", pairs=False, coulomb=False, channels=False)),
                    ("coulomb_only", dict(base, tail_tol=1e9, sheets=None, pairs=False, channels=False, uv=False)),
                    ("coulomb_plus_tails", dict(base, sheets=None, pairs=False, channels=False, uv=False)),
                    ("coulomb_plus_uv", dict(base, tail_tol=1e9, sheets=None, pairs=False, channels=False)),
                    ("coulomb_plus_tails_plus_uv", dict(base, sheets=None, pairs=False, channels=False)),
                    ("coulomb_plus_tails_plus_channels", dict(base, sheets=None, pairs=False, uv=False)),
                    ("coulomb_plus_tails_plus_uv_plus_channels", dict(base, sheets=None, pairs=False)),
                    ("coulomb_plus_tails_plus_uv_plus_channels_plus_pairs", dict(base, sheets=None)),
                    ("all_but_zigzag", dict(base, sheets="twist")),
                    ("all_but_twist", dict(base, sheets="zigzag")),
                    ("no_tails_no_sheets", dict(base, tail_tol=1e9, sheets=None)),
                    ("full_A1_rows", dict(base))):
        r = run_lp(TWO_DERIV, o, "l1", verbose=False, max_rounds=150)
        rec = {"status": r["status"], "norm": r.get("norm_l1"), "x": r.get("x"), "cuts": r["cuts"],
               "worst_rel_eig": r.get("worst_rel_eig"), "max_violation": r.get("max_violation")}
        if r["status"] == "infeasible":
            c = farkas(r["_rows"], TWO_DERIV, o.get("coulomb", True))
            rec["certificate"] = {k: c[k] for k in ("status", "beta", "r_inf", "margin", "support")}
        sub[name] = rec
        log(f"  tol1e-6 subset {name}: {r['status']} norm {r.get('norm_l1')} worst {r.get('worst_rel_eig')} "
            f"x {({k: round(v, 4) for k, v in r['x'].items() if abs(v) > 1e-4} if 'x' in r else '')} "
            f"{('cert ' + str([(l, round(w, 3)) for l, w in rec['certificate']['support'][:5]])) if 'certificate' in rec else ''}")
    PSD_TOL = 1e-7
    R["subsets_tol1e-6"] = sub
    return R


# ============================================================ main
def main():
    stages = ARGV[1:] or ["all"]
    if "all" in stages:
        stages = ["lp", "uvq", "tail", "extra", "extra2", "extra3"]
    out = {}
    if os.path.exists(OUT_JSON):
        try:
            out = json.load(open(OUT_JSON))
        except Exception:
            out = {}
    for s in stages:
        log(f"=== stage {s}")
        if s == "lp":
            out["lp"] = stage_lp()
        elif s == "uvq":
            out["uvq"] = stage_uvq()
        elif s == "tail":
            out["tail"] = stage_tail()
        elif s == "extra":
            out["extra"] = stage_extra()
        elif s == "extra2":
            out["extra2"] = stage_extra2()
        elif s == "extra3":
            out["extra3"] = stage_extra3()
        out["_meta"] = {"date": time.strftime("%Y-%m-%d %H:%M"), "elapsed_s": time.time() - T0,
                        "rows_json": os.path.basename(ROWS_JSON), "wall_s_unused": ROWS.get("wall_s")}
        json.dump(out, open(OUT_JSON, "w"), indent=1, default=float)
        log(f"wrote {OUT_JSON}")


if __name__ == "__main__":
    main()
