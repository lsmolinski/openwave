"""M5.32 R14-D adversarial audit: the P250 bridge on the 4x4 field.

Independent refutation attempt of the eight R14-D claims (D1..D8) with the auditor's
own sympy for V4 and iota on diagonal states, own grid scans with an own strict
8-neighbor local-minimum criterion, own continuous multistart minimization on an
EXTENDED box (the producer's [-1, 1.5]^2 box is a hypothesis, not a fact), own exact
asymptotics along the split, own 1D and 2D order-of-transition analysis (the inf of
r(s) = DV'(s) / (c iota(s)) decides the order), a mutation per gate, and a kink-tension
integral checked with two normalizations.  The producer's script
(m5_32_r14_d_bridge.py) and its result JSON were never opened.

EQUATIONS FIRST (all per unit volume, g = 8, delta = 3/10, W1 = 0.000724023879)
  N = M eta = diag(-m0, m1, m2, m3), the uniform states here diag(g, 1, m2, m3)
  V4(m2, m3) = W1 sum_{p=1..4} (m2^p + m3^p - delta^p)^2     [the (-g)^p + 1 cancel]
  f(x) = (x + g)(x - 1),  iota(m2, m3) = [f(m2) f(m3)]^2 (m2 - m3)^2
  V_eff = V4 - c omega^2 iota,   E_J = V4 + J^2 / (4 c iota),   omega = J / (2 c iota)
  split: m2 = delta/2 + s/2, m3 = delta/2 - s/2;  D': V' = V4 + mu s^2
  order of the ticking transition: omega_onset^2 = inf_{s>0} r(s),
     r(s) = [V'(s) - V'(0)] / (c iota(s)),  r(0+) = mu_eff / (c F0) = omega_inst^2,
     F0 = f(delta/2)^4, mu_eff = mu + V4''(0)/2;  the inf at s -> 0+ is a continuous
     (second-order) onset, an interior inf is a first-order crossing with a barrier
  kink tension: E[s(x)] = int [c kappa(s) s'^2 + DV(s)] dx,  kappa = (f(m2)^4 + f(m3)^4)/8
     Bogomolny:  sigma = int_0^{s*} 2 sqrt(c kappa DV) ds   (the brief writes
     sqrt(2 kappa DV): that is the same number only if kappa carries a 1/4, see report)
  thin wall: R = 2 sigma / p, p = V_eff(0) - V_eff(s*)

Usage: python3 m5_32_r14_d_audit.py
Writes data/m5_32_r14_d_audit.json.
"""

import sys
ARGV = list(sys.argv)                     # captured BEFORE any import

import importlib.util
import json
import os
import time

import numpy as np
import sympy as sp
from scipy.integrate import quad
from scipy.optimize import brentq, minimize, minimize_scalar

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.dirname(HERE)
DATA = os.path.join(RES, "data")
OUT_JSON = os.path.join(DATA, "m5_32_r14_d_audit.json")
T0 = time.time()

G = 8.0
DELTA = 0.3
W1 = 0.000724023879                       # own constant; cross-checked against B3.W1 below


def log(m):
    print(f"[{time.time() - T0:7.1f}s] {m}", flush=True)


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


# ============================================================ own sympy
g_s, d_s = sp.Rational(8), sp.Rational(3, 10)
x2, x3, s_s, t_s = sp.symbols("m2 m3 s t", real=True)


def f_sym(x):
    return (x + g_s) * (x - 1)


# exact rational polynomials; W1 multiplies P4 at evaluation time
P4u = sum((((-g_s) ** p + 1 + x2 ** p + x3 ** p)
           - ((-g_s) ** p + 1 + d_s ** p)) ** 2 for p in range(1, 5))   # unexpanded (numerics)
IOTAu = (f_sym(x2) * f_sym(x3)) ** 2 * (x2 - x3) ** 2
P4 = sp.expand(P4u)                                                     # expanded (degrees)
IOTA = sp.expand(IOTAu)
SPL = {x2: d_s / 2 + s_s / 2, x3: d_s / 2 - s_s / 2}
P4su, IOTAsu = P4u.subs(SPL), IOTAu.subs(SPL)
P4s = sp.expand(P4su)
IOTAs = sp.expand(IOTAsu)
KAPPAs = (f_sym(x2) ** 4 + f_sym(x3) ** 4).subs(SPL) / 8

_P4n = sp.lambdify((x2, x3), P4u, "numpy")
_IOn = sp.lambdify((x2, x3), IOTAu, "numpy")
_P4g = [sp.lambdify((x2, x3), sp.diff(P4u, v), "numpy") for v in (x2, x3)]
_IOg = [sp.lambdify((x2, x3), sp.diff(IOTAu, v), "numpy") for v in (x2, x3)]
_P4H = [[sp.lambdify((x2, x3), sp.diff(P4u, u, v), "numpy") for v in (x2, x3)] for u in (x2, x3)]
_IOH = [[sp.lambdify((x2, x3), sp.diff(IOTAu, u, v), "numpy") for v in (x2, x3)] for u in (x2, x3)]
_P4sn = sp.lambdify(s_s, P4su, "numpy")
_IOsn = sp.lambdify(s_s, IOTAsu, "numpy")
_KAPn = sp.lambdify(s_s, KAPPAs, "numpy")
np.seterr(all="ignore")


def V4n(a, b):
    return W1 * _P4n(a, b)


def iota_n(a, b):
    return _IOn(a, b)


def f_n(x):
    return (x + G) * (x - 1.0)


def veff(a, b, om, c=1.0):
    return V4n(a, b) - c * om ** 2 * iota_n(a, b)


def veff_grad(a, b, om, c=1.0):
    return np.array([W1 * _P4g[0](a, b) - c * om ** 2 * _IOg[0](a, b),
                     W1 * _P4g[1](a, b) - c * om ** 2 * _IOg[1](a, b)])


def veff_hess(a, b, om, c=1.0):
    return np.array([[W1 * _P4H[i][j](a, b) - c * om ** 2 * _IOH[i][j](a, b)
                      for j in range(2)] for i in range(2)])


def V4s(s):
    return W1 * _P4sn(s)


def iota_s(s):
    return _IOsn(s)


# ============================================================ helpers
def grid_local_minima(Z):
    """strict local minima against the 8 neighbors, interior points only."""
    C = Z[1:-1, 1:-1]
    ok = np.ones_like(C, dtype=bool)
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            if di == 0 and dj == 0:
                continue
            ok &= C < Z[1 + di:Z.shape[0] - 1 + di, 1 + dj:Z.shape[1] - 1 + dj]
    ii, jj = np.nonzero(ok)
    return ii + 1, jj + 1


def continuous_minima(fun, jac, hess, starts, box, gtol=1e-11, dedupe=1e-4):
    """multistart BFGS on fun/W1-scaled objective; keep converged interior points with a
    positive-definite Hessian; dedupe."""
    found = []
    for x0 in starts:
        r = minimize(fun, x0, jac=jac, method="BFGS",
                     options={"gtol": gtol, "maxiter": 2000})
        x = r.x
        if not (box[0] < x[0] < box[1] and box[0] < x[1] < box[1]):
            continue
        if np.linalg.norm(jac(x)) > 1e-8:
            continue
        ev = np.linalg.eigvalsh(hess(x))
        if ev.min() <= 0:
            continue
        if all(np.linalg.norm(x - y) > dedupe for y, _ in found):
            found.append((x.copy(), float(fun(x))))
    found.sort(key=lambda p: p[1])
    return found


def veff_min_set(om, c, box, nstart=25, gtol=1e-11):
    fun = lambda x: veff(x[0], x[1], om, c) / W1
    jac = lambda x: veff_grad(x[0], x[1], om, c) / W1
    hess = lambda x: veff_hess(x[0], x[1], om, c) / W1
    ax = np.linspace(box[0] + 0.05, box[1] - 0.05, nstart)
    starts = [np.array([a, b]) for a in ax for b in ax]
    return continuous_minima(fun, jac, hess, starts, box, gtol=gtol)


def scan_local_minima_1d(F, s_lo, s_hi, n=24001):
    """local minima of F on (s_lo, s_hi]: dense scan + golden refinement."""
    ss = np.linspace(s_lo, s_hi, n)
    Fv = F(ss)
    out = []
    for k in range(1, n - 1):
        if Fv[k] < Fv[k - 1] and Fv[k] < Fv[k + 1]:
            r = minimize_scalar(F, bounds=(ss[k - 1], ss[k + 1]), method="bounded",
                                options={"xatol": 1e-12})
            out.append((float(r.x), float(r.fun)))
    return out



def newton_min(x0, om, c=1.0, step_cap=0.05):
    """polish a local minimum of V_eff at fixed omega by damped, step-capped Newton (strictly
    local: never leaves the basin of x0).  Convergence = Newton decrement below 1e-7 (scale
    free), the line search tolerates roundoff.  Returns (x, min Hessian eigenvalue / W1) or
    None if the point is not a strict local minimum."""
    fun = lambda x: veff(x[0], x[1], om, c) / W1
    jac = lambda x: veff_grad(x[0], x[1], om, c) / W1
    hess = lambda x: veff_hess(x[0], x[1], om, c) / W1
    x = np.asarray(x0, float).copy()
    dec = np.inf
    for _ in range(300):
        g, H = jac(x), hess(x)
        if not np.all(np.isfinite(H)):
            return None
        ev, U = np.linalg.eigh(H)
        emax = max(1.0, np.abs(ev).max())
        evm = np.where(ev > 1e-12 * emax, ev, np.abs(ev) + 1e-6 * emax)
        dx = U @ ((U.T @ g) / evm)
        dec = float(np.sqrt(abs(g @ dx)))
        if dec < 1e-7 and ev.min() > 0:
            break
        n = np.linalg.norm(dx)
        if n > step_cap:
            dx *= step_cap / n
        f0 = fun(x)
        t = 1.0
        while t > 1e-4 and not (fun(x - t * dx) <= f0 + 1e-12 * abs(f0)):
            t *= 0.5
        x = x - t * dx
    g, H = jac(x), hess(x)
    ev = np.linalg.eigvalsh(H)
    if not (dec < 1e-6) or ev.min() <= 0 or np.linalg.norm(x) > 1e3:
        return None
    return x, float(ev.min())


def continue_branch(x_seed, om_seed, direction, om_lo=1e-7, om_hi=1.0, nsteps=600):
    """follow a minimum family in omega from (x_seed, om_seed); direction +1 up, -1 down.
    Returns the list of points and the refined end (fold) omega if the branch is lost."""
    oms = np.geomspace(om_seed, om_hi if direction > 0 else om_lo, nsteps)
    x = np.asarray(x_seed, float)
    pts, lost_at = [], None
    for om in oms:
        r = newton_min(x, om)
        if r is None or np.linalg.norm(r[0] - x) > 0.25:
            lost_at = float(om)
            break
        x, evmin = r
        pts.append({"omega": float(om), "m2": float(x[0]), "m3": float(x[1]),
                    "Veff": float(veff(x[0], x[1], om)), "hess_min_eig_over_W1": evmin})
    fold = None
    if lost_at is not None and pts:
        lo, hi = pts[-1]["omega"], lost_at
        xl = np.array([pts[-1]["m2"], pts[-1]["m3"]])
        for _ in range(40):
            mid = np.sqrt(lo * hi)
            r = newton_min(xl, mid)
            if r is None or np.linalg.norm(r[0] - xl) > 0.25:
                hi = mid
            else:
                lo, xl = mid, r[0]
        fold = {"omega": float(lo), "m2": float(xl[0]), "m3": float(xl[1]),
                "hess_min_eig_over_W1": float(newton_min(xl, lo)[1])}
    return pts, fold


# ============================================================ D1
def d1_formulas(B3, R14T):
    log("D1: registry K_P_h and e_parts V4 on uniform lattice states vs own sympy")
    res = {"cases": [], "e_u_max": 0.0, "W1_stack": float(B3.W1), "W1_own": W1,
           "W1_match": bool(abs(B3.W1 - W1) < 1e-15)}
    cfg = B3.base_cfg(s=-1, n=6, L=9.0)      # s = -1: the code branch with M00 = +g
    vol = cfg["n"] ** 3 * cfg["h"] ** 3
    G1 = np.zeros((4, 4))
    G1[2, 3], G1[3, 2] = -1.0, 1.0
    worst = 0.0
    for (a, b) in [(0.3, 0.0), (0.5, -0.2), (0.15, 0.15), (-0.7, 1.3)]:
        Md = np.diag([G, 1.0, a, b])
        M = np.broadcast_to(Md, (cfg["n"],) * 3 + (4, 4)).copy()
        a0 = G1 @ M - M @ G1
        assert abs(a0[0, 0, 0, 2, 3] - (a - b)) < 1e-14
        kin = R14T.kin_energy("K_P_h", M, a0, cfg) / vol
        e_u, e_v = B3.e_parts(M, cfg)
        e_v = e_v / vol
        io, v4 = iota_n(a, b), V4n(a, b)
        # own closed forms, no sympy: direct arithmetic
        io_direct = (f_n(a) * f_n(b)) ** 2 * (a - b) ** 2
        v4_direct = W1 * sum((a ** p + b ** p - DELTA ** p) ** 2 for p in range(1, 5))
        d_io = abs(kin - io) / max(abs(io), 1.0)
        d_v4 = abs(e_v - v4) / max(abs(v4), 1.0)
        worst = max(worst, d_io, d_v4)
        # mutation: literal roots (+g, 1) instead of (-g, 1)
        kin_lit = float(np.sum(R14T.kp_h_kin(a0, M, R14T.ROOTS_LITERAL))) / cfg["n"] ** 3
        res["cases"].append({"m2": a, "m3": b, "kin_registry": kin, "iota_sympy": io,
                             "iota_direct": io_direct, "V4_eparts": e_v, "V4_sympy": v4,
                             "V4_direct": v4_direct, "rel_iota": d_io, "rel_V4": d_v4,
                             "kin_literal_roots_mutant": kin_lit,
                             "mutant_rel": abs(kin_lit - io) / max(io, 1e-30)})
        res["e_u_max"] = max(res["e_u_max"], float(abs(e_u)))
        log(f"   ({a:5.2f},{b:5.2f}) kin={kin:.9g} iota={io:.9g} rel={d_io:.1e} | "
            f"V4={e_v:.9g} own={v4:.9g} rel={d_v4:.1e} | literal-root mutant kin={kin_lit:.6g}")
    # sign-convention check: the code branch s = +1 with M00 = -g gives the SAME V4
    cfg_p = B3.base_cfg(s=1, n=4, L=6.0)
    Md = np.diag([-G, 1.0, 0.5, -0.2])
    M = np.broadcast_to(Md, (4, 4, 4, 4, 4)).copy()
    e_v_p = B3.e_parts(M, cfg_p)[1] / (4 ** 3 * cfg_p["h"] ** 3)
    res["sign_branch_check"] = {"V4_s_plus_1_M00_minus_g": float(e_v_p), "V4_own": V4n(0.5, -0.2),
                                "rel": abs(e_v_p - V4n(0.5, -0.2)) / V4n(0.5, -0.2)}
    # wrong-branch mutation: M00 = +g fed to the s = +1 cfg (C_p of the other sign)
    Md = np.diag([G, 1.0, 0.5, -0.2])
    M = np.broadcast_to(Md, (4, 4, 4, 4, 4)).copy()
    e_v_w = B3.e_parts(M, cfg_p)[1] / (4 ** 3 * cfg_p["h"] ** 3)
    res["sign_branch_mutant"] = {"V4_wrong_branch": float(e_v_w),
                                 "rel_vs_own": abs(e_v_w - V4n(0.5, -0.2)) / V4n(0.5, -0.2)}
    res["worst_rel"] = worst
    res["verdict"] = "CONFIRMED" if worst < 1e-9 and res["e_u_max"] < 1e-12 else "REFUTED"
    log(f"   worst rel = {worst:.2e}; e_u max = {res['e_u_max']:.1e}; "
        f"wrong-sign-branch V4 rel = {res['sign_branch_mutant']['rel_vs_own']:.3g}  -> {res['verdict']}")
    return res


# ============================================================ D2
def d2_ticks():
    log("D2: the vacuum ticks and the degenerate state")
    io_vac = IOTA.subs({x2: d_s, x3: 0})
    p4_deg = P4.subs({x2: d_s / 2, x3: d_s / 2})
    io_deg = IOTA.subs({x2: d_s / 2, x3: d_s / 2})
    v4_deg = W1 * float(p4_deg)
    cross = v4_deg / float(io_vac)
    res = {"iota_vac_exact": str(io_vac), "iota_vac": float(io_vac),
           "V4_deg_over_W1_exact": str(p4_deg), "V4_deg": v4_deg, "iota_deg": float(io_deg),
           "V4_vac": V4n(DELTA, 0.0), "omega2_cross_c1": cross, "omega_cross_c1": cross ** 0.5,
           "omega2_cross_c3": cross / 3.0,
           "mutant_delta_0.31": {"iota_vac": iota_n(0.31, 0.0),
                                 "note": "iota at (0.31, 0) is not 194.435: the number is delta-specific"}}
    ok = (abs(res["iota_vac"] - 194.435) < 5e-4 and abs(v4_deg / 1.799e-6 - 1) < 1e-3
          and res["iota_deg"] == 0 and abs(cross / 9.25e-9 - 1) < 2e-3)
    res["verdict"] = "CONFIRMED" if ok else "REFUTED"
    log(f"   iota(delta,0) = {io_vac} = {float(io_vac):.6f}; V4(deg) = {v4_deg:.4e}; "
        f"iota(deg) = {float(io_deg)}; omega^2_cross = {cross:.4e} -> {res['verdict']}")
    return res


# ============================================================ D3
def d3_plane():
    log("D3: local minima of V_eff on the (m2, m3) plane, grid + extended box + multistart")
    res = {"grid": {}, "extended": {}, "continuous": {}, "omega_ladder": []}
    ax = np.linspace(-1.0, 1.5, 251)
    A, B = np.meshgrid(ax, ax, indexing="ij")
    axe = np.linspace(-3.0, 3.5, 651)
    Ae, Be = np.meshgrid(axe, axe, indexing="ij")
    P4A, IOA = _P4n(A, B), _IOn(A, B)
    P4E, IOE = _P4n(Ae, Be), _IOn(Ae, Be)
    omegas = [0.0, 1e-5, 3e-5, 1e-4, 2e-4, 3e-4, 5e-4, 1e-3, 3e-3, 1e-2, 1e-1]
    for om in omegas:
        Z = W1 * P4A - om ** 2 * IOA
        ii, jj = grid_local_minima(Z)
        mins = sorted([(round(float(ax[i]), 3), round(float(ax[j]), 3), float(Z[i, j]))
                       for i, j in zip(ii, jj)])
        Ze = W1 * P4E - om ** 2 * IOE
        ie, je = grid_local_minima(Ze)
        mins_e = sorted([(round(float(axe[i]), 3), round(float(axe[j]), 3), float(Ze[i, j]))
                         for i, j in zip(ie, je)])
        cm = veff_min_set(om, 1.0, (-3.0, 3.5))
        cont = [(round(float(x[0]), 5), round(float(x[1]), 5), float(v * W1)) for x, v in cm]
        res["grid"][str(om)] = mins
        res["extended"][str(om)] = mins_e
        res["continuous"][str(om)] = cont
        log(f"   omega={om:<7g} grid251: {len(mins)} min {mins[:4]} | ext651: {len(mins_e)} "
            f"{mins_e[:4]} | continuous: {len(cont)} {cont[:4]}")
    # Newton continuation of the minimum families in omega
    ptsA, foldA = continue_branch([0.0, 0.3], 1e-6, +1)
    res["family_A_vacuum"] = {"points": ptsA[::15], "fold": foldA, "n": len(ptsA),
                              "last": ptsA[-1] if ptsA else None}
    exits = [b for b in ptsA if not (-1 <= b["m2"] <= 1.5 and -1 <= b["m3"] <= 1.5)]
    res["omega_branch_leaves_producer_box"] = exits[0]["omega"] if exits else None
    res["omega_saddle_node"] = foldA["omega"] if foldA else None
    log(f"   family A (the vacuum): tracked over {len(ptsA)} steps to omega = {ptsA[-1]['omega']:.4e} at "
        f"({ptsA[-1]['m2']:.3f},{ptsA[-1]['m3']:.3f}); leaves the producer box at omega = "
        f"{res['omega_branch_leaves_producer_box']}; fold: {foldA}")
    # family B: seeded from the multistart hit near (-2.34, 3.08) at omega = 1.58e-3, both directions
    seedB = newton_min([-2.34, 3.08], 1.58e-3)
    if seedB is not None:
        upB, _ = continue_branch(seedB[0], 1.58e-3, +1)
        dnB, foldB = continue_branch(seedB[0], 1.58e-3, -1)
        res["family_B"] = {"seed": [float(seedB[0][0]), float(seedB[0][1])], "up": upB[::30], "down": dnB[::15],
                           "birth_fold": foldB, "last_up": upB[-1] if upB else None,
                           "lowest_omega_reached": dnB[-1]["omega"] if dnB else None}
        log(f"   family B: seed {seedB[0]}; reaches down to omega = {res['family_B']['lowest_omega_reached']}, "
            f"birth fold {foldB}; up to {upB[-1] if upB else None}")
    else:
        res["family_B"] = None
        log("   family B: the seed did not polish to a minimum")
    # family A': the (-1.61, -0.18) hit at 1.58e-3: same branch as A or a different one?
    seedC = newton_min([-1.61, -0.18], 1.58e-3)
    if seedC is not None:
        dnC, foldC = continue_branch(seedC[0], 1.58e-3, -1)
        upC, _ = continue_branch(seedC[0], 1.58e-3, +1)
        res["family_C"] = {"seed": [float(seedC[0][0]), float(seedC[0][1])], "down": dnC[::15], "birth_fold": foldC,
                           "lowest_omega_reached": dnC[-1]["omega"] if dnC else None, "last_up": upC[-1] if upC else None}
        log(f"   family C (m2<0, m3~0 at 1.58e-3): reaches down to omega = {res['family_C']['lowest_omega_reached']}, "
            f"birth fold {foldC}; up to {upC[-1] if upC else None}")
    res["omega_ladder"] = ptsA[::40]
    branch = ptsA
    # the count of NON-vacuum minima over the ladder, on the extended box, continuous
    other = 0
    for om in np.geomspace(1e-6, 1e-1, 26):
        cm = veff_min_set(om, 1.0, (-6.0, 6.0), nstart=13)
        for x, v in cm:
            onb = any(abs(b["omega"] / om - 1) < 0.2 and
                      min(np.hypot(x[0] - b["m2"], x[1] - b["m3"]), np.hypot(x[1] - b["m2"], x[0] - b["m3"])) < 0.35
                      for b in branch)
            if not onb:
                other += 1
                log(f"   NON-VACUUM MINIMUM at omega={om:.3g}: {x}, Veff={v * W1:.3e}")
    res["non_vacuum_minima_count"] = other
    # verdict
    g0 = res["grid"]["0.0"]
    g4 = res["grid"]["0.0001"]
    ok0 = len(g0) == 2 and any(abs(p[0] - 0.3) < 0.011 and abs(p[1]) < 0.011 for p in g0)
    ok4 = len(g4) == 2 and any(abs(p[0] + 0.04) < 0.011 and abs(p[1] - 0.33) < 0.011 for p in g4)
    ok3 = all(len(res["grid"][str(om)]) == 0 for om in (1e-3, 3e-3, 1e-2, 1e-1))
    okc = all(len(res["continuous"][str(om)]) == 0 for om in (1e-3, 3e-3, 1e-2, 1e-1))
    hidden = {str(om): res["extended"][str(om)] for om in (1e-3, 3e-3, 1e-2, 1e-1) if res["extended"][str(om)]}
    fb = res.get("family_B")
    fa_end = res["omega_saddle_node"] if res["omega_saddle_node"] else 1.0
    fb_birth = fb["birth_fold"]["omega"] if (fb and fb["birth_fold"]) else (fb["lowest_omega_reached"] if fb else None)
    res["second_family_coexists_with_vacuum"] = bool(fb_birth is not None and fb_birth < fa_end)
    res["checks"] = {"omega0_two_minima_at_vacuum": ok0, "omega1e-4_at_(-0.04,0.33)": ok4,
                     "omega>=1e-3_no_grid_minimum_LITERAL": ok3,
                     "omega>=1e-3_no_minimum_in_extended_box": okc,
                     "no_non_vacuum_minimum": other == 0}
    res["box_hides_deformed_vacuum"] = hidden
    res["verdict"] = ("REFUTED" if (not (ok0 and ok4 and ok3) or other > 0) else
                      ("QUALIFIED" if hidden else "CONFIRMED"))
    # mutation: c < 0 (a wrong-sign inertia) creates a NEW minimum away from the vacuum? and
    # a smaller box hides the disappearance: report the grid count at omega=5e-4 on [-0.2,0.6]
    Zm = W1 * P4A + (3e-4) ** 2 * IOA
    im, jm = grid_local_minima(Zm)
    res["mutant_c_minus1_omega3e-4_grid_minima"] = sorted([(round(float(ax[i]), 3), round(float(ax[j]), 3))
                                                           for i, j in zip(im, jm)])
    log(f"   checks {res['checks']} -> {res['verdict']}; c=-1 mutant minima: "
        f"{res['mutant_c_minus1_omega3e-4_grid_minima']}")
    return res


# ============================================================ D4
def d4_unbounded():
    log("D4: unboundedness of V_eff along the split, exact asymptotics")
    pi_ = sp.Poly(IOTAs, s_s)
    pv = sp.Poly(P4s, s_s)
    res = {"deg_iota": int(pi_.degree()), "lc_iota": str(pi_.LC()), "deg_P4": int(pv.degree()),
           "lc_V4_over_W1": str(pv.LC()), "ratio_leading": f"s^2 / (4 W1) = {1 / (4 * W1):.6g} s^2",
           "rows": []}
    for s in (1.0, 2.0, 4.0, 8.0, 16.0, 64.0):
        io, v4 = iota_s(s), V4s(s)
        res["rows"].append({"s": s, "iota": io, "V4": v4, "ratio": io / v4,
                            "ratio_leading_only": s ** 2 / (4 * W1)})
        log(f"   s={s:4g}: iota={io:.6g} V4={v4:.6g} iota/V4={io / v4:.4g} (leading-only {s ** 2 / (4 * W1):.4g})")
    v8 = V4s(8.0) - 1e-4 * iota_s(8.0)
    res["Veff_s8_omega1e-2_c1"] = float(v8)
    # first zero crossing of V_eff along the split for omega = 1e-2 and 1e-3
    # the runaway threshold along the split: the LAST local maximum of V_eff(s) (the top of the
    # last barrier before the monotone descent to -infinity) and its height above the ticking vacuum
    for om in (1e-2, 1e-3, 1e-4, 1e-5):
        F = lambda s: V4s(s) - om ** 2 * iota_s(s)
        ss = np.geomspace(0.05, 1e4, 400000)
        v = F(ss)
        k = np.nonzero((v[1:-1] > v[:-2]) & (v[1:-1] > v[2:]))[0]
        if len(k):
            kl = k[-1] + 1
            rr = minimize_scalar(lambda s: -F(s), bounds=(ss[kl - 1], ss[kl + 1]), method="bounded")
            s_top, v_top = float(rr.x), float(-rr.fun)
            vac = min(F(np.linspace(0.05, 0.6, 5001)))
            res[f"runaway_last_barrier_omega{om:g}"] = {"s_top": s_top, "Veff_top": v_top,
                                                        "height_above_ticking_vacuum": v_top - float(vac)}
        else:
            res[f"runaway_last_barrier_omega{om:g}"] = None
    res["iota_zeros_on_split"] = [float(z) for z in sorted(sp.solve(sp.Eq(f_sym(x2).subs(SPL) * f_sym(x3).subs(SPL), 0), s_s))]
    # generic rays m3 = k m2: degree of iota in t is 10 unless k in {0, 1}
    rays = {}
    for k in (sp.Rational(-1), sp.Rational(1, 2), sp.Rational(-3), sp.Rational(0), sp.Rational(1)):
        pi_k = sp.Poly(sp.expand(IOTA.subs({x2: t_s, x3: k * t_s})), t_s)
        pv_k = sp.Poly(sp.expand(P4.subs({x2: t_s, x3: k * t_s})), t_s)
        rays[str(k)] = {"deg_iota": str(pi_k.degree()), "deg_P4": str(pv_k.degree()),
                        "unbounded": bool(pi_k.degree() > pv_k.degree())}
    res["rays"] = rays
    # mutations: omega = 0 or c = -1 restore boundedness (V4 >= 0); a linear f (degree-6 iota) too
    f_lin = sp.expand(((x2 - 1) * (x3 - 1)) ** 2 * (x2 - x3) ** 2)
    res["mutant_linear_f_deg"] = int(sp.Poly(sp.expand(f_lin.subs(SPL)), s_s).degree())
    res["mutant_c_minus1_min_on_split"] = float(min(V4s(s) + 1e-4 * iota_s(s) for s in np.linspace(-10, 10, 20001)))
    ok = (res["deg_iota"] == 10 and res["deg_P4"] == 8 and abs(v8 / -3.6e3 - 1) < 0.03
          and abs(res["rows"][0]["ratio"] / 4.8e6 - 1) < 0.02 and abs(res["rows"][1]["ratio"] / 1.9e5 - 1) < 0.03
          and abs(res["rows"][2]["ratio"] / 7.9e5 - 1) < 0.02 and abs(res["rows"][3]["ratio"] / 2.0e5 - 1) < 0.03)
    res["verdict"] = "CONFIRMED" if ok else "REFUTED"
    log(f"   deg iota={res['deg_iota']} LC={res['lc_iota']}; deg V4={res['deg_P4']} LC=W1*{res['lc_V4_over_W1']}; "
        f"Veff(8,1e-2)={v8:.4g}; rays {rays} -> {res['verdict']}")
    return res


# ============================================================ D5
def d5_fixed_j():
    log("D5: fixed-J minima, grid + continuous + extended box")
    res = {"rows": []}
    ax = np.linspace(-1.0, 1.5, 251)
    A, B = np.meshgrid(ax, ax, indexing="ij")
    P4A, IOA = _P4n(A, B), _IOn(A, B)
    c = 1.0
    for J in (1e-3, 1e-2, 1e-1, 1.0, 10.0, 100.0):
        with np.errstate(divide="ignore", invalid="ignore"):
            E = W1 * P4A + J ** 2 / (4 * c * IOA)
        E[~np.isfinite(E)] = np.inf
        E[IOA < 1e-14] = np.inf
        k = np.unravel_index(np.argmin(E), E.shape)
        gm = (float(ax[k[0]]), float(ax[k[1]]), float(E[k]))
        # continuous from the grid minimum
        fun = lambda x: (V4n(x[0], x[1]) + J ** 2 / (4 * c * iota_n(x[0], x[1]))) / W1
        r = minimize(fun, np.array(gm[:2]), method="Nelder-Mead",
                     options={"xatol": 1e-10, "fatol": 1e-14, "maxiter": 20000})
        xc = r.x
        # extended-box multistart (Nelder-Mead; E_J is singular on the diagonal)
        best = (None, np.inf)
        allmin = []
        axs = np.linspace(-4.9, 5.4, 22)
        for a in axs:
            for b in axs:
                if abs(a - b) < 0.2:
                    continue
                rr = minimize(fun, np.array([a, b]), method="Nelder-Mead",
                              options={"xatol": 1e-8, "fatol": 1e-12, "maxiter": 4000})
                if not np.isfinite(rr.fun) or np.linalg.norm(rr.x) > 8:
                    continue
                if rr.fun < best[1]:
                    best = (rr.x.copy(), float(rr.fun))
                xs = rr.x if rr.x[0] <= rr.x[1] else rr.x[::-1]     # fold the mirror
                if all(np.linalg.norm(xs - y) > 1e-3 for y, _ in allmin):
                    allmin.append((xs.copy(), float(rr.fun)))
        allmin.sort(key=lambda q: q[1])
        def ej_hess_pd(x, h=1e-4):
            H = np.zeros((2, 2))
            for i in range(2):
                for j in range(2):
                    e_i, e_j = np.eye(2)[i] * h, np.eye(2)[j] * h
                    H[i, j] = (fun(x + e_i + e_j) - fun(x + e_i - e_j) - fun(x - e_i + e_j) + fun(x - e_i - e_j)) / (4 * h * h)
            return bool(np.linalg.eigvalsh(H).min() > 0)
        allmin = [(q[0], q[1], ej_hess_pd(q[0])) for q in allmin]
        om = J / (2 * c * iota_n(*xc))
        row = {"J": J, "grid_min": gm[:2], "grid_E": gm[2], "cont_min": [float(xc[0]), float(xc[1])],
               "cont_E": float(r.fun * W1), "omega": float(om),
               "ext_global_min": [float(best[0][0]), float(best[0][1])], "ext_global_E": best[1] * W1,
               "ext_agrees_with_grid": bool(min(np.linalg.norm(best[0] - xc),
                                                np.linalg.norm(best[0][::-1] - xc)) < 1e-3),
               "all_local_minima_ext_box_mirror_folded": [[float(q[0][0]), float(q[0][1]), q[1] * W1, q[2]] for q in allmin[:8]]}
        res["rows"].append(row)
        log(f"   J={J:<6g} grid {gm[:2]} -> cont ({xc[0]:.4f},{xc[1]:.4f}) E={r.fun * W1:.4e} "
            f"omega={om:.3e} | extended global ({best[0][0]:.4f},{best[0][1]:.4f}) E={best[1] * W1:.4e} "
            f"agrees={row['ext_agrees_with_grid']} | all minima (mirror folded): "
            f"{[(round(q[0][0], 3), round(q[0][1], 3), q[2]) for q in allmin[:6]]}")
    rows = {r_["J"]: r_ for r_ in res["rows"]}
    near = lambda p, q, tol=0.011: (abs(p[0] - q[0]) < tol and abs(p[1] - q[1]) < tol) or \
        (abs(p[0] - q[1]) < tol and abs(p[1] - q[0]) < tol)
    checks = {"J1e-3_at_vac": near(rows[1e-3]["grid_min"], (0.0, 0.3)),
              "J1e-2_at_vac": near(rows[1e-2]["grid_min"], (0.0, 0.3)),
              "J1e-1_at_(-0.07,0.35)": near(rows[1e-1]["grid_min"], (-0.07, 0.35)),
              "J1_at_(-0.36,0.41)": near(rows[1.0]["grid_min"], (-0.36, 0.41)),
              "omegas": [abs(rows[J]["omega"] / w - 1) < 0.05 for J, w in
                         ((1e-3, 2.6e-6), (1e-2, 2.6e-5), (1e-1, 1.3e-4), (1.0, 3.2e-4))],
              "extended_box_agrees_J<=1": all(rows[J]["ext_agrees_with_grid"] for J in (1e-3, 1e-2, 1e-1, 1.0))}
    res["checks"] = checks
    # boundedness: E_J >= V4 >= 0 and E_J -> +inf at iota -> 0 and at |m| -> inf (deg 8 vs J^2/iota -> 0)
    res["bounded_below"] = True
    # mutation: c = -1 gives E_J = V4 - J^2/(4 iota), unbounded near the diagonal
    res["mutant_c_minus1_E_at_(0.15,0.1500001)"] = float(V4n(0.15, 0.1500001) - 1.0 / (4 * iota_n(0.15, 0.1500001)))
    ok = all(v if not isinstance(v, list) else all(v) for v in checks.values())
    res["verdict"] = "CONFIRMED" if ok else "REFUTED"
    log(f"   checks {checks} -> {res['verdict']}")
    return res


# ============================================================ D6
def d6_dprime():
    log("D6: D' split stiffness at the degenerate point")
    d2 = sp.diff(P4s, s_s, 2).subs(s_s, 0)
    v4pp = W1 * float(d2)
    mu_loc = -v4pp / 2
    res = {"V4pp_over_W1_exact": str(d2), "V4pp_0": v4pp, "mu_local_threshold": mu_loc, "rows": []}
    ss = np.linspace(0, 1.2, 12001)
    v4 = V4s(ss)
    # the mu below which s = 0 stops being the GLOBAL min on [0, 1.2]
    with np.errstate(divide="ignore", invalid="ignore"):
        q = (v4[0] - v4[1:]) / ss[1:] ** 2
    kq = np.argmax(q)
    res["mu_global_threshold"] = float(q[kq])
    res["mu_global_threshold_at_s"] = float(ss[1 + kq])
    for mu in (1e-4, 1e-3, 1e-2, 5e-5, 3e-5):
        Vp = v4 + mu * ss ** 2
        k = np.argmin(Vp)
        loc = scan_local_minima_1d(lambda s: V4s(s) + mu * s ** 2, 0.0, 1.2, 12001)
        row = {"mu": mu, "Vp_0": float(Vp[0]), "argmin_s": float(ss[k]), "min_Vp": float(Vp[k]),
               "s0_is_global_min": bool(k == 0), "Vpp_0": v4pp + 2 * mu,
               "interior_local_minima": loc}
        res["rows"].append(row)
        log(f"   mu={mu:g}: V''(0)={v4pp + 2 * mu:+.3e}, argmin on [0,1.2] s={ss[k]:.3f}, "
            f"s=0 global: {k == 0}; interior local minima {loc}")
    ok = (abs(v4pp / -7.9e-5 - 1) < 0.02 and abs(mu_loc / 4e-5 - 1) < 0.02
          and all(r_["s0_is_global_min"] for r_ in res["rows"][:3]))
    res["verdict"] = "CONFIRMED" if ok else "REFUTED"
    log(f"   V4''(0) = W1 * {d2} = {v4pp:.4e}; mu_loc = {mu_loc:.4e}; mu_glob = "
        f"{res['mu_global_threshold']:.4e} at s = {res['mu_global_threshold_at_s']:.3f} -> {res['verdict']}")
    return res


# ============================================================ D7
def order_analysis(Vp, dVp_num, c, s_hi, tag):
    """Vp(s): the static potential along the split (numpy-vectorized).  Returns the two
    thresholds, the inf of r(s) on (0, s_hi], and a ladder tracking the split minimum."""
    F0 = f_n(DELTA / 2) ** 4
    ss = np.linspace(1e-4, s_hi, 60001)
    dV = Vp(ss) - Vp(0.0)
    r = dV / (c * iota_s(ss))
    k = int(np.argmin(r))
    # refine the interior minimum of r if it is interior
    rf = lambda s: (Vp(s) - Vp(0.0)) / (c * iota_s(s))
    if 0 < k < len(ss) - 1:
        rr = minimize_scalar(rf, bounds=(ss[k - 1], ss[k + 1]), method="bounded", options={"xatol": 1e-12})
        s_arg, r_min = float(rr.x), float(rr.fun)
    else:
        s_arg, r_min = float(ss[k]), float(r[k])
    # r(0+) by finite second derivatives: (Vp''(0)/2) / (c F0)
    h = 1e-4
    vpp0 = (Vp(h) - 2 * Vp(0.0) + Vp(-h)) / h ** 2
    r0 = (vpp0 / 2) / (c * F0)
    out = {"tag": tag, "c": c, "s_hi": s_hi, "omega2_inst_r0plus": float(r0),
           "r_inf": r_min, "r_argmin": s_arg, "r_at_small_s_1e-4": float(r[0]),
           "onset_is_interior": bool(s_arg > 5e-3 and r_min < r0 * (1 - 1e-6)),
           "ladder": []}
    om2_on = min(r_min, r0)
    out["omega2_onset"] = float(om2_on)
    for fac in (0.5, 0.8, 0.9, 0.95, 0.99, 1.0, 1.01, 1.05, 1.1, 1.2, 1.5, 2.0, 4.0):
        om2 = om2_on * fac
        Ve = lambda s: Vp(s) - c * om2 * iota_s(s)
        loc = scan_local_minima_1d(Ve, 0.0, s_hi, 24001)
        row = {"fac": fac, "omega2": float(om2), "minima": []}
        for smin, vmin in loc:
            sb = np.linspace(0, smin, 4001)
            bar = float(np.max(Ve(sb) - Ve(0.0)))
            row["minima"].append({"s": smin, "dV": float(vmin - Ve(0.0)), "barrier": bar})
        out["ladder"].append(row)
    return out


def sigma_of(Vp, c, om2, s_star):
    """kink tension between s = 0 and s = s_star at omega^2 = om2, two normalizations."""
    dVe = lambda s: max(Vp(s) - c * om2 * iota_s(s) - (Vp(0.0)), 0.0)
    bog = quad(lambda s: 2.0 * np.sqrt(c * _KAPn(s) * dVe(s)), 0.0, s_star, limit=200)[0]
    brief = quad(lambda s: np.sqrt(2.0 * _KAPn(s) * dVe(s)), 0.0, s_star, limit=200)[0]
    return float(bog), float(brief)


def d7_order():
    log("D7: order of the ticking transition in D' (1D split), all nine (mu, c)")
    F0 = f_n(DELTA / 2) ** 4
    v4pp = W1 * float(sp.diff(P4s, s_s, 2).subs(s_s, 0))
    res = {"F0": F0, "V4pp_0": v4pp, "rows": []}
    io_vac = iota_n(DELTA, 0.0)
    for mu in (1e-4, 1e-3, 1e-2):
        for c in (0.3, 1.0, 3.0):
            mu_eff = mu + v4pp / 2
            om2_inst = mu_eff / (c * F0)
            om2_cross = ((0.0 + mu * DELTA ** 2) - V4s(0.0)) / (c * io_vac)
            Vp = lambda s, mu=mu: V4s(s) + mu * s ** 2
            oa = order_analysis(Vp, None, c, 1.2, f"mu={mu:g},c={c:g}")
            oa17 = order_analysis(Vp, None, c, 1.69, f"mu={mu:g},c={c:g},s_hi=1.69")
            # the first omega^2 (above onset) where a split minimum lies below the exterior
            first = next((row for row in oa["ladder"] if any(m["dV"] < 0 for m in row["minima"])), None)
            row = {"mu": mu, "c": c, "mu_eff": mu_eff, "omega2_inst": om2_inst, "omega2_cross_delta": om2_cross,
                   "inst_over_cross": om2_inst / om2_cross, "percent_earlier": 100 * (1 - om2_inst / om2_cross),
                   "r_inf_on_(0,1.2]": oa["r_inf"], "r_argmin": oa["r_argmin"], "r0plus_fd": oa["omega2_inst_r0plus"],
                   "onset_interior_1.2": oa["onset_is_interior"], "onset_interior_1.69": oa17["onset_is_interior"],
                   "r_inf_on_(0,1.69]": oa17["r_inf"], "r_argmin_1.69": oa17["r_argmin"],
                   "ladder": oa["ladder"],
                   "first_dip": first}
            res["rows"].append(row)
            fd = "none" if first is None else f"fac={first['fac']} s={first['minima'][0]['s']:.4f} barrier={first['minima'][0]['barrier']:.2e}"
            log(f"   mu={mu:g} c={c:g}: inst={om2_inst:.3e} cross(delta)={om2_cross:.3e} "
                f"({row['percent_earlier']:.1f}% earlier); inf r on (0,1.2] = {oa['r_inf']:.3e} at s={oa['r_argmin']:.4f} "
                f"(r(0+)={oa['omega2_inst_r0plus']:.3e}); interior onset: {oa['onset_is_interior']}; first dip {fd}")
    # the runaway window: for mu = 1e-3, c = 1, where does r(s) fall below r(0+) again at large s?
    mu, c = 1e-3, 1.0
    Vp = lambda s: V4s(s) + mu * s ** 2
    r0 = (mu + v4pp / 2) / (c * F0)
    ss = np.geomspace(16.4, 1e4, 200000)
    rr = (Vp(ss) - Vp(0.0)) / (c * iota_s(ss))
    below = np.nonzero(rr < r0)[0]
    res["runaway_s_where_r_below_r0plus_mu1e-3_c1"] = float(ss[below[0]]) if len(below) else None
    res["runaway_s_leading_estimate"] = float(np.sqrt(4 * W1 / r0))
    # a mu ladder: is the onset ever interior with a UNIFORM quadratic penalty?
    lad = []
    for mu in np.geomspace(4.1e-5, 1e2, 30):
        Vp = lambda s, mu=mu: V4s(s) + mu * s ** 2
        oa = order_analysis(Vp, None, 1.0, 1.2, "")
        lad.append({"mu": float(mu), "onset_interior": oa["onset_is_interior"], "r_argmin": oa["r_argmin"],
                    "r_inf_over_r0plus": oa["r_inf"] / oa["omega2_inst_r0plus"]})
    res["mu_ladder_uniform_penalty"] = lad
    res["uniform_penalty_ever_first_order"] = any(l_["onset_interior"] for l_ in lad)
    # the delta-well comparison in closed form (why a uniform stiffness cannot do it)
    res["uniform_penalty_algebra"] = ("r(delta)/r(0+) = [(mu delta^2 - V4(0)) / iota(delta,0)] / [(mu + V4''(0)/2) / F0]; "
                                      "with F0/iota(delta,0) = %.5f and V4(0)/delta^2 = %.4e < -V4''(0)/2 = %.4e the ratio "
                                      "exceeds 1 for every mu > -V4''(0)/2" % (F0 / io_vac, V4s(0.0) / DELTA ** 2, -v4pp / 2))
    # ---- MUTANT: a LOCALIZED split stiffness at the degenerate point
    mu_l, w, nu, c = 1e-3, 0.1, 1e-4, 1.0
    Vm = lambda s: V4s(s) + mu_l * s ** 2 * np.exp(-s ** 2 / w ** 2) + nu * s ** 2
    oam = order_analysis(Vm, None, c, 1.2, "mutant local stiffness")
    om2_c = oam["r_inf"]
    s_c = oam["r_argmin"]
    Ve_c = lambda s: Vm(s) - c * om2_c * iota_s(s)
    sb = np.linspace(0, s_c, 4001)
    barrier_c = float(np.max(Ve_c(sb) - Ve_c(0.0)))
    sig_bog, sig_brief = sigma_of(Vm, c, om2_c, s_c)
    om2_b = 1.5 * om2_c
    Ve_b = lambda s: Vm(s) - c * om2_b * iota_s(s)
    loc_b = scan_local_minima_1d(Ve_b, 0.0, 1.2, 24001)
    s_b = min(loc_b, key=lambda p: p[1])[0] if loc_b else None
    p_b = float(Ve_b(0.0) - Ve_b(s_b)) if s_b else None
    sig_b = sigma_of(Vm, c, om2_b, s_b)[0] if s_b else None
    R_b = 2 * sig_bog / p_b if s_b else None            # thin wall: sigma at degeneracy
    ss = np.linspace(0, 1.2, 12001)
    res["mutant"] = {"form": "V4 + mu_l s^2 exp(-s^2/w^2) + nu s^2", "mu_l": mu_l, "w": w, "nu": nu, "c": c,
                     "Vpp_0": float((Vm(1e-4) - 2 * Vm(0) + Vm(-1e-4)) / 1e-8),
                     "exterior_global_min_omega0_on_[0,1.2]": bool(np.argmin(Vm(ss)) == 0),
                     "omega2_inst_r0plus": oam["omega2_inst_r0plus"], "omega2_cross_true": om2_c, "s_star": s_c,
                     "onset_interior": oam["onset_is_interior"], "barrier_at_crossing": barrier_c,
                     "sigma_bogomolny": sig_bog, "sigma_brief_formula": sig_brief,
                     "sigma_ratio": sig_bog / sig_brief if sig_brief else None,
                     "bag_at_1.5x": {"omega2": om2_b, "s_star": s_b, "p": p_b, "sigma_clipped_at_1.5x": sig_b,
                                     "sigma_c": sig_bog, "R": R_b},
                     "ladder": oam["ladder"]}
    log(f"   MUTANT local stiffness: V''(0)={res['mutant']['Vpp_0']:+.3e}, inst r(0+)={oam['omega2_inst_r0plus']:.3e}, "
        f"true crossing {om2_c:.3e} at s*={s_c:.4f} (interior: {oam['onset_is_interior']}), barrier {barrier_c:.3e}, "
        f"sigma={sig_bog:.4e} (brief formula {sig_brief:.4e}); bag at 1.5x: p={p_b}, R={R_b}")
    ok = (all(not r_["onset_interior_1.2"] for r_ in res["rows"])
          and all(r_["first_dip"] is not None and r_["first_dip"]["minima"][0]["barrier"] < 1e-12 for r_ in res["rows"])
          and all(0.05 < r_["percent_earlier"] / 100 < 0.32 for r_ in res["rows"])
          and oam["onset_is_interior"] and barrier_c > 0 and sig_bog > 0)
    res["verdict"] = "CONFIRMED" if ok else "REFUTED"
    log(f"   -> {res['verdict']}")
    return res


# ============================================================ D7 at plane level (not claimed)
def bogo_sigma_path(P, om2, c, mu, ext_val, nseg=1500):
    """Bogomolny tension along a piecewise-linear path P (k x 2) in the (m2, m3) plane with the
    K_P^h static metric on diagonal states, (c/2)(f(m2)^4 dm2^2 + f(m3)^4 dm3^2), and
    DV = V' - c om2 iota - ext_val clipped at 0.  Returns (sigma, max DV on the path)."""
    tot, bar = 0.0, -np.inf
    for i in range(len(P) - 1):
        a_, b_ = np.asarray(P[i]), np.asarray(P[i + 1])
        lam = np.linspace(0, 1, nseg + 1)
        mid = 0.5 * (lam[1:] + lam[:-1])
        pts = a_[None, :] + mid[:, None] * (b_ - a_)[None, :]
        dm = (b_ - a_) / nseg
        dV = V4n(pts[:, 0], pts[:, 1]) + mu * (pts[:, 0] - pts[:, 1]) ** 2 - c * om2 * iota_n(pts[:, 0], pts[:, 1]) - ext_val
        bar = max(bar, float(dV.max()))
        ds2 = (c / 2) * (f_n(pts[:, 0]) ** 4 * dm[0] ** 2 + f_n(pts[:, 1]) ** 4 * dm[1] ** 2)
        tot += float(np.sum(2 * np.sqrt(np.maximum(dV, 0) * ds2)))
    return tot, bar


def plane_min(x0, om2, c, mu):
    """strict local minimum of V_eff' = V4 + mu (m2 - m3)^2 - c om2 iota near x0 (trust-exact +
    PD check); returns (x, eigs) or None."""
    J2 = np.array([[1.0, -1.0], [-1.0, 1.0]])
    fun = lambda x: (V4n(x[0], x[1]) + mu * (x[0] - x[1]) ** 2 - c * om2 * iota_n(x[0], x[1])) / W1
    jac = lambda x: (veff_grad(x[0], x[1], np.sqrt(om2), c) + 2 * mu * (x[0] - x[1]) * np.array([1.0, -1.0])) / W1
    hess = lambda x: (veff_hess(x[0], x[1], np.sqrt(om2), c) + 2 * mu * J2) / W1
    r = minimize(fun, np.asarray(x0, float), jac=jac, hess=hess, method="trust-exact",
                 options={"gtol": 1e-12, "maxiter": 500})
    ev = np.linalg.eigvalsh(hess(r.x))
    if np.linalg.norm(jac(r.x)) > 1e-7 * max(1.0, np.abs(hess(r.x)).max()) or ev.min() <= 0:
        return None
    return r.x, ev * W1


def plane_case(mu, c, t_star, ax, A, B, P4A, IOA):
    Vp = lambda a_, b_: V4n(a_, b_) + mu * (a_ - b_) ** 2
    ext = float(Vp(t_star, t_star))
    Fs = f_n(t_star) ** 4
    vss = (Vp(t_star + 5e-5, t_star - 5e-5) - 2 * ext + Vp(t_star - 5e-5, t_star + 5e-5)) / 1e-8
    om2_inst = (vss / 2) / (c * Fs)
    with np.errstate(divide="ignore", invalid="ignore"):
        R = (W1 * P4A + mu * (A - B) ** 2 - ext) / (c * IOA)
    R[np.abs(A - B) < 0.02] = np.inf
    k = np.unravel_index(np.argmin(R), R.shape)
    fun = lambda x: (Vp(x[0], x[1]) - ext) / (c * iota_n(x[0], x[1]))
    rr = minimize(fun, np.array([ax[k[0]], ax[k[1]]]), method="Nelder-Mead",
                  options={"xatol": 1e-10, "fatol": 1e-16, "maxiter": 20000})
    row = {"mu": mu, "c": c, "exterior": [t_star, t_star], "V_ext": ext, "omega2_inst_plane": float(om2_inst),
           "plane_inf_r": float(rr.fun), "plane_argmin": [float(rr.x[0]), float(rr.x[1])],
           "grid_argmin": [float(ax[k[0]]), float(ax[k[1]])], "grid_argmin_on_box_edge": bool(k[0] in (0, len(ax) - 1) or k[1] in (0, len(ax) - 1)),
           "onset_interior_plane": bool(rr.fun < om2_inst * (1 - 1e-4) and abs(rr.x[0] - rr.x[1]) > 0.02)}
    if row["onset_interior_plane"]:
        om2_c = float(rr.fun)
        ext_min = plane_min([t_star, t_star], om2_c, c, mu)
        int_min = plane_min(rr.x, om2_c, c, mu)
        row["exterior_is_local_min_at_crossing"] = ext_min is not None
        row["interior_is_local_min_at_crossing"] = int_min is not None
        if int_min is not None and ext_min is not None:
            xi = int_min[0]
            row["interior_state"] = [float(xi[0]), float(xi[1])]
            row["interior_split"] = float(xi[0] - xi[1])
            row["interior_sum_minus_delta"] = float(xi[0] + xi[1] - DELTA)
            row["interior_hess_eigs"] = [float(e) for e in int_min[1]]
            row["exterior_hess_eigs"] = [float(e) for e in ext_min[1]]
            row["omega_c"] = float(np.sqrt(om2_c))
            P0 = [np.array(ext_min[0]), xi]
            sig_line, bar_line = bogo_sigma_path(P0, om2_c, c, mu, ext)
            # 3-waypoint path optimization of the Bogomolny tension
            w0 = np.array([ext_min[0] + (xi - ext_min[0]) * q for q in (0.25, 0.5, 0.75)]).ravel()
            obj = lambda w: bogo_sigma_path([ext_min[0]] + list(w.reshape(3, 2)) + [xi], om2_c, c, mu, ext, nseg=400)[0]
            ro = minimize(obj, w0, method="Nelder-Mead", options={"xatol": 1e-6, "fatol": 1e-12, "maxiter": 6000})
            Popt = [ext_min[0]] + list(ro.x.reshape(3, 2)) + [xi]
            sig_opt, bar_opt = bogo_sigma_path(Popt, om2_c, c, mu, ext)
            row["sigma_straight_path"] = sig_line
            row["barrier_straight_path"] = bar_line
            row["sigma_optimized_path"] = sig_opt
            row["barrier_optimized_path"] = bar_opt
            row["optimized_waypoints"] = [[float(q[0]), float(q[1])] for q in Popt]
            # thin-wall bag R = 2 sigma_c / p with sigma_c the tension at degeneracy (the crossing) and
            # p the pressure at omega^2 = 1.5 x crossing and at 0.95 x the exterior instability (top of
            # the coexistence window)
            row["coexistence_window_omega2"] = [om2_c, float(om2_inst)]
            row["bags"] = {}
            for lab, om2_b in (("1.5x_crossing", 1.5 * om2_c), ("0.95x_inst", 0.95 * float(om2_inst))):
                ib = plane_min(xi, om2_b, c, mu)
                eb = plane_min(ext_min[0], om2_b, c, mu)
                if ib is not None and eb is not None:
                    Vi = V4n(*ib[0]) + mu * (ib[0][0] - ib[0][1]) ** 2 - c * om2_b * iota_n(*ib[0])
                    Ve = V4n(*eb[0]) + mu * (eb[0][0] - eb[0][1]) ** 2 - c * om2_b * iota_n(*eb[0])
                    p_b = float(Ve - Vi)
                    row["bags"][lab] = {"omega2": om2_b, "interior": [float(ib[0][0]), float(ib[0][1])], "p": p_b,
                                        "sigma_c": sig_opt, "R": 2 * sig_opt / p_b if p_b > 0 else None,
                                        "both_states_local_minima": True}
                else:
                    row["bags"][lab] = {"omega2": om2_b, "exterior_still_min": eb is not None, "interior_still_min": ib is not None}
    return row


def d7_plane(mu_list=(1e-4, 1e-3, 1e-2), c=1.0):
    log("D7-plane: the D' exterior in the full (m2, m3) plane (the degenerate point is not stationary)")
    res = {"rows": []}
    fd = lambda t: V4n(t, t)
    rt = minimize_scalar(fd, bounds=(-0.5, 1.0), method="bounded", options={"xatol": 1e-12})
    t_star = float(rt.x)
    gdeg = veff_grad(DELTA / 2, DELTA / 2, 0.0)
    res["grad_V4_at_degenerate"] = [float(gdeg[0]), float(gdeg[1])]
    res["t_star"] = t_star
    res["V4_at_t_star"] = float(rt.fun)
    res["V4_at_degenerate"] = V4n(DELTA / 2, DELTA / 2)
    log(f"   grad V4 at (delta/2, delta/2) = {gdeg}; diagonal minimum at t* = {t_star:.6f}, "
        f"V4 = {rt.fun:.4e} (vs {V4n(DELTA / 2, DELTA / 2):.4e} at delta/2)")
    ax = np.linspace(-1.5, 1.5, 601)
    A, B = np.meshgrid(ax, ax, indexing="ij")
    P4A, IOA = _P4n(A, B), _IOn(A, B)
    for mu in mu_list:
        for cc in ((0.3, 1.0, 3.0) if mu == 1e-3 else (c,)):
            row = plane_case(mu, cc, t_star, ax, A, B, P4A, IOA)
            res["rows"].append(row)
            extra = ""
            if row.get("interior_state"):
                extra = (f"; interior state {row['interior_state']} (split {row['interior_split']:.3f}, sum-delta "
                         f"{row['interior_sum_minus_delta']:.3f}), both local minima at omega_c={row['omega_c']:.3e}; "
                         f"barrier (straight/opt) {row['barrier_straight_path']:.2e}/{row['barrier_optimized_path']:.2e}; "
                         f"sigma (straight/opt) {row['sigma_straight_path']:.3e}/{row['sigma_optimized_path']:.3e}; "
                         f"bags {row.get('bags')}")
            log(f"   mu={mu:g} c={cc:g}: inst(plane)={row['omega2_inst_plane']:.3e}; inf r = {row['plane_inf_r']:.3e} at "
                f"{row['plane_argmin']} (grid on edge: {row['grid_argmin_on_box_edge']}); first order: "
                f"{row['onset_interior_plane']}{extra}")
    # the mu threshold of the plane-level first-order onset
    lad = []
    for mu in np.geomspace(1e-4, 1e-2, 17):
        row = plane_case(float(mu), 1.0, t_star, ax, A, B, P4A, IOA)
        lad.append({"mu": float(mu), "first_order": row["onset_interior_plane"], "inf_r_over_inst": row["plane_inf_r"] / row["omega2_inst_plane"],
                    "argmin": row["plane_argmin"]})
    res["mu_ladder"] = lad
    fo = [l_["mu"] for l_ in lad if l_["first_order"]]
    res["mu_first_order_threshold_between"] = [max([l_["mu"] for l_ in lad if not l_["first_order"]] or [None]), min(fo) if fo else None]
    log(f"   plane-level first order for mu in {[round(m_, 6) for m_ in fo]}; threshold between {res['mu_first_order_threshold_between']}")
    res["first_order_found"] = bool(any(r_.get("interior_state") and r_["barrier_optimized_path"] > 0 and r_["sigma_optimized_path"] > 0
                                        for r_ in res["rows"]))
    return res


# ============================================================ main
def main():
    B3 = _load("m5_21_3_a_4d", "m5_21_3_a_4d.py")
    R14T = _load("m5_32_r14_terms", "m5_32_r14_terms.py")
    out = {"argv": ARGV, "conventions": {"g": G, "delta": DELTA, "W1": W1,
                                         "V4": "W1 sum_p (m2^p + m3^p - delta^p)^2",
                                         "iota": "[f(m2) f(m3)]^2 (m2 - m3)^2, f = (x + g)(x - 1)",
                                         "P4_exact": str(P4), "IOTA_exact": str(IOTA)}}
    out["D1"] = d1_formulas(B3, R14T)
    out["D2"] = d2_ticks()
    out["D3"] = d3_plane()
    out["D4"] = d4_unbounded()
    out["D5"] = d5_fixed_j()
    out["D6"] = d6_dprime()
    out["D7"] = d7_order()
    out["D7_plane"] = d7_plane()
    if out["D7"]["verdict"] == "CONFIRMED" and out["D7_plane"]["first_order_found"]:
        out["D7"]["verdict"] = "QUALIFIED"
        out["D7"]["qualification"] = ("true on the 1D split line as defined (sum m2 + m3 frozen at delta, exterior at delta/2); "
                                      "in the (m2, m3) plane the exterior sits at t* and the crossing state relaxes the sum: "
                                      "first order with a barrier for mu = 1e-3 and 1e-2 (see D7_plane)")
    d8ok = (out["D7"]["verdict"] != "REFUTED" and not out["D7"]["uniform_penalty_ever_first_order"]
            and out["D7"]["mutant"]["onset_interior"] and not out["D7_plane"]["first_order_found"])
    out["D8"] = {"verdict": "CONFIRMED" if d8ok else "REFUTED",
                 "basis": ("REFUTED at plane level: with the (2,3) eigenvalue penalty mu >= ~1e-3 the exterior (t*, t*) and a "
                           "rotating interior off the split line are both local minima at the crossing omega with a barrier "
                           "and a nonzero Bogomolny tension (D7_plane rows); the 1D split-line statement holds"
                           if out["D7_plane"]["first_order_found"] else
                           "D7 second order for all nine (mu, c) and the whole uniform-mu ladder; the localized stiffness "
                           "mutant is first order with a barrier, a nonzero sigma and a finite bag radius"),
                 "plane_level": {k: v for k, v in out["D7_plane"].items() if k != "rows"}}
    out["runtime_s"] = time.time() - T0
    with open(OUT_JSON, "w") as fh:
        json.dump(out, fh, indent=1, default=float)
    log(f"wrote {OUT_JSON}")
    for k in ("D1", "D2", "D3", "D4", "D5", "D6", "D7", "D8"):
        log(f"   {k}: {out[k]['verdict']}")


if __name__ == "__main__":
    main()
