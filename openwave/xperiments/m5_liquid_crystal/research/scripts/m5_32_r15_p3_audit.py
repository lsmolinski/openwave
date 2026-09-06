"""M5.32 R15-P-iii ADVERSARIAL AUDIT: the reduced planar functional of
L_P on the degenerate vacuum (independent implementation; the producer's
scripts and data were NOT opened).

Own tools: sympy for the densities + exact gradients (lambdified), a
LINK-BASED trapezoid discretization of the 1D wall functional (fields on
nodes, derivatives on links, midpoint fields on links; different from the
producer's per-cell sym stencil), L-BFGS-B relaxation, the Beltrami first
integral as the discretization-independent check, and the certified 4D
stack (m5_21_3_a_4d.py: e_parts / kin_of / inner_eta / a_fields, G1 from
m5_21_8_b_lattice.py) for the slab check.

Claims audited (W1..W5) + the off-diagonal (2,3) mutation.

Modes: all | w1 | w2 | w3 | w4 | mut
Out: ../data/m5_32_r15_p3_audit.json
"""
from __future__ import annotations

import sys
ARGV = list(sys.argv)                      # captured before any import
import os
os.environ.setdefault("OMP_NUM_THREADS", "2")
import importlib.util
import json
import time

import numpy as np
import sympy as sp
from scipy.optimize import minimize

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
OUT = os.path.join(DATA, "m5_32_r15_p3_audit.json")

G, DELTA, W1 = 8.0, 0.3, 0.000724023879
MUS = [1e-3, 1e-2, 1e-1]
CS = [0.1, 1.0, 10.0]
RNG = np.random.default_rng(20260905)

# --------------- certified stack (allowed) ---------------
_spec = importlib.util.spec_from_file_location(
    "ins4", os.path.join(HERE, "m5_21_3_a_4d.py"))
INS4 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(INS4)
G1 = np.zeros((4, 4)); G1[2, 3] = -1.0; G1[3, 2] = 1.0   # = m5_21_8_b G1


# --------------- reduced densities (own) ---------------
def v4dd(m2, m3):
    v = 0.0
    for p in range(1, 5):
        v = v + (m2 ** p + m3 ** p - 2 * DELTA ** p) ** 2
    return W1 * v


def v4dd_grad(m2, m3):
    g2 = 0.0; g3 = 0.0
    for p in range(1, 5):
        r = m2 ** p + m3 ** p - 2 * DELTA ** p
        g2 = g2 + 2 * r * p * m2 ** (p - 1)
        g3 = g3 + 2 * r * p * m3 ** (p - 1)
    return W1 * g2, W1 * g3


def veff(m2, m3, eps):
    """V_eff = V4dd + (mu - omega^2 c) s^2 = V4dd - eps s^2."""
    return v4dd(m2, m3) - eps * (m2 - m3) ** 2


def veff_grad(m2, m3, eps):
    g2, g3 = v4dd_grad(m2, m3)
    s = m2 - m3
    return np.array([g2 - 2 * eps * s, g3 + 2 * eps * s])


def min_veff(eps, nstart=60, box=6.0):
    """multistart BFGS on the unrestricted plane; returns the lowest
    minimum and the list of distinct minima."""
    starts = RNG.uniform(-box, box, size=(nstart, 2))
    starts = np.vstack([starts, [[DELTA, DELTA], [2, -1], [-1, 2],
                                 [3, -2.5], [-2.5, 3]]])
    mins = []
    for x0 in starts:
        r = minimize(lambda x: veff(x[0], x[1], eps), x0,
                     jac=lambda x: veff_grad(x[0], x[1], eps),
                     method="BFGS", options={"gtol": 1e-12})
        mins.append((r.fun, r.x[0], r.x[1]))
    mins.sort()
    distinct = []
    for f, a, b in mins:
        if all(abs(a - d[1]) + abs(b - d[2]) > 1e-4 for d in distinct):
            distinct.append((float(f), float(a), float(b)))
    return distinct


# --------------- sympy: the 3-field (m2, m3, q) planar densities -------
def build_sym():
    m2, m3, q, d2, d3, dq = sp.symbols("m2 m3 q d2 d3 dq", real=True)
    c, mu, om2 = sp.symbols("c mu om2", positive=True)
    B = sp.Matrix([[m2, q], [q, m3]])
    Bp = sp.Matrix([[d2, dq], [dq, d3]])
    J = sp.Matrix([[0, -1], [1, 0]])                  # G1 restricted
    trp = [(B ** p).trace() for p in range(1, 5)]
    V4 = W1 * sum((trp[p - 1] - 2 * DELTA ** p) ** 2 for p in range(1, 5))
    S2_inv = (m2 - m3) ** 2 + 4 * q ** 2               # eigenvalue split^2
    S2_diag = (m2 - m3) ** 2                           # diagonal-entry split^2
    a0 = J * B - B * J                                 # (2,3) clock G1 M - M G1
    Fz = a0 * Bp - Bp * a0                             # comm_eta on the block
    I1 = 4 * sum(Fz[i, j] ** 2 for i in range(2) for j in range(2))
    KPin = sp.Rational(1, 2) * c * sum(a0[i, j] ** 2 for i in range(2)
                                       for j in range(2))   # = c (s^2 + 4q^2)
    KPst = sp.Rational(1, 2) * c * (d2 ** 2 + d3 ** 2 + 2 * dq ** 2)
    out = {}
    for tag, S2 in (("inv", S2_inv), ("diag", S2_diag)):
        Vnode = V4 + mu * S2 - om2 * KPin
        Llink = KPst - om2 * I1
        args_n = (m2, m3, q, c, mu, om2)
        args_l = (m2, m3, q, d2, d3, dq, c, om2)
        out[tag] = dict(
            V=sp.lambdify(args_n, Vnode, "numpy"),
            dV=sp.lambdify(args_n, [sp.diff(Vnode, v) for v in (m2, m3, q)],
                           "numpy"),
            L=sp.lambdify(args_l, Llink, "numpy"),
            dL=sp.lambdify(args_l, [sp.diff(Llink, v)
                                    for v in (m2, m3, q, d2, d3, dq)], "numpy"),
        )
    out["I1_expr"] = sp.simplify(I1)
    out["KPin_expr"] = sp.simplify(KPin)
    out["I1_diag"] = sp.simplify(I1.subs({q: 0, dq: 0}))
    out["V4_diag"] = sp.expand(V4.subs(q, 0))
    return out


SYM = None


def sym():
    global SYM
    if SYM is None:
        SYM = build_sym()
    return SYM


# --------------- the 1D wall functional (link discretization) --------
class Wall:
    """nodes 0..N (fields m2, m3, q), ends pinned; F = h [sum_i w_i V_i
    + sum_links L_link(mid, diff/h)]; free variables = interior nodes;
    q frozen at 0 when offdiag=False."""

    def __init__(self, tag, c, mu, om2, L, N, left, right, offdiag):
        self.f = sym()[tag]
        self.c, self.mu, self.om2 = c, mu, om2
        self.L, self.N, self.h = L, N, L / N
        self.left, self.right = np.array(left, float), np.array(right, float)
        self.offdiag = offdiag
        self.z = np.linspace(-L / 2, L / 2, N + 1)

    def full(self, x):
        nf = 3 if self.offdiag else 2
        Y = np.zeros((self.N + 1, 3))
        Y[0, :] = self.left; Y[-1, :] = self.right
        Y[1:-1, :nf] = x.reshape(self.N - 1, nf)
        return Y

    def energy_parts(self, Y):
        h = self.h
        Vn = self.f["V"](Y[:, 0], Y[:, 1], Y[:, 2], self.c, self.mu, self.om2)
        w = np.ones(self.N + 1); w[0] = w[-1] = 0.5
        mid = 0.5 * (Y[1:] + Y[:-1]); dif = (Y[1:] - Y[:-1]) / h
        Ll = self.f["L"](mid[:, 0], mid[:, 1], mid[:, 2],
                         dif[:, 0], dif[:, 1], dif[:, 2], self.c, self.om2)
        return h * np.sum(w * Vn), h * np.sum(Ll)

    def energy(self, Y):
        a, b = self.energy_parts(Y)
        return a + b

    def fun_grad(self, x):
        Y = self.full(x); h = self.h
        Vn = self.f["V"](Y[:, 0], Y[:, 1], Y[:, 2], self.c, self.mu, self.om2)
        dV = self.f["dV"](Y[:, 0], Y[:, 1], Y[:, 2], self.c, self.mu, self.om2)
        w = np.ones(self.N + 1); w[0] = w[-1] = 0.5
        mid = 0.5 * (Y[1:] + Y[:-1]); dif = (Y[1:] - Y[:-1]) / h
        la = (mid[:, 0], mid[:, 1], mid[:, 2], dif[:, 0], dif[:, 1], dif[:, 2],
              self.c, self.om2)
        Ll = self.f["L"](*la)
        dL = self.f["dL"](*la)
        F = h * (np.sum(w * Vn) + np.sum(Ll))
        Gn = np.zeros((self.N + 1, 3))
        for k in range(3):
            Gn[:, k] += h * w * np.broadcast_to(dV[k], Vn.shape)
            gm = 0.5 * h * np.broadcast_to(dL[k], Ll.shape)
            gd = np.broadcast_to(dL[3 + k], Ll.shape)      # h * (1/h)
            Gn[:-1, k] += gm - gd
            Gn[1:, k] += gm + gd
        nf = 3 if self.offdiag else 2
        return F, Gn[1:-1, :nf].ravel()

    def uniform_ref(self):
        """F of the uniform left state on the same line."""
        Vl = self.f["V"](*self.left, self.c, self.mu, self.om2)
        return self.L * float(Vl)

    def relax(self, Y0, maxiter=20000):
        nf = 3 if self.offdiag else 2
        x0 = Y0[1:-1, :nf].ravel()
        r = minimize(self.fun_grad, x0, jac=True, method="L-BFGS-B",
                     options={"maxiter": maxiter, "maxfun": 4 * maxiter,
                              "ftol": 1e-15, "gtol": 1e-10, "maxcor": 30})
        return self.full(r.x), r

    def beltrami(self, Y):
        """T_link - V_mid on each link (autonomous Lagrangian => const)."""
        mid = 0.5 * (Y[1:] + Y[:-1]); dif = (Y[1:] - Y[:-1]) / self.h
        T = self.f["L"](mid[:, 0], mid[:, 1], mid[:, 2],
                        dif[:, 0], dif[:, 1], dif[:, 2], self.c, self.om2)
        V = self.f["V"](mid[:, 0], mid[:, 1], mid[:, 2], self.c, self.mu,
                        self.om2)
        return T - V


def width_10_90(z, s):
    smax = s[-1]; smin = s[0]
    lo = smin + 0.1 * (smax - smin); hi = smin + 0.9 * (smax - smin)
    zl = np.interp(lo, s, z) if s[-1] > s[0] else np.interp(lo, s[::-1], z[::-1])
    zh = np.interp(hi, s, z) if s[-1] > s[0] else np.interp(hi, s[::-1], z[::-1])
    return abs(zh - zl)


# ============================ W1 ============================
def run_w1():
    t0 = time.time()
    res = {}
    # (a) analytic: V4dd = 0 => power sums 1..2 match => multiset {delta,delta}
    m2, m3 = sp.symbols("m2 m3", real=True)
    sols = sp.solve([m2 + m3 - 2 * DELTA, m2 ** 2 + m3 ** 2 - 2 * DELTA ** 2],
                    [m2, m3], dict=True)
    res["zero_set_p12"] = [{str(k): float(v) for k, v in s.items()} for s in sols]
    # quartic-flat coefficient along the split at fixed sum
    s = sp.symbols("s", real=True)
    V4s = sym()["V4_diag"].subs({sp.Symbol("m2", real=True): DELTA + s / 2,
                                 sp.Symbol("m3", real=True): DELTA - s / 2})
    ser = sp.Poly(sp.expand(V4s), s)
    res["V4_along_split_poly_coeffs_s^k"] = {
        str(k[0]): float(v) for k, v in zip(ser.monoms(), ser.coeffs())}
    # (b) wide random search for V4dd < 0 and for local minima of V4dd
    X = RNG.uniform(-12, 12, size=(400000, 2))
    vals = v4dd(X[:, 0], X[:, 1])
    res["random_min_V4dd_wide"] = float(vals.min())
    res["random_argmin"] = X[np.argmin(vals)].tolist()
    res["V4dd_local_minima_(multistart BFGS)"] = min_veff(0.0, nstart=300, box=8.0)[:6]
    # (c) the nine points at 0.999 omega_c^2: eps = mu (0.999 - 1) < 0
    pts = {}
    for mu in MUS:
        for c in CS:
            omc2 = mu / c
            row = {}
            for ratio in (0.999, 1.0, 1.001):
                om2 = ratio * omc2
                eps = om2 * c - mu
                d = min_veff(eps, nstart=80, box=6.0)
                row[str(ratio)] = {"min_Veff": d[0][0], "argmin": d[0][1:],
                                   "n_distinct_minima": len(d),
                                   "any_nonvacuum_min_le_0":
                                   bool(any(f <= 1e-15 and abs(a - DELTA)
                                            + abs(b - DELTA) > 1e-3
                                            for f, a, b in d))}
            # grid check on [-2,2]^2 at 0.999
            gx = np.linspace(-2, 2, 801)
            A, Bm = np.meshgrid(gx, gx, indexing="ij")
            eps = 0.999 * omc2 * c - mu
            Vg = veff(A, Bm, eps)
            row["grid_min_0.999"] = float(Vg.min())
            row["grid_argmin_0.999"] = [float(A.ravel()[Vg.argmin()]),
                                        float(Bm.ravel()[Vg.argmin()])]
            pts[f"mu={mu:g},c={c:g}"] = row
    res["nine_points"] = pts
    res["runtime_s"] = time.time() - t0
    return res


# ============================ W2 ============================
def run_w2():
    t0 = time.time()
    res = {}
    for mu in MUS:
        row = {}
        for ratio in (1.01, 2.0):
            eps = mu * (ratio - 1.0)          # (omega^2 c - mu) at omega^2 = ratio mu/c
            d = min_veff(eps, nstart=120, box=8.0)
            f, a, b = d[0]
            # mirror partner must be degenerate
            row[str(ratio)] = {"eps": eps, "Veff_min": f, "a": a, "b": b,
                               "split": abs(a - b), "sum": a + b,
                               "inside_box_2": bool(abs(a) <= 2 and abs(b) <= 2),
                               "n_distinct_minima": len(d),
                               "all_minima": d[:6]}
        res[f"mu={mu:g}"] = row
    # c-independence: V_eff only via eps = omega^2 c - mu; verify numerically
    chk = []
    for c in CS:
        mu = 1e-2; om2 = 2.0 * mu / c
        eps = om2 * c - mu
        chk.append(min_veff(eps, nstart=40)[0])
    res["c_independence_mu=1e-2_ratio2"] = chk
    res["runtime_s"] = time.time() - t0
    return res


# ============================ W3 ============================
def wall_point(mu, c, ratio=2.0, N=600, Lfac=30.0, tag="diag", offdiag=False,
               init="ising", maxiter=20000, Lover=None):
    om2 = ratio * mu / c
    eps = om2 * c - mu
    f, a, b = min_veff(eps, nstart=40)[0]
    ell = np.sqrt(c / mu)
    L = Lover if Lover is not None else Lfac * ell
    w = Wall(tag, c, mu, om2, L, N, (a, b, 0.0), (b, a, 0.0), offdiag)
    z = w.z
    s0 = a - b
    Y0 = np.zeros((N + 1, 3))
    if init == "ising":
        sz = -s0 * np.tanh(z / (0.5 * ell))
        Y0[:, 0] = (a + b) / 2 + sz / 2
        Y0[:, 1] = (a + b) / 2 - sz / 2
        if offdiag:
            Y0[:, 2] = 0.05 * s0 / np.cosh(z / (0.5 * ell))   # seed q
    elif init == "rot":
        th = 0.5 * np.pi * (0.5 + 0.5 * np.tanh(z / (0.2 * L)))  # 0 -> pi/2
        # B = R(th) diag(a,b) R(th)^T
        Y0[:, 0] = a * np.cos(th) ** 2 + b * np.sin(th) ** 2
        Y0[:, 1] = a * np.sin(th) ** 2 + b * np.cos(th) ** 2
        Y0[:, 2] = (a - b) * np.sin(th) * np.cos(th)
    Y0[0] = w.left; Y0[-1] = w.right
    F0 = w.energy(Y0)
    Y, r = w.relax(Y0, maxiter=maxiter)
    F = w.energy(Y)
    ref = w.uniform_ref()
    s = Y[:, 0] - Y[:, 1]
    S = np.sqrt(s ** 2 + 4 * Y[:, 2] ** 2)
    belt = w.beltrami(Y)
    out = {"mu": mu, "c": c, "ratio": ratio, "om2": om2, "eps": eps,
           "a": a, "b": b, "s_star": s0, "L": L, "N": N, "h": w.h,
           "F_init": F0, "F_final": F, "F_uniform": ref,
           "tension": F - ref, "tension_over_sqrt_c_mu": (F - ref) / np.sqrt(c * mu),
           "width_10_90_s": width_10_90(z, s), "width_over_ell": width_10_90(z, s) / ell,
           "max_abs_s": float(np.max(np.abs(s))), "max_S_inv": float(np.max(S)),
           "max_abs_q": float(np.max(np.abs(Y[:, 2]))),
           "coef_sprime2_at_sstar_(c/4-8om2 s*^2)": c / 4 - 8 * om2 * s0 ** 2,
           "producer_criterion_16om2s2<c": bool(16 * om2 * s0 ** 2 < c),
           "my_criterion_32om2s2<c": bool(32 * om2 * s0 ** 2 < c),
           "beltrami_spread": float(np.std(belt[5:-5])),
           "beltrami_mean": float(np.mean(belt[5:-5])),
           "beltrami_expected(-V_uniform)": float(-ref / L),
           "grad_norm_final": float(np.max(np.abs(r.jac))),
           "nit": int(r.nit), "success": bool(r.success), "msg": str(r.message)}
    return out, (z, Y)


def run_w3():
    t0 = time.time()
    res = {"well_defined": {}, "runaway": {}, "threshold_mutation": {}}
    # the 4 producer-well-defined points, two resolutions
    for mu, c in ((1e-3, 1.0), (1e-3, 10.0), (1e-2, 10.0), (1e-1, 10.0)):
        rows = {}
        for N in (300, 600):
            o, _ = wall_point(mu, c, N=N)
            rows[f"N={N}"] = o
        res["well_defined"][f"mu={mu:g},c={c:g}"] = rows
    # two runaway points, capped iterations: energy must fall below the uniform
    # reference by far more than any wall tension, |s| growing
    for mu, c in ((1e-2, 1.0), (1e-3, 0.1)):
        o, (z, Y) = wall_point(mu, c, N=300, maxiter=400)
        res["runaway"][f"mu={mu:g},c={c:g}"] = o
    # continuum runaway family at mu=1e-2, c=1: s = s* + A sin(k z) on a window
    mu, c = 1e-2, 1.0
    om2 = 2 * mu / c; eps = om2 * c - mu
    f, a, b = min_veff(eps, nstart=40)[0]
    s0 = a - b; u0 = a + b
    fam = []
    zz = np.linspace(-5, 5, 200001); dz = zz[1] - zz[0]
    for k in (1.0, 3.0, 10.0, 30.0, 100.0):
        A = 0.3 * s0
        s = s0 + A * np.sin(k * zz) * np.exp(-(zz / 2) ** 4)
        m2 = (u0 + s) / 2; m3 = (u0 - s) / 2
        ds = np.gradient(s, dz)
        dens = (c / 2) * (0.5 * ds ** 2) + veff(m2, m3, eps) - om2 * 8 * s ** 2 * ds ** 2
        Fk = np.sum(dens) * dz - veff(a, b, eps) * (zz[-1] - zz[0])
        fam.append({"k": k, "F_minus_uniform": float(Fk)})
    res["continuum_family_mu=1e-2,c=1"] = {"s_star": s0, "coef": c / 4 - 8 * om2 * s0 ** 2,
                                          "rows": fam}
    # threshold mutation: mu=1e-1, ratio 2, c chosen between the two criteria
    mu = 1e-1
    eps = mu
    f, a, b = min_veff(eps, nstart=40)[0]
    s0 = a - b
    # producer: 32 mu s*^2 < c^2 ; mine: 64 mu s*^2 < c^2
    c_lo = np.sqrt(32 * mu * s0 ** 2); c_hi = np.sqrt(64 * mu * s0 ** 2)
    c_mid = 0.5 * (c_lo + c_hi)
    res["threshold_mutation"]["s_star"] = s0
    res["threshold_mutation"]["c_producer_threshold"] = float(c_lo)
    res["threshold_mutation"]["c_my_threshold"] = float(c_hi)
    for c in (c_mid, 1.05 * c_hi):
        o, _ = wall_point(mu, c, N=300, maxiter=600)
        res["threshold_mutation"][f"c={c:.4f}"] = o
    res["runtime_s"] = time.time() - t0
    return res


# ============================ W4 ============================
def slab_check(nz=48, h=0.5, offdiag=False):
    """profile on a 4 x 4 x nz slab; certified stack vs own reduced sums."""
    cfg = INS4.base_cfg(n=nz, L=nz * h)
    assert abs(cfg["h"] - h) < 1e-15
    z = (np.arange(nz) - (nz - 1) / 2) * h
    s = 1.2 * np.tanh(z / 3.0) + 0.1 * np.sin(z)
    u = 2 * DELTA + 0.2 * np.exp(-(z / 4) ** 2)
    m2 = (u + s) / 2; m3 = (u - s) / 2
    q = 0.3 / np.cosh(z / 2.5) if offdiag else np.zeros(nz)
    M = np.zeros((4, 4, nz, 4, 4))
    M[..., 0, 0] = G; M[..., 1, 1] = 1.0
    M[..., 2, 2] = m2; M[..., 3, 3] = m3
    M[..., 2, 3] = q; M[..., 3, 2] = q
    c = 1.0
    h3 = h ** 3; nxy = 16
    # certified stack
    e_u, _ = INS4.e_parts(M, cfg)
    a0 = G1 @ M - M @ G1
    kin = INS4.kin_of(M, a0, cfg)
    # K_P^23 via the stack's inner_eta on the projected z-derivative
    P = np.diag([0.0, 0.0, 1.0, 1.0])
    kp_st_stack = 0.0
    for br, (A, wt) in INS4.a_fields(M, cfg).items():
        Om = P @ A[2] @ P @ INS4.ETA
        kp_st_stack += wt * (c / 2) * np.sum(INS4.inner_eta(Om, Om))
    kp_st_stack *= h3
    kp_in_stack = h3 * (c / 2) * np.sum(INS4.inner_eta(a0 @ INS4.ETA, a0 @ INS4.ETA))
    # V4dd via the stack's M eta powers? No: own route = eigenvalues of N
    N = M @ INS4.ETA
    lam = np.linalg.eigvals(N).real
    v_stack = 0.0
    for p in range(1, 5):
        Cp = (-G) ** p + 1 + 2 * DELTA ** p
        v_stack = v_stack + (np.sum(lam ** p, axis=-1) - Cp) ** 2
    v_eig = h3 * W1 * np.sum(v_stack)
    # own reduced sums (sym stencil = 1/2 fwd + 1/2 bwd per cell, one-sided
    # rows zero as in d1)
    def dfwd(f):
        o = np.zeros_like(f); o[:-1] = (f[1:] - f[:-1]) / h; return o

    def dbwd(f):
        o = np.zeros_like(f); o[1:] = (f[1:] - f[:-1]) / h; return o
    red = {}
    if not offdiag:
        red["V4dd"] = nxy * h3 * np.sum(v4dd(m2, m3))
        red["split"] = nxy * h3 * np.sum(s ** 2)
        red["KP_static"] = nxy * h3 * (c / 2) * sum(
            0.5 * np.sum(d(m2) ** 2 + d(m3) ** 2) for d in (dfwd, dbwd))
        red["KP_inertia"] = nxy * h3 * c * np.sum(s ** 2)
        red["I1_inertia"] = nxy * h3 * sum(
            0.5 * np.sum(8 * s ** 2 * d(s) ** 2) for d in (dfwd, dbwd))
    else:
        f = sym()["inv"]
        S2 = s ** 2 + 4 * q ** 2
        red["V4dd"] = nxy * h3 * np.sum(
            f["V"](m2, m3, q, c, 0.0, 0.0))            # mu=0, om2=0 => V4dd only
        red["split_inv"] = nxy * h3 * np.sum(S2)
        red["KP_static"] = nxy * h3 * (c / 2) * sum(
            0.5 * np.sum(d(m2) ** 2 + d(m3) ** 2 + 2 * d(q) ** 2) for d in (dfwd, dbwd))
        red["KP_inertia"] = nxy * h3 * c * np.sum(S2)
        I1 = sp.lambdify(sp.symbols("m2 m3 q d2 d3 dq", real=True), sym()["I1_expr"], "numpy")
        red["I1_inertia"] = nxy * h3 * sum(
            0.5 * np.sum(I1(m2, m3, q, d(m2), d(m3), d(q))) for d in (dfwd, dbwd))
    stack = {"E_u": float(e_u), "V4dd_eig": float(v_eig), "KP_static": float(kp_st_stack),
             "KP_inertia": float(kp_in_stack), "I1_inertia": float(kin)}
    rel = {"V4dd": abs(stack["V4dd_eig"] - red["V4dd"]) / red["V4dd"],
           "KP_static": abs(stack["KP_static"] - red["KP_static"]) / red["KP_static"],
           "KP_inertia": abs(stack["KP_inertia"] - red["KP_inertia"]) / red["KP_inertia"],
           "I1_inertia": abs(stack["I1_inertia"] - red["I1_inertia"]) / red["I1_inertia"]}
    return {"nz": nz, "h": h, "offdiag": offdiag, "stack": stack,
            "reduced": {k: float(v) for k, v in red.items()}, "rel_err": rel,
            "kappa_P_measured": float(kp_in_stack / (nxy * h3 * c * np.sum(
                s ** 2 + 4 * q ** 2)))}


def run_w4():
    t0 = time.time()
    res = {"diag_h0.5": slab_check(48, 0.5, False),
           "diag_h1.5": slab_check(32, 1.5, False),
           "offdiag_h0.5": slab_check(48, 0.5, True)}
    res["sympy_I1_diag"] = str(sym()["I1_diag"])
    res["sympy_KPin"] = str(sym()["KPin_expr"])
    res["runtime_s"] = time.time() - t0
    return res


# ============================ MUTATION ============================
def run_mut():
    t0 = time.time()
    res = {}
    # (1) uniform sector with q free, under both readings of the split term
    uni = {}
    for tag in ("inv", "diag"):
        f = sym()[tag]
        rows = {}
        for mu, c in ((1e-3, 1.0), (1e-1, 10.0)):
            for ratio in (0.5, 0.999, 2.0):
                om2 = ratio * mu / c
                best = None
                for x0 in RNG.uniform(-3, 3, size=(60, 3)):
                    r = minimize(lambda x: float(f["V"](x[0], x[1], x[2], c, mu, om2)),
                                 x0, jac=lambda x: np.array(f["dV"](x[0], x[1], x[2], c, mu, om2), float),
                                 method="BFGS", options={"gtol": 1e-12})
                    if best is None or r.fun < best[0]:
                        best = (float(r.fun), r.x.tolist())
                # diagonal-only minimum for comparison
                eps = om2 * c - mu
                d0 = min_veff(eps, nstart=30)[0]
                rows[f"mu={mu:g},c={c:g},ratio={ratio}"] = {
                    "min_with_q": best[0], "argmin_(m2,m3,q)": best[1],
                    "S_inv": float(np.sqrt((best[1][0] - best[1][1]) ** 2 + 4 * best[1][2] ** 2)),
                    "min_diag_only": d0[0], "lower_than_diag": bool(best[0] < d0[0] - 1e-12)}
        uni[tag] = rows
    res["uniform_with_q"] = uni
    # (2) wall with q free at two well-defined points, both readings, two inits
    walls = {}
    for mu, c in ((1e-3, 1.0), (1e-1, 10.0)):
        for tag in ("inv", "diag"):
            for init in ("ising", "rot"):
                o, (z, Y) = wall_point(mu, c, N=400, tag=tag, offdiag=True, init=init,
                                       maxiter=6000)
                # analytic rotation-wall estimates
                ell = np.sqrt(c / mu); s0 = o["s_star"]
                o["rot_wall_tension_diag_reading_2s*^2sqrt(c mu)"] = 2 * s0 ** 2 * np.sqrt(c * mu)
                o["rot_twist_tension_inv_reading_c s*^2 (pi/2)^2/L"] = c * s0 ** 2 * (np.pi / 2) ** 2 / o["L"]
                walls[f"mu={mu:g},c={c:g},{tag},{init}"] = o
        # inv reading: tension must scale ~ 1/L (Goldstone twist), check 2 L
        for Lf in (15.0, 60.0):
            o, _ = wall_point(mu, c, N=400, tag="inv", offdiag=True, init="rot",
                              maxiter=6000, Lfac=Lf)
            walls[f"mu={mu:g},c={c:g},inv,rot,Lfac={Lf}"] = {
                "L": o["L"], "tension": o["tension"],
                "c s*^2 (pi/2)^2 / L": c * o["s_star"] ** 2 * (np.pi / 2) ** 2 / o["L"]}
    res["walls_with_q"] = walls
    # (3) diagonal Ising wall reference at the same points (q frozen)
    ref = {}
    for mu, c in ((1e-3, 1.0), (1e-1, 10.0)):
        o, _ = wall_point(mu, c, N=400)
        ref[f"mu={mu:g},c={c:g}"] = {"tension": o["tension"], "width_over_ell": o["width_over_ell"]}
    res["ising_reference"] = ref
    res["runtime_s"] = time.time() - t0
    return res


# ============================ main ============================
def main():
    mode = ARGV[1] if len(ARGV) > 1 else "all"
    res = {}
    if os.path.exists(OUT):
        with open(OUT) as fh:
            res = json.load(fh)
    todo = ["w1", "w2", "w3", "w4", "mut"] if mode == "all" else [mode]
    for m in todo:
        t0 = time.time()
        res[m] = {"w1": run_w1, "w2": run_w2, "w3": run_w3, "w4": run_w4,
                  "mut": run_mut}[m]()
        print(f"[{m}] done in {time.time() - t0:.1f}s", flush=True)
        with open(OUT, "w") as fh:
            json.dump(res, fh, indent=1, default=float)
    print(json.dumps(res[todo[-1]], indent=1, default=float)[:6000])


if __name__ == "__main__":
    main()
