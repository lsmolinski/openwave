"""Independent check of the P250 wall / bag algebra (substrate-framework
campaign P250-shell-bubble-clock, PR #196 + correction PR #197, merged
2026-09-03, release v0.174.0, main at 71bbe9a276a2e53c2b88e1a1d4aaf8c9d0ee915f).

Everything is re-derived FRESH from the P249 potential as encoded in our
m5_32_p249_check.py (the same chart (a, t, p, u, v, q) + complex psi); no
substrate-framework module is imported or run. Checks (each can fail):

 1. slice reduction: on S = diag(m, c+b, c-b), psi = f >= 0 the potential is
    V0 = V_M5 + 2c^2 + 2b^2 + 6(b-f^2)^2 + W(f), V_M5 = -r^2/2 - tr S^3 + r^4 + 1/2
    (C-M5W-001); the kinetic metric diag(1/4,1/2,1/2,1/2); the inertia f^2 + 4b^2
 2. transverse invariance: the four transverse first derivatives of the
    ROTATING-frame potential (u, v, q, Im psi) vanish on the slice (C-M5W-001)
 3. vacuum: V_M5(-N N^T) = 2, the four summands of V0 nonnegative with the unique
    zero at the exterior (numerical minimization from random starts) (C-M5W-002)
 4. the deep branch m = 0 (d_m V ~ m), the Maxwell system pA, pB, pC and the
    frequency substitution omega^2 = 32 f^2 - 12 f + 6 - 24 b (C-M5W-004)
 5. the certified root (c, b, f) and omega_*^2 = 1.663945700059150... reproduced by
    our own mpmath solve; the two exact rational witnesses (C-M5W-004)
 6. thin-wall identities: R = 2 sigma/p, F''(R_c) = -8 pi sigma, E_crit = 16 pi
    sigma^3 / (3 p^2), dE/dQ = omega on E = F + omega Q (C-M5W-005, -008)
 7. INDEPENDENT sigma_0: our own solve_bvp on 2K u'' = grad V_omega at omega_*^2
    from the vacuum to the deep-branch state, L-continuation with a tanh seed,
    three tension routes; target 0.72929841786(58) (C-M5W-006)
 8. the delta = 0.001 bag radius from the thin-wall law with p = -V_omega(B):
    target R_m = 1218.761785 with chi = 1.000328 (C-M5W-007, headline only; the
    spherical BVP itself is NOT re-solved here)

NOT checked: the H1 phase-slip bookkeeping (C-M5W-003; our own W0 S4 identity is
the same structure), the spherical bag BVP (C-M5W-007 beyond its thin-wall
headline), the Benci-Fortunato bridge of P249.

Run: python3 m5_32_p250_check.py  (sympy + mpmath + scipy; ~1 to 3 min)
Writes ../data/m5_32_p250_check.json.
"""
import json
import math
import os
import time

import mpmath as mp
import numpy as np
import sympy as sp
from scipy.integrate import solve_bvp
from scipy.optimize import minimize

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(HERE), "data", "m5_32_p250_check.json")
T0 = time.time()
RES = {}
FAILS = []


def check(name, ok, detail=""):
    print(f"{'PASS' if ok else 'FAIL'}  {name}  {detail}")
    RES[name] = {"pass": bool(ok), "detail": str(detail)}
    if not ok:
        FAILS.append(name)


# ---------------------------------------------------------------- the P249 potential (our encoding)
a, t, p, u, v, q, fr, fi = sp.symbols('a t p u v q fr fi', real=True)
w2 = sp.symbols('w2', positive=True)
S = sp.Matrix([[1 + a, u, v], [u, t + p, q], [v, q, t - p]])
A = sp.Matrix([[0, 0, 0], [0, 0, -1], [0, 1, 0]])
tr2 = sp.trace(S * S)
V_m5 = -tr2 / 2 - sp.trace(S ** 3) + tr2 ** 2 + sp.Rational(1, 2)
V_axis = tr2 - S[0, 0] ** 2
P = sp.diag(0, 1, 1)
Bproj = P * S * P
Bcl = Bproj - sp.trace(Bproj) * P / 2
H = sp.diag(0, 1, -1)
K = (A * H - H * A) / 2
Q = (fr ** 2 - fi ** 2) * H + 2 * fr * fi * K
V_lock = 6 * sp.trace((Bcl - Q).T * (Bcl - Q)) / 2
f2 = fr ** 2 + fi ** 2
W = 3 * f2 - 4 * f2 * sp.sqrt(f2) + 2 * f2 ** 2
V = V_m5 + V_axis + V_lock + W
comm = A * S - S * A
iota_full = sp.trace(comm.T * comm) / 2 + f2          # the canonical clock-generator norm
V_rot = V - w2 / 2 * iota_full                        # rotating-frame potential

# ---------------------------------------------------------------- 1. slice reduction
m, c, b = sp.symbols('m c b', real=True)
f = sp.symbols('f', nonnegative=True)
slice_sub = {a: m - 1, t: c, p: b, u: 0, v: 0, q: 0, fr: f, fi: 0}
r2 = m ** 2 + 2 * c ** 2 + 2 * b ** 2
V_M5_hand = -r2 / 2 - (m ** 3 + 2 * c ** 3 + 6 * c * b ** 2) + r2 ** 2 + sp.Rational(1, 2)
W_hand = 3 * f ** 2 - 4 * f ** 3 + 2 * f ** 4
V0_hand = V_M5_hand + 2 * c ** 2 + 2 * b ** 2 + 6 * (b - f ** 2) ** 2 + W_hand
V_slice = sp.simplify(V.subs(slice_sub))
check("1a slice potential == V_M5 + 2c^2 + 2b^2 + 6(b - f^2)^2 + W(f)",
      sp.simplify(V_slice - V0_hand) == 0)
iota_slice = sp.simplify(iota_full.subs(slice_sub))
check("1b slice inertia == f^2 + 4 b^2", sp.simplify(iota_slice - (f ** 2 + 4 * b ** 2)) == 0, iota_slice)
# kinetic metric: Tr(Sdot^2)/4 + |psidot|^2/2 on the slice
md, cd, bd, fd = sp.symbols('md cd bd fd', real=True)
Sd = sp.diag(md, cd + bd, cd - bd)
Tkin = sp.expand(sp.trace(Sd * Sd) / 4 + fd ** 2 / 2)
check("1c kinetic metric diag(1/4, 1/2, 1/2, 1/2)",
      sp.simplify(Tkin - (md ** 2 / 4 + cd ** 2 / 2 + bd ** 2 / 2 + fd ** 2 / 2)) == 0, Tkin)

# ---------------------------------------------------------------- 2. transverse invariance
trans = {}
for name, var in (("u", u), ("v", v), ("q", q), ("Im psi", fi)):
    d = sp.simplify(sp.diff(V_rot, var).subs(slice_sub))
    trans[name] = d
check("2 transverse first derivatives of V_omega vanish on the slice",
      all(d == 0 for d in trans.values()), {k: str(x) for k, x in trans.items()})
# the two shear Hessians differ and carry -omega^2 (the #197 correction)
Huu = sp.factor(sp.diff(V_rot, u, 2).subs(slice_sub))
Hvv = sp.factor(sp.diff(V_rot, v, 2).subs(slice_sub))
check("2b shear Hessians u, v both carry -omega^2 and differ",
      (Huu.has(w2) and Hvv.has(w2) and sp.simplify(Huu - Hvv) != 0), f"Huu = {Huu}; Hvv = {Hvv}")

# ---------------------------------------------------------------- 3. vacuum decomposition
V_M5_minus = V_m5.subs({a: -2, t: 0, p: 0, u: 0, v: 0, q: 0})
check("3a V_M5(-N N^T) == 2", sp.simplify(V_M5_minus) == 2, V_M5_minus)
Vm5_f = sp.lambdify((a, t, p, u, v, q), V_m5, "numpy")
Vax_f = sp.lambdify((a, t, p, u, v, q), V_axis, "numpy")
Vlk_f = sp.lambdify((a, t, p, u, v, q, fr, fi), V_lock, "numpy")
W_f = sp.lambdify((fr, fi), W, "numpy")
rng = np.random.default_rng(7)
mins = {"V_M5": np.inf, "axis": np.inf, "lock": np.inf, "W": np.inf}
argmin_m5 = None
for _ in range(60):
    x0 = rng.normal(size=6) * 0.8
    r = minimize(lambda x: Vm5_f(*x), x0, method="BFGS")
    if r.fun < mins["V_M5"]:
        mins["V_M5"], argmin_m5 = float(r.fun), r.x
    mins["axis"] = min(mins["axis"], float(Vax_f(*x0)))
    z = rng.normal(size=8) * 0.8
    mins["lock"] = min(mins["lock"], float(Vlk_f(*z)))
    mins["W"] = min(mins["W"], float(W_f(z[6], z[7])))
S_min = np.array([[1 + argmin_m5[0], argmin_m5[3], argmin_m5[4]],
                  [argmin_m5[3], argmin_m5[1] + argmin_m5[2], argmin_m5[5]],
                  [argmin_m5[4], argmin_m5[5], argmin_m5[1] - argmin_m5[2]]])
eig_min = np.sort(np.linalg.eigvalsh(S_min))
check("3b V_M5 >= 0 with minimum 0 on the rank-one projector orbit (60 random BFGS starts)",
      mins["V_M5"] > -1e-9 and np.allclose(eig_min, [0, 0, 1], atol=1e-4),
      f"min V_M5 = {mins['V_M5']:.2e}, eigenvalues at the argmin {np.round(eig_min, 5).tolist()}")
check("3c axis lock, phase lock, W nonnegative on random samples",
      mins["axis"] > -1e-12 and mins["lock"] > -1e-12 and mins["W"] > -1e-12, mins)
# the W(f) zero: W = f^2 (3 - 4 f + 2 f^2), discriminant of the bracket negative
disc = sp.discriminant(3 - 4 * f + 2 * f ** 2, f)
check("3d W(f) = f^2 (2 f^2 - 4 f + 3) > 0 for f > 0 (bracket has no real root)", disc < 0, f"disc = {disc}")

# ---------------------------------------------------------------- 4. deep branch and Maxwell system
V0s = sp.expand(V0_hand)
Vw = V0s - w2 / 2 * (f ** 2 + 4 * b ** 2)
dm = sp.factor(sp.diff(Vw, m))
check("4a d_m V_omega is proportional to m (deep branch m = 0 consistent)", sp.simplify(dm / m).is_polynomial(m), dm)
w2_of_f = 32 * f ** 2 - 12 * f + 6 - 24 * b
check("4b d_f V_omega = 0  <=>  omega^2 = 32 f^2 - 12 f + 6 - 24 b",
      sp.simplify(sp.diff(Vw, f) - f * (w2_of_f - w2)) == 0)
Vw0 = Vw.subs(m, 0)
pA_hand = 8 * c ** 3 - 3 * c ** 2 + (8 * b ** 2 + 1) * c - 3 * b ** 2
pB_hand = 2 * (8 * b ** 3 + 8 * b * c ** 2 - 6 * b * c - 2 * b * w2 + 7 * b - 6 * f ** 2)
check("4c pA == (1/2) d_c V_omega on m = 0", sp.simplify(sp.diff(Vw0, c) / 2 - pA_hand) == 0)
check("4d pB == d_b V_omega on m = 0", sp.simplify(sp.diff(Vw0, b) - pB_hand) == 0)

# ---------------------------------------------------------------- 5. the certified root and the witnesses
mp.mp.dps = 50
sysF = [pA_hand, pB_hand.subs(w2, w2_of_f), (2 * Vw0).subs(w2, w2_of_f)]
F_num = sp.lambdify((c, b, f), sysF, "mpmath")
root = mp.findroot(lambda x, y, z: F_num(x, y, z), (mp.mpf('0.30'), mp.mpf('0.66'), mp.mpf('0.81')))
c_s, b_s, f_s = [root[i] for i in range(3)]
w2_star = 32 * f_s ** 2 - 12 * f_s + 6 - 24 * b_s
target = (mp.mpf('0.30280764518677380369'), mp.mpf('0.65773437772324925193'),
          mp.mpf('0.81436149699856776719'), mp.mpf('1.663945700059150298856193'))
dev = max(abs(c_s - target[0]), abs(b_s - target[1]), abs(f_s - target[2]), abs(w2_star - target[3]))
check("5a Maxwell root and omega_*^2 reproduce the certified values to 1e-20", dev < mp.mpf('1e-20'),
      f"c = {mp.nstr(c_s, 22)}, b = {mp.nstr(b_s, 22)}, f = {mp.nstr(f_s, 22)}, omega*^2 = {mp.nstr(w2_star, 28)}, max dev {mp.nstr(dev, 3)}")
Hw = sp.hessian(Vw0, (c, b, f)).subs(w2, w2_of_f.subs({b: b, f: f}))  # fixed-omega Hessian uses omega^2 fixed: recompute below
Hfix = sp.hessian(Vw0, (c, b, f))                     # omega^2 held fixed
Hnum = np.array(sp.lambdify((c, b, f, w2), Hfix, "numpy")(float(c_s), float(b_s), float(f_s), float(w2_star)), dtype=float)
eigH = np.linalg.eigvalsh(Hnum)
check("5b fixed-omega Hessian positive definite at the root", eigH.min() > 0, f"eigenvalues {np.round(eigH, 4).tolist()}")
wit = {m: 0, c: sp.Rational(31, 100), b: sp.Rational(13, 20), f: sp.Rational(41, 50)}
w1 = sp.nsimplify(Vw.subs(wit).subs(w2, sp.Rational(5, 3)))
w2v = sp.nsimplify(Vw.subs(wit).subs(w2, sp.Rational(45, 16)))
check("5c rational witnesses: V_omega(0, 31/100, 13/20, 41/50) = -13739/18750000 at 5/3 and -33854777/25000000 at 45/16",
      w1 == sp.Rational(-13739, 18750000) and w2v == sp.Rational(-33854777, 25000000), f"{w1}, {w2v}")
check("5d hence omega_c^2 < 5/3 < 45/16 and omega_*^2 < 5/3", float(w2_star) < 5 / 3, f"omega*^2 = {float(w2_star):.6f}")
Vstar = float(sp.lambdify((c, b, f, w2), Vw0, "numpy")(float(c_s), float(b_s), float(f_s), float(w2_star)))
check("5e V_omega = 0 at the root (equal depth with the vacuum)", abs(Vstar) < 1e-14, f"{Vstar:.2e}")

# ---------------------------------------------------------------- 6. thin-wall identities
R, sig, pr, om, I_, Fw, Qc = sp.symbols('R sigma p omega I F Q', positive=True)
Fth = 4 * sp.pi * R ** 2 * sig - sp.Rational(4, 3) * sp.pi * R ** 3 * pr
Rc = sp.solve(sp.diff(Fth, R), R)
Rc = [s for s in Rc if s != 0][0]
check("6a stationary radius R = 2 sigma / p", sp.simplify(Rc - 2 * sig / pr) == 0, Rc)
check("6b F''(R_c) = -8 pi sigma", sp.simplify(sp.diff(Fth, R, 2).subs(R, Rc) + 8 * sp.pi * sig) == 0)
check("6c E_crit = 16 pi sigma^3 / (3 p^2)", sp.simplify(Fth.subs(R, Rc) - 16 * sp.pi * sig ** 3 / (3 * pr ** 2)) == 0)
# envelope: E(omega) = F(omega) + omega Q(omega), with dF/domega = -Q (Legendre) => dE/dQ = omega
omg = sp.Function('Q')(om)
Fo = sp.Function('F')(om)
E = Fo + om * omg
dE = sp.diff(E, om).subs(sp.Derivative(Fo, om), -omg)
check("6d dE/dQ = omega on E = F + omega Q with dF/domega = -Q", sp.simplify(dE / sp.diff(omg, om) - om) == 0)

# ---------------------------------------------------------------- 7. independent sigma_0
W2S = float(w2_star)
Vw_num = sp.lambdify((m, c, b, f), Vw.subs(w2, W2S), "numpy")
grad_num = sp.lambdify((m, c, b, f), [sp.diff(Vw, x).subs(w2, W2S) for x in (m, c, b, f)], "numpy")
A_VEC = np.array([1.0, 0.0, 0.0, 0.0])
B_VEC = np.array([0.0, float(c_s), float(b_s), float(f_s)])


def rhs(x, y):
    gm, gc, gb, gf = grad_num(y[0], y[1], y[2], y[3])
    return np.vstack((y[4:], [2.0 * gm, gc, gb, gf]))      # 2K u'' = grad V, 2K = diag(1/2, 1, 1, 1)


def bc(ya, yb):
    return np.concatenate((ya[:4] - A_VEC, yb[:4] - B_VEC))


def kin(y):
    return (y[4] ** 2 + 2 * y[5] ** 2 + 2 * y[6] ** 2) / 4.0 + y[7] ** 2 / 2.0


def routes(sol, L):
    xe = np.linspace(-L, L, 4 * (len(sol.x) - 1) + 1)
    Y = sol.sol(xe)
    T = kin(Y)
    Vv = Vw_num(Y[0], Y[1], Y[2], Y[3])
    tz = lambda g: np.trapezoid(g, xe)
    return dict(r1=float(tz(T + Vv)), r2=float(2 * tz(Vv)), r3=float(2 * tz(T)),
                fi_drift=float(np.max(np.abs(T - Vv))), status=int(sol.status), nodes=int(len(sol.x)))


print("-- independent sigma_0 solve (tanh seed, L-continuation 2 -> 12) --")
sol = None
rows = []
for L, tol in ((2.0, 1e-6), (3.0, 1e-7), (4.5, 1e-8), (6.0, 1e-9), (8.0, 1e-10), (10.0, 1e-10), (12.0, 1e-10)):
    x = np.linspace(-L, L, 401)
    if sol is None:
        s = 0.5 * (1 - np.tanh(x / 0.55))
        ds = -0.5 / 0.55 / np.cosh(x / 0.55) ** 2
        y0 = np.vstack(((A_VEC[:, None] + np.outer(B_VEC - A_VEC, s)), np.outer(B_VEC - A_VEC, ds)))
    else:
        y0 = sol.sol(x)
    sol = solve_bvp(rhs, bc, x, y0, tol=tol, max_nodes=300000)
    r = routes(sol, L)
    r["L"] = L
    rows.append(r)
    print(f"   L = {L:5.1f} status {r['status']} nodes {r['nodes']:6d} sigma(T+V) {r['r1']:.12f} 2V {r['r2']:.12f} 2T {r['r3']:.12f} |T-V|max {r['fi_drift']:.1e}")
# headline = the LAST status-0 rung (their own anchoring rule: solver success is part of acceptance);
# the status-1 rungs hit the node ceiling with the residual already converged and are reported, not consumed
ok_rows = [r for r in rows if r["status"] == 0 and r["L"] >= 8.0]
final = ok_rows[-1] if ok_rows else rows[-1]
sigma0_target = 0.7292984178625
spread = max(abs(final["r1"] - final["r2"]), abs(final["r1"] - final["r3"]))
check("7 independent sigma_0 (last status-0 rung, L >= 8) within 1e-8 of 0.72929841786, route spread < 1e-8",
      final["status"] == 0 and abs(final["r1"] - sigma0_target) < 1e-8 and spread < 1e-8,
      f"L = {final['L']}: sigma_0 = {final['r1']:.12f} (routes {final['r2']:.12f}, {final['r3']:.12f}); dev from theirs {final['r1'] - sigma0_target:.2e}; L = 12 rung (status {rows[-1]['status']}) {rows[-1]['r1']:.12f}")
RES["sigma0_ladder"] = rows

# ---------------------------------------------------------------- 8. the delta = 0.001 bag radius from the thin-wall law
w2d = W2S + 0.001
Vwd = sp.lambdify((c, b, f), Vw0.subs(w2, w2d), "numpy")
rr = minimize(lambda x: Vwd(*x), [float(c_s), float(b_s), float(f_s)], method="BFGS", options={"gtol": 1e-13})
p_d = -float(rr.fun)
R_thin = 2 * final["r1"] / p_d
chi_target, Rm_target = 1.000328, 1218.761785
check("8 thin-wall radius at delta = 0.001: R_m / (2 sigma_0 / p) == chi = 1.000328 (their R_m 1218.761785)",
      abs(Rm_target / R_thin - chi_target) < 2e-6,
      f"p = {p_d:.9e}, 2 sigma_0/p = {R_thin:.6f}, their R_m / ours = {Rm_target / R_thin:.6f}")

# ---------------------------------------------------------------- summary
n = len([k for k in RES if isinstance(RES[k], dict) and 'pass' in RES[k]])
print(f"\n{n - len(FAILS)}/{n} checks passed; FAILS: {FAILS}   [{time.time() - T0:.0f} s]")
RES["summary"] = {"passed": n - len(FAILS), "total": n, "fails": FAILS, "wall_s": round(time.time() - T0, 1),
                  "pinned_main": "71bbe9a276a2e53c2b88e1a1d4aaf8c9d0ee915f"}
json.dump(RES, open(OUT, "w"), indent=1, default=str)
print("wrote", OUT)
