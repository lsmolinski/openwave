"""M5.32 R15-V-a (the floor witness on homogeneous jets) and R15-H (the tilt channel), symbolic.

Equations first.  eta = diag(-1, 1, 1, 1); M(x) 4x4 symmetric contravariant; A_mu = d_mu M;
F_mu nu = A_mu eta A_nu - A_nu eta A_mu;  <F, F>_G = tr(G F G F^T) for an internal metric G.
    I1   = sum_{mu<nu} eta^mu eta^nu <F_mu nu, F_mu nu>_eta          (the certified term, L_cert = -4 I1 - V4)
    I1_h = sum_{mu<nu} eta^mu eta^nu <F_mu nu, F_mu nu>_h,  h = eta + 2 (eta u)(eta u)^T,
           u the timelike unit eigenvector of N = M eta (u = Q e0 for M = Q D Q^T, Q in SO(1,3))
E-orientation (static): E = +4 x I1 (resp. I1_h) at omega = 0.

R15-V-a.  The jet of the author's floor witness at one point: the vacuum D = diag(g, 1, delta, delta4)
dressed by a boost along axis a with rapidity chi (gradient b along axis s) and a (1,2) twist by
psi (gradient k along axis t = 3), TWIST INSIDE:  M = L_a(chi) R_12(psi) D R_12^T L_a^T,
TWIST AFTER:  M = R_12(psi) L_a(chi) D L_a^T R_12^T.  The only curvature is F_st = [b d_chi M, k d_psi M]_eta
(s != t), and the static E-density of each form is U_G = 4 <F_st, F_st>_G = b^2 k^2 c_G(chi, psi).
Claims (pre-registered, ledger 6.4): V1 c_eta < 0 (twist inside, the author's -2 b^2 k^2 g^2 (delta-1)^2 at
chi -> 0); V2 c_h > 0 (the author's +2 b^2 k^2 g^2 (delta-1)^2); V3 twist after: c_eta == c_h; V4 the exact
symbolic coefficients.

R15-H.  The rotating split background M = R_23(omega t) R_12(theta(t, z)) D_s R_12^T R_23^T,
D_s = diag(g, 1, delta + s, delta - s), the tilt theta INSIDE the rotating frame.  A_0 = omega d_tau M + theta_t d_theta M,
A_z = theta_z d_theta M; only F_0z is nonzero.  Lagrangian densities (T - V orientation):
    L(-4 I1)   = +4 <F_0z, F_0z>_eta            (eta^0 eta^z = -1)
    L(-4 I1_h) = +4 <F_0z, F_0z>_h              (u = e0 here: h = 1, the Frobenius form)
    L(K_P^23)  = (c_P / 2) [tr(Om_0^T eta Om_0 eta) - tr(Om_z^T eta Om_z eta)],  Om_mu = P23 A_mu eta P23,
                 P23 = Q diag(0, 0, 1, 1) Q^T (= the exact Lagrange projector, Q the spatial rotation)
    L(reg)     = w [tr(A_0 G A_0 G) - tr(A_z G A_z G)],  G in {eta, h} (h = 1 here); kappa_2 := tr(Y G Y G),
                 Y = d_theta M at theta = 0.
    V4^dd and mu (lambda_2 - lambda_3)^2 depend on the spectrum only: constant in theta, dropped.
Quadratic form in (theta, theta_t, theta_z):  L_2 = alpha theta_t^2 + gamma theta_z^2 + eps theta^2 + cross terms.
Plane wave e^{i(k z - Omega t)}: Omega^2 = -(gamma k^2 + eps) / alpha.  Hyperbolic iff alpha > 0 and gamma < 0.
Claims: H1 alpha = 0 for every curvature term and for K_P^23 (only the regulator supplies alpha = w kappa_2);
H2 gamma(-4 I1) = 4 omega^2 |C|^2 with |C|^2 = 4 s^2 (delta + s - 1)^2 (the author's <F_0z, F_0z> = 4 k^2 omega^2 s^2 (delta+s-1)^2);
H3 the hyperbolicity inequality in our normalization (compare the author's w kappa_2 > 16 omega^2 s^2);
H4 K_P^23 has no theta_t^2 and no theta_z^2 term, and its STATIC density is exactly zero on any (1,2) twist sheet
theta = psi(z) of any amplitude (the projector kills the (1,2) block), so the sheet is a free direction of K_P^23.

usage: python3 m5_32_r15_vh_symbolic.py va | h
"""
import sys
ARGS = list(sys.argv[1:])
import os, json, time
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.dirname(HERE)
DATA = os.path.join(RES, "data")
T0 = time.time()


def log(m):
    print(f"[{time.time() - T0:8.1f}s] {m}", flush=True)


ETA = sp.diag(-1, 1, 1, 1)
E0 = sp.Matrix([1, 0, 0, 0])


def boost(a, chi):
    L = sp.eye(4)
    c, s = sp.cosh(chi), sp.sinh(chi)
    L[0, 0] = c; L[a, a] = c; L[0, a] = s; L[a, 0] = s
    return L


def rot(p, q, psi):
    R = sp.eye(4)
    c, s = sp.cos(psi), sp.sin(psi)
    R[p, p] = c; R[q, q] = c; R[p, q] = -s; R[q, p] = s
    return R


def comm(A, B):
    return A * ETA * B - B * ETA * A


def inner(F, G, Gm):
    return (Gm * F * Gm * G.T).trace()


def hcov(u):
    v = ETA * u
    return ETA + 2 * v * v.T


# ------------------------------------------------ R15-V-a
def va_mode():
    g, d, d4, b, k = sp.symbols("g delta delta4 b k", positive=True)
    chi, psi = sp.symbols("chi psi", real=True)
    D = sp.diag(g, 1, d, d4)
    out = {"definitions": __doc__.split("R15-H.")[0], "author_coefficient": "-2 b^2 k^2 g^2 (delta-1)^2 (eta), +2 b^2 k^2 g^2 (delta-1)^2 (h)",
           "pairs": []}
    t = 3
    for order in ("inside", "after"):
        for a in (1, 2, 3):
            La = boost(a, chi)
            R = rot(1, 2, psi)
            Q = La * R if order == "inside" else R * La
            M = Q * D * Q.T
            u = Q * E0
            H = hcov(u)
            dchi, dpsi = sp.diff(M, chi), sp.diff(M, psi)
            for s_ax in (1, 2, 3):
                if s_ax == t:
                    out["pairs"].append({"order": order, "boost_axis": a, "grad_axis": s_ax, "twist_axis": t,
                                         "c_eta": "0 (s = t: one gradient direction, F = 0)", "c_h": "0"})
                    continue
                F = comm(b * dchi, k * dpsi)
                Ue = sp.expand(4 * inner(F, F, ETA))
                Uh = sp.expand(4 * inner(F, F, H))
                ce0 = sp.factor(sp.simplify((Ue / (b * b * k * k)).subs({chi: 0, psi: 0})))
                ch0 = sp.factor(sp.simplify((Uh / (b * b * k * k)).subs({chi: 0, psi: 0})))
                nums = {}
                for d4v, lab in ((sp.Rational(3, 10), "degenerate"), (0, "certified")):
                    for chiv in (0, sp.Rational(1, 2)):
                        for psiv in (0, sp.pi / 4):
                            sub = {g: 8, d: sp.Rational(3, 10), d4: d4v, chi: chiv, psi: psiv, b: 1, k: 1}
                            nums[f"{lab}_chi{float(chiv):g}_psi{float(psiv):.3f}"] = [float(Ue.subs(sub)), float(Uh.subs(sub))]
                rec = {"order": order, "boost_axis": a, "grad_axis": s_ax, "twist_axis": t,
                       "c_eta_chi0_psi0": str(ce0), "c_h_chi0_psi0": str(ch0),
                       "c_eta_minus_c_h_general": str(sp.factor(sp.simplify((Ue - Uh) / (b * b * k * k)))),
                       "numeric_[U_eta,U_h]_b1_k1": nums}
                out["pairs"].append(rec)
                log(f"V-a {order:6s} boost {a} grad {s_ax} twist(1,2) along {t}: c_eta(0,0) = {ce0}   c_h(0,0) = {ch0}")
                log(f"      numbers g8 d0.3: " + ", ".join(f"{kk}: {v[0]:+.3f}/{v[1]:+.3f}" for kk, v in nums.items()))
    json.dump(out, open(os.path.join(DATA, "m5_32_r15_va_jets.json"), "w"), indent=1)
    log("wrote data/m5_32_r15_va_jets.json")


# ------------------------------------------------ R15-H
def h_mode():
    g, d, s, om, cP, w = sp.symbols("g delta s omega c_P w", positive=True)
    tau, th, tht, thz, eps = sp.symbols("tau theta theta_t theta_z epsilon", real=True)
    Ds = sp.diag(g, 1, d + s, d - s)
    Q = rot(2, 3, tau) * rot(1, 2, th)
    M = Q * Ds * Q.T
    A0 = om * sp.diff(M, tau) + tht * sp.diff(M, th)
    Az = thz * sp.diff(M, th)
    F0z = comm(A0, Az)
    P23 = Q * sp.diag(0, 0, 1, 1) * Q.T
    Om0 = P23 * A0 * ETA * P23
    Omz = P23 * Az * ETA * P23
    ONE = sp.eye(4)
    dens = {
        "-4 I1": 4 * inner(F0z, F0z, ETA),
        "-4 I1_h (h = 1 on this background)": 4 * inner(F0z, F0z, ONE),
        "K_P^23": (cP / 2) * ((Om0.T * ETA * Om0 * ETA).trace() - (Omz.T * ETA * Omz * ETA).trace()),
        "regulator G = eta": w * ((A0 * ETA * A0 * ETA).trace() - (Az * ETA * Az * ETA).trace()),
        "regulator G = h (= 1)": w * ((A0 * A0).trace() - (Az * Az).trace()),
    }
    Y = sp.diff(M, th).subs({th: 0, tau: 0})
    kappa2 = {"eta": sp.factor((Y * ETA * Y * ETA).trace()), "h": sp.factor((Y * Y).trace())}
    C = comm(sp.diff(M, tau), sp.diff(M, th)).subs({th: 0})
    C2 = sp.factor(sp.simplify(inner(C, C, ETA)))
    out = {"definitions": "R15-H." + __doc__.split("R15-H.")[1], "kappa_2": {k_: str(v) for k_, v in kappa2.items()},
           "|C|^2 = <[d_tau M, d_theta M]_eta, same>_eta at theta 0": str(C2), "terms": {}}
    log(f"H: kappa_2(eta) = {kappa2['eta']}, kappa_2(h) = {kappa2['h']}, |C|^2 = {C2} (the author: 4 s^2 (delta+s-1)^2)")
    def taylor2(Ld):
        """Taylor coefficients in the tilt amplitude by differentiation (sympy's series chokes on the
        projector trig products): c_n = (1/n!) d^n L / d eps^n at eps = 0."""
        Ls = Ld.subs({th: eps * th, tht: eps * tht, thz: eps * thz})
        cs = {}
        cur = Ls
        for order in (0, 1, 2):
            cs[order] = sp.simplify(sp.expand(cur.subs(eps, 0) / sp.factorial(order)))
            cur = sp.diff(cur, eps)
        return cs

    for name, Ld in dens.items():
        coeffs = taylor2(Ld)
        L2 = sp.Poly(sp.expand(coeffs[2]), th, tht, thz)
        q = {}
        for mon, coef in zip(L2.monoms(), L2.coeffs()):
            key = "*".join(f"{v}^{e}" if e > 1 else f"{v}" for v, e in zip(("theta", "theta_t", "theta_z"), mon) if e)
            q[key] = str(sp.factor(sp.simplify(coef)))
        tau_dep = any(sp.simplify(sp.diff(c, tau)) != 0 for c in coeffs.values())
        rec = {"L_0 (background)": str(sp.factor(coeffs[0])), "L_1": str(sp.factor(coeffs[1])), "L_2 quadratic form": q,
               "tau_dependent": bool(tau_dep)}
        out["terms"][name] = rec
        log(f"H term {name}: L0 = {rec['L_0 (background)']}; L1 = {rec['L_1']}; L2 = {q}; tau-dependent {tau_dep}")
    # the hyperbolicity inequality: alpha from the regulators only, gamma from -4 I1 (K_P adds no theta_z^2)
    def coef_of(Ld, mon):
        return sp.factor(sp.Poly(sp.expand(taylor2(Ld)[2]), th, tht, thz).coeff_monomial(mon))
    alpha_e = coef_of(dens["regulator G = eta"], tht ** 2)
    gam_e = coef_of(dens["regulator G = eta"], thz ** 2)
    gam_i1 = coef_of(dens["-4 I1"], thz ** 2)
    gam_kp = coef_of(dens["K_P^23"], thz ** 2)
    gam_i1 = gam_i1 + gam_kp
    cond = sp.simplify(-(gam_e + gam_i1))
    out["hyperbolicity"] = {"alpha_total": str(alpha_e), "gamma_total": str(sp.factor(gam_e + gam_i1)),
                            "hyperbolic_iff": f"alpha > 0 and gamma < 0  <=>  {sp.factor(cond)} > 0",
                            "author": "w kappa_2 > 16 omega^2 s^2"}
    log(f"H hyperbolicity: alpha = {alpha_e}, gamma = {sp.factor(gam_e + gam_i1)}  ->  hyperbolic iff {sp.factor(cond)} > 0")
    # H4 exact: the static K_P^23 density on a finite (1,2) twist sheet psi(z) (omega = 0)
    psi, dpsi = sp.symbols("psi psi_z", real=True)
    Qs = rot(1, 2, psi)
    Ms = Qs * Ds * Qs.T
    Azs = dpsi * sp.diff(Ms, psi)
    P23s = Qs * sp.diag(0, 0, 1, 1) * Qs.T
    Omzs = P23s * Azs * ETA * P23s
    kp_static_sheet = sp.simplify((Omzs.T * ETA * Omzs * ETA).trace())
    i1_static_sheet = sp.simplify(inner(comm(Azs, Azs), comm(Azs, Azs), ETA))
    out["H4_static_on_finite_twist_sheet"] = {"K_P^23 static density": str(kp_static_sheet), "I1 static density": str(i1_static_sheet)}
    log(f"H4: static K_P^23 density on the finite (1,2) twist sheet = {kp_static_sheet}; I1 static = {i1_static_sheet}")
    json.dump(out, open(os.path.join(DATA, "m5_32_r15_h_tilt.json"), "w"), indent=1)
    log("wrote data/m5_32_r15_h_tilt.json")


if __name__ == "__main__":
    {"va": va_mode, "h": h_mode}[ARGS[0]]()
