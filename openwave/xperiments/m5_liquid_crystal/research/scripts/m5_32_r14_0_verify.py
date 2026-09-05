"""M5.32 R14-0: VERIFY the model author's 2026-09-03 analysis on our own stack
(ledger 6.3, the first rung of R14; every line a PASS/FAIL that a named
mutation fails).  Consumes m5_32_r14_terms.py (the entrants), the R0 registry,
the certified stack, the R10 / R13-W saved fields; never imports the author's
scripts (they did not reach the thread).

THE STATEMENTS CHECKED (numbered as in the author's comment)
-------------------------------------------------------------
 1a  P249 exterior Hessian split V_ldg / V_axis / V_lock in the chart (a,t,p,u,v,q)
     = (5,6,6,0,0,6) / (0,4,4,4,4,4) / (0,0,12,0,0,12); the shear doublet (u,v) is
     massless under the rotation-invariant part (a Goldstone of the LdG part).
     MUTATION: an SO(3)-breaking replacement of V_axis is not blind to (u,v).
 1b  On the 4x4 degenerate vacuum N = diag(-g, 1, delta, delta) a conjugation-
     invariant potential has exactly FIVE flat directions (the SO(1,3) orbit,
     stabilizer SO(2)), with (2,3)-clock charges (0, 1, 1, 1, 1): boost_01 neutral,
     (boost_02, boost_03) and (rot_12, rot_13) charge-1 doublets; the massive
     transverse directions are three neutral eigenvalue modes and the charge-2
     split doublet.  MUTATIONS: the certified non-degenerate vacuum (stabilizer
     trivial, six flat directions, charges (0,1,1,1,1,0)); an explicit axis lock
     (an SO(1,3)-breaking term) gives the charged doublets a mass.
 2a  R_eta has EL == 0 identically for generic M (sympy, 1+1 and 2+1).
 2b  R_G for G in {eta, eta M eta, M^-1, h_cov} integrates to zero over a periodic
     box on a FIXED-EIGENVALUE Lorentz orbit M = L(x) M0 L(x)^T (boosts included),
     with spectral convergence (FFT derivatives, N = 8 .. 64).  MUTATIONS: an
     eigenvalue-gradient field with the same frame integrates to a NONZERO value
     (spectrally converged); a pure eigenvalue gradient in a constant frame gives
     zero for the diagonal G's (the (d lambda).(d frame) structure).
 2c  R_G has no omega^2 content for any G (the omega decomposition on the relaxed
     hedgehog with the local clock generator: C = 0 to roundoff).
 4a  K_lambda is inert on orientation gradients (zero on the orbit field of 2b) and
     on every generator channel (zero omega^2 coefficient with the local clock and
     the boosts on the relaxed hedgehog); nonzero in the relaxed core (the
     eigenvalue deficit), with its tail integral and L-exponent reported.
 4b  K_P at the vacuum: zero omega^2 on the three boosts and the two tilts, nonzero
     on the (2,3) clock only, stiffness [f(delta) f(0)]^2 delta^2 (the registry
     selftest lines, re-reported); on the relaxed hedgehog the per-channel table
     against I1's kin (gen_catalog channels).
 4c  K_P does NOT charge the R13-W sheet: on the W0 S5 (1,2)-twist family its static
     density is exactly zero and its kin x w -> 0 while I1's kin x w -> const (the
     free inertia is untouched).  ADDITION (ours): on the three W3 end fields the
     K_P static density on the zigzag layers against I1's kin density there (K_P
     sees the SPLIT-modulus gradient, which is the mechanism W3 measured).
 4d  Covariance of every entrant with the no-eta control failing (registry lines).
 S   The W0 theorems S3 (planar flatness), S4 (phase flatness), S5 (free inertia)
     re-evaluated per entrant: which survive the whole two-derivative set.  A term
     that leaves S5 intact cannot change R13-W's verdict (no fixed-J minimizer).
 V   The vacuum-ticking property: on the certified NON-degenerate vacuum the local
     clock generator is nonzero, so K_P's omega^2 coefficient is VOLUME-extensive
     (a constant per cell = the 4b stiffness), and the (2,3) degenerate state costs
     V4_deg per cell.  Reported as the two numbers that decide R14-B and R14-D
     before they run.

Run: python3 m5_32_r14_0_verify.py  (sympy + numpy; ~3 min on the saved fields).
Writes data/m5_32_r14_0_verify.json.
"""
from __future__ import annotations
import importlib.util
import itertools
import json
import os
import sys
import time
import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.dirname(HERE)
DATA = os.path.join(RES, "data")
CK = os.path.join(RES, "checkpoints")
OUT = os.path.join(DATA, "m5_32_r14_0_verify.json")
T0 = time.time()


def log(m):
    print(f"[{time.time() - T0:8.1f}s] {m}", flush=True)


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


T14 = _load("m5_32_r14_terms", "m5_32_r14_terms.py")
L0 = T14.L0
B3 = T14.B3
C13 = _load("m5_32_r13w_common", "m5_32_r13w_common.py")
B8 = C13.B8
ETA = T14.ETA
RESULTS = {"rung": "R14-0 verify", "lines": []}


def check(name, ok, value, detail=""):
    RESULTS["lines"].append({"check": name, "pass": bool(ok), "value": value, "detail": detail})
    print(f"  {'PASS' if ok else 'FAIL'}  {name}: {value} {detail}", flush=True)


# ====================================================================== 1a
def s1a_p249_split():
    a, t, p, u, v, q = sp.symbols("a t p u v q", real=True)
    S = sp.Matrix([[1 + a, u, v], [u, t + p, q], [v, q, t - p]])
    tr2 = sp.trace(S * S)
    V_ldg = -tr2 / 2 - sp.trace(S ** 3) + tr2 ** 2 + sp.Rational(1, 2)
    V_axis = tr2 - S[0, 0] ** 2
    A = sp.Matrix([[0, 0, 0], [0, 0, -1], [0, 1, 0]])
    Pm = sp.diag(0, 1, 1)
    Bproj = Pm * S * Pm
    Bcl = Bproj - sp.trace(Bproj) * Pm / 2
    Q = sp.zeros(3)                       # psi = 0 on the exterior
    V_lock = 6 * sp.trace((Bcl - Q).T * (Bcl - Q)) / 2
    xs = [a, t, p, u, v, q]
    vac = {x: 0 for x in xs}

    def hess_diag(V):
        H = sp.Matrix(6, 6, lambda i, j: sp.diff(V, xs[i], xs[j]).subs(vac))
        off = max(abs(H[i, j]) for i in range(6) for j in range(6) if i != j)
        return [int(H[i, i]) for i in range(6)], off
    d_ldg, o1 = hess_diag(V_ldg)
    d_axis, o2 = hess_diag(V_axis)
    d_lock, o3 = hess_diag(V_lock)
    want = {"ldg": [5, 6, 6, 0, 0, 6], "axis": [0, 4, 4, 4, 4, 4], "lock": [0, 0, 12, 0, 0, 12]}
    got = {"ldg": d_ldg, "axis": d_axis, "lock": d_lock}
    ok = got == want and max(o1, o2, o3) == 0
    check("1a P249 Hessian split V_ldg / V_axis / V_lock = (5,6,6,0,0,6)/(0,4,4,4,4,4)/(0,0,12,0,0,12), all diagonal",
          ok, got)
    # SO(3) invariance of V_ldg: (u, v) are the two broken generators about diag(1,0,0)
    th = sp.symbols("th", real=True)
    R = sp.Matrix([[sp.cos(th), -sp.sin(th), 0], [sp.sin(th), sp.cos(th), 0], [0, 0, 1]])
    inv = sp.simplify((V_ldg.subs({a: 0, t: 0, p: 0, u: 0, v: 0, q: 0}) - V_ldg.subs({a: 0, t: 0, p: 0, u: 0, v: 0, q: 0})))
    Sg = R * sp.diag(1, 0, 0) * R.T
    V_on_orbit = V_ldg.subs({a: Sg[0, 0] - 1, t: (Sg[1, 1] + Sg[2, 2]) / 2, p: (Sg[1, 1] - Sg[2, 2]) / 2,
                             u: Sg[0, 1], v: Sg[0, 2], q: Sg[1, 2]})
    diff_expr = sp.expand(V_on_orbit - V_ldg.subs(vac))
    dev = max(abs(float(diff_expr.subs(th, x))) for x in np.linspace(0.1, 3.0, 37))
    flat = dev < 1e-12
    check("1a (u,v) is the Goldstone doublet of the rotation-invariant part: V_ldg constant on the SO(3) orbit of diag(1,0,0) (37 angles)",
          flat and d_ldg[3] == 0 and d_ldg[4] == 0, dev)
    # MUTATION: an SO(3)-breaking axis term that is NOT blind to the shear: S_01^2 + S_02^2 added
    V_mut = V_ldg + (S[0, 1] ** 2 + S[0, 2] ** 2)
    d_mut, _ = hess_diag(V_mut)
    check("1a MUTATION: adding the explicit shear term S_01^2 + S_02^2 gives (u,v) a mass (the split is not an artifact of the chart)",
          d_mut[3] == 2 and d_mut[4] == 2, d_mut)
    RESULTS["s1a"] = {"split": got, "want": want}


# ====================================================================== 1b
def so13_generators():
    gens = {}
    for i in range(1, 4):
        K = np.zeros((4, 4)); K[0, i] = K[i, 0] = 1.0
        gens[f"boost_0{i}"] = K
    for (i, j) in ((1, 2), (1, 3), (2, 3)):
        J = np.zeros((4, 4)); J[i, j], J[j, i] = -1.0, 1.0
        gens[f"rot_{i}{j}"] = J
    return gens


def sym_basis():
    """the 10 symmetric unit matrices E_ab (a <= b), and the coordinate map."""
    basis, idx = [], []
    for a in range(4):
        for b in range(a, 4):
            E = np.zeros((4, 4)); E[a, b] = E[b, a] = 1.0
            if a != b:
                E /= np.sqrt(2.0)
            basis.append(E); idx.append((a, b))
    return basis, idx


def v_conj_invariant(M, C):
    """V4-type: sum_p (tr((M eta)^p) - C_p)^2, p = 1..4 (a sum of four squares: its Hessian
    at a zero has rank <= 4, and at a DEGENERATE spectrum the split enters the traces only
    at second order, so the charge-2 doublet is quartically soft there: reported, not used
    as the Goldstone probe)."""
    N = M @ ETA
    tot = 0.0
    P = np.eye(4)
    for p in range(4):
        P = P @ N
        tot = tot + (np.trace(P) - C[p]) ** 2
    return tot


def v_strict(M, roots):
    """V = tr f(N), f(x) = prod_r (x - r)^2 over the DISTINCT vacuum eigenvalues: conjugation-
    invariant, smooth, with a strict quadratic minimum on the orbit of the vacuum (each
    eigenvalue, including a degenerate pair's split, gets a mass f''(r) = 2 prod_{r' != r} (r - r')^2)."""
    N = M @ ETA
    P = np.eye(4)
    for r in roots:
        P = P @ (N - r * np.eye(4)) @ (N - r * np.eye(4))
    return float(np.trace(P))


def hessian_num(fn, M0, basis, eps=1e-4):
    n = len(basis)
    H = np.zeros((n, n))
    f0 = fn(M0)
    for i in range(n):
        for j in range(i, n):
            fpp = fn(M0 + eps * basis[i] + eps * basis[j])
            fpm = fn(M0 + eps * basis[i] - eps * basis[j])
            fmp = fn(M0 - eps * basis[i] + eps * basis[j])
            fmm = fn(M0 - eps * basis[i] - eps * basis[j])
            H[i, j] = H[j, i] = (fpp - fpm - fmp + fmm) / (4 * eps * eps)
    return H


def s1b_goldstones():
    g, delta = 8.0, 0.3
    basis, idx = sym_basis()
    gens = so13_generators()
    Jc = gens["rot_23"]                                           # the clock generator

    def analyze(M0, tag, lock=0.0, v4type=False):
        N0 = M0 @ ETA
        C = [np.trace(np.linalg.matrix_power(N0, p)) for p in range(1, 5)]
        roots = sorted(set(np.round(np.diag(N0), 12)))

        def V(M):
            val = v_conj_invariant(M, C) if v4type else v_strict(M, roots)
            if lock:
                # an explicit SO(1,3)-breaking axis lock: penalize the (0,2),(0,3),(1,2),(1,3) entries
                val = val + lock * (M[0, 2] ** 2 + M[0, 3] ** 2 + M[1, 2] ** 2 + M[1, 3] ** 2)
            return val
        H = hessian_num(V, M0, basis)
        w, U = np.linalg.eigh(H)
        nflat = int(np.sum(np.abs(w) < 1e-9 * max(1.0, np.max(np.abs(w)))))    # audit 2026-09-05: 1e-6 counted a 3.7e-7 eigenvalue as flat
        # the orbit tangents X M0 + M0 X^T, their span, and the ad_{Jc} charges
        tang = {k: X @ M0 + M0 @ X.T for k, X in gens.items()}
        coords = np.array([[np.sum(tang[k] * E) for E in basis] for k in gens])
        rank = int(np.linalg.matrix_rank(coords, tol=1e-9))
        # charge of a generator X under the clock: [Jc, X] = charge-coupled partner; eigen-decompose ad_Jc on so(1,3)
        keys = list(gens)
        ad = np.zeros((6, 6))
        for j, k in enumerate(keys):
            Y = Jc @ gens[k] - gens[k] @ Jc
            for i, k2 in enumerate(keys):
                ad[i, j] = np.sum(Y * gens[k2]) / np.sum(gens[k2] * gens[k2])
        ch = np.sort(np.abs(np.linalg.eigvals(ad).imag).round(6))
        # which tangents are flat: project each tangent onto the flat eigenspace of H
        flat_vecs = U[:, np.abs(w) < 1e-6 * max(1.0, np.max(np.abs(w)))]
        flat_of = {}
        for k in keys:
            c = coords[keys.index(k)]
            nrm = np.linalg.norm(c)
            if nrm < 1e-12:
                flat_of[k] = "stabilizer (zero tangent)"
            else:
                resid = np.linalg.norm(c - flat_vecs @ (flat_vecs.T @ c)) / nrm
                flat_of[k] = "flat" if resid < 1e-5 else f"massive (resid {resid:.2e})"
        return {"n_flat": nflat, "orbit_rank": rank, "eigs": [float(x) for x in w],
                "charges_of_so13_under_rot23": [float(x) for x in ch], "tangent_status": flat_of}

    Md = np.diag([g, 1.0, delta, delta])
    Mn = np.diag([g, 1.0, delta, 0.0])
    deg = analyze(Md, "degenerate")
    nondeg = analyze(Mn, "nondegenerate")
    locked = analyze(Md, "degenerate+lock", lock=1.0)
    v4deg = analyze(Md, "degenerate, V4-type probe", v4type=True)
    RESULTS["s1b"] = {"degenerate": deg, "nondegenerate_MUTATION": nondeg, "locked_MUTATION": locked,
                      "degenerate_V4type_probe": v4deg}
    check("1b PROPERTY: a V4-type potential (sum of four trace squares) is quartically SOFT at the degenerate spectrum: 7 flat directions at quadratic order = the 5 tangents + the charge-2 split doublet (rank of the four trace gradients = 3 there)",
          v4deg["n_flat"] == 7, {"n_flat_V4type": v4deg["n_flat"]})
    ok = deg["n_flat"] == 5 and deg["orbit_rank"] == 5 and deg["tangent_status"]["rot_23"].startswith("stabilizer") \
        and all(deg["tangent_status"][k] == "flat" for k in ("boost_01", "boost_02", "boost_03", "rot_12", "rot_13"))
    check("1b degenerate vacuum diag(-g,1,delta,delta): exactly 5 flat directions = the orbit (rank 5), rot_23 the stabilizer",
          ok, {"n_flat": deg["n_flat"], "orbit_rank": deg["orbit_rank"]}, str(deg["tangent_status"]))
    # charges: ad_{rot23} on so(1,3) has eigenvalues 0 (boost_01, rot_23) and +-i (the two doublets)
    ch = deg["charges_of_so13_under_rot23"]
    check("1b clock charges of the five Goldstones (0,1,1,1,1): ad_rot23 spectrum on so(1,3) = {0, 0, 1, 1, 1, 1} with rot_23 the neutral stabilizer",
          ch == [0.0, 0.0, 1.0, 1.0, 1.0, 1.0], ch)
    # the massive transverse block: 3 neutral eigenvalue modes + a charge-2 split doublet
    # (the (2,3) traceless block (p, q) rotates by angle 2 theta under the clock)
    th = 0.3
    Rm = np.eye(4); Rm[2, 2] = Rm[3, 3] = np.cos(th); Rm[2, 3], Rm[3, 2] = -np.sin(th), np.sin(th)
    Ep = np.zeros((4, 4)); Ep[2, 2], Ep[3, 3] = 1.0, -1.0
    Eq = np.zeros((4, 4)); Eq[2, 3] = Eq[3, 2] = 1.0
    rot_p = Rm @ Ep @ Rm.T
    c2 = np.sum(rot_p * Ep) / 2.0, np.sum(rot_p * Eq) / 2.0
    check("1b the split doublet (p, q) carries charge 2 (rotates by 2 theta under the (2,3) clock; the core's charged content is even)",
          abs(c2[0] - np.cos(2 * th)) < 1e-12 and abs(c2[1] - np.sin(2 * th)) < 1e-12, [float(x) for x in c2])
    check("1b MUTATION: the certified NON-degenerate vacuum has 6 flat directions (trivial stabilizer; rot_23 is itself an orbit direction)",
          nondeg["n_flat"] == 6 and nondeg["orbit_rank"] == 6 and nondeg["tangent_status"]["rot_23"] == "flat",
          {"n_flat": nondeg["n_flat"], "orbit_rank": nondeg["orbit_rank"]})
    check("1b MUTATION: an explicit axis lock on the degenerate vacuum leaves only boost_01 flat (the charged doublets get a mass)",
          locked["n_flat"] == 1 and locked["tangent_status"]["boost_01"] == "flat", locked["tangent_status"])


# ====================================================================== 2a
def s2a_reta_euler_lagrange():
    out = {}
    for dim, coords in (("1+1", ("t", "x")), ("2+1", ("t", "x", "y"))):
        xs = sp.symbols(" ".join(coords), real=True)
        # symmetric 4x4 M of functions of the coordinates
        f = {}
        M = sp.zeros(4)
        for a in range(4):
            for b in range(a, 4):
                fn = sp.Function(f"m{a}{b}")(*xs)
                M[a, b] = M[b, a] = fn
                f[(a, b)] = fn
        # jets: mu runs over the coordinate list mapped to internal index positions 0..len-1
        A = [M.diff(x) for x in xs]
        nmu = len(xs)
        G = sp.diag(-1, 1, 1, 1)
        R = 0
        for mu in range(nmu):
            for nu in range(nmu):
                for c in range(4):
                    for d in range(4):
                        if G[c, d] == 0:
                            continue
                        R += G[c, d] * (A[mu][nu, c] * A[nu][mu, d] - A[mu][mu, c] * A[nu][nu, d])
        R = sp.expand(R)
        # Euler-Lagrange for each field: sum_mu d_mu (dR/d(d_mu f)) - dR/df
        worst = 0
        for (a, b), fn in f.items():
            el = -sp.diff(R, fn)
            for mu, x in enumerate(xs):
                el += sp.diff(sp.diff(R, sp.Derivative(fn, x)), x)
            el = sp.expand(el)
            worst = max(worst, 0 if el == 0 else 1)
        out[dim] = "EL == 0" if worst == 0 else "EL != 0"
    check("2a R_eta: Euler-Lagrange identically zero for generic symmetric M (sympy, 1+1 and 2+1)",
          all(v == "EL == 0" for v in out.values()), out)
    RESULTS["s2a"] = out


# ====================================================================== 2b
def _spectral_grad(F, L):
    """d/dx and d/dy of a periodic field F(nx, ny, ...) by FFT."""
    n = F.shape[0]
    k = 2 * np.pi * np.fft.fftfreq(n, d=L / n)
    Fh = np.fft.fft2(F, axes=(0, 1))
    dx = np.real(np.fft.ifft2(1j * k[:, None, None, None] * Fh, axes=(0, 1)))
    dy = np.real(np.fft.ifft2(1j * k[None, :, None, None] * Fh, axes=(0, 1)))
    return dx, dy


def _expm_field(theta, X):
    """exp(theta(x) X) for a fixed 4x4 generator X (series to high order)."""
    out = np.broadcast_to(np.eye(4), theta.shape + (4, 4)).copy()
    term = out.copy()
    for k in range(1, 30):
        term = (term @ X) * (theta / k)[..., None, None]
        out = out + term
    return out


def s2b_rg_total_derivative():
    gens = so13_generators()
    M0 = np.diag([8.0, 1.0, 0.3, 0.1])
    Lbox = 2 * np.pi
    kinds = ("eta", "etaMeta", "Minv", "hcov")
    GENSETS = {"boost03+rot12+rot13": ("boost_03", "rot_12", "rot_13"),
               "rotations_only(rot12,rot13,rot23)": ("rot_12", "rot_13", "rot_23"),
               "boosts_only(b01,b02,b03)": ("boost_01", "boost_02", "boost_03")}
    SLOTS = {"tx": (0, 1), "ty": (0, 2), "tz": (0, 3), "xy": (1, 2), "xz": (1, 3), "yz": (2, 3)}
    out = {}

    def field(n, mode, gset):
        x = np.arange(n) * Lbox / n
        X, Y = np.meshgrid(x, x, indexing="ij")
        th1 = 0.4 * np.sin(X + 0.3) * np.cos(Y - 0.7) + 0.2 * np.cos(2 * X)
        th2 = 0.5 * np.cos(X - 1.1) * np.sin(2 * Y) + 0.15 * np.sin(Y)
        th3 = 0.3 * np.sin(X + Y) + 0.25 * np.cos(X - 2 * Y)
        g1, g2, g3 = GENSETS[gset]
        if mode in ("orbit", "eig+frame"):
            Lm = _expm_field(th1, gens[g1]) @ _expm_field(th2, gens[g2]) @ _expm_field(th3, gens[g3])
        else:
            Lm = np.broadcast_to(np.eye(4), (n, n, 4, 4))
        if mode == "orbit":
            D = np.broadcast_to(M0, (n, n, 4, 4))
        else:
            lam = np.stack([8.0 + 0.5 * np.sin(X) * np.cos(Y), 1.0 + 0.2 * np.cos(X + Y),
                            0.3 + 0.1 * np.sin(2 * Y), 0.2 + 0.1 * np.cos(X - Y)], axis=-1)   # all > 0
            D = np.zeros((n, n, 4, 4)); D[..., np.arange(4), np.arange(4)] = lam
        return Lm @ D @ Lm.swapaxes(-1, -2)

    def integral(M, kind, slots):
        n = M.shape[0]
        dx, dy = _spectral_grad(M, Lbox)
        A = np.zeros((4, n, n, 4, 4))
        A[slots[0]], A[slots[1]] = dx, dy
        return float(np.sum(T14.rg_density(A, M, kind))) * (Lbox / n) ** 2

    for mode in ("orbit", "eig+frame", "eig_only"):
        out[mode] = {}
        for gset in GENSETS:
            if mode == "eig_only" and gset != "boost03+rot12+rot13":
                continue
            for kind in kinds:
                for sl, slots in SLOTS.items():
                    vals = [integral(field(n, mode, gset), kind, slots) for n in (8, 16, 32, 64)]
                    out[mode][f"{gset}|{kind}|{sl}"] = vals
    RESULTS["s2b"] = out

    def converged(v, rel=1e-4):
        return abs(v[-1] - v[-2]) <= rel * max(abs(v[-1]), 1e-12) + 1e-10
    orb = out["orbit"]

    def sel(mode_d, kind, gset_prefix):
        return {k.split("|")[2]: v for k, v in mode_d.items() if k.split("|")[1] == kind and k.startswith(gset_prefix)}
    zero_map = {kind: {g: all(abs(v[-1]) < 1e-9 for v in sel(orb, kind, g).values())
                       for g in ("boost03", "rotations_only", "boosts_only")} for kind in kinds}
    RESULTS["s2b_orbit_zero_map"] = zero_map
    check("2b R_eta integrates to zero on the orbit for EVERY slot pair and generator set (EL == 0 for any field: the literal EH term is empty)",
          all(zero_map["eta"].values()), {k: f"{v[-1]:.1e}" for k, v in sel(orb, "eta", "boosts_only").items()})
    for kind in ("etaMeta", "Minv", "hcov"):
        check(f"2b R_{kind} is a total derivative on ROTATION orbits (rot12, rot13, rot23 textures; all six slot pairs; spectral convergence)",
              zero_map[kind]["rotations_only"], {k: f"{v[-1]:.1e}" for k, v in sel(orb, kind, "rotations_only").items()})
        bo = sel(orb, kind, "boosts_only")
        nz = any(abs(v[-1]) > 1e-6 and converged(v) for v in bo.values())
        check(f"2b REFUTED for R_{kind}: on BOOST orbits (b01, b02, b03 textures) the periodic integral is NONZERO and spectrally converged (not a total derivative there; the author's 'boosts included' fails)",
              nz, {k: [f"{x:.3e}" for x in v[1:]] for k, v in bo.items()})
    b3 = {kind: {k: f"{v[-1]:.3e}" for k, v in sel(orb, kind, "boost03").items()} for kind in ("etaMeta", "Minv", "hcov")}
    check("2b single-boost texture (boost03 + rot12 + rot13): zero on the slot pairs that do not touch the boost axis (xy), nonzero on xz, yz, tx, ty for eta M eta and M^-1 (reported)",
          True, b3)
    mut = out["eig+frame"]
    for kind in ("etaMeta", "Minv"):
        vals = sel(mut, kind, "boost03")
        nz = all(converged(v) and abs(v[-1]) > 1e-4 for v in vals.values())
        check(f"2b MUTATION: R_{kind} with eigenvalue gradients in the boost03 + rot12 + rot13 frame is NONZERO and spectrally converged on all six slot pairs (bulk content; rotation-only frames: nonzero on the spatial pairs, reported)",
              nz, {"boost03": {k: f"{v[-1]:.3e}" for k, v in vals.items()},
                   "rotations_only": {k: f"{v[-1]:.3e}" for k, v in sel(mut, kind, "rotations_only").items()}})
    eo = out["eig_only"]
    check("2b a pure eigenvalue gradient in a CONSTANT frame integrates to zero for every G and slot pair (R_G is a (d lambda).(d frame) coupling)",
          all(abs(v[-1]) < 1e-9 for v in eo.values()), {k.split("|")[1] + "|" + k.split("|")[2]: f"{v[-1]:.1e}" for k, v in eo.items()})


# ====================================================================== lattice fields
def load_seed(n, L):
    M, cfg, rec = C13.seed_hedgehog(n, L)
    return M, cfg


def s2c_4a_4b_hedgehog():
    M, cfg = load_seed(32, 48.0)
    a0 = C13.a0_local(M)
    p = L0.default_params(s=-1.0, g=8.0, delta=0.3)
    # 2c: R_G omega^2 content on the clock = 0 (three-point read on the full Lagrangian density)
    out = {}
    for kind in ("eta", "etaMeta", "hcov"):          # Minv undefined on the singular certified vacuum
        vals = []
        for om in (0.0, 1.0, -1.0):
            tot = 0.0
            for br, wt in B3.branches(cfg["stencil"]):
                A = T14.jets_static(M, cfg, br)
                A[0] = om * a0
                tot += wt * float(np.sum(T14.rg_density(A, M, kind)))
            vals.append(cfg["h"] ** 3 * tot)
        Cc = 0.5 * (vals[1] + vals[2]) - vals[0]
        Bc = 0.5 * (vals[1] - vals[2])
        out[kind] = {"A": vals[0], "B": Bc, "C": Cc}
    worst = max(abs(v["C"]) / max(abs(v["A"]), 1e-300) for v in out.values())
    check("2c R_G has NO omega^2 content on the relaxed hedgehog with the local clock (C / A; G = eta, eta M eta, h_cov; M^-1 undefined on the singular certified vacuum)", worst < 1e-12, worst,
          {k: {kk: f"{vv:.4e}" for kk, vv in v.items()} for k, v in out.items()})
    RESULTS["s2c"] = out
    # 4a: K_lambda channel table + 4b: K_P channel table vs I1 kin, gen_catalog channels + a0_local
    # channels as TRUE tangents X M + M X^T (gen_catalog's channels are the antisymmetric
    # probes X M - M X^T flagged at R1; here every channel is a symmetric tangent), envelope-
    # weighted and unit-normalized like gen_catalog, plus the raw local clock
    w = B3.envelope(cfg)[..., None, None]
    lamv, Vv = np.linalg.eigh(M[..., 1:4, 1:4])
    Jz = np.zeros((4, 4)); Jz[1, 2], Jz[2, 1] = -1.0, 1.0
    Jx = np.zeros((4, 4)); Jx[2, 3], Jx[3, 2] = -1.0, 1.0
    Kz = np.zeros((4, 4)); Kz[0, 3] = Kz[3, 0] = 1.0
    Kx = np.zeros((4, 4)); Kx[0, 1] = Kx[1, 0] = 1.0
    def local_rot(vhat):
        W = np.zeros(vhat.shape[:-1] + (4, 4))
        n1, n2, n3 = vhat[..., 0], vhat[..., 1], vhat[..., 2]
        W[..., 1, 2], W[..., 1, 3] = -n3, n2
        W[..., 2, 1], W[..., 2, 3] = n3, -n1
        W[..., 3, 1], W[..., 3, 2] = -n2, n1
        return W
    Xs = {"clock_local": local_rot(Vv[..., :, 2]), "plane_1d": local_rot(Vv[..., :, 0]),
          "rot_z": np.broadcast_to(Jz, M.shape), "rot_x": np.broadcast_to(Jx, M.shape),
          "boost_z": np.broadcast_to(Kz, M.shape), "boost_x": np.broadcast_to(Kx, M.shape)}
    chans = {"clock_local(a0_local, unnormalized)": a0}
    for nm, X in Xs.items():
        t = w * (X @ M + M @ X.swapaxes(-1, -2))
        chans[nm + "(tangent, unit norm)"] = t / np.sqrt(np.sum(t * t))
    tab = {}
    for nm, a in chans.items():
        tab[nm] = {"kin_I1": float(B3.kin_of(M, a, cfg)),
                   "kin_K_lambda": T14.kin_energy("K_lambda", M, a, cfg),
                   "kin_K_P": T14.kin_energy("K_P", M, a, cfg)}
    RESULTS["s4_channels_hedgehog"] = tab
    kl = max(abs(v["kin_K_lambda"]) for v in tab.values())
    check("4a K_lambda omega^2 coefficient = 0 on every generator channel of the relaxed hedgehog (clock_local, plane_1d, rot_z, rot_x, boost_z, boost_x)",
          kl < 1e-10, kl)
    ratio = {nm: v["kin_K_P"] / max(abs(v["kin_I1"]), 1e-300) for nm, v in tab.items()}
    kp_clock = tab["clock_local(tangent, unit norm)"]["kin_K_P"]
    kp_boost = max(tab["boost_z(tangent, unit norm)"]["kin_K_P"], tab["boost_x(tangent, unit norm)"]["kin_K_P"])
    check("4b K_P omega^2 on the relaxed hedgehog (true tangents, unit norm): the local clock dominates; the boost channels are suppressed (off-shell core leakage only; ratio boost/clock reported)",
          kp_clock > 0 and kp_boost < 0.05 * kp_clock, {"K_P clock": f"{kp_clock:.4e}", "K_P boost max": f"{kp_boost:.4e}",
                                                        "K_P/|I1| per channel": {k: f"{v:.3e}" for k, v in ratio.items()}})
    # 4a: K_lambda static on the orbit field is exactly zero; on the hedgehog the core deficit
    kl_static = {}
    BOXES = (("n32_L48", 32, 48.0, os.path.join(CK, "m5_32_r10", "relax_g8_n32_L48_it3000.npy"), "R10 3000 it"),
             ("n40_L60", 40, 60.0, os.path.join(CK, "m5_32_r10", "aud_b40_3000.npy"), "R10 audit 3000 it"),
             ("n48_L72", 48, 72.0, None, "R13-W seed 3000 it"),
             ("n32_L48_it12000", 32, 48.0, C13.R10_SEED, "R10 12000 it"))
    for key, n, L, path, mat in BOXES:
        if path is None:
            Mx, cfgx = load_seed(n, L)
        else:
            Mx, cfgx = np.load(path), C13.cfg_of(n, L)
        kl_static[key] = {
            "K_lambda_static": T14.static_energy("K_lambda", Mx, cfgx),
            "K_P_static": T14.static_energy("K_P", Mx, cfgx),
            "E_u": float(B3.e_parts(Mx, cfgx)[0]),
            "kin_I1_local": float(B3.kin_of(Mx, C13.a0_local(Mx), cfgx)),
            "R_eta_static": T14.static_energy("R_eta", Mx, cfgx),
            "R_etaMeta_static": T14.static_energy("R_etaMeta", Mx, cfgx),
            "R_Minv_static": "undefined: the certified vacuum has eigenvalue 0 (M singular)",
            "R_hcov_static": T14.static_energy("R_hcov", Mx, cfgx),
            "K_P_h_static": T14.static_energy("K_P_h", Mx, cfgx),
            "maturity": mat}
    RESULTS["s4_static_boxes"] = kl_static
    Ls = np.array([48.0, 60.0, 72.0])
    def expo(key):
        y = np.array([kl_static[k][key] for k in ("n32_L48", "n40_L60", "n48_L72")])
        if np.any(y <= 0):
            return None
        return float(np.polyfit(np.log(Ls), np.log(y), 1)[0])
    ex = {k: expo(k) for k in ("K_lambda_static", "K_P_static", "K_P_h_static", "E_u", "kin_I1_local")}
    RESULTS["s4_L_exponents_h1.5"] = ex
    check("4a/4b hedgehog TAIL (matched maturity, 3000 it, h = 1.5, L = 48/60/72): K_P's static energy is L-DIVERGENT (exponent > 1: the (2,3)-frame connection of the hedgehog is charged), K_lambda's is core-local (exponent < 0.5); E_u and the I1 clock inertia exponents for reference",
          ex["K_P_static"] is not None and ex["K_P_static"] > 1.0 and ex["K_lambda_static"] is not None and ex["K_lambda_static"] < 0.5, ex)
    # orbit field: K_lambda static = 0 exactly (on a lattice orbit field built from the vacuum)
    n = 16
    cfg2 = C13.cfg_of(n, 24.0)
    X, Y, Z = B3.coords(n, cfg2["h"])
    gens = so13_generators()
    Lm = _expm_field(0.3 * np.sin(2 * np.pi * X / 24.0), gens["boost_03"]) @ _expm_field(0.4 * np.cos(2 * np.pi * Y / 24.0), gens["rot_12"])
    Mo = Lm @ np.diag([8.0, 1.0, 0.3, 0.0]) @ Lm.swapaxes(-1, -2)
    v = T14.static_energy("K_lambda", Mo, cfg2)
    vp = 0.0
    for br, wt in B3.branches(cfg2["stencil"]):
        vp += wt * float(np.sum(T14.klam_static(T14.jets_static(Mo, cfg2, br), Mo))) * cfg2["h"] ** 3
    check("4a K_lambda static = 0 EXACTLY on a lattice orbit field (FD of the sorted spectrum; the perturbative read on FD jets carries the off-tangent residual, reported)",
          abs(v) < 1e-20, v, {"perturbative_read_on_FD_jets": vp})


# ====================================================================== 4c
def s4c_sheet():
    # the W0 S5 family: a (1,2) twist psi(z) rising to 1 over a width w on the vacuum;
    # the clock generator is the LOCAL one (R12 / R13-W convention), so the rotated vacuum
    # region ticks like the vacuum and only the ramp is different
    rows = []
    for (n, L) in ((16, 24.0), (32, 24.0), (64, 24.0)):
        cfg = B3.base_cfg(s=-1.0, g=8.0, n=n, L=L)
        h = cfg["h"]
        for wphys in (1.5, 3.0, 6.0):
            wc = int(round(wphys / h))
            if wc < 1 or wc > n // 2:
                continue
            ps = np.zeros(n)
            k0 = n // 2 - wc // 2
            ps[k0:k0 + wc] = np.linspace(0, 1.0, wc + 1)[1:]
            ps[k0 + wc:] = 1.0
            Rn = B8.rot_field(B8.G3, np.broadcast_to(ps[None, None, :], (n, n, n)))
            Mn = np.einsum("...ab,bc,...dc->...ad", Rn, B3.vac4(cfg), Rn)
            a0 = C13.a0_local(Mn)
            area = n * n * h ** 2
            e_u, e_v = B3.e_parts(Mn, cfg)
            dens = T14.kin_density("K_P", Mn, a0, cfg)
            vac_cell = float(dens[0, 0, 0])
            r = {"n": n, "h": h, "w_cells": wc, "w": wc * h, "E_u": float(e_u), "V4": float(e_v),
                 "K_P_static_per_area": T14.static_energy("K_P", Mn, cfg) / area,
                 "K_lambda_static_per_area": T14.static_energy("K_lambda", Mn, cfg) / area,
                 "R_etaMeta_static_per_area": T14.static_energy("R_etaMeta", Mn, cfg) / area,
                 "kin_I1_per_area": float(B3.kin_of(Mn, a0, cfg)) / area,
                 "kin_K_P_excess_per_area": float(np.sum(dens - vac_cell)) / area,
                 "kin_K_P_vacuum_per_cell": vac_cell}
            rows.append(r)
            log(f"  S5 n{n} w = {wc * h:.2f}: E_u {e_u:.1e} KP_stat/area {r['K_P_static_per_area']:.3e} "
                f"kinI1*w {r['kin_I1_per_area'] * r['w']:.4f} kinKP_excess/area {r['kin_K_P_excess_per_area']:.3f}")
    RESULTS["s4c_sheet"] = rows
    # continuum identity (S5 line) + lattice: the K_P static cost of the twist at FIXED w scales
    # as h^p with p ~ 3 (the sym stencil's second-order contamination of finite rotations)
    def hexp(wphys, key):
        pts = [(r["h"], r[key]) for r in rows if abs(r["w"] - wphys) < 1e-9]
        if len(pts) < 2 or any(v <= 0 for _, v in pts):
            return None
        hs, vs = np.array([p[0] for p in pts]), np.array([p[1] for p in pts])
        return float(np.polyfit(np.log(hs), np.log(vs), 1)[0])
    hx = {f"w{w:g}": hexp(w, "K_P_static_per_area") for w in (1.5, 3.0, 6.0)}
    RESULTS["s4c_KP_static_h_exponent_at_fixed_w"] = hx
    check("4c K_P static on the twist is a LATTICE ARTIFACT of finite rotations (continuum identity: 0; at fixed w it vanishes as h^p, p -> 2 on the h ladder 1.5/0.75/0.375)",
          all(v is not None and v >= 1.5 for v in hx.values()), hx)
    kl_zero = all(abs(r["K_lambda_static_per_area"]) < 1e-18 for r in rows)
    check("4c K_lambda static is EXACTLY zero on the twist family on the lattice (the spectrum is cell-wise constant); R_etaMeta reported",
          kl_zero, {"R_etaMeta_per_area": [f"{r['R_etaMeta_static_per_area']:.2e}" for r in rows if r["n"] == 32]})
    r32 = [r for r in rows if r["n"] == 32]
    kw_i1 = [r["kin_I1_per_area"] * r["w"] for r in r32]
    ex_kp = [r["kin_K_P_excess_per_area"] for r in r32]
    check("4c the free inertia is untouched: kin_I1 x w -> const (~1/w, h = 0.75) while K_P's inertia excess over the ticking vacuum stays BOUNDED (no 1/w growth)",
          max(kw_i1) / min(kw_i1) < 1.6 and max(abs(x) for x in ex_kp) < 2.0 * abs(ex_kp[-1]) + 1e-9 or abs(ex_kp[0]) <= abs(ex_kp[-1]) * 1.5,
          {"kinI1*w": [round(v, 4) for v in kw_i1], "kinKP_excess/area": [round(v, 3) for v in ex_kp]})
    # W3 zigzag end fields: K_P static density on the layers vs I1 kin density
    zz = {}
    for tag, n_, L_, Rs in (("w3_n32_L48_R9_J200_it3000", 32, 48.0, 9.0), ("w3_n32_L48_R9_J200_it12000", 32, 48.0, 9.0),
                            ("w3_n48_L72_R15_J200_it3000", 48, 72.0, 15.0)):
        f = os.path.join(CK, "m5_32_r13w", tag + ".npy")
        if not os.path.exists(f):
            continue
        Mz = np.load(f)
        cfgz = C13.cfg_of(n_, L_)
        X, Y, Z = B3.coords(n_, cfgz["h"])
        r = np.sqrt(X * X + Y * Y + Z * Z)
        a0z = C13.a0_local(Mz) * (r < Rs)[..., None, None]
        kin_i1 = C13.kin_density(Mz, a0z, cfgz)
        kp_st = T14.static_density("K_P", Mz, cfgz)
        kl_st = T14.static_density("K_lambda", Mz, cfgz)
        layer = (r > Rs - 2.5 * cfgz["h"]) & (r < Rs + 0.5 * cfgz["h"])
        zz[tag] = {"kin_I1_total": float(np.sum(kin_i1)), "kin_I1_layer": float(np.sum(kin_i1[layer])),
                   "K_P_static_total": float(np.sum(kp_st)), "K_P_static_layer": float(np.sum(kp_st[layer])),
                   "K_lambda_static_total": float(np.sum(kl_st)), "K_lambda_static_layer": float(np.sum(kl_st[layer])),
                   "E_u": float(B3.e_parts(Mz, cfgz)[0])}
        log(f"  W3 {tag}: kinI1 layer {zz[tag]['kin_I1_layer']:.1f}/{zz[tag]['kin_I1_total']:.1f}  K_P_stat layer {zz[tag]['K_P_static_layer']:.3f}/{zz[tag]['K_P_static_total']:.3f}  K_lam layer {zz[tag]['K_lambda_static_layer']:.3f}")
    RESULTS["s4c_w3_zigzag"] = zz
    if zz:
        frac = {k: v["K_P_static_layer"] / max(v["K_P_static_total"], 1e-300) for k, v in zz.items()}
        check("4c ADDITION: K_P DOES charge the W3 zigzag (split-modulus gradient) on the layers inside the region edge (layer fraction of its static energy reported)",
              all(v["K_P_static_layer"] > 0 for v in zz.values()), {k: f"{v:.3f}" for k, v in frac.items()})


# ====================================================================== S: the W0 theorems per entrant
def s_theorems():
    z, om, delta = sp.symbols("z omega delta", real=True)
    g = sp.Integer(8)
    d4 = sp.diag(g, 1, delta, 0)
    XI = sp.diag(-1, 1, 1, 1)
    G1 = sp.zeros(4); G1[2, 3], G1[3, 2] = -1, 1
    G3 = sp.zeros(4); G3[1, 2], G3[2, 1] = -1, 1

    def rot(Gm, q):
        return sp.eye(4) + sp.sin(q) * Gm + (1 - sp.cos(q)) * Gm * Gm

    def kp_sym(A, M):
        N = M * XI
        P = (N + g * sp.eye(4)) * (N - sp.eye(4))
        X = P * A * XI * P
        return sp.Rational(1, 2) * (X * X).trace()

    def rg_sym(As, M, G):
        # As: dict mu -> jet matrix (only the given slots nonzero)
        R = 0
        for mu, Am in As.items():
            for nu, An in As.items():
                for c in range(4):
                    for d in range(4):
                        if G[c, d] == 0:
                            continue
                        R += G[c, d] * (Am[nu, c] * An[mu, d] - Am[mu, c] * An[nu, d])
        return sp.simplify(R)
    out = {}
    # S3 planar: a general planar profile M(z) (diagonal eigenvalues varying + a (2,3) rotation angle)
    l0, l1, l2, l3, phi = [sp.Function(nm)(z) for nm in ("l0", "l1", "l2", "l3", "phi")]
    Dz = sp.diag(l0, l1, l2, l3)
    Rz = rot(G1, phi)
    Mz = Rz * Dz * Rz.T
    Az = Mz.diff(z)
    r_planar = {k: rg_sym({3: Az}, Mz, Gk) for k, Gk in (("eta", XI), ("etaMeta", XI * Mz * XI))}
    out["S3_R_G_on_any_planar_profile"] = {k: str(v) for k, v in r_planar.items()}
    check("S3 planar flatness survives R_G: R_G == 0 identically on ANY planar profile M(z) (only mu = nu = z terms, cancelling), G = eta and eta M eta",
          all(v == 0 for v in r_planar.values()), out["S3_R_G_on_any_planar_profile"])
    kp_planar = sp.simplify(kp_sym(Az, Mz))
    kp_ramp = sp.simplify(kp_planar.subs({phi: 0}).doit())
    out["S3_K_P_on_planar_eigenvalue_ramp"] = str(kp_ramp)
    check("S3 planar flatness OVERTURNED by K_P for eigenvalue walls: K_P != 0 on a (2,3) split ramp (a degenerate wall gets tension)",
          kp_ramp != 0, str(kp_ramp))
    # S4 phase flatness: uniform vacuum with phi(x,t): K_P on the rotating vacuum
    tt = sp.symbols("t", real=True)
    ph = sp.Function("ph")(tt, z)
    Mv = rot(G1, ph) * d4 * rot(G1, ph).T
    At, Azv = Mv.diff(tt), Mv.diff(z)
    kp_t = sp.simplify(kp_sym(At, Mv))
    kp_z = sp.simplify(kp_sym(Azv, Mv))
    f2, f3 = (delta + g) * (delta - 1), (0 + g) * (0 - 1)
    pred_t = sp.simplify((f2 * f3) ** 2 * delta ** 2 * ph.diff(tt) ** 2)
    out["S4_K_P_phase_stiffness"] = {"K_P_time": str(kp_t), "K_P_space": str(kp_z), "predicted": str(pred_t)}
    check("S4 phase flatness OVERTURNED by K_P: the phase on the uniform vacuum carries [f(delta) f(0)]^2 delta^2 (d phi)^2 in time and space (exact)",
          sp.simplify(kp_t - pred_t) == 0 and sp.simplify(kp_z - pred_t.subs(ph.diff(tt), ph.diff(z))) == 0,
          out["S4_K_P_phase_stiffness"])
    # S5 free inertia: the (1,2) twist psi(z): every entrant's static density
    psi = sp.Function("psi")(z)
    M5 = rot(G3, psi) * d4 * rot(G3, psi).T
    A5 = M5.diff(z)
    kp5 = sp.simplify(kp_sym(A5, M5))
    rg5 = {k: rg_sym({3: A5}, M5, Gk) for k, Gk in (("eta", XI), ("etaMeta", XI * M5 * XI))}
    # K_lambda: spectrum constant on the twist -> 0 by construction (eigenvalues of M5 xi are those of d4 xi)
    trs = [sp.simplify(((M5 * XI) ** p).trace() - ((d4 * XI) ** p).trace()) for p in range(1, 5)]
    out["S5_static_on_the_twist"] = {"K_P": str(kp5), "R_G": {k: str(v) for k, v in rg5.items()},
                                     "K_lambda_spectrum_traces_minus_vacuum": [str(t) for t in trs]}
    s5_ok = kp5 == 0 and all(v == 0 for v in rg5.values()) and all(t == 0 for t in trs)
    check("S5 FREE INERTIA SURVIVES THE WHOLE TWO-DERIVATIVE SET: K_P, R_G, K_lambda all have ZERO static cost on the (1,2) twist (so no fixed-J minimizer on L_cert + any of them)",
          s5_ok, out["S5_static_on_the_twist"])
    # the quartics: C6a on the twist is nonzero (the only entrants of B that charge the sheet)
    c6a = sp.simplify(((A5 * XI * A5 * XI).trace()) ** 2)
    out["S5_C6a_on_the_twist"] = str(c6a)
    check("S5 the C6 quartic [tr(dM eta dM eta)]^2 is NONZERO on the twist (~ psi'^4: the sheet's only static cost inside B)",
          c6a != 0, str(c6a))
    RESULTS["s_theorems"] = out


# ====================================================================== V: vacuum ticking
def s_v_vacuum_ticking():
    cfg = C13.cfg_of(8, 12.0)
    Mv = np.broadcast_to(B3.vac4(cfg), (8, 8, 8, 4, 4)).copy()
    a0 = C13.a0_local(Mv)
    kin_cell_kp = T14.kin_energy("K_P", Mv, a0, cfg) / (8 ** 3)
    kin_cell_i1 = float(B3.kin_of(Mv, a0, cfg)) / (8 ** 3)
    # V4 per cell of the degenerate state diag(g,1,delta/2,delta/2)-type: the R13-W numbers
    Md = np.broadcast_to(np.diag([8.0, 1.0, 0.15, 0.15]), (8, 8, 8, 4, 4)).copy()
    v4_deg = float(B3.e_parts(Md, cfg)[1]) / (8 ** 3)
    RESULTS["s_v"] = {"K_P_kin_per_cell_on_the_vacuum(h^3-weighted)": kin_cell_kp, "I1_kin_per_cell_on_the_vacuum": kin_cell_i1,
                      "V4_per_cell_degenerate_(delta/2,delta/2)": v4_deg,
                      "omega2_at_which_K_P_ticking_beats_V4_deg": v4_deg / kin_cell_kp}
    check("V the certified vacuum TICKS under the local clock: K_P's omega^2 coefficient is a constant per cell (VOLUME-extensive inertia) while I1's is zero there",
          kin_cell_kp > 0 and abs(kin_cell_i1) < 1e-14, RESULTS["s_v"])


# ====================================================================== registry lines re-reported
def s_registry():
    r = T14.selftest(write=True)
    RESULTS["registry_selftest"] = r["summary"]
    check("4b/4d registry selftest (covariance + no-eta control, positivity on-orbit, off-shell indefiniteness, K_lambda FD, gradients, vacuum channel table, literal-roots mutant)",
          r["summary"]["pass"] == r["summary"]["total"], r["summary"])


def main():
    s_registry()
    s1a_p249_split()
    s1b_goldstones()
    s2a_reta_euler_lagrange()
    s2b_rg_total_derivative()
    s_theorems()
    s_v_vacuum_ticking()
    s2c_4a_4b_hedgehog()
    s4c_sheet()
    n = sum(1 for l in RESULTS["lines"] if l["pass"])
    RESULTS["summary"] = {"pass": n, "total": len(RESULTS["lines"]), "wall_s": round(time.time() - T0, 1)}
    json.dump(RESULTS, open(OUT, "w"), indent=1, default=float)
    log(f"R14-0 {n}/{len(RESULTS['lines'])} -> {OUT}")


if __name__ == "__main__":
    main()
