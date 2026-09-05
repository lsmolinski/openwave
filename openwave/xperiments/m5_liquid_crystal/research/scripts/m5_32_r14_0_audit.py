"""M5.32 R14-0 adversarial audit: an INDEPENDENT re-derivation of the 14 R14-0 claims.

Independent of the producer's code (m5_32_r14_terms.py, m5_32_r14_0_verify.py and
their JSON were never opened): own finite-difference stencils, own sympy, own
eta-orthonormal eigen-frame, own periodic-box textures, own twist family and h ladder.
The certified stack (m5_21_3_a_4d.py) is imported ONLY to cross-check the audit's own
stencil (E_u, kin) and for coords / cfg; m5_32_r13w_common.a0_local is cross-checked
against the audit's own local generator, never used as the primary instrument.

Conventions (the R14-0 brief): M symmetric 4x4, eta = diag(-1,1,1,1), N = M eta,
vacuum M = diag(8, 1, 0.3, 0) (N spectrum (-8, 1, 0.3, 0)); A_mu = d_mu M; sym
stencil = the average of the fwd and bwd branch densities; h^3 weights.

Entrants (E-densities):
  K_lambda  (1/2) sum_a [sum_i (d_i lam_a)^2 + om^2 (d_t lam_a)^2], lam_a eigenvalues of N
  R_G       sum_{mu != nu} G_cd [A_mu[nu,c] A_nu[mu,d] - A_mu[mu,c] A_nu[nu,d]]
  K_P       (1/2)[sum_i tr((P A_i eta P)^2) + om^2 tr((P a0 eta P)^2)], P = (N + g)(N - 1)
  K_P^h     the Frobenius norm of Om = P A eta P in the eta-orthonormal eigenbasis of N

Stages: c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 c11 c12 c13 c14 | all
Out: ../data/m5_32_r14_0_audit.json (merged per stage)
"""
from __future__ import annotations

import importlib.util
import json
import os
import sys
import time

import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.dirname(HERE)
DATA = os.path.join(RES, "data")
CK10 = os.path.join(RES, "checkpoints", "m5_32_r10")
CK13 = os.path.join(RES, "checkpoints", "m5_32_r13w")
OUT = os.path.join(DATA, "m5_32_r14_0_audit.json")
T0 = time.time()

ETA = np.diag([-1.0, 1.0, 1.0, 1.0])
I4 = np.eye(4)
G_, DELTA = 8.0, 0.3
M_VAC = np.diag([8.0, 1.0, 0.3, 0.0])
M_DEG = np.diag([8.0, 1.0, 0.3, 0.3])
W1 = 0.000724023879          # the certified V4 weight (m5_21_3_a_4d.W1)


def log(m):
    print(f"[{time.time() - T0:7.1f}s] {m}", flush=True)


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


INS4 = _load("m5_21_3_a_4d", "m5_21_3_a_4d.py")          # cross-checks only
R13W = _load("m5_32_r13w_common", "m5_32_r13w_common.py")  # cfg_of + a0_local cross-check


def save(stage, rec):
    db = json.load(open(OUT)) if os.path.exists(OUT) else {}
    db[stage] = rec
    db["_meta"] = {"script": "m5_32_r14_0_audit.py", "python": sys.version.split()[0],
                   "numpy": np.__version__, "sympy": sp.__version__,
                   "updated": time.strftime("%Y-%m-%d %H:%M:%S")}
    json.dump(db, open(OUT, "w"), indent=1, default=float)


# ============================================================ algebra
def gen_rot(i, j):
    """rotation generator in the (i, j) plane: G[i, j] = -1, G[j, i] = +1 (G1 = gen_rot(2, 3))."""
    X = np.zeros((4, 4)); X[i, j], X[j, i] = -1.0, 1.0
    return X


def gen_boost(i):
    X = np.zeros((4, 4)); X[0, i] = X[i, 0] = 1.0
    return X


GENS = {"boost_01": gen_boost(1), "boost_02": gen_boost(2), "boost_03": gen_boost(3),
        "rot_12": gen_rot(1, 2), "rot_13": gen_rot(1, 3), "rot_23": gen_rot(2, 3)}
G1, G3 = GENS["rot_23"], GENS["rot_12"]


def tangent(X, M):
    """orbit tangent X M + M X^T (rotations: X M - M X)."""
    return X @ M + M @ X.T


def rot_of(X, th):
    th = np.asarray(th, dtype=float)
    return (I4 + np.sin(th)[..., None, None] * X + (1.0 - np.cos(th))[..., None, None] * (X @ X))


def boost_of(K, b):
    b = np.asarray(b, dtype=float)
    return (I4 + np.sinh(b)[..., None, None] * K + (np.cosh(b) - 1.0)[..., None, None] * (K @ K))


def conj(L, M):
    """L M L^T per cell."""
    return L @ M @ np.swapaxes(L, -1, -2)


def sym_basis():
    """the 10 symmetric unit directions E_ab (a <= b), unit Frobenius norm."""
    B, lab = [], []
    for a in range(4):
        for b in range(a, 4):
            E = np.zeros((4, 4)); E[a, b] = E[b, a] = 1.0
            B.append(E / np.sqrt(np.sum(E * E))); lab.append(f"M{a}{b}")
    return B, lab


def eig_frame(M):
    """eigen-decomposition of N = M eta per cell, eigenvalues ascending (index 0 = the
    timelike one), columns normalized to v^T eta v = -1 (timelike) / +1 (spacelike)."""
    N = M @ ETA
    w, V = np.linalg.eig(N)
    idx = np.argsort(w.real, axis=-1)
    w = np.take_along_axis(w, idx, -1)
    V = np.take_along_axis(V, idx[..., None, :], -1)
    imag = float(np.max(np.abs(w.imag))) if w.size else 0.0
    w, V = w.real, V.real
    s = np.einsum("...ia,ij,...ja->...a", V, ETA, V)
    V = V / np.sqrt(np.abs(s))[..., None, :]
    return w, V, np.sign(s), imag


def timelike_u(M):
    w, V, s, _ = eig_frame(M)
    return V[..., :, 0], s[..., 0]


def G_of(name, M):
    if name == "eta":
        return np.broadcast_to(ETA, M.shape)
    if name == "eta_M_eta":
        return ETA @ M @ ETA
    if name == "Minv":
        return np.linalg.inv(M)
    if name == "h_cov":
        u, s = timelike_u(M)
        assert np.all(s < 0), "timelike eigenvector not timelike"
        eu = (ETA @ u[..., None])
        return ETA + 2.0 * eu @ np.swapaxes(eu, -1, -2)
    if name == "identity":              # non-covariant control
        return np.broadcast_to(I4, M.shape)
    raise ValueError(name)


# ============================================================ own stencils
def _sl(nd, ax, s):
    t = [slice(None)] * nd; t[ax] = s; return tuple(t)


def dfwd(f, ax, h):
    out = np.zeros_like(f)
    out[_sl(f.ndim, ax, slice(0, -1))] = (f[_sl(f.ndim, ax, slice(1, None))]
                                          - f[_sl(f.ndim, ax, slice(0, -1))]) / h
    return out


def dbwd(f, ax, h):
    out = np.zeros_like(f)
    out[_sl(f.ndim, ax, slice(1, None))] = (f[_sl(f.ndim, ax, slice(1, None))]
                                            - f[_sl(f.ndim, ax, slice(0, -1))]) / h
    return out


def dctr(f, ax, h):
    return 0.5 * (dfwd(f, ax, h) + dbwd(f, ax, h))


def jets(M, h, naxes=3, st="sym"):
    """[(A_list, weight)]; A_list[i] = d_i M along array axis i."""
    if st == "sym":
        return [([dfwd(M, ax, h) for ax in range(naxes)], 0.5),
                ([dbwd(M, ax, h) for ax in range(naxes)], 0.5)]
    if st == "ctr":
        return [([dctr(M, ax, h) for ax in range(naxes)], 1.0)]
    raise ValueError(st)


def eta_inner(F, Gm):
    return np.einsum("...ab,...ab,a,b->...", F, Gm, np.diag(ETA), np.diag(ETA))


def comm_eta(A, B):
    return A @ ETA @ B - B @ ETA @ A


# ============================================================ term densities (per cell, unweighted)
def rg_density(A, Gm):
    """R_G = sum_{mu != nu} G_cd [A_mu[nu,c] A_nu[mu,d] - A_mu[mu,c] A_nu[nu,d]]; A: dict slot -> jet."""
    dens = 0.0
    for mu in A:
        for nu in A:
            if mu == nu:
                continue
            t1 = np.einsum("...cd,...c,...d->...", Gm, A[mu][..., nu, :], A[nu][..., mu, :])
            t2 = np.einsum("...cd,...c,...d->...", Gm, A[mu][..., mu, :], A[nu][..., nu, :])
            dens = dens + (t1 - t2)
    return dens


def i1_density(Alist):
    d = 0.0
    n = len(Alist)
    for i in range(n):
        for j in range(i + 1, n):
            F = comm_eta(Alist[i], Alist[j])
            d = d + 4.0 * eta_inner(F, F)
    return d


def kin_density(Alist, a0):
    d = 0.0
    for A in Alist:
        F = comm_eta(a0, A)
        d = d + 4.0 * eta_inner(F, F)
    return d


def P_of(M, p1=-G_, p2=1.0):
    N = M @ ETA
    return (N - p1 * I4) @ (N - p2 * I4)


def kp_plain_density(Alist, P):
    d = 0.0
    for A in Alist:
        Om = P @ A @ ETA @ P
        d = d + 0.5 * np.einsum("...ab,...ba->...", Om, Om)
    return d


def kp_h_density(Alist, M, P=None, variant="eigenbasis"):
    """H-adjoint K_P: variant 'eigenbasis' = ||V^-1 Om V||_F^2 (the brief's parenthetical),
    'stated' = tr(Om H Om^T H^-1) (the brief's formula), 'covariant' = tr(Om^T H Om H^-1)."""
    P = P_of(M) if P is None else P
    w, V, s, _ = eig_frame(M)
    u = V[..., :, 0]
    eu = ETA @ u[..., None]
    H = ETA + 2.0 * eu @ np.swapaxes(eu, -1, -2)
    Hi = ETA + 2.0 * u[..., None] @ u[..., None, :]
    Vi = np.linalg.inv(V)
    d = 0.0
    for A in Alist:
        Om = P @ A @ ETA @ P
        if variant == "eigenbasis":
            Op = Vi @ Om @ V
            d = d + 0.5 * np.einsum("...ab,...ab->...", Op, Op)
        elif variant == "stated":
            d = d + 0.5 * np.einsum("...aa->...", Om @ H @ np.swapaxes(Om, -1, -2) @ Hi)
        elif variant == "covariant":
            d = d + 0.5 * np.einsum("...aa->...", np.swapaxes(Om, -1, -2) @ H @ Om @ Hi)
        else:
            raise ValueError(variant)
    return d


def lam_sorted(M):
    w = np.linalg.eigvals(M @ ETA)
    imag = float(np.max(np.abs(w.imag)))
    return np.sort(w.real, axis=-1), imag


def klam_static_density(M, h, naxes=3, st="sym"):
    lam, imag = lam_sorted(M)
    d = 0.0
    for Al, wt in jets(lam, h, naxes, st):
        for A in Al:
            d = d + wt * 0.5 * np.sum(A * A, axis=-1)
    return d, imag


def klam_omega2_density(M, a0):
    """(1/2) sum_a (d_t lam_a)^2 with d_t lam_a = (v_a^T eta a0 eta v_a)/(v_a^T eta v_a)."""
    w, V, s, _ = eig_frame(M)
    num = np.einsum("...ia,ij,...jk,kl,...la->...a", V, ETA, a0, ETA, V)
    dl = num / s
    return 0.5 * np.sum(dl * dl, axis=-1)


def a0_local_own(M):
    """J M - M J, J = [n]_x the so(3) generator about the leading spatial eigenvector n."""
    w, V = np.linalg.eigh(M[..., 1:, 1:])
    n = V[..., :, -1]
    J = np.zeros(M.shape)
    J[..., 1, 2], J[..., 1, 3] = -n[..., 2], n[..., 1]
    J[..., 2, 1], J[..., 2, 3] = n[..., 2], -n[..., 0]
    J[..., 3, 1], J[..., 3, 2] = -n[..., 1], n[..., 0]
    return J @ M - M @ J


def a0_lowest_own(M):
    w, V = np.linalg.eigh(M[..., 1:, 1:])
    n = V[..., :, 0]
    J = np.zeros(M.shape)
    J[..., 1, 2], J[..., 1, 3] = -n[..., 2], n[..., 1]
    J[..., 2, 1], J[..., 2, 3] = n[..., 2], -n[..., 0]
    J[..., 3, 1], J[..., 3, 2] = -n[..., 1], n[..., 0]
    return J @ M - M @ J


def coords3(n, h):
    x = (np.arange(n) - (n - 1) / 2.0) * h
    return np.meshgrid(x, x, x, indexing="ij")


FIELDS = {
    "n32_L48_it3000": (os.path.join(CK10, "relax_g8_n32_L48_it3000.npy"), 32, 48.0),
    "n32_L48_it12000": (os.path.join(CK10, "relax_g8_n32_L48_it12000.npy"), 32, 48.0),
    "n40_L60_it3000": (os.path.join(CK10, "aud_b40_3000.npy"), 40, 60.0),
    "n48_L72_it3000": (os.path.join(CK13, "seed_n48_L72_it3000.npy"), 48, 72.0),
    "w3_n32_R9_it3000": (os.path.join(CK13, "w3_n32_L48_R9_J200_it3000.npy"), 32, 48.0),
    "w3_n32_R9_it12000": (os.path.join(CK13, "w3_n32_L48_R9_J200_it12000.npy"), 32, 48.0),
    "w3_n48_R15_it3000": (os.path.join(CK13, "w3_n48_L72_R15_J200_it3000.npy"), 48, 72.0),
}


def load_field(key):
    path, n, L = FIELDS[key]
    return np.load(path), n, L, L / n


# ============================================================ C1: the P249 exterior Hessian
def stage_c1():
    a, t, p, u, v, q = X = sp.symbols("a t p u v q", real=True)
    S = sp.Matrix([[1 + a, u, v], [u, t + p, q], [v, q, t - p]])
    tr2, tr3 = (S * S).trace(), (S * S * S).trace()
    V_ldg = -tr2 / 2 - tr3 + tr2 ** 2 + sp.Rational(1, 2)
    V_axis = tr2 - S[0, 0] ** 2
    P3 = sp.diag(0, 1, 1)
    B = P3 * S * P3 - (P3 * S * P3).trace() * P3 / 2
    V_lock = 6 * (B.T * B).trace() / 2
    at0 = {x: 0 for x in X}

    def hess(Vp):
        return sp.hessian(Vp, X).subs(at0)

    def grad(Vp):
        return [sp.diff(Vp, x).subs(at0) for x in X]

    H = {k: hess(Vp) for k, Vp in (("ldg", V_ldg), ("axis", V_axis), ("lock", V_lock))}
    expect = {"ldg": sp.diag(5, 6, 6, 0, 0, 6), "axis": sp.diag(0, 4, 4, 4, 4, 4),
              "lock": sp.diag(0, 0, 12, 0, 0, 12)}
    rec = {"hessians": {k: [[int(H[k][i, j]) for j in range(6)] for i in range(6)] for k in H},
           "hessian_match": {k: bool(H[k] == expect[k]) for k in H},
           "grad_ldg_at_origin": [str(g) for g in grad(V_ldg)],
           "V_ldg_at_origin": str(V_ldg.subs(at0))}
    # Goldstone: V_ldg on the SO(3) orbit of diag(1,0,0) (three Euler angles, symbolic)
    al, be, ga = sp.symbols("alpha beta gamma", real=True)
    Rz = lambda th: sp.Matrix([[sp.cos(th), -sp.sin(th), 0], [sp.sin(th), sp.cos(th), 0], [0, 0, 1]])
    Ry = lambda th: sp.Matrix([[sp.cos(th), 0, sp.sin(th)], [0, 1, 0], [-sp.sin(th), 0, sp.cos(th)]])
    Rm = Rz(al) * Ry(be) * Rz(ga)
    S_orb = Rm * sp.diag(1, 0, 0) * Rm.T
    tr2o, tr3o = (S_orb * S_orb).trace(), (S_orb ** 3).trace()
    V_orb = sp.simplify(-tr2o / 2 - tr3o + tr2o ** 2 + sp.Rational(1, 2))
    rec["V_ldg_on_rotation_orbit"] = str(V_orb)
    # tangents [X, S0] of the orbit in the chart: X_01 -> u, X_02 -> v, X_12 -> 0
    S0 = sp.diag(1, 0, 0)
    tang = {}
    for (i, j) in ((0, 1), (0, 2), (1, 2)):
        Xg = sp.zeros(3, 3); Xg[i, j], Xg[j, i] = -1, 1
        T = Xg * S0 - S0 * Xg
        # chart read: a = T00, t = (T11+T22)/2, p = (T11-T22)/2, u = T01, v = T02, q = T12
        vec = sp.Matrix([T[0, 0], (T[1, 1] + T[2, 2]) / 2, (T[1, 1] - T[2, 2]) / 2, T[0, 1], T[0, 2], T[1, 2]])
        tang[f"X{i}{j}"] = {"chart": [int(c) for c in vec], "H_ldg_times": [int(c) for c in (H["ldg"] * vec)]}
    rec["orbit_tangents"] = tang
    null_ldg = H["ldg"].nullspace()
    rec["H_ldg_nullspace"] = [[int(c) for c in nv] for nv in null_ldg]
    # MUTATION: tr(S^3) coefficient doubled: origin is no longer critical, Hessian changes
    V_mut = -tr2 / 2 - 2 * tr3 + tr2 ** 2 + sp.Rational(1, 2)
    Hm = hess(V_mut)
    rec["mutation_tr3_doubled"] = {"grad_at_origin": [str(g) for g in grad(V_mut)],
                                   "hessian_diag": [int(Hm[i, i]) for i in range(6)],
                                   "matches_claim": bool(Hm == expect["ldg"])}
    ok = all(rec["hessian_match"].values()) and V_orb == 0 and all(
        all(c == 0 for c in tang[k]["H_ldg_times"]) for k in ("X01", "X02"))
    rec["verdict"] = "CONFIRMED" if ok else "REFUTED"
    log(f"C1 {rec['verdict']} hessians {rec['hessian_match']} V_orbit={V_orb} mutant_matches={rec['mutation_tr3_doubled']['matches_claim']}")
    save("c1", rec)


# ============================================================ C2: Goldstone count on the degenerate vacuum
def poly_trace_hessian(coefs, M):
    """exact Hessian (10x10) of V = sum_k c_k tr(N^k), N = M eta, on the symmetric basis:
    d^2 tr N^k [E_I, E_J] = k sum_{j=0}^{k-2} tr(N^j E_I eta N^{k-2-j} E_J eta)."""
    B, lab = sym_basis()
    N = M @ ETA
    kmax = len(coefs) - 1
    pw = [I4]
    for _ in range(kmax):
        pw.append(pw[-1] @ N)
    H = np.zeros((10, 10))
    grad = np.zeros(10)
    for k, c in enumerate(coefs):
        if c == 0 or k < 1:
            continue
        for I_, EI in enumerate(B):
            grad[I_] += c * k * np.trace(pw[k - 1] @ EI @ ETA)
            if k < 2:
                continue
            for J_, EJ in enumerate(B):
                acc = 0.0
                for j in range(k - 1):
                    acc += np.trace(pw[j] @ EI @ ETA @ pw[k - 2 - j] @ EJ @ ETA)
                H[I_, J_] += c * k * acc
    return H, grad, lab


def fd_hessian(Vfun, M, eps=1e-3):
    B, _ = sym_basis()
    H = np.zeros((10, 10))
    for i, Ei in enumerate(B):
        for j, Ej in enumerate(B):
            f = lambda s, t: Vfun(M + s * eps * Ei + t * eps * Ej)
            H[i, j] = (f(1, 1) - f(1, -1) - f(-1, 1) + f(-1, -1)) / (4 * eps * eps)
    return H


def null_count(H, rtol=1e-8):
    w = np.linalg.eigvalsh(0.5 * (H + H.T))
    scale = float(np.max(np.abs(w)))
    return int(np.sum(np.abs(w) < rtol * scale)), w, scale


def charge_of(T, X=G1, nphi=64):
    """dominant Fourier mode of phi -> <T, R(phi) T R(phi)^T>_F (the (2,3)-clock charge)."""
    phis = 2 * np.pi * np.arange(nphi) / nphi
    R = rot_of(X, phis)
    c = np.einsum("ab,...ab->...", T, conj(R, T))
    spec = np.abs(np.fft.rfft(c))
    spec[0] = spec[0] if np.max(spec[1:]) < 1e-12 * max(spec[0], 1e-300) else 0.0
    return int(np.argmax(spec)), float(np.max(spec) / max(np.sum(spec), 1e-300))


def stage_c2():
    rec = {}
    P = np.polynomial.polynomial

    def f_coefs(roots_mult, extra=(1.0,)):
        c = np.array([1.0])
        for r, m in roots_mult:
            for _ in range(m):
                c = P.polymul(c, [-r, 1.0])
        c = P.polymul(c, np.array(extra, dtype=float))
        return c

    cases = {
        "f_producer_deg": (f_coefs([(-G_, 2), (1.0, 2), (DELTA, 2)]), M_DEG),
        "f_own_x2p3_deg": (f_coefs([(-G_, 2), (1.0, 2), (DELTA, 2)], extra=(3.0, 0.0, 1.0)), M_DEG),
        "f_own_shifted_deg": (f_coefs([(-G_, 2), (1.0, 2), (DELTA, 2)], extra=(2.0, 0.1)), M_DEG),
        "MUT_nondegenerate_vac": (f_coefs([(-G_, 2), (1.0, 2), (DELTA, 2), (0.0, 2)]), M_VAC),
        "MUT_simple_root_delta": (f_coefs([(-G_, 2), (1.0, 2), (DELTA, 1)]), M_DEG),
    }
    for key, (c, M) in cases.items():
        H, g, lab = poly_trace_hessian(c, M)
        Vfun = lambda Mx, c=c: float(np.trace(sum(ck * np.linalg.matrix_power(Mx @ ETA, k) for k, ck in enumerate(c))))
        Hfd = fd_hessian(Vfun, M)
        nz, w, sc = null_count(H)
        rec[key] = {"grad_max": float(np.max(np.abs(g))), "n_flat": nz,
                    "hess_eigs": [float(x) for x in w], "fd_vs_exact_rel": float(np.max(np.abs(H - Hfd)) / sc),
                    "hess_psd": bool(np.min(w) > -1e-8 * sc)}
        log(f"C2 {key}: grad_max {rec[key]['grad_max']:.2e} flat {nz} psd {rec[key]['hess_psd']} fd_rel {rec[key]['fd_vs_exact_rel']:.1e}")
    # null space vs orbit tangents (producer f, degenerate vacuum)
    H, g, lab = poly_trace_hessian(cases["f_producer_deg"][0], M_DEG)
    w, U = np.linalg.eigh(H)
    nullU = U[:, np.abs(w) < 1e-8 * np.max(np.abs(w))]
    B, _ = sym_basis()
    Tvecs, names = [], []
    for nm, X in GENS.items():
        T = tangent(X, M_DEG)
        vec = np.array([np.sum(T * E) for E in B])
        Tvecs.append(vec); names.append(nm)
        rec.setdefault("tangent_norms_deg", {})[nm] = float(np.linalg.norm(vec))
    Tm = np.array(Tvecs).T
    rank_T = int(np.linalg.matrix_rank(Tm, tol=1e-10))
    # subspace agreement: project null space onto tangent span
    Q, _ = np.linalg.qr(Tm[:, [i for i, nm in enumerate(names) if nm != "rot_23"]])
    resid = float(np.linalg.norm(nullU - Q @ (Q.T @ nullU)))
    rec["null_vs_tangents"] = {"tangent_rank": rank_T, "null_dim": int(nullU.shape[1]),
                               "null_minus_projection_norm": resid}
    # charges under the (2,3) clock
    ch = {}
    for nm, X in GENS.items():
        T = tangent(X, M_DEG)
        if np.linalg.norm(T) < 1e-12:
            ch[nm] = "null (stabilizer)"
        else:
            ch[nm] = charge_of(T)[0]
    Ep = np.diag([0.0, 0.0, 1.0, -1.0]); Eq = np.zeros((4, 4)); Eq[2, 3] = Eq[3, 2] = 1.0
    ch["p_split"] = charge_of(Ep)[0]; ch["q_split"] = charge_of(Eq)[0]
    # ad_{rot23} spectrum on so(1,3)
    basis = list(GENS.values())
    ad = np.zeros((6, 6))
    for j, Y in enumerate(basis):
        C = G1 @ Y - Y @ G1
        coef = np.linalg.lstsq(np.array([b.ravel() for b in basis]).T, C.ravel(), rcond=None)[0]
        ad[:, j] = coef
    ad_eigs = np.sort(np.abs(np.linalg.eigvals(ad)))
    rec["charges"] = ch
    rec["ad_rot23_abs_spectrum"] = [float(x) for x in ad_eigs]
    # lock mutation: + (M02^2 + M03^2 + M12^2 + M13^2)
    Hl = H.copy()
    for i, lb in enumerate(lab):
        if lb in ("M02", "M03", "M12", "M13"):
            Hl[i, i] += 2.0 * 1.0 * 2.0  # E is unit-normalized (1/sqrt2 per entry): M_ab^2 = (E_ab.M)^2/2 -> d2 = 2*(1/2)*2 ... computed below exactly
    # exact: lock(M) = sum M_ab^2 over the four pairs; along unit E_ab, M_ab = s/sqrt(2): lock = s^2/2 -> second derivative 1
    Hl = H.copy()
    for i, lb in enumerate(lab):
        if lb in ("M02", "M03", "M12", "M13"):
            Hl[i, i] += 1.0
    nzl, wl, _ = null_count(Hl)
    wl_, Ul = np.linalg.eigh(Hl)
    nv = Ul[:, np.abs(wl_) < 1e-8 * np.max(np.abs(wl_))]
    Tb = tangent(GENS["boost_01"], M_DEG); Tb = np.array([np.sum(Tb * E) for E in B]); Tb /= np.linalg.norm(Tb)
    rec["MUT_lock"] = {"n_flat": nzl, "flat_overlap_with_boost01": float(np.abs(Tb @ nv).max()) if nv.size else 0.0}
    # V4-type potential at the degenerate spectrum, C_p matched: Hessian = 2 J^T J
    for pmax in (4, 6):
        N = M_DEG @ ETA
        Cp = [np.trace(np.linalg.matrix_power(N, p)) for p in range(1, pmax + 1)]
        V4 = lambda Mx, Cp=Cp, pmax=pmax: float(sum((np.trace(np.linalg.matrix_power(Mx @ ETA, p)) - Cp[p - 1]) ** 2 for p in range(1, pmax + 1)))
        Hv = fd_hessian(V4, M_DEG, eps=1e-3)
        # exact J^T J
        Jm = np.zeros((pmax, 10))
        for p in range(1, pmax + 1):
            for i, E in enumerate(B):
                Jm[p - 1, i] = p * np.trace(np.linalg.matrix_power(N, p - 1) @ E @ ETA)
        Hx = 2.0 * Jm.T @ Jm
        sv = np.linalg.svd(Jm, compute_uv=False)
        rank_J = int(np.sum(sv > 1e-9 * sv[0]))
        wv = np.linalg.eigvalsh(Hx); scv = float(np.max(np.abs(wv)))
        nzf, wvf, _ = null_count(Hv, rtol=1e-6)
        rec[f"V4type_pmax{pmax}"] = {"n_flat_exact": 10 - rank_J, "n_flat_fd_rtol1e-6": nzf, "rank_J": rank_J,
                                     "singular_values_J": [float(x) for x in sv],
                                     "hess_eigs_exact": [float(x) for x in wv], "hess_eigs_fd": [float(x) for x in wvf],
                                     "fd_vs_exact_rel": float(np.max(np.abs(Hv - Hx)) / scv)}
        log(f"C2 V4-type pmax {pmax}: flat exact {10 - rank_J} (fd rtol 1e-6: {nzf}), rank J {rank_J}, sv(J) {sv}")
    # non-invariant control: ||M - M0||_F^2 has 0 flat directions
    rec["MUT_noninvariant_frobenius"] = {"n_flat": null_count(2.0 * np.eye(10))[0]}
    main_ok = (rec["f_producer_deg"]["n_flat"] == 5 and rec["f_own_x2p3_deg"]["n_flat"] == 5
               and rec["MUT_nondegenerate_vac"]["n_flat"] == 6 and rec["MUT_lock"]["n_flat"] == 1
               and resid < 1e-8 and sorted(int(ch[k]) for k in ("boost_01", "boost_02", "boost_03", "rot_12", "rot_13")) == [0, 1, 1, 1, 1]
               and ch["p_split"] == 2 and ch["q_split"] == 2)
    v4_ok = rec["V4type_pmax4"]["n_flat_exact"] == 8
    rec["verdict"] = ("CONFIRMED" if main_ok and v4_ok else "QUALIFIED" if main_ok else "REFUTED")
    rec["verdict_note"] = (f"5/6/1 and charges confirmed; V4-type flat count at quadratic order = "
                           f"{rec['V4type_pmax4']['n_flat_exact']} (producer 8); Hessian eigs {rec['V4type_pmax4']['hess_eigs_exact']}")
    log(f"C2 {rec['verdict']}: {rec['verdict_note']}; charges {ch}")
    save("c2", rec)


# ============================================================ C3 / C11 sympy: R_G Euler-Lagrange, planar profiles, twist
def sym_field(coords):
    fs = {}
    M = sp.zeros(4, 4)
    for a in range(4):
        for b in range(a, 4):
            f = sp.Function(f"m{a}{b}")(*coords)
            fs[(a, b)] = f
            M[a, b] = f; M[b, a] = f
    return M, fs


def rg_sym(M, coords, dslots, Gm):
    """R_G with derivative index mu in dslots (mu also the raw row index), coords[k] the coordinate of dslots[k]."""
    A = {mu: M.diff(coords[k]) for k, mu in enumerate(dslots)}
    R = 0
    for mu in dslots:
        for nu in dslots:
            if mu == nu:
                continue
            for c in range(4):
                for d in range(4):
                    if Gm[c, d] == 0:
                        continue
                    R += Gm[c, d] * (A[mu][nu, c] * A[nu][mu, d] - A[mu][mu, c] * A[nu][nu, d])
    return sp.expand(R)


def euler_lagrange(Lag, fs, coords):
    out = {}
    for key, f in fs.items():
        dL_df = sp.diff(Lag, f)
        tot = dL_df
        for x in coords:
            tot -= sp.diff(sp.diff(Lag, sp.Derivative(f, x)), x)
        out[key] = sp.expand(tot)
    return out


def stage_c3():
    rec = {}
    XI = sp.diag(-1, 1, 1, 1)
    for label, coords, slots in (("1+1", sp.symbols("t x", real=True), (0, 1)),
                                 ("2+1", sp.symbols("t x y", real=True), (0, 1, 2))):
        M, fs = sym_field(coords)
        R = rg_sym(M, coords, slots, XI)
        EL = euler_lagrange(R, fs, coords)
        nz = {f"m{a}{b}": str(e) for (a, b), e in EL.items() if e != 0}
        rec[label] = {"n_terms_density": len(R.args), "EL_nonzero_components": nz, "EL_all_zero": len(nz) == 0}
        # MUTATIONS: (i) plus sign between the two terms; (ii) G = eta M eta
        A = {mu: M.diff(coords[k]) for k, mu in enumerate(slots)}
        Rp = 0
        for mu in slots:
            for nu in slots:
                if mu == nu:
                    continue
                for c in range(4):
                    Rp += XI[c, c] * (A[mu][nu, c] * A[nu][mu, c] + A[mu][mu, c] * A[nu][nu, c])
        ELp = euler_lagrange(sp.expand(Rp), fs, coords)
        rec[label]["MUT_plus_sign_EL_nonzero_count"] = sum(1 for e in ELp.values() if e != 0)
        if label == "1+1":
            Rm = rg_sym(M, coords, slots, XI * M * XI)
            ELm = euler_lagrange(Rm, fs, coords)
            rec[label]["MUT_G_etaMeta_EL_nonzero_count"] = sum(1 for e in ELm.values() if e != 0)
        log(f"C3 {label}: EL all zero {rec[label]['EL_all_zero']}; mutant plus-sign nonzero {rec[label]['MUT_plus_sign_EL_nonzero_count']}")
    rec["verdict"] = "CONFIRMED" if rec["1+1"]["EL_all_zero"] and rec["2+1"]["EL_all_zero"] else "REFUTED"
    save("c3", rec)


def stage_c11():
    rec = {}
    XI = sp.diag(-1, 1, 1, 1)
    g, d = sp.Integer(8), sp.Rational(3, 10)
    z, t = sp.symbols("z t", real=True)
    f = lambda x: (x + g) * (x - 1)
    # S3: planar profile M(z), generic 10 functions: R_G with the single derivative slot 3
    M, fs = sym_field((z,))
    R_eta = rg_sym(M, (z,), (3,), XI)
    R_eMe = rg_sym(M, (z,), (3,), XI * M * XI)
    rec["S3_planar_R_eta"] = str(R_eta); rec["S3_planar_R_etaMeta"] = str(R_eMe)
    # S3 mutation: a 2D profile M(y, z) is NOT identically zero
    y = sp.symbols("y", real=True)
    M2, fs2 = sym_field((y, z))
    R2 = rg_sym(M2, (y, z), (2, 3), XI)
    rec["S3_MUT_2D_profile_R_eta_nonzero"] = bool(R2 != 0)
    # S3: K_P on a planar (2,3)-split ramp M = diag(8, 1, d/2 + s(z), d/2 - s(z))
    s = sp.Function("s")(z)
    Ms = sp.diag(g, 1, d / 2 + s, d / 2 - s)
    Ns = Ms * XI
    Pm = (Ns + g * sp.eye(4)) * (Ns - sp.eye(4))
    Om = Pm * Ms.diff(z) * XI * Pm
    KP_split = sp.simplify((Om * Om).trace() / 2)
    rec["S3_KP_split_ramp"] = str(KP_split)
    # S4: uniform vacuum with a phase phi(t, z), R = exp(phi G1)
    phi = sp.Function("phi")(t, z)
    G1s = sp.zeros(4, 4); G1s[2, 3], G1s[3, 2] = -1, 1
    Rm = sp.eye(4) + sp.sin(phi) * G1s + (1 - sp.cos(phi)) * G1s * G1s
    d4 = sp.diag(g, 1, d, 0)
    Mp = Rm * d4 * Rm.T
    Np = Mp * XI
    Pp = (Np + g * sp.eye(4)) * (Np - sp.eye(4))
    out = {}
    for nm, x in (("time", t), ("space", z)):
        Omp = Pp * Mp.diff(x) * XI * Pp
        dens = sp.simplify((Omp * Omp).trace() / 2)
        coef = sp.simplify(dens / sp.Derivative(phi, x) ** 2)
        out[nm] = str(coef)
    target = sp.simplify((f(d) * f(0)) ** 2 * d ** 2)
    rec["S4_KP_coefficient"] = out
    rec["S4_target_[f(d)f(0)]^2 d^2"] = str(target)
    rec["S4_match"] = all(sp.simplify(sp.sympify(v) - target) == 0 for v in out.values())
    # S4 mutation: literal roots (+g, 1) give a different coefficient
    Pl = (Np - g * sp.eye(4)) * (Np - sp.eye(4))
    Oml = Pl * Mp.diff(z) * XI * Pl
    rec["S4_MUT_literal_roots_coef"] = str(sp.simplify((Oml * Oml).trace() / 2 / sp.Derivative(phi, z) ** 2))
    # S5: the (1,2) twist psi(z), R = exp(psi G3): K_P, R_G, K_lambda static = 0; C6 quartic = 4 (d-1)^4 psi'^4
    psi = sp.Function("psi")(z)
    G3s = sp.zeros(4, 4); G3s[1, 2], G3s[2, 1] = -1, 1
    Rt = sp.eye(4) + sp.sin(psi) * G3s + (1 - sp.cos(psi)) * G3s * G3s
    Mt = Rt * d4 * Rt.T
    Nt = Mt * XI
    Pt = (Nt + g * sp.eye(4)) * (Nt - sp.eye(4))
    dMt = Mt.diff(z)
    Omt = Pt * dMt * XI * Pt
    KP_t = sp.simplify((Omt * Omt).trace() / 2)
    RG_t = rg_sym(Mt, (z,), (3,), XI)
    lam_t = [sp.simplify(e) for e in Nt.eigenvals().keys()]
    Klam_t = sp.simplify(sum(sp.diff(e, z) ** 2 for e in lam_t) / 2)
    quart = sp.simplify(((dMt * XI * dMt * XI).trace()) ** 2)
    rec["S5"] = {"KP_static": str(KP_t), "RG_eta": str(RG_t), "Klambda_static": str(Klam_t),
                 "eigenvalues": [str(e) for e in lam_t], "C6_quartic": str(quart),
                 "C6_quartic_target": str(sp.expand(4 * (d - 1) ** 4 * sp.Derivative(psi, z) ** 4)),
                 "C6_match": bool(sp.simplify(quart - 4 * (d - 1) ** 4 * sp.Derivative(psi, z) ** 4) == 0)}
    # S5 mutation: the (2,3) twist (G1) has nonzero K_P static cost
    Rm1 = sp.eye(4) + sp.sin(psi) * G1s + (1 - sp.cos(psi)) * G1s * G1s
    Mm = Rm1 * d4 * Rm1.T; Nm = Mm * XI
    Pmm = (Nm + g * sp.eye(4)) * (Nm - sp.eye(4))
    Omm = Pmm * Mm.diff(z) * XI * Pmm
    rec["S5_MUT_23_twist_KP_static"] = str(sp.simplify((Omm * Omm).trace() / 2))
    ok = (R_eta == 0 and R_eMe == 0 and KP_split != 0 and rec["S4_match"] and KP_t == 0 and RG_t == 0
          and Klam_t == 0 and rec["S5"]["C6_match"])
    rec["verdict"] = "CONFIRMED" if ok else "REFUTED"
    log(f"C11 {rec['verdict']}: S3 planar R=({R_eta},{R_eMe}) KP_split={KP_split}; S4 {out} target {target}; "
        f"S5 KP={KP_t} RG={RG_t} Klam={Klam_t} quartic={quart}")
    save("c11", rec)


# ============================================================ C4 / C5(2D): periodic-box integrals of R_G
SLOT_PAIRS = {"(t,x)": (0, 1), "(t,y)": (0, 2), "(t,z)": (0, 3), "(x,y)": (1, 2), "(x,z)": (1, 3), "(y,z)": (2, 3)}


def box_grid(N):
    x = 2 * np.pi * np.arange(N) / N
    return np.meshgrid(x, x, indexing="ij")


def spec_d(F, axis):
    N = F.shape[axis]
    k = np.fft.fftfreq(N, d=1.0 / N)
    Fk = np.fft.fft(F, axis=axis)
    shape = [1] * F.ndim; shape[axis] = N
    return np.real(np.fft.ifft(1j * k.reshape(shape) * Fk, axis=axis))


def fd4_d(F, axis, N):
    h = 2 * np.pi / N
    r = lambda s: np.roll(F, -s, axis=axis)
    return (-r(2) + 8 * r(1) - 8 * r(-1) + r(-2)) / (12 * h)


def textures(N, seed, kind):
    X, Y = box_grid(N)
    if kind == "det":
        return [0.8 * np.sin(X) + 0.3 * np.cos(2 * Y), 0.6 * np.cos(X - Y) + 0.2 * np.sin(3 * X), 0.7 * np.sin(Y) * np.cos(X)]
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(3):
        f = np.zeros_like(X)
        for k1 in range(-2, 3):
            for k2 in range(-2, 3):
                if k1 == 0 and k2 == 0:
                    continue
                f += rng.normal() * 0.25 / (1 + k1 * k1 + k2 * k2) * np.cos(k1 * X + k2 * Y + rng.uniform(0, 2 * np.pi))
        out.append(f)
    return out


def L_of(orbit, th):
    th1, th2, th3 = th
    if orbit == "rot":          # order (12, 23, 13): different from a (12, 13, 23) ordering
        return rot_of(GENS["rot_12"], th1) @ rot_of(GENS["rot_23"], th2) @ rot_of(GENS["rot_13"], th3)
    if orbit == "boost":        # order (03, 01, 02)
        return boost_of(GENS["boost_03"], th1) @ boost_of(GENS["boost_01"], th2) @ boost_of(GENS["boost_02"], th3)
    if orbit == "mixed_K03":    # a single boost K_03 with rot_12, rot_13
        return boost_of(GENS["boost_03"], th1) @ rot_of(GENS["rot_12"], th2) @ rot_of(GENS["rot_13"], th3)
    if orbit == "boost_single_K01":
        return boost_of(GENS["boost_01"], th1)
    raise ValueError(orbit)


def rg_box_integral(M, slots, Gname, deriv="spec", N=None, sign=-1.0):
    N = M.shape[0] if N is None else N
    d = (lambda F, ax: spec_d(F, ax)) if deriv == "spec" else (lambda F, ax: fd4_d(F, ax, N))
    A = {slots[0]: d(M, 0), slots[1]: d(M, 1)}
    Gm = G_of(Gname, M)
    if sign > 0:   # wrong-sign mutant: plus between the two terms
        dens = 0.0
        for mu in A:
            for nu in A:
                if mu == nu:
                    continue
                t1 = np.einsum("...cd,...c,...d->...", Gm, A[mu][..., nu, :], A[nu][..., mu, :])
                t2 = np.einsum("...cd,...c,...d->...", Gm, A[mu][..., mu, :], A[nu][..., nu, :])
                dens = dens + (t1 + t2)
    else:
        dens = rg_density(A, Gm)
    dA = (2 * np.pi / N) ** 2
    return float(np.sum(dens) * dA), float(np.max(np.abs(dens)))


def stage_c4():
    rec = {"M0": [8.0, 1.0, 0.3, 0.1], "box": "2 pi periodic, spectral derivatives"}
    M0 = np.diag([8.0, 1.0, 0.3, 0.1])
    Gnames = ("eta", "eta_M_eta", "Minv", "h_cov")
    Ns = (8, 16, 32, 64)
    table = {}
    for orbit in ("rot", "boost", "mixed_K03"):
        for tex in ("det", "rnd1", "rnd2"):
            for sl, slots in SLOT_PAIRS.items():
                key = f"{orbit}|{tex}|{sl}"
                row = {}
                for Gn in Gnames:
                    vals = []
                    for N in Ns:
                        th = textures(N, {"det": 0, "rnd1": 1, "rnd2": 2}[tex], "det" if tex == "det" else "rnd")
                        L = L_of(orbit, th)
                        M = conj(L, M0)
                        vals.append(rg_box_integral(M, slots, Gn)[0])
                    row[Gn] = vals
                table[key] = row
    rec["integrals_by_N"] = table
    # summary: rotation orbits zero? boosts nonzero and converged?
    summ = {}
    for orbit in ("rot", "boost", "mixed_K03"):
        for Gn in Gnames:
            mx, mn_boost, conv = 0.0, np.inf, 0.0
            per_slot = {}
            for key, row in table.items():
                o, tex, sl = key.split("|")
                if o != orbit:
                    continue
                v = row[Gn]
                per_slot.setdefault(sl, []).append(v[-1])
                mx = max(mx, abs(v[-1]))
                conv = max(conv, abs(v[-1] - v[-2]) / max(abs(v[-1]), 1e-12))
            summ[f"{orbit}|{Gn}"] = {"max_abs_N64": mx, "max_rel_change_32_to_64": conv,
                                     "N64_by_slot_det_rnd1_rnd2": per_slot}
    rec["summary"] = summ
    # pointwise: is the rotation-orbit density zero pointwise or only integrated?
    th = textures(64, 0, "det")
    M = conj(L_of("rot", th), M0)
    rec["rot_det_pointwise_max_density_xy"] = {Gn: rg_box_integral(M, (1, 2), Gn)[1] for Gn in Gnames}
    # finite-difference variant on a fine grid (boost det, slots (x,y) and (x,z))
    fdv = {}
    for N in (128, 256):
        th = textures(N, 0, "det")
        M = conj(L_of("boost", th), M0)
        fdv[N] = {sl: {Gn: rg_box_integral(M, SLOT_PAIRS[sl], Gn, deriv="fd4", N=N)[0] for Gn in Gnames}
                  for sl in ("(x,y)", "(x,z)")}
    rec["boost_det_fd4"] = fdv
    rec["boost_det_spec_N64"] = {sl: {Gn: table[f"boost|det|{sl}"][Gn][-1] for Gn in Gnames} for sl in ("(x,y)", "(x,z)")}
    # MUTATIONS: wrong sign -> nonzero on rotation orbits for G = eta; non-covariant G = identity on rotation orbit
    th = textures(64, 0, "det")
    Mr = conj(L_of("rot", th), M0)
    rec["MUT_plus_sign_rot_det_xy_eta"] = rg_box_integral(Mr, (1, 2), "eta", sign=+1.0)[0]
    rec["MUT_identity_G_rot_det"] = {sl: rg_box_integral(Mr, SLOT_PAIRS[sl], "identity")[0] for sl in SLOT_PAIRS}
    Mb = conj(L_of("boost", th), M0)
    rec["MUT_identity_G_boost_det"] = {sl: rg_box_integral(Mb, SLOT_PAIRS[sl], "identity")[0] for sl in SLOT_PAIRS}
    # amplitude scan on boosts (x,y): does the integral scale like amplitude^2 or higher?
    amp = {}
    for a in (0.25, 0.5, 1.0, 2.0):
        th = [a * f for f in textures(64, 0, "det")]
        Mb = conj(L_of("boost", th), M0)
        amp[a] = {Gn: rg_box_integral(Mb, (1, 2), Gn)[0] for Gn in Gnames}
    rec["boost_amplitude_scan_xy"] = amp
    # ordering study: which spatial slot pair vanishes depends on the OUTERMOST boost factor
    import itertools
    Ks = {"01": GENS["boost_01"], "02": GENS["boost_02"], "03": GENS["boost_03"]}
    th = textures(32, 1, "rnd")
    order_rows = {}
    for order in itertools.permutations(["01", "02", "03"]):
        L = boost_of(Ks[order[0]], th[0]) @ boost_of(Ks[order[1]], th[1]) @ boost_of(Ks[order[2]], th[2])
        Mo = conj(L, M0)
        order_rows["-".join(order)] = {sl: rg_box_integral(Mo, s2, "eta_M_eta")[0] for sl, s2 in SLOT_PAIRS.items()}
    rec["ordering_study_rnd1_N32_etaMeta"] = order_rows
    zero_pairs = {o: [sl for sl, v in row.items() if abs(v) < 1e-9] for o, row in order_rows.items()}
    rec["ordering_study_zero_slots"] = zero_pairs
    single = {}
    for k, K in Ks.items():
        Mo = conj(boost_of(K, th[0]), M0)
        single[k] = {Gn: max(abs(rg_box_integral(Mo, s2, Gn)[0]) for s2 in SLOT_PAIRS.values()) for Gn in Gnames}
    rec["single_boost_texture_max_abs_all_slots"] = single
    # the producer's ordering (01, 02, 03) on the det texture at N = 64
    L = boost_of(Ks["01"], textures(64, 0, "det")[0]) @ boost_of(Ks["02"], textures(64, 0, "det")[1]) @ boost_of(Ks["03"], textures(64, 0, "det")[2])
    Mo = conj(L, M0)
    rec["producer_ordering_det_N64"] = {Gn: {sl: rg_box_integral(Mo, s2, Gn)[0] for sl, s2 in SLOT_PAIRS.items()} for Gn in Gnames}
    rot_zero = all(summ[f"rot|{Gn}"]["max_abs_N64"] < 1e-9 for Gn in Gnames)
    eta_zero = all(summ[f"{o}|eta"]["max_abs_N64"] < 1e-9 for o in ("rot", "boost", "mixed_K03"))
    # nonzero + converged: for each G, every slot with |I(64)| > 1e-3 must have |I(64) - I(32)| / |I(64)| < 1e-6
    conv = {}
    for Gn in ("eta_M_eta", "Minv", "h_cov"):
        worst, nnz = 0.0, 0
        for key, row in table.items():
            if not key.startswith("boost|"):
                continue
            v = row[Gn]
            if abs(v[-1]) > 1e-3:
                nnz += 1
                worst = max(worst, abs(v[-1] - v[-2]) / abs(v[-1]))
        conv[Gn] = {"n_nonzero_slots_of_18": nnz, "worst_rel_change_32_to_64": worst}
    rec["boost_convergence"] = conv
    boost_nonzero = all(conv[Gn]["n_nonzero_slots_of_18"] >= 6 and conv[Gn]["worst_rel_change_32_to_64"] < 1e-5 for Gn in conv)
    mixed = {sl: max(abs(table[f"mixed_K03|{tex}|{sl}"][Gn][-1]) for tex in ("det", "rnd1", "rnd2") for Gn in ("eta_M_eta", "Minv", "h_cov"))
             for sl in SLOT_PAIRS}
    rec["mixed_K03_max_abs_by_slot"] = mixed
    mixed_ok = mixed["(x,y)"] < 1e-9 and all(mixed[sl] > 1e-4 for sl in ("(x,z)", "(y,z)", "(t,x)", "(t,y)"))
    rec["checks"] = {"rot_zero": rot_zero, "eta_zero_everywhere": eta_zero, "boost_nonzero_converged": boost_nonzero, "mixed_pattern": mixed_ok,
                     "mutant_plus_sign_nonzero": abs(rec["MUT_plus_sign_rot_det_xy_eta"]) > 1e-6}
    rec["verdict"] = "QUALIFIED" if all(rec["checks"].values()) else "REFUTED"
    rec["verdict_note"] = ("qualitative content confirmed with own textures/orderings; the producer's slot values are texture-specific "
                           "and not reproducible from the brief; found: the spatial slot pair orthogonal to the OUTERMOST boost factor "
                           "integrates to zero, and a single-plane boost texture integrates to zero on every slot for every G")
    log(f"C4 {rec['verdict']} checks {rec['checks']}; boost det N64 {rec['boost_det_spec_N64']}; fd4 {fdv[256]}")
    log(f"C4 mixed K03 by slot {mixed}; rot pointwise {rec['rot_det_pointwise_max_density_xy']}")
    log(f"C4 ordering zero slots {zero_pairs}; single boost {single}; convergence {conv}")
    save("c4", rec)


# ============================================================ C5: no omega^2 content of R_G on the hedgehog
def rg_lattice(M, h, a0, omega, Gname, st="sym"):
    Gm = G_of(Gname, M)
    tot = 0.0
    for Al, wt in jets(M, h, 3, st):
        A = {0: omega * a0, 1: Al[0], 2: Al[1], 3: Al[2]}
        tot += wt * float(np.sum(rg_density(A, Gm)))
    return h ** 3 * tot


def stage_c5():
    rec = {}
    M, n, L, h = load_field("n32_L48_it12000")
    a0 = a0_local_own(M)
    rec["a0_local_vs_r13w_maxdiff"] = float(np.max(np.abs(a0 - R13W.a0_local(M))))
    for Gn in ("eta", "eta_M_eta", "h_cov"):
        Lp, Lm, L0 = (rg_lattice(M, h, a0, om, Gn) for om in (1.0, -1.0, 0.0))
        C, Bc = 0.5 * (Lp + Lm) - L0, 0.5 * (Lp - Lm)
        L0c = rg_lattice(M, h, a0, 0.0, Gn, st="ctr")
        rec[Gn] = {"A_static_sym": L0, "A_static_ctr": L0c, "B_linear": Bc, "C_omega2": C}
        log(f"C5 {Gn}: A {L0:.4f} (ctr {L0c:.4f}) B {Bc:.4e} C {C:.3e}")
    # positive control: the same three-point read on the certified clock inertia has C = kin != 0
    kin_own = 0.0
    for Al, wt in jets(M, h):
        kin_own += wt * float(np.sum(kin_density(Al, a0)))
    kin_own *= h ** 3
    Lk = lambda om: sum(wt * float(np.sum(i1_density([om * a0] + Al))) for Al, wt in jets(M, h)) * h ** 3
    Ck = 0.5 * (Lk(1.0) + Lk(-1.0)) - Lk(0.0)
    rec["control_I1_three_point_C"] = Ck
    rec["control_kin_own"] = kin_own
    rec["control_kin_certified"] = float(INS4.kin_of(M, a0, R13W.cfg_of(n, L)))
    # MUTATION: a generator with an explicit omega^2 carrier: A_0 = omega a0 + omega^2 b0 (must give C != 0)
    rng = np.random.default_rng(7)
    b0 = rng.normal(size=M.shape); b0 = 0.5 * (b0 + np.swapaxes(b0, -1, -2))
    Lq = lambda om: sum(wt * float(np.sum(rg_density({0: om * a0 + om * om * b0, 1: Al[0], 2: Al[1], 3: Al[2]}, G_of("eta", M))))
                        for Al, wt in jets(M, h)) * h ** 3
    rec["MUT_quadratic_carrier_C"] = 0.5 * (Lq(1.0) + Lq(-1.0)) - Lq(0.0)
    # block structure of the relaxed hedgehog (explains R_eta == R_hcov and B == 0)
    rec["hedgehog_max_abs_M0i"] = float(np.max(np.abs(M[..., 0, 1:])))
    rec["hedgehog_M00_range"] = [float(M[..., 0, 0].min()), float(M[..., 0, 0].max())]
    # the linear-in-omega coefficient B on generator tangents of the hedgehog and of a boost-dressed state
    Bt = {}
    cfg = R13W.cfg_of(n, L)
    Md = R13W.B8.dressed(cfg, 0.12)
    for lab, Mx in (("hedgehog", M), ("dressed_m0.12", Md)):
        for nm in ("boost_01", "rot_12", "rot_23", "local_clock"):
            ax = a0_local_own(Mx) if nm == "local_clock" else tangent(GENS[nm], Mx)
            Lp, Lm, L0 = (rg_lattice(Mx, h, ax, om, "eta_M_eta") for om in (1.0, -1.0, 0.0))
            Bt[f"{lab}|{nm}"] = {"A": L0, "B": 0.5 * (Lp - Lm), "C": 0.5 * (Lp + Lm) - L0}
    rec["linear_term_B_etaMeta"] = Bt
    # 2D sub-claims: constant frame + eigenvalue gradients; eigenvalue gradients + frame texture
    N = 64
    X, Y = box_grid(N)
    lam = np.stack([8.0 + 0.4 * np.sin(X) * np.cos(Y), 1.0 + 0.1 * np.cos(2 * X), 0.3 + 0.05 * np.sin(X + Y), 0.1 + 0.03 * np.cos(Y)], -1)
    Md = np.zeros((N, N, 4, 4))
    for a in range(4):
        Md[..., a, a] = lam[..., a]
    L0 = boost_of(GENS["boost_01"], 0.4) @ rot_of(GENS["rot_12"], 0.7)
    Mc = conj(np.broadcast_to(L0, Md.shape), Md)
    rec["constframe_eigengrad"] = {Gn: {sl: rg_box_integral(Mc, SLOT_PAIRS[sl], Gn) for sl in ("(x,y)", "(t,z)")}
                                   for Gn in ("eta", "eta_M_eta", "Minv", "h_cov")}
    tex = {}
    for N in (32, 64):
        X, Y = box_grid(N)
        lam = np.stack([8.0 + 0.4 * np.sin(X) * np.cos(Y), 1.0 + 0.1 * np.cos(2 * X), 0.3 + 0.05 * np.sin(X + Y), 0.1 + 0.03 * np.cos(Y)], -1)
        Md = np.zeros((N, N, 4, 4))
        for a in range(4):
            Md[..., a, a] = lam[..., a]
        Lr = L_of("rot", textures(N, 0, "det"))
        Mt = conj(Lr, Md)
        tex[N] = {Gn: {sl: rg_box_integral(Mt, SLOT_PAIRS[sl], Gn)[0] for sl in ("(x,y)", "(x,z)")} for Gn in ("eta", "eta_M_eta", "Minv", "h_cov")}
    rec["eigengrad_plus_rot_texture"] = tex
    cf_zero = all(abs(v[0]) < 1e-9 for Gd in rec["constframe_eigengrad"].values() for v in Gd.values())
    tex_nonzero = all(abs(tex[64][Gn]["(x,y)"]) > 1e-4 and abs(tex[64][Gn]["(x,y)"] - tex[32][Gn]["(x,y)"]) < 1e-6 for Gn in ("eta_M_eta", "Minv"))
    rec["checks"] = {"C_zero_all_G": all(abs(rec[Gn]["C_omega2"]) < 1e-8 for Gn in ("eta", "eta_M_eta", "h_cov")),
                     "control_detects_omega2": abs(Ck - kin_own) < 1e-6 * kin_own,
                     "mutant_nonzero": abs(rec["MUT_quadratic_carrier_C"]) > 1e-6,
                     "constframe_zero": cf_zero, "texture_nonzero_converged": tex_nonzero}
    rec["verdict"] = "CONFIRMED" if all(rec["checks"].values()) else "QUALIFIED"
    log(f"C5 {rec['verdict']} {rec['checks']} control C {Ck:.3f} kin {kin_own:.3f}; texture {tex[64]}")
    save("c5", rec)


# ============================================================ C6 / C13: K_lambda, the triple, exponents
def field_energies(key):
    M, n, L, h = load_field(key)
    cfg = R13W.cfg_of(n, L)
    eu = ev = kp = kpc = klam = kin = 0.0
    a0 = a0_local_own(M)
    P = P_of(M)
    for Al, wt in jets(M, h):
        eu += wt * float(np.sum(i1_density(Al)))
        kp += wt * float(np.sum(kp_plain_density(Al, P)))
        kin += wt * float(np.sum(kin_density(Al, a0)))
    kpc = float(np.sum(kp_plain_density(jets(M, h, st="ctr")[0][0], P)))
    kl_d, imag = klam_static_density(M, h)
    klam = float(np.sum(kl_d))
    e_u_cert, v4_cert = INS4.e_parts(M, cfg)
    kin_cert = float(INS4.kin_of(M, R13W.a0_local(M), cfg))
    return {"n": n, "L": L, "h": h, "E_u_own": h ** 3 * eu, "E_u_cert": float(e_u_cert), "V4_cert": float(v4_cert),
            "KP_static_sym": h ** 3 * kp, "KP_static_ctr": h ** 3 * kpc, "kin_I1_local_own": h ** 3 * kin, "kin_I1_local_cert": kin_cert,
            "Klambda_static": h ** 3 * klam, "eig_max_imag": imag}


def exponents(vals, Ls):
    lv, lL = np.log(np.array(vals)), np.log(np.array(Ls))
    return {"pair_48_60": float((lv[1] - lv[0]) / (lL[1] - lL[0])), "pair_60_72": float((lv[2] - lv[1]) / (lL[2] - lL[1])),
            "lsq_3pt": float(np.polyfit(lL, lv, 1)[0])}


def stage_triple():
    rec = {}
    for key in ("n32_L48_it3000", "n40_L60_it3000", "n48_L72_it3000"):
        rec[key] = field_energies(key)
        log(f"triple {key}: {rec[key]}")
    Ls = [48.0, 60.0, 72.0]
    keys = ("n32_L48_it3000", "n40_L60_it3000", "n48_L72_it3000")
    rec["exponents"] = {q: exponents([rec[k][q] for k in keys], Ls)
                        for q in ("E_u_own", "KP_static_sym", "KP_static_ctr", "kin_I1_local_own", "Klambda_static")}
    # hypothesis: the producer's n32 member is the 12000-it field (seed_hedgehog(32, 48) returns it)
    rec["n32_L48_it12000"] = field_energies("n32_L48_it12000")
    keys2 = ("n32_L48_it12000", "n40_L60_it3000", "n48_L72_it3000")
    rec["exponents_producer_triple_it12000_n32"] = {q: exponents([rec[k][q] for k in keys2], Ls)
                                                    for q in ("E_u_own", "KP_static_sym", "kin_I1_local_own", "Klambda_static")}
    log(f"triple exps (listed 3000-it n32): {rec['exponents']}")
    log(f"triple exps (12000-it n32): {rec['exponents_producer_triple_it12000_n32']}")
    save("triple", rec)
    return rec


def stage_c6():
    rec = {}
    trip = json.load(open(OUT)).get("triple") or stage_triple()
    M, n, L, h = load_field("n32_L48_it12000")
    # omega^2 coefficient of K_lambda on the generator channels
    chans = {"clock_local": a0_local_own(M), "plane_1d": a0_lowest_own(M),
             "rot_z": tangent(GENS["rot_12"], M), "rot_x": tangent(GENS["rot_23"], M),
             "boost_z": tangent(GENS["boost_03"], M), "boost_x": tangent(GENS["boost_01"], M)}
    om2 = {}
    for nm, a0 in chans.items():
        om2[nm] = float(h ** 3 * np.sum(klam_omega2_density(M, a0)))
        # cross-check by finite differences of the sorted eigenvalues along the flow
        eps = 1e-4
        lp, _ = lam_sorted(M + eps * a0); lm, _ = lam_sorted(M - eps * a0)
        om2[nm + "_fd"] = float(h ** 3 * 0.5 * np.sum(((lp - lm) / (2 * eps)) ** 2))
    rec["omega2_coefficients"] = om2
    # MUTATION: a non-tangent direction (the z-jet itself) carries eigenvalue motion
    b0 = dctr(M, 2, h)
    rec["MUT_nontangent_zjet_omega2"] = float(h ** 3 * np.sum(klam_omega2_density(M, b0)))
    # static = 0 on cell-wise constant spectrum: own 3D orbit field with random textures
    n3, h3 = 16, 1.5
    X, Y, Z = coords3(n3, h3)
    th = [0.5 * np.sin(2 * np.pi * X / 24) * np.cos(2 * np.pi * Y / 24), 0.4 * np.cos(2 * np.pi * (Y + Z) / 24), 0.3 * np.sin(2 * np.pi * Z / 24)]
    Lr = L_of("mixed_K03", th)
    Morb = conj(Lr, M_VAC)
    kl_d, imag = klam_static_density(Morb, h3)
    rec["static_on_orbit_field"] = {"Klambda": float(h3 ** 3 * np.sum(kl_d)), "eig_max_imag": imag,
                                    "E_u_of_that_field": float(h3 ** 3 * sum(wt * np.sum(i1_density(Al)) for Al, wt in jets(Morb, h3)))}
    # MUTATION: eigenvalue-textured field has nonzero static K_lambda
    Mev = Morb.copy()
    lamt = np.stack([8 + 0.2 * np.sin(2 * np.pi * X / 24), np.ones_like(X), 0.3 + 0.05 * np.cos(2 * np.pi * Y / 24), 0.1 * np.sin(2 * np.pi * Z / 24) ** 2], -1)
    Md = np.zeros_like(Morb)
    for a in range(4):
        Md[..., a, a] = lamt[..., a]
    Mev = conj(Lr, Md)
    rec["MUT_eigen_textured_static"] = float(h3 ** 3 * np.sum(klam_static_density(Mev, h3)[0]))
    rec["triple_Klambda"] = {k: trip[k]["Klambda_static"] for k in ("n32_L48_it3000", "n40_L60_it3000", "n48_L72_it3000")}
    rec["triple_exponents_Klambda"] = trip["exponents"]["Klambda_static"]
    ok_om = all(abs(v) < 1e-20 for k, v in om2.items() if not k.endswith("_fd"))
    ok_static = rec["static_on_orbit_field"]["Klambda"] < 1e-20
    vals = list(rec["triple_Klambda"].values())
    rec["Klambda_n32_it12000"] = trip["n32_L48_it12000"]["Klambda_static"]
    vals_p = [rec["Klambda_n32_it12000"], vals[1], vals[2]]
    ok_vals = all(abs(v - p) / p < 0.02 for v, p in zip(vals_p, (9.18, 6.38, 6.71)))
    rec["producer_values_reproduced_with_it12000_n32"] = ok_vals
    rec["triple_exponents_Klambda_producer_triple"] = trip["exponents_producer_triple_it12000_n32"]["Klambda_static"]
    rec["checks"] = {"omega2_zero_all_channels": ok_om, "static_zero_on_orbit": ok_static, "triple_values_match": ok_vals,
                     "mutants_nonzero": rec["MUT_nontangent_zjet_omega2"] > 1e-6 and rec["MUT_eigen_textured_static"] > 1e-6}
    rec["verdict"] = "QUALIFIED" if all(rec["checks"].values()) else "REFUTED"
    rec["verdict_note"] = ("omega^2 = 0 and static = 0 confirmed; the 9.18 / 6.38 / 6.71 triple reproduces only with the "
                           "12000-it n32 field (the listed 3000-it n32 gives %.3f), so the triple is maturity-mismatched; "
                           "with the true 3000-it triple K_lambda = %.3f / %.3f / %.3f, exponents %s" % (vals[0], vals[0], vals[1], vals[2], rec["triple_exponents_Klambda"]))
    log(f"C6 {rec['verdict']} {rec['checks']} om2 {om2} triple {vals} exps {rec['triple_exponents_Klambda']}")
    save("c6", rec)


def stage_c13():
    trip = json.load(open(OUT)).get("triple") or stage_triple()
    keys = ("n32_L48_it3000", "n40_L60_it3000", "n48_L72_it3000")
    rec = {"values": {q: [trip[k][q] for k in keys] for q in ("E_u_own", "E_u_cert", "KP_static_sym", "KP_static_ctr", "kin_I1_local_own", "kin_I1_local_cert", "Klambda_static")},
           "exponents": trip["exponents"],
           "producer": {"KP": [95381, 151565, 193653], "KP_exp": 1.34, "E_u_exp": 0.22, "kin_exp": 1.23, "Klam_exp": 0.26}}
    kp = rec["values"]["KP_static_sym"]
    rec["KP_match_1pct_listed_triple"] = all(abs(a - b) / b < 0.01 for a, b in zip(kp, rec["producer"]["KP"]))
    kp2 = [trip["n32_L48_it12000"]["KP_static_sym"], kp[1], kp[2]]
    rec["KP_match_1pct_it12000_n32"] = all(abs(a - b) / b < 0.01 for a, b in zip(kp2, rec["producer"]["KP"]))
    rec["values_n32_it12000"] = {q: trip["n32_L48_it12000"][q] for q in ("E_u_own", "KP_static_sym", "kin_I1_local_own", "Klambda_static")}
    rec["exponents_producer_triple_it12000_n32"] = trip["exponents_producer_triple_it12000_n32"]
    rec["verdict"] = "QUALIFIED" if rec["KP_match_1pct_it12000_n32"] else "REFUTED"
    log(f"C13 {rec['verdict']} KP {kp} exps {trip['exponents']}")
    save("c13", rec)


# ============================================================ C7 / C12: K_P omega^2 at the vacuum
def stage_c7():
    rec = {}
    f = lambda x: (x + G_) * (x - 1.0)
    target = (f(DELTA) * f(0.0)) ** 2 * DELTA ** 2
    rec["target_[f(d)f(0)]^2 d^2"] = target
    for roots, lab in (((-G_, 1.0), "intended"), ((G_, 1.0), "literal")):
        P = P_of(M_VAC, *roots)
        row = {}
        for nm, X in GENS.items():
            a0 = tangent(X, M_VAC)
            row[nm] = float(kp_plain_density([a0], P))
        rec[lab] = row
    # local clock generator at the vacuum equals the rot_23 tangent
    a0l = a0_local_own(M_VAC[None])[0]
    rec["local_clock_equals_rot23_tangent"] = float(np.max(np.abs(a0l - tangent(GENS["rot_23"], M_VAC))))
    # MUTATION: degenerate vacuum -> the (2,3) clock tangent vanishes, coefficient 0
    rec["MUT_degenerate_vacuum_rot23"] = float(kp_plain_density([tangent(GENS["rot_23"], M_DEG)], P_of(M_DEG)))
    # MUTATION: wrong internal metric (identity instead of eta in Om = P a0 eta P) on the boost_03 literal case
    P = P_of(M_VAC, G_, 1.0)
    a0 = tangent(GENS["boost_03"], M_VAC)
    Om_id = P @ a0 @ P
    rec["MUT_identity_metric_literal_boost03"] = float(0.5 * np.trace(Om_id @ Om_id))
    ok = (abs(rec["intended"]["rot_23"] - target) < 1e-9 * target
          and all(abs(rec["intended"][k]) < 1e-12 for k in ("boost_01", "boost_02", "boost_03", "rot_12", "rot_13"))
          and abs(abs(rec["literal"]["boost_03"]) - 8.5e7) / 8.5e7 < 0.01)
    rec["verdict"] = "CONFIRMED" if ok else "REFUTED"
    rec["note"] = f"literal boost_03 coefficient is {rec['literal']['boost_03']:.4e} (sign matters: NEGATIVE)"
    log(f"C7 {rec['verdict']} intended {rec['intended']} literal {rec['literal']} target {target:.4f}")
    save("c7", rec)


def stage_c12():
    h = 1.5
    om2_kp = (((DELTA + G_) * (DELTA - 1.0) * (0.0 + G_) * (0.0 - 1.0)) ** 2 * DELTA ** 2) * h ** 3
    # numeric on a small uniform vacuum grid with the local generator
    n = 6
    M = np.broadcast_to(M_VAC, (n, n, n, 4, 4)).copy()
    a0 = a0_local_own(M)
    P = P_of(M)
    kp_cell = float(h ** 3 * kp_plain_density([a0], P)[n // 2, n // 2, n // 2])
    kin_cell = float(h ** 3 * sum(wt * kin_density(Al, a0)[n // 2, n // 2, n // 2] for Al, wt in jets(M, h)))
    # V4 cost per cell of the degenerate state diag(8, 1, 0.15, 0.15) against the certified targets
    Cp = [(-G_) ** p + 1.0 + DELTA ** p for p in range(1, 5)]
    Nd = np.diag([8.0, 1.0, 0.15, 0.15]) @ ETA
    tp = [np.trace(np.linalg.matrix_power(Nd, p)) for p in range(1, 5)]
    v4_cell = h ** 3 * W1 * sum((a - b) ** 2 for a, b in zip(tp, Cp))
    cfg = R13W.cfg_of(6, 9.0)
    v4_cert = float(INS4.e_parts(np.broadcast_to(np.diag([8.0, 1.0, 0.15, 0.15]), (6, 6, 6, 4, 4)).copy(), cfg)[1]) / 6 ** 3
    rec = {"KP_omega2_per_cell_analytic": om2_kp, "KP_omega2_per_cell_numeric": kp_cell, "I1_omega2_per_cell": kin_cell,
           "V4_degenerate_per_cell": v4_cell, "V4_degenerate_per_cell_cert": v4_cert, "crossover_omega2": v4_cell / om2_kp,
           "MUT_degenerate_targets_V4": h ** 3 * W1 * sum((a - b) ** 2 for a, b in zip(tp, [(-G_) ** p + 1.0 + 2 * 0.15 ** p for p in range(1, 5)]))}
    ok = abs(om2_kp - 656.22) < 0.01 and abs(kp_cell - om2_kp) < 1e-9 and kin_cell == 0.0 and abs(v4_cell - 6.073e-6) / 6.073e-6 < 1e-3 \
        and abs(rec["crossover_omega2"] - 9.25e-9) / 9.25e-9 < 1e-2
    rec["verdict"] = "CONFIRMED" if ok else "REFUTED"
    log(f"C12 {rec['verdict']} {rec}")
    save("c12", rec)


# ============================================================ C8: indefiniteness of plain K_P; H-adjoint
def stage_c8():
    rng = np.random.default_rng(20260904)
    nprobe = 2000
    Ms, As = [], []
    for _ in range(nprobe):
        R = rng.uniform(-1, 1, (4, 4)); R = 0.5 * (R + R.T)
        Ms.append(M_VAC + 0.3 * R)
        Ai = []
        for i in range(3):
            J = rng.normal(size=(4, 4)); Ai.append(0.5 * (J + J.T))
        As.append(Ai)
    Ms = np.array(Ms); As = np.array(As)          # (n,4,4), (n,3,4,4)
    Al = [As[:, i] for i in range(3)]
    P = P_of(Ms)
    plain = kp_plain_density(Al, P)
    w, V, s, imag = eig_frame(Ms)
    real_ok = imag < 1e-9 and bool(np.all(s[:, 0] < 0) and np.all(s[:, 1:] > 0))
    h_eig = kp_h_density(Al, Ms, P, "eigenbasis")
    h_stated = kp_h_density(Al, Ms, P, "stated")
    h_cov = kp_h_density(Al, Ms, P, "covariant")
    rec = {"nprobe": nprobe, "plain_min": float(plain.min()), "plain_max": float(plain.max()), "plain_neg_fraction": float(np.mean(plain < 0)),
           "frame_real_and_signed": real_ok, "eig_max_imag": imag,
           "h_eigenbasis_min": float(h_eig.min()), "h_stated_min": float(h_stated.min()), "h_covariant_min": float(h_cov.min()),
           "stated_vs_eigenbasis_max_rel": float(np.max(np.abs(h_stated - h_eig) / np.maximum(h_eig, 1e-12))),
           "covariant_vs_eigenbasis_max_rel": float(np.max(np.abs(h_cov - h_eig) / np.maximum(h_eig, 1e-12)))}
    # on the vacuum orbit: random Lorentz L, tangent jets and random jets
    orb = {}
    for jets_kind in ("tangent", "random"):
        Ms2, As2 = [], []
        for _ in range(500):
            th = rng.uniform(-0.6, 0.6, 3); ph = rng.uniform(-1, 1, 3)
            L = L_of("boost", th)[()] @ L_of("rot", ph)[()]
            Ms2.append(L @ M_VAC @ L.T)
            Ai = []
            for i in range(3):
                if jets_kind == "tangent":
                    X = sum(rng.normal() * Xg for Xg in GENS.values())
                    Ai.append(L @ tangent(X, M_VAC) @ L.T)
                else:
                    J = rng.normal(size=(4, 4)); Ai.append(0.5 * (J + J.T))
            As2.append(Ai)
        Ms2 = np.array(Ms2); As2 = np.array(As2)
        Al2 = [As2[:, i] for i in range(3)]
        P2 = P_of(Ms2)
        pl = kp_plain_density(Al2, P2); he = kp_h_density(Al2, Ms2, P2, "eigenbasis"); hs = kp_h_density(Al2, Ms2, P2, "stated"); hc = kp_h_density(Al2, Ms2, P2, "covariant")
        orb[jets_kind] = {"plain_vs_eigenbasis_max_rel": float(np.max(np.abs(pl - he) / np.maximum(np.abs(he), 1e-12))),
                          "stated_vs_eigenbasis_max_rel": float(np.max(np.abs(hs - he) / np.maximum(np.abs(he), 1e-12))),
                          "covariant_vs_eigenbasis_max_rel": float(np.max(np.abs(hc - he) / np.maximum(np.abs(he), 1e-12))),
                          "plain_min": float(pl.min()), "stated_min": float(hs.min())}
    rec["vacuum_orbit"] = orb
    # single-point invariance test of the three H forms under a random Lorentz transformation, off the orbit
    Mx = Ms[0]; Ax = [As[0, i] for i in range(3)]
    th = rng.uniform(-0.7, 0.7, 3); L = L_of("boost", th)[()]
    inv = {}
    for var in ("eigenbasis", "stated", "covariant"):
        v0 = float(kp_h_density(Ax, Mx[None], None, var)[0])
        v1 = float(kp_h_density([L @ A @ L.T for A in Ax], (L @ Mx @ L.T)[None], None, var)[0])
        inv[var] = {"lab": v0, "boosted": v1, "rel_change": abs(v1 - v0) / max(abs(v0), 1e-12)}
    inv["plain"] = {"lab": float(kp_plain_density(Ax, P_of(Mx))), "boosted": float(kp_plain_density([L @ A @ L.T for A in Ax], P_of(L @ Mx @ L.T)))}
    rec["lorentz_invariance_probe"] = inv
    # MUTATIONS: H -> eta (plain) already negative; H <-> H^-1 swap = the 'stated' variant; sign flip of the 2uu term
    u = V[:, :, 0]
    eu = ETA @ u[..., None]
    Hm = ETA - 2.0 * eu @ np.swapaxes(eu, -1, -2); Hmi = np.linalg.inv(Hm)
    dm = 0.0
    for A in Al:
        Om = P @ A @ ETA @ P
        dm = dm + 0.5 * np.einsum("...aa->...", np.swapaxes(Om, -1, -2) @ Hm @ Om @ Hmi)
    rec["MUT_flipped_uu_sign_min"] = float(dm.min())
    ok = plain.min() < 0 and h_eig.min() >= -1e-12 and real_ok and orb["tangent"]["plain_vs_eigenbasis_max_rel"] < 1e-8
    rec["verdict"] = "CONFIRMED" if ok else "REFUTED"
    rec["note"] = ("the brief's written formula tr(Om H Om^T H^-1) is NOT the eigenbasis Frobenius norm off the vacuum orbit "
                   f"(max rel dev {rec['stated_vs_eigenbasis_max_rel']:.3g}); tr(Om^T H Om H^-1) is "
                   f"(max rel dev {rec['covariant_vs_eigenbasis_max_rel']:.3g})")
    log(f"C8 {rec['verdict']} plain min {plain.min():.4e} h_eig min {h_eig.min():.3e} stated min {h_stated.min():.3e} cov min {h_cov.min():.3e}")
    log(f"C8 {rec['note']}; orbit {orb}; invariance {inv}")
    save("c8", rec)


# ============================================================ C9: the (1,2) twist family on the lattice
def twist_line(n, L, w, center=0.0):
    h = L / n
    z = (np.arange(n) - (n - 1) / 2.0) * h
    psi = np.clip((z - center + w / 2.0) / w, 0.0, 1.0)
    R = rot_of(G3, psi)
    return conj(R, M_VAC), h, z, psi


def line_energies(M, h):
    P = P_of(M)
    a0 = a0_local_own(M)
    kp = kin = eu = 0.0
    for Al, wt in jets(M, h, naxes=1):
        kp += wt * float(np.sum(kp_plain_density(Al, P)))
        kin += wt * float(np.sum(kin_density(Al, a0)))
        eu += wt * float(np.sum(i1_density(Al)))   # zero: only one spatial jet
    return h * kp, h * kin


def stage_c9():
    rec = {"L": 24.0, "ns": [16, 32, 64, 128]}
    for center_lab, center in (("centered", 0.0), ("offset_h32/2", 0.375)):
        out = {}
        for w in (1.5, 3.0, 6.0):
            rows = {}
            for n in rec["ns"]:
                M, h, z, psi = twist_line(n, 24.0, w, center)
                kp, kin = line_energies(M, h)
                rows[n] = {"h": h, "KP_static_per_area": kp, "kin_I1_local_per_area_times_w": kin * w}
            ps = {}
            ns = rec["ns"]
            for i in range(len(ns) - 1):
                a, b = rows[ns[i]]["KP_static_per_area"], rows[ns[i + 1]]["KP_static_per_area"]
                ps[f"{ns[i]}->{ns[i+1]}"] = float(np.log(a / b) / np.log(2.0)) if a > 0 and b > 0 else None
            e3 = [rows[nn]["KP_static_per_area"] for nn in (16, 32, 64)]
            ps["lsq_16_32_64"] = float(-np.polyfit(np.log([1.5, 0.75, 0.375]), np.log(e3), 1)[0] * -1.0)
            out[str(w)] = {"rows": rows, "p_exponents": ps}
            log(f"C9 {center_lab} w {w}: KP/area {[round(rows[n]['KP_static_per_area'], 6) for n in ns]} p {ps} kin*w(h=0.75) {rows[32]['kin_I1_local_per_area_times_w']:.4f}")
        rec[center_lab] = out
    # MUTATION: an eigenvalue ramp (the (2,3) split, off the orbit) does not vanish with h
    mut = {}
    for n in rec["ns"]:
        h = 24.0 / n
        z = (np.arange(n) - (n - 1) / 2.0) * h
        s = 0.15 * np.clip((z + 1.5) / 3.0, 0, 1)
        M = np.zeros((n, 4, 4)); M[:, 0, 0] = 8.0; M[:, 1, 1] = 1.0; M[:, 2, 2] = 0.15 + s; M[:, 3, 3] = 0.15 - s
        mut[n] = line_energies(M, h)[0]
    rec["MUT_split_ramp_KP_per_area"] = mut
    c = rec["offset_h32/2"]
    kin_w = [c[str(w)]["rows"][32]["kin_I1_local_per_area_times_w"] for w in (1.5, 3.0, 6.0)]
    rec["kin_w_h075_offset"] = kin_w
    rec["kin_w_h075_centered"] = [rec["centered"][str(w)]["rows"][32]["kin_I1_local_per_area_times_w"] for w in (1.5, 3.0, 6.0)]
    rec["producer_kin_w"] = [0.3244, 0.3455, 0.3510]
    rec["producer_p"] = [1.53, 1.89, 1.97]
    vanish = all(c[str(w)]["rows"][128]["KP_static_per_area"] < c[str(w)]["rows"][16]["KP_static_per_area"] / 8 for w in (1.5, 3.0, 6.0))
    rec["checks"] = {"KP_vanishes_with_h": vanish, "kin_w_within_3pct": all(abs(a - b) / b < 0.03 for a, b in zip(kin_w, rec["producer_kin_w"])),
                     "mutant_persists": mut[128] > 0.5 * mut[16]}
    rec["verdict"] = "CONFIRMED" if all(rec["checks"].values()) else "QUALIFIED"
    log(f"C9 {rec['verdict']} {rec['checks']} kin*w {kin_w} mutant {mut}")
    save("c9", rec)


# ============================================================ C10: zigzag concentration on the W3 end fields
def stage_c10():
    rec = {}
    for key, Rs in (("w3_n32_R9_it3000", 9.0), ("w3_n32_R9_it12000", 9.0), ("w3_n48_R15_it3000", 15.0)):
        M, n, L, h = load_field(key)
        X, Y, Z = coords3(n, h)
        r = np.sqrt(X * X + Y * Y + Z * Z)
        region = r < Rs
        shell = (r > Rs - 2.5 * h) & (r < Rs + 0.5 * h)
        outer = (r > Rs + 0.5 * h) & (r < Rs + 3.5 * h)
        P = P_of(M)
        a0 = a0_local_own(M) * region[..., None, None]
        kp_d = np.zeros(M.shape[:3]); kin_d = np.zeros(M.shape[:3])
        for Al, wt in jets(M, h):
            kp_d += wt * kp_plain_density(Al, P)
            kin_d += wt * kin_density(Al, a0)
        kp_d *= h ** 3; kin_d *= h ** 3
        cfg = R13W.cfg_of(n, L)
        kin_cert = float(INS4.kin_of(M, R13W.a0_local(M) * region[..., None, None], cfg))
        rec[key] = {"R_s": Rs, "h": h, "n_shell_cells": int(shell.sum()),
                    "KP_total": float(kp_d.sum()), "KP_shell": float(kp_d[shell].sum()), "KP_shell_frac": float(kp_d[shell].sum() / kp_d.sum()),
                    "KP_outer_frac": float(kp_d[outer].sum() / kp_d.sum()), "KP_min_density": float(kp_d.min()),
                    "kin_total": float(kin_d.sum()), "kin_shell": float(kin_d[shell].sum()), "kin_cert_total": kin_cert,
                    "kin_outer": float(kin_d[outer].sum())}
        log(f"C10 {key}: {rec[key]}")
    prod = {"w3_n32_R9_it3000": (0.879, 2695.5, 2699.3), "w3_n32_R9_it12000": (0.974, 10408.7, 10411.2), "w3_n48_R15_it3000": (0.745, 1233.0, 1290.0)}
    rec["producer"] = prod
    ok = all(abs(rec[k]["KP_shell_frac"] - p[0]) < 0.01 and abs(rec[k]["kin_shell"] - p[1]) / p[1] < 0.01 and abs(rec[k]["kin_total"] - p[2]) / p[2] < 0.01
             for k, p in prod.items())
    rec["verdict"] = "CONFIRMED" if ok else "QUALIFIED"
    save("c10", rec)


# ============================================================ C14
def stage_c14():
    rec = {"det_M_vac": float(np.linalg.det(M_VAC)), "eigs_M_vac": [float(x) for x in np.linalg.eigvalsh(M_VAC)],
           "det_M0_C4": float(np.linalg.det(np.diag([8.0, 1.0, 0.3, 0.1])))}
    try:
        np.linalg.inv(M_VAC); rec["inv_raises"] = False
    except np.linalg.LinAlgError:
        rec["inv_raises"] = True
    # on the hedgehog: the pinned shell is the exact vacuum -> singular cells
    M, n, L, h = load_field("n32_L48_it12000")
    dets = np.linalg.det(M)
    rec["hedgehog_min_abs_det"] = float(np.min(np.abs(dets)))
    rec["hedgehog_singular_cells_frac"] = float(np.mean(np.abs(dets) < 1e-12))
    rec["verdict"] = "CONFIRMED" if rec["inv_raises"] and rec["det_M_vac"] == 0.0 else "REFUTED"
    log(f"C14 {rec}")
    save("c14", rec)


STAGES = {"c1": stage_c1, "c2": stage_c2, "c3": stage_c3, "c4": stage_c4, "c5": stage_c5, "triple": stage_triple,
          "c6": stage_c6, "c7": stage_c7, "c8": stage_c8, "c9": stage_c9, "c10": stage_c10, "c11": stage_c11,
          "c12": stage_c12, "c13": stage_c13, "c14": stage_c14}

if __name__ == "__main__":
    args = sys.argv[1:] or ["all"]
    if args == ["all"]:
        args = ["c1", "c2", "c3", "c11", "c7", "c12", "c14", "c8", "c4", "c5", "triple", "c6", "c13", "c9", "c10"]
    for a in args:
        log(f"=== stage {a}")
        STAGES[a]()
    log("done")
