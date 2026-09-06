"""M5.32 R15 shared instruments (ledger 6.4): the DEGENERATE vacuum, the split
stiffness, the exact-projector kinetic term K_P^23, the I1^h read, and a FIRE
with callables.  Every certified piece is consumed from the certified stack
(m5_21_3_a_4d.py, read-only; two of its module-level functions are PATCHED on
every loaded instance: c4_of reads the trace targets from cfg["cp"] when present,
vac4 reads the fourth eigenvalue from cfg["lam4"] when present).

EQUATIONS FIRST
---------------
Conventions (the run's locked table, R14): M real symmetric 4x4 per cell, eta =
diag(-1, 1, 1, 1), N = M eta; code branch s = -1: the certified vacuum
M = diag(g, 1, delta, 0) has the N-spectrum (-g, 1, delta, 0).  The author's
09-05 object lives on the DEGENERATE vacuum M = diag(g, 1, delta, delta),
N-spectrum (-g, 1, delta, delta), reached by the trace targets
    C_p = (-g)^p + 1 + 2 delta^p                              (cfg["cp"])
in the certified form V4 = W1 sum_p (tr N^p - C_p)^2 (INS4.e_parts / grad, patched).

The split stiffness, through symmetric functions of the pair (no pair eigenvectors,
smooth at the degenerate point):
    s = tr N - lambda_g - lambda_1,   p = det N / (lambda_g lambda_1)   (= lambda_2 lambda_3),
    (lambda_2 - lambda_3)^2 = s^2 - 4 p,        E_split = mu h^3 sum_cells (s^2 - 4 p),
    lambda_g, lambda_1 the two ISOLATED eigenvalues (the most negative and the largest of the
    sorted spectrum), read per cell and Newton-polished on the characteristic polynomial
    chi(x) = x^4 - e1 x^3 + e2 x^2 - e3 x + e4 (holomorphic in M: complex-step exact).

The exact spectral projector onto the (2,3) eigenplane and the kinetic term:
    P_g = (N - lambda_1)(N^2 - s N + p) / [(lambda_g - lambda_1)(lambda_g^2 - s lambda_g + p)]
    P_1 = (N - lambda_g)(N^2 - s N + p) / [(lambda_1 - lambda_g)(lambda_1^2 - s lambda_1 + p)]
    P23 = I - P_g - P_1          (= the author's (N - g)(N - 1)/[(l23 - g)(l23 - 1)] at l2 = l3)
    Om_mu = P23 A_mu eta P23,    K_P^23 E-density = (1/2)[sum_i tr(Om_i^T eta Om_i eta)
                                                    + omega^2 tr(Om_0^T eta Om_0 eta)]
    THEOREM (checked in the selftest): tr(Om^T H Om H^-1) = tr(Om^T eta Om eta) for
    H = eta + 2 (eta u)(eta u)^T, because P23 u = 0 and (eta u)^T P23 = 0: on the
    projected block the author's H-adjoint form and the eta form coincide.
Gradient (exact): dE = tr(Y dOm), Y = eta Om^T eta; dOm = dP X P + P dX P + P X dP,
X = A eta;  dE/dA = (eta P Y P)^T;  through the projector, with Z = X P Y + Y P X and
    dP23 = sum_{j in {g,1}} [R_j dN P_j + P_j dN R_j],
    R_j  = P23 [(s - 2 lambda_j) I - (N - lambda_j)] P23 / (p - s lambda_j + lambda_j^2)
(the inverse of (N - lambda_j) restricted to the pair block, Cayley-Hamilton; smooth
at the degenerate point):  dE/dN = W^T, W = sum_j (P_j Z R_j + R_j Z P_j);
dN = dM eta so dE/dM = (eta W)^T, symmetrized.  Split: d(s^2 - 4p) = tr(W_s dN),
    W_s = 2 s P23 - 4 [ (tr N) I - N - lambda_1 P_g - lambda_g P_1 - s (P_g + P_1) - (lambda_g + lambda_1) P23 ].

L_P (the author's object under our reading) = E_u(-4 I1, certified) + V4^dd + mu split + c_P K_P^23,
fixed J:  E_J = E_stat + J^2 / (4 kin_tot),  kin_tot = kin_I1 + c_P kin_KP23 (a0 = the local clock,
frozen inside the gradient, the R13-W convention).

Selftests (python3 m5_32_r15_common.py): the projector identities (P23^2 = P23, P23 u = 0,
rank 2, the author's formula at the degenerate point), the H-adjoint equality, the vacuum
facts (V4^dd = 0, a0_local = 0, both inertias 0 on the degenerate vacuum), complex-step
and central-difference gradient gates for K_P^23 static, K_P^23 kinetic and the split term,
covariance under global boosts and rotations with the no-eta control FAILING, and the
E_u == 4 x (the I1 Lagrangian read at omega 0) identity that fixes the E-orientation of the I1^h read.
"""
from __future__ import annotations
import importlib.util
import json
import os
import sys
import time
import types

import numpy as np

ARGV = list(sys.argv)
HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.dirname(HERE)
DATA, PLOTS = os.path.join(RES, "data"), os.path.join(RES, "plots")
CK = os.path.join(RES, "checkpoints", "m5_32_r15")
os.makedirs(CK, exist_ok=True)
T0 = time.time()
ETA = np.diag([-1.0, 1.0, 1.0, 1.0])
EYE = np.eye(4)
G, DELTA, S, W1 = 8.0, 0.3, -1.0, 0.000724023879


def log(m):
    print(f"[{time.time() - T0:8.1f}s] {m}", flush=True)


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


INS4 = _load("m5_21_3_a_4d", "m5_21_3_a_4d.py")
C13 = _load("m5_32_r13w_common", "m5_32_r13w_common.py")
T14 = _load("m5_32_r14_terms", "m5_32_r14_terms.py")
EXT = _load("m5_32_terms_ext", "m5_32_terms_ext.py")
LAG = C13.LAG
B8 = C13.B8


# ------------------------------------------------ the patch: trace targets + fourth eigenvalue from cfg
def _patch_stack(root_mods):
    seen = set()
    patched = []

    def visit(m, depth):
        if id(m) in seen or depth > 4:
            return
        seen.add(id(m))
        f = getattr(m, "__file__", "") or ""
        if f.endswith("m5_21_3_a_4d.py") and not getattr(m, "_R15_PATCHED", False):
            orig_c4, orig_vac = m.c4_of, m.vac4
            m.c4_of = (lambda cfg, _o=orig_c4: tuple(cfg["cp"]) if "cp" in cfg else _o(cfg))
            m.vac4 = (lambda cfg, _o=orig_vac: (np.diag([-cfg["sg"], 1.0, cfg["delta"], cfg["lam4"]])
                                                if "lam4" in cfg else _o(cfg)))
            m._R15_PATCHED = True
            patched.append(m)
        for v in list(vars(m).values()):
            if isinstance(v, types.ModuleType):
                visit(v, depth + 1)
    for r in root_mods:
        visit(r, 0)
    return patched


PATCHED = _patch_stack([INS4, C13, T14, EXT, LAG, B8])


def cp_dd(g=G, delta=DELTA):
    return tuple((-g) ** p + 1.0 + 2.0 * delta ** p for p in range(1, 5))


def cfg_dd(n, L, mu=0.0, cP=0.0, g=G, delta=DELTA):
    """the degenerate-vacuum configuration: trace targets of diag(g, 1, delta, delta), the
    fourth eigenvalue delta, the split stiffness mu and the K_P^23 coefficient c_P."""
    cfg = INS4.base_cfg(s=S, g=g, n=n, L=float(L), delta=delta)
    cfg["cp"] = cp_dd(g, delta)
    cfg["lam4"] = delta
    cfg["mu"] = float(mu)
    cfg["cP"] = float(cP)
    return cfg


def cfg_cert(n, L, g=G, delta=DELTA):
    return INS4.base_cfg(s=S, g=g, n=n, L=float(L), delta=delta)


# ------------------------------------------------ spectrum parts and projectors
def spectrum_parts(M):
    """per cell: lambda_g (most negative), lambda_1 (largest), s = l2 + l3, p = l2 l3 of the
    pair (the two middle eigenvalues of N = M eta sorted by real part).  Complex-step safe:
    the two isolated roots are polished by Newton steps on the characteristic polynomial
    (a holomorphic map of the traces; LAPACK's eigvals is only the starting point), and
    s = e1 - lg - l1 and p = det N / (lg l1) (holomorphic, no near-degenerate product)."""
    N = M @ ETA
    lam = np.linalg.eigvals(N)
    order = np.argsort(lam.real, axis=-1)
    lam = np.take_along_axis(lam, order, axis=-1)
    lg, l1 = lam[..., 0].real, lam[..., 3].real
    if np.iscomplexobj(M):
        lg, l1 = lg.astype(complex), l1.astype(complex)
    N2 = N @ N
    t1, t2 = _tr(N), _tr(N2)
    t3 = _tr(N2 @ N)
    e1 = t1
    e2 = (t1 * t1 - t2) / 2.0
    e3 = (t1 ** 3 - 3.0 * t1 * t2 + 2.0 * t3) / 6.0
    e4 = np.linalg.det(N)                      # no cancellation (the trace form loses 4 digits)

    def polish(x):
        for _ in range(4):
            chi = (((x - e1) * x + e2) * x - e3) * x + e4
            dchi = ((4.0 * x - 3.0 * e1) * x + 2.0 * e2) * x - e3
            x = x - chi / dchi
        return x
    lg, l1 = polish(lg), polish(l1)
    s = e1 - lg - l1
    p = e4 / (lg * l1)                        # l2 l3 = det N / (lg l1), exact to roundoff
    return N, lg, l1, s, p


def projectors(M):
    """(N, P_g, P_1, P23, R_g, R_1, lg, l1, s, p) per cell, all matrix polynomials in N."""
    N, lg, l1, s, p = spectrum_parts(M)
    I = np.broadcast_to(EYE, N.shape)
    Q = N @ N - s[..., None, None] * N + p[..., None, None] * I          # (N - l2)(N - l3)
    Pg = ((N - l1[..., None, None] * I) @ Q) / ((lg - l1) * (lg * lg - s * lg + p))[..., None, None]
    P1 = ((N - lg[..., None, None] * I) @ Q) / ((l1 - lg) * (l1 * l1 - s * l1 + p))[..., None, None]
    P23 = I - Pg - P1
    Rs = []
    for lj in (lg, l1):
        B = (s - 2.0 * lj)[..., None, None] * I - (N - lj[..., None, None] * I)
        Rs.append(P23 @ B @ P23 / (p - s * lj + lj * lj)[..., None, None])
    return N, Pg, P1, P23, Rs[0], Rs[1], lg, l1, s, p


def _tr(A):
    return np.einsum("...aa->...", A)


def kp23_cells(A_list, M, need_grad=True):
    """per cell: E = (1/2) sum_A tr(Om^T eta Om eta), Om = P23 (A eta) P23 for each A in
    A_list (held fixed); returns (E_cells, dE/dM_cells at fixed A, [dE/dA]).  The M-gradient
    is through the projector only (the jets are chained by the caller)."""
    N, Pg, P1, P, Rg, R1, lg, l1, s, p = projectors(M)
    E = np.zeros(M.shape[:-2], dtype=M.dtype)
    Zsum = np.zeros_like(M)
    dA_out = []
    for A in A_list:
        X = A @ ETA
        Om = P @ X @ P
        Y = ETA @ np.swapaxes(Om, -1, -2) @ ETA
        E = E + 0.5 * _tr(Y @ Om)
        if need_grad:
            dA_out.append(np.swapaxes(ETA @ P @ Y @ P, -1, -2))
            Zsum = Zsum + X @ P @ Y + Y @ P @ X
    if not need_grad:
        return E, None, None
    W = Pg @ Zsum @ Rg + Rg @ Zsum @ Pg + P1 @ Zsum @ R1 + R1 @ Zsum @ P1
    Gm = np.swapaxes(ETA @ W, -1, -2)
    Gm = 0.5 * (Gm + np.swapaxes(Gm, -1, -2))
    return E, Gm, dA_out


def kp23_energy_grad(M, cfg, a0=None, need_grad=True):
    """K_P^23 on the certified stencil: (E_stat, G_stat, kin, G_kin), h^3-weighted; kin is the
    omega^2 coefficient for A_0 = a0 (frozen a0 in G_kin)."""
    h = cfg["h"]
    h3 = h ** 3
    E_stat, G_stat = 0.0, (np.zeros_like(M) if need_grad else None)
    for br, wt in INS4.branches(cfg["stencil"]):
        A = [INS4.d1(M, ax, h, br) for ax in range(3)]
        Ec, Gc, dA = kp23_cells(A, M, need_grad)
        E_stat = E_stat + wt * np.sum(Ec)
        if need_grad:
            G_stat += wt * Gc
            for ax in range(3):
                G_stat += wt * INS4.d1_adj(dA[ax], ax, h, br)
    kin, G_kin = None, None
    if a0 is not None:
        Ek, Gk, _ = kp23_cells([a0], M, need_grad)
        kin = h3 * np.sum(Ek)
        G_kin = (h3 * Gk) if need_grad else None
    return h3 * E_stat, (h3 * G_stat if need_grad else None), kin, G_kin


def kp23_static_density(M, cfg):
    h = cfg["h"]
    d = 0.0
    for br, wt in INS4.branches(cfg["stencil"]):
        A = [INS4.d1(M, ax, h, br) for ax in range(3)]
        Ec, _, _ = kp23_cells(A, M, need_grad=False)
        d = d + wt * Ec
    return h ** 3 * d


def split_cells(M, need_grad=True):
    """(lambda_2 - lambda_3)^2 = s^2 - 4 p per cell and its M-gradient."""
    N, Pg, P1, P23, Rg, R1, lg, l1, s, p = projectors(M)
    val = s * s - 4.0 * p
    if not need_grad:
        return val, None
    I = np.broadcast_to(EYE, N.shape)
    t1 = _tr(N)
    Ws = (2.0 * s)[..., None, None] * P23 - 4.0 * (t1[..., None, None] * I - N - l1[..., None, None] * Pg
                                                   - lg[..., None, None] * P1 - s[..., None, None] * (Pg + P1)
                                                   - (lg + l1)[..., None, None] * P23)
    Gm = np.swapaxes(ETA @ Ws, -1, -2)
    return val, 0.5 * (Gm + np.swapaxes(Gm, -1, -2))


def split_energy_grad(M, cfg, need_grad=True):
    h3 = cfg["h"] ** 3
    v, g = split_cells(M, need_grad)
    return cfg["mu"] * h3 * np.sum(v), (cfg["mu"] * h3 * g if need_grad else None)


def sorted_spectrum(M):
    lam = np.linalg.eigvals(M @ ETA).real
    return np.sort(lam, axis=-1)


# ------------------------------------------------ the full L_P energy and gradient
def lp_parts(M, cfg, a0=None):
    """dict of the h^3-weighted parts: E_u, V4, split, KP (c_P-weighted), E_stat, kin_I1, kin_KP, kin_tot."""
    e_u, e_v = INS4.e_parts(M, cfg)
    es, _ = split_energy_grad(M, cfg, need_grad=False)
    ek, _, kk, _ = kp23_energy_grad(M, cfg, a0, need_grad=False)
    num = (lambda x: complex(x)) if np.iscomplexobj(M) else (lambda x: float(x))
    out = {"E_u": num(e_u), "V4": num(e_v), "split": num(es), "KP": num(cfg["cP"] * ek),
           "KP_static_raw": num(ek)}
    out["E_stat"] = out["E_u"] + out["V4"] + out["split"] + out["KP"]
    if a0 is not None:
        ki = num(INS4.kin_of(M, a0, cfg))
        out["kin_I1"], out["kin_KP_raw"] = ki, num(kk)
        out["kin_tot"] = ki + cfg["cP"] * num(kk)
    return out


def lp_grad(M, cfg):
    Gt = INS4.grad(M, cfg)
    if cfg["mu"] != 0.0:
        Gt = Gt + split_energy_grad(M, cfg)[1]
    if cfg["cP"] != 0.0:
        Gt = Gt + cfg["cP"] * kp23_energy_grad(M, cfg)[1]
    return Gt


def lp_kin_grad(M, a0, cfg):
    Gk = INS4.kin_grad(M, a0, cfg)
    if cfg["cP"] != 0.0:
        Gk = Gk + cfg["cP"] * kp23_energy_grad(M, cfg, a0)[3]
    return Gk


def fire_lp(M0, cfg, free_mask, max_iter, J=None, a0_of=None, log_every=100, tag="", f_tol=1e-6,
            plateau=(2000, 1e-10), dt0=0.01, dt_max=0.1, diag=None, ck_path=None, ck_every=500):
    """C13.fire_proj verbatim in its FIRE logic on the L_P energy (no projection hook);
    J with a0_of adds J^2 / (4 kin_tot), a0 refreshed each step and frozen in the gradient;
    ck_path: the field is saved there every ck_every iterations (resume-complete)."""
    M = M0.copy()
    free = free_mask[..., None, None].astype(float)
    v = np.zeros_like(M)
    dt, alpha, n_up = dt0, 0.1, 0
    hist = []

    def parts(Mx):
        a0 = a0_of(Mx) if J is not None else None
        pp = lp_parts(Mx, cfg, a0)
        if J is not None:
            pp["E_J"] = pp["E_stat"] + J * J / (4.0 * pp["kin_tot"])
            pp["omega"] = J / (2.0 * pp["kin_tot"])
        return pp

    def tot_grad(Mx):
        Gt = lp_grad(Mx, cfg)
        if J is not None:
            a0 = a0_of(Mx)
            k = float(INS4.kin_of(Mx, a0, cfg)) + cfg["cP"] * float(kp23_energy_grad(Mx, cfg, a0, need_grad=False)[2])
            Gt = Gt - (J * J / (4.0 * k * k)) * lp_kin_grad(Mx, a0, cfg)
        return Gt

    F = -tot_grad(M) * free
    t0 = time.time()
    stop = "max_iter"
    it = 0
    for it in range(1, max_iter + 1):
        Pw = float(np.sum(F * v))
        if Pw > 0.0:
            n_up += 1
            vn = np.sqrt(np.sum(v * v))
            fn = np.sqrt(np.sum(F * F))
            v = (1 - alpha) * v + alpha * (F / max(fn, 1e-300)) * vn
            if n_up > 5:
                dt = min(dt * 1.1, dt_max)
                alpha *= 0.99
        else:
            v[:] = 0.0
            dt *= 0.5
            alpha = 0.1
            n_up = 0
        v += dt * F
        M += dt * v
        F = -tot_grad(M) * free
        fmax = float(np.max(np.abs(F)))
        if not np.isfinite(fmax):
            stop = "non-finite"
            break
        if it % log_every == 0 or it == max_iter:
            pp = parts(M)
            row = {"it": it, "fmax": fmax, "dt": dt}
            row.update(pp)
            if diag is not None:
                row.update(diag(M))
            hist.append(row)
            key = "E_J" if J is not None else "E_stat"
            extra = f" kin {pp['kin_tot']:10.4f} om {pp['omega']:.5f}" if J is not None else ""
            print(f"  {tag} it {it:6d} {key} {pp[key]:14.6f} E_u {pp['E_u']:10.5f} V4 {pp['V4']:.4e} split {pp['split']:.4e} "
                  f"KP {pp['KP']:.4e} fmax {fmax:.3e}{extra} [{time.time() - t0:.0f}s]", flush=True)
            back = max(1, plateau[0] // max(log_every, 1))
            if len(hist) > back and abs(row[key] - hist[-1 - back][key]) < plateau[1]:
                stop = "plateau"
                break
            if ck_path is not None and it % ck_every == 0:
                np.save(ck_path, M)
        if fmax < f_tol:
            stop = "f_tol"
            break
    return M, {"stop": stop, "trace": hist, "wall_s": round(time.time() - t0, 1), "iters": it}


# ------------------------------------------------ seeds and reads
def seed_uniaxial(cfg):
    """the radial hedgehog of the degenerate vacuum: the leading eigenvector (eigenvalue 1)
    radial, the (2,3) pair degenerate (nothing to wind); B8.dressed at zero rapidity on the
    patched vacuum diag(g, 1, delta, delta)."""
    return B8.dressed(cfg, 0.0)


def i1h_static(M, cfg):
    """the static energy of the author's kinetic term -4 I1^h read on a field, E-orientation:
    L = -4 I1^h gives E_static = +4 x the Lagrangian read at omega = 0 (the same orientation as
    E_u = +4 x the I1 read, checked in the selftest)."""
    T = EXT.REGISTRY_EXT["I1_h"]
    p = LAG.default_params(s=S, g=cfg["g"], delta=cfg["delta"])
    return 4.0 * float(LAG.term_lagrangian(T, M, cfg, p))


def i1_static_registry(M, cfg):
    T = LAG.REGISTRY["I1"]
    p = LAG.default_params(s=S, g=cfg["g"], delta=cfg["delta"])
    return 4.0 * float(LAG.term_lagrangian(T, M, cfg, p))


def h_form_check(M, A):
    """tr(Om^T H Om H^-1) with H from the timelike eigenvector, against tr(Om^T eta Om eta)."""
    lam, V, sig = T14.eig_N(M)
    k = np.argmin(sig, axis=-1)
    u = np.take_along_axis(V, k[..., None, None], axis=-1)[..., :, 0]          # (..., 4) timelike eigenvector
    eu = np.einsum("ab,...b->...a", ETA, u)
    H = ETA + 2.0 * eu[..., :, None] * eu[..., None, :]
    Hi = ETA + 2.0 * u[..., :, None] * u[..., None, :]
    _, _, _, P, _, _, _, _, _, _ = projectors(M)
    Om = P @ (A @ ETA) @ P
    OmT = np.swapaxes(Om, -1, -2)
    return _tr(OmT @ H @ Om @ Hi), _tr(OmT @ ETA @ Om @ ETA), np.max(np.abs(np.einsum("...ab,...b->...a", P, u)))


# ------------------------------------------------ selftests
def _random_field(rng, n, cfg, amp=0.25):
    """a smooth random symmetric field around the degenerate vacuum (the pair split, tilted,
    boosted): vac + amp * smooth noise, symmetrized."""
    X = rng.normal(size=(n, n, n, 4, 4))
    for ax in range(3):
        X = 0.5 * (X + np.roll(X, 1, axis=ax))
        X = 0.5 * (X + np.roll(X, -1, axis=ax))
    X = 0.5 * (X + np.swapaxes(X, -1, -2))
    return INS4.vac4(cfg)[None, None, None] + amp * X


def selftest(write=True):
    res, lines = {}, []
    rng = np.random.default_rng(1505)
    n = 6
    cfg = cfg_dd(n, 9.0, mu=1e-2, cP=1.0)
    M = _random_field(rng, n, cfg)

    def check(name, ok, val):
        res[name] = {"ok": bool(ok), "value": val}
        lines.append(f"{'PASS' if ok else 'FAIL'} {name}: {val}")
        log(lines[-1])

    # 1. projector identities
    N, Pg, P1, P, Rg, R1, lg, l1, s, p = projectors(M)
    check("P23 idempotent", np.max(np.abs(P @ P - P)) < 1e-9, float(np.max(np.abs(P @ P - P))))
    check("P23 rank 2 (trace)", np.max(np.abs(_tr(P) - 2.0)) < 1e-9, float(np.max(np.abs(_tr(P) - 2.0))))
    check("P23 commutes with N", np.max(np.abs(P @ N - N @ P)) < 1e-8, float(np.max(np.abs(P @ N - N @ P))))
    # the author's formula at an exactly degenerate point
    Md = INS4.vac4(cfg)[None, None, None] + 0.0 * M
    Nd = Md @ ETA
    auth = (Nd - (-G) * EYE) @ (Nd - 1.0 * EYE) / ((DELTA + G) * (DELTA - 1.0))
    _, _, _, Pd, _, _, _, _, _, _ = projectors(Md)
    check("P23 == the author's formula at lambda_2 = lambda_3 (vacuum)", np.max(np.abs(Pd - auth)) < 1e-12, float(np.max(np.abs(Pd - auth))))
    # 2. the H-adjoint equality and P23 u = 0
    A = rng.normal(size=M.shape); A = 0.5 * (A + np.swapaxes(A, -1, -2))
    hf, ef, pu = h_form_check(M, A)
    check("H-adjoint form == eta form on the projected block", np.max(np.abs(hf - ef)) < 1e-10 * max(1.0, float(np.max(np.abs(ef)))), float(np.max(np.abs(hf - ef))))
    check("P23 u = 0 (timelike eigenvector killed)", pu < 1e-10, float(pu))
    # 3. vacuum facts
    Mv = np.broadcast_to(INS4.vac4(cfg), M.shape).copy()
    ev = INS4.e_parts(Mv, cfg)
    a0v = C13.a0_local(Mv)
    check("V4^dd = 0 and E_u = 0 on the degenerate vacuum", abs(ev[0]) < 1e-14 and abs(ev[1]) < 1e-14, [float(ev[0]), float(ev[1])])
    check("a0_local = 0 on the degenerate vacuum (the (2,3) rotation is a symmetry)", np.max(np.abs(a0v)) < 1e-12, float(np.max(np.abs(a0v))))
    a0G1 = C13.a0_G1(Mv)
    check("kin_I1 = 0 and kin_KP23 = 0 on the degenerate vacuum for the G1 clock", abs(INS4.kin_of(Mv, a0G1, cfg)) < 1e-14 and abs(kp23_energy_grad(Mv, cfg, a0G1, need_grad=False)[2]) < 1e-14,
          [float(INS4.kin_of(Mv, a0G1, cfg)), float(kp23_energy_grad(Mv, cfg, a0G1, need_grad=False)[2])])
    sv, _ = split_cells(Mv, need_grad=False)
    check("split = 0 on the degenerate vacuum", np.max(np.abs(sv)) < 1e-14, float(np.max(np.abs(sv))))
    # 4. gradient gates: complex step + central differences on a random direction
    D = rng.normal(size=M.shape); D = 0.5 * (D + np.swapaxes(D, -1, -2))
    a0 = C13.a0_local(M)
    for name, fn_e, fn_g in (
        ("K_P^23 static", lambda X: kp23_energy_grad(X, cfg, need_grad=False)[0], lambda X: kp23_energy_grad(X, cfg)[1]),
        ("K_P^23 kinetic (frozen a0)", lambda X: kp23_energy_grad(X, cfg, a0, need_grad=False)[2], lambda X: kp23_energy_grad(X, cfg, a0)[3]),
        ("split term", lambda X: split_energy_grad(X, cfg, need_grad=False)[0], lambda X: split_energy_grad(X, cfg)[1]),
        ("L_P total static", lambda X: lp_parts(X, cfg)["E_stat"], lambda X: lp_grad(X, cfg)),
    ):
        g = fn_g(M)
        an = float(np.sum(g * D))
        eps = 1e-6
        fd = (float(fn_e(M + eps * D)) - float(fn_e(M - eps * D))) / (2 * eps)
        cs = None
        try:
            cs = float(np.imag(fn_e(M + 1e-20j * D)) / 1e-20)
        except Exception as e:                                                          # noqa: BLE001
            cs = f"complex step unavailable: {e!r}"
        rel_fd = abs(an - fd) / max(abs(fd), 1e-300)
        rel_cs = (abs(an - cs) / max(abs(cs), 1e-300)) if isinstance(cs, float) else None
        ok = rel_fd < 1e-6 and (rel_cs is None or rel_cs < 1e-8)
        check(f"gradient gate {name}: analytic vs central FD vs complex step", ok, {"analytic": an, "fd": fd, "cs": cs, "rel_fd": rel_fd, "rel_cs": rel_cs})
    # 5. covariance under global Lorentz maps (constant maps commute with the stencil)
    def lorentz(kind):
        if kind == "boost":
            K = np.zeros((4, 4)); K[0, 2] = K[2, 0] = 1.0
            b = 0.3
            return EYE + np.sinh(b) * K + (np.cosh(b) - 1.0) * (K @ K)
        Gm = np.zeros((4, 4)); Gm[1, 3], Gm[3, 1] = -1.0, 1.0
        q = 0.7
        return EYE + np.sin(q) * Gm + (1 - np.cos(q)) * (Gm @ Gm)
    for kind in ("boost", "rotation"):
        Lm = lorentz(kind)
        ML = Lm @ M @ Lm.T
        e0 = kp23_energy_grad(M, cfg, need_grad=False)[0]
        e1 = kp23_energy_grad(ML, cfg, need_grad=False)[0]
        s0 = split_energy_grad(M, cfg, need_grad=False)[0]
        s1 = split_energy_grad(ML, cfg, need_grad=False)[0]
        check(f"covariance of K_P^23 and the split term under a global {kind}", abs(e1 - e0) < 1e-9 * abs(e0) and abs(s1 - s0) < 1e-9 * abs(s0), [float(e0), float(e1), float(s0), float(s1)])
    # the no-eta control: the plain Frobenius form tr(Om^T Om) must FAIL under a boost
    Lm = lorentz("boost")
    ML = Lm @ M @ Lm.T

    def frob(X):
        tot = 0.0
        for br, wt in INS4.branches(cfg["stencil"]):
            A = [INS4.d1(X, ax, cfg["h"], br) for ax in range(3)]
            _, _, _, P, _, _, _, _, _, _ = projectors(X)
            for Ai in A:
                Om = P @ (Ai @ ETA) @ P
                tot += wt * 0.5 * np.sum(Om * Om)
        return tot
    f0, f1 = frob(M), frob(ML)
    check("no-eta control (plain Frobenius) FAILS covariance under the boost", abs(f1 - f0) > 1e-3 * abs(f0), [float(f0), float(f1)])
    # 6. the sign of the registry read: E_u == +4 x the I1 Lagrangian read (E-orientation)
    eu = INS4.e_parts(M, cfg)[0]
    ir = i1_static_registry(M, cfg)
    check("E_u == 4 I1 (registry read at omega 0, fixes the I1^h E-orientation)", abs(eu - ir) < 1e-9 * abs(eu), [float(eu), float(ir)])
    ih = i1h_static(M, cfg)
    res["I1h_static_on_random_field"] = float(ih)
    log(f"the static energy of the kinetic term -4 I1^h (= 4 x the I1^h read) on the random field: {ih:.6f} (E_u = 4 I1: {eu:.6f})")
    # 7. the stencil patch reached every stack instance
    check("the certified stack patched on every loaded instance", len(PATCHED) >= 1 and all(getattr(m, "_R15_PATCHED", False) for m in PATCHED), len(PATCHED))
    res["n_pass"] = sum(1 for v in res.values() if isinstance(v, dict) and v.get("ok"))
    res["n_total"] = sum(1 for v in res.values() if isinstance(v, dict) and "ok" in v)
    log(f"selftest {res['n_pass']}/{res['n_total']}")
    if write:
        json.dump(res, open(os.path.join(DATA, "m5_32_r15_common_selftest.json"), "w"), indent=1, default=float)
    return res


if __name__ == "__main__":
    selftest()
