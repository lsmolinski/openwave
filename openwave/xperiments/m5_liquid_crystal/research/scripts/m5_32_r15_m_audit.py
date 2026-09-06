"""M5.32 R15-M adversarial audit: the degenerate-vacuum object L_P on the lattice.

Independent implementation from the DEFINITIONS only (the producer's scripts and
summaries were not opened).  Consumes: the certified stack `m5_21_3_a_4d.py`
(base_cfg, coords, pin_shell, e_parts (E_u only), a_fields, comm_eta, inner_eta,
kin_of, sym4, d1, branches), `m5_21_8_b_lattice.py` (G1, G2, G3, rot_field),
`m5_32_r13w_common.py` (a0_local, kin_density, seed_hedgehog), and the END FIELDS
`checkpoints/m5_32_r15/m_hedgehog/relax_n{n}_L{L}_mu{mu}_cP{cP}.npy`.

Own objects: V4^dd (trace targets of the degenerate spectrum (-g, 1, delta, delta)),
SPLIT, the (2,3) spectral projector (eigendecomposition route, cross-checked against
the Lagrange-product route), K_P^23 (static and clock), the E_stat density, shells,
tail fits, drift, the M-a Hessian (finite differences AND the Gauss-Newton closed
form), and a 200-step FIRE mutation with an OWN analytic gradient (complex-step
gated) of E_u + V4^dd.

Modes: all | hess | reads | mutate | ref | embed
Out: ../data/m5_32_r15_m_audit.json
"""
from __future__ import annotations

import sys
ARGV = sys.argv[1:]                    # captured before any import side effect
import os
os.environ.setdefault("OMP_NUM_THREADS", "2")
import importlib.util
import json
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.dirname(HERE)
DATA = os.path.join(RES, "data")
CK = os.path.join(RES, "checkpoints", "m5_32_r15", "m_hedgehog")
OUT = os.path.join(DATA, "m5_32_r15_m_audit.json")
T0 = time.time()


def log(m):
    print(f"[{time.time() - T0:7.1f}s] {m}", flush=True)


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


INS4 = _load("m5_21_3_a_4d", "m5_21_3_a_4d.py")
B8 = _load("m5_21_8_b_lattice", "m5_21_8_b_lattice.py")
R13 = _load("m5_32_r13w_common", "m5_32_r13w_common.py")

G, DELTA, W1 = 8.0, 0.3, 0.000724023879
ETA = np.diag([-1.0, 1.0, 1.0, 1.0])
SGN = np.diag(ETA)                                   # (-1, 1, 1, 1)
D_DEG = np.diag([G, 1.0, DELTA, DELTA])              # the degenerate vacuum M
C_DD = tuple((-G) ** p + 1.0 + 2.0 * DELTA ** p for p in range(1, 5))
VAC_SPEC = np.array([-G, DELTA, DELTA, 1.0])         # sorted N-spectrum


def cfg_of(n, L):
    return INS4.base_cfg(s=-1.0, g=G, n=n, L=float(L), delta=DELTA)


# ============================================================ own densities
def inner(F, Gm):
    """<F, G>_eta = tr(eta F eta G^T) = sum_ab eta_a eta_b F_ab G_ab (complex-safe)."""
    return np.einsum("...ab,...ab,a,b->...", F, Gm, SGN, SGN)


def traces(M):
    N = M @ ETA
    P = N
    t = [np.einsum("...kk->...", P)]
    for _ in range(3):
        P = P @ N
        t.append(np.einsum("...kk->...", P))
    return t


def v4dd_density(M):
    """W1 sum_p (tr N^p - C_p)^2 per cell (unweighted)."""
    t = traces(M)
    return W1 * sum((t[p] - C_DD[p]) ** 2 for p in range(4))


def spectral(M):
    """eigendecomposition of N = M eta per cell, sorted by real part.
    Returns lam (complex, sorted), P23 (real part), max |Im lam|, max |Im P23|."""
    N = M @ ETA
    w, V = np.linalg.eig(N)
    order = np.argsort(w.real, axis=-1)
    w = np.take_along_axis(w, order, axis=-1)
    V = np.take_along_axis(V, order[..., None, :], axis=-1)
    Vinv = np.linalg.inv(V)
    P23 = V[..., :, 1:3] @ Vinv[..., 1:3, :]
    return w, P23.real, float(np.max(np.abs(w.imag))), float(np.max(np.abs(P23.imag)))


def p23_lagrange(M, lam):
    """the definition's route: P23 = I - P_g - P_1 with Lagrange projectors."""
    N = M @ ETA
    I = np.broadcast_to(np.eye(4), N.shape)
    l = lam.real
    lg, l2, l3, l1 = l[..., 0], l[..., 1], l[..., 2], l[..., 3]

    def proj(lk, others):
        P = I.copy()
        for lo in others:
            P = P @ (N - lo[..., None, None] * np.eye(4)) / (lk - lo)[..., None, None]
        return P
    Pg = proj(lg, (l1, l2, l3))
    P1 = proj(l1, (lg, l2, l3))
    return I - Pg - P1


def split_density(M, lam=None):
    if lam is None:
        lam = spectral(M)[0]
    d = lam[..., 2] - lam[..., 1]
    return np.abs(d) ** 2


def eu_density(M, cfg):
    """4 sum_{i<j} <F_ij, F_ij>_eta, branch-weighted (unweighted by h^3)."""
    dens = np.zeros(M.shape[:3], dtype=M.dtype)
    for br, wt in INS4.branches(cfg["stencil"]):
        A = [INS4.d1(M, ax, cfg["h"], br) for ax in range(3)]
        for i in range(3):
            for j in range(i + 1, 3):
                F = A[i] @ ETA @ A[j] - A[j] @ ETA @ A[i]
                dens = dens + wt * 4.0 * inner(F, F)
    return dens


def kp23_density(M, cfg, P23):
    """(1/2) sum_i tr(Om_i^T eta Om_i eta), Om_i = P23 A_i eta P23 (unweighted)."""
    dens = np.zeros(M.shape[:3])
    for br, wt in INS4.branches(cfg["stencil"]):
        for ax in range(3):
            A = INS4.d1(M, ax, cfg["h"], br)
            Om = P23 @ A @ ETA @ P23
            dens += wt * 0.5 * inner(Om, Om)
    return dens


def kp_clock(M, cfg, P23, a0):
    Om = P23 @ a0 @ ETA @ P23
    return cfg["h"] ** 3 * float(np.sum(0.5 * inner(Om, Om)))


# ============================================================ reads
def shells(dens_h3, r, L, nsh=12):
    """energy per radial shell (12 shells to L/2) and the shell-mean density."""
    w = (L / 2.0) / nsh
    edges = np.arange(nsh + 1) * w
    cen = 0.5 * (edges[1:] + edges[:-1])
    esh, dmean = [], []
    for k in range(nsh):
        m = (r >= edges[k]) & (r < edges[k + 1])
        esh.append(float(np.sum(dens_h3[m])))
        dmean.append(float(np.mean(dens_h3[m])) if m.any() else float("nan"))
    return cen, np.array(esh), np.array(dmean)


def loglog_slope(x, y):
    m = (x > 0) & (y > 0) & np.isfinite(y)
    if m.sum() < 3:
        return float("nan")
    return float(np.polyfit(np.log(x[m]), np.log(y[m]), 1)[0])


def tail_fit(dens, r, L, lo, hi, nb=None):
    """log-log slope of the shell-mean density over lo <= r <= hi (unit-width bins)."""
    nb = nb or max(4, int(round(hi - lo)))
    edges = np.linspace(lo, hi, nb + 1)
    xs, ys = [], []
    for k in range(nb):
        m = (r >= edges[k]) & (r < edges[k + 1])
        if m.sum() > 0:
            xs.append(0.5 * (edges[k] + edges[k + 1]))
            ys.append(float(np.mean(dens[m])))
    return loglog_slope(np.array(xs), np.array(ys))


def reads_on(M, cfg, mu, cP, tag=""):
    n, L, h = cfg["n"], cfg["L"], cfg["h"]
    h3 = h ** 3
    X, Y, Z = INS4.coords(n, h)
    r = np.sqrt(X * X + Y * Y + Z * Z)
    lam, P23, im_lam, im_P = spectral(M)
    # projector cross-check against the Lagrange route
    PL = p23_lagrange(M, lam)
    dP = float(np.max(np.abs(PL - P23)))
    # idempotency / rank
    idem = float(np.max(np.abs(P23 @ P23 - P23)))
    rank = float(np.mean(np.einsum("...kk->...", P23)))
    d_u = eu_density(M, cfg)
    d_v = v4dd_density(M)
    d_s = split_density(M, lam)
    d_k = kp23_density(M, cfg, P23)
    E_u = h3 * float(np.sum(d_u))
    V4 = h3 * float(np.sum(d_v))
    SPLIT = h3 * float(np.sum(d_s))
    KP = h3 * float(np.sum(d_k))
    E_stat = E_u + V4 + mu * SPLIT + cP * KP
    eu_cert = float(INS4.e_parts(M, cfg)[0])
    d_stat = h3 * (d_u + d_v + mu * d_s + cP * d_k)
    frac = {f"r<L/{k}": float(np.sum(d_stat[r < L / k]) / E_stat) for k in (4, 8)}
    cen, esh, dmean = shells(d_stat, r, L)
    outer = cen > L / 4
    slope_shell = loglog_slope(cen[outer], esh[outer])
    slope_shell_alt = loglog_slope(cen[cen > L / 3], esh[cen > L / 3])
    win = (6.0, L / 2 - 3.0)
    win_alt = (L / 6.0, L / 2 - 1.0)
    win_in = (3.0, L / 4.0)
    tails = {}
    for nm, dd in (("KP23", d_k), ("E_u", d_u), ("V4dd", d_v), ("SPLIT", d_s), ("E_stat", d_stat)):
        tails[nm] = {"producer_window": tail_fit(dd, r, L, *win),
                     "alt_window": tail_fit(dd, r, L, *win_alt),
                     "inner_window": tail_fit(dd, r, L, *win_in)}
    ext = (r > 0.35 * L) & (r < 0.45 * L)
    drift = float(np.max(np.abs(lam.real[ext] - VAC_SPEC)))
    drift_all_interior = float(np.max(np.abs(lam.real[(r > 0.25 * L) & (r < 0.5 * L)] - VAC_SPEC)))
    split = np.sqrt(d_s)
    ic = np.unravel_index(int(np.argmax(split)), split.shape)
    max_split = float(split.max())
    split_ext = float(split[ext].max())
    # inertias with the local clock (unnormalized, as a0_local returns it) and unit-normalized
    a0 = R13.a0_local(M)
    a0n = float(np.sqrt(np.sum(a0 * a0)))
    kin_I1 = float(INS4.kin_of(M, a0, cfg))
    kin_KP = kp_clock(M, cfg, P23, a0)
    a0u = a0 / max(a0n, 1e-300)
    kin_I1_unit = float(INS4.kin_of(M, a0u, cfg))
    kin_KP_unit = kp_clock(M, cfg, P23, a0u)
    kd = R13.kin_density(M, a0, cfg)
    kdens_ext = float(np.max(np.abs(kd[ext])))
    Om0 = P23 @ a0 @ ETA @ P23
    kp_ext = float(np.max(np.abs(0.5 * inner(Om0, Om0)[ext])))
    a0_ext = float(np.max(np.abs(a0[ext])))
    # the vacuum-shell identity: a0_local on an exactly degenerate cell is zero
    off = float(np.max(np.abs(M[..., 0, 1:])))
    row = {"tag": tag, "n": n, "L": L, "mu": mu, "cP": cP,
           "E_u": E_u, "E_u_certified_e_parts": eu_cert, "V4dd": V4, "SPLIT": SPLIT, "KP23": KP,
           "E_stat": E_stat, "frac": frac,
           "shell_centers": cen.tolist(), "shell_energy": esh.tolist(),
           "outer_slope_r>L/4": slope_shell, "outer_slope_r>L/3": slope_shell_alt,
           "tail_exponents": tails, "tail_windows": {"producer": win, "alt": win_alt, "inner": win_in},
           "drift_0.35-0.45L": drift, "drift_0.25-0.5L": drift_all_interior,
           "max_split": max_split, "max_split_cell_r": float(r[ic]), "split_exterior_max": split_ext,
           "kin_I1_a0local": kin_I1, "kin_KP23_a0local": kin_KP, "a0_local_norm": a0n,
           "kin_I1_unit_a0": kin_I1_unit, "kin_KP23_unit_a0": kin_KP_unit,
           "kin_density_ext_max": kdens_ext, "kp_clock_density_ext_max": kp_ext, "a0_ext_max": a0_ext,
           "P23_lagrange_vs_eig_maxdiff": dP, "P23_idempotency": idem, "P23_mean_trace": rank,
           "max_imag_lambda": im_lam, "max_imag_P23": im_P, "offblock_max": off}
    log(f"{tag}: E_u {E_u:.4f} (cert {eu_cert:.4f}) V4dd {V4:.4f} SPLIT {SPLIT:.3e} KP23 {KP:.4f} "
        f"E_stat {E_stat:.4f} frac {frac} slope {slope_shell:.3f}/{slope_shell_alt:.3f} "
        f"tails KP {tails['KP23']['producer_window']:.2f}/{tails['KP23']['alt_window']:.2f} "
        f"Eu {tails['E_u']['producer_window']:.2f}/{tails['E_u']['alt_window']:.2f} "
        f"drift {drift:.2e} maxsplit {max_split:.4f}@r={r[ic]:.1f} kinI1 {kin_I1:.3e} kinKP {kin_KP:.3e} "
        f"unit: {kin_I1_unit:.3e}/{kin_KP_unit:.3e} dP {dP:.1e} imlam {im_lam:.1e}")
    return row


# ============================================================ Hessian at the vacuum
def sym_basis():
    """Frobenius-orthonormal basis of the 10-dim symmetric 4x4 space."""
    B = []
    for k in range(4):
        E = np.zeros((4, 4)); E[k, k] = 1.0; B.append(E)
    for k in range(4):
        for l in range(k + 1, 4):
            E = np.zeros((4, 4)); E[k, l] = E[l, k] = 1.0 / np.sqrt(2.0); B.append(E)
    return B


def cell_energy(M, mu):
    """V4^dd + mu SPLIT on a single cell (M 4x4)."""
    v = float(v4dd_density(M[None]).real[0])
    if mu:
        lam = np.sort_complex(np.linalg.eigvals(M @ ETA))
        lam = lam[np.argsort(lam.real)]
        v += mu * float(np.abs(lam[2] - lam[1]) ** 2)
    return v


def hessian_fd(mu, t=1e-3):
    B = sym_basis()
    H = np.zeros((10, 10))
    e0 = cell_energy(D_DEG, mu)
    for i in range(10):
        for j in range(i, 10):
            if i == j:
                H[i, i] = (cell_energy(D_DEG + t * B[i], mu) - 2 * e0 + cell_energy(D_DEG - t * B[i], mu)) / t ** 2
            else:
                H[i, j] = H[j, i] = (cell_energy(D_DEG + t * (B[i] + B[j]), mu)
                                     - cell_energy(D_DEG + t * (B[i] - B[j]), mu)
                                     - cell_energy(D_DEG - t * (B[i] - B[j]), mu)
                                     + cell_energy(D_DEG - t * (B[i] + B[j]), mu)) / (4 * t ** 2)
    return H


def hessian_gauss_newton():
    """exact at the zero of the residual: H = 2 W1 sum_p J_p J_p^T,
    J_p = d tr(N^p)/dM = p sym(eta N^{p-1})^T, N = d eta diagonal."""
    B = sym_basis()
    lam = np.diag(D_DEG @ ETA)
    H = np.zeros((10, 10))
    for p in range(1, 5):
        Jm = p * INS4.sym4((ETA @ np.diag(lam ** (p - 1))).T)
        Jv = np.array([np.sum(Jm * Bi) for Bi in B])
        H += 2 * W1 * np.outer(Jv, Jv)
    return H


def hess_mode():
    names = ["E00", "E11", "E22", "E33", "b01", "b02", "b03", "r12", "r13", "s23"]
    out = {}
    Hgn = hessian_gauss_newton()
    for mu in (0.0, 0.01):
        rows = {}
        for t in (1e-3, 1e-4):
            H = hessian_fd(mu, t)
            w, V = np.linalg.eigh(H)
            thr = 1e-6
            null = int(np.sum(np.abs(w) < thr))
            rows[f"t{t:g}"] = {"eigs": w.tolist(), "null_count": null, "threshold": thr,
                               "nonzero": [float(x) for x in w if abs(x) >= thr]}
            if t == 1e-3:
                # split stiffness along the two traceless (2,3) directions
                d22m33 = (np.zeros((4, 4))); d22m33[2, 2], d22m33[3, 3] = 1 / np.sqrt(2), -1 / np.sqrt(2)
                s23 = np.zeros((4, 4)); s23[2, 3] = s23[3, 2] = 1 / np.sqrt(2)
                stiff = {}
                for nm, Dm in (("E22-E33", d22m33), ("sym23", s23)):
                    stiff[nm] = float((cell_energy(D_DEG + t * Dm, mu) - 2 * cell_energy(D_DEG, mu)
                                       + cell_energy(D_DEG - t * Dm, mu)) / t ** 2)
                rows["split_stiffness"] = stiff
                rows["gauss_newton_maxdiff"] = float(np.max(np.abs(H - Hgn))) if mu == 0 else None
                # which basis directions are null: project the null eigenvectors on the basis
                nullvec = V[:, np.abs(w) < thr]
                rows["null_span_weights"] = {names[k]: float(np.sum(nullvec[k] ** 2)) for k in range(10)}
        out[f"mu{mu:g}"] = rows
        log(f"hess mu {mu}: null {rows['t0.001']['null_count']} nonzero {rows['t0.001']['nonzero']} "
            f"stiff {rows['split_stiffness']}")
    wg = np.linalg.eigvalsh(Hgn)
    out["gauss_newton_eigs"] = wg.tolist()
    log(f"gauss-newton eigs: {wg}")
    return out


# ============================================================ mutation: own gradient + FIRE
def d1_adj(g, ax, h, st):
    """adjoint of INS4.d1 (verified by the inner-product test below)."""
    out = np.zeros_like(g)
    sl = [slice(None)] * g.ndim

    def at(i):
        s = list(sl); s[ax] = i; return tuple(s)
    if st == "fwd":
        out[at(slice(1, None))] += g[at(slice(0, -1))] / h
        out[at(slice(0, -1))] -= g[at(slice(0, -1))] / h
    elif st == "bwd":
        out[at(slice(1, None))] += g[at(slice(1, None))] / h
        out[at(slice(0, -1))] -= g[at(slice(1, None))] / h
    else:
        raise ValueError(st)
    return out


def energy_uv(M, cfg):
    h3 = cfg["h"] ** 3
    return h3 * np.sum(eu_density(M, cfg)) + h3 * np.sum(v4dd_density(M))


def grad_uv(M, cfg):
    """own gradient of E_u + V4^dd wrt symmetric M.
    <F,F>_eta with F = Ai eta Aj - Aj eta Ai:  d<F,F>/dAi = 2 [eta F eta (eta Aj)^T ... ] derived as
    d<F,F> = 2 <F, dF>_eta, <F, X eta Y>_eta = tr(eta F eta Y^T eta X^T) -> d/dX = eta Y eta F^T eta."""
    h3 = cfg["h"] ** 3
    h = cfg["h"]
    Gt = np.zeros_like(M)
    for br, wt in INS4.branches(cfg["stencil"]):
        A = [INS4.d1(M, ax, h, br) for ax in range(3)]
        dA = [np.zeros_like(M) for _ in range(3)]
        for i in range(3):
            for j in range(i + 1, 3):
                F = A[i] @ ETA @ A[j] - A[j] @ ETA @ A[i]
                # <F, Ai eta Aj>: d/dAi = eta Aj eta F^T eta ; d/dAj = (eta Ai eta)^T ... use generic:
                # <F, X eta Y>_eta : d/dX = eta Y eta F^T eta ; d/dY = eta F^T eta X eta ... careful:
                # tr(eta F eta Y^T eta X^T) -> d/dY = (eta X^T eta F eta)^T = eta F^T eta X eta
                # F is antisymmetric: F^T = -F.
                FT = F.swapaxes(-1, -2)
                dX = ETA @ A[j] @ ETA @ FT @ ETA           # d<F, Ai eta Aj>/dAi
                dY = ETA @ FT @ ETA @ A[i] @ ETA            # d<F, Ai eta Aj>/dAj
                dX2 = ETA @ A[i] @ ETA @ FT @ ETA           # d<F, Aj eta Ai>/dAj
                dY2 = ETA @ FT @ ETA @ A[j] @ ETA           # d<F, Aj eta Ai>/dAi
                dA[i] += 8.0 * (dX - dY2)
                dA[j] += 8.0 * (dY - dX2)
        for ax in range(3):
            Gt += wt * d1_adj(dA[ax], ax, h, br)
    N = M @ ETA
    pows = [np.broadcast_to(np.eye(4), M.shape).copy()]
    for _ in range(3):
        pows.append(pows[-1] @ N)
    t = [np.einsum("...kk->...", P @ N) for P in pows]
    GV = np.zeros_like(M)
    for p in range(1, 5):
        coef = 2.0 * W1 * (t[p - 1] - C_DD[p - 1]) * p
        GV += coef[..., None, None] * (ETA @ pows[p - 1]).swapaxes(-1, -2)
    return h3 * INS4.sym4(Gt) + h3 * INS4.sym4(GV)


def gates():
    rng = np.random.default_rng(7)
    out = {}
    # adjoint test
    f = rng.normal(size=(9, 9, 9, 4, 4)); g = rng.normal(size=f.shape)
    errs = []
    for st in ("fwd", "bwd"):
        for ax in range(3):
            a = np.sum(INS4.d1(f, ax, 0.7, st) * g); b = np.sum(f * d1_adj(g, ax, 0.7, st))
            errs.append(abs(a - b) / abs(a))
    out["adjoint_relerr"] = float(max(errs))
    # complex-step gradient gate (own energy is complex-safe: no abs, no eig)
    cfg = cfg_of(9, 13.5)
    M = INS4.sym4(rng.normal(size=(9, 9, 9, 4, 4)))
    Gm = grad_uv(M, cfg)
    errs = []
    for _ in range(4):
        V = INS4.sym4(rng.normal(size=M.shape))
        de_an = float(np.sum(Gm * V))
        de_cs = energy_uv(M + 1e-30j * V, cfg).imag / 1e-30
        errs.append(abs(de_cs - de_an) / abs(de_cs))
    out["complex_step_relerr"] = float(max(errs))
    # E_u agreement with the certified stack on the same random field
    out["E_u_vs_certified_relerr"] = float(abs(cfg["h"] ** 3 * np.sum(eu_density(M, cfg)) - INS4.e_parts(M, cfg)[0])
                                           / INS4.e_parts(M, cfg)[0])
    # V4dd zero on the degenerate vacuum, positive on the certified vacuum
    Mv = np.broadcast_to(D_DEG, (3, 3, 3, 4, 4)).copy()
    out["V4dd_on_degenerate_vacuum"] = float(np.max(v4dd_density(Mv)))
    out["V4dd_on_certified_vacuum"] = float(v4dd_density(np.diag([8.0, 1.0, 0.3, 0.0])[None])[0])
    # a0_local vanishes on the degenerate vacuum (rotated)
    R = B8.rot_field(B8.G3, np.array([[[0.7]]])) @ B8.rot_field(B8.G2, np.array([[[0.4]]]))
    Mr = R @ D_DEG @ R.swapaxes(-1, -2)
    out["a0_local_on_degenerate_vacuum"] = float(np.max(np.abs(R13.a0_local(Mr))))
    out["pass"] = bool(out["adjoint_relerr"] < 1e-12 and out["complex_step_relerr"] < 1e-8
                       and out["E_u_vs_certified_relerr"] < 1e-12 and out["V4dd_on_degenerate_vacuum"] < 1e-20)
    log(f"gates: {out}")
    return out


def fire_own(M0, cfg, free, steps, dt0=0.01, dt_max=0.1, log_every=50):
    M = M0.copy(); v = np.zeros_like(M)
    fr = free[..., None, None].astype(float)
    dt, alpha, nup = dt0, 0.1, 0
    F = -grad_uv(M, cfg) * fr
    hist = [{"it": 0, "E": float(energy_uv(M, cfg)), "fmax": float(np.max(np.abs(F)))}]
    for it in range(1, steps + 1):
        P = float(np.sum(F * v))
        if P > 0:
            nup += 1
            vn, fn = np.sqrt(np.sum(v * v)), np.sqrt(np.sum(F * F))
            v = (1 - alpha) * v + alpha * F / max(fn, 1e-300) * vn
            if nup > 5:
                dt = min(dt * 1.1, dt_max); alpha *= 0.99
        else:
            v[:] = 0; dt *= 0.5; alpha = 0.1; nup = 0
        v += dt * F
        M += dt * v
        F = -grad_uv(M, cfg) * fr
        if it % log_every == 0 or it == steps:
            E = float(energy_uv(M, cfg)); fm = float(np.max(np.abs(F)))
            hist.append({"it": it, "E": E, "fmax": fm, "dt": dt})
            log(f"  fire it {it} E {E:.6f} fmax {fm:.3e} dt {dt:.4f}")
    return M, hist


def mutate_mode():
    """(a) sampled gradient components of the FULL E_stat (FD on own energy) at the producer's
    n32 mu0 cP1 field; (b) a 200-step own FIRE on E_u + V4^dd from the n32 mu0 cP0 field;
    (c) 200 own FIRE steps from the degenerate-vacuum seed on n16 L24 (the small-grid replica)."""
    out = {"gates": gates()}
    # (a) stationarity samples: FD directional derivative of E_stat = E_u + V4dd + cP KP23 (mu = 0)
    n, L = 32, 48.0
    cfg = cfg_of(n, L)
    M = np.load(os.path.join(CK, f"relax_n{n}_L{L:g}_mu0_cP1.npy"))
    h3 = cfg["h"] ** 3

    def e_stat(Mx, cP=1.0):
        lam, P23, _, _ = spectral(Mx)
        return h3 * (np.sum(eu_density(Mx, cfg)) + np.sum(v4dd_density(Mx)) + cP * np.sum(kp23_density(Mx, cfg, P23)))
    rng = np.random.default_rng(11)
    B = sym_basis()
    samples = []
    c = n // 2
    sites = [(c, c, c), (c + 1, c, c), (c + 2, c + 1, c), (c + 4, c, c), (c + 6, c + 2, c + 1),
             (c + 8, c, c), (c + 10, c + 3, c + 2), (c + 12, c, c), (c + 3, c + 3, c + 3), (c - 5, c + 1, c - 2)]
    t = 1e-4
    e0 = e_stat(M)
    for (i, j, k) in sites:
        for b in (0, 1, 2, 4, 8, 9):                 # E00, E11, E22, b01, r13, s23
            V = np.zeros_like(M); V[i, j, k] = B[b]
            g = (e_stat(M + t * V) - e_stat(M - t * V)) / (2 * t)
            samples.append({"site": [i, j, k], "basis": b, "grad": float(g)})
    gmax = max(abs(s["grad"]) for s in samples)
    out["stationarity_n32_mu0_cP1"] = {"E_stat_own": float(e0), "samples": samples, "max_abs_sampled_grad": float(gmax),
                                       "note": "gradient per h^3-weighted energy wrt one Frobenius-unit cell direction"}
    log(f"stationarity n32 mu0 cP1: E_stat {e0:.4f} max sampled |grad| {gmax:.3e}")
    # (b) own 200-step FIRE from the producer's mu0 cP0 field (E_u + V4dd only)
    M0 = np.load(os.path.join(CK, f"relax_n{n}_L{L:g}_mu0_cP0.npy"))
    free = ~INS4.pin_shell(n, cfg["h"])
    r0 = reads_on(M0, cfg, 0.0, 0.0, "pre-mutation n32 mu0 cP0")
    M1, hist = fire_own(M0, cfg, free, 200)
    r1 = reads_on(M1, cfg, 0.0, 0.0, "post-200-FIRE n32 mu0 cP0")
    out["descent_n32_mu0_cP0"] = {"trace": hist, "E_before": hist[0]["E"], "E_after": hist[-1]["E"],
                                  "reads_before": r0, "reads_after": r1,
                                  "field_rel_move": float(np.sqrt(np.sum((M1 - M0) ** 2)) / np.sqrt(np.sum((M0 - D_DEG) ** 2)))}
    # (c) small-grid replica: 200 own FIRE steps from the degenerate seed, n16 L24
    cfg16 = cfg_of(16, 24.0)
    Ms = seed_degenerate(cfg16)
    free16 = ~INS4.pin_shell(16, cfg16["h"])
    rs0 = reads_on(Ms, cfg16, 0.0, 0.0, "seed n16 L24")
    Ms1, hist16 = fire_own(Ms, cfg16, free16, 200, log_every=50)
    rs1 = reads_on(Ms1, cfg16, 0.0, 0.0, "seed+200FIRE n16 L24")
    out["seed_descent_n16_L24"] = {"trace": hist16, "reads_seed": rs0, "reads_after": rs1}
    return out


def seed_degenerate(cfg):
    """B8.dressed(cfg, 0) with the DEGENERATE vacuum: Q = Qh (no boost), M = Qh d Qh^T."""
    n, h = cfg["n"], cfg["h"]
    X, Y, Z = INS4.coords(n, h)
    rho = np.sqrt(X * X + Y * Y)
    phi = np.arctan2(Y, X)
    th = -np.arctan2(Z, rho)
    Qh = np.einsum("...ab,...bc->...ac", B8.rot_field(B8.G3, phi), B8.rot_field(B8.G2, th))
    M = np.einsum("...ab,bc,...dc->...ad", Qh, D_DEG, Qh)
    return INS4.sym4(M)


# ============================================================ certified reference
def ref_mode():
    M, cfg, rec = R13.seed_hedgehog(32, 48.0)
    n, L, h = cfg["n"], cfg["L"], cfg["h"]
    h3 = h ** 3
    X, Y, Z = INS4.coords(n, h)
    r = np.sqrt(X * X + Y * Y + Z * Z)
    d_u = eu_density(M, cfg)
    # certified V4 density (spectrum (-g, 1, delta, 0)): own trace form with the certified targets
    t = traces(M)
    cc = tuple((-G) ** p + 1.0 + DELTA ** p for p in range(1, 5))
    d_v = W1 * sum((t[p] - cc[p]) ** 2 for p in range(4))
    E_u, V4 = h3 * float(np.sum(d_u)), h3 * float(np.sum(d_v))
    eu_c, ev_c = INS4.e_parts(M, cfg)
    d_stat = h3 * (d_u + d_v)
    E = E_u + V4
    cen, esh, dm = shells(d_stat, r, L)
    outer = cen > L / 4
    out = {"source": rec, "E_u": E_u, "V4": V4, "E_u_certified": float(eu_c), "V4_certified": float(ev_c),
           "frac": {f"r<L/{k}": float(np.sum(d_stat[r < L / k]) / E) for k in (4, 8)},
           "shell_centers": cen.tolist(), "shell_energy": esh.tolist(),
           "outer_slope_r>L/4": loglog_slope(cen[outer], esh[outer]),
           "outer_slope_r>L/3": loglog_slope(cen[cen > L / 3], esh[cen > L / 3]),
           "tail_E_u": {"producer_window": tail_fit(d_u, r, L, 6.0, L / 2 - 3), "alt_window": tail_fit(d_u, r, L, L / 6, L / 2 - 1)},
           "tail_E_stat": {"producer_window": tail_fit(d_u + d_v, r, L, 6.0, L / 2 - 3)}}
    log(f"certified ref: E_u {E_u:.4f} V4 {V4:.4f} frac {out['frac']} slope {out['outer_slope_r>L/4']:.3f}/"
        f"{out['outer_slope_r>L/3']:.3f} tail E_u {out['tail_E_u']}")
    return out


# ============================================================ driver
def reads_mode():
    rows = {}
    for n, L in ((32, 48.0), (48, 72.0)):
        cfg = cfg_of(n, L)
        for mu in (0.0, 0.01):
            for cP in (0, 1):
                key = f"n{n}_L{L:g}_mu{mu:g}_cP{cP}"
                p = os.path.join(CK, f"relax_{key}.npy")
                if not os.path.exists(p):
                    log(f"missing {p}"); continue
                rows[key] = reads_on(np.load(p), cfg, mu, float(cP), key)
    # L-exponents between the boxes, p = ln[E(72)/E(48)] / ln 1.5
    expo = {}
    for mu in (0.0, 0.01):
        for cP in (0, 1):
            a, b = rows.get(f"n32_L48_mu{mu:g}_cP{cP}"), rows.get(f"n48_L72_mu{mu:g}_cP{cP}")
            if a and b:
                e = {}
                for q in ("E_u", "KP23", "V4dd", "E_stat", "SPLIT"):
                    if a[q] > 0 and b[q] > 0:
                        e[q] = float(np.log(b[q] / a[q]) / np.log(1.5))
                # a 1/L-tail extrapolation: E(L) = E_inf - c/L
                for q in ("E_u", "KP23"):
                    Einf = (72 * b[q] - 48 * a[q]) / (72 - 48)
                    e[f"{q}_Einf_1/L"] = float(Einf)
                    e[f"{q}_tail_frac_L48"] = float(1 - a[q] / Einf)
                expo[f"mu{mu:g}_cP{cP}"] = e
    # mu-sensitivity: field-level and read-level
    musens = {}
    for n, L in ((32, 48.0), (48, 72.0)):
        for cP in (0, 1):
            a, b = rows.get(f"n{n}_L{L:g}_mu0_cP{cP}"), rows.get(f"n{n}_L{L:g}_mu0.01_cP{cP}")
            if a and b:
                Ma = np.load(os.path.join(CK, f"relax_n{n}_L{L:g}_mu0_cP{cP}.npy"))
                Mb = np.load(os.path.join(CK, f"relax_n{n}_L{L:g}_mu0.01_cP{cP}.npy"))
                musens[f"n{n}_cP{cP}"] = {"field_maxdiff": float(np.max(np.abs(Ma - Mb))),
                                         "reads_absdiff": {q: abs(a[q] - b[q]) for q in ("E_u", "KP23", "V4dd", "max_split", "drift_0.35-0.45L", "outer_slope_r>L/4")}}
    return {"fields": rows, "L_exponents": expo, "mu_sensitivity": musens}


# ============================================================ embedding / seed-texture test
def embed_mode():
    """are the two boxes independent?  (a) n48 center 32^3 block vs the n32 field; (b) each field
    vs the analytic seed texture outside r > L/8 and r > L/4; (c) the seed texture's own energies,
    tail exponents and L-exponent (h = 1.5 in both boxes)."""
    out = {}
    for mu in (0.0, 0.01):
        for cP in (0, 1):
            M32 = np.load(os.path.join(CK, f"relax_n32_L48_mu{mu:g}_cP{cP}.npy"))
            M48 = np.load(os.path.join(CK, f"relax_n48_L72_mu{mu:g}_cP{cP}.npy"))
            blk = M48[8:40, 8:40, 8:40]
            key = f"mu{mu:g}_cP{cP}"
            out[key] = {"n48_center_block_vs_n32_maxdiff": float(np.max(np.abs(blk - M32))),
                        "n48_center_block_vs_n32_maxdiff_r<L/8": None}
            rows = {}
            for n, L, M in ((32, 48.0, M32), (48, 72.0, M48)):
                cfg = cfg_of(n, L)
                S = seed_degenerate(cfg)
                X, Y, Z = INS4.coords(n, cfg["h"])
                r = np.sqrt(X * X + Y * Y + Z * Z)
                dd = np.max(np.abs(M - S), axis=(-1, -2))
                rows[f"n{n}"] = {"maxdiff_vs_seed_all": float(dd.max()),
                                 "maxdiff_vs_seed_r>L/8": float(dd[r > L / 8].max()),
                                 "maxdiff_vs_seed_r>L/4": float(dd[r > L / 4].max()),
                                 "maxdiff_vs_seed_r>0.35L": float(dd[r > 0.35 * L].max()),
                                 "r_where_diff>1e-3": float(r[dd > 1e-3].max()) if (dd > 1e-3).any() else 0.0,
                                 "r_where_diff>1e-6": float(r[dd > 1e-6].max()) if (dd > 1e-6).any() else 0.0}
            if mu == 0.0 and cP == 0:
                X, Y, Z = INS4.coords(32, 1.5)
                r = np.sqrt(X * X + Y * Y + Z * Z)
                dd = np.max(np.abs(blk - M32), axis=(-1, -2))
                out[key]["n48_center_block_vs_n32_maxdiff_r<L/8"] = float(dd[r < 6].max())
            out[key]["vs_seed"] = rows
            log(f"embed {key}: {out[key]}")
    # the seed texture itself
    seed = {}
    for n, L in ((32, 48.0), (48, 72.0)):
        cfg = cfg_of(n, L)
        S = seed_degenerate(cfg)
        seed[f"n{n}"] = reads_on(S, cfg, 0.0, 1.0, f"SEED texture n{n} L{L:g}")
    seed["L_exponent_E_u"] = float(np.log(seed["n48"]["E_u"] / seed["n32"]["E_u"]) / np.log(1.5))
    seed["L_exponent_KP23"] = float(np.log(seed["n48"]["KP23"] / seed["n32"]["KP23"]) / np.log(1.5))
    out["seed_texture"] = seed
    log(f"seed L-exponents: E_u {seed['L_exponent_E_u']:.4f} KP23 {seed['L_exponent_KP23']:.4f}")
    return out


if __name__ == "__main__":
    mode = ARGV[0] if ARGV else "all"
    res = {}
    if os.path.exists(OUT):
        res = json.load(open(OUT))
    if mode in ("all", "hess"):
        res["hessian"] = hess_mode()
    if mode in ("all", "ref"):
        res["certified_reference"] = ref_mode()
    if mode in ("all", "reads"):
        res.update(reads_mode())
    if mode in ("all", "mutate"):
        res["mutation"] = mutate_mode()
    if mode in ("all", "embed"):
        res["embedding"] = embed_mode()
    res["wall_s"] = round(time.time() - T0, 1)
    json.dump(res, open(OUT, "w"), indent=1)
    log(f"wrote {OUT}")
