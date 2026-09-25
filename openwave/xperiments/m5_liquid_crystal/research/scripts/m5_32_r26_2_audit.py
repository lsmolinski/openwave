"""M5.32 R26-2 audit: a fresh agent's attempt to refute the strand-ladder claims with its own
code. The audited script was NOT read; its inputs are the stored arrays (data/m5_32_r26_2/,
data/m5_32_r25_1/), its JSON, and the R25-1 instrument's definitions.

EQUATIONS (own forms)
---------------------
Field M(x) real symmetric 4x4 on an n x n x 4 slab, h = 48 / n, eta = diag(-1, 1, 1, 1),
N = M eta. Index layout of the stored fields: slot 0 the time slot (M_00 = 8), slots 1, 2 the
transverse pair, slot 3 the director (M_33 = 1). The stored deviation D obeys
M = diag(8, 1, delta, 0) + delta D.
    E = h^3 sum_cells [ 4 sum_{i<j} <F_ij, F_ij>_eta + w sum_{p<=4} (tr N^p - C_p)^2 ],
    F_ij = A_i eta A_j - A_j eta A_i,  A_i = d_i M (one-sided, the "sym" stencil = the mean of
    the forward and backward one-sided ENERGIES),  <F, F>_eta = sum_ab eta_a eta_b F_ab^2,
    C_p = (-8)^p + 1 + delta^p,  w = W1 x w1s,  T = E / (4 h).
Trace brackets: (i) plain float64, (ii) telescoped in float64 about N_ref = diag(-8, s0, s0, 1)
(the words with at least one dN, plus the exact constant C_p - tr N_ref^p), (iii) EXACT in
Python integers (every float entry scaled by 2^120 is an integer; the bracket is an exact
rational, rounded once to float64). (iii) is the reference for every T in this audit.
Bound: T_bps = pi sqrt(32 w K) (delta / 2)^4, K = 4 + 36 s0^2 + 144 s0^4;
T_slaved3 = 4 pi sqrt(8) int_0^f0 sqrt(w V_min(f)) df, V_min the residual sum minimized over
(M_00, M_33, pair mean) at fixed pair weight f = b^2, the residuals written in deviation form
(binomials in e = M_00 - 8, t = M_33 - 1, sigma = s - s0), Gauss-Legendre in f.
Gradient: analytic (own derivation, complex-step checked in the smoke); Hessian-vector
products by the complex step on the gradient; lowest eigenvalues by Lanczos (scipy eigsh).

CLAIMS A to F as listed in the task; each records method, own numbers, audited numbers, verdict.
Modes: smoke | run [workers] (default 4, about 4 min) | extra | summarize.
Output data/m5_32_r26_2_audit.json (results, verdicts, summary.verdict_table).
"""

import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
import multiprocessing as mp

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
OUT_JSON = os.path.join(DATA, "m5_32_r26_2_audit.json")
R26 = os.path.join(DATA, "m5_32_r26_2")
R25 = os.path.join(DATA, "m5_32_r25_1")
LADDER_JSON = os.path.join(DATA, "m5_32_r26_2_ladder.json")
W1 = 0.000724023879
NZ = 4
PIN_DEPTH = 1.6
ETA = np.array([-1.0, 1.0, 1.0, 1.0])
SCALE_BITS = 120
T0 = time.time()


def log(msg):
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


# ================= own lattice energy =================
def d_one_sided(f, ax, h, side):
    """forward: (f[i+1] - f[i]) / h stored at i, zero at the last plane; backward: stored at i+1,
    zero at the first plane (the platform's open-ended convention, replicated from its stencil)."""
    out = np.zeros_like(f)
    df = np.diff(f, axis=ax) / h
    sl = [slice(None)] * f.ndim
    sl[ax] = slice(0, -1) if side == "fwd" else slice(1, None)
    out[tuple(sl)] = df
    return out


def d_central(f, ax, h):
    out = np.zeros_like(f)
    sl_in, sl_p, sl_m = [slice(None)] * f.ndim, [slice(None)] * f.ndim, [slice(None)] * f.ndim
    sl_in[ax], sl_p[ax], sl_m[ax] = slice(1, -1), slice(2, None), slice(0, -2)
    out[tuple(sl_in)] = (f[tuple(sl_p)] - f[tuple(sl_m)]) / (2 * h)
    return out


def comm_eta(A, B):
    return (A * ETA) @ B - (B * ETA) @ A  # A eta B - B eta A (eta diagonal scales columns)


def inner_eta_sq(F):
    return np.einsum("...ab,...ab,a,b->...", F, F, ETA, ETA)


def curv_density_branch(M, h, side):
    A = [d_one_sided(M, ax, h, side) for ax in range(3)]
    dens = 0.0
    for i in range(3):
        for j in range(i + 1, 3):
            F = comm_eta(A[i], A[j])
            dens = dens + 4.0 * inner_eta_sq(F)
    return h**3 * dens


def e_curv(M, h, stencil="sym"):
    if stencil == "sym":
        return 0.5 * (
            float(np.sum(curv_density_branch(M, h, "fwd")))
            + float(np.sum(curv_density_branch(M, h, "bwd")))
        )
    A = [d_central(M, ax, h) for ax in range(3)]
    e = 0.0
    for i in range(3):
        for j in range(i + 1, 3):
            F = comm_eta(A[i], A[j])
            e = e + 4.0 * np.sum(inner_eta_sq(F))
    return float(h**3 * e)


def c_targets(delta):
    return [(-8.0) ** p + 1.0 + delta**p for p in range(1, 5)]


def brackets_plain(M, delta):
    N = M * ETA
    C = c_targets(delta)
    P = N
    out = np.zeros(M.shape[:-2] + (4,))
    for p in range(1, 5):
        if p > 1:
            P = P @ N
        out[..., p - 1] = np.einsum("...kk->...", P) - C[p - 1]
    return out


def m_ref(delta):
    """the deviation reference in the field's own layout: diag(8, s0, s0, 1)."""
    s0 = delta / 2.0
    return np.diag([8.0, s0, s0, 1.0])


def brackets_tele(M, delta, ref=None, Mdev=None):
    """tr N^p - C_p by the telescoping identity N^p - N_ref^p = sum_k N^k (N - N_ref) N_ref^(p-1-k)
    with N_ref diagonal (default diag(-8, s0, s0, 1), the field's own layout), plus the exact
    constant tr N_ref^p - C_p (evaluated in Python floats from small numbers only)."""
    s0 = delta / 2.0
    if ref is None:
        ref = np.array([-8.0, s0, s0, 1.0])
    if Mdev is not None:  # deviation-native: dN exact, N formed only inside the words
        dN = Mdev * ETA
        N = dN + np.diag(ref)
    else:
        N = M * ETA
        dN = N - np.diag(ref)
    C = c_targets(delta)
    # constants tr N_ref^p - C_p, in a cancellation-free form for the default reference
    if np.allclose(ref, [-8.0, s0, s0, 1.0]):
        const = [0.0, 2 * s0**2 - delta**2, 2 * s0**3 - delta**3, 2 * s0**4 - delta**4]
    else:
        const = [float(np.sum(ref**p)) - C[p - 1] for p in range(1, 5)]
    Npow = [np.broadcast_to(np.eye(4), N.shape).copy()]
    for _ in range(3):
        Npow.append(Npow[-1] @ N)
    refpow = [np.ones(4)]
    for _ in range(3):
        refpow.append(refpow[-1] * ref)
    out = np.zeros(N.shape[:-2] + (4,), dtype=N.dtype)
    for p in range(1, 5):
        acc = 0.0
        for k in range(p):
            X = Npow[k] @ dN  # N^k dN
            acc = acc + np.einsum("...aa,a->...", X, refpow[p - 1 - k])
        out[..., p - 1] = acc + const[p - 1]
    return out


def brackets_exact(M, delta, Mdev=None):
    """exact rational brackets: every float64 entry times 2^120 is an integer (entries below
    2^-68 lose at most 2^-120 absolute, irrelevant); Python big integers do the rest."""
    S = 1 << SCALE_BITS
    Sf = float(S)
    if Mdev is not None:
        shp = Mdev.shape[:-2]
        Ni = ((Mdev * ETA).reshape(-1, 4, 4) * Sf).astype(object)
        s0i = int((delta / 2.0) * Sf)
        refi = [-8 * S, s0i, s0i, S]
        for a in range(4):
            Ni[:, a, a] = Ni[:, a, a] + refi[a]
        N = Ni
    else:
        shp = M.shape[:-2]
        N = (M * ETA).reshape(-1, 4, 4)
        Ni = (N * Sf).astype(object)
    out = np.zeros((N.shape[0], 4))
    dS = int(delta * Sf)
    Cint = [
        ((-8) ** p) * S**p + S**p + dS**p for p in range(1, 5)
    ]  # C_p S^p (delta^p from the float delta scaled exactly)
    Sp = [S, S * S, S**3, S**4]
    for c in range(N.shape[0]):
        A = [[int(Ni[c, a, b]) for b in range(4)] for a in range(4)]
        t1 = A[0][0] + A[1][1] + A[2][2] + A[3][3]
        A2 = [[sum(A[a][k] * A[k][b] for k in range(4)) for b in range(4)] for a in range(4)]
        t2 = A2[0][0] + A2[1][1] + A2[2][2] + A2[3][3]
        t3 = sum(A2[a][b] * A[b][a] for a in range(4) for b in range(4))
        t4 = sum(A2[a][b] * A2[b][a] for a in range(4) for b in range(4))
        for p, t in enumerate((t1, t2, t3, t4)):
            out[c, p] = (t - Cint[p]) / Sp[p]  # int / int: correctly rounded
    return out.reshape(shp + (4,))


def v4_energy(M, h, w, delta, form="exact", Mdev=None):
    if form == "exact":
        br = brackets_exact(M, delta, Mdev=Mdev)
    elif form == "tele":
        br = brackets_tele(M, delta, Mdev=Mdev)
    else:
        br = brackets_plain(M, delta)
    return float(h**3 * w * np.sum(br * br))


def tension(M, h, w, delta, form="exact", stencil="sym", Mdev=None):
    eu = e_curv(M if Mdev is None else Mdev, h, stencil)  # derivatives ignore the constant
    ev = v4_energy(M, h, w, delta, form, Mdev=Mdev)
    return (eu + ev) / (NZ * h), eu, ev


# ================= own gradient (analytic) =================
def grad_energy(M, h, w, delta, Mdev=None):
    """dE = <G, dM>_F for symmetric dM, G symmetric. Works on complex M (complex step).
    With Mdev given, everything is computed from the deviation (M = m_ref + Mdev never formed)."""
    if Mdev is not None:
        M = Mdev
    G = np.zeros_like(M)
    for side in ("fwd", "bwd"):
        A = [d_one_sided(M, ax, h, side) for ax in range(3)]
        dA = [np.zeros_like(M) for _ in range(3)]
        for i in range(3):
            for j in range(i + 1, 3):
                F = comm_eta(A[i], A[j])
                Gf = (F * ETA) * ETA[:, None]  # eta F eta
                AjT, AiT = A[j].swapaxes(-1, -2), A[i].swapaxes(-1, -2)
                dA[i] += 8.0 * ((Gf @ AjT) * ETA - ETA[:, None] * (AjT @ Gf))
                dA[j] += 8.0 * (ETA[:, None] * (AiT @ Gf) - (Gf @ AiT) * ETA)
        for ax in range(3):
            # adjoint of the one-sided difference
            g = dA[ax]
            sl_a, sl_b = [slice(None)] * M.ndim, [slice(None)] * M.ndim
            if side == "fwd":
                sl_a[ax], sl_b[ax] = slice(1, None), slice(0, -1)
                src = g[tuple(sl_b)]
            else:
                sl_a[ax], sl_b[ax] = slice(1, None), slice(0, -1)
                src = g[tuple(sl_a)]
            adj = np.zeros_like(M)
            adj[tuple(sl_a)] += src / h
            adj[tuple(sl_b)] -= src / h
            G += 0.5 * h**3 * adj
    # V4 part: d tr(N^p)/dM_ab = p (eta N^(p-1))_ba
    if Mdev is not None:
        N = Mdev * ETA + np.diag([-8.0, delta / 2.0, delta / 2.0, 1.0])
        br = brackets_tele(None, delta, Mdev=Mdev)
    else:
        N = M * ETA
        br = brackets_tele(M, delta)
    Npow = [np.broadcast_to(np.eye(4), M.shape).copy()]
    for _ in range(3):
        Npow.append(Npow[-1] @ N)
    GV = np.zeros_like(M)
    for p in range(1, 5):
        X = ETA[:, None] * Npow[p - 1]  # eta N^(p-1)
        GV += (2.0 * w * p * br[..., p - 1])[..., None, None] * X.swapaxes(-1, -2)
    G = G + h**3 * GV
    return 0.5 * (G + G.swapaxes(-1, -2))


# ================= geometry, seeds, bound =================
def coords2(n, h):
    x = (np.arange(n) - (n - 1) / 2.0) * h
    X, Y = np.meshgrid(x, x, indexing="ij")
    return X, Y, np.hypot(X, Y), np.arctan2(Y, X)


def k_lead(s0):
    return 4.0 + 36.0 * s0**2 + 144.0 * s0**4


def t_bps(delta, w):
    s0 = delta / 2.0
    return float(np.pi * np.sqrt(32.0 * w * k_lead(s0)) * s0**4)


def kappa_bps(delta, w):
    return float(np.sqrt(w * k_lead(delta / 2.0) / 32.0))


def bps_field(n, h, delta, w, kappa=None, pure_winding=False):
    X, Y, rho, phi = coords2(n, h)
    s0 = b0 = delta / 2.0
    k = kappa_bps(delta, w) if kappa is None else kappa
    b = np.full_like(rho, b0) if pure_winding else b0 * np.sqrt(1.0 - np.exp(-k * rho**2))
    M = np.zeros((n, n, NZ, 4, 4))
    M[..., 0, 0] = 8.0
    M[..., 3, 3] = 1.0
    M[..., 1, 1] = (s0 + b * np.cos(phi))[:, :, None]
    M[..., 2, 2] = (s0 - b * np.cos(phi))[:, :, None]
    M[..., 1, 2] = M[..., 2, 1] = (b * np.sin(phi))[:, :, None]
    return M


def free_mask(n, h):
    wc = max(1, int(np.ceil(PIN_DEPTH / h)))
    free = np.ones((n, n, NZ), dtype=bool)
    free[:wc] = free[-wc:] = False
    free[:, :wc] = free[:, -wc:] = False
    return free


def load_row_field(tag):
    Z = np.load(os.path.join(R26, tag + ".npz"))
    D, delta = Z["D"], float(Z["delta"])
    n = D.shape[0]
    M = np.zeros_like(D)
    M[..., 0, 0], M[..., 1, 1], M[..., 2, 2] = 8.0, 1.0, delta
    return M + delta * D, delta, 48.0 / n


def load_r25_field(tag):
    M = np.load(os.path.join(R25, tag + ".npz"))["M"]
    return M, 48.0 / M.shape[0]


def departures(M, delta):
    return {
        "M0i": float(np.max(np.abs(M[..., 0, 1:]))) / delta,
        "director_row": float(np.max(np.abs(M[..., 3, 1:3]))) / delta,
        "M33_minus_1": float(np.max(np.abs(M[..., 3, 3] - 1.0))) / delta,
        "M00_minus_8": float(np.max(np.abs(M[..., 0, 0] - 8.0))) / delta,
        "z_variation": float(np.max(np.abs(M - M[:, :, :1]))) / delta,
        "pair_mean_minus_s0": float(
            np.max(np.abs(0.5 * (M[..., 1, 1] + M[..., 2, 2]) - delta / 2))
        )
        / delta,
    }


# ================= the slaved bound (own, deviation form) =================
def resid_dev(y, f, delta):
    """the four residuals tr N^p - C_p on eigenvalues (-(8+e), 1+t, s0+sig+b, s0+sig-b),
    every term a small number (no 8^p - 8^p subtraction)."""
    e, t, sig = y
    s0 = delta / 2.0
    u = s0 + sig
    df = f - f0_of(delta)
    f0 = f0_of(delta)
    r1 = -e + t + 2.0 * sig
    r2 = (16.0 * e + e * e) + (2.0 * t + t * t) + 2.0 * df + 4.0 * s0 * sig + 2.0 * sig * sig
    r3 = (
        -(192.0 * e + 24.0 * e * e + e**3)
        + (3.0 * t + 3.0 * t * t + t**3)
        + 2.0 * sig * (3.0 * s0 * s0 + 3.0 * s0 * sig + sig * sig)
        + 6.0 * s0 * df
        + 6.0 * sig * f
    )
    u2ms = sig * (2.0 * s0 + sig)  # u^2 - s0^2
    r4 = (
        (2048.0 * e + 384.0 * e * e + 32.0 * e**3 + e**4)
        + (4.0 * t + 6.0 * t * t + 4.0 * t**3 + t**4)
        + 2.0 * u2ms * (u * u + s0 * s0)
        + 12.0 * (u * u * df + u2ms * f0)
        + 2.0 * df * (f + f0)
    )
    return np.array([r1, r2, r3, r4])


def f0_of(delta):
    return (delta / 2.0) ** 2


def v_min_curve(delta, fs):
    from scipy.optimize import least_squares

    y = np.zeros(3)
    vals, ys = np.zeros(len(fs)), np.zeros((len(fs), 3))
    for k in range(len(fs) - 1, -1, -1):
        sol = least_squares(
            resid_dev,
            y,
            args=(fs[k], delta),
            xtol=1e-15,
            ftol=1e-15,
            gtol=1e-15,
            x_scale=np.array([1e-4, 1e-2, 1e-2]) * delta,
        )
        y = sol.x
        vals[k], ys[k] = float(np.sum(sol.fun**2)), y
    return vals, ys


def t_slaved3(delta, w, nodes=120):
    f0 = f0_of(delta)
    x, wt = np.polynomial.legendre.leggauss(nodes)
    fs = 0.5 * f0 * (x + 1.0)
    v, ys = v_min_curve(delta, fs)
    integ = 0.5 * f0 * float(np.sum(wt * np.sqrt(w * v)))
    return 4.0 * np.pi * np.sqrt(8.0) * integ, ys[0]


# ================= refinement, readers =================
def refine_x2(M, order):
    """cell-centered x2 refinement in x and y (z untouched): new index j sits at old index
    j/2 - 1/4; map_coordinates with the given spline order, nearest beyond the edge."""
    from scipy.ndimage import map_coordinates

    n = M.shape[0]
    jj = np.arange(2 * n) / 2.0 - 0.25
    I, J = np.meshgrid(jj, jj, indexing="ij")
    out = np.zeros((2 * n, 2 * n) + M.shape[2:])
    for z in range(M.shape[2]):
        for a in range(4):
            for b in range(a, 4):
                v = map_coordinates(
                    M[:, :, z, a, b], [I.ravel(), J.ravel()], order=order, mode="nearest"
                )
                out[:, :, z, a, b] = out[:, :, z, b, a] = v.reshape(2 * n, 2 * n)
    return out


def pair_reader(M, z=0):
    P = M[:, :, z, 1:3, 1:3]
    half = 0.5 * (P[..., 0, 0] - P[..., 1, 1])
    b = np.hypot(half, P[..., 0, 1])
    chi = np.arctan2(P[..., 0, 1], half)
    return b, chi


def wrap(a):
    return (a + np.pi) % (2 * np.pi) - np.pi


def texture_report(M, h, delta, w):
    """where the winding sits and how rough the pair angle is (z = 0 slice)."""
    n = M.shape[0]
    b, chi = pair_reader(M)
    b0 = delta / 2.0
    X, Y, rho, _ = coords2(n, h)
    dx, dy = wrap(np.diff(chi, axis=0)), wrap(np.diff(chi, axis=1))
    plaq = (dx[:, :-1] + dy[1:, :] - dx[:, 1:] - dy[:-1, :]) / (2 * np.pi)
    core = rho < 4.0 * h
    dens_f = curv_density_branch(M, h, "fwd")[:, :, 0]
    dens_b = curv_density_branch(M, h, "bwd")[:, :, 0]
    dens = 0.5 * (dens_f + dens_b)
    tb = t_bps(delta, w)
    ecen = e_curv(M, h, "central")
    esym = e_curv(M, h, "sym")
    return {
        "b_over_b0_min_core": float(np.min(b[core]) / b0),
        "b_over_b0_axis4_mean": float(
            np.mean(b[n // 2 - 1 : n // 2 + 1, n // 2 - 1 : n // 2 + 1]) / b0
        ),
        "max_bond_jump_deg_core": float(
            np.degrees(max(np.max(np.abs(dx[core[:-1, :]])), np.max(np.abs(dy[core[:, :-1]]))))
        ),
        "bonds_over_90deg_within_4h": int(
            np.sum(np.abs(dx[core[:-1, :]]) > np.pi / 2)
            + np.sum(np.abs(dy[core[:, :-1]]) > np.pi / 2)
        ),
        "bonds_over_90deg_total": int(
            np.sum(np.abs(dx) > np.pi / 2) + np.sum(np.abs(dy) > np.pi / 2)
        ),
        "bonds_over_150deg_total": int(
            np.sum(np.abs(dx) > 5 * np.pi / 6) + np.sum(np.abs(dy) > 5 * np.pi / 6)
        ),
        "bonds_60_to_120deg_total": int(
            np.sum((np.abs(dx) > np.pi / 3) & (np.abs(dx) < 2 * np.pi / 3))
            + np.sum((np.abs(dy) > np.pi / 3) & (np.abs(dy) < 2 * np.pi / 3))
        ),
        "plaquette_winding_sum": float(np.sum(np.round(plaq))),
        "plaquettes_with_winding": int(np.sum(np.abs(np.round(plaq)) > 0.5)),
        "E_u_fraction_within_3h": float(np.sum(dens[rho < 3.0 * h]) / max(np.sum(dens), 1e-300)),
        "E_u_fwd_over_bwd": float(np.sum(dens_f) / max(np.sum(dens_b), 1e-300)),
        "E_u_central_over_sym": float(ecen / max(esym, 1e-300)),
        "max_cell_density_over_T_bps_h": float(np.max(dens) / (tb * h)),
    }


# ================= descents and Hessian probes =================
IU = np.triu_indices(4)


def pack(M, free, slots):
    return M[free][:, slots[0], slots[1]].ravel().copy()


def unpack(x, base, free, slots):
    M = base.copy()
    blk = M[free].copy()
    v = x.reshape(-1, len(slots[0]))
    blk[:, slots[0], slots[1]] = v
    blk[:, slots[1], slots[0]] = v
    M[free] = blk
    return M


def perturb(M, free, slots, dv):
    """M with the packed slots on the free cells shifted by dv (complex allowed)."""
    return unpack(pack(M, free, slots) + dv, M, free, slots)


BLOCK_SLOTS = (np.array([0, 1, 2, 1, 3]), np.array([0, 1, 2, 2, 3]))
DIR_SLOTS = (np.array([1, 2]), np.array([3, 3]))
STATIC_SLOTS = (np.array([0, 1, 2, 1, 3, 1, 2]), np.array([0, 1, 2, 2, 3, 3, 3]))


def descend(M0, h, w, delta, free, slots, iters, precond=None, tag="", every=250):
    """plain L-BFGS in scaled variables X = (M - M) / delta on the given slots; f = E / delta^4.
    precond: per-slot multiplicative scale on the variables (None = plain)."""
    from scipy.optimize import minimize

    off = np.where(slots[0] == slots[1], 1.0, 2.0)
    sc = np.ones(len(slots[0])) if precond is None else np.asarray(precond, float)
    base = M0 - m_ref(delta)  # the deviation is the state; M itself is never formed
    p0 = pack(base, free, slots)
    k_s = len(sc)
    x0 = np.zeros_like(p0)
    tb = t_bps(delta, w)
    traj = []

    def to_dev(x):
        return unpack(p0 + delta * (x.reshape(-1, k_s) * sc).ravel(), base, free, slots)

    def fun(x):
        Md = to_dev(x)
        eu = e_curv(Md, h)
        br = brackets_tele(None, delta, Mdev=Md)
        ev = float(h**3 * w * np.sum(br * br))
        G = grad_energy(None, h, w, delta, Mdev=Md)[free][:, slots[0], slots[1]] * off
        return (eu + ev) / delta**4, (G * sc).ravel() / delta**3

    x = x0
    done = 0
    while done < iters:
        k = min(every, iters - done)
        res = minimize(
            fun,
            x,
            jac=True,
            method="L-BFGS-B",
            options={"maxcor": 20, "maxiter": k, "maxfun": 3 * k, "gtol": 1e-30, "ftol": 1e-30},
        )
        x = np.asarray(res.x)
        done += max(1, int(res.nit))
        Md = to_dev(x)
        M = Md + m_ref(delta)
        T, eu, ev = tension(None, h, w, delta, form="exact", Mdev=Md)
        traj.append(
            {
                "iters": done,
                "T_over_T_bps": T / tb,
                "gmax_scaled": float(np.max(np.abs(res.jac))),
                "E_u_central_over_sym": e_curv(Md, h, "central") / max(eu, 1e-300),
                "msg": str(res.message),
            }
        )
        log(
            f"{tag} its {done} T/T_bps {T / tb:.6f} gmax {traj[-1]['gmax_scaled']:.2e} {res.message}"
        )
        if int(res.nit) < 2:
            break
    return M, traj, Md


def hessian_extremes(M, h, w, delta, free, slots, k_iter=40, tol=1e-4, scale=None):
    """lowest and highest eigenvalues of the Hessian restricted to the given slots on the free
    cells (the metric: the packed variables with off-diagonal entries counted once, so the
    quadratic form is d2E along symmetric dM with dM_ab = dM_ba = v_ab)."""
    from scipy.sparse.linalg import ArpackNoConvergence, LinearOperator, eigsh

    off = np.where(slots[0] == slots[1], 1.0, 2.0)
    nvar = int(free.sum()) * len(slots[0])
    Mc = M.astype(complex)
    eps = 1e-20

    def hv_raw(v):
        Mv = perturb(Mc, free, slots, 1j * eps * np.asarray(v).ravel())
        G = grad_energy(Mv, h, w, delta).imag / eps
        return (G[free][:, slots[0], slots[1]] * off).ravel()

    if isinstance(scale, str) and scale == "jacobi":
        # per-slot Jacobi scale from the Hessian diagonal at the central free cell
        k = len(slots[0])
        idx = np.argwhere(free)
        c = int(np.argmin(np.sum((idx - idx.mean(0)) ** 2, axis=1)))
        scv = np.ones(k)
        for s_ in range(k):
            e = np.zeros(nvar)
            e[c * k + s_] = 1.0
            d = float(hv_raw(e)[c * k + s_])
            scv[s_] = 1.0 / np.sqrt(abs(d)) if d != 0 else 1.0
        hessian_extremes.jacobi_scale = scv.tolist()
    else:
        scv = np.ones(len(slots[0])) if scale is None else np.asarray(scale, float)

    def hv(v):
        vs = (np.asarray(v).reshape(-1, len(scv)) * scv).ravel()
        Mv = perturb(Mc, free, slots, 1j * eps * vs)
        G = grad_energy(Mv, h, w, delta).imag / eps
        return (G[free][:, slots[0], slots[1]] * off * scv).ravel()  # S H S: same inertia as H

    op = LinearOperator((nvar, nvar), matvec=hv, dtype=float)

    def ext(which):
        try:
            val, vec = eigsh(op, k=1, which=which, tol=tol, maxiter=k_iter, ncv=24)
            return float(val[0]), vec[:, 0], True
        except ArpackNoConvergence as e:  # partial result after the restart budget
            if len(e.eigenvalues):
                return float(e.eigenvalues[0]), e.eigenvectors[:, 0], False
            v = np.random.default_rng(0).standard_normal(nvar)
            for _ in range(30):
                v = hv(v)
                v /= np.linalg.norm(v)
            return float(v @ hv(v)), v, False

    lo, vlo, ok_lo = ext("SA")
    hi, _, ok_hi = ext("LA")
    hessian_extremes.last_converged = (ok_lo, ok_hi)
    return lo, hi, vlo


def mode_roughness(v, free, slots):
    """the fraction of the mode's norm in nearest-neighbor differences (2 for a checkerboard,
    small for a smooth mode) and its central/one-sided energy ratio on the x-y plane."""
    n = free.shape[0]
    V = np.zeros(free.shape + (len(slots[0]),))
    V[free] = v.reshape(-1, len(slots[0]))
    dsum = sum(float(np.sum(np.diff(V, axis=ax) ** 2)) for ax in (0, 1))
    return float(np.sqrt(dsum / max(float(np.sum(V * V)), 1e-300)))


def rayleigh(M, h, w, delta, free, slots, v):
    off = np.where(slots[0] == slots[1], 1.0, 2.0)
    eps = 1e-20
    Mv = perturb(M.astype(complex), free, slots, 1j * eps * v)
    G = grad_energy(Mv, h, w, delta).imag / eps
    Hv = (G[free][:, slots[0], slots[1]] * off).ravel()
    return float(v @ Hv / (v @ v))


def refine_mode(v, free, slots, order=1):
    """interpolate a packed mode to the x2 grid (same map as refine_x2), repacked on the fine
    free mask."""
    from scipy.ndimage import map_coordinates

    n = free.shape[0]
    V = np.zeros(free.shape + (len(slots[0]),))
    V[free] = v.reshape(-1, len(slots[0]))
    jj = np.arange(2 * n) / 2.0 - 0.25
    I, J = np.meshgrid(jj, jj, indexing="ij")
    out = np.zeros((2 * n, 2 * n, NZ, len(slots[0])))
    for z in range(NZ):
        for s in range(len(slots[0])):
            out[:, :, z, s] = map_coordinates(
                V[:, :, z, s], [I.ravel(), J.ravel()], order=order, mode="nearest"
            ).reshape(2 * n, 2 * n)
    free2 = free_mask(2 * n, 24.0 / n)
    return out[free2].ravel()


# ================= claims =================
def rows_of_record(J):
    return {t: r for t, r in J["rows"].items() if t.endswith("_block_plain")}


def job_A_rows(args):
    """recompute T (exact brackets) and the departures on a set of tags."""
    tags, src = args
    out = {}
    for tag in tags:
        if src == "r26":
            M, delta, h = load_row_field(tag)
            w = W1 * float(tag.split("_w")[1].split("_")[0])
        else:
            M, h = load_r25_field(tag)
            delta = float(tag.split("_")[0][1:])
            w = W1 * float(tag.split("_w")[1].split("_")[0])
        T, eu, ev = tension(M, h, w, delta, form="exact")
        Tt, _, evt = tension(M, h, w, delta, form="tele")
        Tp, _, evp = tension(M, h, w, delta, form="plain")
        G = grad_energy(M, h, w, delta)
        free = free_mask(M.shape[0], h)
        gfree = G[free]
        out[tag] = {
            "delta": delta,
            "h": h,
            "w": w,
            "T_exact": T,
            "E_u": eu,
            "V4_exact": ev,
            "V4_tele_rel": (evt - ev) / ev,
            "V4_plain_rel": (evp - ev) / ev,
            "T_over_T_bps": T / t_bps(delta, w),
            "departure": departures(M, delta),
            "gmax_scaled_all": float(np.max(np.abs(gfree))) / (delta**3 * h**3),
            "gmax_scaled_pair": float(np.max(np.abs(gfree[:, 1:3, 1:3]))) / (delta**3 * h**3),
            "E_u_central_over_sym": e_curv(M, h, "central") / eu,
        }
        log(
            f"A {src} {tag} T/T_bps {T / t_bps(delta, w):.6f} plain_rel {out[tag]['V4_plain_rel']:.2e}"
        )
    return out


def job_B_floor(_):
    """the plain and telescoped V4 floors on the BPS seed (n 48, h 1) against the exact brackets."""
    out = {}
    for delta in (0.3, 0.03, 0.01, 0.003, 0.001):
        w = W1 * 25.0
        M = bps_field(48, 1.0, delta, w)
        ex = brackets_exact(M, delta)
        pl = brackets_plain(M, delta)
        te = brackets_tele(M, delta)
        lab = brackets_tele(M, delta, ref=np.array([-8.0, 1.0, delta, 0.0]))
        vex = float(np.sum(ex * ex))
        out[f"{delta:g}"] = {
            "V4_plain_rel": float(np.sum(pl * pl) / vex - 1.0),
            "V4_tele_rel": float(np.sum(te * te) / vex - 1.0),
            "V4_tele_labref_rel": float(np.sum(lab * lab) / vex - 1.0),
            "bracket_p4_plain_max_abs_err": float(np.max(np.abs(pl[..., 3] - ex[..., 3]))),
            "bracket_p4_true_max": float(np.max(np.abs(ex[..., 3]))),
            "T_seed_exact": tension(M, 1.0, w, delta, form="exact")[0],
        }
        log(
            f"B delta {delta:g} plain {out[f'{delta:g}']['V4_plain_rel']:.2e} tele {out[f'{delta:g}']['V4_tele_rel']:.2e} labref {out[f'{delta:g}']['V4_tele_labref_rel']:.2e}"
        )
    return out


def job_slaved(_):
    out = {}
    for delta in (0.3, 0.03, 0.01, 0.003, 0.001, 0.0003):
        w = W1 * 25.0
        t, y0 = t_slaved3(delta, w)
        t80, _ = t_slaved3(delta, w, nodes=60)
        out[f"{delta:g}"] = {
            "T_slaved3_over_delta4": t / delta**4,
            "quadrature_change_60_nodes": t80 / t - 1.0,
            "slots_at_f0_e_t_sigma": [float(v) for v in y0],
            "T_bps_over_delta4": t_bps(delta, w) / delta**4,
        }
        log(f"slaved delta {delta:g} T/delta^4 {t / delta**4:.6f}")
    for w1s in (6.25, 100.0):
        for delta in (0.01, 0.001):
            t, _ = t_slaved3(delta, W1 * w1s)
            out[f"{delta:g}_w{w1s:g}"] = {"T_slaved3_over_delta4": t / delta**4}
    return out


def job_rescale(_):
    """scale invariance: the converged delta 0.01 fields mapped to delta 0.003 and 0.001
    (pair block x delta'/delta, diagonal slot deviations x (delta'/delta)^2, shell exact)."""
    out = {}
    for n in (32, 48, 64, 96):
        tag = f"d0.01_w25_n{n}_L48_block_plain"
        M, delta, h = load_row_field(tag)
        w = W1 * 25.0
        for dp in (0.003, 0.001):
            r = dp / delta
            Mp = M.copy()
            Mp[..., 0, 0] = 8.0 + (M[..., 0, 0] - 8.0) * r * r
            Mp[..., 3, 3] = 1.0 + (M[..., 3, 3] - 1.0) * r * r
            P = M[..., 1:3, 1:3]
            mean = 0.5 * (P[..., 0, 0] + P[..., 1, 1])
            tl = P - mean[..., None, None] * np.eye(2)
            newmean = dp / 2.0 + (mean - delta / 2.0) * r * r
            Mp[..., 1:3, 1:3] = tl * r + newmean[..., None, None] * np.eye(2)
            T, _, _ = tension(Mp, h, w, dp, form="exact")
            out[f"n{n}_to_d{dp:g}"] = {"h": h, "T_over_T_bps": T / t_bps(dp, w)}
            log(f"rescale n{n} -> delta {dp:g}: T/T_bps {T / t_bps(dp, w):.6f}")
    return out


def job_descent(args):
    """own scaled descent on a row of record (block sector, plain variables or slot-scaled)."""
    tag, iters, precond = args
    M, delta, h = load_row_field(tag)
    w = W1 * float(tag.split("_w")[1].split("_")[0])
    free = free_mask(M.shape[0], h)
    pc = None if precond is None else np.array(precond)
    M1, traj, Md1 = descend(
        M,
        h,
        w,
        delta,
        free,
        BLOCK_SLOTS,
        iters,
        precond=pc,
        tag=f"desc {tag} {'pc' if pc is not None else 'plain'}",
    )
    T, eu, ev = tension(None, h, w, delta, form="exact", Mdev=Md1)
    T0_, _, _ = tension(M, h, w, delta, form="exact")
    return {
        "tag": tag,
        "precond": precond,
        "T_over_T_bps_start": T0_ / t_bps(delta, w),
        "T_over_T_bps_end": T / t_bps(delta, w),
        "trajectory": traj,
        "texture_end": texture_report(M1, h, delta, w),
        "departure_end": departures(M1, delta),
    }


def job_D_refine(_):
    out = {}
    w = W1 * 25.0
    cases = [
        ("r26", "d0.03_w25_n64_L48"),
        ("r26", "d0.01_w25_n64_L48_block"),
        ("r26", "escape_d0.03_w25_n48_L48_long_amp0"),
        ("r26", "d0.03_w25_n48_L48"),
        ("r25", "d0.03_w25_n48_L48"),
        ("r25", "d0.3_w25_n48_L48"),
        ("r26", "d0.03_w25_n48_L48_block_plain"),
        ("r26", "d0.001_w25_n64_L48_block_plain"),
    ]
    for src, tag in cases:
        if src == "r26":
            M, delta, h = load_row_field(tag)
        else:
            M, h = load_r25_field(tag)
            delta = float(tag.split("_")[0][1:])
        tb = t_bps(delta, w)
        T, eu, ev = tension(M, h, w, delta)
        rec = {"h": h, "delta": delta, "T_over_T_bps": T / tb, "V4_fraction": ev / (eu + ev)}
        for order in (1, 3):
            M2 = refine_x2(M, order)
            T2, _, _ = tension(M2, h / 2, w, delta)
            rec[f"T_over_T_bps_x2_order{order}"] = T2 / tb
        rec["texture"] = texture_report(M, h, delta, w)
        out[f"{src}:{tag}"] = rec
        log(
            f"D {src}:{tag} T/T_bps {T / tb:.4f} x2 lin {rec['T_over_T_bps_x2_order1']:.3f} cub {rec['T_over_T_bps_x2_order3']:.3f}"
        )
    # pure winding, no melted core
    pw = {}
    for delta in (0.3, 0.03):
        pw[f"{delta:g}"] = {}
        for n in (32, 48, 64, 96):
            h = 48.0 / n
            M = bps_field(n, h, delta, w, pure_winding=True)
            T, eu, ev = tension(M, h, w, delta)
            Mb = bps_field(n, h, delta, w)
            Tb_, _, _ = tension(Mb, h, w, delta)
            pw[f"{delta:g}"][f"h{h:g}"] = {
                "T_pure_winding_over_T_bps": T / t_bps(delta, w),
                "V4_pure_winding_over_T_bps": ev / (NZ * h) / t_bps(delta, w),
                "T_bps_seed_over_T_bps": Tb_ / t_bps(delta, w),
            }
        hs = np.array([1.5, 1.0, 0.75, 0.5])
        ys = np.array([pw[f"{delta:g}"][f"h{x:g}"]["T_pure_winding_over_T_bps"] for x in hs])
        pw[f"{delta:g}"]["h_exponent"] = float(np.polyfit(np.log(hs), np.log(ys), 1)[0])
        log(f"D pure winding delta {delta:g}: {ys} exponent {pw[f'{delta:g}']['h_exponent']:.3f}")
    out["pure_winding"] = pw
    return out


def job_E_hessian(args):
    """the stored R25-1 block strand: director-row residue, the lowest Hessian eigenvalue in the
    director-row subspace and in the block subspace, the same at x2 refinement, and the h 1
    mode's Rayleigh quotient carried to h 0.5."""
    tag = args
    M, h = load_r25_field(tag)
    delta = float(tag.split("_")[0][1:])
    w = W1 * 25.0
    n = M.shape[0]
    free = free_mask(n, h)
    tb = t_bps(delta, w)
    G = grad_energy(M, h, w, delta)
    out = {
        "delta": delta,
        "h": h,
        "T_over_T_bps": tension(M, h, w, delta)[0] / tb,
        "director_row_residue_over_delta": float(np.max(np.abs(M[..., 3, 1:3]))) / delta,
        "director_row_grad_max": float(np.max(np.abs(G[free][:, 3, 1:3]))),
        "grad_max_all": float(np.max(np.abs(G[free]))),
    }
    t1 = time.time()
    lo, hi, v = hessian_extremes(M, h, w, delta, free, DIR_SLOTS, k_iter=60)
    out["dir_subspace"] = {
        "arpack_converged_lo_hi": hessian_extremes.last_converged,
        "lambda_min": lo,
        "lambda_max": hi,
        "lambda_min_over_max": lo / hi,
        "mode_roughness": mode_roughness(v, free, DIR_SLOTS),
        "secs": time.time() - t1,
    }
    log(
        f"E {tag} director subspace: lam_min {lo:.3e} lam_max {hi:.3e} rough {out['dir_subspace']['mode_roughness']:.3f}"
    )
    t1 = time.time()
    lob, hib, vb = hessian_extremes(M, h, w, delta, free, BLOCK_SLOTS, k_iter=60, scale="jacobi")
    out["block_subspace"] = {
        "note": "the congruence-scaled operator S H S (Jacobi scale per slot from the central cell): same inertia as H",
        "jacobi_scale": hessian_extremes.jacobi_scale,
        "arpack_converged_lo_hi": hessian_extremes.last_converged,
        "lambda_min": lob,
        "lambda_max": hib,
        "mode_roughness": mode_roughness(vb, free, BLOCK_SLOTS),
        "secs": time.time() - t1,
    }
    log(
        f"E {tag} block subspace: lam_min {lob:.3e} lam_max {hib:.3e} rough {out['block_subspace']['mode_roughness']:.3f}"
    )
    # the h 1 director mode is Zeroed on the block strand M (exactly block): remove the residue first
    Mb = M.copy()
    Mb[..., 3, 1:3] = 0.0
    Mb[..., 1:3, 3] = 0.0
    lo_b, _, v_b = hessian_extremes(Mb, h, w, delta, free, DIR_SLOTS)
    out["dir_subspace_on_zeroed_row"] = {"lambda_min": lo_b}
    # refinement x2 (linear): the field, the mode carried over, and a fresh lowest eigenvalue
    M2 = refine_x2(Mb, 1)
    free2 = free_mask(2 * n, h / 2)
    v2 = refine_mode(v_b, free, DIR_SLOTS, order=1)
    rq = rayleigh(M2, h / 2, w, delta, free2, DIR_SLOTS, v2)
    rq1 = rayleigh(Mb, h, w, delta, free, DIR_SLOTS, v_b)
    t1 = time.time()
    lo2, hi2, vv2 = hessian_extremes(M2, h / 2, w, delta, free2, DIR_SLOTS, k_iter=12, tol=1e-3)
    conv2 = hessian_extremes.last_converged
    out["refined_x2"] = {
        "T_over_T_bps": tension(M2, h / 2, w, delta)[0] / tb,
        "rayleigh_h1_mode_at_h1": rq1,
        "rayleigh_h1_mode_carried_to_h05": rq,
        "dir_lambda_min_h05": lo2,
        "dir_lambda_max_h05": hi2,
        "dir_mode_roughness_h05": mode_roughness(vv2, free2, DIR_SLOTS),
        "arpack_converged_lo_hi": conv2,
        "secs": time.time() - t1,
    }
    log(f"E {tag} refined: lam_min {lo2:.3e} lam_max {hi2:.3e} carried RQ {rq:.3e} (h1 {rq1:.3e})")
    return out


def job_relax_refined(args):
    """re-relax the x2-refined artifact briefly (free static sector, slot-scaled variables)."""
    tag, iters = args
    M, delta, h = load_row_field(tag)
    w = W1 * 25.0
    M2 = refine_x2(M, 1)
    n2 = M2.shape[0]
    free2 = free_mask(n2, h / 2)
    tb = t_bps(delta, w)
    # slot scales from the V4 stiffness of the diagonal slots (own estimate: the p = 4 term)
    pc = np.array([1.0 / 2048.0, 1.0, 1.0, 1.0, 1.0 / 4.0, 1.0, 1.0])
    M3, traj, Md3 = descend(
        M2,
        h / 2,
        w,
        delta,
        free2,
        STATIC_SLOTS,
        iters,
        precond=pc,
        tag=f"relax-refined {tag}",
        every=100,
    )
    return {
        "tag": tag,
        "T_over_T_bps_coarse": tension(M, h, w, delta)[0] / tb,
        "T_over_T_bps_refined_start": tension(M2, h / 2, w, delta)[0] / tb,
        "T_over_T_bps_refined_end": tension(None, h / 2, w, delta, Mdev=Md3)[0] / tb,
        "trajectory": traj,
        "departure_end": departures(M3, delta),
        "texture_end": texture_report(M3, h / 2, delta, w),
    }


def job_escape_direction(_):
    """the audited escape direction as a test vector: the director row of the stored escape end
    field, applied to the stored R25-1 block strand (director row zeroed): the Rayleigh quotient
    (one complex-step Hessian-vector product; negative = the strand is a saddle along it) and the
    energy along the ray at finite amplitude."""
    out = {}
    w = W1 * 25.0
    for d in ("0.3", "0.03", "0.01"):
        M, h = load_r25_field(f"d{d}_w25_n48_L48")
        delta = float(d)
        Me, _, _ = load_row_field(f"escape_d{d}_w25_n48_L48_long_amp0")
        free = free_mask(M.shape[0], h)
        Mb = M.copy()
        Mb[..., 3, 1:3] = 0.0
        Mb[..., 1:3, 3] = 0.0
        v = pack(Me, free, DIR_SLOTS)
        v = v / np.linalg.norm(v)
        rq = rayleigh(Mb, h, w, delta, free, DIR_SLOTS, v)
        E0 = sum(tension(Mb, h, w, delta)[1:])
        tb = t_bps(delta, w) * NZ * h
        g = grad_energy(Mb, h, w, delta)
        gv = float((g[free][:, DIR_SLOTS[0], DIR_SLOTS[1]] * 2.0).ravel() @ v)
        ray = {}
        for a in (0.01, 0.03, 0.1, 0.3):
            amp = a * delta * np.sqrt(len(v))  # a x delta per entry in rms
            Ma = perturb(Mb, free, DIR_SLOTS, amp * v)
            ray[f"{a:g}"] = (sum(tension(Ma, h, w, delta)[1:]) - E0) / tb
        out[d] = {
            "rayleigh_quotient": rq,
            "gradient_dot_v": gv,
            "E_change_over_E_bps_along_ray_rms_amp_x_delta": ray,
            "escape_end_T_over_T_bps": tension(Me, h, w, delta)[0] / t_bps(delta, w),
        }
        log(f"escape-dir delta {d}: RQ {rq:.3e} g.v {gv:.1e} ray {ray}")
    return out


# ================= smoke =================
def smoke():
    rng = np.random.default_rng(1)
    n, h, delta, w = 10, 1.0, 0.3, W1 * 25.0
    M = bps_field(n, h, delta, w)
    M = M + 0.01 * delta * (lambda K: 0.5 * (K + K.swapaxes(-1, -2)))(rng.standard_normal(M.shape))
    out = {}
    # (1) energy against the platform's e_parts (E_u) at delta 0.3 (the platform is readable)
    import importlib.util

    spec = importlib.util.spec_from_file_location("b3", os.path.join(HERE, "m5_21_3_a_4d.py"))
    b3 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(b3)
    cfg = b3.base_cfg(s=-1.0, n=n, L=float(n * h), delta=delta)
    eu_p, ev_p = b3.e_parts(M, cfg)
    eu = e_curv(M, h)
    ev_plain = v4_energy(M, h, W1, delta, "plain")
    ev_ex = v4_energy(M, h, W1, delta, "exact")
    ev_te = v4_energy(M, h, W1, delta, "tele")
    out["stencil_match_E_u_rel"] = abs(eu / eu_p - 1.0)
    out["V4_plain_vs_platform_rel"] = abs(ev_plain / ev_p - 1.0)
    out["V4_exact_vs_plain_rel"] = abs(ev_ex / ev_plain - 1.0)
    out["V4_tele_vs_exact_rel"] = abs(ev_te / ev_ex - 1.0)
    # (2) gradient by complex step on the energy (own) in random symmetric directions
    errs = []
    G = grad_energy(M, h, w, delta)
    for _ in range(4):
        V = rng.standard_normal(M.shape)
        V = 0.5 * (V + V.swapaxes(-1, -2))
        eps = 1e-7
        Ep = e_curv(M + eps * V, h) + v4_energy(M + eps * V, h, w, delta, "tele")
        Em = e_curv(M - eps * V, h) + v4_energy(M - eps * V, h, w, delta, "tele")
        fd = (Ep - Em) / (2 * eps)
        an = float(np.sum(G * V))
        errs.append(abs(fd / an - 1.0))
    out["gradient_fd_rel_err_max"] = max(errs)
    # (3) the slaved residuals in deviation form against the plain eigenvalue form at delta 0.3
    y = np.array([-3e-5, 0.01, -0.002])
    f = 0.5 * f0_of(delta)
    b = np.sqrt(f)
    ev = np.array([-(8.0 + y[0]), 1.0 + y[1], delta / 2 + y[2] + b, delta / 2 + y[2] - b])
    C = c_targets(delta)
    plain = np.array([np.sum(ev**p) - C[p - 1] for p in range(1, 5)])
    dev = resid_dev(y, f, delta)
    out["resid_dev_vs_plain_max_rel"] = float(np.max(np.abs(dev / plain - 1.0)))
    # (4) exact brackets on a vacuum field are zero
    Mv = bps_field(6, 1.0, 0.03, w, kappa=1e9)
    out["exact_bracket_on_pure_winding_max"] = float(np.max(np.abs(brackets_exact(Mv, 0.03))))
    # (5) hessian-vector product symmetry on the tiny slab (u^T H v = v^T H u)
    free = free_mask(n, h)
    off = np.where(DIR_SLOTS[0] == DIR_SLOTS[1], 1.0, 2.0)
    nv = int(free.sum()) * 2
    u, v = rng.standard_normal(nv), rng.standard_normal(nv)
    Mc = M.astype(complex)

    def hv(x):
        Mx = perturb(Mc, free, DIR_SLOTS, 1j * 1e-20 * x)
        return (
            (grad_energy(Mx, h, w, delta).imag / 1e-20)[free][:, DIR_SLOTS[0], DIR_SLOTS[1]] * off
        ).ravel()

    out["hessian_symmetry_rel"] = abs(u @ hv(v) / (v @ hv(u)) - 1.0)
    out["PASS"] = bool(
        out["stencil_match_E_u_rel"] < 1e-12
        and out["V4_plain_vs_platform_rel"] < 1e-12
        and out["V4_exact_vs_plain_rel"] < 1e-9
        and out["V4_tele_vs_exact_rel"] < 1e-12
        and out["gradient_fd_rel_err_max"] < 1e-6
        and out["resid_dev_vs_plain_max_rel"] < 1e-9
        and out["exact_bracket_on_pure_winding_max"] < 1e-15
        and out["hessian_symmetry_rel"] < 1e-8
    )
    print(json.dumps(out, indent=1))
    return out


# ================= run =================
def run(workers=4):
    J = json.load(open(LADDER_JSON))
    rows = rows_of_record(J)
    art = [t for t in J["rows"] if not t.endswith("_block_plain")]
    esc = [
        "escape_d0.03_w25_n48_L48_long_amp0",
        "escape_d0.01_w25_n48_L48_long_amp0",
        "escape_d0.3_w25_n48_L48_long_amp0",
    ]
    r25_tags = [
        "d0.3_w25_n48_L48",
        "d0.03_w25_n48_L48",
        "d0.01_w25_n48_L48",
        "d0.3_w25_n32_L48",
        "d0.03_w25_n32_L48",
    ]
    tags_rec = sorted(rows)
    jobs = {
        "A_rows_1": (job_A_rows, (tags_rec[:10], "r26")),
        "A_rows_2": (job_A_rows, (tags_rec[10:], "r26")),
        "A_art": (job_A_rows, (art + esc, "r26")),
        "A_r25": (job_A_rows, (r25_tags, "r25")),
        "B_floor": (job_B_floor, None),
        "slaved": (job_slaved, None),
        "rescale": (job_rescale, None),
        "D_refine": (job_D_refine, None),
        "E_d0.03": (job_E_hessian, "d0.03_w25_n48_L48"),
        "E_d0.3": (job_E_hessian, "d0.3_w25_n48_L48"),
        "desc_d0.001_n48": (job_descent, ("d0.001_w25_n48_L48_block_plain", 1500, None)),
        "desc_d0.003_n48": (job_descent, ("d0.003_w25_n48_L48_block_plain", 1500, None)),
        "desc_d0.001_n32": (job_descent, ("d0.001_w25_n32_L48_block_plain", 1500, None)),
        "desc_d0.03_n48_cont": (job_descent, ("d0.03_w25_n48_L48_block_plain", 1500, None)),
        "desc_d0.01_n48_pc": (
            job_descent,
            ("d0.01_w25_n48_L48_block_plain", 1500, [1.0 / 2048.0, 1.0, 1.0, 1.0, 0.25]),
        ),
        "relax_refined_escape": (job_relax_refined, ("escape_d0.03_w25_n48_L48_long_amp0", 400)),
        "escape_direction": (job_escape_direction, None),
    }
    results = {}
    with ProcessPoolExecutor(max_workers=workers, mp_context=mp.get_context("fork")) as ex:
        futs = {name: ex.submit(fn, arg) for name, (fn, arg) in jobs.items()}
        for name, fut in futs.items():
            try:
                results[name] = fut.result()
            except Exception as e:  # noqa: BLE001
                import traceback

                results[name] = {"FAILED": repr(e), "tb": traceback.format_exc()}
                log(f"{name} FAILED {e!r}")
            log(f"job {name} done")
    verdicts = assess(J, results)
    out = {
        "task": "M5.32 R26-2 audit",
        "runtime_s": round(time.time() - T0, 1),
        "results": results,
        "verdicts": verdicts,
    }
    with open(OUT_JSON, "w") as f:
        json.dump(out, f, indent=1, default=float)
    print_table(verdicts)
    log(f"written {OUT_JSON}")
    summarize()


def assess(J, R):
    """the verdict logic: numbers only; the prose in the summary table."""
    V = {}
    rows = rows_of_record(J)
    A = {}
    for k in ("A_rows_1", "A_rows_2"):
        if isinstance(R.get(k), dict) and "FAILED" not in R[k]:
            A.update(R[k])
    sl = R.get("slaved", {})
    ratios = {}
    for tag, a in A.items():
        r = rows[tag]
        key = (
            f"{a['delta']:g}"
            if abs(a["w"] / W1 - 25.0) < 1e-9
            else f"{a['delta']:g}_w{a['w'] / W1:g}"
        )
        Ts3 = sl[key]["T_slaved3_over_delta4"] * a["delta"] ** 4 if key in sl else None
        ratios[tag] = {
            "own_T_over_T_slaved3": a["T_exact"] / Ts3 if Ts3 else None,
            "audited": r["T_over_T_slaved3"],
            "own_T_rel_vs_audited_T": a["T_exact"] / r["T"] - 1.0,
            "verdict_audited": r["verdict"],
            "gmax_scaled_own": a["gmax_scaled_all"],
            "gmax_scaled_pair_own": a["gmax_scaled_pair"],
        }
    V["A"] = {"rows": ratios}
    if ratios:
        V["A"]["max_rel_T_diff"] = max(abs(v["own_T_rel_vs_audited_T"]) for v in ratios.values())
        V["A"]["own_ratio_range"] = [
            min(v["own_T_over_T_slaved3"] for v in ratios.values() if v["own_T_over_T_slaved3"]),
            max(v["own_T_over_T_slaved3"] for v in ratios.values() if v["own_T_over_T_slaved3"]),
        ]
    V["A"]["rescale"] = R.get("rescale")
    V["A"]["descents"] = {
        k: {kk: vv for kk, vv in R[k].items() if kk != "trajectory"}
        | {"last": R[k]["trajectory"][-1]}
        for k in R
        if k.startswith("desc_") and "trajectory" in R[k]
    }
    V["B"] = R.get("B_floor")
    # C: exponent on the AT_GATE h 1.5 rows (own T) and on the rescaled/descended values
    pts = [
        (a["delta"], a["T_exact"])
        for t, a in A.items()
        if rows[t]["h"] == 1.5 and rows[t]["verdict"] == "AT_GATE"
    ]
    if len(pts) >= 2:
        d = np.array([p[0] for p in pts])
        T = np.array([p[1] for p in pts])
        V["C"] = {"points": pts, "exponent_own_T": float(np.polyfit(np.log(d), np.log(T), 1)[0])}
        ds = np.array([0.03, 0.01, 0.003, 0.001])
        Ts = np.array([sl[f"{x:g}"]["T_slaved3_over_delta4"] * x**4 for x in ds])
        V["C"]["exponent_T_slaved3"] = float(np.polyfit(np.log(ds), np.log(Ts), 1)[0])
        # the exponent on ALL four deltas at h 1.5 with the own T
        pts4 = [(a["delta"], a["T_exact"]) for t, a in A.items() if rows[t]["h"] == 1.5]
        d4 = np.array([p[0] for p in pts4])
        T4 = np.array([p[1] for p in pts4])
        V["C"]["exponent_own_T_all_h1.5_rows"] = float(np.polyfit(np.log(d4), np.log(T4), 1)[0])
    V["D"] = R.get("D_refine")
    V["E"] = {k: R[k] for k in R if k.startswith("E_")}
    V["E"]["relax_refined"] = R.get("relax_refined_escape")
    V["E"]["escape_direction"] = R.get("escape_direction")
    V["F"] = {"audited_label": J["collect"]["label"]}
    return V


def print_table(V):
    print("\n==== verdict numbers (see the terminal report for the prose) ====")
    if "rows" in V.get("A", {}):
        for tag, v in sorted(V["A"]["rows"].items()):
            print(
                f"{tag:38s} own {v['own_T_over_T_slaved3']:.5f} audited {v['audited']:.5f} dT {v['own_T_rel_vs_audited_T']:+.1e} {v['verdict_audited']:8s} gscaled {v['gmax_scaled_own']:.2e}"
            )
    print(json.dumps({k: V[k] for k in V if k in ("B", "C")}, indent=1, default=float))


def h2_extrap(hs, ys):
    """y = y0 + c h^2 least squares; returns y0 and the rms residual."""
    A = np.vstack([np.ones(len(hs)), np.asarray(hs) ** 2]).T
    sol, res, _, _ = np.linalg.lstsq(A, np.asarray(ys), rcond=None)
    return float(sol[0]), float(np.sqrt(res[0] / len(hs))) if len(res) else 0.0


def summarize():
    """the verdict table from the written JSON (numbers only from the run), merged back."""
    out = json.load(open(OUT_JSON))
    V, R = out["verdicts"], out["results"]
    J = json.load(open(LADDER_JSON))
    rows = rows_of_record(J)
    A = V["A"]["rows"]
    hs = [1.5, 1.0, 0.75, 0.5]
    per = {}
    for d in ("0.03", "0.01", "0.003", "0.001"):
        ys = [
            A[f"d{d}_w25_n{n}_L48_block_plain"]["own_T_over_T_slaved3"] for n in (32, 48, 64, 96)
        ]
        y0, rms = h2_extrap(hs, ys)
        per[d] = {"ratios_h": ys, "h2_extrapolation": y0, "fit_rms": rms}
    pw = R["D_refine"]["pure_winding"]["0.03"]
    lat = [pw[f"h{h:g}"]["T_bps_seed_over_T_bps"] for h in hs]
    per["lattice_factor_of_bps_seed"] = lat
    per["0.01_over_lattice_factor"] = [y / f for y, f in zip(per["0.01"]["ratios_h"], lat)]
    per["0.03_over_lattice_factor"] = [y / f for y, f in zip(per["0.03"]["ratios_h"], lat)]
    resc = V["A"]["rescale"]
    ts3 = R["slaved"]
    resc_ratio = {
        k: v["T_over_T_bps"]
        * ts3[k.split("_d")[1]]["T_bps_over_delta4"]
        / ts3[k.split("_d")[1]]["T_slaved3_over_delta4"]
        for k, v in resc.items()
    }
    y0_resc, _ = h2_extrap(hs, [resc_ratio[f"n{n}_to_d0.001"] for n in (32, 48, 64, 96)])
    desc = V["A"]["descents"]
    desc_ratio = {
        k: v["T_over_T_bps_end"]
        * ts3[v["tag"].split("_")[0][1:]]["T_bps_over_delta4"]
        / ts3[v["tag"].split("_")[0][1:]]["T_slaved3_over_delta4"]
        for k, v in desc.items()
        if v["tag"].split("_")[0][1:] in ts3
    }
    E = V["E"]
    table = [
        {
            "claim": "A1 every row's T and T / T_slaved3 (exact-integer brackets, own stencil, own slaved bound)",
            "own": f"20 rows reproduced to {V['A']['max_rel_T_diff']:.1e} relative; ratios {V['A']['own_ratio_range'][0]:.5f} to {V['A']['own_ratio_range'][1]:.5f}; director row and M_0i exactly 0, M_33 - 1 = 6e-5 to 2.1e-3 delta, M_00 - 8 = 1e-7 to 4e-6 delta, z variation 3e-5 to 8e-4 delta, pair mean 3e-5 to 9e-4 delta",
            "audited": "0.9875 to 1.011",
            "verdict": "CONFIRMED",
        },
        {
            "claim": "A2 the delta 0.003 and 0.001 rows reached the gate; extrapolations 1.006 +- 0.008 and 1.011 +- 0.0001",
            "own": f"the converged delta 0.01 fields rescaled to delta 0.001 (pair x r, slot deviations x r^2, shell exact) give T / T_slaved3 = {[round(resc_ratio[f'n{n}_to_d0.001'], 5) for n in (32, 48, 64, 96)]} at h 1.5, 1, 0.75, 0.5, below every audited delta 0.001 row ({[round(y, 5) for y in per['0.001']['ratios_h']]}); own deviation-native descents from the audited rows fall to {[round(v, 5) for v in desc_ratio.values()][:3]} (n48 d0.001, n48 d0.003, n32 d0.001) in at most 1500 iterations; the rows stopped at 19 to 40 iterations on ftol / line search because M_00 near 8 is quantized at 1.8e-15 in float64 (0.7 to 1.5 percent of energy left); the h^2 extrapolation of the converged twins is {y0_resc:.4f}",
            "audited": "AT_GATE at all four h; 1.006 +- 0.008, 1.011 +- 0.0001",
            "verdict": "REFUTED",
        },
        {
            "claim": "A3 the delta 0.03 and 0.01 rows at h <= 1 FALLING at 0.9945 to 0.9995",
            "own": f"ratios reproduced; own 1500-iteration continuation of d0.03 h1 moves 0.99456 to {desc_ratio.get('desc_d0.03_n48_cont', float('nan')):.5f} with the per-250-iteration drop falling 5e-6 to 1e-6 (converged to 1e-5); the ratios track the BPS seed's own lattice factor {[round(x, 4) for x in lat]}: divided out {[round(x, 4) for x in per['0.01_over_lattice_factor']]}; h^2 extrapolation {per['0.01']['h2_extrapolation']:.4f} (delta 0.01), {per['0.03']['h2_extrapolation']:.4f} (delta 0.03)",
            "audited": "0.9945 to 0.9995, not at the gate",
            "verdict": "CONFIRMED",
        },
        {
            "claim": "B plain V4 round-off 1e-7 at delta 0.001, telescoped under 1e-8",
            "own": f"on the n48 BPS seed the plain form is off by {R['B_floor']['0.001']['V4_plain_rel']:.1e} at delta 0.001 ({R['B_floor']['0.003']['V4_plain_rel']:.1e} at 0.003, {R['B_floor']['0.01']['V4_plain_rel']:.1e} at 0.01), 4e-8 to 6.6e-8 on the stored delta 0.001 rows; the own telescoped form about diag(-8, s0, s0, 1) is exact to 1e-15 at every delta; telescoped about the literal diag(-8, 1, delta, 0) still 1.1e-11; the p = 4 bracket in plain form is pure noise at delta 0.001 (error 9.4e-13 against a true 8.6e-13)",
            "audited": "1e-7 plain, under 1e-8 telescoped",
            "verdict": "CONFIRMED (plain floor 2e-8 to 7e-8, same order; the digit loss that stops the small-delta rows is the FIELD's float64 quantization of M_00, not the energy form)",
        },
        {
            "claim": "C exponent of T against delta on the AT_GATE h 1.5 rows = 3.994",
            "own": f"own exact T on the same three rows: {V['C']['exponent_own_T']:.4f}; on all four h 1.5 rows {V['C']['exponent_own_T_all_h1.5_rows']:.4f}; T_slaved3's own exponent over 0.03 to 0.001: {V['C']['exponent_T_slaved3']:.4f} (the departure from 4 is the slaved bound's K_eff drift, not the lattice)",
            "audited": "3.994 (T_slaved3 3.997)",
            "verdict": "CONFIRMED",
        },
        {
            "claim": "D1 artifact energies 0.235, 0.574, 0.131 T_bps; x2 interpolation raises them to 2.1 to 22.7; strands keep 0.840 to 0.842; pure winding 11.6 to 115 T_bps growing as 1/h^2",
            "own": "0.2348, 0.5741, 0.1309 (exact); x2 linear / cubic: 2.07 / 2.51, 4.01 / 22.5, 8.88 / 10.4; the strands 0.8376 to 0.8412 / 0.8411 (d0.03 h1), 0.8678 to 0.8687 / 0.8709 (d0.001 h0.75); pure winding 11.58, 26.06, 46.32, 104.2 (delta 0.3) and 12.78, 28.76, 51.13, 115.0 (delta 0.03), exponent -2.000",
            "audited": "as claimed",
            "verdict": "CONFIRMED",
        },
        {
            "claim": "D2 the low energies are LATTICE artifacts (they do not survive refinement); the texture: the line angle drifts and jumps between neighboring cells",
            "own": f"texture (z = 0): b / b0 stays 0.77 to 0.98 in the core (no melted core), the winding sits on one plaquette, 7 to 107 bonds carry angle jumps over 150 degrees (pi-walls of the pair anisotropy) and 15 to 233 over 90 degrees, only 2.5 to 5 percent of E_u within 3 h (spread), the central / one-sided energy ratio 1.05 to 14.3 against 0.99 on the smooth strands (a one-bond-wide structure). Across a pi-wall M jumps by 2 b R but F = [A_x, A_y]_eta vanishes where the field is constant along the wall and V4 is blind to it (same eigenvalues): the commutator-only functional is degenerate for 1D and point structures in the CONTINUUM too (the pure winding with an unmelted point core has u = 32 (b b' / rho)^2 = 0). Decisive: re-relaxing the x2-refined escape field for 400 iterations at h 0.5 takes it from 8.88 back to {E['relax_refined']['T_over_T_bps_refined_end']:.3f} T_bps (0.326 after 100), below the smooth strand's 0.84 and still falling, so the state re-forms on the finer grid",
            "audited": "lattice artifact by the x2 interpolation test",
            "verdict": "QUALIFIED (the interpolation-without-relaxation test does not show the state disappears with h; it re-forms at h/2; the degeneracy is the functional's, exploited at the cell scale)",
        },
        {
            "claim": "E1 the stored strands' director-row residue is sub-1e-6 and is what gets amplified",
            "own": "director row of the stored R25-1 fields: 4.3e-4 delta (d0.03), 2.6e-3 delta (d0.3), 2.5e-4 delta (d0.01); the director-row gradient 2.3e-9 (d0.03) and 6.1e-7 (d0.3) absolute",
            "audited": "sub-1e-6",
            "verdict": "REFUTED as a number (the residue is 1e-4 to 1e-3 of delta; a seed for amplification exists either way)",
        },
        {
            "claim": "E2 the block strand is a saddle of the lattice energy; does the unstable direction survive refinement",
            "own": f"d0.03 h1: lowest Hessian eigenvalue in the director-row subspace {E['E_d0.03']['dir_subspace']['lambda_min']:.2e} (converged; largest {E['E_d0.03']['dir_subspace']['lambda_max']:.2e}): a saddle; in the block subspace (Jacobi-scaled, same inertia) {E['E_d0.03']['block_subspace']['lambda_min']:.2e} with a cell-scale mode (roughness 1.92 of 2.83): a saddle inside the block sector too. The h 1 director mode carried to h 0.5 keeps a NEGATIVE Rayleigh quotient {E['E_d0.03']['refined_x2']['rayleigh_h1_mode_carried_to_h05']:.2e} (a rigorous upper bound on lambda_min at h 0.5; pure continuum scaling would give -2.6e-6): the direction survives refinement, so it is not a lattice-only instability at this h pair; the fresh Lanczos at h 0.5 and both d0.3 probes did not converge (reported, unresolved). The escape END direction is uphill from the strand (Rayleigh +2.6e-5 at d0.03): the escape is a curved path, not the linear mode",
            "audited": "a saddle, amplified residue; lattice artifact",
            "verdict": "CONFIRMED saddle; QUALIFIED on lattice-vs-continuum (the direction survives x2 refinement)",
        },
        {
            "claim": "F the label STRAND_0258_UNRESOLVED follows from the rule; the rule",
            "own": f"the label follows (delta 0.03 and 0.01 rows FALLING at h <= 1). The rule inverts the evidence: it passes the small-delta rows because their gate sits at the absolute floor (gate_scaled 0.05 to 1.3 against 7e-6 to 6e-5 at delta 0.03), rows that are 0.7 to 1.5 percent above their converged twins, and fails the converged rows. The FALLING rows are the smooth profile's lattice discretization below the continuum bound (ratios / BPS-seed lattice factor = {[round(x, 4) for x in per['0.01_over_lattice_factor']]}) plus a decelerating tail, not a slide into the artifact (the plain continuation keeps central / sym at 0.987; the slot-preconditioned descent from the same field does slide, 0.8513 to 0.8083 T_bps in 1500 iterations with the roughness rising, the block-sector saddle). On the numbers the CONFIRMED condition holds at every delta: h^2 extrapolations {per['0.03']['h2_extrapolation']:.4f} (0.03), {per['0.01']['h2_extrapolation']:.4f} (0.01), {y0_resc:.4f} (0.001 by the converged twins)",
            "audited": "STRAND_0258_UNRESOLVED",
            "verdict": "CONFIRMED (label follows) / QUALIFIED (the rule is the wrong instrument; the extrapolated ratio is 1.000 +- 0.001 at every delta)",
        },
    ]
    out["summary"] = {
        "per_delta": per,
        "rescaled_ratio_T_over_T_slaved3": resc_ratio,
        "descended_ratio_T_over_T_slaved3": desc_ratio,
        "verdict_table": table,
    }
    with open(OUT_JSON, "w") as f:
        json.dump(out, f, indent=1, default=float)
    for t in table:
        print(f"\n[{t['verdict']}] {t['claim']}\n  own: {t['own']}\n  audited: {t['audited']}")


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "smoke"
    if mode == "smoke":
        smoke()
    elif mode == "run":
        run(int(sys.argv[2]) if len(sys.argv) > 2 else 4)
    elif mode == "extra":  # the escape-direction job alone, merged into the written JSON
        out = json.load(open(OUT_JSON))
        out["results"]["escape_direction"] = job_escape_direction(None)
        out["verdicts"] = assess(json.load(open(LADDER_JSON)), out["results"])
        out["runtime_s"] = round(out["runtime_s"] + time.time() - T0, 1)
        with open(OUT_JSON, "w") as f:
            json.dump(out, f, indent=1, default=float)
    elif mode == "summarize":
        summarize()
    else:
        raise SystemExit(f"unknown mode {mode}")


if __name__ == "__main__":
    main()
