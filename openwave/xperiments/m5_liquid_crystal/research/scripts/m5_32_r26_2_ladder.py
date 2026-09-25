"""M5.32 R26-2: the strand's delta ladder below 0.03 on a split-residual potential, in scaled
variables, against the author's small-delta constant 0.258 (T_slaved3 / delta^4 as delta -> 0).

EQUATIONS FIRST
---------------
The R25-1 slab (nz 4, the z stencil open-ended, the x-y shell pinned at depth 1.6), the
certified biaxial vacuum M0 = diag(8, 1, delta, 0), N = M eta, and the static energy
    E = 4 h^3 sum_{i<j} <F_ij, F_ij>_eta + V4,   V4 = w sum_{p<=4} (tr N^p - C_p)^2.
THE DEVIATION FORM. Write M = M0 + delta D (D the scaled deviation, the static sector
D_0i = 0, D_00 free). The curvature term is homogeneous of degree 4 in the first
derivatives, so E_u[M] = delta^4 E_u[D] with no cancellation. The trace brackets are
telescoped, eps = N - N0 = delta D eta:
    tr N^p - C_p = tr N^p - tr N0^p = sum_{k=0}^{p-1} tr(N^k eps N0^{p-1-k}),
every term of order eps, so the 8^p of the time slot is never formed and subtracted (R26-0
check d: the plain form's relative error 4e-7 at delta 0.001, the telescoped 8e-9). The
gradient of V4 uses the same brackets: dV4/dM = 2 w sum_p p B_p sym(eta (N)^{p-1})^T,
dE/dD = delta dE/dM on the spatial entries and on D_00.
THE SCALING. The objective is f(D) = E / T_ref, T_ref = T_bps(delta, w) nz h, so every row's
optimizer sees O(1) numbers; the gate is stated in the M-units of R25-1
    fmax_M < max(1e-6 T_bps / h, 10 x floor),
floor = max abs(g_plain - g_dev) on the BPS seed, the plain form's round-off on the profile
itself (an upper bound on the deviation form's floor, whose own round-off is the smaller by
the gain of check d), converted to the scaled units by delta / T_ref.
THE READOUT. T = E / (nz h); T / delta^4; T / T_bps; T / T_slaved3 where T_slaved3 is the
Bogomolny integral 4 pi sqrt(8) int_0^{f0} sqrt(w V_min(f)) df with the three diagonal slots
(M_00, M_33, the pair mean s) relaxed per f (the R25-1 audit's three-slot potential),
computed HERE in the deviation form (the binomial brackets in e = M_00 - 8 and t = M_33 - 1)
so it does not round off at small delta; its small-delta limit T_slaved3 / delta^4 is the
constant the author quotes as 0.258. Per delta the h -> 0 extrapolation of T / T_slaved3
from h 1.5, 1, 0.75, 0.5 (a fit T0 + a h^2, with the three-term fit's difference as the
error); the exponent of T against delta at the finest h (least squares on the logs; report
018 item 1: the trace-power potential's law is delta^4 from the block identity).

THE ROWS: delta in {0.03, 0.01, 0.003, 0.001} x h in {1.5, 1, 0.75, 0.5} at w1s 25, L 48
(n 32, 48, 64, 96), plus w1s 6.25 and 100 at h 1 for delta 0.01 and 0.001. The seed: the
BPS profile plus the energy-sized kick of R25-1 (0.25 T_bps on the six spatial entries,
M_00 unkicked, the static sector). FIRE 200 (short: R25-1 found it does not move the strand)
then L-BFGS 4000 in chunks of 500 (at the run FIRE was dropped and the variables preconditioned, see slot_scales).

PRE-REGISTERED LABELS
---------------------
    STRAND_0258_CONFIRMED   at every delta the extrapolated T / T_slaved3 within 1 percent of 1
                            with the gate reached on the rows used
    STRAND_0258_REFUTED     off by more than 3 percent at two deltas with the gate reached
    STRAND_0258_UNRESOLVED  otherwise (a gate not reached, or the extrapolation's error above
                            1 percent)
    the exponent is reported, never a label.

Modes: smoke | run [workers] | collect | jobs. Output: data/m5_32_r26_2_ladder.json, arrays in
data/m5_32_r26_2/ (local, kept). Regenerate: run 12 about one hour; smoke about a minute.
"""

import importlib.util
import json
import multiprocessing as mp
import os
import sys
import time
import zlib
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
OUT_JSON = os.path.join(DATA, "m5_32_r26_2_ladder.json")
OUT_NPZ = os.path.join(DATA, "m5_32_r26_2")
NZ = 4
PIN_DEPTH = 1.6
KICK_EXCESS = 0.25
FIRE_ITERS = 0  # dropped at the run, see slot_scales
LBFGS_ITERS = 4000
LBFGS_CHUNK = 500
GATE_REL = 1e-6
FLOOR_FACTOR = 10.0
ETA = np.diag([-1.0, 1.0, 1.0, 1.0])
T0 = time.time()


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


R25_1 = _load("m5_32_r25_1_strand", "m5_32_r25_1_strand.py")
B3, R0, W1 = R25_1.B3, R25_1.R0, R25_1.W1
IU4 = np.triu_indices(4)
OFF4 = np.where(IU4[0] == IU4[1], 1.0, 2.0)


def log(msg):
    print(f"[{time.time() - T0:8.1f}s] {msg}", flush=True)


# ================= the bound, the slaved potential in deviation form =================
def K_lead(s0):
    return 4.0 + 36.0 * s0**2 + 144.0 * s0**4


def T_bps(delta, w):
    s0 = delta / 2.0
    return float(np.pi * np.sqrt(32.0 * w * K_lead(s0)) * s0**4)


def kappa_bps(delta, w):
    return float(np.sqrt(w * K_lead(delta / 2.0) / 32.0))


def _comb(p, k):
    from math import comb

    return comb(p, k)


def resid_dev(y, b, delta):
    """the four trace brackets with M_00 = 8 + e, M_33 = 1 + t, the pair (s + b, s - b),
    every bracket written as a difference of small numbers (no 8^p formed)."""
    e, t, s = y
    out = []
    for p in range(1, 5):
        b0 = sum(_comb(p, k) * (-8.0) ** (p - k) * (-e) ** k for k in range(1, p + 1))
        b3 = sum(_comb(p, k) * t**k for k in range(1, p + 1))
        bp = (s + b) ** p + (s - b) ** p - delta**p
        out.append(b0 + b3 + bp)
    return np.array(out)


def v_slaved3_curve(delta, fs):
    """min over (e, t, s) of sum_p resid_p^2 at each f = b^2, continued from f0 downward."""
    from scipy.optimize import least_squares

    s0 = delta / 2.0
    y = np.array([0.0, 0.0, s0])
    v = np.zeros(len(fs))
    ys = np.zeros((len(fs), 3))
    for k in range(len(fs) - 1, -1, -1):
        b = np.sqrt(max(fs[k], 0.0))
        sol = least_squares(
            resid_dev,
            y,
            args=(b, delta),
            xtol=1e-15,
            ftol=1e-15,
            gtol=1e-15,
            x_scale=np.array([delta**2 / 512.0 + 1e-300, delta**2 + 1e-300, delta + 1e-300]),
        )
        y = sol.x
        v[k] = float(np.sum(sol.fun**2))
        ys[k] = y
    return v, ys


def T_slaved3(delta, w, nodes=80):
    f0 = (delta / 2.0) ** 2
    x, wt = np.polynomial.legendre.leggauss(nodes)
    fs = 0.5 * f0 * (x + 1.0)
    v, ys = v_slaved3_curve(delta, fs)
    integ = 0.5 * f0 * np.sum(wt * np.sqrt(w * v))
    return float(4.0 * np.pi * np.sqrt(8.0) * integ), {
        "axis_e_t_s_minus_s0": [float(ys[0][0]), float(ys[0][1]), float(ys[0][2] - delta / 2)],
        "nodes": nodes,
    }


# ================= the slab =================
def slab_cfg(n, L, delta):
    return R25_1.slab_cfg(n, L, delta)


def M_vac(delta):
    return np.diag([8.0, 1.0, delta, 0.0])


def bps_D(n, nz, h, delta, w):
    """the BPS seed as the scaled deviation D = (M - M0) / delta."""
    M = R25_1.bps_field(n, nz, h, delta, w)
    return (M - M_vac(delta)) / delta


def free_mask_slab(n, nz, h):
    return R25_1.free_mask_slab(n, nz, h)


# ---------- the energy and gradient in the deviation form ----------
def e_curv_D(D, cfg):
    """delta^4 E_u[D] / delta^4: the curvature energy of the deviation (h^3-weighted)."""
    eu, _ = B3.e_parts(D, cfg)
    return float(eu)


def grad_curv(D, cfg):
    """the curvature part of B3.grad (its V4 part left out), on the deviation field."""
    h3 = cfg["h"] ** 3
    G = np.zeros_like(D)
    for br, wt in B3.branches(cfg["stencil"]):
        A = [B3.d1(D, ax, cfg["h"], br) for ax in range(3)]
        dA = [np.zeros_like(D) for _ in range(3)]
        for i in range(3):
            for j in range(i + 1, 3):
                F = B3.comm_eta(A[i], A[j])
                WF = 8.0 * (ETA @ F @ ETA)
                dA[i] += WF @ (ETA @ A[j]).swapaxes(-1, -2) - (A[j] @ ETA).swapaxes(-1, -2) @ WF
                dA[j] += (A[i] @ ETA).swapaxes(-1, -2) @ WF - WF @ (ETA @ A[i]).swapaxes(-1, -2)
        for ax in range(3):
            G += wt * B3.d1_adj(dA[ax], ax, cfg["h"], br)
    return h3 * G


def v4_dev(D, cfg, delta, w, need_grad=True):
    """V4 and dV4/dM by the telescoped brackets on M = M0 + delta D (h^3-weighted)."""
    h3 = cfg["h"] ** 3
    N0 = M_vac(delta) @ ETA
    eps = delta * (D @ ETA)
    N = N0 + eps
    N0p = [np.eye(4)]
    for k in range(1, 4):
        N0p.append(N0p[-1] @ N0)
    Np = [np.broadcast_to(np.eye(4), D.shape).copy()]
    for k in range(1, 4):
        Np.append(Np[-1] @ N)
    B = []
    for p in range(1, 5):
        b = 0.0
        for k in range(p):
            b = b + np.einsum("...ij,...jk,ki->...", Np[k], eps, N0p[p - 1 - k])
        B.append(b)
    V = h3 * w * sum(b * b for b in B)
    Vs = np.sum(V)
    Vs = Vs if np.iscomplexobj(Vs) else float(Vs)
    if not need_grad:
        return Vs, None
    GV = np.zeros_like(D)
    for p in range(1, 5):
        coef = 2.0 * w * B[p - 1] * p
        X = ETA @ Np[p - 1]
        GV += coef[..., None, None] * X.swapaxes(-1, -2)
    return Vs, h3 * B3.sym4(GV)


def energy_grad_M(D, cfg, delta, w, need_grad=True):
    """E[M0 + delta D] and dE/dM in M-units, the static sector (time row zeroed, M_00 free)."""
    eu = delta**4 * e_curv_D(D, cfg)
    ev, gv = v4_dev(D, cfg, delta, w, need_grad)
    if not need_grad:
        return eu, ev, None
    G = delta**4 * grad_curv(D, cfg) / delta + gv  # dE_u/dM = delta^3 grad_curv(D)
    G = R25_1.static_sector(G)
    return eu, ev, G


def energy_grad_plain(D, cfg, delta, w):
    """R25-1's plain evaluation of the same field, for the floor."""
    M = M_vac(delta) + delta * D
    eu, ev = R25_1.energy_parts(M, cfg, w)
    G = R25_1.grad_w(M, cfg, w)
    return eu, ev, G


# ================= optimizers on the scaled objective =================
BLOCK_MASK = np.ones((4, 4))
BLOCK_MASK[3, 1:3] = 0.0
BLOCK_MASK[1:3, 3] = 0.0


class Scaled:
    def __init__(self, cfg, delta, w, free, T_ref, block=False):
        self.cfg, self.delta, self.w, self.free, self.T_ref = cfg, delta, w, free, T_ref
        self.fr = free[..., None, None].astype(float)
        # RUN-TIME ADDITION (2026-09-25 16:40 UTC): the block sector as a constraint. The unconstrained
        # descent leaves the block sector into a SUB-CELL texture whose energy does not survive a
        # grid refinement (refine_test: 0.23 -> 2.1 to 2.5 T_bps, 0.16 -> 9.6 to 10.9), while the
        # block strand keeps its value (0.840 -> 0.842); the director row is frozen at zero here
        # (its gradient masked, so L-BFGS never moves it), which is the class the Bogomolny value
        # belongs to
        if block:
            self.fr = self.fr * BLOCK_MASK
        self.block = block

    def fun_grad(self, D):
        eu, ev, G = energy_grad_M(D, self.cfg, self.delta, self.w)
        # dE/dD = delta dE/dM; scaled by 1 / T_ref
        return (eu + ev) / self.T_ref, (self.delta / self.T_ref) * G * self.fr

    def energy(self, D):
        eu, ev, _ = energy_grad_M(D, self.cfg, self.delta, self.w, need_grad=False)
        return eu + ev


def fire(S, D0, iters, dt0=0.02, dt_max=0.2, tag=""):
    D = D0.copy()
    v = np.zeros_like(D)
    dt, alpha, n_up = dt0, 0.1, 0
    f, g = S.fun_grad(D)
    F = -g
    E0 = f
    for it in range(1, iters + 1):
        P = float(np.sum(F * v))
        if P > 0:
            n_up += 1
            vn, fn = np.sqrt(np.sum(v * v)), np.sqrt(np.sum(F * F))
            v = (1 - alpha) * v + alpha * F * (vn / max(fn, 1e-300))
            if n_up > 5:
                dt, alpha = min(dt * 1.1, dt_max), alpha * 0.99
        else:
            v[:] = 0.0
            dt, alpha, n_up = dt * 0.5, 0.1, 0
        v = v + dt * F
        D = B3.sym4(D + dt * v * S.fr)
        f, g = S.fun_grad(D)
        F = -g
    return D, {"f_start": E0, "f_end": f, "gmax_end": float(np.max(np.abs(F))), "iters": iters}


FREE_IJ = [(0, 0), (1, 1), (1, 2), (1, 3), (2, 2), (2, 3), (3, 3)]  # the static sector's entries
FREE_I = np.array([a for a, b in FREE_IJ])
FREE_J = np.array([b for a, b in FREE_IJ])
FREE_OFF = np.where(FREE_I == FREE_J, 1.0, 2.0)
RT2 = np.sqrt(2.0)


def slot_scales(cfg, delta, w, T_ref):
    """the preconditioner: the diagonal slots' V4 stiffnesses in the scaled units (per cell, M-unit
    stiffness x delta^2 / T_ref): M_00 (e), M_33 (t), the pair mean u = (D_11 + D_22) / sqrt 2;
    the pair anisotropy v = (D_11 - D_22) / sqrt 2 and the off-diagonals are the soft channels and
    keep scale 1. RUN-TIME DEVIATION (2026-09-25): in the scaled variables every diagonal slot's
    stiffness is delta-independent while the strand's is O(1), so the raw condition number is
    1e11 at delta 0.001 (M_00) and 5e5 (M_33, the pair mean); FIRE overflowed on the first row
    and is dropped, L-BFGS runs on these rescaled variables (a diagonal preconditioner, the
    functional untouched)."""
    h3 = cfg["h"] ** 3
    s0 = delta / 2.0
    k_e = 2.0 * w * sum((p * 8.0 ** (p - 1)) ** 2 for p in range(1, 5)) * h3
    k_t = 2.0 * w * sum(p**2 for p in range(1, 5)) * h3
    k_u = 2.0 * w * sum((2.0 * p * s0 ** (p - 1)) ** 2 for p in range(1, 5)) * h3 / 2.0
    sc = np.ones(7)
    for k, val in ((0, k_e), (6, k_t), (1, k_u)):
        sc[k] = np.sqrt(max(val * delta**2 / T_ref, 1.0))
    return sc, {"k_e": k_e, "k_t": k_t, "k_u": k_u}


def pack(D, free, sc):
    y = D[free][:, FREE_I, FREE_J].copy()
    u = (y[:, 1] + y[:, 4]) / RT2
    v = (y[:, 1] - y[:, 4]) / RT2
    y[:, 1], y[:, 4] = u, v
    return (y * sc).ravel()


def unpack(x, D_base, free, sc):
    y = x.reshape(-1, 7) / sc
    d11 = (y[:, 1] + y[:, 4]) / RT2
    d22 = (y[:, 1] - y[:, 4]) / RT2
    D = D_base.copy()
    blk = D[free].copy()
    vals = y.copy()
    vals[:, 1], vals[:, 4] = d11, d22
    blk[:, FREE_I, FREE_J] = vals
    blk[:, FREE_J, FREE_I] = vals
    D[free] = blk
    return D


def grad_pack(g, free, sc):
    gy = g[free][:, FREE_I, FREE_J] * FREE_OFF
    gu = (gy[:, 1] + gy[:, 4]) / RT2
    gv = (gy[:, 1] - gy[:, 4]) / RT2
    gy[:, 1], gy[:, 4] = gu, gv
    return (gy / sc).ravel()


def lbfgs(S, D0, gate_scaled, max_iter, sc, chunk=LBFGS_CHUNK, tag=""):
    from scipy.optimize import minimize

    D = D0.copy()
    done, chunks, verdict = 0, [], "FALLING"
    f_start = S.fun_grad(D)[0]
    while done < max_iter:
        base = D.copy()

        def fun(x, base=base):
            Dx = unpack(x, base, S.free, sc)
            f, g = S.fun_grad(Dx)
            return f, grad_pack(g, S.free, sc)

        k = min(chunk, max_iter - done)
        res = minimize(
            fun,
            pack(D, S.free, sc),
            jac=True,
            method="L-BFGS-B",
            options={"maxcor": 20, "maxiter": k, "maxfun": 3 * k, "gtol": 1e-16, "ftol": 1e-18},
        )
        D = unpack(np.asarray(res.x), base, S.free, sc)
        its = int(res.nit)
        done += max(1, its)
        f, g = S.fun_grad(D)
        gmax = float(np.max(np.abs(g[S.free])))
        chunks.append(
            {"iters": done, "f": f, "gmax": gmax, "chunk_iters": its, "scipy": str(res.message)}
        )
        log(f"{tag} LBFGS its {done} f {f:.9e} gmax {gmax:.2e} (gate {gate_scaled:.1e})")
        if gmax < gate_scaled:
            verdict = "AT_GATE"
            break
        if its < 3:
            verdict = "LINE_SEARCH_STALL"
            break
    return D, {
        "f_start": f_start,
        "f_end": f,
        "gmax_end": gmax,
        "iters": done,
        "verdict": verdict,
        "chunks": chunks,
    }


# ================= jobs =================
def job_tag(j):
    return (
        f"d{j['delta']:g}_w{j['w1s']:g}_n{j['n']}_L{j['L']:g}"
        + ("_block" if j.get("sector") == "block" else "")
        + ("_plain" if j.get("precond") == "plain" else "")
    )


def jobs_all(sector="block", precond="plain"):
    jobs = []
    for delta in (0.03, 0.01, 0.003, 0.001):
        for n in (32, 48, 64, 96):
            jobs.append(
                {
                    "delta": delta,
                    "w1s": 25.0,
                    "n": n,
                    "L": 48.0,
                    "arm": "h",
                    "sector": sector,
                    "precond": precond,
                }
            )
    for delta in (0.01, 0.001):
        for w1s in (6.25, 100.0):
            jobs.append(
                {
                    "delta": delta,
                    "w1s": w1s,
                    "n": 48,
                    "L": 48.0,
                    "arm": "w",
                    "sector": sector,
                    "precond": precond,
                }
            )
    return jobs


def departure(M, delta):
    return R25_1.departure(M, delta)


def run_job(j, fire_iters=FIRE_ITERS, lbfgs_iters=LBFGS_ITERS):
    t0 = time.time()
    tag = job_tag(j)
    n, L, delta, w = j["n"], j["L"], j["delta"], W1 * j["w1s"]
    cfg = slab_cfg(n, L, delta)
    h = cfg["h"]
    row = dict(j, tag=tag, h=h, w=w)
    try:
        free = free_mask_slab(n, NZ, h)
        Tb = T_bps(delta, w)
        T_ref = Tb * NZ * h
        D0 = bps_D(n, NZ, h, delta, w)
        # the floor on the BPS profile: plain against deviation, in M-units
        eu_d, ev_d, g_d = energy_grad_M(D0, cfg, delta, w)
        eu_p, ev_p, g_p = energy_grad_plain(D0, cfg, delta, w)
        floor = float(np.max(np.abs((g_p - g_d)[free])))
        gate_M = max(GATE_REL * Tb / h, FLOOR_FACTOR * floor)
        gate_scaled = gate_M * delta / T_ref
        row.update(
            T_bps=Tb,
            T_ref=T_ref,
            kappa_bps=kappa_bps(delta, w),
            floor_M=floor,
            gate_M=gate_M,
            gate_scaled=gate_scaled,
            gate_at_floor=bool(gate_M == FLOOR_FACTOR * floor),
            E_seed_dev=(eu_d + ev_d),
            E_seed_plain=(eu_p + ev_p),
            seed_rel_diff=abs((eu_d + ev_d) - (eu_p + ev_p)) / max(abs(eu_p + ev_p), 1e-300),
            gmax_seed_M=float(np.max(np.abs(g_d[free]))),
        )
        row["T_seed_bps"] = (eu_d + ev_d) / (NZ * h)
        block = j.get("sector") == "block"
        S = Scaled(cfg, delta, w, free, T_ref, block=block)
        # the kick (R25-1: sized to KICK_EXCESS x T_bps above the seed, six spatial entries)
        rng = np.random.default_rng(zlib.crc32(tag.encode()))
        kick = 0.02 * rng.standard_normal(D0.shape)
        kick[..., 0, :] = 0.0
        kick[..., :, 0] = 0.0
        kick = B3.sym4(kick) * free[..., None, None]
        if block:
            kick = kick * BLOCK_MASK
        for _ in range(3):
            exc = S.energy(D0 + kick) / (NZ * h) - row["T_seed_bps"]
            if exc <= 0:
                break
            kick = kick * np.sqrt(KICK_EXCESS * Tb / exc)
        D0 = D0 + kick
        row["T_seed_kicked"] = S.energy(D0) / (NZ * h)
        os.makedirs(OUT_NPZ, exist_ok=True)
        sc, ks = slot_scales(cfg, delta, w, T_ref)
        if j.get("precond") == "plain":
            # RUN-TIME DEVIATION 3 (2026-09-25 16:55 UTC): the preconditioned L-BFGS finds sub-cell textures
            # below the bound in EVERY sector (refine_test: the block-sector row at delta 0.01, h 0.75 rises
            # from 0.57 to 4.4 and 22.7 T_bps under a x2 refinement); R25-1's plain L-BFGS stayed on the
            # smooth block strand (0.840 -> 0.842 under refinement). The ladder therefore runs R25-1's
            # instrument as planned (plain variables, the deviation-form energy and floor its only change)
            sc = np.ones(7)
        row["preconditioner"] = {
            "scales_e_u_v_t": [float(sc[0]), float(sc[1]), float(sc[4]), float(sc[6])],
            **ks,
        }
        D = D0
        if fire_iters > 0:
            D, fr = fire(S, D0, fire_iters, tag=tag)
            row["fire"] = fr
        D, lb = lbfgs(S, D, gate_scaled, lbfgs_iters, sc, tag=tag)
        row["lbfgs"] = lb
        eu, ev, G = energy_grad_M(D, cfg, delta, w)
        M = M_vac(delta) + delta * D
        T = (eu + ev) / (NZ * h)
        Ts3, s3info = T_slaved3(delta, w)
        row.update(
            E_u_per_len=eu / (NZ * h),
            V4_per_len=ev / (NZ * h),
            T=T,
            T_over_delta4=T / delta**4,
            T_over_T_bps=T / Tb,
            T_slaved3=Ts3,
            T_slaved3_over_delta4=Ts3 / delta**4,
            T_over_T_slaved3=T / Ts3,
            slaved3_axis=s3info,
            fmax_M_end=float(np.max(np.abs(G[free]))),
            departure=departure(M, delta),
            core_radius=R25_1.core_radius(M, cfg, delta),
            core_radius_bps=float(np.sqrt(np.log(2.0) / kappa_bps(delta, w))),
        )
        row["departure_max"] = R25_1.dep_max(row["departure"])
        row["verdict"] = lb["verdict"]
        row["status"] = "OK" if np.all(np.isfinite(D)) else "NONFINITE"
        np.savez_compressed(
            os.path.join(OUT_NPZ, tag + ".npz"), D=D.astype(np.float64), delta=delta
        )
    except Exception as e:  # noqa: BLE001
        import traceback

        row.update(status="FAILED", stop=repr(e), traceback=traceback.format_exc())
    row["wall_s"] = round(time.time() - t0, 1)
    log(
        f"DONE {tag} {row.get('status')} {row.get('verdict')} T/Ts3 {row.get('T_over_T_slaved3')} "
        f"T/d4 {row.get('T_over_delta4')} floor {row.get('floor_M')} gate_at_floor {row.get('gate_at_floor')} "
        f"wall {row['wall_s']}"
    )
    return row


def load_json():
    if os.path.exists(OUT_JSON):
        with open(OUT_JSON) as f:
            return json.load(f)
    return {"task": "M5.32 R26-2", "rows": {}, "collect": {}}


def save_json(J):
    tmp = OUT_JSON + ".tmp"
    with open(tmp, "w") as f:
        json.dump(J, f, indent=1, default=str)
    os.replace(tmp, OUT_JSON)


def run_pool(jobs, workers):
    workers = min(int(workers), 12)
    rows = load_json()["rows"]
    pending = [j for j in jobs if rows.get(job_tag(j), {}).get("status") != "OK"]
    pending.sort(key=lambda j: -j["n"])
    log(f"pool: {len(pending)} jobs, {workers} workers")
    with ProcessPoolExecutor(max_workers=workers, mp_context=mp.get_context("spawn")) as ex:
        futs = [ex.submit(run_job, j) for j in pending]
        for fut in as_completed(futs):
            row = fut.result()
            Jn = load_json()
            Jn["rows"][row["tag"]] = row
            save_json(Jn)
    log("pool done")


# ================= collect =================
def extrapolate(hs, ys):
    """T0 from y = T0 + a h^2 (two-term) with the three-term (+ b h^3) difference as the error."""
    hs, ys = np.asarray(hs, float), np.asarray(ys, float)
    if len(hs) < 3:
        return None, None
    A2 = np.stack([np.ones_like(hs), hs**2], 1)
    c2, *_ = np.linalg.lstsq(A2, ys, rcond=None)
    A3 = np.stack([np.ones_like(hs), hs**2, hs**3], 1)
    c3, *_ = np.linalg.lstsq(A3, ys, rcond=None)
    return float(c2[0]), float(abs(c3[0] - c2[0]))


def collect(rows=None, sector="block", precond="plain"):
    J = load_json()
    rows = J["rows"] if rows is None else rows
    rows = {t: r for t, r in rows.items() if r.get("precond", "scaled") == precond}
    out = {"per_delta": {}, "exponent": {}, "w_arm": {}, "sector": sector, "precond": precond}
    n_ok, n_conf, n_ref = 0, 0, 0
    for delta in (0.03, 0.01, 0.003, 0.001):
        rs = [
            r
            for r in rows.values()
            if r.get("status") == "OK"
            and abs(r["delta"] - delta) < 1e-12
            and r["w1s"] == 25.0
            and r.get("sector") == sector
        ]
        rs.sort(key=lambda r: -r["h"])
        gate_ok = [r for r in rs if r["verdict"] == "AT_GATE"]
        pts = [(r["h"], r["T_over_T_slaved3"], r["verdict"], r["gate_at_floor"]) for r in rs]
        T0v, err = (
            extrapolate([r["h"] for r in gate_ok], [r["T_over_T_slaved3"] for r in gate_ok])
            if len(gate_ok) >= 3
            else (None, None)
        )
        rec = {
            "rows_h_ratio_verdict_gatefloor": pts,
            "n_at_gate": len(gate_ok),
            "extrapolated_T_over_T_slaved3": T0v,
            "extrapolation_error": err,
            "T_slaved3_over_delta4": rs[0]["T_slaved3_over_delta4"] if rs else None,
            "finest_h_T_over_delta4": rs[-1]["T_over_delta4"] if rs else None,
        }
        if T0v is not None and len(gate_ok) == len(rs) == 4 and err is not None and err < 0.01:
            n_ok += 1
            if abs(T0v - 1.0) <= 0.01:
                n_conf += 1
            elif abs(T0v - 1.0) > 0.03:
                n_ref += 1
        out["per_delta"][str(delta)] = rec
    if n_ok == 4 and n_conf == 4:
        label = "STRAND_0258_CONFIRMED"
    elif n_ref >= 2:
        label = "STRAND_0258_REFUTED"
    else:
        label = "STRAND_0258_UNRESOLVED"
    out["label"] = label
    # the exponent at each h (rows at the gate)
    for n in (32, 48, 64, 96):
        rs = [
            r
            for r in rows.values()
            if r.get("status") == "OK"
            and r["n"] == n
            and r["w1s"] == 25.0
            and r["verdict"] == "AT_GATE"
            and r.get("sector") == sector
        ]
        if len(rs) >= 3:
            p = np.polyfit(np.log([r["delta"] for r in rs]), np.log([r["T"] for r in rs]), 1)
            ps = np.polyfit(
                np.log([r["delta"] for r in rs]), np.log([r["T_slaved3"] for r in rs]), 1
            )
            out["exponent"][f"n{n}"] = {
                "T_exponent": float(p[0]),
                "T_slaved3_exponent": float(ps[0]),
                "n_points": len(rs),
            }
    for delta in (0.01, 0.001):
        rs = [
            r
            for r in rows.values()
            if r.get("status") == "OK"
            and abs(r["delta"] - delta) < 1e-12
            and r["n"] == 48
            and r.get("sector") == sector
        ]
        out["w_arm"][str(delta)] = {
            f"w{r['w1s']:g}": (r["T_over_T_slaved3"], r["verdict"]) for r in rs
        }
    out["small_delta_constant"] = {
        str(d): T_slaved3(d, W1 * 25.0)[0] / d**4 for d in (0.03, 0.01, 0.003, 0.001, 0.0003)
    }
    out["small_delta_constant"]["bps_leading"] = float(
        np.pi * np.sqrt(32.0 * W1 * 25.0 * 4.0) / 16.0
    )
    J["collect"] = out
    save_json(J)
    return out


def smoke():
    """the deviation form against the plain one and a complex step; T_slaved3 against the R25-1
    audit's three-slot integral where the plain form is exact; a short n 16 row; the labels."""
    out = {}
    cfg = slab_cfg(16, 24.0, 0.03)
    n, h = 16, cfg["h"]
    w = W1 * 25.0
    D0 = bps_D(n, NZ, h, 0.03, w)
    rng = np.random.default_rng(3)
    D0 = D0 + 0.05 * B3.sym4(rng.standard_normal(D0.shape))
    D0[..., 0, 1:] = 0.0
    D0[..., 1:, 0] = 0.0
    eu_d, ev_d, g_d = energy_grad_M(D0, cfg, 0.03, w)
    eu_p, ev_p, g_p = energy_grad_plain(D0, cfg, 0.03, w)
    out["form_vs_plain_d0.03"] = {
        "E_rel": abs((eu_d + ev_d) - (eu_p + ev_p)) / abs(eu_p + ev_p),
        "V4_rel": abs(ev_d - ev_p) / abs(ev_p),
        "grad_max_abs_diff": float(np.max(np.abs(g_d - g_p))),
        "grad_max_abs": float(np.max(np.abs(g_p))),
    }
    # complex step on the deviation-form energy (the static sector, M_00 free)
    V = B3.sym4(rng.standard_normal(D0.shape))
    V[..., 0, 1:] = 0.0
    V[..., 1:, 0] = 0.0
    eps = 1e-20
    Dc = D0.astype(complex) + 1j * eps * V
    euc = B3.e_parts(Dc, cfg)[0] * 0.03**4
    evc, _ = v4_dev(Dc, cfg, 0.03, w, need_grad=False)
    dE_cs = float(np.imag(euc + evc) / eps)  # dE/dD . V
    dE_an = float(np.sum(g_d * V) * 0.03)  # dE/dM . (delta V)
    out["complex_step"] = {
        "cs": dE_cs,
        "analytic": dE_an,
        "rel": abs(dE_cs - dE_an) / max(abs(dE_cs), 1e-300),
    }
    # T_slaved3 against the audit's plain three-slot integral
    AU = _load("m5_32_r25_1_audit", "m5_32_r25_1_audit.py")
    ts = {}
    for d in (0.3, 0.03):
        mine = T_slaved3(d, w)[0]
        theirs = AU.t_slaved(d, w, nslot=3, nodes=80)[0]
        ts[str(d)] = {"mine": mine, "audit": theirs, "rel": abs(mine / theirs - 1.0)}
    out["T_slaved3_vs_audit"] = ts
    out["small_delta_constant"] = {
        str(d): T_slaved3(d, w)[0] / d**4 for d in (0.03, 0.003, 0.001, 0.0001)
    }
    # a short row
    t = time.time()
    global OUT_NPZ
    keep = OUT_NPZ
    OUT_NPZ = keep + "_smoke"
    r = run_job(
        {"delta": 0.001, "w1s": 25.0, "n": 16, "L": 24.0, "arm": "smoke"},
        fire_iters=0,
        lbfgs_iters=150,
    )
    OUT_NPZ = keep
    out["row_n16_d0.001"] = {
        k: r.get(k)
        for k in (
            "status",
            "verdict",
            "T_over_T_slaved3",
            "T_over_delta4",
            "floor_M",
            "gate_M",
            "gate_at_floor",
            "seed_rel_diff",
            "departure_max",
            "stop",
        )
    }
    out["row_wall_s"] = round(time.time() - t, 1)

    # labels on synthetic rows
    def syn(ratio, verdict="AT_GATE"):
        rows = {}
        for d in (0.03, 0.01, 0.003, 0.001):
            for n in (32, 48, 64, 96):
                h = 48.0 / n
                rows[f"d{d}_n{n}"] = {
                    "status": "OK",
                    "delta": d,
                    "w1s": 25.0,
                    "n": n,
                    "h": h,
                    "verdict": verdict,
                    "sector": "block",
                    "precond": "plain",
                    "gate_at_floor": False,
                    "T_over_T_slaved3": ratio + 0.002 * h**2,
                    "T_slaved3_over_delta4": 0.258,
                    "T_over_delta4": 0.258 * ratio,
                    "T": 0.258 * ratio * d**4,
                    "T_slaved3": 0.258 * d**4,
                }
        return rows

    keep_json = OUT_JSON
    globals()["OUT_JSON"] = keep_json.replace(".json", "_smoke_rows.json")
    lab = {
        "conf": collect(syn(1.0))["label"],
        "ref": collect(syn(1.05))["label"],
        "unres": collect(syn(1.0, "FALLING"))["label"],
        "exp": collect(syn(1.0))["exponent"].get("n96", {}).get("T_exponent"),
    }
    if os.path.exists(OUT_JSON):
        os.remove(OUT_JSON)
    globals()["OUT_JSON"] = keep_json
    out["labels"] = lab
    out["labels_ok"] = (
        lab["conf"] == "STRAND_0258_CONFIRMED"
        and lab["ref"] == "STRAND_0258_REFUTED"
        and lab["unres"] == "STRAND_0258_UNRESOLVED"
        and abs(lab["exp"] - 4.0) < 1e-6
    )
    tags = [job_tag(j) for j in jobs_all()]
    out["wiring"] = {"jobs": len(tags), "unique": len(set(tags)) == len(tags)}
    out["PASS"] = bool(
        out["form_vs_plain_d0.03"]["E_rel"] < 1e-9
        and out["complex_step"]["rel"] < 1e-8
        and all(v["rel"] < 1e-6 for v in ts.values())
        and r.get("status") == "OK"
        and out["labels_ok"]
        and out["wiring"]["unique"]
    )
    with open(OUT_JSON.replace(".json", "_smoke.json"), "w") as f:
        json.dump(out, f, indent=1, default=str)
    print(json.dumps(out, indent=1, default=str))
    return out


def escape_test(tag="d0.03_w25_n48_L48", amps=(0.0, 1e-4, 1e-2), iters=1500, suffix=""):
    """RUN-TIME ADDITION (2026-09-25 15:55 UTC): the first ladder rows left the block sector through
    the director row (0.46 delta at delta 0.01, h 0.75) and fell to half the bound. The block sector
    is an invariant manifold of the flow (a block-diagonal field has a block-diagonal gradient), so
    the question is whether R25-1's converged block-sector strand is a saddle (a tiny director-row
    kick grows) or metastable (a finite kick is needed). This seeds R25-1's stored field with a
    director-row kick of amplitude amp x delta (M units, the free cells, the (1,3) and (2,3)
    entries) and descends in the preconditioned variables, recording the energy and the
    departure per chunk."""
    r25 = os.path.join(DATA, "m5_32_r25_1", tag + ".npz")
    M = np.load(r25)["M"]
    parts = tag.split("_")
    delta = float(parts[0][1:])
    w = W1 * float(parts[1][1:])
    n, L = int(parts[2][1:]), float(parts[3][1:])
    cfg = slab_cfg(n, L, delta)
    h = cfg["h"]
    free = free_mask_slab(n, NZ, h)
    Tb = T_bps(delta, w)
    T_ref = Tb * NZ * h
    S = Scaled(cfg, delta, w, free, T_ref)
    sc, ks = slot_scales(cfg, delta, w, T_ref)
    D_base = (M - M_vac(delta)) / delta
    eu, ev, G = energy_grad_M(D_base, cfg, delta, w)
    D_bps = bps_D(n, NZ, h, delta, w)
    _, _, g_d = energy_grad_M(D_bps, cfg, delta, w)
    _, _, g_p = energy_grad_plain(D_bps, cfg, delta, w)
    floor = float(np.max(np.abs((g_p - g_d)[free])))
    gate_scaled = max(GATE_REL * Tb / h, FLOOR_FACTOR * floor) * delta / T_ref
    out = {
        "tag": tag,
        "delta": delta,
        "h": h,
        "T_r25": (eu + ev) / (NZ * h),
        "T_r25_over_T_bps": (eu + ev) / (NZ * h) / Tb,
        "T_slaved3": T_slaved3(delta, w)[0],
        "fmax_M_r25": float(np.max(np.abs(G[free]))),
        "director_row_grad_max": float(np.max(np.abs(G[..., 3, 1:3]))),
        "departure_r25": departure(M, delta),
        "runs": {},
    }
    Ts3 = out["T_slaved3"]
    for amp in amps:
        rng = np.random.default_rng(11)
        D = D_base.copy()
        if amp > 0:
            k = np.zeros_like(D)
            k[..., 3, 1] = k[..., 1, 3] = amp * rng.standard_normal(D.shape[:3])
            k[..., 3, 2] = k[..., 2, 3] = amp * rng.standard_normal(D.shape[:3])
            D = D + k * free[..., None, None]
        traj = []
        done = 0
        Dc = D
        while done < iters:
            Dc, lb = lbfgs(S, Dc, gate_scaled, 250, sc, chunk=250, tag=f"escape amp {amp:g}")
            done += lb["iters"]
            Mc = M_vac(delta) + delta * Dc
            dep = departure(Mc, delta)
            traj.append(
                {
                    "iters": done,
                    "T_over_T_bps": lb["f_end"],
                    "T_over_T_slaved3": lb["f_end"] * Tb / Ts3,
                    "director_row": dep["director_row"],
                    "M33": dep["M33_minus_1"],
                    "z_var": dep["z_variation"],
                }
            )
            if lb["verdict"] in ("LINE_SEARCH_STALL", "AT_GATE"):
                break
        out["runs"][f"amp{amp:g}"] = {
            "start_T_over_T_bps": S.fun_grad(D)[0],
            "trajectory": traj,
            "verdict": lb["verdict"],
            "gate_scaled": gate_scaled,
        }
        if amp == amps[-1] or len(amps) == 1:
            np.savez_compressed(
                os.path.join(OUT_NPZ, f"escape_{tag}{suffix}_amp{amp:g}.npz"), D=Dc, delta=delta
            )
        log(
            f"escape amp {amp:g}: T/T_bps {traj[-1]['T_over_T_bps']:.4f} director_row {traj[-1]['director_row']:.3e} at {done} its"
        )
    with open(OUT_JSON.replace(".json", f"_escape_{tag}{suffix}.json"), "w") as f:
        json.dump(out, f, indent=1, default=str)
    print(json.dumps({k: v for k, v in out.items() if k != "runs"}, indent=1, default=str))
    for a, r in out["runs"].items():
        print(
            a,
            [
                (t["iters"], round(t["T_over_T_bps"], 4), "%.1e" % t["director_row"])
                for t in r["trajectory"]
            ],
        )
    return out


def profile_of(D, delta, cfg):
    """the escaped texture: per radius (the x axis at z = 0, the two rows about y = 0), the director's
    tilt from z, the three spatial eigenvalues, the pair gap and the anisotropy angle; the z variation
    and the winding of the transverse line on the outer loop."""
    M = M_vac(delta) + delta * D
    n, h = cfg["n"], cfg["h"]
    lam, V = np.linalg.eigh(M[..., 1:, 1:])
    d = V[..., :, 2]
    tilt = np.degrees(np.arccos(np.clip(np.abs(d[..., 2]), 0.0, 1.0)))
    x = (np.arange(n) - (n - 1) / 2.0) * h
    c = n // 2
    k0 = M.shape[2] // 2
    rows = []
    for i in range(c, n):
        cells = [(i, c, k0), (i, c - 1, k0)]
        rows.append(
            {
                "x": float(x[i]),
                "tilt_deg": float(np.mean([tilt[q] for q in cells])),
                "lam": [float(np.mean([lam[q][a] for q in cells])) for a in range(3)],
                "pair_gap_over_delta": float(
                    np.mean([(lam[q][1] - lam[q][0]) for q in cells]) / delta
                ),
                "dir_gap": float(np.mean([(lam[q][2] - lam[q][1]) for q in cells])),
                "M33_minus_1_over_delta": float(
                    np.mean([M[q][3, 3] - 1.0 for q in cells]) / delta
                ),
                "M00_minus_8_over_delta": float(
                    np.mean([M[q][0, 0] - 8.0 for q in cells]) / delta
                ),
                "director_row_over_delta": float(
                    np.max([np.abs(M[q][3, 1:3]).max() for q in cells]) / delta
                ),
            }
        )
    # the tilt map on the z = 0 plane: max and its radius, the azimuthal symmetry of the tilt
    tz = tilt[:, :, k0]
    X, Y = np.meshgrid(x, x, indexing="ij")
    rho = np.sqrt(X * X + Y * Y)
    imax = np.unravel_index(int(np.argmax(tz)), tz.shape)
    # the winding of the transverse line (the middle eigenvector) on the loop rho in [12, 13.5] at z = 0
    e1 = V[..., :, 1][:, :, k0]
    ring = (rho > 12.0) & (rho < 13.5)
    ph = np.arctan2(Y[ring], X[ring])
    order = np.argsort(ph)
    th = np.arctan2(e1[ring][order][:, 1], e1[ring][order][:, 0])
    dth = np.diff(np.concatenate([th, th[:1]]))
    dth = dth - np.pi * np.rint(dth / np.pi)
    # where the energy sits: the curvature and potential densities by radius, per stencil branch
    h3 = h**3
    Tb = T_bps(delta, W1 * 25.0)
    dens = {}
    for br, wt in (("fwd", 0.5), ("bwd", 0.5)):
        A = [B3.d1(M, ax, h, br) for ax in range(3)]
        e = 0.0
        for i in range(3):
            for j in range(i + 1, 3):
                F = B3.comm_eta(A[i], A[j])
                e = e + 4.0 * B3.inner_eta(F, F)
        dens[br] = h3 * e
    eu = 0.5 * (dens["fwd"] + dens["bwd"])
    N0 = M_vac(delta) @ ETA
    eps = M @ ETA - N0
    Nn = N0 + eps
    Np = [np.broadcast_to(np.eye(4), M.shape).copy()]
    for k in range(1, 4):
        Np.append(Np[-1] @ Nn)
    N0p = [np.eye(4)]
    for k in range(1, 4):
        N0p.append(N0p[-1] @ N0)
    v = 0.0
    for pp in range(1, 5):
        b = 0.0
        for k in range(pp):
            b = b + np.einsum("...ij,...jk,ki->...", Np[k], eps, N0p[pp - 1 - k])
        v = v + b * b
    ev = h3 * (W1 * 25.0) * v
    rho3 = np.broadcast_to(rho[:, :, None], M.shape[:3])
    nz_len = M.shape[2] * h
    by_r = {}
    for lo, hi in ((0.0, 3.0), (3.0, 6.0), (6.0, 12.0), (12.0, 1e9)):
        m = (rho3 >= lo) & (rho3 < hi)
        by_r[f"rho_{lo:g}_{hi if hi < 1e8 else 'wall':}"] = {
            "E_u_over_T_bps": float(eu[m].sum() / nz_len / Tb),
            "V4_over_T_bps": float(ev[m].sum() / nz_len / Tb),
            "fwd_minus_bwd_over_sum": float(
                (dens["fwd"][m].sum() - dens["bwd"][m].sum())
                / max(dens["fwd"][m].sum() + dens["bwd"][m].sum(), 1e-300)
            ),
        }
    return {
        "axis_profile": rows,
        "energy_by_radius": by_r,
        "T_over_T_bps": float((eu.sum() + ev.sum()) / nz_len / Tb),
        "tilt_max_deg": float(tz[imax]),
        "tilt_max_rho": float(rho[imax]),
        "tilt_max_xy": [float(X[imax]), float(Y[imax])],
        "tilt_mean_rho_lt_6": float(tz[rho < 6.0].mean()),
        "tilt_azimuthal_spread_at_rho_3": float(tz[(rho > 2.5) & (rho < 3.5)].std()),
        "z_variation_over_delta": float(np.max(np.abs(M - M[:, :, :1])) / delta),
        "outer_loop_half_turns": float(np.sum(dth) / np.pi),
        "pair_gap_min_over_delta": float((lam[..., 1] - lam[..., 0]).min() / delta),
        "dir_gap_min": float((lam[..., 2] - lam[..., 1]).min()),
    }


def profiles(paths):
    out = {}
    for f in paths:
        Z = np.load(f)
        b = os.path.basename(f)
        if "D" in Z.files:
            D, delta = Z["D"], float(Z["delta"])
        else:
            delta = float(b.split("_")[0][1:])
            D = (Z["M"] - M_vac(delta)) / delta
        n = D.shape[0]
        L = float(b.split("_L")[1].split("_")[0].split(".npz")[0]) if "_L" in b else 48.0
        cfg = slab_cfg(n, L, delta)
        pr = profile_of(D, delta, cfg)
        out[b] = pr
        print(
            b,
            "tilt max %.2f deg at rho %.1f, mean(rho<6) %.2f, spread %.2f, z_var %.1e, outer loop %.1f half-turns, gap min %.3f, dir gap min %.3f"
            % (
                pr["tilt_max_deg"],
                pr["tilt_max_rho"],
                pr["tilt_mean_rho_lt_6"],
                pr["tilt_azimuthal_spread_at_rho_3"],
                pr["z_variation_over_delta"],
                pr["outer_loop_half_turns"],
                pr["pair_gap_min_over_delta"],
                pr["dir_gap_min"],
            ),
        )
        print(
            "   energy by radius (E_u, V4)/T_bps:",
            {
                k: (round(v["E_u_over_T_bps"], 3), round(v["V4_over_T_bps"], 3))
                for k, v in pr["energy_by_radius"].items()
            },
            "T/T_bps",
            round(pr["T_over_T_bps"], 4),
        )
        print(
            "   x: tilt / pair_gap / dir_row:",
            [
                (
                    r["x"],
                    round(r["tilt_deg"], 2),
                    round(r["pair_gap_over_delta"], 3),
                    round(r["director_row_over_delta"], 3),
                )
                for r in pr["axis_profile"][:8]
            ],
        )
    with open(OUT_JSON.replace(".json", "_profiles.json"), "w") as f:
        json.dump(out, f, indent=1, default=str)
    return out


def winding_residue():
    """RUN-TIME ADDITION (2026-09-25 16:20 UTC). The escaped rows keep the pair anisotropic through the
    axis and carry their energy as a spread residue with a negligible potential. The curvature of the
    certified action is a commutator of two first derivatives, zero on any one-coordinate texture, so
    the PURE WINDING (b = b0 everywhere, the line singular on the axis, no melted core) costs nothing
    in the continuum; the lattice pays a discretization residue (R25-0 check b). This evaluates that
    residue on the slab at four spacings and four deltas, no descent, and beside it the BPS field's
    energy: if the residue falls with h, the tension of the half-line under this action has no floor
    above the residue and the block-sector Bogomolny value is the energy of a smooth-core SADDLE.
    """
    out = {}
    for delta in (0.3, 0.03, 0.01, 0.001):
        w = W1 * 25.0
        Tb = T_bps(delta, w)
        rec = {}
        for n in (32, 48, 64, 96):
            cfg = slab_cfg(n, 48.0, delta)
            h = cfg["h"]
            Mpw = R25_1.bps_field(n, NZ, h, delta, w, kappa=1e12)
            Dpw = (Mpw - M_vac(delta)) / delta
            eu, ev, _ = energy_grad_M(Dpw, cfg, delta, w, need_grad=False)
            Mb = R25_1.bps_field(n, NZ, h, delta, w)
            Db = (Mb - M_vac(delta)) / delta
            eub, evb, _ = energy_grad_M(Db, cfg, delta, w, need_grad=False)
            # the residue by radius
            X, Y, rho, phi = R25_1.slab_coords(n, h)
            A = [B3.d1(Mpw, ax, h, "fwd") for ax in range(3)]
            e = 0.0
            for i in range(3):
                for j in range(i + 1, 3):
                    F = B3.comm_eta(A[i], A[j])
                    e = e + 4.0 * B3.inner_eta(F, F)
            dens = h**3 * e
            rho3 = np.broadcast_to(rho[:, :, None], dens.shape)
            byr = {
                f"rho_lt_{R:g}": float(dens[rho3 < R].sum() / (NZ * h) / Tb)
                for R in (1.5, 3.0, 6.0, 12.0)
            }
            rec[f"h{h:g}"] = {
                "T_pure_winding_over_T_bps": (eu + ev) / (NZ * h) / Tb,
                "V4_pure_winding_over_T_bps": ev / (NZ * h) / Tb,
                "T_bps_seed_over_T_bps": (eub + evb) / (NZ * h) / Tb,
                "residue_within_rho_fwd_branch": byr,
            }
            log(
                f"delta {delta:g} h {h:g}: pure winding T/T_bps {rec[f'h{h:g}']['T_pure_winding_over_T_bps']:.4f} (V4 {ev / (NZ * h) / Tb:.1e}), BPS seed {rec[f'h{h:g}']['T_bps_seed_over_T_bps']:.4f}, within rho 3: {byr['rho_lt_3']:.4f}"
            )
        hs = sorted(rec, key=lambda k: float(k[1:]))
        ys = [rec[k]["T_pure_winding_over_T_bps"] for k in hs]
        xs = [float(k[1:]) for k in hs]
        pf = np.polyfit(np.log(xs), np.log(ys), 1)
        rec["h_exponent_of_residue"] = float(pf[0])
        out[str(delta)] = rec
    with open(OUT_JSON.replace(".json", "_winding_residue.json"), "w") as f:
        json.dump(out, f, indent=1, default=str)
    for d, rec in out.items():
        print(
            d,
            "h exponent",
            round(rec["h_exponent_of_residue"], 3),
            {
                k: round(v["T_pure_winding_over_T_bps"], 4)
                for k, v in rec.items()
                if k.startswith("h")
            },
        )
    return out


def refine_test(paths, factor=2):
    """RUN-TIME ADDITION (2026-09-25 16:35 UTC): is the escaped texture a sub-cell arrangement (a
    lattice artifact) or a continuum texture? Interpolate the field to a grid `factor` times finer
    (cubic in x and y, z-invariant) and re-evaluate the energy: a continuum texture keeps its
    energy (the block strand of R25-1 is the control), a sub-cell one changes it by O(1)."""
    from scipy.ndimage import zoom

    out = {}
    for f in paths:
        Z = np.load(f)
        b = os.path.basename(f)
        if "D" in Z.files:
            D, delta = Z["D"], float(Z["delta"])
        else:
            delta = float(b.split("_")[0][1:])
            D = (Z["M"] - M_vac(delta)) / delta
        n = D.shape[0]
        w = W1 * 25.0
        cfg = slab_cfg(n, 48.0, delta)
        h = cfg["h"]
        Tb = T_bps(delta, w)
        eu, ev, _ = energy_grad_M(D, cfg, delta, w, need_grad=False)
        rec = {"h": h, "T_over_T_bps": (eu + ev) / (NZ * h) / Tb}
        for order in (1, 3):
            Dz = D[:, :, :1]
            Df = np.zeros((n * factor, n * factor, NZ, 4, 4))
            for a in range(4):
                for c in range(a, 4):
                    v = zoom(Dz[:, :, 0, a, c], factor, order=order, mode="nearest")
                    Df[:, :, :, a, c] = v[:, :, None]
                    Df[:, :, :, c, a] = v[:, :, None]
            cff = slab_cfg(n * factor, 48.0, delta)
            euf, evf, _ = energy_grad_M(Df, cff, delta, w, need_grad=False)
            rec[f"T_over_T_bps_refined_x{factor}_order{order}"] = (
                (euf + evf) / (NZ * cff["h"]) / Tb
            )
        out[b] = rec
        print(b, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in rec.items()})
    with open(OUT_JSON.replace(".json", "_refine.json"), "w") as f:
        json.dump(out, f, indent=1, default=str)
    return out


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "smoke"
    if mode == "refine":
        refine_test(sys.argv[2:])
        return
    if mode == "winding_residue":
        winding_residue()
        return
    if mode == "profile":
        import glob

        paths = sys.argv[2:] or sorted(glob.glob(os.path.join(OUT_NPZ, "*.npz")))
        profiles(paths)
        return
    if mode == "escape":
        escape_test(sys.argv[2] if len(sys.argv) > 2 else "d0.03_w25_n48_L48")
        return
    if mode == "escape_long":
        # the zero-kick leg only, to the gate or the cap, the end field kept
        escape_test(
            sys.argv[2],
            amps=(0.0,),
            iters=int(sys.argv[3]) if len(sys.argv) > 3 else 8000,
            suffix="_long",
        )
        return
    if mode == "smoke":
        smoke()
    elif mode == "run":
        run_pool(jobs_all("block", "plain"), sys.argv[2] if len(sys.argv) > 2 else 12)
    elif mode == "run_precond":
        run_pool(jobs_all("block", "scaled"), sys.argv[2] if len(sys.argv) > 2 else 12)
    elif mode == "run_free":
        run_pool(jobs_all("free", "scaled"), sys.argv[2] if len(sys.argv) > 2 else 12)
    elif mode == "collect":
        c = collect()
        print(json.dumps(c, indent=1, default=str))
    elif mode == "jobs":
        for j in jobs_all():
            print(job_tag(j), j)
    else:
        raise SystemExit(f"unknown mode {mode}")


if __name__ == "__main__":
    main()
