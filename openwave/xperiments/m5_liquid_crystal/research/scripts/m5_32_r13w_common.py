"""M5.32 R13-W shared instruments: the certified stack (m5_21_3_a_4d.py) plus
(a) the degeneracy projection of the (2,3) eigenpair, (b) a FIRE with a
projection hook and an optional FIXED-J term, (c) the local clock generator,
(d) the seeds.  Every energy is the certified E = E_u + V4 (INS4.e_parts) and
every inertia is INS4.kin_of (cross-checked against the registry in W2).

Conventions (the R10 to R12 protocol): cfg = base_cfg(s = -1, g = 8, delta = 0.3),
vacuum d4 = diag(8, 1, 0.3, 0), sym stencil, pinned vacuum shell depth 1.6 at the
box edge only, FIRE dt0 = 0.01, dt_max = 0.1 (the R13-W packet: 12000 iterations
or fmax < 1e-6).

FIXED J (new here; the record's fixed-J states were family minima, never a free
field relaxation): E_J[M] = E_stat[M] + J^2 / (4 kin[M]), kin = INS4.kin_of(M, a0(M)),
a0 refreshed from M every step and held FROZEN inside the gradient (the stack's
velocity-field convention, INS4.kin_grad):
    dE_J/dM = dE_stat/dM - (J^2 / (4 kin^2)) dkin/dM|_{a0 frozen},
    omega = J / (2 kin)   (the Legendre partner of the frozen-a0 read).
"""
from __future__ import annotations

import importlib.util
import json
import os
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.dirname(HERE)
DATA, PLOTS = os.path.join(RES, "data"), os.path.join(RES, "plots")
CK = os.path.join(RES, "checkpoints", "m5_32_r13w")
os.makedirs(CK, exist_ok=True)
T0 = time.time()


def log(m):
    print(f"[{time.time() - T0:8.1f}s] {m}", flush=True)


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


INS4 = _load("m5_21_3_a_4d", "m5_21_3_a_4d.py")
B8 = _load("m5_21_8_b_lattice", "m5_21_8_b_lattice.py")
LAG = _load("m5_32_lagrangian", "m5_32_lagrangian.py")
G, DELTA, S = 8.0, 0.3, -1.0
G1 = B8.G1                     # (2,3) rotation: the clock
G3 = B8.G3                     # (1,2) rotation: the non-commuting twist


def cfg_of(n, L, g=G, delta=DELTA):
    return INS4.base_cfg(s=S, g=g, n=n, L=float(L), delta=delta)


def rot(Gm, q):
    """exp(q G) for a spatial rotation generator (G^3 = -G); q scalar or array."""
    q = np.asarray(q, dtype=float)
    return (np.eye(4) + np.sin(q)[..., None, None] * Gm
            + (1.0 - np.cos(q))[..., None, None] * (Gm @ Gm))


# ------------------------------------------------ degeneracy instruments
def spatial_eigs(M):
    """eigenvalues of the spatial 3x3 block, ascending: (d3, d2, d1) -> returned as (d1, d2, d3)."""
    w = np.linalg.eigvalsh(M[..., 1:, 1:])
    return w[..., ::-1]


def gap23(M):
    """d2 - d3 of the spatial block (the wall diagnostic; 0.3 in the vacuum, 0 on a degenerate wall)."""
    w = np.linalg.eigvalsh(M[..., 1:, 1:])
    return w[..., 1] - w[..., 0]


def degenerate_project(M, mask):
    """on the cells of `mask`, replace the two smallest spatial eigenvalues by their mean
    (the (2,3) block projected onto its trace part, in the local eigenframe)."""
    out = M.copy()
    sub = M[mask][..., 1:, 1:]
    w, V = np.linalg.eigh(sub)
    m = 0.5 * (w[..., 0] + w[..., 1])
    w = w.copy(); w[..., 0] = m; w[..., 1] = m
    sub2 = np.einsum("...ik,...k,...jk->...ij", V, w, V)
    blk = M[mask].copy()
    blk[..., 1:, 1:] = sub2
    out[mask] = blk
    return out


# ------------------------------------------------ the clock generator
def a0_local(M):
    """a0 = J M - M J, J the rotation generator about the local LEADING spatial
    eigenvector (m5_32_r12_a_ring.a0_local; equals B8.a0_unit on the hedgehog up
    to the eigenvector sign, R12 H12-a).  On the vacuum this is G1 M - M G1."""
    w, V = np.linalg.eigh(M[..., 1:, 1:])
    n1 = V[..., :, -1]
    J = np.zeros(M.shape)
    J[..., 1, 2], J[..., 2, 1] = -n1[..., 2], n1[..., 2]
    J[..., 1, 3], J[..., 3, 1] = n1[..., 1], -n1[..., 1]
    J[..., 2, 3], J[..., 3, 2] = -n1[..., 0], n1[..., 0]
    return J @ M - M @ J


def a0_G1(M):
    return G1 @ M - M @ G1


def kin_density(M, a0, cfg):
    """per-cell h^3-weighted kinetic density (sum = INS4.kin_of)."""
    h3 = cfg["h"] ** 3
    dens = np.zeros(M.shape[:3])
    for br, (A, wt) in INS4.a_fields(M, cfg).items():
        for i in range(3):
            F = INS4.comm_eta(a0, A[i])
            dens += wt * 4.0 * INS4.inner_eta(F, F)
    return h3 * dens


def kin_registry(M, cfg, a0):
    """-4 x the omega^2 coefficient of I1 through the registry (the R10 kin_c2 measure)."""
    p = LAG.default_params(s=S, g=cfg["g"], delta=cfg["delta"])
    _, _, C = LAG.omega_decompose(LAG.REGISTRY["I1"], M, cfg, p, a0)
    return -4.0 * float(C)


# ------------------------------------------------ FIRE with hooks
def fire_proj(M0, cfg, free_mask, max_iter, project=None, J=None, a0_of=None,
              log_every=500, tag="", f_tol=1e-6, plateau=(2000, 1e-10),
              dt0=0.01, dt_max=0.1, diag=None):
    """INS4.fire verbatim in its FIRE logic, plus: `project(M) -> M` applied after
    every step (the degeneracy constraint); `J` with `a0_of(M)` adds the fixed-J
    term J^2/(4 kin) (a0 refreshed each step, frozen in the gradient); `diag(M)`
    adds a dict to every trace row."""
    M = M0.copy()
    if project is not None:
        M = project(M)
    free = free_mask[..., None, None].astype(float)
    v = np.zeros_like(M)
    dt, alpha, n_up = dt0, 0.1, 0
    hist = []

    def parts(Mx):
        e_u, e_v = INS4.e_parts(Mx, cfg)
        if J is None:
            return float(e_u + e_v), float(e_u), float(e_v), None, None
        a0 = a0_of(Mx)
        k = float(INS4.kin_of(Mx, a0, cfg))
        return float(e_u + e_v + J * J / (4.0 * k)), float(e_u), float(e_v), k, a0

    def tot_grad(Mx):
        Gt = INS4.grad(Mx, cfg)
        if J is not None:
            a0 = a0_of(Mx)
            k = float(INS4.kin_of(Mx, a0, cfg))
            Gt = Gt - (J * J / (4.0 * k * k)) * INS4.kin_grad(Mx, a0, cfg)
        return Gt

    def force(Mx):
        Fx = -tot_grad(Mx) * free
        if project is not None:
            # audit J7 (2026-09-02): fmax and P = F.v must see the force in the
            # constraint TANGENT; linearized projection, s small
            s = 1e-6
            Fx = (project(Mx + s * Fx) - Mx) / s
        return Fx
    F = force(M)
    t0 = time.time()
    stop = "max_iter"
    for it in range(1, max_iter + 1):
        P = float(np.sum(F * v))
        if P > 0.0:
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
        if project is not None:
            M = project(M)
        F = force(M)
        fmax = float(np.max(np.abs(F)))
        if not np.isfinite(fmax):
            stop = "non-finite"
            break
        if it % log_every == 0 or it == max_iter:
            E, e_u, e_v, k, _ = parts(M)
            row = {"it": it, "E": E, "E_u": e_u, "V4": e_v, "fmax": fmax, "dt": dt}
            if k is not None:
                row["kin"] = k
                row["omega"] = J / (2.0 * k)
            if diag is not None:
                row.update(diag(M))
            hist.append(row)
            extra = f" kin {k:10.4f} om {row['omega']:.5f}" if k is not None else ""
            print(f"  {tag} it {it:6d} E {E:14.6f} E_u {e_u:10.5f} V4 {e_v:.4e} fmax {fmax:.3e}{extra} "
                  f"[{time.time() - t0:.0f}s]", flush=True)
            # audit J1 (2026-09-02): back must be >= 1, or the row compares with
            # itself and every run stops "plateau" at its first log line when
            # log_every > plateau[0] (the same latent defect exists in INS4.fire)
            back = max(1, plateau[0] // max(log_every, 1))
            if len(hist) > back and abs(E - hist[-1 - back]["E"]) < plateau[1]:
                stop = "plateau"
                break
        if fmax < f_tol:
            stop = "f_tol"
            break
    return M, {"stop": stop, "trace": hist, "wall_s": round(time.time() - t0, 1), "iters": it}


# ------------------------------------------------ seeds
R10_SEED = os.path.join(RES, "checkpoints", "m5_32_r10", "relax_g8_n32_L48_it12000.npy")


def seed_hedgehog(n, L, maxit=3000):
    """the R10 protocol: B8.dressed(cfg, 0) relaxed statically (pin_shell 1.6,
    dt0 0.01, dt_max 0.1); cached.  n = 32, L = 48 reuses R10's 12000-iteration end state."""
    cfg = cfg_of(n, L)
    if n == 32 and abs(L - 48.0) < 1e-9 and os.path.exists(R10_SEED):
        return np.load(R10_SEED), cfg, {"source": "m5_32_r10 relax_g8_n32_L48_it12000.npy"}
    key = f"seed_n{n}_L{L:g}_it{maxit}"
    npy, js = os.path.join(CK, key + ".npy"), os.path.join(CK, key + ".json")
    if os.path.exists(npy):
        return np.load(npy), cfg, json.load(open(js))
    M0 = B8.dressed(cfg, 0.0)
    free = ~INS4.pin_shell(n, cfg["h"])
    M, info = fire_proj(M0, cfg, free, maxit, tag=key, log_every=500)
    e0, e1 = INS4.e_parts(M0, cfg), INS4.e_parts(M, cfg)
    rec = {"n": n, "L": L, "h": cfg["h"], "maxit": maxit, "stop": info["stop"], "wall_s": info["wall_s"],
           "E_u_start": float(e0[0]), "E_u_end": float(e1[0]), "V4_start": float(e0[1]), "V4_end": float(e1[1])}
    np.save(npy, M)
    json.dump(rec, open(js, "w"), indent=1)
    return M, cfg, rec


if __name__ == "__main__":
    import sys
    n, L = int(sys.argv[1]), float(sys.argv[2])
    it = int(sys.argv[3]) if len(sys.argv) > 3 else 3000
    M, cfg, rec = seed_hedgehog(n, L, it)
    log(f"seed n{n} L{L:g}: {rec}")
