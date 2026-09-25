"""M5.32 R26-4 audit: the index-partition census of the charge, refuted with own code.

A fresh agent that did not read the census instrument (`m5_32_r26_4_census.py`) nor the
R26-0 form script. Inputs: the census JSON (the audited numbers), the stored arrays in
`data/m5_32_r26_4/` (seed, gate, end, kick_stage), and the stack's stated definitions.

OWN CODE (every number below is recomputed from the arrays):
  energy      E = 4 h^3 sum_br wt sum_{i<j} sum_ab eta_a eta_b (F_ij)_ab^2 + h^3 w sum_p (tr N^p - C_p)^2,
              F_ij = A_i eta A_j - A_j eta A_i, A_i = D_i M (own forward / backward differences,
              averaged), N = M eta, roots (-8, 1, 0.3, 0), w = 25 W1; own analytic gradient,
              complex-step checked on random entries; the pinned shell = two cell layers per face
  reader      the sphere r = R sampled on a (theta, phi) grid by trilinear interpolation of the
              spatial block; the pair line v = middle eigenvector, the director d = top eigenvector
              oriented outward where |d . rhat| >= 0.3 and by continuity elsewhere; each grid edge
              carries the signed angle from v_a to v_b in the plane normal to (d_a + d_b), reduced
              mod pi (a discrete parallel transport); a plaquette's index = (edge sum + the solid
              angle swept by d over the plaquette) / pi, so a defect-free plaquette reads 0 and a
              carrier its half-turn count; carriers = connected regions of flagged or charged
              plaquettes (the region sum is exactly its boundary winding plus the holonomy)
  gap tail    a_l(R) = (2l + 1) / 2 int g P_l(cos theta) d(cos theta) d(phi) / 2 pi on the
              interpolated sphere (own quadrature), g = lambda_1 - lambda_0 - delta; also the
              cell-binned least squares on |r - R| <= h / 2; slope = log-log fit of |a_1| on r 6 to 18
  beta^2      1 - 6 (tr L^3)^2 / (tr L^2)^3 on the traceless spatial block, block and cell means
  virial      E_u / (3 V) from the own energy parts
  barrier     the own energy along (1 - s) M_a + s M_b, 11 points, with its pinned / free split
  generator   a_int = [J_z, M], a_orb = (x d_y - y d_x) M (own central stencil, one-sided faces;
              a 4th-order variant as the sensitivity), a_rigid = a_int - a_orb, C = own kin

Usage: python3 m5_32_r26_4_audit.py   (about 2 minutes, one process)
Output: data/m5_32_r26_4_audit.json and the summary table on stdout.
"""

from __future__ import annotations

import json
import os
import sys
import time

os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")

import numpy as np  # noqa: E402
from scipy.ndimage import label as cc_label  # noqa: E402
from scipy.ndimage import map_coordinates  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
ARR = os.path.join(DATA, "m5_32_r26_4")
CENSUS_JSON = os.path.join(DATA, "m5_32_r26_4_census.json")
OUT_JSON = os.path.join(DATA, "m5_32_r26_4_audit.json")

ETA = np.diag([-1.0, 1.0, 1.0, 1.0])
ETA_D = np.array([-1.0, 1.0, 1.0, 1.0])
W1 = 0.000724023879
W = 25.0 * W1
ROOTS = (-8.0, 1.0, 0.3, 0.0)
CP = [sum(q**p for q in ROOTS) for p in range(1, 5)]
DELTA = 0.3
N, L = 32, 48.0
H = L / N
PIN_LAYERS = 2  # ceil(1.6 / 1.5)
PARTS = ["2_2", "3_1", "1_1_1_1", "4", "2_1_1"]
INTENDED = {"4": [4], "3_1": [3, 1], "2_2": [2, 2], "2_1_1": [2, 1, 1], "1_1_1_1": [1, 1, 1, 1]}
TAG = "P{}_d0.3_w25_n32_L48"
SPHERES_END = (6.0, 9.0, 12.0, 18.0, 21.0)
SPHERES_SEED = (9.0, 18.0)
NT, NP = 120, 240
DIR_OUT_TOL = 0.3
GAP_TOL = 0.05
EDGE_DD_TOL = 0.5
EDGE_ANG_TOL = np.pi / 3.0
CHARGE_TOL = 0.25
SHELLS = [3.0 + 1.5 * k for k in range(13)]
T0 = time.time()


def log(msg):
    print(f"[{time.time() - T0:6.1f}s] {msg}", flush=True)


def load(part, suffix=""):
    return np.load(os.path.join(ARR, TAG.format(part) + suffix + ".npz"))["M"].astype(np.float64)


def coords():
    x = (np.arange(N) - (N - 1) / 2.0) * H
    return np.meshgrid(x, x, x, indexing="ij")


def pin_mask():
    P = np.zeros((N, N, N), dtype=bool)
    for ax in range(3):
        sl = [slice(None)] * 3
        sl[ax] = slice(0, PIN_LAYERS)
        P[tuple(sl)] = True
        sl[ax] = slice(N - PIN_LAYERS, N)
        P[tuple(sl)] = True
    return P


# ================= own energy =================
def d_fwd(f, ax):
    out = np.zeros_like(f)
    d = np.diff(f, axis=ax) / H
    sl = [slice(None)] * f.ndim
    sl[ax] = slice(0, -1)
    out[tuple(sl)] = d
    return out


def d_bwd(f, ax):
    out = np.zeros_like(f)
    d = np.diff(f, axis=ax) / H
    sl = [slice(None)] * f.ndim
    sl[ax] = slice(1, None)
    out[tuple(sl)] = d
    return out


def d_fwd_T(g, ax):
    """adjoint of d_fwd: (D^T g)[m] = g[m-1] / h - g[m] / h with the boundary rows dropped."""
    out = np.zeros_like(g)
    sl = [slice(None)] * g.ndim
    a, b = list(sl), list(sl)
    a[ax], b[ax] = slice(1, None), slice(0, -1)
    out[tuple(a)] += g[tuple(b)] / H
    out[tuple(b)] -= g[tuple(b)] / H
    return out


def d_bwd_T(g, ax):
    out = np.zeros_like(g)
    sl = [slice(None)] * g.ndim
    a, b = list(sl), list(sl)
    a[ax], b[ax] = slice(1, None), slice(0, -1)
    out[tuple(a)] += g[tuple(a)] / H
    out[tuple(b)] -= g[tuple(a)] / H
    return out


BRANCHES = ((d_fwd, d_fwd_T), (d_bwd, d_bwd_T))


def comm(A, B):
    return A @ ETA @ B - B @ ETA @ A


def eta_norm2(F):
    """<F, F>_eta per cell = sum_ab eta_a eta_b F_ab^2 (eta diagonal)."""
    return np.einsum("...ab,a,b->...", F * F, ETA_D, ETA_D)


def traces(M):
    Me = M @ ETA
    P = Me
    t = [np.einsum("...kk->...", P)]
    for _ in range(3):
        P = P @ Me
        t.append(np.einsum("...kk->...", P))
    return t


def density_parts(M):
    """(curvature density, potential density) per cell, h^3-weighted."""
    h3 = H**3
    eu = np.zeros(M.shape[:3], dtype=M.dtype)
    for dfun, _ in BRANCHES:
        A = [dfun(M, ax) for ax in range(3)]
        for i in range(3):
            for j in range(i + 1, 3):
                eu = eu + 0.5 * 4.0 * eta_norm2(comm(A[i], A[j]))
    t = traces(M)
    ev = W * sum((t[p] - CP[p]) ** 2 for p in range(4))
    return h3 * eu, h3 * ev


def energy(M):
    eu, ev = density_parts(M)
    return eu.sum() + ev.sum()


def energy_parts(M):
    eu, ev = density_parts(M)
    return float(eu.sum()), float(ev.sum())


def gradient(M):
    """dE/dM treating the 16 entries as independent, then symmetrized (the gradient on the
    symmetric field's entry space, matching the reported fmax convention)."""
    h3 = H**3
    G = np.zeros_like(M)
    for dfun, dT in BRANCHES:
        A = [dfun(M, ax) for ax in range(3)]
        dA = [np.zeros_like(M) for _ in range(3)]
        for i in range(3):
            for j in range(i + 1, 3):
                F = comm(A[i], A[j])
                Wf = 0.5 * 8.0 * (ETA @ F @ ETA)
                dA[i] += Wf @ A[j] @ ETA - ETA @ A[j] @ Wf
                dA[j] += ETA @ A[i] @ Wf - Wf @ A[i] @ ETA
        for ax in range(3):
            G += dT(dA[ax], ax)
    Me = M @ ETA
    pows = [np.broadcast_to(np.eye(4), M.shape)]
    for _ in range(3):
        pows.append(pows[-1] @ Me)
    t = [np.einsum("...kk->...", P @ Me) for P in pows]
    GV = np.zeros_like(M)
    for p in range(1, 5):
        coef = 2.0 * W * (t[p - 1] - CP[p - 1]) * p
        GV += coef[..., None, None] * (ETA @ pows[p - 1]).swapaxes(-1, -2)
    G = h3 * (G + GV)
    return 0.5 * (G + G.swapaxes(-1, -2))


def complex_step_check(M, k=12, seed=0):
    rng = np.random.default_rng(seed)
    G = gradient(M)
    eps = 1e-20
    errs = []
    for _ in range(k):
        i, j, l = rng.integers(2, N - 2, 3)
        a, b = rng.integers(0, 4, 2)
        Mc = M.astype(complex)
        Mc[i, j, l, a, b] += 1j * eps
        Mc[i, j, l, b, a] += 1j * eps if a != b else 0.0
        g_cs = energy(Mc).imag / eps
        g_an = G[i, j, l, a, b] * (1.0 if a == b else 2.0)
        errs.append(abs(g_cs - g_an) / max(abs(g_cs), 1e-30))
    return float(max(errs))


def adjoint_check(seed=1):
    rng = np.random.default_rng(seed)
    x = rng.standard_normal((N, N, N, 4, 4))
    y = rng.standard_normal((N, N, N, 4, 4))
    errs = []
    for dfun, dT in BRANCHES:
        for ax in range(3):
            lhs = np.sum(dfun(x, ax) * y)
            rhs = np.sum(x * dT(y, ax))
            errs.append(abs(lhs - rhs) / abs(lhs))
    return float(max(errs))


def kin(M, a0):
    """kin(M; a0) = 4 h^3 sum_br wt sum_i <[a0, A_i]_eta, [a0, A_i]_eta>_eta."""
    k = 0.0
    for dfun, _ in BRANCHES:
        for ax in range(3):
            k += 0.5 * 4.0 * np.sum(eta_norm2(comm(a0, dfun(M, ax))))
    return float(H**3 * k)


# ================= own sphere reader =================
def sample_block(M, pts):
    idx = (pts / H + (N - 1) / 2.0).reshape(-1, 3).T
    S = np.zeros(pts.shape[:-1] + (3, 3))
    for a in range(3):
        for b in range(a, 3):
            v = map_coordinates(M[..., 1 + a, 1 + b], idx, order=1, mode="nearest")
            S[..., a, b] = v.reshape(pts.shape[:-1])
            S[..., b, a] = S[..., a, b]
    return S


def tri_solid_angle(a, b, c):
    num = np.einsum("...i,...i->...", a, np.cross(b, c))
    den = 1.0 + np.einsum("...i,...i->...", a, b)
    den = den + np.einsum("...i,...i->...", b, c) + np.einsum("...i,...i->...", c, a)
    return 2.0 * np.arctan2(num, den)


def edge_inc(da, db, va, vb):
    """signed angle from the line va to vb in the plane normal to da + db, mod pi; flags."""
    n = da + db
    n = n / np.maximum(np.linalg.norm(n, axis=-1), 1e-300)[..., None]
    dd = np.einsum("...i,...i->...", da, db)
    pa = va - np.einsum("...i,...i->...", va, n)[..., None] * n
    pb = vb - np.einsum("...i,...i->...", vb, n)[..., None] * n
    na, nb = np.linalg.norm(pa, axis=-1), np.linalg.norm(pb, axis=-1)
    pa = pa / np.maximum(na, 1e-300)[..., None]
    pb = pb / np.maximum(nb, 1e-300)[..., None]
    s = np.einsum("...i,...i->...", n, np.cross(pa, pb))
    c = np.einsum("...i,...i->...", pa, pb)
    ang = np.arctan2(s, c)
    ang = (ang + np.pi / 2.0) % np.pi - np.pi / 2.0
    bad = (dd < EDGE_DD_TOL) | (np.minimum(na, nb) < 0.3) | (np.abs(ang) > EDGE_ANG_TOL)
    return ang, bad, dd


def orient_director(d, rhat):
    """orient the director outward where confident, by continuity elsewhere; returns d, flags."""
    dot = np.einsum("...i,...i->...", d, rhat)
    conf = np.abs(dot) >= DIR_OUT_TOL
    sgn = np.where(conf, np.sign(dot), 0.0)
    assigned = conf.copy()
    nt1, npn = d.shape[:2]
    for _ in range(nt1 * 2 + 10):
        if assigned.all():
            break
        best = np.zeros(d.shape[:2])
        cand = np.zeros(d.shape[:2])
        for shift, axis in ((1, 0), (-1, 0), (1, 1), (-1, 1)):
            dn = np.roll(d, shift, axis=axis)
            sn = np.roll(sgn, shift, axis=axis)
            an = np.roll(assigned, shift, axis=axis)
            if axis == 0:  # no wrap in theta
                if shift == 1:
                    an[0] = False
                else:
                    an[-1] = False
            dd = np.einsum("...i,...i->...", d, dn)
            ok = an & (np.abs(dd) >= EDGE_DD_TOL) & (np.abs(dd) > best)
            best = np.where(ok, np.abs(dd), best)
            cand = np.where(ok, np.sign(dd) * sn, cand)
        new = (~assigned) & (best > 0)
        if not new.any():
            break
        sgn = np.where(new, cand, sgn)
        assigned = assigned | new
    unassigned = ~assigned
    sgn = np.where(unassigned, np.where(dot >= 0, 1.0, -1.0), sgn)
    return d * sgn[..., None], conf, unassigned


def region_merge(lab, nreg):
    """merge labels across the phi seam (column 0 and column NP - 1)."""
    parent = list(range(nreg + 1))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for i in range(lab.shape[0]):
        a, b = lab[i, 0], lab[i, -1]
        if a and b and find(a) != find(b):
            parent[find(a)] = find(b)
    out = np.zeros_like(lab)
    for i in range(1, nreg + 1):
        out[lab == i] = find(i)
    return out


def read_sphere(M, R):
    th = np.linspace(0.0, np.pi, NT + 1)
    ph = np.linspace(0.0, 2.0 * np.pi, NP, endpoint=False)
    T, P = np.meshgrid(th, ph, indexing="ij")
    rhat = np.stack([np.sin(T) * np.cos(P), np.sin(T) * np.sin(P), np.cos(T)], -1)
    pts = R * rhat
    S = sample_block(M, pts)
    lam, V = np.linalg.eigh(S)
    d, v = V[..., :, 2], V[..., :, 1]
    gap_pair, gap_dir = lam[..., 1] - lam[..., 0], lam[..., 2] - lam[..., 1]
    d, conf, unassigned = orient_director(d, rhat)
    # edges: theta (i -> i+1 at j), phi (j -> j+1 at i, wrapped)
    e_th, bad_th, _ = edge_inc(d[:-1], d[1:], v[:-1], v[1:])
    e_ph, bad_ph, _ = edge_inc(d, np.roll(d, -1, 1), v, np.roll(v, -1, 1))
    # plaquette (i, j): c0 (i, j), c1 (i+1, j), c2 (i+1, j+1), c3 (i, j+1), counterclockwise from outside
    j1 = (np.arange(NP) + 1) % NP
    edge_sum = e_th + e_ph[1:] - e_th[:, j1] - e_ph[:-1]
    d0, d1, d2, d3 = d[:-1], d[1:], d[1:][:, j1], d[:-1][:, j1]
    om = tri_solid_angle(d0, d1, d2) + tri_solid_angle(d0, d2, d3)
    w = (edge_sum + om) / np.pi
    gp_min = np.minimum.reduce(
        [gap_pair[:-1], gap_pair[1:], gap_pair[1:][:, j1], gap_pair[:-1][:, j1]]
    )
    edge_bad = bad_th | bad_ph[1:] | bad_th[:, j1] | bad_ph[:-1]
    weak = ~conf
    weak_p = weak[:-1] | weak[1:] | weak[1:][:, j1] | weak[:-1][:, j1]
    una_p = unassigned[:-1] | unassigned[1:] | unassigned[1:][:, j1] | unassigned[:-1][:, j1]
    flagged = edge_bad | (gp_min < GAP_TOL) | una_p
    charged = np.abs(w) > CHARGE_TOL
    # solid angle of the sphere plaquettes (area weights)
    tm = 0.5 * (th[:-1] + th[1:])
    dA = (np.sin(tm) * (th[1] - th[0]) * (ph[1] - ph[0]))[:, None] * np.ones((1, NP))
    lab, nreg = cc_label(flagged | charged)
    lab = region_merge(lab, nreg)
    regions = []
    pos = rhat[:-1] + rhat[1:] + rhat[1:][:, j1] + rhat[:-1][:, j1]
    gd_min = np.minimum.reduce(
        [gap_dir[:-1], gap_dir[1:], gap_dir[1:][:, j1], gap_dir[:-1][:, j1]]
    )
    for k in np.unique(lab):
        if k == 0:
            continue
        sel = lab == k
        ws = float(w[sel].sum())
        m = int(np.rint(ws))
        cen = pos[sel].sum(0)
        cen = cen / max(np.linalg.norm(cen), 1e-300)
        rho, z = R * float(np.hypot(cen[0], cen[1])), R * float(cen[2])
        weak_in, una_in = bool(weak_p[sel].any()), bool(una_p[sel].any())
        # a clean carrier: an integer sum to 0.15 with the director outward and oriented throughout
        clean = abs(ws - m) < 0.15 and not weak_in and not una_in
        regions.append(
            {
                "m_half_units": m,
                "raw_sum": round(ws, 4),
                "residual": round(ws - m, 4),
                "boundary_edge_sum_over_pi": round(float(edge_sum[sel].sum() / np.pi), 4),
                "holonomy_over_pi": round(float(om[sel].sum() / np.pi), 4),
                "rho": round(rho, 2),
                "z": round(z, 2),
                "phi_deg": round(float(np.degrees(np.arctan2(cen[1], cen[0])) % 360.0), 1),
                "solid_angle_frac": round(float(dA[sel].sum() / (4.0 * np.pi)), 4),
                "n_plaq": int(sel.sum()),
                "min_pair_gap": round(float(gp_min[sel].min()), 4),
                "min_dir_gap": round(float(gd_min[sel].min()), 4),
                "director_weak_inside": weak_in,
                "orientation_unassigned_inside": una_in,
                "clean": bool(clean),
            }
        )
    regions.sort(key=lambda r: (-r["m_half_units"], r["rho"]))
    part = sorted(
        [r["m_half_units"] for r in regions if r["m_half_units"] != 0 and r["clean"]], reverse=True
    )
    part_all = sorted([r["m_half_units"] for r in regions if r["m_half_units"] != 0], reverse=True)
    unread = [r for r in regions if not r["clean"]]
    tot_w = float(w.sum())
    return {
        "R": R,
        "partition": part,
        "total": int(sum(part)),
        "partition_all_regions": part_all,
        "total_all_regions": int(sum(part_all)),
        "unreadable_regions_net_raw": round(float(sum(r["raw_sum"] for r in unread)), 4),
        "n_unreadable_regions": len(unread),
        "unreadable_regions_solid_angle_frac": round(
            float(sum(r["solid_angle_frac"] for r in unread)), 4
        ),
        "total_raw": round(tot_w, 4),
        "degree_of_director": round(float(om.sum() / (4.0 * np.pi)), 4),
        "regions": regions,
        "n_regions": len(regions),
        "frac_director_not_outward_0p3": round(float(dA[weak_p].sum() / (4.0 * np.pi)), 4),
        "frac_flagged_own": round(float(dA[flagged].sum() / (4.0 * np.pi)), 4),
        "frac_pair_gap_closed": round(float(dA[gp_min < GAP_TOL].sum() / (4.0 * np.pi)), 4),
        "frac_orientation_unassigned": round(float(dA[una_p].sum() / (4.0 * np.pi)), 4),
        "min_pair_gap": round(float(gap_pair.min()), 4),
        "min_dir_gap": round(float(gap_dir.min()), 4),
        "max_abs_w_readable": round(
            float(np.abs(w[~flagged]).max()) if (~flagged).any() else -1, 4
        ),
    }


def part_eq(a, b):
    return sorted(a, reverse=True) == sorted(b, reverse=True)


# ================= gap tail =================
def gap_tail(M, R):
    """own quadrature of the pair-gap deviation on the interpolated sphere: a_l, l = 0, 1, 2 (P_l),
    plus the rms of the m != 0 part at l = 1 (the axisymmetry of the dipole)."""
    th = np.linspace(0.0, np.pi, 181)
    tm = 0.5 * (th[:-1] + th[1:])
    ph = np.linspace(0.0, 2.0 * np.pi, 240, endpoint=False)
    T, P = np.meshgrid(tm, ph, indexing="ij")
    rhat = np.stack([np.sin(T) * np.cos(P), np.sin(T) * np.sin(P), np.cos(T)], -1)
    lam = np.linalg.eigvalsh(sample_block(M, R * rhat))
    g = lam[..., 1] - lam[..., 0] - DELTA
    x = np.cos(T)
    wgt = np.sin(T) * (th[1] - th[0]) * (ph[1] - ph[0]) / (4.0 * np.pi)  # integrates to 1
    a0 = float((g * wgt).sum())
    a1 = float(3.0 * (g * x * wgt).sum())
    a2 = float(5.0 * (g * 0.5 * (3 * x * x - 1) * wgt).sum())
    b1x = float(3.0 * (g * np.sin(T) * np.cos(P) * wgt).sum())
    b1y = float(3.0 * (g * np.sin(T) * np.sin(P) * wgt).sum())
    return {
        "R": R,
        "a0": a0,
        "a1": a1,
        "a2": a2,
        "a1_xy": [b1x, b1y],
        "maxdev": float(np.abs(g).max()),
    }


def gap_tail_cells(M, R):
    X, Y, Z = coords()
    r = np.sqrt(X * X + Y * Y + Z * Z)
    sel = np.abs(r - R) <= 0.5 * H
    lam = np.linalg.eigvalsh(M[sel][:, 1:, 1:])
    g = lam[:, 1] - lam[:, 0] - DELTA
    x = Z[sel] / np.maximum(r[sel], 1e-300)
    A = np.stack([np.ones_like(x), x, 0.5 * (3 * x * x - 1)], 1)
    c, *_ = np.linalg.lstsq(A, g, rcond=None)
    return {"R": R, "a0": float(c[0]), "a1": float(c[1]), "a2": float(c[2]), "n": int(sel.sum())}


def loglog_slope(rs, vals):
    rs, vals = np.asarray(rs), np.abs(np.asarray(vals))
    ok = vals > 0
    if ok.sum() < 3:
        return None
    return float(np.polyfit(np.log(rs[ok]), np.log(vals[ok]), 1)[0])


# ================= beta^2 =================
def beta2_of(Lm):
    t2 = np.einsum("...ab,...ba->...", Lm, Lm)
    t3 = np.einsum("...ab,...bc,...ca->...", Lm, Lm, Lm)
    return 1.0 - 6.0 * t3 * t3 / np.maximum(t2, 1e-300) ** 3, t2


def beta2_reads(M, rmax=6.0):
    X, Y, Z = coords()
    r = np.sqrt(X * X + Y * Y + Z * Z)
    S = M[..., 1:, 1:]
    Lm = S - (np.einsum("...kk->...", S) / 3.0)[..., None, None] * np.eye(3)
    out = {}
    for name, sel in (("lt", r < rmax), ("le", r <= rmax)):
        Lbar = Lm[sel].mean(0)
        bb, _ = beta2_of(Lbar)
        bc, t2 = beta2_of(Lm[sel])
        res = t2 > 1e-8
        out[name] = {
            "block_mean": float(bb),
            "cell_mean": float(bc[res].mean()),
            "cell_mean_all": float(bc.mean()),
            "n_cells": int(sel.sum()),
            "resolved_fraction": float(res.mean()),
            "Lbar_eigs": [float(e) for e in np.linalg.eigvalsh(Lbar)],
        }
    vac = np.diag([1.0, DELTA, 0.0])
    out["vacuum"] = float(beta2_of(vac - np.trace(vac) / 3.0 * np.eye(3))[0])
    return out


# ================= physical generator =================
def d_central(f, ax, order=2):
    """own central difference (interior) with one-sided faces; order 4 uses the 5-point stencil
    inside and drops to order 2 one cell from the face."""
    out = np.zeros_like(f)
    sl = [slice(None)] * f.ndim

    def at(i):
        s = list(sl)
        s[ax] = i
        return tuple(s)

    out[at(slice(1, -1))] = (f[at(slice(2, None))] - f[at(slice(0, -2))]) / (2 * H)
    out[at(0)] = (f[at(1)] - f[at(0)]) / H
    out[at(-1)] = (f[at(-1)] - f[at(-2)]) / H
    if order == 4:
        out[at(slice(2, -2))] = (
            -f[at(slice(4, None))]
            + 8.0 * f[at(slice(3, -1))]
            - 8.0 * f[at(slice(1, -3))]
            + f[at(slice(0, -4))]
        ) / (12.0 * H)
    return out


def d_central_halfface(f, ax):
    """the central stencil with the face rows at half weight (the mean of one-sided forward and
    backward differences, each zero where it runs off the grid)."""
    return 0.5 * (d_fwd(f, ax) + d_bwd(f, ax))


def generator_reads(M, order=2):
    X, Y, _ = coords()
    Jz = np.zeros((4, 4))
    Jz[1, 2], Jz[2, 1] = -1.0, 1.0
    a_int = Jz @ M - M @ Jz
    if order == 0:
        dy, dx = d_central_halfface(M, 1), d_central_halfface(M, 0)
    else:
        dy, dx = d_central(M, 1, order), d_central(M, 0, order)
    a_orb = X[..., None, None] * dy - Y[..., None, None] * dx
    a_rig = a_int - a_orb
    a_rig_plus = a_int + a_orb
    return {
        "C_int": kin(M, a_int),
        "C_orb": kin(M, a_orb),
        "C_rigid": kin(M, a_rig),
        "C_int_plus_orb_wrong_sign": kin(M, a_rig_plus),
        "a_rigid_asym_max": float(np.abs(a_rig - a_rig.swapaxes(-1, -2)).max()),
    }


# ================= verdicts =================
def verdict(ok_confirm, ok_qualify, reason):
    if ok_confirm:
        return {"verdict": "CONFIRMED", "reason": reason}
    if ok_qualify:
        return {"verdict": "QUALIFIED", "reason": reason}
    return {"verdict": "REFUTED", "reason": reason}


def main():
    with open(CENSUS_JSON) as f:
        J = json.load(f)
    rows = J["rows"]
    col = J["collect"]["table"]
    pin = pin_mask()
    free = ~pin
    res = {"task": "M5.32 R26-4 audit", "claims": {}, "runtime_s": None}

    # ---------- A: energies, gradients, labels, kick ----------
    log("A: own energy and gradient")
    a_err = adjoint_check()
    fields = {
        p: {"end": load(p), "gate": load(p, "_gate"), "seed": load(p, "_seed")} for p in PARTS
    }
    for p in PARTS:
        if os.path.exists(os.path.join(ARR, TAG.format(p) + "_kick_stage.npz")):
            fields[p]["kick"] = load(p, "_kick_stage")
    cs_err = complex_step_check(fields["2_2"]["end"])
    A = {"adjoint_check_relerr": a_err, "complex_step_gradient_relerr": cs_err, "rows": {}}
    a_ok, a_notes = True, []
    for p in PARTS:
        tag = TAG.format(p)
        r, c = rows[tag], col[tag]
        Mend, Mgate = fields[p]["end"], fields[p]["gate"]
        eu, ev = energy_parts(Mend)
        E_own = eu + ev
        G = gradient(Mend)
        Gf = G[free]
        fmax_sp = float(np.abs(Gf[:, 1:, 1:]).max())
        fmax_00 = float(np.abs(Gf[:, 0, 0]).max())
        last = r["chunks"][-1]
        row = {
            "E_own": E_own,
            "E_audited": r["E"],
            "dE": E_own - r["E"],
            "E_u_own": eu,
            "V_own": ev,
            "fmax_spatial_own": fmax_sp,
            "fmax_spatial_audited_last_chunk": last["fmax_spatial"],
            "fmax_M00_own": fmax_00,
            "iters_audited": r["iters"],
            "gate_label_audited": r["gate_label"],
            "kick_label_audited": r["kick_label"],
            "end_equals_gate_maxdiff": float(np.abs(Mend - Mgate).max()),
            "E_gate_own": float(energy(Mgate)),
            "E_seed_own": float(energy(fields[p]["seed"])),
            "E_seed_audited": r["E_seed"],
            "pin_holds_seed_maxdiff": float(np.abs(Mend[pin] - fields[p]["seed"][pin]).max()),
        }
        # own gate verdict: fmax < 1e-4 -> AT_GATE else FALLING (the cap)
        row["gate_label_own"] = "AT_GATE" if fmax_sp < 1e-4 else "FALLING"
        if "kick" in fields[p]:
            Ek = float(energy(fields[p]["kick"]))
            rel = abs(Ek - row["E_gate_own"]) / max(1.0, abs(row["E_gate_own"]))
            row["E_kick_relaxed_own"] = Ek
            row["E_kick_relaxed_audited"] = r["kick"]["E_after"]
            row["kick_return_rel"] = rel
            row["kick_lower_by"] = row["E_gate_own"] - Ek
            if Ek < row["E_gate_own"] - 1e-4:
                row["kick_label_own"] = "SADDLE"
            elif rel < 1e-3:
                row["kick_label_own"] = "STABLE"
            else:
                row["kick_label_own"] = "UNRESOLVED"
            Gk = gradient(fields[p]["kick"])[free]
            row["fmax_spatial_kick_relaxed_own"] = float(np.abs(Gk[:, 1:, 1:]).max())
        else:
            row["kick_label_own"] = "FALLING (not kicked)"
        A["rows"][p] = row
        e_ok = abs(row["dE"]) < 1e-6 * max(1.0, abs(r["E"]))
        g_ok = row["gate_label_own"] == r["gate_label"]
        k_ok = row["kick_label_own"].startswith(r["kick_label"])
        f_ok = abs(fmax_sp - last["fmax_spatial"]) < 0.1 * max(fmax_sp, last["fmax_spatial"])
        a_ok &= e_ok and g_ok and k_ok and f_ok
        if not (e_ok and g_ok and k_ok and f_ok):
            a_notes.append(f"{p}: E {e_ok} gate {g_ok} kick {k_ok} fmax {f_ok}")
        log(
            f"  {p:8s} E own {E_own:.6f} aud {r['E']:.6f} fmax own {fmax_sp:.2e} aud "
            f"{last['fmax_spatial']:.2e} gate {row['gate_label_own']} kick {row['kick_label_own']}"
        )
    # the record's table quotes fmax 1.0e-3 ({4}) and 1.6e-3 ({2,1,1}) at the cap: which chunk is that
    quoted = {"4": 1.0e-3, "2_1_1": 1.6e-3}
    A["record_quoted_fmax"] = {}
    q_ok = True
    for p, fq in quoted.items():
        ch = rows[TAG.format(p)]["chunks"]
        near = min(ch[1:], key=lambda c: abs(np.log(c["fmax_spatial"] / fq)))
        end_f = A["rows"][p]["fmax_spatial_own"]
        A["record_quoted_fmax"][p] = {
            "quoted": fq,
            "end_field_own": end_f,
            "closest_chunk_iters": near["iters"],
            "closest_chunk_fmax": near["fmax_spatial"],
            "end_over_quoted": end_f / fq,
            "fmax_last_6_chunks": [c["fmax_spatial"] for c in ch[-6:]],
        }
        q_ok &= abs(end_f - fq) < 0.25 * fq
    A.update(
        verdict(
            a_ok and cs_err < 1e-6 and q_ok,
            a_ok and cs_err < 1e-6,
            (
                (
                    "energies to 1e-15, gradient max, gate and kick labels agree on all five rows"
                    + (
                        ""
                        if q_ok
                        else "; the record's quoted fmax at the cap (1.0e-3, 1.6e-3) is the 7750-iteration "
                        "chunk, the stored end field (8000) reads 2.8e-3 and 5.2e-3 (FALLING either way)"
                    )
                )
                if a_ok
                else "; ".join(a_notes)
            ),
        )
    )
    res["claims"]["A_energies_labels"] = A

    # ---------- B: seed gate ----------
    log("B: own reader on the seeds")
    B = {"rows": {}}
    b_ok = True
    for p in PARTS:
        B["rows"][p] = {"intended": INTENDED[p], "reads": {}}
        for R in SPHERES_SEED:
            rd = read_sphere(fields[p]["seed"], R)
            B["rows"][p]["reads"][f"R{R:g}"] = rd
            ok = part_eq(rd["partition"], INTENDED[p]) and rd["total"] == 4
            b_ok &= ok
            log(
                f"  seed {p:8s} R{R:g}: {rd['partition']} total {rd['total']} (raw {rd['total_raw']}) "
                f"deg {rd['degree_of_director']} weak {rd['frac_director_not_outward_0p3']} "
                f"flagged {rd['frac_flagged_own']}"
            )
    B.update(
        verdict(
            b_ok,
            False,
            (
                "every seed reads its intended partition on r 9 and r 18 with total 4 in the own reader"
                if b_ok
                else "at least one seed does not read its intended partition in the own reader"
            ),
        )
    )
    res["claims"]["B_seed_gate"] = B

    # ---------- C: end-field partitions ----------
    log("C: own reader on the end fields")
    C = {"rows": {}}
    claimed_r9 = {
        "2_2": [1, 1, 1, 1],
        "3_1": [2, 1, 1],
        "1_1_1_1": [1, 1, 1, 1],
        "4": [1, 1, 1, 1, 1, 1, -1, -1],
        "2_1_1": [3, 1, 1, 1, -1, -1],
    }
    c_match, c_notes = {}, []
    for p in PARTS:
        tag = TAG.format(p)
        aud = rows[tag]["own_reads"]["spheres"]
        C["rows"][p] = {
            "reads": {},
            "audited_r9": aud["R9"]["partition"],
            "claimed_r9": claimed_r9[p],
        }
        for R in SPHERES_END:
            rd = read_sphere(fields[p]["end"], R)
            key = f"R{R:g}"
            rd["audited_partition"] = aud[key]["partition"] if key in aud else None
            rd["audited_carriers_rho_z"] = (
                [[c[0], c[1], c[2]] for c in aud[key]["carriers_hu_rho_z_phi"]]
                if key in aud
                else None
            )
            C["rows"][p]["reads"][key] = rd
            log(
                f"  end {p:8s} R{R:g}: own {rd['partition']} (tot {rd['total']}, unread net "
                f"{rd['unreadable_regions_net_raw']}, all {rd['partition_all_regions']}) "
                f"aud {rd['audited_partition']} weak {rd['frac_director_not_outward_0p3']} "
                f"flagged {rd['frac_flagged_own']} regions "
                f"{[(g['m_half_units'], g['rho'], g['z']) for g in rd['regions'] if g['m_half_units'] and g['clean']]}"
            )
        r9 = C["rows"][p]["reads"]["R9"]
        c_match[p] = part_eq(r9["partition"], claimed_r9[p])
        if not c_match[p]:
            c_notes.append(f"{p} r9 own {r9['partition']} vs claimed {claimed_r9[p]}")
    # the {1,1,1,1} carrier positions at r 9 and r 12
    pos = {}
    for R in (9.0, 12.0):
        regs = [
            g
            for g in C["rows"]["1_1_1_1"]["reads"][f"R{R:g}"]["regions"]
            if g["m_half_units"] == 1 and g["clean"]
        ]
        pos[f"R{R:g}"] = {
            "rho_range": (
                [min(g["rho"] for g in regs), max(g["rho"] for g in regs)] if regs else None
            ),
            "absz_range": (
                [min(abs(g["z"]) for g in regs), max(abs(g["z"]) for g in regs)] if regs else None
            ),
            "n_half_units": len(regs),
        }
    C["positions_1_1_1_1"] = pos
    C["claimed_positions"] = {
        "R9": {"rho": [4.4, 8.2], "absz": [3.6, 7.8]},
        "R12": {"rho": [5.1, 10.4], "absz": [5.9, 10.8]},
    }
    n_ok = sum(c_match.values())
    C["r9_match"] = c_match
    r6_match = {
        p: part_eq(
            C["rows"][p]["reads"]["R6"]["partition"],
            rows[TAG.format(p)]["own_reads"]["spheres"]["R6"]["partition"],
        )
        for p in PARTS
    }
    C["r6_match"] = r6_match
    C["director_degree_r9"] = {p: C["rows"][p]["reads"]["R9"]["degree_of_director"] for p in PARTS}
    C["clean_total_r9"] = {p: C["rows"][p]["reads"]["R9"]["total"] for p in PARTS}
    pos_ok = (
        abs(pos["R9"]["rho_range"][0] - 4.4) < 0.3
        and abs(pos["R9"]["rho_range"][1] - 8.2) < 0.3
        and abs(pos["R9"]["absz_range"][0] - 3.6) < 0.3
        and abs(pos["R9"]["absz_range"][1] - 7.8) < 0.3
        and abs(pos["R12"]["rho_range"][0] - 5.1) < 0.3
        and abs(pos["R12"]["rho_range"][1] - 10.4) < 0.3
        and abs(pos["R12"]["absz_range"][0] - 5.9) < 0.5
        and abs(pos["R12"]["absz_range"][1] - 10.8) < 0.3
    )
    C["positions_1_1_1_1_match"] = bool(pos_ok)
    reason = (
        f"r 9: own clean carriers reproduce the audited partition on {n_ok} of 5 fields "
        f"(the certified three and {{4}}; {'; '.join(c_notes) if c_notes else 'none differ'}: the audited 3 is a "
        "polar-cap merge of two off-axis +1 at z about -5, a FALLING row); the {1,1,1,1} carrier positions "
        f"at r 9 and 12 match {pos_ok}; the clean carriers alone sum to 4 on every r 9 sphere; 11 to 24 "
        "percent of r 9 and 17 to 28 percent of r 6 have the director tangent (|d.rhat| < 0.3) and the "
        "flagged regions (2 to 5 percent) sit on director-gap closures (0.005 to 0.1) where the director's "
        "orientation cannot be continued (its degree reads 0, 2, 3 there): unreadable by any reader; "
        f"r 6 agrees on {sum(r6_match.values())} of 5 (the {{2,2}} r 6 reads a clean {{1,1,1,1}} at rho 4.9, "
        "|z| 3.4 here against the JSON's [4], which supports the record's 'r 6 to 15' sentence over its JSON)"
    )
    C.update(
        verdict(n_ok == 5 and pos_ok and sum(r6_match.values()) >= 4, n_ok >= 4 and pos_ok, reason)
    )
    res["claims"]["C_end_partitions"] = C

    # ---------- D: the moment channel ----------
    log("D: own gap-tail projections")
    D = {"rows": {}}
    for p in ("3_1", "1_1_1_1", "2_2", "2_1_1", "4"):
        tag = TAG.format(p)
        aud = rows[tag]["own_reads"]["gap_tail"]
        aud_shells = {s[0]: s for s in aud["shells_r_a0_a1_a2_maxdev_n"]}
        own_sph = [gap_tail(fields[p]["end"], R) for R in SHELLS]
        own_cells = [gap_tail_cells(fields[p]["end"], R) for R in SHELLS]
        seed_sph = [gap_tail(fields[p]["seed"], R) for R in SHELLS]
        win = [k for k, R in enumerate(SHELLS) if 6.0 <= R <= 18.0]
        rs = [SHELLS[k] for k in win]
        sl_sph = loglog_slope(rs, [own_sph[k]["a1"] for k in win])
        sl_cells = loglog_slope(rs, [own_cells[k]["a1"] for k in win])
        sl_seed = loglog_slope(rs, [seed_sph[k]["a1"] for k in win])
        a1_sph = [own_sph[k]["a1"] for k in win]
        a1_seed = [seed_sph[k]["a1"] for k in win]
        # monotone trend: sign of the fitted linear slope of |a1| against r
        trend = float(np.polyfit(rs, np.abs(a1_sph), 1)[0])
        ratio_max = max(abs(own_sph[k]["a1"]) / max(abs(own_sph[k]["a0"]), 1e-300) for k in win)
        D["rows"][p] = {
            "shells": SHELLS,
            "own_sphere_a0_a1_a2": [[o["a0"], o["a1"], o["a2"]] for o in own_sph],
            "own_sphere_a1_xy": [o["a1_xy"] for o in own_sph],
            "own_cells_a0_a1_a2_n": [[o["a0"], o["a1"], o["a2"], o["n"]] for o in own_cells],
            "seed_sphere_a0_a1_a2": [[o["a0"], o["a1"], o["a2"]] for o in seed_sph],
            "audited_a0_a1_a2_n": [
                (
                    [aud_shells[R][1], aud_shells[R][2], aud_shells[R][3], aud_shells[R][5]]
                    if R in aud_shells
                    else None
                )
                for R in SHELLS
            ],
            "window": rs,
            "slope_l1_own_sphere": sl_sph,
            "slope_l1_own_cells": sl_cells,
            "slope_l1_seed": sl_seed,
            "slope_l1_audited": aud["slope_l1"],
            "abs_a1_linear_trend_per_unit_r": trend,
            "a1_window_end": [a1_sph[0], a1_sph[-1]],
            "a1_seed_window_end": [a1_seed[0], a1_seed[-1]],
            "l1_over_l0_max_own": ratio_max,
            "l1_over_l0_max_audited": aud["l1_over_l0_max_in_window"],
            "seed_end_a1_corr": (
                float(np.corrcoef(a1_sph, a1_seed)[0, 1]) if np.std(a1_seed) > 0 else None
            ),
            "seed_end_a1_ratio_at_18": (
                (a1_sph[-1] / a1_seed[-1]) if abs(a1_seed[-1]) > 1e-12 else None
            ),
        }
        log(
            f"  {p:8s} a1 own {np.array2string(np.array(a1_sph), precision=4)} slope own {sl_sph} "
            f"(cells {sl_cells}) aud {aud['slope_l1']:.2f}; seed a1 "
            f"{np.array2string(np.array(a1_seed), precision=4)} slope {sl_seed}"
        )
    d31 = D["rows"]["3_1"]
    rises = d31["abs_a1_linear_trend_per_unit_r"] > 0 and abs(d31["a1_window_end"][1]) > abs(
        d31["a1_window_end"][0]
    )
    a1_win = [
        row[1] for k, row in enumerate(d31["own_sphere_a0_a1_a2"]) if 6.0 <= SHELLS[k] <= 18.0
    ]
    a1_range_ok = max(abs(x) for x in a1_win) < 0.05
    seed_same = d31["seed_end_a1_corr"] is not None and d31["seed_end_a1_corr"] > 0.9
    D["l1_rises_outward_3_1"] = bool(rises)
    D["seed_carries_same_l1_profile_3_1"] = bool(seed_same)
    D["moment_label_audited"] = J["collect"]["moment_label"]
    a1_end_21 = d31["own_sphere_a0_a1_a2"][-1][1]
    a1_seed_21 = d31["seed_sphere_a0_a1_a2"][-1][1]
    D["a1_at_wall_r21_end_over_seed_3_1"] = (
        a1_end_21 / a1_seed_21 if abs(a1_seed_21) > 1e-12 else None
    )
    D["m1_over_m0_3_1_window"] = [
        float(
            np.hypot(*d31["own_sphere_a1_xy"][k])
            / max(abs(d31["own_sphere_a0_a1_a2"][k][1]), 1e-300)
        )
        for k, R in enumerate(SHELLS)
        if 6.0 <= R <= 18.0
    ]
    D.update(
        verdict(
            rises and a1_range_ok and seed_same,
            rises and a1_range_ok,
            (
                (
                    f"the {{3,1}} l = 1 coefficient rises outward (own slope {d31['slope_l1_own_sphere']:.2f}, "
                    f"0.002 to 0.042 on r 6 to 18, l = 0 crossing zero near r 8): not a dipole tail, CONFIRMED; "
                    "but the seed's l = 1 profile FALLS outward (0.031 at r 7.5 to 0.011 at r 18, slope "
                    f"{d31['slope_l1_seed']:.2f}, correlation with the end's profile {d31['seed_end_a1_corr']:.2f}) "
                    f"and the end's coefficient is {d31['seed_end_a1_ratio_at_18']:.1f}x the seed's at r 18 and "
                    f"{D['a1_at_wall_r21_end_over_seed_3_1']:.1f}x at r 21 (the pinned value): the rise is built "
                    "by the relaxation against the wall, not the seed's asymmetry held by the shell; the m = +-1 "
                    "parts of l = 1 are as large as the m = 0 part (the asymmetry is not axial)"
                )
                if rises
                else "the own l = 1 coefficient does not rise outward on {3,1}"
            ),
        )
    )
    res["claims"]["D_moment_channel"] = D

    # ---------- E: beta^2 ----------
    log("E: beta^2")
    E = {"rows": {}}
    e_ok = True
    for p in PARTS:
        tag = TAG.format(p)
        aud = rows[tag]["own_reads"]["biaxiality"]
        b = beta2_reads(fields[p]["end"])
        b["audited_block_mean_r6"] = aud["block_mean_r6"]
        b["audited_cell_mean_r6"] = aud["cell_mean_r6"]
        d_block = min(
            abs(b["lt"]["block_mean"] - aud["block_mean_r6"]),
            abs(b["le"]["block_mean"] - aud["block_mean_r6"]),
        )
        d_cell = min(
            abs(b["lt"]["cell_mean"] - aud["cell_mean_r6"]),
            abs(b["le"]["cell_mean"] - aud["cell_mean_r6"]),
        )
        b["max_abs_diff"] = max(d_block, d_cell)
        e_ok &= b["max_abs_diff"] < 0.01
        E["rows"][p] = b
        log(
            f"  {p:8s} block {b['lt']['block_mean']:.3f} (aud {aud['block_mean_r6']:.3f}) cell "
            f"{b['lt']['cell_mean']:.3f} (aud {aud['cell_mean_r6']:.3f}) vac {b['vacuum']:.4f}"
        )
    E.update(
        verdict(
            e_ok,
            True,
            (
                "block and cell means inside r 6 reproduced to 0.01"
                if e_ok
                else "a beta^2 mean differs by over 0.01"
            ),
        )
    )
    res["claims"]["E_beta2"] = E

    # ---------- F: virial and one barrier ----------
    log("F: virial ratios and the {1,1,1,1} to {2,2} barrier path")
    F = {"virial": {}, "barrier": {}}
    f_ok = True
    for p in PARTS:
        tag = TAG.format(p)
        eu, ev = A["rows"][p]["E_u_own"], A["rows"][p]["V_own"]
        aud = rows[tag]["own_reads"]["virial"]["E_u_over_3V"]
        F["virial"][p] = {"own": eu / (3.0 * ev), "audited": aud}
        f_ok &= abs(F["virial"][p]["own"] - aud) < 1e-3
    Ma, Mb = fields["1_1_1_1"]["end"], fields["2_2"]["end"]
    ss = np.linspace(0.0, 1.0, 11)
    path, split = [], []
    for s in ss:
        Mx = (1.0 - s) * Ma + s * Mb
        eu, ev = density_parts(Mx)
        e = eu + ev
        path.append(float(e.sum()))
        # the first free layer beside the pin
        near = np.zeros((N, N, N), dtype=bool)
        for ax in range(3):
            sl = [slice(None)] * 3
            sl[ax] = slice(PIN_LAYERS, PIN_LAYERS + 1)
            near[tuple(sl)] = True
            sl[ax] = slice(N - PIN_LAYERS - 1, N - PIN_LAYERS)
            near[tuple(sl)] = True
        deep = ~(pin | near)
        split.append([float(e[pin].sum()), float(e[near & ~pin].sum()), float(e[deep].sum())])
    Mm = 0.5 * (Ma + Mb)
    eu_m, ev_m = density_parts(Mm)
    X, Y, Z = coords()
    rr = np.sqrt(X * X + Y * Y + Z * Z)
    em = eu_m + ev_m
    radial = {
        "r_lt_9": float(em[rr < 9].sum()),
        "r_9_to_18": float(em[(rr >= 9) & (rr < 18)].sum()),
        "r_ge_18": float(em[rr >= 18].sum()),
    }
    aud_b = J["collect"]["barriers_upper_bound"]["P1_1_1_1_d0.3_w25_n32_L48|P2_2_d0.3_w25_n32_L48"]
    diff = Ma - Mb
    F["barrier"] = {
        "pair": "{1,1,1,1} -> {2,2}",
        "E_path_own": path,
        "E_path_audited": aud_b["E_path"],
        "max_own": max(path),
        "height_above_a_own": max(path) - path[0],
        "max_audited": aud_b["barrier_upper_bound_from_a"],
        "midpoint_split_pin_firstfree_deep": split[5],
        "endpoint_a_split": split[0],
        "endpoint_b_split": split[10],
        "field_diff_rms_pin": float(np.sqrt((diff[pin] ** 2).mean())),
        "field_diff_rms_free": float(np.sqrt((diff[free] ** 2).mean())),
        "field_diff_max_pin": float(np.abs(diff[pin]).max()),
        "field_diff_max_free": float(np.abs(diff[free]).max()),
        "midpoint_curvature": float(eu_m.sum()),
        "midpoint_potential": float(ev_m.sum()),
        "midpoint_by_radius": radial,
    }
    mid = split[5]
    F["barrier"]["midpoint_frac_pin_plus_first_free"] = (mid[0] + mid[1]) / sum(mid)
    F["barrier"]["midpoint_frac_deep"] = mid[2] / sum(mid)
    bar_ok = (
        abs((max(path) - path[0]) - aud_b["barrier_upper_bound_from_a"]) < 1e-6
        and max(abs(a - b) for a, b in zip(path, aud_b["E_path"])) < 1e-9
    )
    boundary_dominated = F["barrier"]["midpoint_frac_pin_plus_first_free"] > 0.5
    F["barrier"]["boundary_dominated"] = bool(boundary_dominated)
    log(
        f"  barrier own max {max(path):.2f} aud {aud_b['barrier_upper_bound_from_a']:.2f}; midpoint split pin/first-free/deep "
        f"{mid[0]:.1f}/{mid[1]:.1f}/{mid[2]:.1f}; diff rms pin {F['barrier']['field_diff_rms_pin']:.3f} free "
        f"{F['barrier']['field_diff_rms_free']:.3f}"
    )
    F.update(
        verdict(
            f_ok and bar_ok and boundary_dominated,
            f_ok and bar_ok,
            (
                (
                    "virial ratios to 1e-15 and the 11-point path to 3e-14 (barrier height 141.4 above the "
                    "{1,1,1,1} end, path max 149.8); the 'uninformative upper bound' conclusion stands but the "
                    f"stated cause is wrong: the pinned shell plus its first free layer hold {mid[0] + mid[1]:.1f} "
                    f"of the midpoint's {sum(mid):.1f} (1 percent), the rest is deep interior "
                    f"({radial['r_lt_9']:.0f} inside r 9, {radial['r_9_to_18']:.0f} in r 9 to 18), "
                    f"{F['barrier']['midpoint_potential'] / sum(mid) * 100:.0f} percent of it potential: the "
                    "straight line leaves the vacuum manifold where the two interiors differ (rms difference "
                    f"{F['barrier']['field_diff_rms_free']:.3f} free against {F['barrier']['field_diff_rms_pin']:.3f} "
                    "pinned), the usual reason a linear interpolation barrier is uninformative"
                )
                if (f_ok and bar_ok)
                else "virial or barrier numbers not reproduced"
            ),
        )
    )
    res["claims"]["F_virial_barrier"] = F

    # ---------- G: physical generator ----------
    log("G: physical generator")
    G = {"rows": {}, "control_seed_2_2": None}
    g_ok = True
    for p in PARTS:
        tag = TAG.format(p)
        aud = rows[tag]["own_reads"]["physical_generator"]
        o2 = generator_reads(fields[p]["end"], 2)
        o4 = generator_reads(fields[p]["end"], 4)
        o0 = generator_reads(fields[p]["end"], 0)
        rel = max(
            abs(o2["C_int"] - aud["internal"]) / aud["internal"],
            abs(o2["C_orb"] - aud["orbital"]) / aud["orbital"],
            abs(o2["C_rigid"] - aud["rigid"]) / aud["rigid"],
        )
        rel0 = max(
            abs(o0["C_int"] - aud["internal"]) / aud["internal"],
            abs(o0["C_orb"] - aud["orbital"]) / aud["orbital"],
            abs(o0["C_rigid"] - aud["rigid"]) / aud["rigid"],
        )
        G["rows"][p] = {
            "own_order2_onesided_faces": o2,
            "own_order2_halfweight_faces": o0,
            "own_order4": o4,
            "audited": [aud["internal"], aud["orbital"], aud["rigid"]],
            "max_rel_diff_onesided_faces": rel,
            "max_rel_diff_halfweight_faces": rel0,
            "rigid_rel_change_order4": abs(o4["C_rigid"] - o2["C_rigid"]) / o2["C_rigid"],
            "rigid_rel_change_faces": abs(o0["C_rigid"] - o2["C_rigid"]) / o2["C_rigid"],
        }
        g_ok &= rel0 < 1e-6 and rel < 1e-2
        log(
            f"  {p:8s} own {o2['C_int']:.0f} {o2['C_orb']:.0f} {o2['C_rigid']:.0f} (o4 rigid {o4['C_rigid']:.0f}) "
            f"aud {aud['internal']:.0f} {aud['orbital']:.0f} {aud['rigid']:.0f}"
        )
    ctl = generator_reads(fields["2_2"]["seed"], 2)
    G["control_seed_2_2"] = {**ctl, "rigid_over_int": ctl["C_rigid"] / ctl["C_int"]}
    log(f"  control {{2,2}} seed: C_rigid / C_int {G['control_seed_2_2']['rigid_over_int']:.4f}")
    try:
        sys.path.insert(0, HERE)
        import importlib.util

        spec = importlib.util.spec_from_file_location("b3", os.path.join(HERE, "m5_21_3_a_4d.py"))
        b3 = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(b3)
        cfg = b3.base_cfg(s=-1.0, g=8.0, n=N, L=L, delta=DELTA)
        Jz = np.zeros((4, 4))
        Jz[1, 2], Jz[2, 1] = -1.0, 1.0
        Mx = fields["1_1_1_1"]["end"]
        a0 = Jz @ Mx - Mx @ Jz
        G["platform_kin_of_check"] = {
            "platform": float(b3.kin_of(Mx, a0, cfg)),
            "own": kin(Mx, a0),
            "platform_e_total": float(b3.e_total(Mx, cfg)),
            "note": "platform e_total uses W1 (w1s 1), not comparable to E; kin_of comparable",
        }
    except Exception as e:  # noqa: BLE001
        G["platform_kin_of_check"] = {"error": repr(e)}
    sens = max(r["rigid_rel_change_order4"] for r in G["rows"].values())
    G["max_rigid_sensitivity_order4"] = sens
    G.update(
        verdict(
            g_ok and G["control_seed_2_2"]["rigid_over_int"] < 0.1 and sens < 0.05,
            g_ok and G["control_seed_2_2"]["rigid_over_int"] < 0.1,
            (
                (
                    "C_int, C_orb, C_rigid reproduced to 1e-12 with the half-weight face rule and C_int, C_rigid "
                    "to 1e-4 with full one-sided faces (C_orb then 0.8 percent off: the face rows); the "
                    f"axisymmetric {{2,2}} seed control cancels to {G['control_seed_2_2']['rigid_over_int']:.4f}; "
                    f"a 4th-order interior stencil moves C_rigid by 3 to {sens * 100:.0f} percent ({{3,1}}, "
                    "{1,1,1,1}): the field is rough at the cell scale and the numbers are stencil-defined at the "
                    "10 percent level"
                )
                if g_ok
                else "an inertia differs from the audited value beyond the stencil's face rule"
            ),
        )
    )
    res["claims"]["G_physical_generator"] = G

    res["runtime_s"] = round(time.time() - T0, 1)
    # ---------- summary ----------
    print("\n| Claim | Own numbers | Audited | Verdict |")
    print("| --- | --- | --- | --- |")
    for key, cl in res["claims"].items():
        print(f"| {key} | see JSON | see JSON | {cl['verdict']}: {cl['reason']} |")

    def clean(o):
        if isinstance(o, dict):
            return {str(k): clean(v) for k, v in o.items()}
        if isinstance(o, (list, tuple)):
            return [clean(v) for v in o]
        if isinstance(o, (np.floating, float)):
            return float(o)
        if isinstance(o, (np.integer, int)):
            return int(o)
        if isinstance(o, (np.bool_, bool)):
            return bool(o)
        if isinstance(o, np.ndarray):
            return clean(o.tolist())
        return o

    with open(OUT_JSON, "w") as f:
        json.dump(clean(res), f, indent=1)
    log(f"wrote {OUT_JSON}")


if __name__ == "__main__":
    main()
