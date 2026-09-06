"""M5.32 R15-P-iv ADVERSARIAL AUDIT: the fixed-J descent on the projector object L_P.

Independent re-implementation from DEFINITIONS only (the producer's scripts and
logs were not opened).  Consumes the certified stack (m5_21_3_a_4d.py: base_cfg,
coords, pin_shell, e_parts (E_u only; V4^dd is rebuilt here with the DEGENERATE
targets), a_fields, comm_eta, inner_eta, kin_of) and B8.G3.  Everything else
(the projector P23, K_P^23, SPLIT, a0_local, densities, E_J, finite differences,
radial profiles, mutations) is written here.

Definitions audited:
  eta = diag(-1,1,1,1), N = M eta, vacuum diag(g,1,delta,delta), g = 8, delta = 0.3
  E_u      = e_parts(M)[0]
  V4^dd    = W1 h^3 sum_cells sum_p (tr N^p - C_p)^2, C_p = (-g)^p + 1 + 2 delta^p
  SPLIT    = h^3 sum (lambda_2 - lambda_3)^2 (middle eigenvalues of N by real part)
  K_P^23   = (1/2) h^3 sum_i tr(Om_i^T eta Om_i eta), Om_i = P23 A_i eta P23,
             P23 = I - P_g - P_1 (Lagrange projectors of the isolated eigenvalues)
  E_stat   = E_u + V4^dd + mu SPLIT + c_P K_P^23,  mu = 1e-2, c_P = 1
  a0       = J M - M J, J = rotation about the leading spatial eigenvector
  kin_I1   = kin_of(M, a0);  kin_KP23 = (1/2) h^3 tr(Om_0^T eta Om_0 eta), Om_0 = P23 a0 eta P23
  E_J      = E_stat + J^2 / (4 (kin_I1 + c_P kin_KP23)),  J = 200

Run: OMP_NUM_THREADS=2 python3 m5_32_r15_p4_audit.py [quick]
Out: ../data/m5_32_r15_p4_audit.json
"""
from __future__ import annotations

import sys
ARGV = list(sys.argv[1:])          # captured BEFORE any import
import importlib.util
import json
import os
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.dirname(HERE)
DATA = os.path.join(RES, "data")
CK = os.path.join(RES, "checkpoints", "m5_32_r15")
T0 = time.time()
QUICK = "quick" in ARGV


def log(m):
    print(f"[{time.time() - T0:7.1f}s] {m}", flush=True)


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


INS4 = _load("m5_21_3_a_4d", "m5_21_3_a_4d.py")
B8 = _load("m5_21_8_b_lattice", "m5_21_8_b_lattice.py")
ETA = INS4.ETA
W1 = INS4.W1
G, DELTA, MU, CP, JFIX = 8.0, 0.3, 1e-2, 1.0, 200.0
I4 = np.eye(4)

SEED = os.path.join(CK, "m_hedgehog", "relax_n32_L48_mu0.01_cP1.npy")
END = os.path.join(CK, "p4_fixedj", "fixedJ_n32_L48_J200.npy")
END48 = os.path.join(CK, "p4_fixedj", "fixedJ_n48_L72_J200.npy")


# ------------------------------------------------------------------ own instruments
def n_eigs(M):
    """eigenvalues of N = M eta per cell, sorted by real part (4 columns)."""
    N = M @ ETA
    w = np.linalg.eigvals(N)
    idx = np.argsort(w.real, axis=-1)
    return np.take_along_axis(w, idx, axis=-1), N


def projector23(M):
    """P23 = I - P_g - P_1 with Lagrange projectors of the isolated eigenvalues
    (smallest real part = g-branch, largest = the 1-branch)."""
    w, N = n_eigs(M)
    wr = w.real
    lg, l2, l3, l1 = wr[..., 0], wr[..., 1], wr[..., 2], wr[..., 3]

    def shifted(lam):
        return N - lam[..., None, None] * I4

    Pg = shifted(l1) @ shifted(l2) @ shifted(l3)
    Pg = Pg / ((lg - l1) * (lg - l2) * (lg - l3))[..., None, None]
    P1 = shifted(lg) @ shifted(l2) @ shifted(l3)
    P1 = P1 / ((l1 - lg) * (l1 - l2) * (l1 - l3))[..., None, None]
    P23 = I4 - Pg - P1
    return P23, w, Pg, P1


def split_density(w, h3):
    """h^3 (lambda_2 - lambda_3)^2 per cell, from the sorted N spectrum."""
    s = w[..., 2].real - w[..., 1].real
    return h3 * s * s, s


def v4dd(M, cfg):
    h3 = cfg["h"] ** 3
    g, d = cfg["g"], cfg["delta"]
    Me = M @ ETA
    P = Me
    vd = 0.0
    for p in range(1, 5):
        if p > 1:
            P = P @ Me
        Cp = (-g) ** p + 1.0 + 2.0 * d ** p
        vd = vd + (np.einsum("...kk->...", P) - Cp) ** 2
    return W1 * h3 * np.sum(vd), W1 * h3 * vd


def kp23_density(M, cfg, P23):
    """(1/2) h^3 sum_i tr(Om_i^T eta Om_i eta), Om_i = P23 A_i eta P23, branch-weighted."""
    h3 = cfg["h"] ** 3
    dens = np.zeros(M.shape[:3])
    for br, (A, wt) in INS4.a_fields(M, cfg).items():
        for i in range(3):
            Om = P23 @ A[i] @ ETA @ P23
            dens += wt * 0.5 * INS4.inner_eta(Om, Om)
    return h3 * dens


def a0_local_own(M):
    """J M - M J, J the rotation generator about the leading spatial eigenvector."""
    w, V = np.linalg.eigh(M[..., 1:, 1:])
    n1 = V[..., :, -1]
    J = np.zeros(M.shape)
    J[..., 1, 2], J[..., 2, 1] = -n1[..., 2], n1[..., 2]
    J[..., 1, 3], J[..., 3, 1] = n1[..., 1], -n1[..., 1]
    J[..., 2, 3], J[..., 3, 2] = -n1[..., 0], n1[..., 0]
    return J @ M - M @ J


def kin_i1_density(M, a0, cfg):
    h3 = cfg["h"] ** 3
    dens = np.zeros(M.shape[:3])
    for br, (A, wt) in INS4.a_fields(M, cfg).items():
        for i in range(3):
            F = INS4.comm_eta(a0, A[i])
            dens += wt * 4.0 * INS4.inner_eta(F, F)
    return h3 * dens


def kin_kp23_density(M, a0, cfg, P23):
    h3 = cfg["h"] ** 3
    Om = P23 @ a0 @ ETA @ P23
    return h3 * 0.5 * INS4.inner_eta(Om, Om)


def reads(M, cfg, want_dens=False):
    """all reads of one field; E_J with a0 refreshed."""
    h3 = cfg["h"] ** 3
    e_u, _ = INS4.e_parts(M, cfg)
    ev4, v4d = v4dd(M, cfg)
    P23, w, Pg, P1 = projector23(M)
    sd, s = split_density(w, h3)
    split = float(np.sum(sd))
    kpd = kp23_density(M, cfg, P23)
    kp = float(np.sum(kpd))
    a0 = a0_local_own(M)
    ki1d = kin_i1_density(M, a0, cfg)
    ki1 = float(np.sum(ki1d))
    kkpd = kin_kp23_density(M, a0, cfg, P23)
    kkp = float(np.sum(kkpd))
    kin_tot = ki1 + CP * kkp
    e_stat = float(e_u + ev4 + MU * split + CP * kp)
    e_fix = JFIX ** 2 / (4.0 * kin_tot)
    out = {"E_u": float(e_u), "V4dd": float(ev4), "SPLIT": split, "muSPLIT": MU * split,
           "KP23": kp, "E_stat": e_stat, "kin_I1": ki1, "kin_KP23": kkp,
           "kin_tot": kin_tot, "E_fixJ": float(e_fix), "E_J": e_stat + float(e_fix),
           "omega": JFIX / (2.0 * kin_tot), "max_split": float(np.max(np.abs(s))),
           "max_imag_eig": float(np.max(np.abs(w.imag)))}
    if want_dens:
        out["_dens"] = {"split": sd, "s": s, "kp23": kpd, "kin_I1": ki1d, "kin_KP23": kkpd,
                        "v4": v4d, "a0": a0, "P23": P23}
    return out


def e_j(M, cfg):
    return reads(M, cfg)["E_J"]


def radial(dens, r, h, L):
    """shell profile: mean density per cell and shell sum, bins of width h."""
    edges = np.arange(0.0, r.max() + h, h)
    idx = np.digitize(r, edges) - 1
    prof = []
    for k in range(len(edges) - 1):
        m = idx == k
        if m.any():
            prof.append({"r_lo": float(edges[k]), "r_hi": float(edges[k + 1]),
                         "cells": int(m.sum()), "sum": float(dens[m].sum()),
                         "mean": float(dens[m].mean()), "max": float(dens[m].max())})
    return prof


def frac_inside(dens, r, rc):
    tot = dens.sum()
    ins = dens[r < rc].sum()
    return float(ins / tot) if tot != 0 else float("nan"), float(np.abs(dens)[r < rc].sum() / np.abs(dens).sum())


def split_rescale(M, f, free):
    """in each free cell's spatial eigenframe, scale the (2,3) eigenvalue split by f
    (the mean of the two and the leading eigenvalue untouched)."""
    out = M.copy()
    sub = M[free][..., 1:, 1:]
    w, V = np.linalg.eigh(sub)
    m = 0.5 * (w[..., 0] + w[..., 1])
    d = 0.5 * (w[..., 1] - w[..., 0])
    w = w.copy()
    w[..., 0] = m - f * d
    w[..., 1] = m + f * d
    sub2 = np.einsum("...ik,...k,...jk->...ij", V, w, V)
    blk = M[free].copy()
    blk[..., 1:, 1:] = sub2
    out[free] = blk
    return out


def rand_dir(rng, shape, free):
    D = rng.standard_normal(shape)
    D = 0.5 * (D + D.swapaxes(-1, -2))
    D = D * free[..., None, None]
    return D / np.sqrt(np.sum(D * D))


def dderiv(M, D, cfg, eps):
    """central first and second derivative of E_J along unit D at step eps."""
    ep = e_j(M + eps * D, cfg)
    em = e_j(M - eps * D, cfg)
    e0 = e_j(M, cfg)
    return (ep - em) / (2 * eps), (ep - 2 * e0 + em) / (eps * eps), e0


# ------------------------------------------------------------------ the crossing section
def crossing_cells(M, gap_tol=0.01):
    """cells where the leading and middle spatial eigenvalues are within gap_tol."""
    w = np.linalg.eigvalsh(M[..., 1:, 1:])
    return (w[..., 2] - w[..., 1]) < gap_tol, w


def crossing_dir(M, mask):
    """per crossing cell: v_mid v_mid^T - v_lead v_lead^T in the spatial block
    (raises the middle eigenvalue, lowers the leading one: eps > 0 crosses)."""
    D = np.zeros_like(M)
    w, V = np.linalg.eigh(M[mask][..., 1:, 1:])
    vm, vl = V[..., :, 1], V[..., :, 2]
    blk = np.zeros((mask.sum(), 4, 4))
    blk[..., 1:, 1:] = np.einsum("...i,...j->...ij", vm, vm) - np.einsum("...i,...j->...ij", vl, vl)
    D[mask] = blk
    return D


def e_j_frozen(M, a0, cfg):
    """E_J with a0 held FIXED (the descent's gradient object), own E_stat."""
    R = reads_static(M, cfg)
    P23, w, _, _ = projector23(M)
    ki1 = float(INS4.kin_of(M, a0, cfg))
    kkp = float(np.sum(kin_kp23_density(M, a0, cfg, P23)))
    kt = ki1 + CP * kkp
    return R["E_stat"] + JFIX ** 2 / (4.0 * kt), kt


def reads_static(M, cfg):
    h3 = cfg["h"] ** 3
    e_u, _ = INS4.e_parts(M, cfg)
    ev4, _ = v4dd(M, cfg)
    P23, w, _, _ = projector23(M)
    sd, s = split_density(w, h3)
    kp = float(np.sum(kp23_density(M, cfg, P23)))
    return {"E_stat": float(e_u + ev4 + MU * np.sum(sd) + CP * kp), "E_u": float(e_u),
            "V4dd": float(ev4), "SPLIT": float(np.sum(sd)), "KP23": kp}


def closest_pair_reads(M, cfg):
    """ALTERNATIVE reading of the object: the '(2,3) plane' = the CLOSEST pair of the
    three spatial eigenvalues (the degenerate pair), the clock = rotation about the
    ISOLATED eigenvector.  Label-free; coincides with the ordered reading wherever
    the isolated eigenvalue is the largest."""
    h3 = cfg["h"] ** 3
    w, V = np.linalg.eigh(M[..., 1:, 1:])
    gaps = np.stack([w[..., 1] - w[..., 0], w[..., 2] - w[..., 1]], axis=-1)
    iso = np.where(gaps[..., 0] < gaps[..., 1], 2, 0)          # index of the isolated eigenvalue
    n_iso = np.take_along_axis(V, iso[..., None, None].repeat(3, -2), axis=-1)[..., 0]
    pair_split = np.where(iso == 2, gaps[..., 0], gaps[..., 1])
    # clock about the isolated eigenvector
    J = np.zeros(M.shape)
    J[..., 1, 2], J[..., 2, 1] = -n_iso[..., 2], n_iso[..., 2]
    J[..., 1, 3], J[..., 3, 1] = n_iso[..., 1], -n_iso[..., 1]
    J[..., 2, 3], J[..., 3, 2] = -n_iso[..., 0], n_iso[..., 0]
    a0 = J @ M - M @ J
    # projector on the closest pair (spatial block, M_0i = 0 verified on both fields)
    Pp = np.zeros(M.shape)
    Pp[..., 1:, 1:] = np.eye(3) - np.einsum("...i,...j->...ij", n_iso, n_iso)
    kpd = kp23_density(M, cfg, Pp)
    ki1 = float(INS4.kin_of(M, a0, cfg))
    kkp = float(np.sum(kin_kp23_density(M, a0, cfg, Pp)))
    split = float(h3 * np.sum(pair_split ** 2))
    e_u, _ = INS4.e_parts(M, cfg)
    ev4, _ = v4dd(M, cfg)
    e_stat = float(e_u + ev4 + MU * split + CP * np.sum(kpd))
    kt = ki1 + CP * kkp
    return {"SPLIT_pair": split, "KP23_pair": float(np.sum(kpd)), "E_stat_pair": e_stat,
            "kin_I1_iso": ki1, "kin_KP23_pair": kkp, "kin_tot": kt,
            "E_J_pair": e_stat + JFIX ** 2 / (4.0 * kt), "max_pair_split": float(pair_split.max()),
            "cells_isolated_is_lowest": int((iso == 0).sum())}


def crossing_section(Ms, Me, cfg, free, r, rng, res):
    n, h = cfg["n"], cfg["h"]
    mask, w = crossing_cells(Me)
    sec = {"n_crossing_cells": int(mask.sum()),
           "crossing_cells": [{"idx": [int(v) for v in ij], "r": float(r[ij]),
                               "w": [float(v) for v in w[ij]], "gap12": float(w[ij][2] - w[ij][1])}
                              for ij in zip(*np.where(mask))]}
    a0e = a0_local_own(Me)
    # (1) the controlled crossing direction, refreshed vs frozen a0
    D = crossing_dir(Me, mask)
    rows = []
    for eps in [-1e-2, -3e-3, -1e-3, -3e-4, -1e-4, -1e-5, 0.0, 1e-6, 1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2]:
        Mx = Me + eps * D
        R = reads(Mx, cfg)
        ef, ktf = e_j_frozen(Mx, a0e, cfg)
        wx = np.linalg.eigvalsh(Mx[mask][..., 1:, 1:])
        rows.append({"eps": eps, "E_J_refreshed": R["E_J"], "E_J_frozen_a0": ef, "E_stat": R["E_stat"],
                     "kin_tot_refreshed": R["kin_tot"], "kin_tot_frozen": ktf, "kin_I1": R["kin_I1"],
                     "n_cells_crossed": int((wx[..., 2] - wx[..., 1] < 0).sum()) if False else
                     int(np.sum(np.linalg.eigvalsh(Mx[mask][..., 1:, 1:])[..., 1] > w[mask][..., 1] + 0.5 * (w[mask][..., 2] - w[mask][..., 1])))})
        log(f"cross eps {eps:+.0e}: E_J refreshed {R['E_J']:.4f} frozen {ef:.4f} E_stat {R['E_stat']:.4f} kin refreshed {R['kin_tot']:.3f} frozen {ktf:.3f}")
    sec["crossing_line"] = rows
    # (2) random directions: crossing cells masked OUT vs ONLY the crossing cells
    for label, fm in (("free_minus_crossing", free & ~mask), ("crossing_only", mask)):
        out = []
        for k in range(4):
            Dk = rand_dir(rng, Me.shape, fm)
            row = {"k": k}
            for eps in (5e-4, 1e-3, 2e-3):
                g1, g2, _ = dderiv(Me, Dk, cfg, eps)
                row[f"dE_ds_eps{eps:g}"], row[f"d2_eps{eps:g}"] = g1, g2
            out.append(row)
            log(f"{label} dir {k}: dE/ds {row['dE_ds_eps0.001']:+.3e} d2 {row['d2_eps0.0005']:.3e}/{row['d2_eps0.001']:.3e}/{row['d2_eps0.002']:.3e}")
        sec[f"random_{label}"] = out
    # (3) the alternative closest-pair reading on both fields
    sec["closest_pair_reading"] = {"seed": closest_pair_reads(Ms, cfg), "end": closest_pair_reads(Me, cfg)}
    log(f"closest-pair reading end: {sec['closest_pair_reading']['end']}")
    # (4) the spectrum along the (1,1,1) diagonal ray from the center, both fields
    c = n // 2
    ray = []
    for k in range(c, n):
        ij = (k, k, k)
        ws, we = np.linalg.eigvalsh(Ms[ij][1:, 1:]), np.linalg.eigvalsh(Me[ij][1:, 1:])
        ray.append({"r": float(r[ij]), "seed_w": [float(v) for v in ws], "end_w": [float(v) for v in we]})
    sec["spectrum_ray_111"] = ray
    ray = []
    for k in range(c, n):
        ij = (k, c, c)
        ws, we = np.linalg.eigvalsh(Ms[ij][1:, 1:]), np.linalg.eigvalsh(Me[ij][1:, 1:])
        ray.append({"r": float(r[ij]), "seed_w": [float(v) for v in ws], "end_w": [float(v) for v in we]})
    sec["spectrum_ray_100"] = ray
    # (5) the 8 innermost cells: which crossed
    inner = [(i, j, k) for i in (c - 1, c) for j in (c - 1, c) for k in (c - 1, c)]
    sec["innermost_8"] = [{"idx": list(ij), "end_w": [float(v) for v in np.linalg.eigvalsh(Me[ij][1:, 1:])],
                           "crossing": bool(mask[ij])} for ij in inner]
    # (6) kin_I1 concentration
    ki = kin_i1_density(Me, a0e, cfg)
    order = np.sort(ki.ravel())[::-1]
    cum = np.cumsum(order) / order.sum()
    sec["kin_I1_concentration"] = {"top6_frac": float(cum[5]), "top8_frac": float(cum[7]),
                                   "top100_frac": float(cum[99]), "in_crossing_cells_frac": float(ki[mask].sum() / ki.sum())}
    res["crossing"] = sec


# ------------------------------------------------------------------ main
def main():
    cfg = INS4.base_cfg(s=-1, g=G, n=32, L=48.0, delta=DELTA)
    n, h, L = cfg["n"], cfg["h"], cfg["L"]
    h3 = h ** 3
    X, Y, Z = INS4.coords(n, h)
    r = np.sqrt(X * X + Y * Y + Z * Z)
    pin = INS4.pin_shell(n, h)
    free = ~pin
    res = {"argv": ARGV, "cfg": {k: v for k, v in cfg.items()}, "files": {"seed": os.path.relpath(SEED, os.path.dirname(os.path.dirname(CK))), "end": os.path.relpath(END, os.path.dirname(os.path.dirname(CK)))},   # research-relative (no machine-local paths in tracked JSON)
           "n48_exists": os.path.exists(END48)}

    Ms = np.load(SEED)
    Me = np.load(END)
    log(f"loaded seed {Ms.shape} end {Me.shape}; n48 exists: {res['n48_exists']}")
    res["field_checks"] = {
        "seed_sym_max": float(np.max(np.abs(Ms - Ms.swapaxes(-1, -2)))),
        "end_sym_max": float(np.max(np.abs(Me - Me.swapaxes(-1, -2)))),
        "end_finite": bool(np.all(np.isfinite(Me))),
        "pinned_shell_cells": int(pin.sum()), "free_cells": int(free.sum()),
        "end_minus_seed_on_pin_max": float(np.max(np.abs((Me - Ms)[pin]))),
        "end_minus_seed_max": float(np.max(np.abs(Me - Ms))),
        "end_minus_seed_fro": float(np.sqrt(np.sum((Me - Ms) ** 2))),
        "seed_M00_center": float(Ms[n // 2, n // 2, n // 2, 0, 0]),
        "vac_corner_diag": [float(v) for v in np.diag(Ms[0, 0, 0])],
    }
    # a0 cross-check against the allowed reference implementation
    try:
        R13 = _load("m5_32_r13w_common", "m5_32_r13w_common.py")
        a_ref, a_own = R13.a0_local(Me), a0_local_own(Me)
        res["field_checks"]["a0_vs_reference_max"] = float(np.max(np.abs(a_ref - a_own)))
        kd_ref = R13.kin_density(Me, a_own, cfg)
        res["field_checks"]["kin_density_vs_reference_max"] = float(
            np.max(np.abs(kd_ref - kin_i1_density(Me, a_own, cfg))))
    except Exception as ex:      # noqa: BLE001
        res["field_checks"]["a0_vs_reference_max"] = f"skipped: {ex}"

    # ---- P1 / P2: the reads
    Rs = reads(Ms, cfg, want_dens=True)
    Re = reads(Me, cfg, want_dens=True)
    ds, de = Rs.pop("_dens"), Re.pop("_dens")
    for R, d, M, name in ((Rs, ds, Ms, "seed"), (Re, de, Me, "end")):
        core = r < L / 4
        ktot = d["kin_I1"] + CP * d["kin_KP23"]
        R["core_mean_split_abs"] = float(np.mean(np.abs(d["s"])[core]))
        R["core_mean_split_signed"] = float(np.mean(d["s"][core]))
        R["core_rms_split"] = float(np.sqrt(np.mean(d["s"][core] ** 2)))
        R["kin_frac_inside_L4"], R["kin_absfrac_inside_L4"] = frac_inside(ktot, r, L / 4)
        R["kin_I1_frac_inside_L4"] = frac_inside(d["kin_I1"], r, L / 4)[0]
        R["kin_KP23_frac_inside_L4"] = frac_inside(d["kin_KP23"], r, L / 4)[0]
        R["split_frac_inside_L4"] = frac_inside(d["split"], r, L / 4)[0]
        R["kin_I1_neg_sum"] = float(d["kin_I1"][d["kin_I1"] < 0].sum())
        R["kin_I1_frac_on_pin"] = float(d["kin_I1"][pin].sum() / d["kin_I1"].sum())
        R["kin_KP23_on_pin"] = float(d["kin_KP23"][pin].sum())
        R["r_of_max_split"] = float(r[np.unravel_index(np.argmax(np.abs(d["s"])), r.shape)])
        R["r_of_max_kin"] = float(r[np.unravel_index(np.argmax(ktot), r.shape)])
        # the (1,2) share
        a12 = B8.G3 @ M - M @ B8.G3
        R["kin_I1_a0_12"] = float(INS4.kin_of(M, a12, cfg))
        R["share_12"] = R["kin_I1_a0_12"] / R["kin_I1"]
        # the split-inertia identity per cell
        ratio = np.where(d["split"] > 1e-12 * d["split"].max(), d["kin_KP23"] / np.maximum(d["split"], 1e-300), np.nan)
        R["identity_kinKP23_over_SPLIT_total"] = R["kin_KP23"] / R["SPLIT"]
        R["identity_cell_ratio_median"] = float(np.nanmedian(ratio))
        R["identity_cell_ratio_min_max"] = [float(np.nanmin(ratio)), float(np.nanmax(ratio))]
        # does kin_I1 follow the split?
        m = d["split"] > 1e-6 * d["split"].max()
        R["kin_I1_over_split_cell_median"] = float(np.median(d["kin_I1"][m] / d["split"][m]))
        R["kin_I1_over_split_cell_p10_p90"] = [float(np.percentile(d["kin_I1"][m] / d["split"][m], 10)),
                                                float(np.percentile(d["kin_I1"][m] / d["split"][m], 90))]
        R["corr_kinI1_split"] = float(np.corrcoef(d["kin_I1"].ravel(), d["split"].ravel())[0, 1])
        # spatial-block split vs N split (the time-space coupling)
        ws = np.linalg.eigvalsh(M[..., 1:, 1:])
        s_sp = ws[..., 1] - ws[..., 0]
        R["max_abs_spatial_split_minus_N_split"] = float(np.max(np.abs(s_sp - np.abs(d["s"]))))
        R["max_abs_M0i"] = float(np.max(np.abs(M[..., 0, 1:])))
        R["radial_split"] = radial(d["split"], r, h, L)
        R["radial_kin_tot"] = radial(ktot, r, h, L)
        R["radial_kin_I1"] = radial(d["kin_I1"], r, h, L)
        R["radial_E_u_proxy_v4"] = radial(d["v4"], r, h, L)
        R["radial_absS_max"] = [float(np.abs(d["s"])[(r >= p["r_lo"]) & (r < p["r_hi"])].max())
                                for p in R["radial_split"]]
        log(f"{name}: E_J {R['E_J']:.4f} E_stat {R['E_stat']:.4f} E_u {R['E_u']:.4f} V4dd {R['V4dd']:.4f} "
            f"muSPLIT {R['muSPLIT']:.3e} KP23 {R['KP23']:.4f} kin_I1 {R['kin_I1']:.4e} "
            f"kin_KP23 {R['kin_KP23']:.4e} omega {R['omega']:.4f} maxsplit {R['max_split']:.4f} "
            f"core|s| {R['core_mean_split_abs']:.4f} frac<L/4 {R['kin_frac_inside_L4']:.4f} share12 {R['share_12']:.3e}")
    res["seed"], res["end"] = Rs, Re

    # ---- P3: stationarity by finite differences (a0 refreshed)
    rng = np.random.default_rng(20260905)
    nd = 4 if QUICK else 12
    dirs = []
    for k in range(nd):
        D = rand_dir(rng, Me.shape, free)
        row = {"k": k}
        for eps in ((1e-3,) if QUICK else (1e-3, 2e-3)):
            g1, g2, e0 = dderiv(Me, D, cfg, eps)
            row[f"dE_ds_eps{eps:g}"] = g1
            row[f"d2E_ds2_eps{eps:g}"] = g2
        dirs.append(row)
        log(f"random dir {k}: dE/ds {row['dE_ds_eps0.001']:+.4e} d2E/ds2 {row['d2E_ds2_eps0.001']:.4e}")
    res["random_dirs"] = dirs
    g1s = np.array([d["dE_ds_eps0.001"] for d in dirs])
    Dfree = int(free.sum()) * 10
    res["grad_norm_estimate"] = {"D_free_components": Dfree,
                                 "|grad| ~ sqrt(D mean(g.d)^2)": float(np.sqrt(Dfree * np.mean(g1s ** 2))),
                                 "dE_ds_min_max": [float(g1s.min()), float(g1s.max())],
                                 "d2_min_max": [float(min(d["d2E_ds2_eps0.001"] for d in dirs)),
                                                float(max(d["d2E_ds2_eps0.001"] for d in dirs))]}
    # same at the seed (is the seed stationary for E_stat alone? for E_J?)
    seed_dirs = []
    for k in range(3 if QUICK else 6):
        D = rand_dir(rng, Ms.shape, free)
        g1, g2, e0 = dderiv(Ms, D, cfg, 1e-3)
        Rp = reads(Ms + 1e-3 * D, cfg); Rm = reads(Ms - 1e-3 * D, cfg)
        seed_dirs.append({"k": k, "dEJ_ds": g1, "d2EJ_ds2": g2,
                          "dEstat_ds": (Rp["E_stat"] - Rm["E_stat"]) / 2e-3,
                          "dkin_tot_ds": (Rp["kin_tot"] - Rm["kin_tot"]) / 2e-3})
    res["seed_random_dirs"] = seed_dirs
    log(f"seed random dirs: dE_J/ds {[f'{d['dEJ_ds']:+.3e}' for d in seed_dirs]}")

    # ---- the segment seed -> end
    Dseg = Me - Ms
    seg_norm = float(np.sqrt(np.sum(Dseg ** 2)))
    ss = [0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 1.0, 1.05, 1.1, 1.25, 1.5] if not QUICK else [0.0, 0.5, 1.0, 1.5]
    seg = []
    for s in ss:
        R = reads(Ms + s * Dseg, cfg)
        seg.append({"s": s, **{k: R[k] for k in ("E_J", "E_stat", "E_u", "V4dd", "muSPLIT", "KP23",
                                                  "kin_I1", "kin_KP23", "E_fixJ", "max_split", "omega")}})
        log(f"segment s={s:5.2f}: E_J {R['E_J']:.4f} E_stat {R['E_stat']:.4f} kin_tot {R['kin_tot']:.4e} maxsplit {R['max_split']:.4f}")
    res["segment"] = {"norm": seg_norm, "points": seg}
    # derivative along the unit segment direction at the end point
    Dhat = Dseg / seg_norm
    g1, g2, _ = dderiv(Me, Dhat, cfg, 1e-3)
    res["segment"]["dEJ_ds_at_end_unit_dir"] = g1
    res["segment"]["d2EJ_ds2_at_end_unit_dir"] = g2
    log(f"unit segment dir at end: dE/ds {g1:+.4e} d2 {g2:.4e}")

    # ---- mutations: rescale the end field's split by f in the spatial eigenframe
    muts = []
    for f in ([0.5, 1.5] if QUICK else [0.0, 0.5, 0.8, 0.9, 1.1, 1.2, 1.5, 2.0]):
        Mf = split_rescale(Me, f, free)
        R = reads(Mf, cfg)
        muts.append({"f": f, **{k: R[k] for k in ("E_J", "E_stat", "E_u", "V4dd", "muSPLIT", "KP23",
                                                  "kin_I1", "kin_KP23", "E_fixJ", "max_split")}})
        log(f"split x{f:4.2f}: E_J {R['E_J']:.4f} E_stat {R['E_stat']:.4f} kin_I1 {R['kin_I1']:.4e} kin_KP23 {R['kin_KP23']:.4e}")
    res["split_rescale"] = muts
    # derivative w.r.t. the split scale at f=1
    Rp, Rm = reads(split_rescale(Me, 1.001, free), cfg), reads(split_rescale(Me, 0.999, free), cfg)
    res["dEJ_df_split_at_1"] = (Rp["E_J"] - Rm["E_J"]) / 0.002
    res["dEstat_df_split_at_1"] = (Rp["E_stat"] - Rm["E_stat"]) / 0.002
    res["dkin_df_split_at_1"] = (Rp["kin_tot"] - Rm["kin_tot"]) / 0.002
    log(f"d E_J / d f(split) at f=1: {res['dEJ_df_split_at_1']:+.4e} (E_stat {res['dEstat_df_split_at_1']:+.4e}, kin {res['dkin_df_split_at_1']:+.4e})")

    # ---- E_J vs J on the end field (how visible is the fixed-J term)
    crossing_section(Ms, Me, cfg, free, r, rng, res)

    res["fixJ_visibility"] = {"E_fixJ_over_E_J": Re["E_fixJ"] / Re["E_J"],
                              "E_fixJ_over_E_stat": Re["E_fixJ"] / Re["E_stat"]}

    if res["n48_exists"]:
        try:
            M48 = np.load(END48)
            cfg48 = INS4.base_cfg(s=-1, g=G, n=48, L=72.0, delta=DELTA)
            R48 = reads(M48, cfg48)
            res["end_n48"] = {"shape": list(M48.shape), "mtime": os.path.getmtime(END48), **R48}
            log(f"n48 end: E_J {R48['E_J']:.4f} E_stat {R48['E_stat']:.4f} kin_tot {R48['kin_tot']:.4e}")
        except Exception as ex:  # noqa: BLE001
            res["end_n48"] = f"failed: {ex}"

    res["runtime_s"] = time.time() - T0
    out = os.path.join(DATA, "m5_32_r15_p4_audit.json")
    with open(out, "w") as fh:
        json.dump(res, fh, indent=1, default=float)
    log(f"wrote {out}")


if __name__ == "__main__":
    main()
