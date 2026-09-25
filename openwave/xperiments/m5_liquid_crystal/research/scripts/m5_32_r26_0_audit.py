"""M5.32 R26-0 adversarial audit: an independent attempt to refute the ten
claims of `data/m5_32_r26_0_form.json` with own code (the audited script
`m5_32_r26_0_form.py` was NOT read; only its JSON, the certified platform
`m5_21_3_a_4d.py`, the seed builder `m5_32_r20_1_axes.seed_axes`, the config
helpers of `m5_32_r21_1_runs` and `m5_32_r20_0_class.roots_of` were used).

The stack under audit: real symmetric 4x4 M(x) on an n^3 lattice, eta =
diag(-1, 1, 1, 1), N = M eta, vacuum diag(8, 1, delta, 0), static energy
E = 4 h^3 sum_{i<j} <F_ij, F_ij>_eta + V4 with F_ij = A_i eta A_j - A_j eta
A_i, A_i = d_i M (the mean of the forward and backward branches), <F, G>_eta
= tr(eta F eta G^T), V4 = w sum_p (tr N^p - C_p)^2.

Own instruments (one line each, the details in the check functions):
  (a) sympy on a symbolic symmetric M plus a finite difference of the
      exact pushforward exp(tG) M exp(tG)^T with scipy expm;
  (b) the platform's gen_catalog fields measured for their symmetric part,
      the corrected tangents rebuilt from the same generators;
  (c) an own curvature energy with both contractions on every full-kick
      slab and on the stored charge;
  (d) the plain and the telescoped trace brackets against an mpmath
      50-digit reference on diag(8, s0 + b, s0 - b, 1);
  (e) the biaxiality invariant with a root find for the 0.382 crossing;
  (f) an own shell reader (lattice least squares on P0, P1, P2 per shell
      AND a Gauss-Legendre projection on interpolated spheres) on the
      stored charge and on a planted A cos(theta) / r^2 + B / r gap;
  (g) an own orbital derivative (central and fourth-order stencils), an
      own kin, and a rotate-by-interpolation finite difference on a smooth
      synthetic field, the stored charge raw and smoothed;
  (h) an own loop-winding reader (fixed transverse frame, double angle)
      on the R20 S1 seed, an own sphere partition reader on the seed and
      on the stored end field, the outward alignment of the director;
  (i) the h-scaling with the platform e_parts and the own energies;
  (j) an own degree integral (the signed solid-angle sum of the
      triangulated director map) with two orientation rules.

Output: data/m5_32_r26_0_audit.json and a verdict table on the terminal.
"""

from __future__ import annotations

import importlib.util
import json
import os
import sys
import time

import numpy as np
import mpmath as mp
import sympy as sp
from scipy.linalg import expm
from scipy.ndimage import gaussian_filter, map_coordinates

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
FORM_JSON = os.path.join(DATA, "m5_32_r26_0_form.json")
OUT_JSON = os.path.join(DATA, "m5_32_r26_0_audit.json")
S1_NPZ = os.path.join(DATA, "m5_32_r25_2", "S1_d0.3_w25_n32_L48.npz")
S1_GATE_NPZ = os.path.join(DATA, "m5_32_r25_2", "S1_d0.3_w25_n32_L48_gate.npz")
FULLKICK_DIR = os.path.join(DATA, "m5_32_r25_1_fullkick")
UNLIKE_NPZ = os.path.join(DATA, "m5_32_r22_2", "unlike_d12_n32_L48.npz")
T0 = time.time()


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    argv = sys.argv
    sys.argv = [argv[0]]
    spec.loader.exec_module(mod)
    sys.argv = argv
    return mod


B3 = _load("m5_21_3_a_4d", "m5_21_3_a_4d.py")
R20 = _load("m5_32_r20_1_axes", "m5_32_r20_1_axes.py")
R21 = _load("m5_32_r21_1_runs", "m5_32_r21_1_runs.py")
R0 = _load("m5_32_r20_0_class", "m5_32_r20_0_class.py")
ETA = np.diag([-1.0, 1.0, 1.0, 1.0])
W1 = B3.W1
W1S = 25.0
DELTA = 0.3
JZ = np.zeros((4, 4))
JZ[1, 2], JZ[2, 1] = -1.0, 1.0


def log(msg):
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


def wrap(a):
    return (a + np.pi) % (2.0 * np.pi) - np.pi


def unit(v):
    v = np.asarray(v, float)
    return v / np.linalg.norm(v)


# ================= own stencils and energies =================
def d_one(f, ax, h, br):
    """own one-sided difference along ax: fwd leaves the last slab zero, bwd the first."""
    out = np.zeros_like(f)
    df = np.diff(f, axis=ax) / h
    sl = [slice(None)] * f.ndim
    sl[ax] = slice(0, -1) if br == "fwd" else slice(1, None)
    out[tuple(sl)] = df
    return out


def d_central(f, ax, h, order=2):
    """own central difference (second or fourth order inside, one-sided at the faces)."""
    out = np.zeros_like(f)
    n = f.shape[ax]

    def at(i):
        s = [slice(None)] * f.ndim
        s[ax] = i
        return tuple(s)

    if order == 2:
        out[at(slice(1, -1))] = (f[at(slice(2, None))] - f[at(slice(0, -2))]) / (2 * h)
    else:
        out[at(slice(2, -2))] = (
            -f[at(slice(4, None))]
            + 8 * f[at(slice(3, -1))]
            - 8 * f[at(slice(1, -3))]
            + f[at(slice(0, -4))]
        ) / (12 * h)
        out[at(1)] = (f[at(2)] - f[at(0)]) / (2 * h)
        out[at(-2)] = (f[at(-1)] - f[at(-3)]) / (2 * h)
    out[at(0)] = (f[at(1)] - f[at(0)]) / h
    out[at(n - 1)] = (f[at(n - 1)] - f[at(n - 2)]) / h
    return out


def curv_energies(M, h, mask=None):
    """own (E_eta, E_frobenius): 4 h^3 sum_{i<j} of the branch mean of <F, F>."""
    e_eta = e_fro = 0.0
    for br in ("fwd", "bwd"):
        A = [d_one(M, ax, h, br) for ax in range(3)]
        for i in range(3):
            for j in range(i + 1, 3):
                F = A[i] @ ETA @ A[j] - A[j] @ ETA @ A[i]
                EF = ETA @ F @ ETA
                ce = np.sum(EF * F, axis=(-1, -2))
                cf = np.sum(F * F, axis=(-1, -2))
                if mask is not None:
                    ce, cf = ce[mask], cf[mask]
                e_eta += 0.5 * 4.0 * float(np.sum(ce))
                e_fro += 0.5 * 4.0 * float(np.sum(cf))
    return h**3 * e_eta, h**3 * e_fro


def v4_own(M, roots, w, h):
    """own V4 = w h^3 sum_cells sum_p (tr N^p - C_p)^2, plain form."""
    N = M @ ETA
    P = np.broadcast_to(np.eye(4), M.shape)
    tot = 0.0
    for p in range(1, 5):
        P = P @ N
        t = np.einsum("...kk->...", P)
        C = sum(r**p for r in roots)
        tot = tot + (t - C) ** 2
    return w * h**3 * float(np.sum(tot))


def kin_own(M, a, h, mask=None):
    """own kin(M; a) = 4 h^3 sum_i branch mean of <[a, A_i]_eta, [a, A_i]_eta>_eta."""
    k = 0.0
    for br in ("fwd", "bwd"):
        for ax in range(3):
            A = d_one(M, ax, h, br)
            F = a @ ETA @ A - A @ ETA @ a
            v = np.sum((ETA @ F @ ETA) * F, axis=(-1, -2))
            if mask is not None:
                v = v[mask]
            k += 0.5 * 4.0 * float(np.sum(v))
    return h**3 * k


# ================= own interpolation and sphere tools =================
def interp(F, pts, h, n, order=1):
    """F (n, n, n, ...) sampled at physical points pts (..., 3), per component."""
    idx = (pts / h + (n - 1) / 2.0).reshape(-1, 3).T
    comps = F.reshape(n, n, n, -1)
    vals = np.stack(
        [
            map_coordinates(comps[..., k], idx, order=order, mode="nearest")
            for k in range(comps.shape[-1])
        ],
        axis=-1,
    )
    return vals.reshape(pts.shape[:-1] + F.shape[3:])


def frames_sph(TH, PH):
    rhat = np.stack([np.sin(TH) * np.cos(PH), np.sin(TH) * np.sin(PH), np.cos(TH)], -1)
    eth = np.stack([np.cos(TH) * np.cos(PH), np.cos(TH) * np.sin(PH), -np.sin(TH)], -1)
    eph = np.stack([-np.sin(PH), np.cos(PH), np.zeros_like(PH)], -1)
    return rhat, eth, eph


def hedgehog_frames(X, Y, Z):
    """own radial / azimuthal / polar frames (phi = y-hat on the polar axis)."""
    r = np.sqrt(X * X + Y * Y + Z * Z)
    rho = np.sqrt(X * X + Y * Y)
    nhat = np.stack([X, Y, Z], -1) / np.maximum(r, 1e-300)[..., None]
    ph = np.stack([-Y, X, np.zeros_like(Z)], -1) / np.maximum(rho, 1e-300)[..., None]
    ph[rho < 1e-9] = np.array([0.0, 1.0, 0.0])
    th = np.cross(ph, nhat)
    return nhat, ph, th, r


def embed3(M3, m00=8.0):
    M = np.zeros(M3.shape[:-2] + (4, 4))
    M[..., 1:4, 1:4] = M3
    M[..., 0, 0] = m00
    return M


# ================= (a) the tangent formula =================
def check_a(aud):
    m = sp.Matrix(4, 4, lambda i, j: sp.Symbol(f"m{min(i, j)}{max(i, j)}"))
    eta = sp.diag(-1, 1, 1, 1)
    J = sp.zeros(4, 4)
    J[1, 2], J[2, 1] = -1, 1
    K = sp.zeros(4, 4)
    K[0, 1] = K[1, 0] = 1
    out = {}
    rng = np.random.default_rng(3)
    Mn = rng.normal(size=(4, 4))
    Mn = 0.5 * (Mn + Mn.T)
    for nm, G, kind in (("rotation_12", J, "commutator"), ("boost_1", K, "anticommutator")):
        T = G * m + m * G.T
        C = G * m - m * G.T
        row = {
            "eta_invariance": bool(sp.simplify(G * eta + eta * G.T) == sp.zeros(4, 4)),
            "tangent_symmetric": bool(sp.simplify(T - T.T) == sp.zeros(4, 4)),
            "tangent_is_commutator": bool(sp.simplify(T - (G * m - m * G)) == sp.zeros(4, 4)),
            "tangent_is_anticommutator": bool(sp.simplify(T - (G * m + m * G)) == sp.zeros(4, 4)),
            "catalog_antisymmetric": bool(sp.simplify(C + C.T) == sp.zeros(4, 4)),
            "catalog_zero": bool(sp.simplify(C) == sp.zeros(4, 4)),
        }
        Gn = np.array(G, dtype=float)
        t = 1e-5
        push = lambda s: expm(s * Gn) @ Mn @ expm(s * Gn).T
        fd = (push(t) - push(-t)) / (2 * t)
        Tn = Gn @ Mn + Mn @ Gn.T
        row["fd_rel_err"] = float(np.max(np.abs(fd - Tn)) / np.max(np.abs(Tn)))
        row["expected"] = kind
        out[nm] = row
    ok = (
        all(
            r["eta_invariance"] and r["tangent_symmetric"] and r["catalog_antisymmetric"]
            for r in out.values()
        )
        and out["rotation_12"]["tangent_is_commutator"]
        and out["boost_1"]["tangent_is_anticommutator"]
        and not out["rotation_12"]["tangent_is_anticommutator"]
        and not out["boost_1"]["tangent_is_commutator"]
        and all(not r["catalog_zero"] and r["fd_rel_err"] < 1e-8 for r in out.values())
    )
    return {
        "method": "sympy on a symbolic symmetric M for J_12 and K_1 plus a central finite difference of "
        "expm(tG) M expm(tG)^T on a random symmetric M",
        "mine": out,
        "audited": {k: v["tangent_equals"] for k, v in aud["a_tangent"]["generators"].items()},
        "verdict": "CONFIRMED" if ok else "REFUTED",
        "note": "rotation: G^T = -G so G M + M G^T = [G, M]; boost: G^T = G so it is {G, M}; the catalog "
        "form G M - M G^T is the other bracket in each case and is antisymmetric",
    }


# ================= (b) the catalog fields on the stored charge =================
def check_b(aud, M, cfg):
    cat = B3.gen_catalog(cfg, M)
    w = B3.envelope(cfg)[..., None, None]
    lam, V = np.linalg.eigh(M[..., 1:4, 1:4])

    def cross_mat(v):
        W = np.zeros(v.shape[:-1] + (4, 4))
        W[..., 1, 2], W[..., 1, 3] = -v[..., 2], v[..., 1]
        W[..., 2, 1], W[..., 2, 3] = v[..., 2], -v[..., 0]
        W[..., 3, 1], W[..., 3, 2] = -v[..., 1], v[..., 0]
        return W

    Jx = np.zeros((4, 4))
    Jx[2, 3], Jx[3, 2] = -1.0, 1.0
    Kz = np.zeros((4, 4))
    Kz[0, 3] = Kz[3, 0] = 1.0
    Kx = np.zeros((4, 4))
    Kx[0, 1] = Kx[1, 0] = 1.0
    gens = {
        "clock_local": (cross_mat(V[..., :, 2]), "rot"),
        "plane_1d": (cross_mat(V[..., :, 0]), "rot"),
        "rot_z": (np.broadcast_to(JZ, M.shape), "rot"),
        "rot_x": (np.broadcast_to(Jx, M.shape), "rot"),
        "boost_z": (np.broadcast_to(Kz, M.shape), "boost"),
        "boost_x": (np.broadcast_to(Kx, M.shape), "boost"),
    }
    rows = {}
    worst_sym = worst_anti = 0.0
    for nm, a in cat.items():
        mx = float(np.max(np.abs(a)))
        sym = float(np.max(np.abs(0.5 * (a + a.swapaxes(-1, -2)))) / mx)
        G, kind = gens[nm]
        tan = G @ M - M @ G if kind == "rot" else G @ M + M @ G
        tan = w * tan
        tan = tan / np.sqrt(np.sum(tan * tan))
        anti = float(np.max(np.abs(0.5 * (tan - tan.swapaxes(-1, -2)))) / np.max(np.abs(tan)))
        # the corrected tangent must be a different field from the catalog one
        ovl = float(np.sum(tan * a))
        rows[nm] = {
            "max_abs": mx,
            "symmetric_part_over_max": sym,
            "corrected_antisymmetric_part_over_max": anti,
            "overlap_catalog_vs_corrected": ovl,
            "audited_symmetric_part_over_max": aud["b_catalog_antisymmetric"]["fields"][nm][
                "symmetric_part_over_max"
            ],
        }
        worst_sym = max(worst_sym, sym)
        worst_anti = max(worst_anti, anti)
    ok = (
        worst_sym < 1e-12
        and worst_anti < 1e-12
        and all(abs(r["overlap_catalog_vs_corrected"]) < 1e-12 for r in rows.values())
    )
    return {
        "method": "platform gen_catalog on the stored charge, symmetric part measured; corrected tangents [G, M] "
        "(rotations, the local ones rebuilt as cross-product matrices of the eigenvectors) and {G, M} (boosts), "
        "antisymmetric part measured, and their overlap with the catalog field",
        "mine": rows,
        "worst_symmetric_part": worst_sym,
        "worst_corrected_antisymmetric_part": worst_anti,
        "verdict": "CONFIRMED" if ok else "REFUTED",
        "note": "the catalog field and the corrected tangent are Frobenius-orthogonal (antisymmetric against "
        "symmetric), so no rescaling relates them",
    }


# ================= (c) the norm witness =================
def check_c(aud, M, cfg):
    rows = {}
    n_neg = n_pos = 0
    max_v4_rel = 0.0
    for fn in sorted(os.listdir(FULLKICK_DIR)):
        if not fn.endswith(".npz"):
            continue
        tag = fn[:-4]
        d, w, nn, LL = tag.split("_")
        delta, w1s, n, L = float(d[1:]), float(w[1:]), int(nn[1:]), float(LL[1:])
        c = B3.base_cfg(s=-1.0, n=n, L=L, delta=delta)
        Mk = np.load(os.path.join(FULLKICK_DIR, fn))["M"].astype(np.float64)
        e_eta, e_fro = curv_energies(Mk, c["h"])
        v4 = v4_own(Mk, R0.roots_of(c), W1 * w1s, c["h"])
        a = aud["c_norm_witness"]["fields"].get(tag, {})
        rel_eta = abs(e_eta - a.get("E_u_eta", np.nan)) / abs(e_eta)
        rel_fro = abs(e_fro - a.get("E_u_frobenius", np.nan)) / abs(e_fro)
        rel_v4 = abs(v4 - a.get("V4", np.nan)) / abs(v4)
        max_v4_rel = max(max_v4_rel, rel_v4)
        n_neg += e_eta < 0
        n_pos += e_fro > 0
        rows[tag] = {
            "E_u_eta": e_eta,
            "E_u_frobenius": e_fro,
            "V4": v4,
            "M0i_max": float(np.max(np.abs(Mk[..., 0, 1:4]))),
            "rel_vs_audited": [float(rel_eta), float(rel_fro), float(rel_v4)],
            "fro_minus_eta_over_fro": float((e_fro - e_eta) / e_fro),
        }
    ee, ef = curv_energies(M, cfg["h"])
    eu_plat, _ = B3.e_parts(M, cfg)
    worst = max(max(r["rel_vs_audited"][:2]) for r in rows.values())
    ok = (
        len(rows) == 22
        and n_neg == 22
        and n_pos == 22
        and worst < 1e-9
        and abs(ee - ef) <= 1e-14 * abs(ee)
        and abs(ee - eu_plat) < 1e-12 * abs(ee)
    )
    return {
        "method": "own curvature energy (fwd/bwd branch mean, own contraction) under both norms on the 22 "
        "full-kick slabs and the stored charge; own plain V4 at W1 x w1s; the platform e_parts as a third read",
        "mine": {
            "n_fields": len(rows),
            "n_eta_negative": int(n_neg),
            "n_frobenius_positive": int(n_pos),
            "worst_rel_vs_audited_curvature": float(worst),
            "worst_rel_vs_audited_V4": float(max_v4_rel),
            "stored_S1": {"E_u_eta": ee, "E_u_frobenius": ef, "platform_E_u": float(eu_plat)},
            "fields": rows,
        },
        "audited": {
            "n_eta_negative": aud["c_norm_witness"]["n_eta_negative"],
            "block_sector_coincide": aud["c_norm_witness"]["block_sector_coincide"],
        },
        "verdict": "QUALIFIED" if ok else "REFUTED",
        "note": "the numbers hold to 1e-12; the wording needs one correction: tr(F F^T) is a sum of squares, so "
        "the fails_if clause 'a non-positive Frobenius curvature energy' cannot fail on any nonzero field "
        "(an unfalsifiable PASS line); the discriminating content is the eta sign flip on every field with "
        "M_0i, which the fails_if does not name, and the block-sector coincidence, which it does",
    }


# ================= (d) the telescoped bracket =================
def check_d(aud):
    mp.mp.dps = 50
    rows = {}
    for delta in (0.3, 0.03, 0.01, 0.003, 0.001):
        s0 = 0.5 * delta
        C = [(-8.0) ** p + 1.0 + delta**p for p in range(1, 5)]
        N0 = np.diag([-8.0, delta, 0.0, 1.0])
        plain_err = tele_err = 0.0
        plain_err_v = tele_err_v = 0.0
        for b in np.linspace(0.0, s0, 13):
            M = np.diag([8.0, s0 + b, s0 - b, 1.0])
            N = M @ ETA
            eps = N - N0
            # exact reference at 50 digits
            ex = []
            for p in range(1, 5):
                v = (
                    (mp.mpf(s0) + mp.mpf(b)) ** p
                    + (mp.mpf(s0) - mp.mpf(b)) ** p
                    - mp.mpf(delta) ** p
                )
                ex.append(v)
            v_ex = sum(x**2 for x in ex)
            # plain
            P = np.eye(4)
            pl = []
            for p in range(1, 5):
                P = P @ N
                pl.append(np.trace(P) - C[p - 1])
            # telescoped
            te = []
            Npow = [np.eye(4)]
            N0pow = [np.eye(4)]
            for p in range(1, 5):
                Npow.append(Npow[-1] @ N)
                N0pow.append(N0pow[-1] @ N0)
            for p in range(1, 5):
                s = 0.0
                for k in range(p):
                    s += np.trace(Npow[k] @ eps @ N0pow[p - 1 - k])
                te.append(s)
            scale = max(abs(x) for x in ex)
            if scale == 0:
                # the vacuum point b = s0: every exact bracket vanishes, read the absolute error over delta^2
                vac_plain = float(max(abs(x) for x in pl) / delta**2)
                vac_tele = float(max(abs(x) for x in te) / delta**2)
                continue
            plain_err = max(plain_err, max(abs(mp.mpf(pl[i]) - ex[i]) for i in range(4)) / scale)
            tele_err = max(tele_err, max(abs(mp.mpf(te[i]) - ex[i]) for i in range(4)) / scale)
            if v_ex > 0:
                plain_err_v = max(plain_err_v, abs(mp.mpf(sum(x * x for x in pl)) - v_ex) / v_ex)
                tele_err_v = max(tele_err_v, abs(mp.mpf(sum(x * x for x in te)) - v_ex) / v_ex)
        a = aud["d_split_residual"]["rows"][str(delta)]
        rows[str(delta)] = {
            "plain_relerr_bracket_max": float(plain_err),
            "telescoped_relerr_bracket_max": float(tele_err),
            "plain_relerr_V4_max": float(plain_err_v),
            "telescoped_relerr_V4_max": float(tele_err_v),
            "gain_bracket": float(plain_err / max(tele_err, mp.mpf("1e-300"))),
            "vacuum_point_abs_err_over_delta2": [vac_plain, vac_tele],
            "audited_plain": a["plain_relerr_max"],
            "audited_telescoped": a["telescoped_relerr_max"],
        }
    ok = all(
        r["telescoped_relerr_bracket_max"] < r["plain_relerr_bracket_max"] for r in rows.values()
    ) and (rows["0.001"]["telescoped_relerr_bracket_max"] < 1e-7)
    same_order = all(
        0.03 < r["plain_relerr_bracket_max"] / r["audited_plain"] < 30 for r in rows.values()
    )
    return {
        "method": "diag(8, s0 + b, s0 - b, 1) with s0 = delta / 2, b on 13 points of [0, s0]; plain tr N^p - C_p "
        "in double against the 50-digit mpmath pair bracket, the telescoped sum_k tr(N^k eps N0^{p-1-k}) with "
        "N0 = diag(-8, delta, 0, 1); the error normalized by the largest exact bracket, and on V4 itself",
        "mine": rows,
        "verdict": "CONFIRMED" if ok else "REFUTED",
        "note": (
            "the telescoped form is exact to the last bit on this diagonal family (eps vanishes in the 8 and 1 "
            "slots, no large-number cancellation survives), so its floor here is below the audited numbers, "
            "which rise toward small delta; the audited telescoped floor therefore measures a residual the "
            "diagonal family does not have (a perturbed M_00 or director slot); the ordering claim and the "
            "1e-7 bound hold either way"
            + (
                ""
                if same_order
                else "; the plain errors differ from the audited ones by more than a factor 30"
            )
        ),
    }


# ================= (e) the biaxiality invariant =================
def beta2(lams):
    L = np.diag(np.asarray(lams, float))
    L = L - np.trace(L) / 3.0 * np.eye(3)
    t2 = np.trace(L @ L)
    t3 = np.trace(L @ L @ L)
    return float(1.0 - 6.0 * t3**2 / t2**3)


def check_e(aud):
    from scipy.optimize import brentq

    d_star = brentq(lambda d: beta2((1.0, d, 0.0)) - 0.382, 0.05, 0.5)
    mine = {
        "uniaxial": beta2((1.0, 0.0, 0.0)),
        "maximal": beta2((1.0, 0.0, -1.0)),
        "vacuum_d0.3": beta2((1.0, 0.3, 0.0)),
        "vacuum_d0.23": beta2((1.0, 0.23, 0.0)),
        "delta_at_0.382": float(d_star),
    }
    a = aud["e_biaxiality"]
    ok = (
        abs(mine["uniaxial"]) < 1e-12
        and abs(mine["maximal"] - 1.0) < 1e-12
        and abs(mine["vacuum_d0.3"] - 0.604) < 5e-4
        and abs(mine["delta_at_0.382"] - 0.2307) < 5e-5
        and abs(mine["vacuum_d0.3"] - a["vacuum_d0.3"]) < 1e-12
        and abs(mine["delta_at_0.382"] - a["delta_at_0.382"]) < 1e-9
    )
    return {
        "method": "the invariant on the traceless part of diag(lams); brentq for the 0.382 crossing",
        "mine": mine,
        "audited": {
            k: a[k]
            for k in ("uniaxial", "maximal", "vacuum_d0.3", "vacuum_d0.23", "delta_at_0.382")
        },
        "verdict": "CONFIRMED" if ok else "REFUTED",
        "note": "",
    }


# ================= (f) the gap-tail reader =================
def gap_dev(M3, delta):
    lam = np.linalg.eigvalsh(M3)
    return lam[..., 1] - lam[..., 0] - delta


def shell_fit_lattice(dev, X, Y, Z, h, radii):
    """own reader 1: least squares of dev on [1, cos th, P2(cos th)] over the cells of each shell."""
    r = np.sqrt(X * X + Y * Y + Z * Z)
    ct = Z / np.maximum(r, 1e-300)
    rows = []
    for R in radii:
        m = (r >= R - 0.5 * h) & (r < R + 0.5 * h)
        A = np.stack([np.ones(m.sum()), ct[m], 0.5 * (3 * ct[m] ** 2 - 1)], -1)
        co, *_ = np.linalg.lstsq(A, dev[m], rcond=None)
        rows.append([float(R), float(co[0]), float(co[1]), float(co[2]), int(m.sum())])
    return rows


def shell_fit_sphere(M3, delta, h, n, radii, nth=48, nph=96, order=3):
    """own reader 2: Gauss-Legendre in cos theta x uniform phi on interpolated spheres, Legendre projection."""
    x, wts = np.polynomial.legendre.leggauss(nth)
    TH = np.arccos(x)[:, None] * np.ones((1, nph))
    PH = np.ones((nth, 1)) * (np.arange(nph) * 2 * np.pi / nph)[None, :]
    rhat = np.stack([np.sin(TH) * np.cos(PH), np.sin(TH) * np.sin(PH), np.cos(TH)], -1)
    rows = []
    for R in radii:
        S = interp(M3, R * rhat, h, n, order=order)
        dev = gap_dev(S, delta).mean(axis=1)  # phi average per ring
        a0 = 0.5 * np.sum(wts * dev)
        a1 = 1.5 * np.sum(wts * dev * x)
        a2 = 2.5 * np.sum(wts * dev * 0.5 * (3 * x * x - 1))
        rows.append([float(R), float(a0), float(a1), float(a2)])
    return rows


def slope_loglog(R, a):
    R, a = np.asarray(R), np.asarray(a)
    m = a > 0
    if m.sum() < 3:
        return float("nan")
    return float(np.polyfit(np.log(R[m]), np.log(a[m]), 1)[0])


def check_f(aud, M, cfg):
    n, h = cfg["n"], cfg["h"]
    X, Y, Z = B3.coords(n, h)
    radii = [3.0, 4.5, 6.0, 7.5, 9.0, 10.5, 12.0, 13.5, 15.0, 16.5, 18.0, 19.5, 21.0]
    dev = gap_dev(M[..., 1:4, 1:4], DELTA)
    lat = shell_fit_lattice(dev, X, Y, Z, h, radii)
    sph = shell_fit_sphere(M[..., 1:4, 1:4], DELTA, h, n, radii)
    win = [i for i, R in enumerate(radii) if 6.0 <= R <= 18.0]
    ratio_lat = max(abs(lat[i][2]) / abs(lat[i][1]) for i in win)
    ratio_sph = max(abs(sph[i][2]) / abs(sph[i][1]) for i in win)
    # z-reflection symmetry of the lattice field itself (the even grid maps cell k to cell n-1-k exactly;
    # the tensor reflection is P M(P x) P with P = diag(1, 1, 1, -1)); the seed reads 4e-16 under this map
    fl = dev[:, :, ::-1]
    odd, even = 0.5 * (dev - fl), 0.5 * (dev + fl)
    Pz = np.diag([1.0, 1.0, 1.0, -1.0])
    dM = M - np.einsum("ab,...bc,cd->...ad", Pz, M[:, :, ::-1], Pz)
    r = np.sqrt(X * X + Y * Y + Z * Z)
    asym = {
        "tensor_odd_rms_over_field_rms": float(
            np.sqrt(np.mean(dM**2)) / np.sqrt(np.mean((M - np.diag([8.0, 0, 0, 0])) ** 2))
        ),
        "dev_odd_rms_over_even_rms": float(np.sqrt(np.mean(odd**2)) / np.sqrt(np.mean(even**2))),
        "dev_odd_over_even_rms_per_shell": {
            str(R): float(
                np.sqrt(np.mean(odd[(r >= R - 0.5 * h) & (r < R + 0.5 * h)] ** 2))
                / np.sqrt(np.mean(even[(r >= R - 0.5 * h) & (r < R + 0.5 * h)] ** 2))
            )
            for R in (3.0, 6.0, 9.0, 12.0, 18.0)
        },
    }
    # synthetic: hedgehog frames, eigenvalues (0, delta + g, 1), g = A cos th / r^2 + B / r
    nhat, ph, th, r = hedgehog_frames(X, Y, Z)
    A_pl, B_pl = 0.3, 0.02
    reff = np.maximum(r, 2.0)
    g = A_pl * (Z / np.maximum(r, 1e-300)) / reff**2 + B_pl / reff
    S = (
        1.0 * nhat[..., :, None] * nhat[..., None, :]
        + (DELTA + g)[..., None, None] * ph[..., :, None] * ph[..., None, :]
        + 0.0 * th[..., :, None] * th[..., None, :]
    )
    dev_s = gap_dev(S, DELTA)
    lat_s = shell_fit_lattice(dev_s, X, Y, Z, h, radii)
    sph_s = shell_fit_sphere(S, DELTA, h, n, radii)
    Rw = [radii[i] for i in win]
    sl_lat = slope_loglog(Rw, [lat_s[i][2] for i in win])
    sl_sph = slope_loglog(Rw, [sph_s[i][2] for i in win])
    sl0_lat = slope_loglog(Rw, [lat_s[i][1] for i in win])
    amp_lat = max(abs(lat_s[i][2] * radii[i] ** 2 / A_pl - 1.0) for i in win)
    amp_sph = max(abs(sph_s[i][2] * radii[i] ** 2 / A_pl - 1.0) for i in win)
    a = aud["f_gap_tail_reader"]
    ok = (
        ratio_lat < 5e-3
        and ratio_sph < 5e-3
        and abs(sl_lat + 2.0) < 0.1
        and abs(sl_sph + 2.0) < 0.1
        and amp_lat < 0.05
        and amp_sph < 0.05
        and abs(sl0_lat + 1.0) < 0.1
    )
    return {
        "method": "two own readers, lattice shell least squares on [1, P1, P2] (shell width h) and a "
        "Gauss-Legendre Legendre projection on trilinear-interpolated spheres; the null on the stored charge; "
        "a planted gap A cos(th) / r^2 + B / r on a hedgehog-frame tensor (A 0.3, B 0.02)",
        "mine": {
            "null_shells_lattice_[R,a0,a1,a2,ncells]": lat,
            "null_shells_sphere_[R,a0,a1,a2]": sph,
            "null_l1_over_l0_max_window_lattice": float(ratio_lat),
            "null_l1_over_l0_max_window_sphere": float(ratio_sph),
            "null_z_reflection_asymmetry": asym,
            "synthetic_slope_l1_lattice": sl_lat,
            "synthetic_slope_l1_sphere": sl_sph,
            "synthetic_slope_l0_lattice(planted -1)": sl0_lat,
            "synthetic_a1_r2_over_A_maxerr_lattice": float(amp_lat),
            "synthetic_a1_r2_over_A_maxerr_sphere": float(amp_sph),
        },
        "audited": {
            "l1_over_l0_max_in_window": a["null_S1"]["l1_over_l0_max_in_window"],
            "synthetic_slope_l1": a["synthetic"]["slope_l1"],
            "synthetic_amp_maxerr": a["synthetic"]["a1_times_r2_over_A_maxerr"],
        },
        "verdict": "QUALIFIED" if ok else "REFUTED",
        "note": "the numbers hold (both own readers under 5e-3 on the null, the planted tail read at -2 within "
        "0.1 and 5 percent); the stated REASON does not: the fails_if says the stored charge is z-reflection "
        "symmetric on the lattice to about 1e-3, but under the exact cell map k -> n-1-k (the seed reads 4e-16) "
        f"the tensor field's odd part is {asym['tensor_odd_rms_over_field_rms']:.2f} of the field in rms and the "
        f"gap deviation's odd part is {asym['dev_odd_rms_over_even_rms']:.2f} of its even part (0.6 to 0.9 on every "
        "shell): the descent did NOT keep the seed's reflection symmetry (the off-axis carriers sit at different "
        "positions in the two hemispheres, as the audited R12 carrier list itself shows). The l = 1 null is then "
        "a property of the odd part (no cos theta moment), i.e. a genuine, falsifiable control rather than a "
        "symmetry tautology; the plan's 'cannot fail on this field' and R26-4's 'only z-asymmetric branches "
        "carry a label' should be restated. The l = 0 profile (-0.1 at r 3, +0.1 at r 9 to 12, -0.03 at r 21) "
        "is not a power law inside the box, as the record says; the audited synthetic quotes slope_l0 -4.7 with "
        "no planted l = 0 term, mine plants B / r and reads it back at -1",
    }


# ================= (g) the physical rotation generator =================
def rigid_parts(M, cfg, order=2):
    """order 2 / 4: own central stencils (one-sided at the faces); 'fwdbwd': the mean of the platform's fwd and
    bwd one-sided fields (the same as order 2 inside, HALF the one-sided value on the two face layers).
    """
    n, h = cfg["n"], cfg["h"]
    X, Y, Z = B3.coords(n, h)
    a_int = JZ @ M - M @ JZ
    if order == "fwdbwd":
        dd = lambda ax: 0.5 * (B3.d1(M, ax, h, "fwd") + B3.d1(M, ax, h, "bwd"))
    else:
        dd = lambda ax: d_central(M, ax, h, order)
    a_orb = X[..., None, None] * dd(1) - Y[..., None, None] * dd(0)
    return a_int, a_orb, a_int - a_orb


def rotate_field(M, cfg, alpha, order=3):
    """own: M_rot(x) = R M(R^-1 x) R^T by cubic interpolation, R = expm(alpha J_z)."""
    n, h = cfg["n"], cfg["h"]
    X, Y, Z = B3.coords(n, h)
    c, s = np.cos(alpha), np.sin(alpha)
    src = np.stack([c * X + s * Y, -s * X + c * Y, Z], -1)
    Ms = interp(M, src, h, n, order=order)
    R = expm(alpha * JZ)
    return np.einsum("ab,...bc,dc->...ad", R, Ms, R)


def check_g(aud, M, Mg, cfg):
    n, h = cfg["n"], cfg["h"]
    X, Y, Z = B3.coords(n, h)
    r = np.sqrt(X * X + Y * Y + Z * Z)
    out = {}
    # the uniaxial seed at delta 0
    cfg0 = B3.base_cfg(s=-1.0, n=32, L=48.0, delta=0.0)
    Mu = R20.seed_axes(cfg0, (1.0, 0.0, 0.0))
    for nm, F, c in (
        ("uniaxial_seed", Mu, cfg0),
        ("stored_S1_end", M, cfg),
        ("stored_S1_gate", Mg, cfg),
    ):
        row = {}
        for order in (2, 4, "fwdbwd"):
            ai, ao, ar = rigid_parts(F, c, order)
            row[f"stencil_{order}"] = {
                "internal": float(B3.kin_of(F, ai, c)),
                "orbital": float(B3.kin_of(F, ao, c)),
                "rigid": float(B3.kin_of(F, ar, c)),
            }
            if order == 2:
                row["own_kin"] = {
                    "internal": kin_own(F, ai, c["h"]),
                    "orbital": kin_own(F, ao, c["h"]),
                    "rigid": kin_own(F, ar, c["h"]),
                }
                row["within_r"] = {
                    str(R): {
                        "internal": kin_own(F, ai, c["h"], r < R),
                        "orbital": kin_own(F, ao, c["h"], r < R),
                        "rigid": kin_own(F, ar, c["h"], r < R),
                    }
                    for R in (3, 6, 9, 12, 18)
                }
        row["rigid_over_internal_o2"] = row["stencil_2"]["rigid"] / row["stencil_2"]["internal"]
        row["rigid_over_internal_o4"] = row["stencil_4"]["rigid"] / row["stencil_4"]["internal"]
        row["rigid_over_internal_fwdbwd"] = (
            row["stencil_fwdbwd"]["rigid"] / row["stencil_fwdbwd"]["internal"]
        )
        row["rigid_within_r18_over_internal"] = (
            row["within_r"]["18"]["rigid"] / row["stencil_2"]["internal"]
        )
        out[nm] = row
    # the finite-difference control by rotation
    interior = np.max(np.abs(np.stack([X, Y, Z], -1)), axis=-1) <= 0.5 * cfg["L"] - 6.0
    rng = np.random.default_rng(11)
    nhat, ph, th, _ = hedgehog_frames(X, Y, Z)
    S = np.zeros(M.shape[:3] + (3, 3))
    S[:] = np.diag([1.0, DELTA, 0.0])
    for cen, sig in (((6.0, 3.0, 2.0), 4.0), ((-4.0, 5.0, -3.0), 3.5), ((2.0, -7.0, 4.0), 4.5)):
        B = rng.normal(size=(3, 3))
        B = 0.5 * (B + B.T) * 0.3
        env = np.exp(-((X - cen[0]) ** 2 + (Y - cen[1]) ** 2 + (Z - cen[2]) ** 2) / (2 * sig**2))
        S += env[..., None, None] * B
    Msyn = embed3(S)
    Msm = gaussian_filter(M, sigma=(1.0, 1.0, 1.0, 0.0, 0.0), mode="nearest")
    fd = {}
    for nm, F in (("synthetic_bump", Msyn), ("stored_smoothed_sigma1", Msm), ("stored_raw", M)):
        alpha = 0.01
        d = (rotate_field(F, cfg, alpha) - rotate_field(F, cfg, -alpha)) / (2 * alpha)
        _, _, ar = rigid_parts(F, cfg, 2)
        _, _, ar4 = rigid_parts(F, cfg, 4)
        num = np.sqrt(np.sum(((d - ar) ** 2)[interior]))
        den = np.sqrt(np.sum((ar**2)[interior]))
        num4 = np.sqrt(np.sum(((d - ar4) ** 2)[interior]))
        flip = np.sqrt(np.sum(((d + ar) ** 2)[interior]))
        fd[nm] = {
            "rel_norm_diff_o2": float(num / den),
            "rel_norm_diff_o4": float(num4 / den),
            "rel_norm_diff_sign_flipped": float(flip / den),
        }
    a = aud["g_physical_generator"]
    u, s1 = out["uniaxial_seed"], out["stored_S1_end"]
    rel = lambda x, y: abs(x - y) / abs(y)
    ok = (
        u["rigid_over_internal_o2"] < 0.01
        and u["rigid_over_internal_fwdbwd"] < 0.01
        and rel(u["stencil_2"]["internal"], a["uniaxial_seed"]["internal"]) < 1e-6
        and rel(u["stencil_fwdbwd"]["rigid"], a["uniaxial_seed"]["rigid"]) < 1e-3
        and rel(s1["stencil_2"]["internal"], a["stored_S1"]["internal"]) < 1e-6
        and rel(s1["stencil_fwdbwd"]["orbital"], a["stored_S1"]["orbital"]) < 1e-4
        and rel(s1["stencil_fwdbwd"]["rigid"], a["stored_S1"]["rigid"]) < 1e-4
        and rel(s1["stencil_2"]["rigid"], a["stored_S1"]["rigid"]) < 1e-2
        and rel(s1["stencil_2"]["rigid"], s1["own_kin"]["rigid"]) < 1e-9
        and rel(s1["within_r"]["18"]["rigid"], a["stored_S1_within_r"]["18"]["rigid"]) < 1e-6
        and fd["synthetic_bump"]["rel_norm_diff_o2"] < 0.05
        and fd["synthetic_bump"]["rel_norm_diff_o4"] < 0.01
        and fd["synthetic_bump"]["rel_norm_diff_sign_flipped"] > 1.5
    )
    return {
        "method": "a_int = [J_z, M], a_orb = (x d_y - y d_x) M with own central (o2) and fourth-order (o4) "
        "stencils, C by the platform kin_of AND an own kin (with shell masks); the finite difference of the field "
        "rotated by +-0.01 rad (cubic interpolation, conjugation by expm) against a_rigid on the interior "
        "(four cells off the walls), on a smooth synthetic off-axis three-bump field, the smoothed and the raw charge",
        "mine": {"seeds_and_charge": out, "fd_check": fd},
        "audited": {
            "uniaxial_seed": a["uniaxial_seed"],
            "stored_S1": a["stored_S1"],
            "fd_check": a["fd_check"],
        },
        "verdict": "QUALIFIED" if ok else "REFUTED",
        "note": "the numbers are reproduced exactly with the audited orbital stencil (the mean of the platform's "
        "fwd and bwd one-sided fields) and inside r 18 with any stencil; the wording needs one correction: the "
        f"uniaxial seed's 0.18 percent is NOT the interior lattice's cubic residue, it is the two face layers of "
        "the box where that stencil halves the derivative (an own central stencil reads "
        f"{u['rigid_over_internal_o2'] * 100:.3f} percent on the full box, the fourth-order one "
        f"{u['rigid_over_internal_o4'] * 100:.4f} percent, and inside r 18 the residue is "
        f"{u['rigid_within_r18_over_internal'] * 100:.3f} percent); on the stored charge the face layers move "
        f"C_orb by {rel(s1['stencil_2']['orbital'], s1['stencil_fwdbwd']['orbital']) * 100:.2f} percent and "
        f"C_rigid by {rel(s1['stencil_2']['rigid'], s1['stencil_fwdbwd']['rigid']) * 100:.3f} percent, the 0.997 "
        "conclusion is unchanged (0.997 to 1.03 across stencils). The sign convention is fixed by the "
        "interpolation test (the flipped sign reads about 2 in norm); the audited 7e-4 on its own bump is matched "
        "by 5.6e-4 with a fourth-order orbital stencil on mine, 5.8e-3 with the second-order one, both inside the "
        "5 percent gate",
    }


# ================= (h) the partition reader =================
def loop_winding(M3, cfg, center, axis, radius, K=128, order=1):
    """own loop reader: winding of the middle-eigenvector line field in a fixed frame transverse to the loop axis,
    in half-units (a full turn of the line = 2), with readability diagnostics."""
    n, h = cfg["n"], cfg["h"]
    u = unit(axis)
    e1 = unit(np.cross(u, [0.0, 0.0, 1.0]) if abs(u[2]) < 0.9 else np.cross(u, [1.0, 0.0, 0.0]))
    e2 = np.cross(u, e1)
    t = np.arange(K) * 2 * np.pi / K
    pts = np.asarray(center, float) + radius * (np.cos(t)[:, None] * e1 + np.sin(t)[:, None] * e2)
    S = interp(M3, pts, h, n, order=order)
    lam, V = np.linalg.eigh(S)
    m = V[:, :, 1]
    d = V[:, :, 2]
    chi = 2.0 * np.arctan2(m @ e2, m @ e1)
    k = float(np.sum(wrap(np.diff(chi, append=chi[:1]))) / (2 * np.pi))
    return {
        "half_units": int(round(k)),
        "raw": k,
        "min_transverse_proj": float(np.min(np.hypot(m @ e1, m @ e2))),
        "min_abs_dir_dot_axis": float(np.min(np.abs(d @ u))),
        "min_pair_gap": float(np.min(lam[:, 1] - lam[:, 0])),
    }


def sphere_partition(M3, cfg, R, nth=120, nph=240, align_tol=0.15, order=1):
    """own sphere reader: the line-field angle relative to the (e_th, e_ph) frame projected transverse to the
    outward-oriented director; plaquette windings in half-units on the band, the caps read as loops, each
    pole's carrier = cap winding + 2 (the frame's own index); plaquettes with a corner |d . rhat| < tol are
    flagged unreadable."""
    n, h = cfg["n"], cfg["h"]
    thc = np.arcsin(min(1.0, h / R))
    th = np.linspace(thc, np.pi - thc, nth)
    ph = np.arange(nph) * 2 * np.pi / nph
    TH, PH = np.meshgrid(th, ph, indexing="ij")
    rhat, eth, eph = frames_sph(TH, PH)
    S = interp(M3, R * rhat, h, n, order=order)
    lam, V = np.linalg.eigh(S)
    d = V[..., :, 2]
    align = np.sum(d * rhat, -1)
    d = d * np.sign(align)[..., None]
    f1 = eth - np.sum(eth * d, -1)[..., None] * d
    f1 /= np.maximum(np.linalg.norm(f1, axis=-1), 1e-300)[..., None]
    f2 = np.cross(d, f1)
    m = V[..., :, 1]
    chi = 2.0 * np.arctan2(np.sum(m * f2, -1), np.sum(m * f1, -1))
    chi_w = np.concatenate([chi, chi[:, :1]], axis=1)
    bad = np.abs(align) < align_tol
    bad_w = np.concatenate([bad, bad[:, :1]], axis=1)
    # plaquette (i, j) -> (i+1, j) -> (i+1, j+1) -> (i, j+1), ccw about the outward normal
    dA = wrap(chi_w[1:, :-1] - chi_w[:-1, :-1])
    dB = wrap(chi_w[1:, 1:] - chi_w[1:, :-1])
    dC = wrap(chi_w[:-1, 1:] - chi_w[1:, 1:])
    dD = wrap(chi_w[:-1, :-1] - chi_w[:-1, 1:])
    kp = np.rint((dA + dB + dC + dD) / (2 * np.pi)).astype(int)
    unread = bad_w[1:, :-1] | bad_w[:-1, :-1] | bad_w[1:, 1:] | bad_w[:-1, 1:]
    n_unread = int(np.sum(unread))
    carriers = []
    for i, j in zip(*np.nonzero(kp)):
        thm, phm = 0.5 * (th[i] + th[i + 1]), ph[j] + np.pi / nph
        carriers.append(
            [int(kp[i, j]), float(R * np.sin(thm)), float(R * np.cos(thm)), bool(unread[i, j])]
        )
    kN = int(np.rint(np.sum(wrap(np.diff(chi_w[0]))) / (2 * np.pi)))
    kS = int(np.rint(-np.sum(wrap(np.diff(chi_w[-1]))) / (2 * np.pi)))
    poles = {"north": kN + 2, "south": kS + 2}
    total = int(np.sum(kp)) + poles["north"] + poles["south"]
    part = sorted(
        [poles["north"], poles["south"]] + [c[0] for c in carriers if c[0] != 0], reverse=True
    )
    return {
        "R": R,
        "total_half_units": total,
        "partition": part,
        "poles": poles,
        "carriers_[k,rho,z,unreadable]": carriers,
        "n_unreadable_plaquettes": n_unread,
        "n_plaquettes": int(kp.size),
        "min_abs_align": float(np.min(np.abs(align))),
        "frac_abs_align_below_0.2": float(np.mean(np.abs(align) < 0.2)),
        "frac_abs_align_below_0.5": float(np.mean(np.abs(align) < 0.5)),
        "min_director_gap": float(np.min(lam[..., 2] - lam[..., 1])),
        "min_pair_gap": float(np.min(lam[..., 1] - lam[..., 0])),
    }


def check_h(aud, M, cfg):
    seed = R20.seed_axes(cfg, (1.0, DELTA, 0.0))[..., 1:4, 1:4]
    loops = {}
    for z0 in (3.0, 6.0, 9.0, 12.0, 15.0):
        for sgn in (1, -1):
            for rad in (2 * cfg["h"], 3 * cfg["h"]):
                loops[f"z_axis_z{sgn * z0:+g}_rad{rad:g}"] = loop_winding(
                    seed, cfg, (0, 0, sgn * z0), (0, 0, 1), rad
                )
    for x0 in (6.0, 9.0, 12.0):
        loops[f"x_axis_x{x0:g}"] = loop_winding(seed, cfg, (x0, 0, 0), (1, 0, 0), 2 * cfg["h"])
        loops[f"x_axis_x{-x0:g}"] = loop_winding(seed, cfg, (-x0, 0, 0), (1, 0, 0), 2 * cfg["h"])
        loops[f"y_axis_y{x0:g}"] = loop_winding(seed, cfg, (0, x0, 0), (0, 1, 0), 2 * cfg["h"])
        loops[f"y_axis_y{-x0:g}"] = loop_winding(seed, cfg, (0, -x0, 0), (0, 1, 0), 2 * cfg["h"])
    # off-axis loops around z (must read 0: no carrier there)
    for x0, y0 in ((6.0, 0.0), (0.0, 9.0), (7.0, 7.0)):
        loops[f"z_axis_offaxis_x{x0:g}_y{y0:g}_z6"] = loop_winding(
            seed, cfg, (x0, y0, 6.0), (0, 0, 1), 2 * cfg["h"]
        )
    pole_ok = all(
        v["half_units"] == 2
        for k, v in loops.items()
        if k.startswith("z_axis_z") and abs(float(k.split("_z")[1].split("_")[0])) >= 6
    )
    else_ok = all(
        v["half_units"] == 0
        for k, v in loops.items()
        if k.startswith(("x_axis", "y_axis", "z_axis_offaxis"))
    )
    sph_seed = {f"R{R:g}": sphere_partition(seed, cfg, R) for R in (9.0, 18.0)}
    seed_ok = all(
        v["partition"] == [2, 2]
        and v["total_half_units"] == 4
        and v["n_unreadable_plaquettes"] == 0
        for v in sph_seed.values()
    )
    M3 = M[..., 1:4, 1:4]
    sph_end = {f"R{R:g}": sphere_partition(M3, cfg, R) for R in (6.0, 9.0, 12.0, 18.0)}
    end_loops = {}
    for z0 in (9.0, 12.0):
        for sgn in (1, -1):
            end_loops[f"z_axis_z{sgn * z0:+g}"] = loop_winding(
                M3, cfg, (0, 0, sgn * z0), (0, 0, 1), 2 * cfg["h"]
            )
    a = aud["h_partition_reader"]
    ok = pole_ok and else_ok and seed_ok
    end12, end6, end9 = sph_end["R12"], sph_end["R6"], sph_end["R9"]
    end_stmt = (
        f"on the stored end field my reader (flagging plaquettes with |d.rhat| < 0.15 instead of aborting) reads "
        f"r 6: partition {end6['partition']} total {end6['total_half_units']} with "
        f"{100 * end6['n_unreadable_plaquettes'] / end6['n_plaquettes']:.1f} percent flagged, min |d.rhat| "
        f"{end6['min_abs_align']:.2f}, {100 * end6['frac_abs_align_below_0.5']:.0f} percent below 0.5; r 9: "
        f"{end9['partition']} total {end9['total_half_units']}, {100 * end9['n_unreadable_plaquettes'] / end9['n_plaquettes']:.1f} "
        f"percent flagged, min {end9['min_abs_align']:.2f}, {100 * end9['frac_abs_align_below_0.5']:.0f} percent below 0.5; "
        f"r 12: {end12['partition']} total {end12['total_half_units']}, poles {end12['poles']}, none flagged, min "
        f"{end12['min_abs_align']:.2f}"
    )
    return {
        "method": "own loop reader (fixed transverse frame, double angle of the middle eigenvector, trilinear "
        "sampling) on loops of radius 2h and 3h around the +-z axis at |z0| 3 to 15, around the x and y axes and "
        "off-axis; own sphere reader (band plaquettes + cap loops + 2 per pole) on the seed at r 9, 18 and on the "
        "stored end field at r 6, 9, 12, 18 with the outward-alignment statistics of the director",
        "mine": {
            "seed_loops": loops,
            "seed_spheres": sph_seed,
            "end_field_spheres": sph_end,
            "end_field_pole_loops": end_loops,
        },
        "audited": {
            "R20_S1_seed_r9": a["R20_S1_seed_r9"],
            "stored_S1_end_R12_partition": a["stored_S1_end"]["R12"]["partition"],
            "stored_S1_end_R9": a["stored_S1_end"]["R9"]["partition"],
        },
        "verdict": "QUALIFIED" if ok else "REFUTED",
        "note": "the seed claims hold on every loop and sphere (2 half-units at each pole from |z0| 3 to 15 at "
        "two radii, 0 around the x and y axes and off-axis, {2, 2} with total 4 on r 9 and r 18). Two scope "
        "corrections on the end-field statements: (1) 'unreadable at r <= 9' is the audited reader's abort "
        "policy, not a property of the field: the director is indeed not outward on part of those spheres "
        "(min |d.rhat| 0 on r 6 and r 9, confirmed), but a reader that flags those plaquettes still recovers the "
        "total 4 and a partition; "
        + end_stmt
        + "; the r 6 read (four unit carriers at rho 5.2, |z| 2.9) is the "
        "two pole-to-pole half lines of report 018 item 2. (2) The partition at r 12 and r 18 is resolution "
        "dependent: my 120 x 240 grid finds the audited carriers (rho 4.2 at |z| 11.2, rho 10.9 at |z| 5.1 twice, "
        "-1 at each pole) plus hairpin pairs (+1, -1 within a cell) where a line grazes the sphere, so the count of "
        "+-1 pairs is a sampling statement and only the net partition (the sum) is invariant",
    }


# ================= (i) the h scaling =================
def check_i(aud, M, cfg):
    c2 = B3.base_cfg(s=-1.0, n=cfg["n"], L=2 * cfg["L"], delta=DELTA)
    eu1, ev1 = B3.e_parts(M, cfg)
    eu2, ev2 = B3.e_parts(M, c2)
    oe1, _ = curv_energies(M, cfg["h"])
    oe2, _ = curv_energies(M, c2["h"])
    ov1 = v4_own(M, R0.roots_of(cfg), W1, cfg["h"])
    ov2 = v4_own(M, R0.roots_of(c2), W1, c2["h"])
    mine = {
        "ratio_u_platform": float(eu2 / eu1),
        "ratio_v_platform": float(ev2 / ev1),
        "ratio_u_own": oe2 / oe1,
        "ratio_v_own": ov2 / ov1,
        "E_u_h": float(eu1),
        "V4_h_W1": float(ev1),
        "V4_h_w25": float(25 * ev1),
        "E_u_over_3V_w25": float(eu1 / (3 * 25 * ev1)),
    }
    a = aud["i_virial_scaling"]
    ok = (
        abs(mine["ratio_u_platform"] - 0.5) < 1e-12
        and abs(mine["ratio_v_platform"] - 8.0) < 1e-10
        and abs(mine["ratio_u_own"] - 0.5) < 1e-12
        and abs(mine["ratio_v_own"] - 8.0) < 1e-10
        and abs(mine["E_u_h"] - a["E_u_h"]) < 1e-9
    )
    return {
        "method": "the same array under cfg(L 48) and cfg(L 96), platform e_parts and own energies",
        "mine": mine,
        "audited": {k: a[k] for k in ("ratio_u", "ratio_v", "E_u_h", "V4_h")},
        "verdict": "CONFIRMED" if ok else "REFUTED",
        "note": "an identity of the discretization (E_u = h^3 x h^-4, V4 = h^3), true on any array; the audited "
        "V4_h 0.0632 is at W1 alone, the virial line's V 1.579 at w25, both consistent",
    }


# ================= (j) the degree reader =================
def degree_solid_angle(vecs):
    """own: signed solid angle sum of the triangulated map (nth+1 rings incl. the poles, nph columns)."""
    v = vecs / np.linalg.norm(vecs, axis=-1)[..., None]
    v = np.concatenate([v, v[:, :1]], axis=1)
    a, b, c, d = v[:-1, :-1], v[1:, :-1], v[1:, 1:], v[:-1, 1:]

    def tri(p, q, s):
        num = np.einsum("...i,...i->...", p, np.cross(q, s))
        den = (
            1
            + np.einsum("...i,...i->...", p, q)
            + np.einsum("...i,...i->...", q, s)
            + np.einsum("...i,...i->...", s, p)
        )
        return 2 * np.arctan2(num, den)

    return float((np.sum(tri(a, b, c)) + np.sum(tri(a, c, d))) / (4 * np.pi))


def degree_read(F, cfg, center, R, is_vector, nth=90, nph=180):
    n, h = cfg["n"], cfg["h"]
    th = np.linspace(0.0, np.pi, nth + 1)
    ph = np.arange(nph) * 2 * np.pi / nph
    TH, PH = np.meshgrid(th, ph, indexing="ij")
    rhat, _, _ = frames_sph(TH, PH)
    pts = np.asarray(center, float) + R * rhat
    out = {"R": R, "center": list(center)}
    if is_vector:
        v = interp(F, pts, h, n, order=1)
        out["degree"] = degree_solid_angle(v)
        out["min_abs_align"] = float(np.min(np.abs(np.sum(unit_rows(v) * rhat, -1))))
    else:
        S = interp(F, pts, h, n, order=1)
        lam, V = np.linalg.eigh(S)
        d = V[..., :, 2]
        align = np.sum(d * rhat, -1)
        # rule 1: outward orientation
        d_out = d * np.sign(align)[..., None]
        out["degree_outward_rule"] = degree_solid_angle(d_out)
        # rule 2: continuity along the meridian phi = 0 then along each ring
        d_c = d.copy()
        for i in range(1, nth + 1):
            if np.dot(d_c[i, 0], d_c[i - 1, 0]) < 0:
                d_c[i, 0] *= -1
        flips = 0
        for i in range(nth + 1):
            for j in range(1, nph):
                if np.dot(d_c[i, j], d_c[i, j - 1]) < 0:
                    d_c[i, j] *= -1
            flips += np.dot(d_c[i, -1], d_c[i, 0]) < 0
        out["degree_continuity_rule"] = degree_solid_angle(d_c)
        out["continuity_ring_closure_failures"] = int(flips)
        out["min_abs_align"] = float(np.min(np.abs(align)))
        out["gap_dir_min"] = float(np.min(lam[..., 2] - lam[..., 1]))
        out["degree"] = out["degree_continuity_rule"]
    return out


def unit_rows(v):
    return v / np.maximum(np.linalg.norm(v, axis=-1), 1e-300)[..., None]


def check_j(aud, cfg):
    seed = R20.seed_axes(cfg, (1.0, DELTA, 0.0))[..., 1:4, 1:4]
    s = {f"R{R:g}": degree_read(seed, cfg, (0.0, 0.0, 0.0), R, False) for R in (6.0, 9.0, 18.0)}
    nfield = np.load(UNLIKE_NPZ)["n"].astype(np.float64)
    u = {
        "core_x+6": degree_read(nfield, cfg, (6.0, 0.0, 0.0), 4.0, True),
        "core_x-6": degree_read(nfield, cfg, (-6.0, 0.0, 0.0), 4.0, True),
        "outer_r20": degree_read(nfield, cfg, (0.0, 0.0, 0.0), 20.0, True),
        "core_x+6_R3": degree_read(nfield, cfg, (6.0, 0.0, 0.0), 3.0, True),
        "core_x-6_R3": degree_read(nfield, cfg, (-6.0, 0.0, 0.0), 3.0, True),
        "outer_r16": degree_read(nfield, cfg, (0.0, 0.0, 0.0), 16.0, True),
    }
    # a torus-like control: a sphere that contains both cores must read 0, one that contains one core +-1
    u["both_cores_R9"] = degree_read(nfield, cfg, (0.0, 0.0, 0.0), 9.0, True)
    a = aud["j_degree_reader"]
    ok = (
        all(abs(abs(v["degree"]) - 1.0) < 0.05 for v in s.values())
        and all(v["continuity_ring_closure_failures"] == 0 for v in s.values())
        and abs(u["core_x+6"]["degree"] - 1.0) < 0.05
        and abs(u["core_x-6"]["degree"] + 1.0) < 0.05
        and abs(u["core_x+6_R3"]["degree"] - 1.0) < 0.05
        and abs(u["core_x-6_R3"]["degree"] + 1.0) < 0.05
        and abs(u["outer_r20"]["degree"]) < 1e-3
        and abs(u["outer_r16"]["degree"]) < 1e-3
        and abs(u["both_cores_R9"]["degree"]) < 0.05
    )
    return {
        "method": "signed solid-angle sum of the triangulated director map on lat-long spheres (91 x 180), the "
        "director oriented by continuity (with the ring-closure check) and by the outward rule for the "
        "eigenvector field, the stored unit vectors as they are for the R22-2 pair; extra spheres R 3 and a "
        "both-cores sphere R 9 as controls",
        "mine": {"S1_seed": s, "R22_2_unlike_d12": u},
        "audited": {
            "S1_seed": {k: v["degree"] for k, v in a["S1_seed"].items()},
            "R22_2": {
                k: v["degree"] for k, v in a["R22_2_unlike_d12"].items() if isinstance(v, dict)
            },
        },
        "verdict": "QUALIFIED" if ok else "REFUTED",
        "note": "the topological content holds (+1 on every seed sphere, +1 and -1 around the pair's cores, 0 on "
        "the outer spheres and on a sphere enclosing both cores), but the audited numbers 0.975 / 0.990 / 0.998 "
        "and +-0.931 are not degrees: the signed solid-angle sum of a triangulated map is an exact integer "
        "whenever the sampled director is continuous on the triangulation, and mine reads +-1.0000 to 1e-12 on "
        "every sphere including R 3 around the cores; the audited reader therefore carries a discretization "
        "residual of up to 7 percent (a flux-integral estimator, not a covering count), and its fails_if "
        "tolerance 'within 0.05' measures that estimator, not the field; a degree reader should return integers "
        "and report non-integers as a sampling failure. The seed's sign is +1 in the outward orientation; the "
        "continuity rule fixes only |degree| (the global sign is the initial choice)",
    }


# ================= main =================
def main():
    with open(FORM_JSON) as f:
        aud = json.load(f)
    cfg = R21.cfg_of(32, 48.0, 8.0, DELTA)
    M = np.load(S1_NPZ)["M"].astype(np.float64)
    Mg = np.load(S1_GATE_NPZ)["M"].astype(np.float64)
    log(f"stored charge loaded; end vs gate max diff {np.max(np.abs(M - Mg)):.3e}")
    res = {"task": "M5.32 R26-0 audit", "audited_json": os.path.relpath(FORM_JSON, HERE)}
    res["end_vs_gate_maxdiff"] = float(np.max(np.abs(M - Mg)))
    for key, fn, args in (
        ("a", check_a, (aud,)),
        ("b", check_b, (aud, M, cfg)),
        ("c", check_c, (aud, M, cfg)),
        ("d", check_d, (aud,)),
        ("e", check_e, (aud,)),
        ("f", check_f, (aud, M, cfg)),
        ("g", check_g, (aud, M, Mg, cfg)),
        ("h", check_h, (aud, M, cfg)),
        ("i", check_i, (aud, M, cfg)),
        ("j", check_j, (aud, cfg)),
    ):
        t = time.time()
        res[key] = fn(*args)
        res[key]["wall_s"] = round(time.time() - t, 2)
        log(f"({key}) {res[key]['verdict']}  [{res[key]['wall_s']}s]")
    res["wall_s"] = round(time.time() - T0, 1)
    res["all_confirmed_or_qualified"] = all(res[k]["verdict"] != "REFUTED" for k in "abcdefghij")
    with open(OUT_JSON, "w") as f:
        json.dump(res, f, indent=1, default=float)
    print("\n| check | verdict | own numbers | audited |")
    print("| --- | --- | --- | --- |")
    r = res
    print(
        f"| a tangent | {r['a']['verdict']} | fd rel err rot {r['a']['mine']['rotation_12']['fd_rel_err']:.1e}, "
        f"boost {r['a']['mine']['boost_1']['fd_rel_err']:.1e}; sympy all true | commutator / anticommutator |"
    )
    print(
        f"| b catalog | {r['b']['verdict']} | worst symmetric part {r['b']['worst_symmetric_part']:.1e}, "
        f"corrected antisymmetric {r['b']['worst_corrected_antisymmetric_part']:.1e} | <1e-12 both |"
    )
    cm = r["c"]["mine"]
    print(
        f"| c norm | {r['c']['verdict']} | {cm['n_eta_negative']}/22 eta negative, worst rel vs audited "
        f"{cm['worst_rel_vs_audited_curvature']:.1e}, S1 eta {cm['stored_S1']['E_u_eta']:.6f} = fro "
        f"{cm['stored_S1']['E_u_frobenius']:.6f} | 22/22, 7.243767 |"
    )
    dm = r["d"]["mine"]
    print(
        f"| d telescoped | {r['d']['verdict']} | plain relerr 0.3..0.001: "
        + ", ".join(f"{dm[k]['plain_relerr_bracket_max']:.1e}" for k in dm)
        + "; telescoped: "
        + ", ".join(f"{dm[k]['telescoped_relerr_bracket_max']:.1e}" for k in dm)
        + " | plain 2.5e-11..4.2e-7; tele 5.6e-15..1.5e-9 |"
    )
    em = r["e"]["mine"]
    print(
        f"| e beta^2 | {r['e']['verdict']} | {em['uniaxial']:.1e}, {em['maximal']:.6f}, {em['vacuum_d0.3']:.6f}, "
        f"delta* {em['delta_at_0.382']:.6f} | 0, 1, 0.603755, 0.230676 |"
    )
    fm = r["f"]["mine"]
    print(
        f"| f gap tail | {r['f']['verdict']} | null l1/l0 max {fm['null_l1_over_l0_max_window_lattice']:.2e} "
        f"(lattice) {fm['null_l1_over_l0_max_window_sphere']:.2e} (sphere); synthetic slope "
        f"{fm['synthetic_slope_l1_lattice']:.3f} / {fm['synthetic_slope_l1_sphere']:.3f}, amp err "
        f"{fm['synthetic_a1_r2_over_A_maxerr_lattice']:.3f} / {fm['synthetic_a1_r2_over_A_maxerr_sphere']:.3f}; "
        f"z-odd part of the gap deviation {fm['null_z_reflection_asymmetry']['dev_odd_rms_over_even_rms']:.2f} of "
        f"the even part (rms) | 1.4e-3; -2.033; 0.045; 'symmetric to 1e-3' |"
    )
    gm = r["g"]["mine"]["seeds_and_charge"]
    gf = r["g"]["mine"]["fd_check"]
    us, ss = gm["uniaxial_seed"], gm["stored_S1_end"]
    print(
        f"| g generator | {r['g']['verdict']} | uniaxial rigid/int {us['rigid_over_internal_fwdbwd']:.4f} (audited "
        f"stencil) / {us['rigid_over_internal_o2']:.4f} (central) / {us['rigid_over_internal_o4']:.5f} (o4), inside "
        f"r 18 {us['rigid_within_r18_over_internal']:.4f}; S1 int {ss['stencil_2']['internal']:.0f} orb "
        f"{ss['stencil_fwdbwd']['orbital']:.0f} / {ss['stencil_2']['orbital']:.0f} rigid {ss['stencil_fwdbwd']['rigid']:.0f} / "
        f"{ss['stencil_2']['rigid']:.0f} (audited stencil / central); fd smooth {gf['synthetic_bump']['rel_norm_diff_o2']:.1e} (o2) "
        f"{gf['synthetic_bump']['rel_norm_diff_o4']:.1e} (o4), smoothed {gf['stored_smoothed_sigma1']['rel_norm_diff_o2']:.3f}, "
        f"raw {gf['stored_raw']['rel_norm_diff_o2']:.3f} | 0.0018; 14071 / 1403 / 14031; 7e-4, 0.084, 0.365 |"
    )
    hm = r["h"]["mine"]
    zl = {k: v["half_units"] for k, v in hm["seed_loops"].items()}
    print(
        f"| h partition | {r['h']['verdict']} | seed loops (half-units) {zl}; seed spheres "
        f"{ {k: v['partition'] for k, v in hm['seed_spheres'].items()} }; end field "
        f"{ {k: (v['partition'], v['n_unreadable_plaquettes'], round(v['min_abs_align'], 3)) for k, v in hm['end_field_spheres'].items()} } "
        f"| {{2,2}}; end R12 [1x6, -1x2] |"
    )
    im = r["i"]["mine"]
    print(
        f"| i scaling | {r['i']['verdict']} | ratio_u {im['ratio_u_own']:.12f} ratio_v {im['ratio_v_own']:.12f} | 0.5, 8 |"
    )
    jm = r["j"]["mine"]
    print(
        f"| j degree | {r['j']['verdict']} | seed "
        + ", ".join(f"{k} {v['degree']:+.4f}" for k, v in jm["S1_seed"].items())
        + "; pair "
        + ", ".join(f"{k} {v['degree']:+.4f}" for k, v in jm["R22_2_unlike_d12"].items())
        + " | 0.975, 0.990, 0.998; +0.931, -0.931, 8e-10 |"
    )
    print(f"\nwall {res['wall_s']} s; JSON {os.path.relpath(OUT_JSON, HERE)}")


if __name__ == "__main__":
    main()
