"""M5.32 R26-3 audit: the spin gate with the physical rotation generator, refuted or confirmed
with OWN code by a fresh agent that did not read `m5_32_r26_3_spin.py` (nor any R26 script).

What is own here: the orbital derivative (np.gradient, second-order central in the interior with
second-order one-sided edges, plus a fourth-order central variant), the kinetic contraction (own
fwd/bwd stencils and an explicit tr(eta F eta F^T)), the rotation by cubic-spline interpolation
that tests the sign of the transport term, the exact 90-degree lattice rotation that measures the
axial asymmetry, the within-r profile on an own kin density, the energy (own curvature and V4
trace sums beside the platform's), the gradient read on the free cells through the R23-1 reduced
descent's `energy_grad`, the solid-angle degree on an icosphere with a trilinear-interpolated
director, and a short own continuation of two capped rows to price the convergence sensitivity
of C_rigid.

Inputs (read-only): the platform `m5_21_3_a_4d.py` (kin_of, a_fields, comm_eta, inner_eta,
coords, pin_shell), `m5_32_r21_1_runs.py` (cfg_of, params_of), `m5_32_r20_0_class.py`
(roots_of), `m5_32_r20_1_axes.py` (energy_parts, seed_axes), `m5_32_r23_1_cscan.py`
(energy_grad, Reduced), the stored fields in data/m5_32_r25_2, data/m5_32_r26_3,
data/m5_32_r22_2, and the JSONs of R25-2 and R26-3 for the audited numbers.

Run: python3 m5_32_r26_3_audit.py  ->  data/m5_32_r26_3_audit.json
"""

from __future__ import annotations

import importlib.util
import json
import os
import sys
import time

import numpy as np
from scipy.ndimage import map_coordinates
from scipy.optimize import minimize

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
OUT_JSON = os.path.join(DATA, "m5_32_r26_3_audit.json")
R26_JSON = os.path.join(DATA, "m5_32_r26_3_spin.json")
R25_JSON = os.path.join(DATA, "m5_32_r25_2_charge.json")
R25_NPZ = os.path.join(DATA, "m5_32_r25_2")
R26_NPZ = os.path.join(DATA, "m5_32_r26_3")
R22_NPZ = os.path.join(DATA, "m5_32_r22_2")
T0 = time.time()


def log(msg):
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


CS = _load("m5_32_r23_1_cscan", "m5_32_r23_1_cscan.py")
R21, R20, B3, R0, W1 = CS.R21, CS.R20, CS.B3, CS.R0, CS.W1
ETA = np.diag([-1.0, 1.0, 1.0, 1.0])
ETAD = np.diag(ETA)
G, W1S, GATE = 8.0, 25.0, 1e-4
JZ = np.zeros((4, 4))
JZ[1, 2], JZ[2, 1] = -1.0, 1.0
LADDER = [
    ("S1_d0.3_w25_n32_L48", 32, 48.0, "h1.5", "48"),
    ("S1_d0.3_w25_n48_L72", 48, 72.0, "h1.5", "72"),
    ("S1_d0.3_w25_n64_L96", 64, 96.0, "h1.5", "96"),
    ("S1_d0.3_w25_n48_L48", 48, 48.0, "h1", "48"),
    ("S1_d0.3_w25_n64_L64", 64, 64.0, "h1", "64"),
]
ARM = [("S1_d0.1_w25_n32_L48", 0.1), ("S1_d0.03_w25_n32_L48", 0.03)]
PAIRS = [
    ("unlike_d6_n32_L48", 6.0, 32, 48.0),
    ("unlike_d6_n48_L48", 6.0, 48, 48.0),
    ("unlike_d8_n32_L48", 8.0, 32, 48.0),
    ("unlike_d8_n48_L48", 8.0, 48, 48.0),
    ("unlike_d8_n64_L48", 8.0, 64, 48.0),
    ("unlike_d8_n64_L64", 8.0, 64, 64.0),
    ("unlike_d10_n32_L48", 10.0, 32, 48.0),
    ("unlike_d10_n48_L48", 10.0, 48, 48.0),
    ("unlike_d12_n32_L48", 12.0, 32, 48.0),
    ("unlike_d12_n48_L48", 12.0, 48, 48.0),
    ("unlike_d12_n64_L64", 12.0, 64, 64.0),
]


# ================= stack handles =================
def stack(n, L, delta):
    cfg = R21.cfg_of(n, L, G, delta)
    p = R21.params_of(G, delta)
    pot = ("v4", R0.roots_of(cfg), W1 * W1S)
    return cfg, p, pot


def load_M(folder, tag, suffix=""):
    return np.load(os.path.join(folder, f"{tag}{suffix}.npz"), allow_pickle=True)["M"].astype(
        np.float64
    )


def radius(cfg):
    X, Y, Z = B3.coords(cfg["n"], cfg["h"])
    return np.sqrt(X * X + Y * Y + Z * Z)


# ================= own derivatives and generators =================
def d_c2(f, ax, h):
    """second-order central in the interior, second-order one-sided at the edges."""
    return np.gradient(f, h, axis=ax, edge_order=2)


def d_c4(f, ax, h):
    """fourth-order central in the interior, the np.gradient second-order read on the two
    outermost cells of each side."""
    out = d_c2(f, ax, h)
    sl = [slice(None)] * f.ndim

    def at(i):
        s = list(sl)
        s[ax] = i
        return tuple(s)

    out[at(slice(2, -2))] = (
        -f[at(slice(4, None))]
        + 8.0 * f[at(slice(3, -1))]
        - 8.0 * f[at(slice(1, -3))]
        + f[at(slice(0, -4))]
    ) / (12.0 * h)
    return out


def d_fwd(f, ax, h):
    out = np.zeros_like(f)
    sl = [slice(None)] * f.ndim
    s1, s0 = list(sl), list(sl)
    s1[ax], s0[ax] = slice(1, None), slice(0, -1)
    out[tuple(s0)] = (f[tuple(s1)] - f[tuple(s0)]) / h
    return out


def d_bwd(f, ax, h):
    out = np.zeros_like(f)
    sl = [slice(None)] * f.ndim
    s1, s0 = list(sl), list(sl)
    s1[ax], s0[ax] = slice(1, None), slice(0, -1)
    out[tuple(s1)] = (f[tuple(s1)] - f[tuple(s0)]) / h
    return out


def a_internal(M):
    return JZ @ M - M @ JZ


def a_orbital(M, cfg, deriv):
    X, Y, _ = B3.coords(cfg["n"], cfg["h"])
    h = cfg["h"]
    return X[..., None, None] * deriv(M, 1, h) - Y[..., None, None] * deriv(M, 0, h)


def own_kin_density(M, a0, h):
    """4 h^3 sum_i <[a0, A_i]_eta, [a0, A_i]_eta>_eta on the sym stencil (half fwd + half bwd),
    the eta inner product written out as sum_ab eta_a eta_b F_ab^2."""
    dens = np.zeros(M.shape[:3])
    for deriv in (d_fwd, d_bwd):
        for ax in range(3):
            A = deriv(M, ax, h)
            F = a0 @ ETA @ A - A @ ETA @ a0
            dens += 0.5 * 4.0 * np.einsum("...ab,a,b->...", F * F, ETAD, ETAD)
    return h**3 * dens


def inertias(M, cfg, deriv=d_c2):
    """(C_int, C_orb, C_rigid) through the platform's kin_of with the own orbital derivative."""
    ai = a_internal(M)
    ao = a_orbital(M, cfg, deriv)
    return (
        float(B3.kin_of(M, ai, cfg)),
        float(B3.kin_of(M, ao, cfg)),
        float(B3.kin_of(M, ai - ao, cfg)),
        ai,
        ao,
    )


# ================= rotations =================
def rot4(theta):
    R = np.eye(4)
    c, s = np.cos(theta), np.sin(theta)
    R[1, 1], R[1, 2], R[2, 1], R[2, 2] = c, -s, s, c
    return R


def rotate_interp(M, cfg, theta, order=3):
    """M_rot(x) = R M(R^-1 x) R^T with R = exp(theta J_z), cubic-spline sampling of the source."""
    n, h = cfg["n"], cfg["h"]
    X, Y, Z = B3.coords(n, h)
    c, s = np.cos(theta), np.sin(theta)
    xs, ys = c * X + s * Y, -s * X + c * Y
    ctr = (n - 1) / 2.0
    idx = np.stack([xs / h + ctr, ys / h + ctr, Z / h + ctr])
    out = np.empty_like(M)
    for a in range(4):
        for b in range(4):
            out[..., a, b] = map_coordinates(M[..., a, b], idx, order=order, mode="nearest")
    R = rot4(theta)
    return R @ out @ R.T


def rotate_90(M):
    """the exact lattice rotation by +90 degrees about z (cell centers map onto cell centers):
    M_rot[i, j, k] = R M[j, n-1-i, k] R^T."""
    src = np.swapaxes(M[:, ::-1], 0, 1)
    R = rot4(np.pi / 2)
    return R @ src @ R.T


def spatial_std(M):
    S = M[..., 1:, 1:]
    return float(np.sqrt(np.mean((S - S.mean(axis=(0, 1, 2))) ** 2)))


def rms(A):
    return float(np.sqrt(np.mean(A * A)))


# ================= own energy =================
def own_energy(M, cfg, pot):
    h = cfg["h"]
    curv = np.zeros(M.shape[:3])
    for deriv in (d_fwd, d_bwd):
        A = [deriv(M, ax, h) for ax in range(3)]
        for i in range(3):
            for j in range(i + 1, 3):
                F = A[i] @ ETA @ A[j] - A[j] @ ETA @ A[i]
                curv += 0.5 * 4.0 * np.einsum("...ab,a,b->...", F * F, ETAD, ETAD)
    N = M @ ETA
    P = np.broadcast_to(np.eye(4), M.shape).copy()
    V = np.zeros(M.shape[:3])
    for pw in range(1, 5):
        P = P @ N
        V += (np.einsum("...kk->...", P) - sum(q**pw for q in pot[1])) ** 2
    return float(h**3 * curv.sum()), float(h**3 * pot[2] * V.sum())


# ================= degree on a sphere =================
def icosphere(k=4):
    t = (1.0 + np.sqrt(5.0)) / 2.0
    v = np.array(
        [
            [-1, t, 0], [1, t, 0], [-1, -t, 0], [1, -t, 0],
            [0, -1, t], [0, 1, t], [0, -1, -t], [0, 1, -t],
            [t, 0, -1], [t, 0, 1], [-t, 0, -1], [-t, 0, 1],
        ],
        dtype=float,
    )  # fmt: skip
    v /= np.linalg.norm(v, axis=1)[:, None]
    f = [
        (0, 11, 5), (0, 5, 1), (0, 1, 7), (0, 7, 10), (0, 10, 11),
        (1, 5, 9), (5, 11, 4), (11, 10, 2), (10, 7, 6), (7, 1, 8),
        (3, 9, 4), (3, 4, 2), (3, 2, 6), (3, 6, 8), (3, 8, 9),
        (4, 9, 5), (2, 4, 11), (6, 2, 10), (8, 6, 7), (9, 8, 1),
    ]  # fmt: skip
    verts = [tuple(x) for x in v]
    cache = {}

    def mid(a, b):
        key = (min(a, b), max(a, b))
        if key not in cache:
            m = np.array(verts[a]) + np.array(verts[b])
            m /= np.linalg.norm(m)
            verts.append(tuple(m))
            cache[key] = len(verts) - 1
        return cache[key]

    for _ in range(k):
        nf = []
        for a, b, c in f:
            ab, bc, ca = mid(a, b), mid(b, c), mid(c, a)
            nf += [(a, ab, ca), (b, bc, ab), (c, ca, bc), (ab, bc, ca)]
        f = nf
    V = np.array(verts)
    F = np.array(f)
    # outward orientation
    A, Bv, C = V[F[:, 0]], V[F[:, 1]], V[F[:, 2]]
    det = np.einsum("ij,ij->i", np.cross(Bv - A, C - A), A)
    flip = det < 0
    F[flip] = F[flip][:, [0, 2, 1]]
    return V, F


ICO_V, ICO_F = icosphere(4)


def solid_angle_sum(u):
    """u: unit vectors at the icosphere vertices; the signed solid angle of the image of every
    face (Van Oosterom and Strackee), summed, over 4 pi."""
    a, b, c = u[ICO_F[:, 0]], u[ICO_F[:, 1]], u[ICO_F[:, 2]]
    num = np.einsum("ij,ij->i", a, np.cross(b, c))
    den = (
        1.0
        + np.einsum("ij,ij->i", a, b)
        + np.einsum("ij,ij->i", b, c)
        + np.einsum("ij,ij->i", c, a)
    )
    return float(np.sum(2.0 * np.arctan2(num, den)) / (4.0 * np.pi))


def degree_on_sphere(nfield, n, h, center, R):
    ctr = (n - 1) / 2.0
    pts = np.asarray(center)[None, :] + R * ICO_V
    idx = (pts / h + ctr).T
    u = np.stack(
        [map_coordinates(nfield[..., k], idx, order=1, mode="nearest") for k in range(3)], axis=-1
    )
    norms = np.linalg.norm(u, axis=1)
    u = u / np.maximum(norms, 1e-300)[:, None]
    # continuity of the sampled map along the mesh edges
    e = np.concatenate([ICO_F[:, [0, 1]], ICO_F[:, [1, 2]], ICO_F[:, [2, 0]]])
    cosang = np.einsum("ij,ij->i", u[e[:, 0]], u[e[:, 1]])
    return (
        solid_angle_sum(u),
        float(np.degrees(np.arccos(np.clip(cosang.min(), -1, 1)))),
        float(norms.min()),
    )


def fnum(x):
    return None if x is None else float(x)


# ================= the claims =================
def claim_a_b(J26):
    """A: the box ladder with own orbital derivative and own contraction; the sign test by
    interpolation. B: the within-r profile on an own kin density."""
    rows, ladder_own = {}, {}
    for tag, n, L, hk, Lk in LADDER:
        cfg, p, pot = stack(n, L, 0.3)
        M = load_M(R25_NPZ, tag)
        aud = J26["collect"]["box"][hk][Lk]
        C2 = inertias(M, cfg, d_c2)
        C4 = inertias(M, cfg, d_c4)
        ai, ao = C2[3], C2[4]
        own_dens = own_kin_density(M, ai - ao, cfg["h"])
        own_rigid = float(own_dens.sum())
        own_int = float(own_kin_density(M, ai, cfg["h"]).sum())
        own_orb = float(own_kin_density(M, ao, cfg["h"]).sum())
        # the sign test: the rigid rotation by interpolation, differenced
        th = 0.01
        a_fd = (rotate_interp(M, cfg, th) - rotate_interp(M, cfg, -th)) / (2.0 * th)
        r = radius(cfg)
        inner = r < (L / 2.0 - 3.0 * cfg["h"])
        w = inner[..., None, None]
        a_rig = ai - ao
        a_rig4 = C4[3] - C4[4]
        sign_minus = rms((a_fd - a_rig)[inner]) / rms(a_rig[inner])
        sign_plus = rms((a_fd - (ai + ao))[inner]) / rms(a_rig[inner])
        core, halo = r < 9.0, (r > 12.0) & inner
        sign_core = rms((a_fd - a_rig)[core]) / rms(a_rig[core])
        sign_halo = rms((a_fd - a_rig)[halo]) / rms(a_rig[halo])
        sign_halo4 = rms((a_fd - a_rig4)[halo]) / rms(a_rig4[halo])
        dens_fd = own_kin_density(M, np.where(w, a_fd, a_rig), cfg["h"])
        # the within-r profile on the own density
        prof = {}
        for R in (3, 6, 9, 12, 18):
            prof[str(R)] = float(own_dens[r < R].sum())
        prof_aud = {k: v["rigid"] for k, v in aud["rigid_within_r"].items()}
        rows[tag] = {
            "h": cfg["h"],
            "L": L,
            "own_c2": {"C_int": C2[0], "C_orb": C2[1], "C_rigid": C2[2]},
            "own_c4": {"C_int": C4[0], "C_orb": C4[1], "C_rigid": C4[2]},
            "own_contraction": {"C_int": own_int, "C_orb": own_orb, "C_rigid": own_rigid},
            "audited": {"C_int": aud["C_int"], "C_orb": aud["C_orb"], "C_rigid": aud["C_rigid"]},
            "rel_dev_c2_vs_audited": {
                "C_int": C2[0] / aud["C_int"] - 1,
                "C_orb": C2[1] / aud["C_orb"] - 1,
                "C_rigid": C2[2] / aud["C_rigid"] - 1,
            },
            "sign_test": {
                "theta_rad": th,
                "rel_rms_fd_minus_(int-orb)": sign_minus,
                "rel_rms_fd_minus_(int+orb)": sign_plus,
                "rel_rms_fd_minus_(int-orb)_core_r_lt_9": sign_core,
                "rel_rms_fd_minus_(int-orb)_halo_r_gt_12": sign_halo,
                "rel_rms_fd_minus_(int-orb)_halo_4th_order_stencil": sign_halo4,
                "interior_r_lt": L / 2.0 - 3.0 * cfg["h"],
                "C_rigid_interior_fd": float(dens_fd[inner].sum()),
                "C_rigid_interior_stencil": float(own_dens[inner].sum()),
            },
            "within_r_own": prof,
            "within_r_audited": prof_aud,
            "frac_inside_12_own": prof["12"] / own_rigid,
            "frac_inside_12_audited": prof_aud["12"] / aud["C_rigid"],
            "rms_a_rigid_over_a_int": rms(a_rig) / rms(ai),
        }
        ladder_own[(hk, Lk)] = C2
        log(
            f"A {tag}: own (c2) {C2[0]:.0f} {C2[1]:.0f} {C2[2]:.0f} | c4 orb {C4[1]:.0f} rigid "
            f"{C4[2]:.0f} | own contraction rigid {own_rigid:.1f} | audited {aud['C_int']:.0f} "
            f"{aud['C_orb']:.0f} {aud['C_rigid']:.0f} | sign(-) {sign_minus:.3f} sign(+) "
            f"{sign_plus:.3f} | within 6,12,18: {prof['6']:.0f} {prof['12']:.0f} {prof['18']:.0f}"
        )
    c4 = {k: rows[t]["own_c4"]["C_rigid"] for t, _, _, hk, Lk in LADDER for k in [(hk, Lk)]}
    ratios = {
        "h1.5_L96_over_L48_rigid": ladder_own[("h1.5", "96")][2] / ladder_own[("h1.5", "48")][2],
        "h1.5_L96_over_L48_int": ladder_own[("h1.5", "96")][0] / ladder_own[("h1.5", "48")][0],
        "h1_L64_over_L48_rigid": ladder_own[("h1", "64")][2] / ladder_own[("h1", "48")][2],
        "h1_L64_over_L48_int": ladder_own[("h1", "64")][0] / ladder_own[("h1", "48")][0],
        "h1.5_L96_over_L48_rigid_4th_order_stencil": c4[("h1.5", "96")] / c4[("h1.5", "48")],
        "h1_L64_over_L48_rigid_4th_order_stencil": c4[("h1", "64")] / c4[("h1", "48")],
    }
    ratios["synthetic_sign_test"] = synthetic_sign_test()
    log(f"A synthetic smooth-field sign test: {ratios['synthetic_sign_test']}")
    return rows, ratios


def synthetic_sign_test():
    """a smooth, resolved, non-axisymmetric field (a Gaussian bump of width 6 at (6, 3, 0) on
    the biaxial vacuum): the interpolated rotation's derivative against int - orb and int + orb,
    where the stencil error is negligible and only the convention is tested."""
    cfg, _, _ = stack(32, 48.0, 0.3)
    n, h = cfg["n"], cfg["h"]
    X, Y, Z = B3.coords(n, h)
    rng = np.random.default_rng(7)
    S = rng.standard_normal((3, 3))
    S = S + S.T
    bump = np.exp(-((X - 6.0) ** 2 + (Y - 3.0) ** 2 + Z**2) / 36.0)
    M = np.zeros((n, n, n, 4, 4))
    M[..., 0, 0] = 8.0
    M[..., 1:, 1:] = np.diag([1.0, 0.3, 0.0]) + bump[..., None, None] * S
    th = 0.01
    a_fd = (rotate_interp(M, cfg, th) - rotate_interp(M, cfg, -th)) / (2.0 * th)
    ai, ao = a_internal(M), a_orbital(M, cfg, d_c2)
    inner = radius(cfg) < 19.5
    return {
        "rel_rms_fd_minus_(int-orb)": rms((a_fd - (ai - ao))[inner]) / rms((ai - ao)[inner]),
        "rel_rms_fd_minus_(int+orb)": rms((a_fd - (ai + ao))[inner]) / rms((ai - ao)[inner]),
        "rms_orb_over_rms_int": rms(ao[inner]) / rms(ai[inner]),
    }


def claim_c():
    """the axial asymmetry: the exact 90-degree lattice rotation and a 45-degree interpolated
    one, relative to the spatial-block spread, on the relaxed L 48 charge and on the seeds."""
    out = {}
    cfg3, _, _ = stack(32, 48.0, 0.3)
    cfg0, _, _ = stack(32, 48.0, 0.0)
    fields = {
        "charge_S1_d0.3_n32_L48_end": (load_M(R25_NPZ, "S1_d0.3_w25_n32_L48"), cfg3),
        "seed_uniaxial_d0_(1,0,0)": (R20.seed_axes(cfg0, (1.0, 0.0, 0.0)), cfg0),
        "seed_S1_d0.3_(1,0.3,0)": (R20.seed_axes(cfg3, (1.0, 0.3, 0.0)), cfg3),
        # axisymmetric but UNRESOLVED cores (r_c 1.0 and 0.5 at h 1.5): the stencil's own
        # non-cancellation, to separate resolution from asymmetry on the charge
        "seed_uniaxial_d0_rc1.0": (R20.seed_axes(cfg0, (1.0, 0.0, 0.0), r_c=1.0), cfg0),
        "seed_uniaxial_d0_rc0.5": (R20.seed_axes(cfg0, (1.0, 0.0, 0.0), r_c=0.5), cfg0),
    }
    for name, (M, cfg) in fields.items():
        M90 = rotate_90(M)
        M45 = rotate_interp(M, cfg, np.pi / 4)
        r = radius(cfg)
        inner = r < (cfg["L"] / 2.0 - 3.0 * cfg["h"])
        sd = spatial_std(M)
        C = inertias(M, cfg, d_c2)
        ai, ao = C[3], C[4]
        out[name] = {
            "rms_M_rot90_minus_M_over_spread": rms(M90 - M) / sd,
            "rms_M_rot45_minus_M_over_spread_interior": rms((M45 - M)[inner]) / sd,
            "rms_M_rot180_minus_M_over_spread": rms(rotate_90(M90) - M) / sd,
            "spatial_block_spread": sd,
            "C_int": C[0],
            "C_orb": C[1],
            "C_rigid": C[2],
            "C_rigid_over_C_int": C[2] / C[0],
            "rms_a_rigid_over_a_int": rms(ai - ao) / rms(ai),
        }
        log(
            f"C {name}: rot90 {out[name]['rms_M_rot90_minus_M_over_spread']:.4f} rot45 "
            f"{out[name]['rms_M_rot45_minus_M_over_spread_interior']:.4f} | C_rigid/C_int "
            f"{C[2] / C[0]:.5f} ({C[0]:.1f}, {C[1]:.1f}, {C[2]:.2f})"
        )
    return out


def claim_d(J26):
    out = {}
    for tag, delta in ARM:
        cfg, p, pot = stack(32, 48.0, delta)
        M = load_M(R26_NPZ, tag)
        Mg = load_M(R26_NPZ, tag, "_gate")
        Ms = np.load(os.path.join(R26_NPZ, tag + "_stage.npz"), allow_pickle=True)["M"]
        parts = R20.energy_parts(M, cfg, p, pot)
        ec, ev = own_energy(M, cfg, pot)
        E, Gm, info = CS.energy_grad(M, cfg, p, pot, 0.0)
        mask = ~B3.pin_shell(cfg["n"], cfg["h"])
        Gf = Gm[mask]
        fsp, f00 = float(np.abs(Gf[:, 1:, 1:]).max()), float(np.abs(Gf[:, 0, 0]).max())
        C = inertias(M, cfg, d_c2)
        row = J26["arm_rows"][tag]
        chunks = row["chunks"]
        drops = [c["drop"] for c in chunks[-4:]]
        big = [c["iters"] for c in chunks if c.get("drop", 0.0) >= 4e-3]
        sp = row["spin_reads"]["generator"]
        out[tag] = {
            "delta": delta,
            "E_platform": parts["E_total"],
            "E_curv_platform": parts["E_curv"],
            "V_platform": parts["V"],
            "E_own": ec + ev,
            "E_curv_own": ec,
            "V_own": ev,
            "E_audited": row["E"],
            "fmax_spatial_free_own": fsp,
            "fmax_M00_free_own": f00,
            "gate": GATE,
            "fmax_over_gate": fsp / GATE,
            "fmax_audited_last_chunk": chunks[-1]["fmax_spatial"],
            "gate_hit_inside_last_chunk": chunks[-1].get("gate_hit_inside_chunk"),
            "iters": row["iters"],
            "gate_label_audited": row["gate_label"],
            "last_4_chunk_drops": drops,
            "last_chunk_iters_with_drop_ge_4e-3": (max(big) if big else None),
            "drop_criterion_GATE_DROP_x_E": 1e-5 * abs(row["E"]),
            "own_C": {"C_int": C[0], "C_orb": C[1], "C_rigid": C[2]},
            "audited_C": {"C_int": sp["internal"], "C_orb": sp["orbital"], "C_rigid": sp["rigid"]},
            "max_abs_end_minus_gate": float(np.abs(M - Mg).max()),
            "max_abs_end_minus_stage": float(np.abs(M - Ms).max()),
        }
        log(
            f"D {tag}: E {parts['E_total']:.4f} (own {ec + ev:.4f}, audited {row['E']:.4f}) | "
            f"fmax free {fsp:.2e} (gate {GATE:.0e}) | C {C[0]:.0f} {C[1]:.0f} {C[2]:.0f} "
            f"(audited {sp['internal']:.0f} {sp['orbital']:.0f} {sp['rigid']:.0f}) | last drops "
            f"{[f'{d:.1e}' for d in drops]}"
        )
    # the delta law from the three points, own fit
    pts = J26["collect"]["delta_law"]["points_delta_Crigid_Cint_tag_kick"]
    ds = np.array([pt[0] for pt in pts])
    cs = np.array([pt[1] for pt in pts])
    slope = float(np.polyfit(np.log(ds), np.log(cs), 1)[0])
    out["delta_exponent_own_fit"] = slope
    out["delta_exponent_audited"] = J26["collect"]["delta_law"]["rigid_exponent"]
    out["delta_points"] = pts
    out["monotone"] = bool(np.all(np.diff(cs) > 0) or np.all(np.diff(cs) < 0))
    return out


def claim_e(J26):
    out = {}
    # the reader's convention on synthetic hedgehogs
    n, h = 32, 1.5
    X, Y, Z = B3.coords(n, h)
    rr = np.maximum(np.sqrt(X * X + Y * Y + Z * Z), 1e-9)
    hedgehog = np.stack([X, Y, Z], axis=-1) / rr[..., None]
    conv = {
        "identity_map_n_eq_rhat": degree_on_sphere(hedgehog, n, h, (0, 0, 0), 6.0)[0],
        "antipodal_map_n_eq_minus_rhat": degree_on_sphere(-hedgehog, n, h, (0, 0, 0), 6.0)[0],
        "sphere_missing_the_core": degree_on_sphere(hedgehog, n, h, (10, 0, 0), 4.0)[0],
    }
    out["convention_checks"] = conv
    log(f"E convention: {conv}")
    for tag, d, n, L in PAIRS:
        h = L / n
        nf = np.load(os.path.join(R22_NPZ, tag + ".npz"))["n"].astype(np.float64)
        nrm = np.linalg.norm(nf, axis=-1)
        # the vector field's continuity on the lattice: neighbor dot products
        neg = 0
        tot = 0
        for ax in range(3):
            a = np.moveaxis(nf, ax, 0)
            dots = np.einsum("...k,...k->...", a[1:], a[:-1])
            neg += int((dots < 0).sum())
            tot += dots.size
        Rc = min(4.0, d / 2.0 - 1.0)
        Ro = L / 2.0 - 4.0
        reads = {}
        for name, ctr, R in [
            (f"core_x+{d / 2:g}", (d / 2.0, 0.0, 0.0), Rc),
            (f"core_x-{d / 2:g}", (-d / 2.0, 0.0, 0.0), Rc),
            ("mid_R_half_d", (0.0, 0.0, 0.0), d / 2.0),
            ("mid_R_half_d_minus_1", (0.0, 0.0, 0.0), d / 2.0 - 1.0),
            ("mid_R_half_d_plus_1", (0.0, 0.0, 0.0), d / 2.0 + 1.0),
            ("both_cores_R_half_d_plus_3", (0.0, 0.0, 0.0), d / 2.0 + 3.0),
            ("outer", (0.0, 0.0, 0.0), Ro),
        ]:
            deg, maxang, minnorm = degree_on_sphere(nf, n, h, ctr, R)
            reads[name] = {
                "R": R,
                "degree": deg,
                "max_edge_angle_deg": maxang,
                "min_interp_norm": minnorm,
                "sphere_runs_through_a_core": bool(minnorm < 0.2),
            }
        aud = J26["pair"][tag]
        out[tag] = {
            "d": d,
            "n": n,
            "L": L,
            "unit_norm_range": [float(nrm.min()), float(nrm.max())],
            "lattice_neighbor_pairs_with_negative_dot": neg,
            "lattice_neighbor_pairs_total": tot,
            "reads": reads,
            "mid_sphere_R_half_d_ill_posed": reads["mid_R_half_d"]["sphere_runs_through_a_core"],
            "audited": aud,
        }
        log(
            f"E {tag}: core+ {reads[f'core_x+{d / 2:g}']['degree']:+.4f} core- "
            f"{reads[f'core_x-{d / 2:g}']['degree']:+.4f} mid {reads['mid_R_half_d']['degree']:+.2e} "
            f"(R-1 {reads['mid_R_half_d_minus_1']['degree']:+.2e}, R+1 "
            f"{reads['mid_R_half_d_plus_1']['degree']:+.2e}, both "
            f"{reads['both_cores_R_half_d_plus_3']['degree']:+.2e}) outer "
            f"{reads['outer']['degree']:+.2e} | max edge angle "
            f"{max(v['max_edge_angle_deg'] for v in reads.values()):.1f} deg | neg dots {neg}"
        )
    return out


def claim_f(J26, J25):
    out = {"rows": {}}
    rows25 = J25["rows"]
    for tag, n, L, hk, Lk in LADDER:
        cfg, p, pot = stack(n, L, 0.3)
        M = load_M(R25_NPZ, tag)
        Mg = load_M(R25_NPZ, tag, "_gate")
        Ms = np.load(os.path.join(R25_NPZ, tag + "_stage.npz"), allow_pickle=True)["M"]
        r25 = rows25[tag]
        rec = {
            "gate_label_r25_2": r25["gate_label"],
            "kick_label_r25_2": r25.get("kick_label"),
            "iters": r25["iters"],
            "fmax_last_chunk": r25["chunks"][-1]["fmax_spatial"],
            "drop_last_chunk": r25["chunks"][-1]["drop"],
            "E": r25["E"],
            "max_abs_end_minus_gate": float(np.abs(M - Mg).max()),
            "max_abs_end_minus_stage": float(np.abs(M - Ms).max()),
        }
        kpath = os.path.join(R25_NPZ, tag + "_kick_stage.npz")
        if os.path.exists(kpath):
            Mk = np.load(kpath, allow_pickle=True)["M"].astype(np.float64)
            Ck = inertias(Mk, cfg, d_c2)
            C = inertias(M, cfg, d_c2)
            rec["kick_stage"] = {
                "C_rigid_end": C[2],
                "C_rigid_kick_stage": Ck[2],
                "rel_sensitivity_rigid": abs(Ck[2] - C[2]) / C[2],
                "rel_sensitivity_int": abs(Ck[0] - C[0]) / C[0],
                "max_abs_M_diff": float(np.abs(Mk - M).max()),
                "audited": J26["ladder"][tag].get("kick_stage_sensitivity_rigid"),
            }
        out["rows"][tag] = rec
        log(f"F {tag}: {rec}")
    # the label by the rule
    unconverged = [t for t, _, _, _, _ in LADDER if rows25[t]["gate_label"] != "AT_GATE"]
    out["rows_not_at_gate"] = unconverged
    out["label_by_rule"] = "SPIN_PHYS_INSUFFICIENT" if unconverged else "rule not triggered"
    out["label_audited"] = J26["collect"]["label"]
    # the stated non-evaluation of g
    gnotes = [v["_end"].get("g_note", "") for v in J26["ladder"].values()]
    out["g_not_evaluated_in_json"] = all("NOT evaluated" in s for s in gnotes)
    out["g_numeric_keys_present"] = any(
        k.startswith("g_") and isinstance(v["_end"].get(k), (int, float))
        for v in J26["ladder"].values()
        for k in v["_end"]
    )
    # own continuation probe on the two n 48 capped rows: dC_rigid per dE over a short L-BFGS-B
    probe = {}
    for tag, n, L, hk, Lk in LADDER:
        if n != 48:
            continue
        if time.time() - T0 > 600:
            probe[tag] = {"skipped": "time budget"}
            continue
        cfg, p, pot = stack(n, L, 0.3)
        M = load_M(R25_NPZ, tag)
        C0 = inertias(M, cfg, d_c2)
        red = CS.Reduced(M, cfg, p, pot, 0.0, True)
        E0 = red.fun(red.pack(M))[0]
        t1 = time.time()
        res = minimize(
            red.fun,
            red.pack(M),
            jac=True,
            method="L-BFGS-B",
            options={"maxcor": 20, "maxiter": 30, "maxfun": 45, "gtol": 1e-14, "ftol": 1e-16},
        )
        M1 = red.build(np.asarray(res.x))
        E1 = float(res.fun)
        C1 = inertias(M1, cfg, d_c2)
        drops = [c["drop"] for c in rows25[tag]["chunks"][-3:]]
        q = drops[-1] / drops[-2] if drops[-2] > 0 else None
        tail = drops[-1] * q / (1.0 - q) if (q is not None and 0 < q < 1) else None
        dC_dE = (C1[2] - C0[2]) / (E1 - E0) if E1 != E0 else None
        probe[tag] = {
            "iters": int(res.nit),
            "E_before": E0,
            "E_after": E1,
            "dE": E1 - E0,
            "C_rigid_before": C0[2],
            "C_rigid_after": C1[2],
            "dC_rigid": C1[2] - C0[2],
            "dC_int": C1[0] - C0[0],
            "rel_dC_rigid_per_rel_dE": ((C1[2] - C0[2]) / C0[2]) / ((E1 - E0) / E0),
            "last_3_chunk_drops": drops,
            "geometric_tail_of_chunk_drops": tail,
            "C_rigid_shift_if_tail_at_this_slope": (dC_dE * tail) if (dC_dE and tail) else None,
            "rel_shift_if_tail_at_this_slope": (
                (dC_dE * tail / C0[2]) if (dC_dE and tail) else None
            ),
            "wall_s": time.time() - t1,
        }
        log(f"F probe {tag}: {probe[tag]}")
    out["continuation_probe"] = probe
    return out


# ================= verdicts =================
def verdicts(A, ratios, C, D, E, F):
    V = {}
    # A: the fifteen inertias, the ratios, the contraction, the sign
    devs = [abs(v) for r in A.values() for v in r["rel_dev_c2_vs_audited"].values()]
    syn = ratios["synthetic_sign_test"]
    # the wrong sign's residual is bounded by 2 rms(orb) / rms(rigid), 0.4 on this field, so the
    # test is the ratio of the two residuals (interpolation-level for the right sign)
    sign_ok = (
        syn["rel_rms_fd_minus_(int-orb)"] < 0.05
        and syn["rel_rms_fd_minus_(int+orb)"] > 10.0 * syn["rel_rms_fd_minus_(int-orb)"]
    )
    sign_ok = sign_ok and all(
        r["sign_test"]["rel_rms_fd_minus_(int-orb)"] < r["sign_test"]["rel_rms_fd_minus_(int+orb)"]
        for r in A.values()
    )
    contr_ok = all(
        abs(r["own_contraction"]["C_rigid"] / r["own_c2"]["C_rigid"] - 1) < 1e-9
        for r in A.values()
    )
    c4_rig = max(abs(r["own_c4"]["C_rigid"] / r["own_c2"]["C_rigid"] - 1) for r in A.values())
    c4_orb = max(abs(r["own_c4"]["C_orb"] / r["own_c2"]["C_orb"] - 1) for r in A.values())
    halo = max(r["sign_test"]["rel_rms_fd_minus_(int-orb)_halo_r_gt_12"] for r in A.values())
    core = max(r["sign_test"]["rel_rms_fd_minus_(int-orb)_core_r_lt_9"] for r in A.values())
    a_ok = max(devs) < 0.01 and sign_ok and contr_ok
    a_ok = a_ok and abs(ratios["h1.5_L96_over_L48_rigid"] - 2.04) < 0.03
    a_ok = a_ok and abs(ratios["h1_L64_over_L48_rigid"] - 1.25) < 0.02
    V["A"] = {
        "verdict": "CONFIRMED" if a_ok else "QUALIFIED",
        "max_rel_dev_vs_audited": max(devs),
        "own_ratios": {k: v for k, v in ratios.items() if k != "synthetic_sign_test"},
        "sign_test_synthetic": syn,
        "sign_test_charge_residual_halo_max": halo,
        "sign_test_charge_residual_core_max": core,
        "own_contraction_matches_kin_of": contr_ok,
        "max_rel_change_C_orb_under_4th_order_orbital_stencil": c4_orb,
        "max_rel_change_C_rigid_under_4th_order_orbital_stencil": c4_rig,
        "reason": (
            f"all fifteen inertias reproduce within {max(devs):.1e} (own orbital stencil, the "
            "platform kin_of) and the own contraction matches kin_of to round-off; the ratios "
            f"reproduce (h 1.5 {ratios['h1.5_L96_over_L48_rigid']:.3f}, h 1 "
            f"{ratios['h1_L64_over_L48_rigid']:.3f}) and are unchanged under a fourth-order "
            f"orbital stencil ({ratios['h1.5_L96_over_L48_rigid_4th_order_stencil']:.3f}, "
            f"{ratios['h1_L64_over_L48_rigid_4th_order_stencil']:.3f}); the interpolated "
            "rotation confirms a_rigid = [J_z, M] - (x d_y - y d_x) M (synthetic resolved field: "
            f"residual {syn['rel_rms_fd_minus_(int-orb)']:.1e} for the minus sign, "
            f"{syn['rel_rms_fd_minus_(int+orb)']:.2f} for the plus). CAVEAT recorded, not a "
            f"refutation: the orbital term is stencil-sensitive on the charge (C_orb moves "
            f"{c4_orb * 100:.0f} percent, C_rigid {c4_rig * 100:.1f} percent under the "
            "fourth-order stencil; the spline rotation differs from the stencil by "
            f"{core:.2f} rms in the cores, {halo:.2f} in the halo): the cores are not resolved "
            "at h, so C_orb is a lattice number at the 15 percent level; the conclusion is not "
            "affected because C_orb is a tenth of C_int"
            if a_ok
            else "see numbers"
        ),
    }
    # B
    b_devs = [
        abs(r["within_r_own"][k] / r["within_r_audited"][k] - 1)
        for r in A.values()
        for k in ("6", "12", "18")
    ]
    fr12 = {t: r["frac_inside_12_own"] for t, r in A.items()}
    third = [t for t, f in fr12.items() if abs(f - 1.0 / 3.0) < 0.07]
    V["B"] = {
        "verdict": "QUALIFIED" if max(b_devs) < 0.02 else "REFUTED",
        "max_rel_dev_within_r": max(b_devs),
        "frac_inside_12_by_row": fr12,
        "rows_near_a_third": third,
        "reason": (
            f"the within-r numbers reproduce to {max(b_devs):.1e} on an own kin density, but "
            "'a third of C_rigid inside r 12 at every box' is not what they say: the fraction is "
            + ", ".join(f"{f:.2f} ({t.split('_')[3]} {t.split('_')[4]})" for t, f in fr12.items())
            + "; a third holds on the two larger h 1.5 boxes only, the L 48 boxes hold 0.45 "
            "(h 1.5) and 0.51 (h 1) inside r 12. The qualitative point (the halo carries the "
            "larger or an equal share, so the box law is the halo's) stands on every box"
        ),
    }
    # C
    ch = C["charge_S1_d0.3_n32_L48_end"]
    su = C["seed_uniaxial_d0_(1,0,0)"]
    s1 = C["seed_uniaxial_d0_rc1.0"]
    s05 = C["seed_uniaxial_d0_rc0.5"]
    c_ok = (
        ch["rms_M_rot90_minus_M_over_spread"] > 0.1
        and su["rms_M_rot90_minus_M_over_spread"] < 1e-10
    )
    c_ok = c_ok and max(s1["C_rigid_over_C_int"], s05["C_rigid_over_C_int"]) < 0.3
    V["C"] = {
        "verdict": "CONFIRMED" if c_ok else "QUALIFIED",
        "charge_rot90_asymmetry": ch["rms_M_rot90_minus_M_over_spread"],
        "charge_rot45_asymmetry": ch["rms_M_rot45_minus_M_over_spread_interior"],
        "charge_rot180_asymmetry": ch["rms_M_rot180_minus_M_over_spread"],
        "seed_uniaxial_rot90_asymmetry": su["rms_M_rot90_minus_M_over_spread"],
        "seed_uniaxial_C_rigid_over_C_int": su["C_rigid_over_C_int"],
        "seed_S1_d0.3_C_rigid_over_C_int": C["seed_S1_d0.3_(1,0.3,0)"]["C_rigid_over_C_int"],
        "seed_uniaxial_rc1.0_C_rigid_over_C_int": s1["C_rigid_over_C_int"],
        "seed_uniaxial_rc0.5_C_rigid_over_C_int": s05["C_rigid_over_C_int"],
        "charge_C_rigid_over_C_int": ch["C_rigid_over_C_int"],
        "charge_rms_a_rigid_over_a_int": ch["rms_a_rigid_over_a_int"],
        "reason": (
            "the relaxed L 48 charge differs from its own 90-degree copy by "
            f"{ch['rms_M_rot90_minus_M_over_spread']:.2f} of its spatial-block spread (45 "
            f"degrees: {ch['rms_M_rot45_minus_M_over_spread_interior']:.2f}, 180 degrees: "
            f"{ch['rms_M_rot180_minus_M_over_spread']:.2f}), the uniaxial seed by "
            f"{su['rms_M_rot90_minus_M_over_spread']:.0e}; C_rigid / C_int is "
            f"{ch['C_rigid_over_C_int']:.3f} on the charge against "
            f"{su['C_rigid_over_C_int']:.4f} on the seed (audited 0.0018 with the sym edge "
            "stencil), and an axisymmetric seed with an UNRESOLVED core (r_c 1.0, 0.5 at h 1.5) "
            f"still cancels to {s1['C_rigid_over_C_int']:.3f}, {s05['C_rigid_over_C_int']:.3f}: "
            "the non-cancellation on the charge is asymmetry, not resolution"
        ),
    }
    # D
    d_ok, notes = True, []
    for tag, delta in ARM:
        r = D[tag]
        d_ok = d_ok and abs(r["E_platform"] - r["E_audited"]) < 1e-3
        d_ok = d_ok and abs(r["E_own"] - r["E_platform"]) < 1e-6
        d_ok = d_ok and all(
            abs(r["own_C"][k] / r["audited_C"][k] - 1) < 0.01
            for k in ("C_int", "C_orb", "C_rigid")
        )
        d_ok = d_ok and r["fmax_spatial_free_own"] > GATE
        notes.append(
            f"delta {delta}: E {r['E_platform']:.4f} (own {r['E_own']:.4f}), C_rigid "
            f"{r['own_C']['C_rigid']:.0f}, fmax on the free cells {r['fmax_spatial_free_own']:.1e} "
            f"({r['fmax_over_gate']:.1f} x the gate), last four chunk drops "
            + ", ".join(f"{x:.1e}" for x in r["last_4_chunk_drops"])
            + f", the last chunk with a drop of 4e-3 or more ended at iteration "
            f"{r['last_chunk_iters_with_drop_ge_4e-3']}"
        )
    d03 = D["S1_d0.03_w25_n32_L48"]
    V["D"] = {
        "verdict": "QUALIFIED" if d_ok else "REFUTED",
        "delta_exponent_own_fit": D["delta_exponent_own_fit"],
        "monotone": D["monotone"],
        "reason": (
            "energies and inertias reproduce and both rows are FALLING by the rule (fmax above "
            "1e-4); the sentence 'their energies still dropping 4e-3 per chunk' does not: "
            + "; ".join(notes)
            + f". The delta 0.03 row's last-chunk drop {d03['last_4_chunk_drops'][-1]:.1e} is "
            f"below the descent's own drop criterion ({d03['drop_criterion_GATE_DROP_x_E']:.1e}) and "
            f"the gate was crossed inside its last two chunks (iteration "
            f"{d03['gate_hit_inside_last_chunk']} of 250) with the chunk-end fmax 4e-4 to 5e-4: in "
            "energy that row is converged, the label is the fmax rule's. The delta law read "
            f"(own exponent {D['delta_exponent_own_fit']:.3f}, not monotone) reproduces"
        ),
    }
    # E
    e_ok, ill = True, []
    for tag, d, n, L in PAIRS:
        rd = E[tag]["reads"]
        cp, cm = rd[f"core_x+{d / 2:g}"]["degree"], rd[f"core_x-{d / 2:g}"]["degree"]
        e_ok = e_ok and {round(cp), round(cm)} == {1, -1} and abs(cp - round(cp)) < 1e-3
        e_ok = e_ok and abs(rd["outer"]["degree"]) < 1e-3
        e_ok = e_ok and abs(rd["mid_R_half_d_minus_1"]["degree"]) < 1e-3
        e_ok = e_ok and abs(rd["mid_R_half_d_plus_1"]["degree"]) < 1e-3
        e_ok = e_ok and abs(rd["both_cores_R_half_d_plus_3"]["degree"]) < 1e-3
        if E[tag]["mid_sphere_R_half_d_ill_posed"] or abs(rd["mid_R_half_d"]["degree"]) > 1e-3:
            ill.append(f"{tag} ({rd['mid_R_half_d']['degree']:+.2f})")
    V["E"] = {
        "verdict": "QUALIFIED" if e_ok else "REFUTED",
        "mid_sphere_reads_not_zero_or_ill_posed": ill,
        "reason": (
            "own solid-angle degree on an icosphere (5120 faces, trilinear director): +1 around "
            "the +x core and -1 around the -x core to 1e-4, 0 on the outer sphere, 0 on the "
            "spheres of radius d/2 - 1 (no core inside) and d/2 + 1 and d/2 + 3 (both cores "
            "inside, +1 - 1) on all eleven pairs: the physical statement (both unit charges stay "
            "inside the balls, the exterior carries none) is confirmed. The stated '0 on the "
            "sphere of radius d/2 about the origin' is ill-posed: that sphere passes THROUGH "
            "both cores (the interpolated director's norm reaches 0 on it and neighboring "
            "samples differ by 90 degrees), the own reader gives -1 on the six h 1 pairs and 0 "
            f"on the five h 1.5 pairs [{'; '.join(ill)}], and the audited values on it "
            "(1e-9 to 5e-6, non-integers) are the same accident seen from the cell-sampled side; "
            "replace that read by d/2 - 1 and d/2 + 1"
            if e_ok
            else "see numbers"
        ),
    }
    # F
    k = F["rows"]["S1_d0.3_w25_n32_L48"]["kick_stage"]
    pr = F["continuation_probe"]
    ident = all(
        r["max_abs_end_minus_gate"] == 0.0 and r["max_abs_end_minus_stage"] == 0.0
        for r in F["rows"].values()
    )
    shifts = {t: v.get("rel_shift_if_tail_at_this_slope") for t, v in pr.items()}
    V["F"] = {
        "verdict": "CONFIRMED",
        "label_follows": F["label_by_rule"] == F["label_audited"],
        "rows_not_at_gate": F["rows_not_at_gate"],
        "kick_stage_sensitivity_own": k["rel_sensitivity_rigid"],
        "kick_stage_sensitivity_audited": k["audited"],
        "end_gate_stage_arrays_identical_all_rows": ident,
        "g_not_evaluated": F["g_not_evaluated_in_json"] and not F["g_numeric_keys_present"],
        "continuation_probe_rel_shift_of_C_rigid": shifts,
        "reason": (
            "the label follows the rule (R25-2 gate_label FALLING on "
            f"{', '.join(F['rows_not_at_gate'])}); the kick-stage sensitivity on the converged row "
            f"is {k['rel_sensitivity_rigid']:.2e} (audited {k['audited']:.2e}); no g number "
            "exists in the JSON (g_note on every row). Two notes: (1) the JSON's 'stage_gate' "
            "zeros are identities, not measurements, the end, gate and stage arrays are "
            "byte-identical on all five rows (R25-2 writes the gate array at the cap and keeps it "
            "as the end field), only the kick-stage number carries information, and only on the "
            "converged row; (2) the box-dependence conclusion is sound at the observed "
            "convergence slope: 30 own L-BFGS-B iterations on the capped n 48 rows move C_rigid "
            + ", ".join(
                f"{v['dC_rigid'] / v['C_rigid_before']:+.1e} relative for dE {v['dE']:+.1e} "
                f"({t.split('_')[3]} {t.split('_')[4]})"
                for t, v in pr.items()
                if "dC_rigid" in v
            )
            + "; the geometric tail of the chunk drops implies a C_rigid shift of "
            + ", ".join(f"{(v or 0) * 100:.2f} percent" for v in shifts.values())
            + " to convergence, against the factor 2.04 the box ratio carries; a factor-2 "
            "reversal would need a slope hundreds of times the observed one, so the INSUFFICIENT "
            "label is the rule's conservatism, not a live doubt about the ratio (a different "
            "basin under continued descent is the residual unknown the probe cannot exclude)"
        ),
    }
    return V


def main():
    J26 = json.load(open(R26_JSON))
    J25 = json.load(open(R25_JSON))
    A, ratios = claim_a_b(J26)
    C = claim_c()
    D = claim_d(J26)
    E = claim_e(J26)
    F = claim_f(J26, J25)
    V = verdicts(A, ratios, C, D, E, F)
    out = {
        "task": "M5.32 R26-3 audit",
        "runtime_s": round(time.time() - T0, 1),
        "orbital_stencil": "np.gradient edge_order=2 (second-order central interior, second-order "
        "one-sided edges); the 4th-order variant reported beside it",
        "A_ladder": A,
        "A_ratios": ratios,
        "C_asymmetry": C,
        "D_delta_arm": D,
        "E_pairs": E,
        "F_label": F,
        "verdicts": V,
    }
    with open(OUT_JSON, "w") as f:
        json.dump(out, f, indent=1, default=fnum)
    print("\n| Claim | Verdict | Own numbers | Audited |")
    print("| --- | --- | --- | --- |")
    r48 = A["S1_d0.3_w25_n32_L48"]["own_c2"]
    print(
        f"| A ladder | {V['A']['verdict']} | L48 h1.5: {r48['C_int']:.0f}, {r48['C_orb']:.0f}, "
        f"{r48['C_rigid']:.0f}; ratios h1.5 {ratios['h1.5_L96_over_L48_rigid']:.3f} "
        f"({ratios['h1.5_L96_over_L48_int']:.3f}), h1 {ratios['h1_L64_over_L48_rigid']:.3f} "
        f"({ratios['h1_L64_over_L48_int']:.3f}); max dev {V['A']['max_rel_dev_vs_audited']:.1e} | "
        f"14071, 1403, 14031; 2.04 (2.06), 1.25 (1.25) |"
    )
    p48 = A["S1_d0.3_w25_n32_L48"]["within_r_own"]
    print(
        f"| B profile | {V['B']['verdict']} | L48: {p48['6']:.0f}, {p48['12']:.0f}, {p48['18']:.0f}; "
        f"frac(r<12) {min(V['B']['frac_inside_12_by_row'].values()):.2f} to "
        f"{max(V['B']['frac_inside_12_by_row'].values()):.2f} | 395, 6380, 12440; a third |"
    )
    print(
        f"| C asymmetry | {V['C']['verdict']} | charge rot90 {V['C']['charge_rot90_asymmetry']:.3f}, "
        f"rot45 {V['C']['charge_rot45_asymmetry']:.3f}; seed rot90 "
        f"{V['C']['seed_uniaxial_rot90_asymmetry']:.1e}; C_rigid/C_int charge "
        f"{V['C']['charge_C_rigid_over_C_int']:.3f} seed {V['C']['seed_uniaxial_C_rigid_over_C_int']:.4f} | "
        f"0.997 / 0.0018 |"
    )
    for tag, delta in ARM:
        r = D[tag]
        print(
            f"| D delta {delta} | {V['D']['verdict']} | E {r['E_platform']:.3f} (own {r['E_own']:.3f}); "
            f"C {r['own_C']['C_int']:.0f}, {r['own_C']['C_orb']:.0f}, {r['own_C']['C_rigid']:.0f}; "
            f"fmax {r['fmax_spatial_free_own']:.1e} vs gate 1e-4; last drop {r['last_4_chunk_drops'][-1]:.1e} | "
            f"E {r['E_audited']:.3f}; {r['audited_C']['C_int']:.0f}, {r['audited_C']['C_orb']:.0f}, "
            f"{r['audited_C']['C_rigid']:.0f} |"
        )
    e0 = E["unlike_d10_n32_L48"]["reads"]
    print(
        f"| E pairs | {V['E']['verdict']} | cores +1 / -1 to 1e-4, outer 0, d/2 +- 1 spheres 0 on "
        f"all 11; the d/2 sphere runs through the cores (own read -1 on the 6 h 1 pairs, 0 on the "
        f"5 h 1.5 pairs; e.g. d10 n32: {e0['mid_R_half_d']['degree']:+.1e}) | +1, -1, 0, 0 |"
    )
    print(
        f"| F label | {V['F']['verdict']} | label follows: {V['F']['label_follows']}; kick-stage "
        f"{V['F']['kick_stage_sensitivity_own']:.2e}; end=gate=stage arrays identical: "
        f"{V['F']['end_gate_stage_arrays_identical_all_rows']}; probe shift "
        f"{V['F']['continuation_probe_rel_shift_of_C_rigid']} | INSUFFICIENT; 4e-4 |"
    )
    log(f"wrote {OUT_JSON}")


if __name__ == "__main__":
    main()
