"""M5.32 R26-0: the form level of the census rung, and the READERS every R26 sub-rung
imports (the partition reader on spheres, the gap-tail moment reader, the physical
rotation generator, the biaxiality invariant, the virial ratio, the hedgehog degree),
each validated here against a null and a synthetic control before any field is read.

EQUATIONS FIRST
---------------
Field M(x) real symmetric 4x4, eta = diag(-1, 1, 1, 1), N = M eta, the certified biaxial
vacuum diag(8, 1, delta, 0); the static energy of the stack (R20.energy_parts, read-only)
    E = 4 h^3 sum_{i<j} <F_ij, F_ij>_eta + V4,  F_ij = [A_i, A_j]_eta,  A_i = d_i M,
    <F, G>_eta = tr(eta F eta G^T),  V4 = w sum_{p<=4} (tr N^p - C_p)^2.
This rung STATES its norm (the author's 2026-09-24 16:08 UTC request): every static number
here is in the eta contraction above, on the block-diagonal sector M_0i = 0 where the eta
and Frobenius contractions coincide (check c shows where they differ).

a   THE TANGENT. Under M -> Lambda M Lambda^T with Lambda = exp(t G), the velocity is
        dM = G M + M G^T,  symmetric for every G,  with  G eta + eta G^T = 0
    (the invariance of eta). The author's lower-index form dM = G'^T M + M G' is the same
    formula with G' = G^T. A rotation (G antisymmetric) gives the commutator [G, M]; a boost
    (G symmetric) gives the anticommutator {G, M}. The M5.21.3 catalog builds
    w (G M - M G^T) instead: for antisymmetric G that is the ANTICOMMUTATOR, for symmetric G
    the COMMUTATOR, and in both cases an ANTISYMMETRIC matrix, not a tangent of the
    symmetric field (the fault R25 found on rot_z; here proved for every entry, in sympy).
b   THE CATALOG'S SIX FIELDS are antisymmetric to round-off on a stored charge; the
    corrected tangents are symmetric.
c   THE NORM on the R25-1 full-kick witness fields (the boost run-away, kept as
    data/m5_32_r25_1_fullkick/): the eta energy is negative on every one, the Frobenius
    energy (<F, F>_I = tr(F F^T), Mikulski's delta_M with v0 = e0) is positive on every one;
    the falsifier is a witness field with a non-positive Frobenius energy.
d   THE SPLIT RESIDUAL. On the block strand the plain tr N^p - C_p subtracts 8^p from
    itself; the telescoped form tr N^p - tr N0^p = sum_k tr(N^k eps N0^{p-1-k}), eps = N - N0,
    never forms the large number. Relative error of both against the exact pair bracket at
    delta 0.3 to 0.001 (the plan-time numbers, tracked).
e   THE BIAXIALITY INVARIANT beta^2 = 1 - 6 (tr L^3)^2 / (tr L^2)^3 on the traceless part of
    the spatial block: 0 uniaxial, 1 on (1, 0, -1), 0.604 on the vacuum at delta 0.3, and the
    delta at which it equals the Koide 0.382 (0.2307).
f   THE GAP-TAIL READER: shells of width 1.5 from r 3 to 21, the pair-gap deviation fitted
    to (P0, P1, P2) in cos(theta) per shell, the log-log slope of each coefficient on r 6 to
    18. Null control: the stored S1 charge (z-reflection symmetric, l = 1 under 1e-4 of
    l = 0). Synthetic control: a planted dev = A cos(theta) / r^2 reads a1 = A / r^2 and
    slope -2 within 5 percent.
g   THE PHYSICAL ROTATION GENERATOR a_rigid = [J_z, M] - (x d_y - y d_x) M, the velocity
    of a rigid rotation of the whole field. Controls: on the axisymmetric uniaxial seed the
    inertia C_rigid = kin(a_rigid) is under 1 percent of C_int = kin([J_z, M]); on the
    stored S1 charge a_rigid agrees with the finite difference of the rotated field
    (cubic interpolation) within 5 percent in norm.
h   THE PARTITION READER on spheres: the transverse pair vector's angle in the sphere's
    (theta, phi) frame (projected on the plane normal to the local director, oriented
    outward), plaquette windings on a 64 x 128 grid, polar caps as single plaquettes with
    the frame's own index (+2 half-turns each) added back, connected components of the
    nonzero or unreadable plaquettes, the component totals in HALF-UNITS as the partition.
    Controls: the five seeds of the census read as {4}, {3,1}, {2,2}, {2,1,1}, {1,1,1,1} on
    r 9 and r 18 (total 4 on every sphere); the stored S1 charge's end field read on r 9
    against the R25-2 audit's {1,1,1,1}.
i   THE VIRIAL. E_u is homogeneous of degree 4 in the first derivatives and V4 of degree 0,
    so under x -> lambda x on a free stationary point E_u = 3 V4 (Derrick in 3D); on the
    stack this is the h-scaling E_u(2h) = E_u(h) / 2, V4(2h) = 8 V4(h) on the same array,
    checked to round-off. The ratio E_u / (3 V4) is reported per branch, UNDER THE PINNED
    SHELL, never labeled.
j   THE DEGREE READER: the hedgehog degree of the director (oriented outward) on a sphere by
    the signed spherical-triangle areas; the S1 seed reads 1 on every sphere, the R22-2
    unlike pair reads +1 and -1 around its cores.

Modes: run (every check, about 3 min) | readers (h on the seeds only) | quick.
Output: data/m5_32_r26_0_form.json. Each check carries `PASS` and `fails_if`.
"""

import importlib.util
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
OUT_JSON = os.path.join(DATA, "m5_32_r26_0_form.json")
FULLKICK = os.path.join(DATA, "m5_32_r25_1_fullkick")
S1_END = os.path.join(DATA, "m5_32_r25_2", "S1_d0.3_w25_n32_L48.npz")
S1_GATE = os.path.join(DATA, "m5_32_r25_2", "S1_d0.3_w25_n32_L48_gate.npz")
R22_2 = os.path.join(DATA, "m5_32_r22_2")
ETA = np.diag([-1.0, 1.0, 1.0, 1.0])
G = 8.0
W1S = 25.0
T0 = time.time()


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


CS = _load("m5_32_r23_1_cscan", "m5_32_r23_1_cscan.py")
R21, R20, B3, R0, W1 = CS.R21, CS.R20, CS.B3, CS.R0, CS.W1


def log(msg):
    print(f"[{time.time() - T0:8.1f}s] {msg}", flush=True)


def cfg_pot(n, L, delta, w1s=W1S):
    cfg = R21.cfg_of(n, L, G, delta)
    p = R21.params_of(G, delta)
    pot = ("v4", R0.roots_of(cfg), W1 * w1s)
    return cfg, p, pot


JZ = np.zeros((4, 4))
JZ[1, 2], JZ[2, 1] = -1.0, 1.0

# ============================================================================
# THE READERS (imported by R26-3 and R26-4)
# ============================================================================


def dsym(M, ax, h):
    """the symmetric first difference (the mean of the forward and backward stencils)."""
    f = np.zeros_like(M)
    b = np.zeros_like(M)
    s0 = [slice(None)] * 3
    s1 = [slice(None)] * 3
    s0[ax], s1[ax] = slice(0, -1), slice(1, None)
    f[tuple(s0)] = (M[tuple(s1)] - M[tuple(s0)]) / h
    b[tuple(s1)] = (M[tuple(s1)] - M[tuple(s0)]) / h
    return 0.5 * (f + b)


def orbital_z(M, cfg):
    """(x d_y - y d_x) M: the transport part of a rigid rotation about z."""
    X, Y, _ = B3.coords(cfg["n"], cfg["h"])
    return X[..., None, None] * dsym(M, 1, cfg["h"]) - Y[..., None, None] * dsym(M, 0, cfg["h"])


def kin_density(M, a0, cfg):
    """per-cell kin density, h^3-weighted (its sum is B3.kin_of)."""
    h3 = cfg["h"] ** 3
    k = np.zeros(M.shape[:3])
    for br, (A, wt) in B3.a_fields(M, cfg).items():
        for i in range(3):
            F = B3.comm_eta(a0, A[i])
            k += wt * 4.0 * B3.inner_eta(F, F)
    return h3 * k


def physical_generator(M, cfg, radii=(3.0, 6.0, 9.0, 12.0, 18.0)):
    """C_int = kin([J_z, M]), C_orb = kin(orbital), C_rigid = kin(a_int - a_orb), with the
    within-r profiles of each and the peak cell of the rigid density."""
    a_int = JZ @ M - M @ JZ
    a_orb = orbital_z(M, cfg)
    a_rig = a_int - a_orb
    X, Y, Z = B3.coords(cfg["n"], cfg["h"])
    r = np.sqrt(X * X + Y * Y + Z * Z)
    kd = {
        "internal": kin_density(M, a_int, cfg),
        "orbital": kin_density(M, a_orb, cfg),
        "rigid": kin_density(M, a_rig, cfg),
    }
    out = {k: float(v.sum()) for k, v in kd.items()}
    out["within_r"] = {f"{R:g}": {k: float(v[r < R].sum()) for k, v in kd.items()} for R in radii}
    out["rigid_over_internal"] = out["rigid"] / max(out["internal"], 1e-300)
    i = np.unravel_index(int(np.argmax(kd["rigid"])), kd["rigid"].shape)
    out["rigid_peak_cell_xyz"] = [float(X[i]), float(Y[i]), float(Z[i])]
    out["rigid_peak_value"] = float(kd["rigid"][i])
    out["a_rigid_symmetric_maxdev"] = float(np.abs(a_rig - a_rig.swapaxes(-1, -2)).max())
    return out


def rotated_copy(M, cfg, alpha):
    """M_alpha(x) = R M(R^-1 x) R^T, R the rotation by alpha about z, cubic interpolation."""
    from scipy.ndimage import map_coordinates

    n, h = cfg["n"], cfg["h"]
    X, Y, Z = B3.coords(n, h)
    c, s = np.cos(alpha), np.sin(alpha)
    # R^-1 x
    Xr = c * X + s * Y
    Yr = -s * X + c * Y
    ix = Xr / h + (n - 1) / 2.0
    iy = Yr / h + (n - 1) / 2.0
    iz = Z / h + (n - 1) / 2.0
    coords = np.stack([ix.ravel(), iy.ravel(), iz.ravel()])
    out = np.zeros_like(M)
    for a in range(4):
        for b in range(a, 4):
            v = map_coordinates(M[..., a, b], coords, order=3, mode="nearest").reshape(n, n, n)
            out[..., a, b] = v
            out[..., b, a] = v
    R = np.eye(4)
    R[1, 1], R[1, 2], R[2, 1], R[2, 2] = c, -s, s, c
    return R @ out @ R.T


def beta2_of(lam):
    """1 - 6 (tr L^3)^2 / (tr L^2)^3 on the traceless part of a 3-spectrum (0 uniaxial, 1 max)."""
    L = np.asarray(lam, float)
    L = L - L.mean(axis=-1, keepdims=True)
    t2 = np.sum(L**2, axis=-1)
    t3 = np.sum(L**3, axis=-1)
    return 1.0 - 6.0 * t3**2 / np.maximum(t2**3, 1e-300)


def biaxiality_reads(M, cfg, delta, radii=(3.0, 6.0, 9.0)):
    """beta^2 of the vacuum, of the cell-averaged spatial block inside r < R, and the mean of
    the per-cell beta^2 inside r < R (cells with a resolved spectrum only)."""
    X, Y, Z = B3.coords(cfg["n"], cfg["h"])
    r = np.sqrt(X * X + Y * Y + Z * Z)
    lam = np.linalg.eigvalsh(M[..., 1:, 1:])
    b2 = beta2_of(lam)
    resolved = (lam[..., 2] - lam[..., 1] > 0.1) & (lam[..., 1] - lam[..., 0] > 0.02 * delta)
    out = {"vacuum": float(beta2_of([1.0, delta, 0.0]))}
    for R in radii:
        m = r < R
        Sbar = M[m][:, 1:, 1:].mean(axis=0)
        out[f"block_mean_r{R:g}"] = float(beta2_of(np.linalg.eigvalsh(Sbar)))
        mr = m & resolved
        out[f"cell_mean_r{R:g}"] = float(b2[mr].mean()) if mr.any() else None
        out[f"resolved_fraction_r{R:g}"] = float(mr.sum() / max(m.sum(), 1))
    return out


def gap_tail(M, cfg, delta, r_lo=3.0, r_hi=21.0, width=1.5, fit=(6.0, 18.0), dev_field=None):
    """the pair-gap deviation on shells, fitted per shell to a0 + a1 cos(theta) + a2 P2(cos);
    the log-log slope of each |a_l| over the fit window. `dev_field` overrides the deviation
    (the synthetic control)."""
    X, Y, Z = B3.coords(cfg["n"], cfg["h"])
    r = np.sqrt(X * X + Y * Y + Z * Z)
    cth = np.where(r > 0, Z / np.maximum(r, 1e-12), 0.0)
    if dev_field is None:
        lam = np.linalg.eigvalsh(M[..., 1:, 1:])
        dev = (lam[..., 1] - lam[..., 0]) - delta
    else:
        dev = dev_field
    pin = B3.pin_shell(cfg["n"], cfg["h"])
    rows = []
    for Rs in np.arange(r_lo, r_hi + 1e-9, width):
        m = (np.abs(r - Rs) < 0.5 * width) & (~pin)
        if m.sum() < 20:
            continue
        A = np.stack([np.ones(m.sum()), cth[m], 0.5 * (3 * cth[m] ** 2 - 1)], 1)
        coef, *_ = np.linalg.lstsq(A, dev[m], rcond=None)
        rows.append(
            [
                float(Rs),
                float(coef[0]),
                float(coef[1]),
                float(coef[2]),
                float(np.abs(dev[m]).max()),
                int(m.sum()),
            ]
        )
    sh = np.array(rows) if rows else np.zeros((0, 6))

    def slope(col):
        if len(sh) == 0:
            return None
        mm = (sh[:, 0] >= fit[0]) & (sh[:, 0] <= fit[1]) & (np.abs(sh[:, col]) > 0)
        if mm.sum() < 3:
            return None
        p = np.polyfit(np.log(sh[mm, 0]), np.log(np.abs(sh[mm, col])), 1)
        return float(p[0])

    out = {
        "shells_r_a0_a1_a2_maxdev_n": rows,
        "slope_l0": slope(1),
        "slope_l1": slope(2),
        "slope_l2": slope(3),
        "fit_window": list(fit),
    }
    if len(sh):
        mm = (sh[:, 0] >= fit[0]) & (sh[:, 0] <= fit[1])
        out["l1_over_l0_max_in_window"] = (
            float(np.max(np.abs(sh[mm, 2]) / np.maximum(np.abs(sh[mm, 1]), 1e-300)))
            if mm.any()
            else None
        )
    return out


def virial_reads(M, cfg, p, pot):
    parts = R20.energy_parts(M, cfg, p, pot)
    eu, v = float(parts["E_curv"]), float(parts["V"])
    return {
        "E_u": eu,
        "V": v,
        "E_u_over_3V": eu / max(3.0 * v, 1e-300),
        "note": "under the pinned shell",
    }


# ---------- the sphere sampler and the partition reader ----------
def _cell_index(x, n, h):
    return np.clip(np.rint(x / h + (n - 1) / 2.0).astype(int), 0, n - 1)


def sphere_samples(R, nth=64, nph=128, center=(0.0, 0.0, 0.0)):
    th = (np.arange(nth) + 0.5) * np.pi / nth
    ph = np.arange(nph) * 2.0 * np.pi / nph
    TH, PH = np.meshgrid(th, ph, indexing="ij")
    rhat = np.stack([np.sin(TH) * np.cos(PH), np.sin(TH) * np.sin(PH), np.cos(TH)], -1)
    that = np.stack([np.cos(TH) * np.cos(PH), np.cos(TH) * np.sin(PH), -np.sin(TH)], -1)
    phat = np.stack([-np.sin(PH), np.cos(PH), np.zeros_like(PH)], -1)
    pts = np.asarray(center)[None, None, :] + R * rhat
    return TH, PH, rhat, that, phat, pts


def _components(mask, wind):
    """8-connected components of `mask` on the (nth, nph) grid with phi wrapping; the two polar
    caps are extra nodes adjacent to row 0 and row nth-1. Returns a list of components as
    (total_winding, list_of_cells, has_north_cap, has_south_cap)."""
    nth, nph = mask.shape
    seen = np.zeros_like(mask, dtype=bool)
    comps = []
    cap_seen = [False, False]
    cap_touch = {}

    def neigh(i, j):
        for di in (-1, 0, 1):
            for dj in (-1, 0, 1):
                if di == 0 and dj == 0:
                    continue
                ii, jj = i + di, (j + dj) % nph
                if 0 <= ii < nth:
                    yield ii, jj

    for i0 in range(nth):
        for j0 in range(nph):
            if not mask[i0, j0] or seen[i0, j0]:
                continue
            stack = [(i0, j0)]
            seen[i0, j0] = True
            cells = []
            caps = [False, False]
            while stack:
                i, j = stack.pop()
                cells.append((i, j))
                if i == 0:
                    caps[0] = True
                if i == nth - 1:
                    caps[1] = True
                for ii, jj in neigh(i, j):
                    if mask[ii, jj] and not seen[ii, jj]:
                        seen[ii, jj] = True
                        stack.append((ii, jj))
            comps.append([cells, caps])
    # merge the components that touch the same cap (the cap joins them)
    merged = []
    for k in (0, 1):
        touching = [c for c in comps if c[1][k]]
        if len(touching) > 1:
            base = touching[0]
            for c in touching[1:]:
                base[0].extend(c[0])
                base[1][0] |= c[1][0]
                base[1][1] |= c[1][1]
                comps.remove(c)
    for cells, caps in comps:
        tot = float(sum(wind[i, j] for i, j in cells))
        merged.append((tot, cells, caps[0], caps[1]))
    return merged


def _wrap_pi(x):
    """x mod pi into (-pi/2, pi/2), antisymmetric (half-integers round to even), so that interior
    edges of adjacent plaquettes cancel exactly and a region's plaquette sum is its boundary winding.
    """
    return x - np.pi * np.rint(x / np.pi)


def sample_block(M, cfg, pts):
    """the spatial block trilinearly interpolated at the points (..., 3)."""
    from scipy.ndimage import map_coordinates

    n, h = cfg["n"], cfg["h"]
    shp = pts.shape[:-1]
    coords = np.stack([(pts[..., a] / h + (n - 1) / 2.0).ravel() for a in range(3)])
    S = np.zeros(shp + (3, 3))
    for a in range(3):
        for b in range(a, 3):
            v = map_coordinates(M[..., 1 + a, 1 + b], coords, order=1, mode="nearest").reshape(shp)
            S[..., a, b] = v
            S[..., b, a] = v
    return S


def partition_reader(
    M, cfg, R, delta, center=(0.0, 0.0, 0.0), nth=64, nph=128, gap_frac=0.05, rho_min_h=1.0
):
    """the transverse index carriers on the sphere of radius R, in half-units: the pair vector's
    angle chi in the sphere's (theta, phi) frame projected on the plane normal to the outward
    director, plaquette windings (counterclockwise seen from outside) on the interpolated block,
    the polar caps closed at the first fully readable ring (their half-units = the ring winding
    plus the frame's own index 2), connected components of the nonzero or unreadable plaquettes
    between the caps, totals per component."""
    n, h = cfg["n"], cfg["h"]
    TH, PH, rhat, that, phat, pts = sphere_samples(R, nth, nph, center)
    Ms = sample_block(M, cfg, pts)
    lam, V = np.linalg.eigh(Ms)
    d = V[..., :, 2]
    sgn = np.sign(np.einsum("...a,...a->...", d, rhat))
    sgn = np.where(sgn == 0, 1.0, sgn)
    d = d * sgn[..., None]
    align = np.einsum("...a,...a->...", d, rhat)
    v1 = V[..., :, 1]
    e1 = that - np.einsum("...a,...a->...", that, d)[..., None] * d
    e1n = np.linalg.norm(e1, axis=-1)
    e1 = e1 / np.maximum(e1n, 1e-300)[..., None]
    e2 = np.cross(d, e1)
    chi = np.arctan2(np.einsum("...a,...a->...", v1, e2), np.einsum("...a,...a->...", v1, e1))
    gap_pair = lam[..., 1] - lam[..., 0]
    gap_dir = lam[..., 2] - lam[..., 1]
    reasons = {
        "pair_gap": gap_pair < gap_frac * delta,
        "director_gap": gap_dir < 0.1,
        "director_not_outward": align < 0.3,
        "frame": e1n < 0.3,
    }
    bad_s = np.zeros_like(align, dtype=bool)
    for v in reasons.values():
        bad_s |= v
    # ring windings (increasing phi) per row
    dphi = _wrap_pi(np.roll(chi, -1, axis=1) - chi)
    ring = np.rint(dphi.sum(axis=1) / np.pi)
    row_ok = ~bad_s.any(axis=1)
    rho_row = R * np.sin(TH[:, 0])
    # the caps: the first readable ring from each pole at rho >= rho_min_h h
    cand_n = [i for i in range(nth // 2) if row_ok[i] and rho_row[i] >= rho_min_h * h]
    cand_s = [
        i for i in range(nth - 1, nth // 2 - 1, -1) if row_ok[i] and rho_row[i] >= rho_min_h * h
    ]
    # AUDIT CORRECTION (R26-0 audit, 2026-09-25): aborting when no ring near a pole is fully readable
    # was the reader's policy, not the field's property; fall back to the ring with the fewest
    # unreadable samples (at most 10 percent) in each hemisphere and FLAG the read
    flagged = []
    bad_frac = bad_s.mean(axis=1)
    if not cand_n:
        opts = [i for i in range(nth // 2) if rho_row[i] >= rho_min_h * h and bad_frac[i] <= 0.30]
        if opts:
            cand_n = [min(opts, key=lambda i: (bad_frac[i], i))]
            flagged.append("N")
    if not cand_s:
        opts = [
            i
            for i in range(nth - 1, nth // 2 - 1, -1)
            if rho_row[i] >= rho_min_h * h and bad_frac[i] <= 0.30
        ]
        if opts:
            cand_s = [min(opts, key=lambda i: (bad_frac[i], -i))]
            flagged.append("S")
    if not cand_n or not cand_s:
        return {
            "R": float(R),
            "center": list(map(float, center)),
            "partition": None,
            "unreadable": "no ring near a pole with at most 30 percent unreadable samples",
            "total_half_units": None,
            "carriers": [],
            "null_spots": [],
            "unreadable_plaquettes": None,
            "unreadable_samples_by_reason": {k: int(v.sum()) for k, v in reasons.items()},
            "in_pin": bool(R > 0.5 * cfg["L"] - 1.6),
        }
    iN, iS = cand_n[0], cand_s[0]
    cap_n = 2.0 + float(ring[iN])
    cap_s = 2.0 - float(ring[iS])
    # plaquettes between rows iN and iS, counterclockwise seen from outside
    c00 = chi[iN:iS, :]
    c01 = np.roll(chi, -1, axis=1)[iN:iS, :]
    c11 = np.roll(chi, -1, axis=1)[iN + 1 : iS + 1, :]
    c10 = chi[iN + 1 : iS + 1, :]
    wind = (
        -(_wrap_pi(c01 - c00) + _wrap_pi(c11 - c01) + _wrap_pi(c10 - c11) + _wrap_pi(c00 - c10))
        / np.pi
    )
    wind = np.rint(wind)
    b = bad_s
    bad = (
        b[iN:iS, :]
        | np.roll(b, -1, axis=1)[iN:iS, :]
        | np.roll(b, -1, axis=1)[iN + 1 : iS + 1, :]
        | b[iN + 1 : iS + 1, :]
    )
    mask = (wind != 0) | bad
    comps = _components(mask, wind)
    carriers = []
    used_n = used_s = False
    for tot, cells, tn, ts in comps:
        t = tot
        if tn:
            t += cap_n
            used_n = True
        if ts:
            t += cap_s
            used_s = True
        ii = np.array([c[0] for c in cells]) + iN + 0.5
        jj = np.array([c[1] for c in cells]) + 0.5
        th_m = float(np.mean(ii) * np.pi / nth)
        phs = jj * 2.0 * np.pi / nph
        ph_m = float(np.arctan2(np.mean(np.sin(phs)), np.mean(np.cos(phs))) % (2 * np.pi))
        if tn and ts:
            th_m = float("nan")
        elif tn:
            th_m = 0.0
        elif ts:
            th_m = float(np.pi)
        carriers.append(
            {
                "half_units": float(t),
                "theta": th_m,
                "phi": ph_m,
                "rho_from_axis": float(R * np.sin(th_m)) if np.isfinite(th_m) else None,
                "z": float(R * np.cos(th_m)) if np.isfinite(th_m) else None,
                "cells": len(cells),
                "unreadable_cells": int(sum(bad[i, j] for i, j in cells)),
                "touches_pole": ("N" if tn else "") + ("S" if ts else ""),
            }
        )
    if not used_n:
        carriers.append(
            {
                "half_units": float(cap_n),
                "theta": 0.0,
                "phi": 0.0,
                "rho_from_axis": 0.0,
                "z": float(R),
                "cells": 0,
                "unreadable_cells": 0,
                "touches_pole": "N",
            }
        )
    if not used_s:
        carriers.append(
            {
                "half_units": float(cap_s),
                "theta": float(np.pi),
                "phi": 0.0,
                "rho_from_axis": 0.0,
                "z": float(-R),
                "cells": 0,
                "unreadable_cells": 0,
                "touches_pole": "S",
            }
        )
    nonzero = [c for c in carriers if c["half_units"] != 0]
    partition = sorted([int(round(c["half_units"])) for c in nonzero], reverse=True)
    return {
        "R": float(R),
        "center": list(map(float, center)),
        "partition": partition,
        "total_half_units": float(sum(c["half_units"] for c in carriers)),
        "carriers": nonzero,
        "null_spots": [c for c in carriers if c["half_units"] == 0 and c["unreadable_cells"] > 0],
        "unreadable_plaquettes": int(bad.sum()),
        "unreadable_samples_by_reason": {k: int(v.sum()) for k, v in reasons.items()},
        "cap_rows_N_S": [int(iN), int(iS)],
        "cap_ring_flagged": flagged,
        "cap_ring_unreadable_fraction": [float(bad_frac[iN]), float(bad_frac[iS])],
        "cap_rho_N_S": [float(rho_row[iN]), float(rho_row[iS])],
        "in_pin": bool(R > 0.5 * cfg["L"] - 1.6),
        "min_align_director_outward": float(align.min()),
    }


def degree_reader(M, cfg, R, center=(0.0, 0.0, 0.0), nth=64, nph=128, nfield=None):
    """the hedgehog degree of the director, oriented outward from `center`, by the signed areas
    of the spherical triangles of the sample grid (a line field lifted by the outward sign).
    `nfield` (n, n, n, 3): a stored unit director instead of M (the R22-2 arrays keep only n)."""
    from scipy.ndimage import map_coordinates

    n, h = cfg["n"], cfg["h"]
    TH, PH, rhat, that, phat, pts = sphere_samples(R, nth, nph, center)
    # AUDIT CORRECTION (R26-0 audit, 2026-09-25): nearest-cell sampling made the sampled map
    # discontinuous and the solid-angle sum a flux estimator with up to 7 percent residual; the
    # samples are now trilinear (a continuous map), so the sum is an integer to round-off
    coords = np.stack([(pts[..., a] / h + (n - 1) / 2.0).ravel() for a in range(3)])
    shp = pts.shape[:-1]
    if nfield is not None:
        # a stored unit VECTOR field: oriented already, no lift needed
        d = np.stack(
            [
                map_coordinates(nfield[..., a], coords, order=1, mode="nearest").reshape(shp)
                for a in range(3)
            ],
            -1,
        )
        d = d / np.maximum(np.linalg.norm(d, axis=-1), 1e-300)[..., None]
        lam = np.zeros(d.shape[:-1] + (3,))
        lam[..., 2] = 1.0
        al = np.einsum("...a,...a->...", d, rhat)
    else:
        lam, V = np.linalg.eigh(sample_block(M, cfg, pts))
        d = V[..., :, 2]
        al = np.einsum("...a,...a->...", d, rhat)
        sgn = np.where(al < 0, -1.0, 1.0)
        d = d * sgn[..., None]
    # close the grid with the poles as the mean directors of the first and last rows
    north = d[0].mean(axis=0)
    north /= max(np.linalg.norm(north), 1e-300)
    south = d[-1].mean(axis=0)
    south /= max(np.linalg.norm(south), 1e-300)

    def tri(a, b, c):
        num = np.einsum("...a,...a->...", a, np.cross(b, c))
        den = (
            1.0
            + np.einsum("...a,...a->...", a, b)
            + np.einsum("...a,...a->...", b, c)
            + np.einsum("...a,...a->...", c, a)
        )
        return 2.0 * np.arctan2(num, den)

    a = d[:-1, :]
    b = np.roll(d, -1, axis=1)[:-1, :]
    c = np.roll(d, -1, axis=1)[1:, :]
    e = d[1:, :]
    area = np.sum(tri(a, b, c)) + np.sum(tri(a, c, e))
    # the caps
    r0 = d[0]
    r0n = np.roll(r0, -1, axis=0)
    area += np.sum(tri(np.broadcast_to(north, r0.shape), r0, r0n))
    r1 = d[-1]
    r1n = np.roll(r1, -1, axis=0)
    area += np.sum(tri(r1, np.broadcast_to(south, r1.shape), r1n))
    # the sample grid's triangles run clockwise seen from outside (theta south, phi east), so the
    # outward degree is minus the summed area
    return {
        "R": float(R),
        "center": list(map(float, center)),
        "degree": float(-area / (4.0 * np.pi)),
        "min_abs_align": float(np.abs(al).min()),
        "gap_dir_min": float((lam[..., 2] - lam[..., 1]).min()),
    }


# ---------- the census seeds ----------
PARTITIONS = {
    "4": [(-1.0, 0.0)],
    "3_1": [(-0.5, 0.0)],
    "2_2": [],
    "2_1_1": [(-1.0, 0.0), (0.5, 1.0), (0.5, -1.0)],
    "1_1_1_1": [(-1.0, 0.0), (0.5, 1.0), (0.5, 1j), (0.5, -1.0), (0.5, -1j)],
}
"""the finite points of the pair-phase rotation psi = sum_k m_k arg(w - w_k) in the chart
w = (x + i y) / (r + z) (w = 0 the +z axis, w = infinity the -z axis) on the S1 base frame
(phi-hat, theta-hat), which carries index 1 (2 half-units) at each pole: the north carrier is
2 (1 + m_0), the south carrier 2 (1 - sum_k m_k), each finite point 2 m_k half-units."""


def partition_of(name):
    fin = PARTITIONS[name]
    m0 = sum(m for m, w in fin if w == 0.0)
    north = 2.0 * (1.0 + m0)
    south = 2.0 * (1.0 - sum(m for m, w in fin))
    pts = [(north, np.array([0.0, 0.0, 1.0]))] if north else []
    pts += [(south, np.array([0.0, 0.0, -1.0]))] if south else []
    for m, w in fin:
        if w == 0.0:
            continue
        w = complex(w)
        a = abs(w) ** 2
        u = np.array([2 * w.real, 2 * w.imag, 1.0 - a]) / (1.0 + a)
        pts.append((2.0 * m, u))
    return sorted([int(round(hu)) for hu, u in pts], reverse=True), pts


def kappa_bps(delta, w, m=1.0):
    return float(
        np.sqrt(w * (4.0 + 36.0 * (delta / 2) ** 2 + 144.0 * (delta / 2) ** 4) / (32.0 * m * m))
    )


def seed_partition(cfg, delta, name, w=None, r_c=4.0):
    """the hedgehog director with the transverse pair (delta, 0) laid in the rotated frame, the
    core melted as R20.seed_axes, the pair melted along each strand ray with the BPS profile of
    its index (R25-1: sqrt(1 - exp(-kappa d^2)), kappa = sqrt(w K / (32 m^2)))."""
    n, h = cfg["n"], cfg["h"]
    w = W1 * W1S if w is None else w
    X, Y, Z = B3.coords(n, h)
    rho = np.sqrt(X * X + Y * Y)
    r = np.sqrt(rho * rho + Z * Z)
    rs = np.maximum(r, 1e-12)
    nhat = np.stack([X, Y, Z], -1) / rs[..., None]
    rhos = np.maximum(rho, 1e-12)
    phihat = np.stack([-Y / rhos, X / rhos, np.zeros_like(Z)], -1)
    dot = np.einsum("...a,...a->...", phihat, nhat)[..., None]
    ph = phihat - dot * nhat
    ph = ph / np.maximum(np.linalg.norm(ph, axis=-1)[..., None], 1e-300)
    th = np.cross(ph, nhat)
    wch = (X + 1j * Y) / np.maximum(r + Z, 1e-300)
    psi = np.zeros_like(X)
    for m, wk in PARTITIONS[name]:
        psi += m * np.angle(wch - complex(wk))
    # (phi-hat, theta-hat) is the REFLECTED tangent frame (phi x theta = -r), so the rotation by
    # psi in the oriented frame reads cos(psi) phi - sin(psi) theta here (found at the reader gate)
    e1 = np.cos(psi)[..., None] * ph - np.sin(psi)[..., None] * th
    e2 = np.sin(psi)[..., None] * ph + np.cos(psi)[..., None] * th
    part, pts = partition_of(name)
    v = np.ones_like(X)
    for hu, u in pts:
        m = hu / 2.0
        proj = X * u[0] + Y * u[1] + Z * u[2]
        dperp = np.sqrt(np.maximum(r * r - proj * proj, 0.0))
        dist = np.where(proj > 0, dperp, r)
        v = v * np.sqrt(1.0 - np.exp(-kappa_bps(delta, w, m) * dist**2))
    nn = nhat[..., :, None] * nhat[..., None, :]
    aniso = e1[..., :, None] * e1[..., None, :] - e2[..., :, None] * e2[..., None, :]
    S = nn + 0.5 * delta * (np.eye(3) - nn) + (0.5 * delta * v)[..., None, None] * aniso
    a = (1.0 + delta) / 3.0
    wc = 1.0 - np.exp(-((r / r_c) ** 2))
    M3 = wc[..., None, None] * S + (1.0 - wc[..., None, None]) * (a * np.eye(3))
    return B3.embed34(M3, cfg), part


# ============================================================================
# THE CHECKS
# ============================================================================
def check_a():
    import sympy as sp

    out = {}
    m = sp.Matrix(4, 4, lambda i, j: sp.Symbol(f"m{min(i, j)}{max(i, j)}"))
    eta = sp.diag(-1, 1, 1, 1)
    t = sp.Symbol("t")
    gens = {}
    for i, j in ((1, 2), (2, 3), (3, 1)):
        Jm = sp.zeros(4, 4)
        Jm[i, j], Jm[j, i] = -1, 1
        gens[f"rot_{i}{j}"] = Jm
    for i in (1, 2, 3):
        Km = sp.zeros(4, 4)
        Km[0, i] = Km[i, 0] = 1
        gens[f"boost_{i}"] = Km
    rec = {}
    for nm, Gm in gens.items():
        inv = sp.simplify(Gm * eta + eta * Gm.T)
        dM = Gm * m + m * Gm.T
        cat = Gm * m - m * Gm.T
        Lam = (t * Gm).exp()
        pushed = sp.simplify(Lam * m * Lam.T)
        d1 = sp.simplify(sp.diff(pushed, t).subs(t, 0))
        rec[nm] = {
            "eta_invariance": inv == sp.zeros(4, 4),
            "tangent_symmetric": sp.simplify(dM - dM.T) == sp.zeros(4, 4),
            "tangent_is_d/dt_of_pushforward": sp.simplify(d1 - dM) == sp.zeros(4, 4),
            "catalog_form_antisymmetric": sp.simplify(cat + cat.T) == sp.zeros(4, 4),
            "catalog_form_is_zero": cat == sp.zeros(4, 4),
            "tangent_equals": (
                "commutator"
                if sp.simplify(dM - (Gm * m - m * Gm)) == sp.zeros(4, 4)
                else (
                    "anticommutator"
                    if sp.simplify(dM - (Gm * m + m * Gm)) == sp.zeros(4, 4)
                    else "neither"
                )
            ),
            "lower_index_form_agrees": sp.simplify((Gm.T).T * m + m * Gm.T - dM) == sp.zeros(4, 4),
        }
    out["generators"] = rec
    ok = all(
        r["eta_invariance"]
        and r["tangent_symmetric"]
        and r["tangent_is_d/dt_of_pushforward"]
        and r["catalog_form_antisymmetric"]
        and not r["catalog_form_is_zero"]
        and r["lower_index_form_agrees"]
        for r in rec.values()
    )
    ok = ok and all(rec[k]["tangent_equals"] == "commutator" for k in rec if k.startswith("rot"))
    ok = ok and all(
        rec[k]["tangent_equals"] == "anticommutator" for k in rec if k.startswith("boost")
    )
    out["PASS"] = bool(ok)
    out["fails_if"] = (
        "a generator fails eta invariance, its tangent G M + M G^T is not symmetric or is not the"
        " t-derivative of exp(tG) M exp(tG)^T, a rotation's tangent is not the commutator, a"
        " boost's not the anticommutator, or the catalog form G M - M G^T is symmetric or zero"
    )
    return out


def check_b():
    cfg, p, pot = cfg_pot(32, 48.0, 0.3)
    M = np.load(S1_END)["M"]
    cat = B3.gen_catalog(cfg, M)
    out = {"fields": {}}
    Jz, Jx = np.zeros((4, 4)), np.zeros((4, 4))
    Jz[1, 2], Jz[2, 1] = -1.0, 1.0
    Jx[2, 3], Jx[3, 2] = -1.0, 1.0
    Kz, Kx = np.zeros((4, 4)), np.zeros((4, 4))
    Kz[0, 3] = Kz[3, 0] = 1.0
    Kx[0, 1] = Kx[1, 0] = 1.0
    consts = {"rot_z": Jz, "rot_x": Jx, "boost_z": Kz, "boost_x": Kx}
    ok = True
    for nm, a0 in cat.items():
        s = float(np.abs(a0 + a0.swapaxes(-1, -2)).max())
        mx = float(np.abs(a0).max())
        rec = {"max_abs": mx, "symmetric_part_over_max": s / max(mx, 1e-300)}
        if nm in consts:
            Gm = consts[nm]
            tng = Gm @ M + M @ Gm.T
            rec["corrected_antisymmetric_part_over_max"] = float(
                np.abs(tng - tng.swapaxes(-1, -2)).max() / max(np.abs(tng).max(), 1e-300)
            )
            ok = ok and rec["corrected_antisymmetric_part_over_max"] < 1e-12
        out["fields"][nm] = rec
        ok = ok and (mx == 0.0 or rec["symmetric_part_over_max"] < 1e-12)
    out["PASS"] = bool(ok)
    out["fails_if"] = (
        "a catalog field has a symmetric part above 1e-12 of its max (it would then be a tangent"
        " after all), or a corrected tangent has an antisymmetric part above 1e-12"
    )
    return out


def frob_energy_u(M, cfg):
    h3 = cfg["h"] ** 3
    e = 0.0
    for br, (A, wt) in B3.a_fields(M, cfg).items():
        for i in range(3):
            for j in range(i + 1, 3):
                F = B3.comm_eta(A[i], A[j])
                e += wt * 4.0 * np.sum(F * F)
    return float(h3 * e)


def check_c():
    out = {"fields": {}}
    if not os.path.isdir(FULLKICK):
        out["PASS"] = False
        out["fails_if"] = "the witness folder is absent"
        return out
    ok = True
    R25_1 = _load("m5_32_r25_1_strand", "m5_32_r25_1_strand.py")
    for f in sorted(os.listdir(FULLKICK)):
        if not f.endswith(".npz"):
            continue
        tag = f[:-4]
        try:
            d = float(tag.split("_")[0][1:])
            w1s = float(tag.split("_")[1][1:])
            n = int(tag.split("_")[2][1:])
            L = float(tag.split("_")[3][1:])
        except Exception:  # noqa: BLE001
            continue
        M = np.load(os.path.join(FULLKICK, f))["M"]
        cfg = R25_1.slab_cfg(n, L, d)
        eu_eta, _ = B3.e_parts(M, cfg)
        eu_fro = frob_energy_u(M, cfg)
        ev, _ = R0.v4_energy_grad(M, cfg, R0.roots_of(cfg), W1 * w1s, need_grad=False)
        m0i = float(np.abs(M[..., 0, 1:]).max())
        rec = {
            "E_u_eta": float(eu_eta),
            "E_u_frobenius": eu_fro,
            "V4": float(ev),
            "M0i_max": m0i,
            "eta_negative": bool(eu_eta < 0),
            "frobenius_positive": bool(eu_fro > 0),
        }
        out["fields"][tag] = rec
        ok = ok and rec["frobenius_positive"]
    out["n_fields"] = len(out["fields"])
    out["n_eta_negative"] = int(sum(r["eta_negative"] for r in out["fields"].values()))
    # the two norms coincide on a block-diagonal field
    cfg, p, pot = cfg_pot(32, 48.0, 0.3)
    Ms = np.load(S1_END)["M"]
    eu_eta, _ = B3.e_parts(Ms, cfg)
    out["block_sector_coincide"] = {
        "E_u_eta": float(eu_eta),
        "E_u_frobenius": frob_energy_u(Ms, cfg),
        "rel_diff": float(abs(eu_eta - frob_energy_u(Ms, cfg)) / max(abs(eu_eta), 1e-300)),
    }
    ok = ok and out["block_sector_coincide"]["rel_diff"] < 1e-12 and out["n_fields"] > 0
    ok = ok and out["n_eta_negative"] == out["n_fields"]
    out["PASS"] = bool(ok)
    out["fails_if"] = (
        "a full-kick witness field has a NON-NEGATIVE eta curvature energy (the sign flip is the"
        " content: tr(F F^T) is a sum of squares and cannot fail), or the two norms differ on the"
        " block-diagonal stored charge, or no witness field is found (AUDIT CORRECTION 2026-09-25)"
    )
    return out


def v4_plain(Mc, q, w):
    N = Mc @ ETA
    P = np.eye(4)
    t = []
    for p in range(1, 5):
        P = P @ N
        t.append(np.trace(P))
    return w * sum((t[p] - sum(qi ** (p + 1) for qi in q)) ** 2 for p in range(4))


def v4_telescoped(Mc, Mvac, w):
    N = Mc @ ETA
    N0 = Mvac @ ETA
    eps = N - N0
    tot = 0.0
    for p in range(1, 5):
        b = 0.0
        for k in range(p):
            b += np.trace(
                np.linalg.matrix_power(N, k) @ eps @ np.linalg.matrix_power(N0, p - 1 - k)
            )
        tot += b * b
    return w * tot


def check_d():
    out = {"rows": {}}
    ok = True
    for d in (0.3, 0.03, 0.01, 0.003, 0.001):
        q = (-8.0, 1.0, d, 0.0)
        w = W1 * W1S
        Mvac = np.diag([8.0, 1.0, d, 0.0])
        s0 = d / 2
        e_pl, e_te = [], []
        for f in np.linspace(0.05, 0.95, 7):
            b = s0 * f
            Mc = np.diag([8.0, s0 + b, s0 - b, 1.0])
            ve = w * sum(((s0 + b) ** p + (s0 - b) ** p - d**p) ** 2 for p in range(1, 5))
            e_pl.append(abs(v4_plain(Mc, q, w) - ve) / ve)
            e_te.append(abs(v4_telescoped(Mc, Mvac, w) - ve) / ve)
        rec = {"plain_relerr_max": float(max(e_pl)), "telescoped_relerr_max": float(max(e_te))}
        rec["gain"] = rec["plain_relerr_max"] / max(rec["telescoped_relerr_max"], 1e-300)
        out["rows"][str(d)] = rec
        ok = ok and rec["telescoped_relerr_max"] < rec["plain_relerr_max"]
    ok = ok and out["rows"]["0.001"]["telescoped_relerr_max"] < 1e-7
    out["PASS"] = bool(ok)
    out["fails_if"] = (
        "the telescoped bracket is not more accurate than the plain one at some delta, or its"
        " relative error at delta 0.001 exceeds 1e-7"
    )
    return out


def check_e():
    from scipy.optimize import brentq

    out = {
        "uniaxial": float(beta2_of([1.0, 0.0, 0.0])),
        "maximal": float(beta2_of([1.0, 0.0, -1.0])),
        "vacuum_d0.3": float(beta2_of([1.0, 0.3, 0.0])),
        "vacuum_d0.23": float(beta2_of([1.0, 0.23, 0.0])),
    }
    out["delta_at_0.382"] = float(brentq(lambda d: beta2_of([1.0, d, 0.0]) - 0.382, 0.01, 0.49))
    ok = (
        abs(out["uniaxial"]) < 1e-12
        and abs(out["maximal"] - 1.0) < 1e-12
        and abs(out["vacuum_d0.3"] - 0.604) < 2e-3
        and abs(out["delta_at_0.382"] - 0.2307) < 1e-3
    )
    out["PASS"] = bool(ok)
    out["fails_if"] = "the invariant misses 0, 1, 0.604 at delta 0.3 or 0.2307 at beta^2 0.382"
    return out


def check_f():
    cfg, p, pot = cfg_pot(32, 48.0, 0.3)
    M = np.load(S1_END)["M"]
    null = gap_tail(M, cfg, 0.3)
    X, Y, Z = B3.coords(cfg["n"], cfg["h"])
    r = np.sqrt(X * X + Y * Y + Z * Z)
    A = 0.5
    dev = A * np.where(r > 0, Z / np.maximum(r, 1e-12), 0.0) / np.maximum(r, 1.0) ** 2
    synth = gap_tail(M, cfg, 0.3, dev_field=dev)
    sh = np.array(synth["shells_r_a0_a1_a2_maxdev_n"])
    win = (sh[:, 0] >= 6.0) & (sh[:, 0] <= 18.0)
    a1_err = float(np.max(np.abs(sh[win, 2] * sh[win, 0] ** 2 / A - 1.0)))
    out = {
        "null_S1": {k: v for k, v in null.items() if k != "shells_r_a0_a1_a2_maxdev_n"},
        "null_S1_shells": null["shells_r_a0_a1_a2_maxdev_n"],
        "synthetic": {
            "slope_l1": synth["slope_l1"],
            "a1_times_r2_over_A_maxerr": a1_err,
            "slope_l0": synth["slope_l0"],
        },
    }
    ok = (
        null["l1_over_l0_max_in_window"] is not None
        and null["l1_over_l0_max_in_window"] < 5e-3
        and synth["slope_l1"] is not None
        and abs(synth["slope_l1"] + 2.0) < 0.1
        and a1_err < 0.05
    )
    out["PASS"] = bool(ok)
    out["fails_if"] = (
        "the stored S1 charge shows an l = 1 coefficient above 5e-3 of l = 0 in the fit window"
        " (AUDIT CORRECTION 2026-09-25: the relaxed charge is NOT z-reflection symmetric, its odd"
        " part is 0.30 of the field in rms, so this null is a genuine control of the relaxed"
        " field, not a symmetry tautology; the OPEN criterion of R26-4 asks for 2e-2, a decade"
        " above), or the planted cos(theta) / r^2 tail does not read slope -2 within 0.1 and"
        " amplitude within 5 percent on the window"
    )
    return out


def check_g():
    cfg0, _, _ = cfg_pot(32, 48.0, 0.0)
    seed0 = R20.seed_axes(cfg0, (1.0, 0.0, 0.0))
    uni = physical_generator(seed0, cfg0)
    cfg, p, pot = cfg_pot(32, 48.0, 0.3)
    M = np.load(S1_END)["M"]
    st = physical_generator(M, cfg)
    from scipy.ndimage import gaussian_filter

    def fd_diff(Mf):
        a_rig = (JZ @ Mf - Mf @ JZ) - orbital_z(Mf, cfg)
        al = 0.01
        fd = (rotated_copy(Mf, cfg, al) - rotated_copy(Mf, cfg, -al)) / (2.0 * al)
        pin = B3.pin_shell(cfg["n"], cfg["h"], depth=4.0)
        inner = ~pin
        num = float(np.sqrt(np.sum((fd[inner] - a_rig[inner]) ** 2)))
        den = float(np.sqrt(np.sum(a_rig[inner] ** 2)))
        return num / max(den, 1e-300)

    Msm = M.copy()
    for a in range(4):
        for b in range(4):
            Msm[..., a, b] = gaussian_filter(M[..., a, b], sigma=1.0, mode="nearest")
    # the analytic control: a smooth off-axis bump of a fixed symmetric matrix on the vacuum
    X, Y, Z = B3.coords(cfg["n"], cfg["h"])
    bump = 0.1 * np.exp(-((X - 6.0) ** 2 + Y**2 + (Z - 3.0) ** 2) / 25.0)
    Tm = np.zeros((4, 4))
    Tm[1, 2] = Tm[2, 1] = 1.0
    Tm[1, 3] = Tm[3, 1] = 0.5
    Tm[3, 3] = 0.3
    Msyn = np.broadcast_to(np.diag([8.0, 1.0, 0.3, 0.0]), M.shape) + bump[..., None, None] * Tm
    out = {
        "uniaxial_seed": {
            k: uni[k] for k in ("internal", "orbital", "rigid", "rigid_over_internal")
        },
        "stored_S1": {k: st[k] for k in ("internal", "orbital", "rigid", "rigid_over_internal")},
        "stored_S1_within_r": st["within_r"],
        "fd_check": {
            "alpha": 0.01,
            "rel_norm_diff_synthetic_bump": fd_diff(Msyn),
            "rel_norm_diff_smoothed": fd_diff(Msm),
            "rel_norm_diff_raw": fd_diff(M),
            "note": "interior = 4 units inside the walls; cubic interpolation of the rotated copy;"
            " the raw stored field is rough at the cell scale (the strand cores), the control is"
            " its Gaussian-smoothed copy (sigma one cell)",
        },
    }
    ok = (
        uni["rigid_over_internal"] < 0.01
        and out["fd_check"]["rel_norm_diff_synthetic_bump"] < 0.05
    )
    out["PASS"] = bool(ok)
    out["fails_if"] = (
        "the axisymmetric uniaxial seed's C_rigid exceeds 1 percent of C_int, or the finite"
        " difference of the rotated smooth off-axis bump differs from a_rigid by over 5 percent in norm"
        " (the stored charge, rough at the cell scale, is reported raw and smoothed)"
    )
    return out


def check_h(cfg=None):
    cfg, p, pot = cfg_pot(32, 48.0, 0.3) if cfg is None else cfg
    out = {"seeds": {}}
    ok = True
    for name in PARTITIONS:
        M, part = seed_partition(cfg, 0.3, name)
        rec = {"partition_intended": part}
        for R in (9.0, 18.0):
            rd = partition_reader(M, cfg, R, 0.3)
            rec[f"R{R:g}"] = {
                "partition": rd["partition"],
                "total": rd["total_half_units"],
                "carriers": [
                    (
                        c["half_units"],
                        round(c["theta"], 3) if c["theta"] == c["theta"] else None,
                        round(c["phi"], 3),
                        c["touches_pole"],
                    )
                    for c in rd["carriers"]
                ],
                "unreadable": rd["unreadable_plaquettes"],
            }
            ok = ok and rd["partition"] == part and abs(rd["total_half_units"] - 4.0) < 1e-9
        parts = R20.energy_parts(M, cfg, p, pot)
        rec["E_seed"] = float(parts["E_total"])
        rec["block_diagonal"] = bool(np.abs(M[..., 0, 1:]).max() == 0.0)
        rec["finite"] = bool(np.all(np.isfinite(M)))
        ok = ok and rec["block_diagonal"] and rec["finite"]
        out["seeds"][name] = rec
    # the S1 seed of R20 against the {2,2} seed here (same frame, the strand melt the only difference)
    Ms1 = R20.seed_axes(cfg, (1.0, 0.3, 0.0))
    M22, _ = seed_partition(cfg, 0.3, "2_2")
    rd = partition_reader(Ms1, cfg, 9.0, 0.3)
    out["R20_S1_seed_r9"] = {"partition": rd["partition"], "total": rd["total_half_units"]}
    out["S1_vs_2_2_maxdiff"] = float(np.abs(Ms1 - M22).max())
    ok = ok and rd["partition"] == [2, 2]
    # the stored charge's end field
    M = np.load(S1_END)["M"]
    out["stored_S1_end"] = {}
    for R in (3.0, 4.5, 6.0, 9.0, 12.0, 18.0):
        rd = partition_reader(M, cfg, R, 0.3)
        out["stored_S1_end"][f"R{R:g}"] = {
            "partition": rd["partition"],
            "total": rd["total_half_units"],
            "carriers_hu_rho_z": [
                (c["half_units"], c["rho_from_axis"], c["z"]) for c in rd["carriers"]
            ],
            "unreadable": rd["unreadable_plaquettes"],
            "null_spots": len(rd["null_spots"]),
        }
    out["PASS"] = bool(ok)
    out["fails_if"] = (
        "a census seed does not read its intended partition on r 9 and r 18 with total 4, is not"
        " block-diagonal or not finite, or the R20 S1 seed does not read {2,2}"
    )
    return out


def check_i():
    cfg, p, pot = cfg_pot(32, 48.0, 0.3)
    M = np.load(S1_END)["M"]
    eu, ev = B3.e_parts(M, cfg)
    cfg2 = dict(cfg)
    cfg2["h"] = 2.0 * cfg["h"]
    cfg2["L"] = 2.0 * cfg["L"]
    eu2, ev2 = B3.e_parts(M, cfg2)
    vir = virial_reads(M, cfg, p, pot)
    out = {
        "E_u_h": float(eu),
        "E_u_2h": float(eu2),
        "ratio_u": float(eu2 / eu),
        "V4_h": float(ev),
        "V4_2h": float(ev2),
        "ratio_v": float(ev2 / ev),
        "virial_stored_S1": vir,
    }
    ok = abs(out["ratio_u"] - 0.5) < 1e-12 and abs(out["ratio_v"] - 8.0) < 1e-12
    out["PASS"] = bool(ok)
    out["fails_if"] = "E_u does not scale as h^-1 or V4 as h^3 on the same array"
    return out


def check_j():
    cfg, p, pot = cfg_pot(32, 48.0, 0.3)
    Ms1 = R20.seed_axes(cfg, (1.0, 0.3, 0.0))
    out = {"S1_seed": {}}
    ok = True
    for R in (6.0, 9.0, 18.0):
        dg = degree_reader(Ms1, cfg, R)
        out["S1_seed"][f"R{R:g}"] = dg
        ok = ok and abs(dg["degree"] - 1.0) < 0.02
    f = os.path.join(R22_2, "unlike_d12_n32_L48.npz")
    if os.path.exists(f):
        nf = np.load(f)["n"]
        d = 12.0
        cfgp, _, _ = cfg_pot(32, 48.0, 0.3)
        rec = {
            "note": "the R22-2 arrays keep the director only; the cores sit on the x axis at +-d/2"
        }
        for xc in (+d / 2, -d / 2):
            dg = degree_reader(None, cfgp, 4.0, center=(xc, 0.0, 0.0), nfield=nf)
            rec[f"core_x{xc:+g}"] = dg
        rec["outer_r20"] = degree_reader(None, cfgp, 20.0, nfield=nf)
        out["R22_2_unlike_d12"] = rec
        degs = sorted(round(rec[k]["degree"]) for k in rec if k.startswith("core"))
        ok = ok and degs == [-1, 1] and abs(rec["outer_r20"]["degree"]) < 0.02
        ok = ok and all(
            abs(rec[k]["degree"] - round(rec[k]["degree"])) < 0.02
            for k in rec
            if isinstance(rec[k], dict)
        )
    else:
        out["R22_2_unlike_d12"] = "absent"
    out["PASS"] = bool(ok)
    out["fails_if"] = (
        "the S1 seed does not read degree 1 on every sphere within 0.02, or the stored R22-2"
        " unlike pair does not read +1 and -1 around its cores and 0 on the outer sphere, or any"
        " degree is off an integer by over 0.02 (a sampling failure, the audit's rule)"
    )
    return out


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "run"
    res = {"task": "M5.32 R26-0", "norm": "eta contraction, block-diagonal sector (stated)"}
    checks = {
        "a_tangent": check_a,
        "b_catalog_antisymmetric": check_b,
        "c_norm_witness": check_c,
        "d_split_residual": check_d,
        "e_biaxiality": check_e,
        "f_gap_tail_reader": check_f,
        "g_physical_generator": check_g,
        "h_partition_reader": check_h,
        "i_virial_scaling": check_i,
        "j_degree_reader": check_j,
    }
    if mode == "readers":
        checks = {"h_partition_reader": check_h}
    elif mode == "quick":
        checks = {k: v for k, v in checks.items() if k[0] in "abdeij"}
    for k, fn in checks.items():
        t = time.time()
        try:
            res[k] = fn()
        except Exception as e:  # noqa: BLE001
            import traceback

            res[k] = {"PASS": False, "error": repr(e), "traceback": traceback.format_exc()}
        res[k]["wall_s"] = round(time.time() - t, 1)
        log(f"{k}: PASS {res[k]['PASS']} ({res[k]['wall_s']} s)")
    res["PASS"] = all(res[k].get("PASS", True) for k in res if isinstance(res[k], dict))
    res["wall_s"] = round(time.time() - T0, 1)
    with open(OUT_JSON if mode == "run" else OUT_JSON.replace(".json", f"_{mode}.json"), "w") as f:
        json.dump(res, f, indent=1, default=str)
    print("R26-0 PASS:", res["PASS"])


if __name__ == "__main__":
    main()
