"""M5.32 R14-A: the coexistence conjecture as a LINEAR PROGRAM over the frozen basis B+
(ledger 6.3, extended on 2026-09-05 by the author's reply: the constant-coefficient
quadratic jet forms T1..T4, the H-adjoint K_P, the Coulomb-form row, the zigzag-sheet row).

EQUATIONS FIRST
---------------
Action  L = sum_k c_k I_k,  I_k(omega) = A_k + B_k omega + C_k omega^2 (+ D_k omega^4 for the
quartics); energy  E = sum_k c_k (C_k omega^2 - A_k) (Legendre; B drops).  The LP works in the
ENERGY orientation  chat_k := -c_k, so that E_static = sum_k chat_k s_k with s_k := A_k (the
static density) and E_kin = omega^2 sum_k chat_k q_k with q_k := -C_k.  The certified action
is chat_I1 = 4 (E_u = 4 h^3 sum I1), chat_V4 = 1 (the potential is not a variable).  Every
basis element is measured on SAVED FIELDS as the pair (s_k, q_k), h^3-weighted on the certified
sym stencil; nothing is fitted.

The basis B+ (24 free coefficients, chat_I1 = 4 fixed):
  C0/C1  I2, I3, I4, I5, I6;  parity-odd E1..E3;  C2 insertions I1_h, J1, J2, Pgrad  (R1's thirteen)
  C4     K_T (R7);  K_lambda, R_{eta M eta}, R_hcov, K_P^h (the 09-03 / 09-05 entrants);
         T1 = eta^{mu nu} tr(A_mu eta A_nu eta),  T2 = eta^{mu nu} tr(A_mu eta) tr(A_nu eta),
         T3 = div_b eta^{bd} div_d with div^b = sum_mu (A_mu)^{mu b},  T4 = sum_{mu nu} (A_mu)^{mu nu} tr(A_nu eta)
         (the covariant constant-coefficient quadratic jet forms; T5 = T3 + R_eta is a total
         derivative away from T3 and is not a separate element; Q_F = T1 - T3 is in the span)
  C5/C6  (F.F)^2 = I1^2 (Q_I1sq), the two non-F quartics C6a, C6b
Excluded and why: R_eta (EL == 0, empty), R_{M^-1} (undefined on the singular certified
vacuum), the plain K_P (indefinite off-shell; K_P^h is the same term on the orbit), the
registry controls I1_frob and I3_mixed_eta (non-covariant / duplicates).

The witness rows (all linear in chat; every number carries its field, box, stencil):
  (i)   TAIL      a_k := the shell-plateau energy per unit r of s_k on the n48 L72 (3000 it)
                  relaxed hedgehog, averaged over r in [12, 21] (the rigid far field): the
                  coefficient of the L-linear divergence.  Row: |sum chat_k a_k| <= eps_tail
                  (eps_tail = 3 x the certified plateau |4 a_I1|, or 1e-3 of the largest |a_k|).
  (ii)  TWIST     S_k(w) := the static energy per area of the (1,2) twist psi(z), 0 -> 1 over w,
        SHEET     on a 4 x 4 x n_z slab (n_z = 256, L_z = 24), w in {1.5, 3, 6}: sum chat_k
                  S_k(w_min) >= eps_sheet and >= 0 at the other w (the static cost must be
                  positive and grow as w -> 0: only then does the fixed-J functional, whose
                  inertia on this sheet is ~ 1/w, have a minimizer).  Zeros certified by the
                  n_z = 128 -> 256 ratio (an artifact drops by 4).
  (iii) ZIGZAG    the same for the (2,3)-split ramp (delta, 0) -> (delta/2, delta/2) over w in a
        SHEET     fixed frame (the mechanism W3 measured).
  (iv)  PAIR      on the R3 arm-i ansatz (n32 L48, rc = 1, m = 0.1, same-sign centers at +-d/2,
                  d in {8, 10, 12, 16}) and on the R3 arm-ii RELAXED lam = 0 pairs (d in {10, 14,
                  18, 24}): attraction = E_int increases with d: sum chat_k [E_k(d_max) - E_k(d_min)] >= eps.
  (v)   COULOMB   the certified like-charge 3x3 pair fields (m5_21_4 ladder it400, d = 12, 18,
                  24, embedded 4x4): repulsion sum chat_k [E_k(12) - E_k(24)] >= eps_C, and the
                  1/d FORM: with rho_ij = 1/d_i - 1/d_j and a free strength K,
                  |Delta_ij - K rho_ij| <= 0.1 K_cert rho_ij for the three pairs, K >= K_cert/2.
  (vi)  POSITIVE  every channel's omega^2 energy coefficient >= 0: the six vacuum tangents,
                  the seven hedgehog tangents (true tangents X M + M X^T, unit norm, plus the
                  raw local clock), the local clock on both sheets.
  (vii) UV/G5    the necessary boundedness condition: for a plane-wave perturbation of the
                  vacuum, M = M0 + eps cos(k z) X (X a random symmetric direction, k = 2 pi m / L_z
                  with m in {1, 4, 16} on a 4 x 4 x 64 slab, eps in {0.03, 0.3, 3}: the quadratic
                  and the quartic regimes), the cycle-averaged static energy density must be >= 0
                  (total derivatives average out; a negative growing gradient energy is unbounded
                  below since V4 is bounded per cell): sum chat_k <s_k> >= 0, 300 samples.  Without
                  this row the LP admits negative coefficients on manifestly positive terms.
Conditioning: null columns (the parity-odd E1..E3 vanish on every row) are removed; per row,
entries below 1e-7 of the row's largest entry are set to zero (roundoff cannot buy feasibility).
Normalization: NONE for the class theorem (the LP is a cone question: feasible or not), and
sum |chat_k| <= 1 for the O(1)-coefficient (parsimony) reading; a ladder {1, 10, 100, 1e4} shows
where the two readings meet.  Solvers: HiGHS dual simplex and interior point.  On infeasibility
the Farkas certificate y >= 0 is obtained from the explicit dual LP and RE-CHECKED IN EXACT
RATIONAL ARITHMETIC on the rounded rows: with r = A^T y and beta = b^T y, the inequality
beta + B ||r||_inf < 0 proves that no chat with sum |chat| <= B satisfies the rows (B = 1e4).
Sub-basis LPs name the class that carries or blocks feasibility; the LP WITHOUT the
Coulomb-form row tests the author's 09-05 prediction (CONE_FEASIBLE via Q_F without it).

Stages: measure (writes data/m5_32_r14_a_rows.json, ~15 min), lp (writes
data/m5_32_r14_a_lp.json), all.
"""
from __future__ import annotations
import importlib.util
import json
import os
import sys
import time
from fractions import Fraction
import numpy as np

ARGV = list(sys.argv)          # captured BEFORE the imports: m5_32_r7_a_kt_form wipes sys.argv at import
HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.dirname(HERE)
DATA = os.path.join(RES, "data")
CK = os.path.join(RES, "checkpoints")
ROWS_JSON = os.path.join(DATA, "m5_32_r14_a_rows.json")
LP_JSON = os.path.join(DATA, "m5_32_r14_a_lp.json")
T0 = time.time()


def log(m):
    print(f"[{time.time() - T0:8.1f}s] {m}", flush=True)


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


T14 = _load("m5_32_r14_terms", "m5_32_r14_terms.py")
L0 = T14.L0
B3 = T14.B3
EXT = _load("m5_32_terms_ext", "m5_32_terms_ext.py")
KT7 = _load("m5_32_r7_a_kt_form", "m5_32_r7_a_kt_form.py")
R8 = _load("m5_32_r8_a_quartics", "m5_32_r8_a_quartics.py")
C13 = _load("m5_32_r13w_common", "m5_32_r13w_common.py")
R3 = _load("m5_32_r3_i_ansatz", "m5_32_r3_i_ansatz.py")
B8 = C13.B8
ETA = T14.ETA
P8 = L0.default_params(s=-1.0, g=8.0, delta=0.3)


# ============================================================ the quadratic jet forms
def _div(A):
    return sum(A[mu][..., mu, :] for mu in range(4))


def d_T1(A):
    Ae = np.einsum("m...ab,bc->m...ac", A, ETA)
    return sum(ETA[m, m] * np.einsum("...ab,...ba->...", Ae[m], Ae[m]) for m in range(4))


def d_T2(A):
    tr = np.einsum("m...ab,ba->m...", A, ETA)
    return sum(ETA[m, m] * tr[m] ** 2 for m in range(4))


def d_T3(A):
    dv = _div(A)
    return np.einsum("...b,bd,...d->...", dv, ETA, dv)


def d_T4(A):
    tr = np.einsum("m...ab,ba->m...", A, ETA)
    return sum(A[mu][..., mu, nu] * tr[nu] for mu in range(4) for nu in range(4))


# ============================================================ the basis
class Elem:
    """name, group, lag(A, M) = the Lagrangian density (A the full 4-jet), degree in omega,
    plus optional overrides for the static and kinetic reads."""

    def __init__(self, name, group, lag, degree=2, s_override=None, q_override=None, definition=""):
        self.name, self.group, self.lag, self.degree = name, group, lag, degree
        self.s_override, self.q_override = s_override, q_override
        self.definition = definition


def build_basis():
    B = []
    for nm in ("I1", "I2", "I3", "I4", "I5", "I6"):
        T = L0.REGISTRY[nm]
        B.append(Elem(nm, "C0C1", (lambda A, M, T=T: T.density(A, M, P8)), definition=T.definition))
    for nm in sorted(EXT.REGISTRY_EXT):
        T = EXT.REGISTRY_EXT[nm]
        grp = "eps" if nm.startswith("E") else "C2"
        B.append(Elem(nm, grp, (lambda A, M, T=T: T.density(A, M, P8)), definition=T.definition))
    B.append(Elem("K_T", "C4", (lambda A, M: KT7.kt_density_np(A, M, P8)),
                  definition="K_T (R7): (1/2) sum_mu eta^mumu [tr(h A_mu h A_mu) - tr(eta A_mu eta A_mu)], h = h_cov"))
    for nm, fn in (("T1", d_T1), ("T2", d_T2), ("T3", d_T3), ("T4", d_T4)):
        B.append(Elem(nm, "C4", (lambda A, M, fn=fn: fn(A)), definition={"T1": "eta^{mu nu} tr(A_mu eta A_nu eta)",
                                                                           "T2": "eta^{mu nu} tr(A_mu eta) tr(A_nu eta)",
                                                                           "T3": "div_b eta^{bd} div_d, div^b = sum_mu (A_mu)^{mu b}",
                                                                           "T4": "sum_{mu nu} (A_mu)^{mu nu} tr(A_nu eta)"}[nm]))
    # entrants (E-orientation): lag = s - kin(A_0)
    for nm in ("K_lambda", "R_etaMeta", "R_hcov", "K_P_h"):
        e = T14.ENTRANTS[nm]
        if nm == "K_lambda":
            B.append(Elem(nm, "C4", (lambda A, M, e=e: e.static_fn(A, M) - e.kin_fn(A[0], M)),
                          s_override=(lambda M, cfg: T14.klam_static_fd(M, cfg)), definition=e.definition))
        elif nm.startswith("R_"):
            kind = nm[2:]
            B.append(Elem(nm, "C4", (lambda A, M, kind=kind: T14.rg_density(A, M, kind)), definition=e.definition))
        else:
            B.append(Elem(nm, "C4", (lambda A, M, e=e: e.static_fn(A, M) - e.kin_fn(A[0], M)), definition=e.definition))
    B.append(Elem("Q_I1sq", "C5", (lambda A, M: R8.d_I1(A) ** 2), degree=4, definition="(F.F)^2 = I1^2 (class C5)"))
    B.append(Elem("C6a", "C6", (lambda A, M: R8.d_C6a(A)), degree=4, definition=R8.QUARTICS["Q_C6a"][0]))
    B.append(Elem("C6b", "C6", (lambda A, M: R8.d_C6b(A)), degree=4, definition=R8.QUARTICS["Q_C6b"][0]))
    return B


BASIS = build_basis()
NAMES = [e.name for e in BASIS]
FIXED = {"I1": 4.0}


# ============================================================ measurement helpers
def jets(M, cfg, br, a0=None, om=0.0):
    A = np.zeros((4,) + M.shape)
    for ax in range(3):
        A[1 + ax] = B3.d1(M, ax, cfg["h"], br)
    if a0 is not None and om != 0.0:
        A[0] = om * a0
    return A


def static_cells(e, M, cfg):
    """per-cell h^3-weighted static E-density s_k (sym stencil)."""
    if e.s_override is not None:
        return e.s_override(M, cfg)
    d = 0.0
    for br, wt in B3.branches(cfg["stencil"]):
        d = d + wt * e.lag(jets(M, cfg, br), M)
    return cfg["h"] ** 3 * d


def kin_total(e, M, cfg, a0):
    """q_k = -C_k, h^3-weighted total (3-point read; 5-point for the quartics)."""
    if e.q_override is not None:
        return e.q_override(M, cfg, a0)
    h3 = cfg["h"] ** 3

    def Lw(om):
        t = 0.0
        for br, wt in B3.branches(cfg["stencil"]):
            t = t + wt * float(np.sum(e.lag(jets(M, cfg, br, a0, om), M)))
        return h3 * t
    if e.degree == 2:
        l0, lp, lm = Lw(0.0), Lw(1.0), Lw(-1.0)
        C = 0.5 * (lp + lm) - l0
        return -C, 0.0
    # degree 4: L(om) = A + B om + C om^2 + D om^4 (+ odd cubic): 5 points
    oms = (-2.0, -1.0, 0.0, 1.0, 2.0)
    vals = np.array([Lw(o) for o in oms])
    V = np.vander(np.array(oms), 5, increasing=True)
    coef = np.linalg.solve(V, vals)                     # A, B, C, cubic, D
    return -float(coef[2]), -float(coef[4])


def measure_field(M, cfg, tag, a0=None, shells=False):
    out = {}
    X = None
    if shells:
        Xg, Yg, Zg = B3.coords(cfg["n"], cfg["h"])
        r = np.sqrt(Xg * Xg + Yg * Yg + Zg * Zg)
        edges = np.arange(0.0, cfg["L"] / 2 * np.sqrt(3) + cfg["h"], cfg["h"])
        idx = np.digitize(r, edges) - 1
    for e in BASIS:
        t0 = time.time()
        sc = static_cells(e, M, cfg)
        rec = {"s": float(np.sum(sc))}
        if shells:
            prof = np.bincount(idx.ravel(), weights=sc.ravel(), minlength=len(edges))[: len(edges) - 1] / cfg["h"]
            rec["shell_per_r"] = [float(v) for v in prof]
        if a0 is not None:
            q, q4 = kin_total(e, M, cfg, a0)
            rec["q"] = float(q)
            if e.degree == 4:
                rec["q4"] = float(q4)
        out[e.name] = rec
        log(f"    {tag:36s} {e.name:10s} s {rec['s']:14.6e}" + (f" q {rec['q']:12.5e}" if a0 is not None else "") + f"  [{time.time() - t0:.1f}s]")
    return out


# ============================================================ the witness fields
def true_tangents(M, cfg):
    w = B3.envelope(cfg)[..., None, None]
    lamv, Vv = np.linalg.eigh(M[..., 1:4, 1:4])
    Jz = np.zeros((4, 4)); Jz[1, 2], Jz[2, 1] = -1.0, 1.0
    Jx = np.zeros((4, 4)); Jx[2, 3], Jx[3, 2] = -1.0, 1.0
    Kz = np.zeros((4, 4)); Kz[0, 3] = Kz[3, 0] = 1.0
    Kx = np.zeros((4, 4)); Kx[0, 1] = Kx[1, 0] = 1.0

    def local_rot(vhat):
        W = np.zeros(vhat.shape[:-1] + (4, 4))
        n1, n2, n3 = vhat[..., 0], vhat[..., 1], vhat[..., 2]
        W[..., 1, 2], W[..., 1, 3] = -n3, n2
        W[..., 2, 1], W[..., 2, 3] = n3, -n1
        W[..., 3, 1], W[..., 3, 2] = -n2, n1
        return W
    Xs = {"clock_local": local_rot(Vv[..., :, 2]), "plane_1d": local_rot(Vv[..., :, 0]),
          "rot_z": np.broadcast_to(Jz, M.shape), "rot_x": np.broadcast_to(Jx, M.shape),
          "boost_z": np.broadcast_to(Kz, M.shape), "boost_x": np.broadcast_to(Kx, M.shape)}
    out = {"a0_local_raw": C13.a0_local(M)}
    for nm, X in Xs.items():
        t = w * (X @ M + M @ X.swapaxes(-1, -2))
        out[nm] = t / np.sqrt(np.sum(t * t))
    return out


def slab(kind, w, nz, Lz=24.0):
    cfg = B3.base_cfg(s=-1.0, g=8.0, n=nz, L=Lz)
    h = cfg["h"]
    cfg = dict(cfg)
    cfg["nx"] = 4
    wc = max(1, int(round(w / h)))
    z = np.zeros(nz)
    k0 = nz // 2 - wc // 2
    z[k0:k0 + wc] = np.linspace(0, 1.0, wc + 1)[1:]
    z[k0 + wc:] = 1.0
    prof = np.broadcast_to(z[None, None, :], (4, 4, nz))
    if kind == "twist":
        Rn = B8.rot_field(B8.G3, prof)
        M = np.einsum("...ab,bc,...dc->...ad", Rn, B3.vac4(cfg), Rn)
    else:
        d = 0.3
        D = np.zeros((4, 4, nz, 4, 4))
        D[..., 0, 0] = 8.0
        D[..., 1, 1] = 1.0
        D[..., 2, 2] = d - d * prof / 2
        D[..., 3, 3] = d * prof / 2
        M = D
    area = 4 * 4 * h * h
    return M, cfg, area, wc * h


def stage_measure():
    rows = {"basis": [{"name": e.name, "group": e.group, "degree": e.degree, "definition": e.definition} for e in BASIS],
            "fixed": FIXED}
    # (i) TAIL on the matched-maturity hedgehog n48 L72 (3000 it) with shells; n32 12000 as the robustness field
    M48, cfg48, _ = C13.seed_hedgehog(48, 72.0)
    rows["tail_n48_L72_it3000"] = measure_field(M48, cfg48, "tail n48 L72 3000it", shells=True)
    rows["tail_n48_L72_it3000"]["_shell_edges_h"] = cfg48["h"]
    M32 = np.load(C13.R10_SEED); cfg32 = C13.cfg_of(32, 48.0)
    rows["tail_n32_L48_it12000"] = measure_field(M32, cfg32, "tail n32 L48 12000it", shells=True)
    rows["tail_n32_L48_it12000"]["_shell_edges_h"] = cfg32["h"]
    json.dump(rows, open(ROWS_JSON, "w"), indent=1); log("checkpoint: tail rows written")
    # (vi) POSITIVITY channels: vacuum tangents (uniform n = 4 box), hedgehog tangents (n32 12000 it)
    cfgv = B3.base_cfg(s=-1.0, g=8.0, n=4, L=6.0)
    Mv = np.broadcast_to(B3.vac4(cfgv), (4, 4, 4, 4, 4)).copy()
    chans = {}
    gens = {}
    for i in range(1, 4):
        K = np.zeros((4, 4)); K[0, i] = K[i, 0] = 1.0; gens[f"vac_boost_0{i}"] = K
    for (i, j) in ((1, 2), (1, 3), (2, 3)):
        J = np.zeros((4, 4)); J[i, j], J[j, i] = -1.0, 1.0; gens[f"vac_rot_{i}{j}"] = J
    for nm, X in gens.items():
        a0 = np.broadcast_to(X @ B3.vac4(cfgv) + B3.vac4(cfgv) @ X.T, Mv.shape)
        chans[nm] = {k: v.get("q", 0.0) for k, v in measure_field(Mv, cfgv, nm, a0=a0).items()}
    tt = true_tangents(M32, cfg32)
    for nm, a0 in tt.items():
        chans["hedgehog_" + nm] = {k: v.get("q", 0.0) for k, v in measure_field(M32, cfg32, "hedgehog " + nm, a0=a0).items()}
    rows["channels"] = chans
    json.dump(rows, open(ROWS_JSON, "w"), indent=1); log("checkpoint: channel rows written")
    # (ii)/(iii) SHEETS: twist and zigzag slabs at n_z = 128 and 256, w in {1.5, 3, 6}
    sheets = {}
    for kind in ("twist", "zigzag"):
        for nz in (128, 256):
            for w in (1.5, 3.0, 6.0):
                M, cfg, area, wreal = slab(kind, w, nz)
                a0 = C13.a0_local(M)
                m = measure_field(M, cfg, f"{kind} nz{nz} w{w:g}", a0=a0)
                sheets[f"{kind}|nz{nz}|w{w:g}"] = {"w": wreal, "h": cfg["h"], "area": area,
                                                    "S_per_area": {k: v["s"] / area for k, v in m.items()},
                                                    "q_per_area": {k: v["q"] / area for k, v in m.items()}}
    rows["sheets"] = sheets
    json.dump(rows, open(ROWS_JSON, "w"), indent=1); log("checkpoint: sheet rows written")
    # (iv) PAIRS: the R3 arm-i ansatz at n32 L48 (rc = 1, m = 0.1), separations d (centers +-d/2)
    cfgp = C13.cfg_of(32, 48.0)
    M0 = B3.vac4(cfgp)
    pairs = {}
    for d in (8.0, 10.0, 12.0, 16.0):
        M = R3.ansatz(M0, [(0.0, 0.0, +d / 2), (0.0, 0.0, -d / 2)], 0.1, 1.0, 32, 48.0)
        pairs[f"ansatz_same_d{d:g}"] = {k: v["s"] for k, v in measure_field(M, cfgp, f"ansatz d{d:g}").items()}
    Ms = R3.ansatz(M0, [(0.0, 0.0, 0.0)], 0.1, 1.0, 32, 48.0)
    pairs["ansatz_single"] = {k: v["s"] for k, v in measure_field(Ms, cfgp, "ansatz single").items()}
    for lam in ("lam0", "lam1"):
        for d in (10, 14, 18, 24):
            f = os.path.join(DATA, "m5_32_r3_ii", f"{lam}_dr1_same_d{d}_n32.npz")
            M = np.load(f)["M"]
            pairs[f"relaxed_{lam}_same_d{d}"] = {k: v["s"] for k, v in measure_field(M, cfgp, f"relaxed {lam} d{d}").items()}
    rows["pairs"] = pairs
    json.dump(rows, open(ROWS_JSON, "w"), indent=1); log("checkpoint: pair rows written")
    # (v) COULOMB: the certified like-charge 3x3 pairs (m5_21_4 ladder it400), embedded 4x4
    coul = {}
    for d in (12, 18, 24):
        M3 = np.load(os.path.join(DATA, f"m5_21_4_lad_same_d{d}_n32_it400.npz"))["M"]
        M = B3.embed34(M3, cfgp)
        coul[f"coulomb_same_d{d}"] = {k: v["s"] for k, v in measure_field(M, cfgp, f"coulomb d{d}").items()}
    rows["coulomb"] = coul
    rows["wall_s"] = round(time.time() - T0, 1)
    json.dump(rows, open(ROWS_JSON, "w"), indent=1); log(f"measure done -> {ROWS_JSON}")
    return rows


def stage_uv(nsamp=300, seed=1414):
    """(vii): cycle-averaged static densities of plane-wave perturbations of the vacuum."""
    rows = json.load(open(ROWS_JSON))
    rng = np.random.default_rng(seed)
    nz, Lz = 64, 24.0
    cfg = B3.base_cfg(s=-1.0, g=8.0, n=nz, L=Lz)
    h = cfg["h"]
    z = (np.arange(nz) - (nz - 1) / 2.0) * h
    M0 = B3.vac4(cfg)
    uv = []
    for i in range(nsamp):
        X = rng.normal(size=(4, 4)); X = 0.5 * (X + X.T); X /= np.sqrt(np.sum(X * X))
        m = int(rng.choice([1, 4, 16]))
        eps = float(rng.choice([0.03, 0.3, 3.0]))
        k = 2 * np.pi * m / Lz
        prof = eps * np.cos(k * z + rng.uniform(0, 2 * np.pi))
        M = M0[None, None, None] + prof[None, None, :, None, None] * X[None, None, None]
        M = np.broadcast_to(M, (4, 4, nz, 4, 4)).copy()
        rec = {"m": m, "eps": eps, "X": [float(v) for v in X.ravel()], "avg": {}}
        for e in BASIS:
            sc = static_cells(e, M, cfg)
            rec["avg"][e.name] = float(np.mean(sc[1:-1, 1:-1, 4:-4]) / h ** 3)     # interior cells, per volume
        uv.append(rec)
        if i % 50 == 0:
            log(f"  uv sample {i}: m {m} eps {eps}")
    rows["uv"] = uv
    json.dump(rows, open(ROWS_JSON, "w"), indent=1)
    log(f"uv rows written ({nsamp})")
    return rows


# ============================================================ the exact UV quadratic forms
UV_DIRS = ((1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (1, 0, 1), (0, 1, 1), (1, 1, 1))


def sym_basis10():
    out = []
    for a in range(4):
        for b in range(a, 4):
            E = np.zeros((4, 4)); E[a, b] = E[b, a] = 1.0
            if a != b:
                E /= np.sqrt(2.0)
            out.append(E)
    return out


def uv_quadratic_forms(n=24, L=24.0, eps=0.01, m=2):
    """for each wave direction (m_x, m_y, m_z) and each two-derivative element k, the 10x10
    matrix Q_k(dir) of the cycle-averaged quadratic-regime static energy density of the plane
    wave M = M0 + eps cos(k.x) X (box-periodic wave, k = 2 pi m (m_x, m_y, m_z) / L; interior
    cells averaged), by polarization on the 10 symmetric basis directions: Q_ij = [q(X_i + X_j)
    - q(X_i) - q(X_j)] / 2 with q(X) = <s>(eps X) / eps^2.  Degree-4 elements are excluded here
    (their exact rows are the random-jet quartic rows)."""
    cfg = B3.base_cfg(s=-1.0, g=8.0, n=n, L=L)
    h = cfg["h"]
    X, Y, Z = B3.coords(n, h)
    M0 = B3.vac4(cfg)
    basis = sym_basis10()
    elems = [e for e in BASIS if e.degree == 2]
    sl = (slice(2, -2), slice(2, -2), slice(2, -2))
    out = {}
    for d in UV_DIRS:
        kx, ky, kz = (2 * np.pi * m / L) * np.array(d, float)
        phase = np.cos(kx * X + ky * Y + kz * Z)

        def q_of(Xm):
            M = M0[None, None, None] + (eps * phase)[..., None, None] * Xm[None, None, None]
            M = np.ascontiguousarray(np.broadcast_to(M, (n, n, n, 4, 4)))
            return np.array([float(np.mean(static_cells(e, M, cfg)[sl]) / h ** 3) / eps ** 2 for e in elems])
        qi = [q_of(Xi) for Xi in basis]
        Q = np.zeros((len(elems), 10, 10))
        for i in range(10):
            Q[:, i, i] = qi[i]
            for j in range(i + 1, 10):
                qij = q_of(basis[i] + basis[j])
                Q[:, i, j] = Q[:, j, i] = 0.5 * (qij - qi[i] - qi[j])
        out["|".join(str(x) for x in d)] = {e.name: Q[k].tolist() for k, e in enumerate(elems)}
        log(f"  uv form direction {d}: done")
    return out


def stage_uv2(nsamp=150, seed=2929):
    """(viii) CROSSED plane waves (added 2026-09-05 after the overnight audit): a single plane wave
    has F = 0 identically, so the single-wave rows are blind to every curvature-built term (I1..I6,
    the C2 insertions, the F-quartics); two waves along two random directions with random symmetric
    amplitudes X, Y on an n = 16 box (L = 24) sample the F-built sector's UV positivity: the cycle-
    averaged static density of EVERY basis element, amplitudes eps in {0.03, 0.3, 3}."""
    rows = json.load(open(ROWS_JSON))
    rng = np.random.default_rng(seed)
    n, L = 16, 24.0
    cfg = B3.base_cfg(s=-1.0, g=8.0, n=n, L=L)
    h = cfg["h"]
    X0, Y0, Z0 = B3.coords(n, h)
    M0 = B3.vac4(cfg)
    uv2 = []
    skipped = 0
    for i in range(nsamp):
        eps = float(rng.choice([0.03, 0.3, 1.0]))          # eps 3 with two waves drives M eta complex somewhere
        ks = []
        for _ in range(2):
            m = rng.integers(1, 4, size=3) * rng.choice([-1, 1], size=3)
            ks.append(2 * np.pi * m / L)
        Xs = []
        for _ in range(2):
            Xm = rng.normal(size=(4, 4)); Xm = 0.5 * (Xm + Xm.T); Xm /= np.sqrt(np.sum(Xm * Xm)); Xs.append(Xm)
        phase = [np.cos(k[0] * X0 + k[1] * Y0 + k[2] * Z0 + rng.uniform(0, 2 * np.pi)) for k in ks]
        M = M0[None, None, None] + eps * (phase[0][..., None, None] * Xs[0] + phase[1][..., None, None] * Xs[1])
        M = np.ascontiguousarray(M)
        rec = {"eps": eps, "k1": ks[0].tolist(), "k2": ks[1].tolist(), "avg": {}}
        try:
            for e in BASIS:
                sc = static_cells(e, M, cfg)
                rec["avg"][e.name] = float(np.mean(sc[2:-2, 2:-2, 2:-2]) / h ** 3)
        except ValueError as ex:                          # a complex spectrum of M eta on this sample
            skipped += 1
            continue
        uv2.append(rec)
        if i % 50 == 0:
            log(f"  uv2 sample {i}: eps {eps}")
    rows["uv2"] = uv2
    json.dump(rows, open(ROWS_JSON, "w"), indent=1)
    log(f"uv2 rows written ({len(uv2)} kept, {skipped} skipped for a complex spectrum)")
    return rows


def stage_uvq():
    rows = json.load(open(ROWS_JSON))
    rows["uvq"] = uv_quadratic_forms()
    json.dump(rows, open(ROWS_JSON, "w"), indent=1)
    log("uv quadratic forms written")
    return rows


def uv_min_eig(rows, x):
    """min eigenvalue of sum_k chat_k Q_k(dir) over the directions; returns per-direction (min eig, eigvec)."""
    res = {}
    for d, Qs in rows["uvq"].items():
        Q = sum(x.get(k, 0.0) * np.array(Qk) for k, Qk in Qs.items())
        w, V = np.linalg.eigh(Q)
        res[d] = (float(w[0]), V[:, 0].tolist())
    return res


# ============================================================ the LP
def build_rows(rows, include_coulomb_form=True, include_zigzag=True, include_relaxed_pair=True, include_uv=True, extra=()):
    """returns (A_ub, b_ub, labels) over x = chat (all basis names, I1 included; I1 is then fixed
    by an equality) plus the Coulomb strength K as the last variable."""
    names = NAMES
    nv = len(names) + 1                      # + K
    kI = names.index("I1")
    A, b, lab = [], [], []

    def row(coefs, rhs, label):
        v = np.zeros(nv)
        for k, c in coefs.items():
            v[names.index(k)] = c
        A.append(v); b.append(rhs); lab.append(label)
        return v
    # (i) tail: plateau per unit r over r in [12, 21] on n48 L72
    t = rows["tail_n48_L72_it3000"]
    h = t["_shell_edges_h"]
    a = {}
    for k in names:
        prof = np.array(t[k]["shell_per_r"])
        rr = (np.arange(len(prof)) + 0.5) * h
        sel = (rr >= 12.0) & (rr <= 21.0)
        a[k] = float(np.mean(prof[sel]))
    # the tolerance is the certified action's OWN residual plateau (its unconverged 3000-it tail),
    # three times over; the earlier relative-to-max clause let K_P^h's giant tail entry set a
    # bound sixteen times looser (found on the first refined vertex, 2026-09-05 01:30 UTC)
    eps_tail = 3.0 * abs(FIXED["I1"] * a["I1"])
    row({k: a[k] for k in names}, eps_tail, "tail: sum chat a_k <= eps")
    row({k: -a[k] for k in names}, eps_tail, "tail: -sum chat a_k <= eps")
    # the SECOND tail field (n32 L48, 12000 it, r in [10, 20]): a cancellation of L-divergent tails
    # that is structural must hold on every relaxed field, not on the one the row was measured on
    # (added 2026-09-05 01:45 UTC after the first refined vertex cancelled two tails of 277 to 0.3
    # with a ratio tuned to the n48 field)
    t2 = rows["tail_n32_L48_it12000"]
    h2 = t2["_shell_edges_h"]
    a2 = {}
    for k in names:
        prof = np.array(t2[k]["shell_per_r"])
        rr = (np.arange(len(prof)) + 0.5) * h2
        sel = (rr >= 10.0) & (rr <= 20.0)
        a2[k] = float(np.mean(prof[sel]))
    eps_tail2 = 3.0 * abs(FIXED["I1"] * a2["I1"])
    row({k: a2[k] for k in names}, eps_tail2, "tail (n32 12000 it): sum chat a_k <= eps")
    row({k: -a2[k] for k in names}, eps_tail2, "tail (n32 12000 it): -sum chat a_k <= eps")
    # (ii) twist sheet (n_z 256), (iii) zigzag sheet
    for kind in (("twist", "zigzag") if include_zigzag else ("twist",)):
        for w in (1.5, 3.0, 6.0):
            S = rows["sheets"][f"{kind}|nz256|w{w:g}"]["S_per_area"]
            S128 = rows["sheets"][f"{kind}|nz128|w{w:g}"]["S_per_area"]
            # certified zeros: an entry whose n_z 128 -> 256 ratio is > 2.5 is an artifact -> set to 0
            Sc = {}
            for k in names:
                ratio = abs(S128[k]) / max(abs(S[k]), 1e-300)
                Sc[k] = 0.0 if (abs(S[k]) < 1e-12 or ratio > 2.5) else S[k]
            eps = 1e-3 * max(abs(v) for v in Sc.values()) if w == 1.5 else 0.0
            row({k: -Sc[k] for k in names}, -eps, f"{kind} sheet w={w:g}: sum chat S_k >= {'eps' if eps else '0'}")
    # (iv) dressed pairs
    P = rows["pairs"]
    row({k: -(P["ansatz_same_d16"][k] - P["ansatz_same_d8"][k]) for k in names},
        -1e-3 * max(abs(P["ansatz_same_d16"][k] - P["ansatz_same_d8"][k]) for k in names),
        "dressed pair (ansatz): E(16) - E(8) >= eps (attraction)")
    if include_relaxed_pair:
        row({k: -(P["relaxed_lam0_same_d24"][k] - P["relaxed_lam0_same_d10"][k]) for k in names},
            -1e-3 * max(abs(P["relaxed_lam0_same_d24"][k] - P["relaxed_lam0_same_d10"][k]) for k in names),
            "dressed pair (relaxed lam0): E(24) - E(10) >= eps (attraction)")
    # (v) Coulomb
    Cq = rows["coulomb"]
    d12, d18, d24 = (Cq[f"coulomb_same_d{d}"] for d in (12, 18, 24))
    Kcert = FIXED["I1"] * (d12["I1"] - d24["I1"]) / (1 / 12 - 1 / 24)
    row({k: -(d12[k] - d24[k]) for k in names}, -0.25 * FIXED["I1"] * (d12["I1"] - d24["I1"]),
        "coulomb: E(12) - E(24) >= 1/4 of the certified repulsion")
    if include_coulomb_form:
        for (da, db, Ea, Eb) in ((12, 18, d12, d18), (12, 24, d12, d24), (18, 24, d18, d24)):
            rho = 1 / da - 1 / db
            v = row({k: (Ea[k] - Eb[k]) for k in names}, 0.1 * Kcert * rho, f"coulomb form d{da},{db}: Delta - K rho <= 0.1 Kcert rho")
            v[-1] = -rho
            v2 = row({k: -(Ea[k] - Eb[k]) for k in names}, 0.1 * Kcert * rho, f"coulomb form d{da},{db}: -(Delta - K rho) <= 0.1 Kcert rho")
            v2[-1] = rho
        v = row({}, -0.5 * Kcert, "coulomb strength K >= Kcert/2")
        v[-1] = -1.0
    # (vi) positivity on every channel (+ the local clock on both sheets)
    for ch, q in rows["channels"].items():
        row({k: -q[k] for k in names}, 0.0, f"positivity: {ch}")
    for kind in ("twist", "zigzag"):
        q = rows["sheets"][f"{kind}|nz256|w3"]["q_per_area"]
        row({k: -q[k] for k in names}, 0.0, f"positivity: local clock on the {kind} sheet")
    # (vii) UV boundedness samples: the quadratic regime by plane-wave averages (total derivatives
    #       average out), the quartic regime EXACTLY: the degree-4 elements are homogeneous in the
    #       jets, so sum chat_k s_k(A) >= 0 on random jets is the necessary condition on the
    #       highest-degree form (a plane wave of amplitude eps at wavenumber k has jet ~ eps k, and
    #       the quartic dominates once eps k is large; V4 stays bounded per cell)
    if include_uv and "uv" in rows:
        for i, rec in enumerate(rows["uv"]):
            row({k: -rec["avg"][k] for k in names}, 0.0, f"uv: plane wave m={rec['m']} eps={rec['eps']} #{i}")
        if "uv2" in rows:
            for i, rec in enumerate(rows["uv2"]):
                row({k: -rec["avg"][k] for k in names}, 0.0, f"uv2: crossed waves eps={rec['eps']} #{i}")
        rng = np.random.default_rng(2828)
        M0 = np.broadcast_to(B3.vac4(B3.base_cfg(s=-1.0, g=8.0)), (150, 4, 4)).copy()
        Aq = rng.normal(size=(4, 150, 4, 4)); Aq = 0.5 * (Aq + Aq.swapaxes(-1, -2)); Aq[0] = 0.0
        quart = [e for e in BASIS if e.degree == 4]
        vals = {e.name: e.lag(Aq, M0) for e in quart}
        for i in range(150):
            row({e.name: -float(vals[e.name][i]) for e in quart}, 0.0, f"uv quartic form on random jet #{i}")
    for coefs, rhs, label in extra:
        row(coefs, rhs, label)
    A, b = np.array(A), np.array(b)
    # conditioning: per-row relative noise floor
    for i in range(A.shape[0]):
        mx = np.max(np.abs(A[i, :len(names)]))
        if mx > 0:
            A[i, :len(names)][np.abs(A[i, :len(names)]) < 1e-7 * mx] = 0.0
    return A, b, lab, Kcert, a


def solve(A, b, lab, active, norm=None, method="highs-ds", objective=None):
    """active: list of basis names allowed nonzero (others fixed 0); I1 fixed at 4.
    Variables: x (chat for all names), K, and for the norm bound the split x = p - m."""
    from scipy.optimize import linprog
    names = NAMES
    nv = len(names) + 1
    bounds = []
    for k in names:
        if k in FIXED:
            bounds.append((FIXED[k], FIXED[k]))
        elif k in active:
            bounds.append((None, None))
        else:
            bounds.append((0.0, 0.0))
    bounds.append((None, None))            # K
    A_ub, b_ub = A.copy(), b.copy()
    if norm is not None:
        # sum |x_k| <= norm over the active free variables via auxiliary t_k >= |x_k|
        nfree = [i for i, k in enumerate(names) if k in active and k not in FIXED]
        ext = np.zeros((A_ub.shape[0], len(nfree)))
        A_ub = np.hstack([A_ub, ext])
        rows_new, b_new = [], []
        for j, i in enumerate(nfree):
            r1 = np.zeros(nv + len(nfree)); r1[i] = 1.0; r1[nv + j] = -1.0
            r2 = np.zeros(nv + len(nfree)); r2[i] = -1.0; r2[nv + j] = -1.0
            rows_new += [r1, r2]; b_new += [0.0, 0.0]
        rn = np.zeros(nv + len(nfree)); rn[nv:] = 1.0
        rows_new.append(rn); b_new.append(norm)
        A_ub = np.vstack([A_ub] + [np.array(rows_new)])
        b_ub = np.concatenate([b_ub, np.array(b_new)])
        bounds += [(0.0, None)] * len(nfree)
    c = np.zeros(A_ub.shape[1])
    if objective == "min_norm":
        c[nv:] = 1.0                          # sum of the |x_k| auxiliaries
    elif objective is not None:
        for k, v in objective.items():
            c[names.index(k)] = v
    res = linprog(c, A_ub=A_ub, b_ub=b_ub, bounds=bounds, method=method)
    out = {"status": int(res.status), "message": str(res.message), "success": bool(res.success)}
    if res.success:
        x = res.x[:len(names)]
        out["x"] = {k: float(x[i]) for i, k in enumerate(names) if abs(x[i]) > 1e-12 or k in FIXED}
        out["K"] = float(res.x[len(names)])
        slack = b_ub[:len(lab)] - A_ub[:len(lab), :] @ res.x
        out["min_slack"] = float(np.min(slack))
        out["norm_l1"] = float(np.sum(np.abs(x[[i for i, k in enumerate(names) if k not in FIXED]])))
        out["tight_rows"] = [lab[i] for i in np.where(slack < 1e-9 * max(1.0, np.max(np.abs(b_ub))))[0]]
    return out


def farkas(A, b, lab, active, B=1e4, method="highs-ds"):
    """the explicit dual: find y >= 0 with (A^T y)_j = 0 on the FREE active variables, and the
    fixed variables folded into b; certificate valid for sum |x| <= B: beta + B ||r||_inf < 0,
    verified in exact rational arithmetic on the rounded rows."""
    from scipy.optimize import linprog
    names = NAMES
    nv = len(names) + 1
    free = [i for i, k in enumerate(names) if k in active and k not in FIXED] + [nv - 1]
    fixed_vec = np.zeros(nv)
    for k, v in FIXED.items():
        fixed_vec[names.index(k)] = v
    bb = b - A @ fixed_vec                                # rows: A_free x_free <= bb
    Af = A[:, free]
    m = A.shape[0]
    # minimize bb^T y  s.t.  Af^T y = 0, y >= 0, sum y = 1  (infeasible primal <=> optimum < 0)
    A_eq = np.vstack([Af.T, np.ones((1, m))])
    b_eq = np.concatenate([np.zeros(Af.shape[1]), [1.0]])
    res = linprog(bb, A_eq=A_eq, b_eq=b_eq, bounds=[(0.0, None)] * m, method=method)
    out = {"dual_status": int(res.status), "dual_success": bool(res.success)}
    if not res.success:
        return out
    y = res.x
    out["beta_float"] = float(bb @ y)
    out["r_inf_float"] = float(np.max(np.abs(Af.T @ y)))
    # exact check on rounded data
    def fr(x, den=10 ** 12):
        return Fraction(x).limit_denominator(den)
    yq = [fr(v) for v in y]
    beta = sum(yq[i] * fr(bb[i]) for i in range(m))
    r = [sum(yq[i] * fr(Af[i, j]) for i in range(m)) for j in range(Af.shape[1])]
    rinf = max(abs(v) for v in r)
    margin = beta + Fraction(int(B)) * rinf
    out["exact"] = {"beta": float(beta), "r_inf": float(rinf), "B": B, "beta_plus_B_rinf": float(margin),
                    "certificate_valid": bool(margin < 0)}
    out["y_support"] = [(lab[i], float(y[i])) for i in np.argsort(-y)[:12] if y[i] > 1e-9]
    return out


def stage_refine(max_rounds=60, tol=-1e-8):
    """cutting planes: solve the min-norm LP, check the exact UV quadratic form of the solution in
    every direction, add the violating eigen-directions as rows, repeat; ends CONE_FEASIBLE (the
    min-norm vertex has a PSD quadratic form in all seven directions) or CLASS_INFEASIBLE."""
    rows = json.load(open(ROWS_JSON))
    if "uvq" not in rows:
        rows = stage_uvq()
    A0, _, _, _, _ = build_rows(rows)
    rowmax = np.max(np.abs(A0[:, :len(NAMES)]), axis=1)
    rel = np.abs(A0[:, :len(NAMES)]) / np.maximum(rowmax, 1e-300)[:, None]
    null_cols = [k for i, k in enumerate(NAMES) if np.max(rel[:, i]) < 1e-6 and k not in FIXED]
    active = [k for k in NAMES if k not in null_cols]
    extra, history = [], []
    verdict = None
    for rnd in range(max_rounds):
        A, b, lab, Kcert, a = build_rows(rows, extra=tuple(extra))
        sol = solve(A, b, lab, active, norm=1e6, method="highs-ds", objective="min_norm")
        if not sol["success"]:
            sol_ipm = solve(A, b, lab, active, norm=1e6, method="highs-ipm", objective="min_norm")
            cert = farkas(A, b, lab, active, B=1e4)
            history.append({"round": rnd, "status": "INFEASIBLE", "ipm_agrees": not sol_ipm["success"], "farkas": cert, "n_extra_rows": len(extra)})
            verdict = "CLASS_INFEASIBLE" + (" (exact certificate)" if cert.get("exact", {}).get("certificate_valid") else " (solver only)")
            log(f"  refine round {rnd}: INFEASIBLE with {len(extra)} cutting planes; certificate {cert.get('exact', {}).get('certificate_valid')}")
            break
        x = sol["x"]
        me = uv_min_eig(rows, x)
        worst = min(v[0] for v in me.values())
        history.append({"round": rnd, "status": "feasible", "norm_l1": sol["norm_l1"], "x": x, "uv_min_eig_per_direction": {d: v[0] for d, v in me.items()}, "n_extra_rows": len(extra)})
        log(f"  refine round {rnd}: feasible, norm {sol['norm_l1']:.3f}, worst UV eigenvalue {worst:.3e}; x = { {k: round(v, 3) for k, v in x.items()} }")
        if worst >= tol:
            verdict = ("CONE_FEASIBLE (min-norm vertex UV-bounded in all seven directions)" if sol["norm_l1"] <= 100
                       else f"MARGINAL: feasible only at sum|chat| = {sol['norm_l1']:.0f} (outside any O(1) parsimony bound; the linear-response rows are not valid there), UV form PSD to {worst:.1e}")
            break
        for d, (w, v) in me.items():
            if w < tol:
                v = np.array(v)
                coefs = {k: -float(v @ np.array(Qk) @ v) for k, Qk in rows["uvq"][d].items()}
                extra.append((coefs, 0.0, f"uv cut dir {d} round {rnd}"))
    else:
        verdict = "UNDECIDED (cutting planes did not converge)"
    out = {"verdict": verdict, "history": history, "null_columns_removed": null_cols, "n_cutting_planes": len(extra)}
    # the same loop on the sub-bases that were feasible, and on the no-Coulomb-form variant
    out["variants"] = {}
    for vname, kw, sub, nb in (("two-derivative only", {}, SUBSETS["two-derivative only (K_T, K_lambda, R's, K_P_h, T1..T4)"], 1e6),
                               ("two-derivative + quartics (no F-built four-derivative terms)", {}, SUBSETS["two-derivative only (K_T, K_lambda, R's, K_P_h, T1..T4)"] + ["Q_I1sq", "C6a", "C6b"], 1e6),
                               ("without the Coulomb-form row", {"include_coulomb_form": False}, None, 1e6),
                               ("without the zigzag row", {"include_zigzag": False}, None, 1e6),
                               ("without the UV rows (cutting planes only)", {"include_uv": False}, None, 1e6),
                               ("B+ with sum|chat| <= 30", {}, None, 30.0),
                               ("B+ with sum|chat| <= 100", {}, None, 100.0),
                               ("B+ with sum|chat| <= 300", {}, None, 300.0),
                               ("B+ with sum|chat| <= 1000", {}, None, 1000.0)):
        act = [k for k in (sub or NAMES) if k not in null_cols]
        ex, hist = [], []
        vd = "UNDECIDED"
        for rnd in range(max_rounds):
            A, b, lab, Kcert, a = build_rows(rows, extra=tuple(ex), **kw)
            sol = solve(A, b, lab, act, norm=nb, method="highs-ds", objective="min_norm")
            if not sol["success"]:
                cert = farkas(A, b, lab, act, B=(1e4 if nb >= 1e6 else nb)) if nb >= 1e6 else {"note": "norm-bounded: infeasible within the bound"}
                vd = "CLASS_INFEASIBLE" + (" (exact certificate)" if cert.get("exact", {}).get("certificate_valid") else "")
                hist.append({"round": rnd, "status": "INFEASIBLE", "farkas": cert}); break
            me = uv_min_eig(rows, sol["x"])
            worst = min(v[0] for v in me.values())
            hist.append({"round": rnd, "norm_l1": sol["norm_l1"], "worst_uv_eig": worst, "x": sol["x"]})
            if worst >= tol:
                vd = "CONE_FEASIBLE"; break
            for d, (w, v) in me.items():
                if w < tol:
                    v = np.array(v)
                    ex.append(({k: -float(v @ np.array(Qk) @ v) for k, Qk in rows["uvq"][d].items()}, 0.0, f"cut {d} r{rnd}"))
        out["variants"][vname] = {"verdict": vd, "history": hist, "n_cutting_planes": len(ex)}
        log(f"  variant {vname}: {vd} after {len(ex)} cuts")
    json.dump(out, open(os.path.join(DATA, "m5_32_r14_a_refine.json"), "w"), indent=1, default=float)
    log(f"refine done: {verdict}")
    return out


SUBSETS = {
    "B+ (all)": NAMES,
    "two-derivative only (K_T, K_lambda, R's, K_P_h, T1..T4)": ["I1", "K_T", "K_lambda", "R_etaMeta", "R_hcov", "K_P_h", "T1", "T2", "T3", "T4"],
    "four-derivative F-built only (R1's thirteen + Q_I1sq)": ["I1", "I2", "I3", "I4", "I5", "I6", "I1_h", "J1", "J2", "Pgrad", "Q_I1sq"] + [n for n in NAMES if n.startswith("E")],
    "non-F quartics only (C6a, C6b)": ["I1", "C6a", "C6b"],
    "Q_F direction only (T1, T3)": ["I1", "T1", "T3"],
    "R1's thirteen + entrants, no T-forms, no quartics": ["I1", "I2", "I3", "I4", "I5", "I6", "I1_h", "J1", "J2", "Pgrad", "K_T", "K_lambda", "R_etaMeta", "R_hcov", "K_P_h"] + [n for n in NAMES if n.startswith("E")],
}


def stage_lp():
    rows = json.load(open(ROWS_JSON))
    A0, _, _, _, _ = build_rows(rows)
    # a column is null if it vanishes on every row, or NEAR-null if on every row its entry is
    # below 1e-6 of that row's largest entry (roundoff-level content cannot carry feasibility)
    rowmax = np.max(np.abs(A0[:, :len(NAMES)]), axis=1)
    rel = np.abs(A0[:, :len(NAMES)]) / np.maximum(rowmax, 1e-300)[:, None]
    null_cols = [k for i, k in enumerate(NAMES) if np.max(rel[:, i]) < 1e-6 and k not in FIXED]
    colmax = {k: float(np.max(rel[:, i])) for i, k in enumerate(NAMES)}
    out = {"basis": NAMES, "fixed": FIXED, "null_columns_removed": null_cols, "column_max_relative_entry": colmax, "variants": {}}
    log(f"null / near-null columns (removed): {null_cols}")
    for sub in SUBSETS:
        SUBSETS[sub] = [k for k in SUBSETS[sub] if k not in null_cols]
    for variant, kw in (("frozen rows + Coulomb-form + zigzag + UV (the 09-05 extended packet)", {}),
                        ("WITHOUT the UV/G5 row (the packet as first frozen)", {"include_uv": False}),
                        ("WITHOUT the Coulomb-form row (the author's prediction: feasible via Q_F)", {"include_coulomb_form": False}),
                        ("WITHOUT the Coulomb-form AND the UV rows", {"include_coulomb_form": False, "include_uv": False}),
                        ("WITHOUT the zigzag row", {"include_zigzag": False}),
                        ("WITHOUT the relaxed-pair row", {"include_relaxed_pair": False})):
        A, b, lab, Kcert, a = build_rows(rows, **kw)
        vres = {"n_rows": len(lab), "rows": lab, "Kcert": Kcert, "tail_a": a, "subsets": {}}
        for sub, active in SUBSETS.items():
            sres = {}
            for norm in (None, 1.0, 10.0, 100.0, 1e4):
                key = "cone" if norm is None else f"norm<={norm:g}"
                r_ds = solve(A, b, lab, active, norm=norm, method="highs-ds")
                r_ipm = solve(A, b, lab, active, norm=norm, method="highs-ipm")
                sres[key] = {"highs-ds": r_ds, "highs-ipm": r_ipm,
                             "agree": bool(r_ds["success"] == r_ipm["success"])}
                if not r_ds["success"] and norm is None:
                    sres[key]["farkas"] = farkas(A, b, lab, active)
                elif r_ds["success"] and norm is None:
                    # the primary read: the MINIMUM-NORM feasible action (min sum |chat| over the
                    # rows), plus two extreme vertices at a norm bound of 1e4
                    verts = {"min_norm": solve(A, b, lab, active, norm=1e6, method="highs-ds", objective="min_norm")}
                    for onm, obj in (("max_attraction", {k: (rows["pairs"]["ansatz_same_d16"][k] - rows["pairs"]["ansatz_same_d8"][k]) * -1 for k in NAMES}),
                                     ("max_twist_margin", {k: -rows["sheets"]["twist|nz256|w1.5"]["S_per_area"][k] for k in NAMES})):
                        verts[onm] = solve(A, b, lab, active, norm=1e4, method="highs-ds", objective=obj)
                    sres[key]["vertices"] = verts
            vres["subsets"][sub] = sres
            st = {k: ("FEASIBLE" if v["highs-ds"]["success"] else "INFEASIBLE") for k, v in sres.items()}
            log(f"  [{variant[:38]:38s}] {sub[:52]:52s} {st}")
        out["variants"][variant] = vres
    # verdict on the primary variant
    prim = out["variants"]["frozen rows + Coulomb-form + zigzag + UV (the 09-05 extended packet)"]["subsets"]["B+ (all)"]
    cone = prim["cone"]
    if cone["highs-ds"]["success"]:
        out["verdict"] = "CONE_FEASIBLE"
    else:
        cert = cone.get("farkas", {}).get("exact", {}).get("certificate_valid", False)
        out["verdict"] = "CLASS_INFEASIBLE" + (" (exact certificate)" if cert else " (solver only; certificate not validated)")
    out["wall_s"] = round(time.time() - T0, 1)
    json.dump(out, open(LP_JSON, "w"), indent=1, default=float)
    log(f"LP done: {out['verdict']} -> {LP_JSON}")
    return out


if __name__ == "__main__":
    stage = ARGV[1] if len(ARGV) > 1 else "all"
    print("STAGE", repr(stage), ARGV, flush=True)
    if stage in ("measure", "all"):
        stage_measure()
    if stage in ("uv", "all"):
        stage_uv()
    if stage in ("lp", "all"):
        stage_lp()
    if stage == "uvq":
        stage_uvq()
    if stage == "uv2":
        stage_uv2()
    if stage in ("refine", "all"):
        stage_refine()
