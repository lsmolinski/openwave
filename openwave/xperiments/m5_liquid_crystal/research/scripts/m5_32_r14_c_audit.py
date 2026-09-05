"""M5.32 R14-C adversarial audit: the scalar-channel Newton route (R_G on the dressed
pair, K_lambda exchange).  Independent implementation on the SAVED fields; the
producer's script (m5_32_r14_c_newton.py) and its json are never opened.

EQUATIONS FIRST (the brief's definitions, re-implemented here)
--------------------------------------------------------------
Jets:      A_i = d_i M, i = 1..3, the certified sym stencil = average of the fwd and bwd
           one-sided branches of the DENSITY (own one-sided differences, checked against
           the certified B3.d1).
R_G:       density = sum_{i,j} G_cd [ (A_i)[j,c] (A_j)[i,d] - (A_i)[i,c] (A_j)[j,d] ]
           (raw entries; the derivative index of one jet against a raw internal row of
           the other by delta; G covariant (0,2)); E_RG = h^3 sum_br wt sum_cells.
           Own contraction = explicit loops over (i, j) with a per-cell 4x4 G contraction.
E_c:       E_lambda[M] + c_R E_RG[M], E_lambda = m5_32_r2_b_bounded.energy_grad(M, cfg, lam)[0]
           (at lambda = 0 cross-checked against the certified B3.e_parts sum).
E_int(d):  E(pair at d) - 2 E(single), all three under the same (lambda, c_R).
K_lambda:  dl_a(x) = lambda_a(x) - lambda_a^vac, lambda sorted DESCENDING per cell
           (vacuum (1, 0.3, 0, -8)); static E = (1/2) h^3 sum_i sum_a sum_cells
           (nearest-neighbour difference / h)^2 (the fwd and bwd branches carry the same
           set of squared differences, so sym == either branch).
Exchange:  L = (1/2)(d phi)^2 - (1/2) m_s^2 phi^2 + J phi, phi = dl_a:
           J = (-nabla^2 + m_s^2) dl_a; E_int(d) = -sum_x J(x) dl_a(x - d e_z) h^3.
           Two independent kernels here: (K1) real-space 7-point periodic Laplacian +
           np.roll shift; (K2) Parseval: -sum_k |dl_k|^2 (k^2 + m_s^2) cos(k_z d) h^3 / n^3.
"""
from __future__ import annotations
import sys
ARGV = list(sys.argv)
import hashlib
import importlib.util
import json
import os
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.dirname(HERE)
DATA = os.path.join(RES, "data")
CK = os.path.join(RES, "checkpoints")
T0 = time.time()
ETA = np.diag([-1.0, 1.0, 1.0, 1.0])


def log(m):
    print(f"[{time.time() - T0:7.1f}s] {m}", flush=True)


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


B3 = _load("m5_21_3_a_4d", "m5_21_3_a_4d.py")
CM = _load("m5_32_r13w_common", "m5_32_r13w_common.py")
R2B = _load("m5_32_r2_b_bounded", "m5_32_r2_b_bounded.py")
R3I = _load("m5_32_r3_i_ansatz", "m5_32_r3_i_ansatz.py")
R14 = _load("m5_32_r14_terms", "m5_32_r14_terms.py")

CFG = CM.cfg_of(32, 48.0)               # g = 8: the hedgehogs (m5_32_r10), the embedded 3x3, the ansatz
CFG32 = CM.cfg_of(32, 48.0, g=32.0)     # g = 32: the R3 relaxed pairs + the R14-C end fields (M_00 = 32 at the
                                        # corner; the brief's "cfg_of(32, 48.0)" (g = 8) gives V4 = 8.7e13 on them)
H = CFG["h"]
H3 = H ** 3
VAC = B3.vac4(CFG)                      # diag(8, 1, 0.3, 0)
VAC32 = B3.vac4(CFG32)                  # diag(32, 1, 0.3, 0)
LAM_VAC = np.array([1.0, 0.3, 0.0, -8.0])   # sorted descending spectrum of N = M eta

FORBIDDEN = ("m5_32_r14_c_newton.py", "m5_32_r14_c_newton.json")
RESULTS = {"argv": ARGV, "forbidden_never_opened": list(FORBIDDEN)}


def sha(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()[:16]


# ============================================================ own stencils
def d1_own(f, ax, h, br):
    """one-sided difference, fwd placed on cells 0..n-2, bwd on 1..n-1 (the certified
    placement), zero elsewhere."""
    out = np.zeros_like(f)
    df = np.diff(f, axis=ax) / h
    sl = [slice(None)] * f.ndim
    sl[ax] = slice(0, -1) if br == "fwd" else slice(1, None)
    out[tuple(sl)] = df
    return out


def jets_own(M, br):
    return [d1_own(M, ax, H, br) for ax in range(3)]


# ============================================================ own R_G
def G_tensor(kind, M):
    if kind == "eta":
        return np.broadcast_to(ETA, M.shape)
    if kind == "etaMeta":
        return ETA @ M @ ETA
    if kind == "neg_etaMeta":
        return -(ETA @ M @ ETA)
    if kind == "flipM":                                 # M-dependence flipped about the vacuum
        return ETA @ (2.0 * VAC - M) @ ETA
    if kind == "hcov":
        u = R2B.tl_eig(M)[0]
        w = u @ ETA
        return ETA + 2.0 * w[..., :, None] * w[..., None, :]
    raise ValueError(kind)


def rg_density_own(A, G, plus_mutant=False):
    """sum_{i,j} G_cd [A_i[j,c] A_j[i,d] - A_i[i,c] A_j[j,d]], A = [A_1, A_2, A_3]
    with raw entries; internal row index j <-> derivative index j (delta).
    plus_mutant: the relative sign of the two terms flipped (a wrong-contraction
    control).  NOTE: contracting the column instead of the row is NOT a mutation
    (M symmetric => A_i symmetric), checked and dropped."""
    dens = np.zeros(A[0].shape[:-2])
    div = np.zeros(A[0].shape[:-2] + (4,))
    for i in range(3):
        div += A[i][..., 1 + i, :]
    for i in range(3):
        for j in range(3):
            ai = A[i][..., 1 + j, :]                 # (A_i)^{j c}
            aj = A[j][..., 1 + i, :]                 # (A_j)^{i d}
            dens += np.einsum("...cd,...c,...d->...", G, ai, aj)
    sgn = 1.0 if plus_mutant else -1.0
    dens += sgn * np.einsum("...cd,...c,...d->...", G, div, div)
    return dens


def e_rg_own(M, kind, plus_mutant=False, h=None):
    G = G_tensor(kind, M)
    hh = h or H
    tot = 0.0
    for br in ("fwd", "bwd"):
        A = [d1_own(M, ax, hh, br) for ax in range(3)]
        tot += 0.5 * np.sum(rg_density_own(A, G, plus_mutant))
    return float(hh ** 3 * tot)


def e_lambda(M, lam, cfg=None):
    e, _, info = R2B.energy_grad(M, cfg or CFG32, lam)
    if not info.get("ok", True):
        raise RuntimeError(f"energy_grad locus failure {info}")
    return float(e)


# ============================================================ own K_lambda
def sorted_spec_own(M):
    lam = np.linalg.eigvals(M @ ETA)
    if np.max(np.abs(lam.imag)) > 1e-8:
        raise ValueError("complex spectrum")
    return np.sort(lam.real, axis=-1)[..., ::-1]


def klam_static_own(M):
    lam = sorted_spec_own(M)
    tot = 0.0
    for ax in range(3):
        d = np.diff(lam, axis=ax) / H
        tot += 0.5 * np.sum(d * d)
    return float(H3 * tot)


# ============================================================ loaders
def load_r3(lam, name):
    p = os.path.join(DATA, "m5_32_r3_ii", f"lam{lam}_dr1_{name}_n32.npz")
    return np.load(p)["M"].astype(np.float64), p


def load_end(lam, cR, name):
    tag = {1.0: "+1", -1.0: "-1", 0.3: "+0.3", -0.3: "-0.3"}[cR]
    p = os.path.join(CK, "m5_32_r14", "c_relax", f"lam{lam}_cR{tag}_{name}.npy")
    return np.load(p).astype(np.float64), p


def load_hedgehog(it):
    p = os.path.join(CK, "m5_32_r10", f"relax_g8_n32_L48_it{it}.npy")
    return np.load(p).astype(np.float64), p


# ============================================================ stage: sanity
def stage_sanity():
    out = {}
    M, _ = load_r3(0, "same_d10")
    # own stencil == certified stencil
    for br in ("fwd", "bwd"):
        for ax in range(3):
            assert np.allclose(d1_own(M, ax, H, br), B3.d1(M, ax, H, br))
    out["stencil_matches_certified"] = True
    # own R_G contraction against the registry read (allowed cross-check)
    for kind in ("eta", "etaMeta", "hcov"):
        mine = e_rg_own(M, kind)
        reg = R14.static_energy(f"R_{kind}", M, CFG32)
        out[f"E_RG_{kind}_own"] = mine
        out[f"E_RG_{kind}_registry"] = reg
        out[f"E_RG_{kind}_reldiff"] = abs(mine - reg) / max(abs(reg), 1e-12)
        log(f"sanity R_{kind}: own {mine:.6f} registry {reg:.6f}")
    # E_lambda at lam = 0 == certified e_parts sum
    e0 = e_lambda(M, 0.0)
    eu, ev = B3.e_parts(M, CFG32)
    out["corner_cell_M00"] = float(M[0, 0, 0, 0, 0])
    out["V4_with_g8_cfg"] = float(B3.e_parts(M, CFG)[1])
    out["V4_with_g32_cfg"] = float(ev)
    out["E_lambda0_vs_e_parts"] = [e0, eu + ev]
    log(f"sanity E_lambda(0) {e0:.6f} vs e_parts {eu + ev:.6f}")
    # K_lambda own vs registry
    kl = klam_static_own(M)
    kr = R14.static_energy("K_lambda", M, CFG32)
    out["K_lambda_own_vs_registry"] = [kl, kr]
    log(f"sanity K_lambda own {kl:.6f} registry {kr:.6f}")
    # R_eta is a total derivative: on a smooth periodic random field it must vanish
    rng = np.random.default_rng(7)
    n = 24
    k = np.fft.fftfreq(n, d=1.0 / n)
    KX, KY, KZ = np.meshgrid(k, k, k, indexing="ij")
    filt = np.exp(-(KX ** 2 + KY ** 2 + KZ ** 2) / 8.0)
    F = np.zeros((n, n, n, 4, 4))
    for a in range(4):
        for b in range(a, 4):
            f = np.fft.ifftn(np.fft.fftn(rng.standard_normal((n, n, n))) * filt).real
            F[..., a, b] = f
            F[..., b, a] = f
    Mp = VAC + 0.5 * F
    hp = 2 * np.pi / n
    # periodic spectral jets (own) for the total-derivative test
    Ap = []
    for ax in range(3):
        kk = [KX, KY, KZ][ax]
        Ap.append(np.fft.ifftn(1j * kk[..., None, None] * np.fft.fftn(Mp, axes=(0, 1, 2)),
                               axes=(0, 1, 2)).real)
    r_eta = float(np.sum(rg_density_own(Ap, np.broadcast_to(ETA, Mp.shape))))
    r_eme = float(np.sum(rg_density_own(Ap, ETA @ Mp @ ETA)))
    scale = float(np.sum(np.abs(rg_density_own(Ap, np.broadcast_to(ETA, Mp.shape)))))
    out["periodic_total_derivative"] = {"R_eta_sum": r_eta, "R_eta_abs_sum": scale,
                                        "R_etaMeta_sum": r_eme}
    log(f"sanity periodic: sum R_eta {r_eta:.3e} (abs {scale:.3e}), sum R_etaMeta {r_eme:.3e}")
    RESULTS["sanity"] = out


# ============================================================ C1 / C2
def stage_c12():
    out = {}
    for lam in (0, 1):
        rec = {}
        fields = {nm: load_r3(lam, nm)[0] for nm in ("same_d10", "same_d18", "single_d0")}
        rec["files"] = {nm: sha(load_r3(lam, nm)[1]) for nm in fields}
        el = {nm: e_lambda(fields[nm], float(lam)) for nm in fields}
        rec["E_lambda"] = el
        rec["slope_certified"] = el["same_d18"] - el["same_d10"]
        rec["E_RG"] = {}
        rec["slope_RG"] = {}
        for kind in ("etaMeta", "eta", "hcov", "neg_etaMeta", "flipM"):
            erg = {nm: e_rg_own(fields[nm], kind) for nm in fields}
            rec["E_RG"][kind] = erg
            rec["slope_RG"][kind] = erg["same_d18"] - erg["same_d10"]
        erg_t = {nm: e_rg_own(fields[nm], "etaMeta", plus_mutant=True) for nm in fields}
        rec["E_RG"]["etaMeta_plus_sign_mutant"] = erg_t
        rec["slope_RG"]["etaMeta_plus_sign_mutant"] = erg_t["same_d18"] - erg_t["same_d10"]
        rec["dEint"] = {}
        for kind in ("etaMeta", "eta", "hcov", "neg_etaMeta", "flipM", "etaMeta_plus_sign_mutant"):
            rec["dEint"][kind] = {}
            for cR in (1.0, -1.0, 0.3, -0.3, 0.0):
                eint = {nm: el[nm] + cR * rec["E_RG"][kind][nm] for nm in fields}
                e10 = eint["same_d10"] - 2 * eint["single_d0"]
                e18 = eint["same_d18"] - 2 * eint["single_d0"]
                rec["dEint"][kind][str(cR)] = {"Eint10": e10, "Eint18": e18, "d18_minus_d10": e18 - e10}
        log(f"C{1 + lam}: lam {lam} certified slope {rec['slope_certified']:.2f} "
            f"RG slope etaMeta {rec['slope_RG']['etaMeta']:.2f} eta {rec['slope_RG']['eta']:.2f} "
            f"hcov {rec['slope_RG']['hcov']:.2f} flipM {rec['slope_RG']['flipM']:.2f}")
        for cR in ("1.0", "-1.0"):
            log(f"   cR {cR}: dEint(18-10) etaMeta {rec['dEint']['etaMeta'][cR]['d18_minus_d10']:.2f}")
        out[f"lam{lam}"] = rec
    RESULTS["c1c2"] = out


# ============================================================ C3
def stage_c3():
    out = {}
    for lam in (0, 1):
        seeds = {nm: load_r3(lam, nm)[0] for nm in ("same_d10", "same_d18", "single_d0")}
        for cR in (1.0, -1.0, 0.3, -0.3):
            rec = {}
            ends = {}
            for nm in seeds:
                Me, p = load_end(lam, cR, nm)
                ends[nm] = Me
                dM = Me - seeds[nm]
                rec[nm] = {
                    "file": os.path.basename(p), "sha": sha(p),
                    "rel_frob": float(np.linalg.norm(dM) / np.linalg.norm(seeds[nm] - VAC32)),
                    "max_abs_dM": float(np.max(np.abs(dM))),
                    "max_abs_M0i_seed": float(np.max(np.abs(seeds[nm][..., 0, 1:]))),
                    "max_abs_M0i_end": float(np.max(np.abs(Me[..., 0, 1:]))),
                    "E_RG_seed": e_rg_own(seeds[nm], "etaMeta"),
                    "E_RG_end": e_rg_own(Me, "etaMeta"),
                    "E_lambda_seed": e_lambda(seeds[nm], float(lam)),
                    "E_lambda_end": e_lambda(Me, float(lam)),
                }
                rec[nm]["E_c_seed"] = rec[nm]["E_lambda_seed"] + cR * rec[nm]["E_RG_seed"]
                rec[nm]["E_c_end"] = rec[nm]["E_lambda_end"] + cR * rec[nm]["E_RG_end"]
                rec[nm]["E_c_drop"] = rec[nm]["E_c_seed"] - rec[nm]["E_c_end"]
                # the producer's descent record (allowed): its own bookkeeping of the same numbers
                jp = p.replace(".npy", ".json")
                J = json.load(open(jp))
                rec[nm]["producer_json"] = {k: J.get(k) for k in
                                            ("E_seed_total", "E_end_total", "E_RG_seed", "E_RG_end")}
                rec[nm]["producer_json"]["accepted"] = J.get("descent", {}).get("accepted")
                rec[nm]["producer_json"]["verdict"] = J.get("descent", {}).get("verdict")
                rec[nm]["producer_json"]["E_drop"] = J.get("descent", {}).get("E_drop")
                rec[nm]["producer_json"]["last_quarter_dE"] = J.get("descent", {}).get("last_quarter_dE")
                rec[nm]["json_vs_own_seed_total"] = abs(J["E_seed_total"] - rec[nm]["E_c_seed"])
                rec[nm]["json_vs_own_end_total"] = abs(J["E_end_total"] - rec[nm]["E_c_end"])
                rec[nm]["json_vs_own_E_RG_seed"] = abs(J["E_RG_seed"] - rec[nm]["E_RG_seed"])
                rec[nm]["json_vs_own_E_RG_end"] = abs(J["E_RG_end"] - rec[nm]["E_RG_end"])
                # the effective R_G coefficient the producer's totals imply, two independent solves:
                # (i) from E_seed_total = E_lambda + x E_RG on the seed; (ii) from the logged E_drop
                rec[nm]["x_eff_from_seed_total"] = (J["E_seed_total"] - rec[nm]["E_lambda_seed"]) / rec[nm]["E_RG_seed"]
                dEl = rec[nm]["E_lambda_seed"] - rec[nm]["E_lambda_end"]
                dEr = rec[nm]["E_RG_seed"] - rec[nm]["E_RG_end"]
                ed = J.get("descent", {}).get("E_drop")
                rec[nm]["x_eff_from_E_drop"] = ((ed - dEl) / dEr) if (ed is not None and abs(dEr) > 1e-9) else None
                rec[nm]["dE_lambda_drop"] = dEl
                rec[nm]["dE_RG_drop"] = dEr
            e_end = {nm: rec[nm]["E_c_end"] for nm in seeds}
            e_seed = {nm: rec[nm]["E_c_seed"] for nm in seeds}
            rec["dEint_end"] = (e_end["same_d18"] - e_end["same_d10"])
            rec["dEint_seed"] = (e_seed["same_d18"] - e_seed["same_d10"])
            rec["Eint_end"] = {d: e_end[f"same_d{d}"] - 2 * e_end["single_d0"] for d in (10, 18)}
            rec["Eint_seed"] = {d: e_seed[f"same_d{d}"] - 2 * e_seed["single_d0"] for d in (10, 18)}
            # per-field E_c drops: does the heal move d10 and d18 by the same amount?
            rec["drop_d10_minus_d18"] = rec["same_d10"]["E_c_drop"] - rec["same_d18"]["E_c_drop"]
            rec["x_eff_seed_total"] = {nm: rec[nm]["x_eff_from_seed_total"] for nm in seeds}
            rec["x_eff_E_drop"] = {nm: rec[nm]["x_eff_from_E_drop"] for nm in seeds}
            log(f"C3 lam {lam} cR {cR:+}: x_eff(seed_total) {np.round(list(rec['x_eff_seed_total'].values()), 3).tolist()} "
                f"x_eff(E_drop) {[None if v is None else round(v, 3) for v in rec['x_eff_E_drop'].values()]}")
            log(f"C3 lam {lam} cR {cR:+}: dEint(18-10) seed {rec['dEint_seed']:.1f} end {rec['dEint_end']:.1f}; "
                f"drops d10 {rec['same_d10']['E_c_drop']:.1f} d18 {rec['same_d18']['E_c_drop']:.1f} "
                f"single {rec['single_d0']['E_c_drop']:.1f}; rel_frob d10 {rec['same_d10']['rel_frob']:.2e}")
            out[f"lam{lam}_cR{cR:+}"] = rec
    RESULTS["c3"] = out


# ============================================================ C4
def stage_c4():
    out = {}
    for d in (12, 18, 24):
        p = os.path.join(DATA, f"m5_21_4_lad_same_d{d}_n32_it400.npz")
        M3 = np.load(p)["M"].astype(np.float64)
        M4 = B3.embed34(M3, CFG)                       # M_00 = -sg = 8, M_0i = 0
        assert abs(M4[0, 0, 0, 0, 0] - 8.0) < 1e-12
        rec = {"sha": sha(p), "M00": float(M4[0, 0, 0, 0, 0])}
        for kind in ("etaMeta", "eta", "hcov"):
            rec[f"E_RG_{kind}"] = e_rg_own(M4, kind)
        rec["E_RG_etaMeta_plus_sign_mutant"] = e_rg_own(M4, "etaMeta", plus_mutant=True)
        eu, ev = B3.e_parts(M4, CFG)
        rec["E_u_4I1"] = float(eu)
        rec["E_V4"] = float(ev)
        out[f"d{d}"] = rec
        log(f"C4 d{d}: R_etaMeta {rec['E_RG_etaMeta']:.4f} R_eta {rec['E_RG_eta']:.4f} "
            f"R_hcov {rec['E_RG_hcov']:.4f} 4I1 {eu:.4f} V4 {ev:.4f}")
    out["variants"] = {}
    for suf in ("", "_it120"):
        v = {}
        for d in (12, 24):
            M3 = np.load(os.path.join(DATA, f"m5_21_4_lad_same_d{d}_n32{suf}.npz"))["M"].astype(np.float64)
            M4 = B3.embed34(M3, CFG)
            v[d] = {"etaMeta": e_rg_own(M4, "etaMeta"), "eta": e_rg_own(M4, "eta"), "4I1": float(B3.e_parts(M4, CFG)[0])}
        out["variants"][suf or "_plain"] = {k: v[12][k] - v[24][k] for k in v[12]}
        log(f"C4 variant {suf or '_plain'} E(12)-E(24): {out['variants'][suf or '_plain']}")
    out["diff_12_24"] = {k: out["d12"][k] - out["d24"][k] for k in out["d12"] if k.startswith("E_")}
    out["diff_12_18"] = {k: out["d12"][k] - out["d18"][k] for k in out["d12"] if k.startswith("E_")}
    log(f"C4 E(12)-E(24): {out['diff_12_24']}")
    RESULTS["c4"] = out


# ============================================================ C5
def fit_power(ds, es, p):
    """least squares E = a + b d^-p; returns a, b, R^2."""
    ds, es = np.asarray(ds, float), np.asarray(es, float)
    X = np.stack([np.ones_like(ds), ds ** (-p)], axis=1)
    coef, *_ = np.linalg.lstsq(X, es, rcond=None)
    pred = X @ coef
    ss = np.sum((es - es.mean()) ** 2)
    r2 = 1.0 - np.sum((es - pred) ** 2) / ss if ss > 0 else float("nan")
    return float(coef[0]), float(coef[1]), float(r2)


def fit_free_power(ds, es):
    """3-parameter fit E = a + b d^-p by a 1-d scan over p (nested linear lstsq)."""
    best = None
    for p in np.linspace(0.2, 8.0, 391):
        a, b, r2 = fit_power(ds, es, p)
        if best is None or r2 > best[3]:
            best = (float(p), a, b, r2)
    return best


def stage_c5():
    out = {}
    n, L = 32, 48.0
    for axis, vac, cfg in (("z", VAC, CFG), ("x", VAC, CFG), ("z_g32", VAC32, CFG32)):
        rec = {"E": {}, "g": float(cfg["g"])}
        ds = [8, 10, 12, 16, 20, 24]
        for d in ds:
            c1 = [0.0, 0.0, 0.0]
            c2 = [0.0, 0.0, 0.0]
            k = {"x": 0, "z": 2, "z_g32": 2}[axis]
            c1[k], c2[k] = d / 2.0, -d / 2.0
            M = R3I.ansatz(vac, [c1, c2], 0.1, 1.0, n, L)
            e = {"R_etaMeta": e_rg_own(M, "etaMeta"), "R_hcov": e_rg_own(M, "hcov"),
                 "R_eta": e_rg_own(M, "eta")}
            eu, ev = B3.e_parts(M, cfg)
            e["4I1"] = float(eu)
            e["K_P_h"] = float(R14.static_energy("K_P_h", M, cfg))
            e["K_lambda"] = klam_static_own(M)
            rec["E"][str(d)] = e
            log(f"C5 {axis} d{d}: " + " ".join(f"{k} {v:.4f}" for k, v in e.items()))
        # single at the origin, for the bookkeeping of E_int
        Ms = R3I.ansatz(vac, [[0.0, 0.0, 0.0]], 0.1, 1.0, n, L)
        rec["E_single"] = {"R_etaMeta": e_rg_own(Ms, "etaMeta"), "R_hcov": e_rg_own(Ms, "hcov"),
                           "4I1": float(B3.e_parts(Ms, cfg)[0]),
                           "K_P_h": float(R14.static_energy("K_P_h", Ms, cfg))}
        rec["fits"] = {}
        for key in ("R_etaMeta", "R_hcov", "4I1", "K_P_h"):
            E = {d: rec["E"][str(d)][key] for d in ds}
            f = {}
            # the producer's form: local exponent of E(d) - E(16) between 8 and 12
            y8, y12 = E[8] - E[16], E[12] - E[16]
            f["local_exp_8_12_ref16"] = (float(np.log(y8 / y12) / np.log(12.0 / 8.0))
                                         if y8 * y12 > 0 else None)
            f["sign_E8_E16"] = float(np.sign(y8))
            f["sign_E12_E16"] = float(np.sign(y12))
            # 1/d and 1/d^2 fits on E(d) - E(16), d in (8, 10, 12) plus the reference
            dd = [8, 10, 12, 16]
            yy = [E[d] - E[16] for d in dd]
            for p in (1.0, 2.0):
                a, b, r2 = fit_power(dd, yy, p)
                f[f"R2_p{p:g}_ref16_d8to16"] = r2
            # own: free-exponent 3-parameter fit on the full shell, no reference subtraction
            for shell in ((8, 10, 12, 16), (8, 10, 12, 16, 20), (8, 10, 12, 16, 20, 24), (10, 12, 16, 20, 24)):
                pbest, a, b, r2 = fit_free_power(list(shell), [E[d] for d in shell])
                f[f"free_fit_{'_'.join(map(str, shell))}"] = {"p": pbest, "a": a, "b": b, "R2": r2}
            # pairwise local exponents of E(d) - E(24) (the reference in the far corner)
            loc = {}
            for d1_, d2_ in ((8, 10), (10, 12), (12, 16), (16, 20)):
                y1, y2 = E[d1_] - E[24], E[d2_] - E[24]
                loc[f"{d1_}_{d2_}"] = (float(np.log(y1 / y2) / np.log(d2_ / d1_)) if y1 * y2 > 0 else None)
            f["local_exp_ref24"] = loc
            dd_ = [8, 10, 12, 16, 20, 24]
            f["slope_per_unit_d"] = {f"{a}_{b}": (E[b] - E[a]) / (b - a) for a, b in zip(dd_[:-1], dd_[1:])}
            f["E_minus_2single"] = {str(d): E[d] - 2 * rec["E_single"][key] for d in dd_}
            rec["fits"][key] = f
            log(f"C5 {axis} {key}: local(8,12|ref16) {f['local_exp_8_12_ref16']} "
                f"R2(1/d) {f['R2_p1_ref16_d8to16']:.3f} R2(1/d2) {f['R2_p2_ref16_d8to16']:.3f} "
                f"free p(8..16) {f['free_fit_8_10_12_16']['p']:.2f} free p(8..24) {f['free_fit_8_10_12_16_20_24']['p']:.2f}")
        out[axis] = rec
    RESULTS["c5"] = out


def stage_c5h():
    """h-refinement of the ansatz reads (n = 32, 48, 64 at L = 48, g = 8): are the R_G
    pair numbers converged in h, and how big is the total-derivative G = eta piece?"""
    out = {}
    L = 48.0
    for n in (32, 48, 64):
        h = L / n
        cfg = CM.cfg_of(n, L)
        rec = {}
        for tag, centers in (("single", [[0, 0, 0]]), ("d8", [[0, 0, 4.0], [0, 0, -4.0]]),
                             ("d16", [[0, 0, 8.0], [0, 0, -8.0]])):
            M = R3I.ansatz(VAC, centers, 0.1, 1.0, n, L)
            rec[tag] = {"etaMeta": e_rg_own(M, "etaMeta", h=h), "hcov": e_rg_own(M, "hcov", h=h),
                        "eta": e_rg_own(M, "eta", h=h), "4I1": float(B3.e_parts(M, cfg)[0])}
        rec["Eint"] = {d: {k: rec[d][k] - 2 * rec["single"][k] for k in rec["single"]} for d in ("d8", "d16")}
        rec["Eint8_minus_Eint16"] = {k: rec["Eint"]["d8"][k] - rec["Eint"]["d16"][k] for k in rec["single"]}
        log(f"C5h n{n}: single {rec['single']} | Eint8-Eint16 {rec['Eint8_minus_Eint16']}")
        out[f"n{n}"] = rec
    RESULTS["c5h"] = out


def stage_c3grad():
    """which R_G coefficient did each heal actually descend under?  On the end field,
    x* = argmin_x |grad E_lambda + x grad E_RG|^2 over the free cells; and the seed-to-end
    displacement's alignment with -(grad E_lambda + x grad E_RG)(seed)."""
    out = {}
    free = ~B3.pin_shell(32, H)
    for lam, cR, nm in ((0, 1.0, "same_d10"), (0, 1.0, "single_d0"), (0, -1.0, "same_d10"),
                        (0, 0.3, "same_d18"), (1, -1.0, "same_d10"), (1, -1.0, "same_d18"),
                        (1, 1.0, "same_d10"), (1, -0.3, "same_d10")):
        seed = load_r3(lam, nm)[0]
        end, p = load_end(lam, cR, nm)
        rec = {"label_cR": cR}
        for which, M in (("end", end), ("seed", seed)):
            gl = R2B.energy_grad(M, CFG32, float(lam))[1]
            gr = R14.rg_grad(M, CFG32, "etaMeta")
            gl = gl[free]
            gr = gr[free]
            xs = -float(np.sum(gl * gr) / np.sum(gr * gr))
            res = {x: float(np.max(np.abs(gl + x * gr))) for x in (cR, xs)}
            rec[which] = {"x_star": xs, "fmax_at_label": res[cR], "fmax_at_xstar": res[xs],
                          "fmax_ratio_label_over_xstar": res[cR] / res[xs]}
            if which == "seed":
                dlt = (end - seed)[free]
                cos = {}
                for x in sorted({cR, xs, 0.0, 0.7, 1.4, -1.7, -0.7, 0.4}):
                    g = -(gl + x * gr)
                    cos[f"{x:+.2f}"] = float(np.sum(dlt * g) / np.sqrt(np.sum(dlt * dlt) * np.sum(g * g)))
                xbest = max(cos, key=cos.get)
                rec["cos_displacement_vs_minus_grad"] = cos
                rec["x_best_by_cos"] = float(xbest)
        log(f"C3grad lam{lam} cR{cR:+} {nm}: x*(end) {rec['end']['x_star']:.3f} "
            f"fmax label/x* {rec['end']['fmax_ratio_label_over_xstar']:.2f}; "
            f"x*(seed) {rec['seed']['x_star']:.3f}; cos best x {rec['x_best_by_cos']:+.2f} "
            f"({rec['cos_displacement_vs_minus_grad']})")
        out[f"lam{lam}_cR{cR:+}_{nm}"] = rec
    RESULTS["c3grad"] = out


# ============================================================ C6
def stage_c6():
    out = {}
    # (a) ansatz: pure boost orbit, constant spectrum
    M = R3I.ansatz(VAC, [[0, 0, 4.0], [0, 0, -4.0]], 0.1, 1.0, 32, 48.0)
    lam = sorted_spec_own(M)
    out["ansatz_spectrum_spread"] = float(np.max(np.abs(lam - LAM_VAC)))
    M32 = R3I.ansatz(VAC32, [[0, 0, 4.0], [0, 0, -4.0]], 0.1, 1.0, 32, 48.0)
    out["ansatz_g32_K_lambda"] = klam_static_own(M32)
    out["ansatz_K_lambda"] = klam_static_own(M)
    log(f"C6 ansatz: K_lambda {out['ansatz_K_lambda']:.3e}, max |lam - vac| {out['ansatz_spectrum_spread']:.3e}")
    # (b) relaxed lambda = 0 pairs
    E = {}
    for d in (10, 14, 18, 24):
        Mp, p = load_r3(0, f"same_d{d}")
        E[d] = klam_static_own(Mp)
    Es = klam_static_own(load_r3(0, "single_d0")[0])
    out["K_lambda_pairs"] = {str(d): E[d] for d in E}
    out["K_lambda_single"] = Es
    out["E_minus_E24"] = {str(d): E[d] - E[24] for d in (10, 14, 18)}
    out["E_minus_2single"] = {str(d): E[d] - 2 * Es for d in (10, 14, 18, 24)}
    out["local_exp"] = {}
    for d1_, d2_ in ((10, 14), (14, 18)):
        y1, y2 = E[d1_] - E[24], E[d2_] - E[24]
        out["local_exp"][f"{d1_}_{d2_}"] = float(np.log(y1 / y2) / np.log(d2_ / d1_)) if y1 * y2 > 0 else None
    # mutation: the UNSORTED spectrum (eig order as returned) must not give the same numbers
    Eu = {}
    for d in (10, 24):
        Mp, _ = load_r3(0, f"same_d{d}")
        lamu = np.linalg.eigvals(Mp @ ETA).real
        tot = 0.0
        for ax in range(3):
            dd = np.diff(lamu, axis=ax) / H
            tot += 0.5 * np.sum(dd * dd)
        Eu[d] = float(H3 * tot)
    out["mutant_unsorted_E10_minus_E24"] = Eu[10] - Eu[24]
    log(f"C6 pairs: E(d)-E(24) {out['E_minus_E24']}, exps {out['local_exp']}, unsorted mutant {out['mutant_unsorted_E10_minus_E24']:.3e}")
    RESULTS["c6"] = out


# ============================================================ C7 / C9
def shell_profile(f, r, edges):
    prof = []
    for a, b in zip(edges[:-1], edges[1:]):
        m = (r >= a) & (r < b)
        prof.append((0.5 * (a + b), float(np.mean(f[m])) if np.any(m) else float("nan"), int(np.sum(m))))
    return prof


def stage_c79():
    out = {}
    X, Y, Z = B3.coords(32, H)
    r = np.sqrt(X * X + Y * Y + Z * Z)
    for it in (12000, 3000):
        M, p = load_hedgehog(it)
        lam = sorted_spec_own(M)
        dl = lam - LAM_VAC
        rec = {"sha": sha(p)}
        rec["monopole_h3"] = [float(np.sum(dl[..., a]) * H3) for a in range(4)]
        rec["dipole_z_h3"] = [float(np.sum(dl[..., a] * Z) * H3) for a in range(4)]
        rec["dipole_x_h3"] = [float(np.sum(dl[..., a] * X) * H3) for a in range(4)]
        ic = np.unravel_index(np.argmin(r), r.shape)
        rec["core_cell_r"] = float(r[ic])
        rec["core_values"] = [float(dl[ic + (a,)]) for a in range(4)]
        rec["extreme_values"] = [[float(np.min(dl[..., a])), float(np.max(dl[..., a]))] for a in range(4)]
        # far exponents: shell averages r in [6, 20], log-log slope; two binnings
        rec["far_exp"] = {}
        for width in (1.0, 2.0):
            edges = np.arange(6.0, 20.0 + 1e-9, width)
            ex = []
            for a in range(4):
                prof = shell_profile(dl[..., a], r, edges)
                rr = np.array([q[0] for q in prof])
                vv = np.array([q[1] for q in prof])
                ok = np.abs(vv) > 0
                if np.all(np.sign(vv[ok]) == np.sign(vv[ok][0])):
                    sl = np.polyfit(np.log(rr[ok]), np.log(np.abs(vv[ok])), 1)[0]
                    ex.append(float(sl))
                else:
                    ex.append(None)
            rec["far_exp"][f"bin{width:g}"] = ex
        # lambda_2 has a sign change inside [6, 20]: the inner window [6, 13] and the radius
        edges = np.arange(6.0, 13.0 + 1e-9, 1.0)
        rec["far_exp_6_13"] = []
        rec["sign_change_r"] = []
        for a in range(4):
            prof = shell_profile(dl[..., a], r, edges)
            rr = np.array([q[0] for q in prof]); vv = np.array([q[1] for q in prof])
            ok = np.abs(vv) > 0
            if np.all(np.sign(vv[ok]) == np.sign(vv[ok][0])):
                rec["far_exp_6_13"].append(float(np.polyfit(np.log(rr[ok]), np.log(np.abs(vv[ok])), 1)[0]))
            else:
                rec["far_exp_6_13"].append(None)
            prof2 = shell_profile(dl[..., a], r, np.arange(0.0, 24.0 + 1e-9, 1.0))
            vv2 = np.array([q[1] for q in prof2]); rr2 = np.array([q[0] for q in prof2])
            keep = np.isfinite(vv2) & (np.abs(vv2) > 1e-12)
            vv2, rr2 = vv2[keep], rr2[keep]
            sc = np.where(np.sign(vv2[1:]) != np.sign(vv2[:-1]))[0]
            rec["sign_change_r"].append(float(rr2[sc[0] + 1]) if len(sc) else None)
        # radial profile (bin 2) for the record
        edges = np.arange(0.0, 24.0 + 1e-9, 2.0)
        rec["profile_bin2"] = {f"lam{a + 1}": shell_profile(dl[..., a], r, edges) for a in range(4)}
        # trace deficit (C9) = tr(N) - tr(N_vac): linear in M, check both routes
        tr = np.sum(dl, axis=-1)
        trM = np.einsum("...aa->...", M @ ETA) - np.sum(LAM_VAC)
        rec["trace_deficit_max"] = float(np.max(tr))
        rec["trace_deficit_min"] = float(np.min(tr))
        im = np.unravel_index(np.argmax(tr), tr.shape)
        rec["trace_deficit_argmax_r"] = float(r[im])
        rec["trace_deficit_core"] = float(tr[ic])
        rec["trace_route_check"] = float(np.max(np.abs(tr - trM)))
        rec["trace_deficit_monopole_h3"] = float(np.sum(tr) * H3)
        # mutation: an UNSORTED spectrum destroys the per-channel monopoles (control)
        lamu = np.linalg.eigvals(M @ ETA).real
        rec["mutant_unsorted_monopoles_h3"] = [float(np.sum(lamu[..., a] - LAM_VAC[a]) * H3) for a in range(4)]
        log(f"C7/C9 it{it}: monopoles {np.round(rec['monopole_h3'], 1).tolist()} dipoles_z "
            f"{np.round(rec['dipole_z_h3'], 3).tolist()} core {np.round(rec['core_values'], 3).tolist()} "
            f"far_exp(bin1) {rec['far_exp']['bin1']} far_exp(bin2) {rec['far_exp']['bin2']} "
            f"trace max {rec['trace_deficit_max']:.4f} at r {rec['trace_deficit_argmax_r']:.2f}")
        out[f"it{it}"] = rec
    RESULTS["c7c9"] = out


# ============================================================ C8
def lap7(f, h):
    out = -6.0 * f
    for ax in range(3):
        out += np.roll(f, 1, axis=ax) + np.roll(f, -1, axis=ax)
    return out / (h * h)


def stage_c8():
    out = {}
    M, p = load_hedgehog(12000)
    lam = sorted_spec_own(M)
    dl = lam - LAM_VAC
    n = 32
    k = 2 * np.pi * np.fft.fftfreq(n, d=H)
    KX, KY, KZ = np.meshgrid(k, k, k, indexing="ij")
    K2 = KX ** 2 + KY ** 2 + KZ ** 2
    ds = [6, 9, 12, 18, 24]
    rng = np.random.default_rng(11)
    for a in range(4):
        phi = dl[..., a]
        ph_k = np.fft.fftn(phi)
        rec = {}
        for ms in (0.1, 0.3, 1.0):
            # kernel K1: real-space 7-point Laplacian source + roll shift
            J = -lap7(phi, H) + ms * ms * phi
            e1 = []
            e2 = []
            e3 = []
            e_rand = []
            signs = rng.choice([-1.0, 1.0], size=K2.shape)
            for d in ds:
                sh = int(round(d / H))
                assert abs(sh * H - d) < 1e-9
                e1.append(float(-np.sum(J * np.roll(phi, sh, axis=2)) * H3))
                # kernel K2: Parseval with the spectral k^2 (independent of K1's stencil)
                val = -np.sum(np.abs(ph_k) ** 2 * (K2 + ms * ms) * np.cos(KZ * d)) * H3 / n ** 3
                e2.append(float(val.real))
                # K3: spectral source J = (-lap + m^2) phi with the SPECTRAL laplacian, roll shift
                Jk = np.fft.ifftn((K2 + ms * ms) * ph_k).real
                e3.append(float(-np.sum(Jk * np.roll(phi, sh, axis=2)) * H3))
                # mutation: random-sign kernel
                vr = -np.sum(np.abs(ph_k) ** 2 * (K2 + ms * ms) * signs * np.cos(KZ * d)) * H3 / n ** 3
                e_rand.append(float(vr.real))
            mono = lambda v: bool(np.all(np.diff(v) > 0))
            surv = 0
            for sd in range(20):
                sg = np.random.default_rng(100 + sd).choice([-1.0, 1.0], size=K2.shape)
                vals = [float((-np.sum(np.abs(ph_k) ** 2 * (K2 + ms * ms) * sg * np.cos(KZ * d)) * H3 / n ** 3).real) for d in ds]
                surv += int(mono(vals) and np.all(np.array(vals) < 0))
            sgx = np.random.default_rng(5).choice([-1.0, 1.0], size=phi.shape)
            e_src = [float(-np.sum((sgx * J) * np.roll(phi, int(round(d / H)), axis=2)) * H3) for d in ds]
            rec[f"ms{ms:g}"] = {
                "d": ds, "E_int_K1_lap7": e1, "E_int_K2_parseval": e2, "E_int_K3_spectral_source": e3,
                "all_negative_K1": bool(np.all(np.array(e1) < 0)),
                "rising_K1": mono(e1), "rising_K2": mono(e2), "rising_K3": mono(e3),
                "mutant_random_sign": e_rand, "mutant_rising": mono(e_rand),
                "mutant_all_negative": bool(np.all(np.array(e_rand) < 0)),
                "mutant_kernel_20seeds_negative_and_rising_survivors": surv,
                "mutant_random_sign_source": e_src,
                "mutant_source_negative_and_rising": bool(mono(e_src) and np.all(np.array(e_src) < 0)),
            }
            log(f"C8 lam{a + 1} ms {ms:g}: K1 {np.round(e1, 2).tolist()} K3 {np.round(e3, 2).tolist()} "
                f"rising {mono(e1)}/{mono(e3)} | mutant rising {mono(e_rand)}")
        out[f"lam{a + 1}"] = rec
    # also: the d-dependence of the like-source exchange for a POINT-like source of the same
    # monopole (Yukawa e^{-m d}/d reference), to see whether the shape is source-dominated
    RESULTS["c8"] = out


# ============================================================ main
STAGES = {"sanity": stage_sanity, "c12": stage_c12, "c3": stage_c3, "c3grad": stage_c3grad,
          "c4": stage_c4, "c5": stage_c5, "c5h": stage_c5h, "c6": stage_c6, "c79": stage_c79,
          "c8": stage_c8}


def main():
    want = ARGV[1:] or list(STAGES)
    for s in want:
        log(f"=== stage {s}")
        STAGES[s]()
    outp = os.path.join(DATA, "m5_32_r14_c_audit.json")
    prev = {}
    if os.path.exists(outp) and len(want) < len(STAGES):
        prev = json.load(open(outp))
    prev.update(RESULTS)
    prev["wall_s"] = time.time() - T0
    json.dump(prev, open(outp, "w"), indent=1, default=float)
    log(f"wrote {outp}")


if __name__ == "__main__":
    main()
