"""M5.32 R14-B adversarial audit: the fixed-J continuation on L_cert + c K_P^h.

Independent re-measurement of the eight R14-B claims on the SAVED end fields
(checkpoints/m5_32_r14/b_fixedj/<tag>.npy) with the auditor's own code:
own one-sided stencil, own eta-orthonormal eigenbasis of N = M eta, own K_P^h
static and omega^2 densities (Frobenius norm of f_a f_b B_ab), own I1 inertia
density, own E_u density, own finite-difference gate on the assembled fixed-J
force, own 50-step FIRE continuation, own (2,3)-gap and radial statistics, and
the three mutations named by the brief (plain K_P, literal roots, global G1).
The producer's script and result file (m5_32_r14_b_fixedj.py / .json) are never
opened; the certified stack (m5_21_3_a_4d.py), the term module
(m5_32_r14_terms.py) and the R13-W helpers are used only as the allowed
cross-check and for the gradient pieces the brief allows.

EQUATIONS (the brief's definitions)
  E_J[M] = E_u + V4 + c K_P^h_stat[M] + J^2 / (4 kin_tot),
  kin_tot = kin_I1(M; a0) + c kin_KP(M; a0),  a0 = a0_local(M),  omega = J / (2 kin_tot)
  K_P^h_stat density = (1/2) sum_i sum_ab f_a^2 f_b^2 (B_i)_ab^2,  B_i = S V^T eta A_i eta V,
  f_a = (lambda_a + g)(lambda_a - 1),  (lambda, V, S) the eta-orthonormal eigensystem of N = M eta
  kin_KP density = (1/2) sum_ab f_a^2 f_b^2 (B_0)_ab^2 with A_0 = a0
  kin_I1 density = 4 sum_i <[a0, A_i]_eta, [a0, A_i]_eta>_eta
  every density h^3-weighted and branch-averaged over the sym stencil (fwd, bwd) x 1/2.

Usage: python3 m5_32_r14_b_audit.py [steps_n32] [steps_n48]  (defaults 50, 50; both <= 100)
Out:   ../data/m5_32_r14_b_audit.json
"""
from __future__ import annotations

import sys

ARGV = list(sys.argv)                       # captured before any import

import importlib.util
import json
import os
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.dirname(HERE)
DATA = os.path.join(RES, "data")
CKB = os.path.join(RES, "checkpoints", "m5_32_r14", "b_fixedj")
CKW = os.path.join(RES, "checkpoints", "m5_32_r13w")
SEED32 = os.path.join(RES, "checkpoints", "m5_32_r10", "relax_g8_n32_L48_it12000.npy")
T0 = time.time()


def log(m):
    print(f"[{time.time() - T0:8.1f}s] {m}", flush=True)


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


TERMS = _load("m5_32_r14_terms", "m5_32_r14_terms.py")      # allowed: kp_h_energy_grad, roots_for, kp_h_static
COMMON = _load("m5_32_r13w_common", "m5_32_r13w_common.py")  # allowed: a0_local, cfg_of, gap23, kin_density
STACK = COMMON.INS4                                          # allowed: e_parts, grad, kin_of, kin_grad, pin_shell

ETA = np.diag([-1.0, 1.0, 1.0, 1.0])
ETA_D = np.diag(ETA)
G = 8.0
ROOTS = (-G, 1.0)                 # (lambda_t, 1): the intended roots in this convention
ROOTS_LITERAL = (G, 1.0)          # the literal reading (mutation)
G1 = np.zeros((4, 4)); G1[2, 3], G1[3, 2] = -1.0, 1.0        # the (2,3) rotation: the global clock

FIELDS = [("n32_c0.3_J200", 32, 48.0, 0.3, 200.0, False),
          ("n32_c1_J200", 32, 48.0, 1.0, 200.0, False),
          ("n32_c3_J200", 32, 48.0, 3.0, 200.0, False),
          ("n32_c1_J50", 32, 48.0, 1.0, 50.0, False),
          ("n32_c1_J800", 32, 48.0, 1.0, 800.0, False),
          ("n32_c1_J200_frz", 32, 48.0, 1.0, 200.0, True),
          ("n48L48_c1_J200", 48, 48.0, 1.0, 200.0, False),
          ("n48L72_c1_J200", 48, 72.0, 1.0, 200.0, False)]


# ============================================================ own building blocks
def own_d1(f, ax, h, br):
    """one-sided difference along axis ax (fwd: cells 0..n-2, bwd: cells 1..n-1, else 0)."""
    out = np.zeros_like(f)
    sl = [slice(None)] * f.ndim
    lo, hi = list(sl), list(sl)
    if br == "fwd":
        lo[ax], hi[ax] = slice(0, -1), slice(1, None)
        out[tuple(lo)] = (f[tuple(hi)] - f[tuple(lo)]) / h
    else:
        lo[ax], hi[ax] = slice(0, -1), slice(1, None)
        out[tuple(hi)] = (f[tuple(hi)] - f[tuple(lo)]) / h
    return out


def own_jets(M, h, br):
    return [own_d1(M, ax, h, br) for ax in range(3)]


def own_eig(M):
    """eta-orthonormal eigensystem of N = M eta: lam (...,4), V (...,4,4) columns, sig (...,4)."""
    lam, V = np.linalg.eig(M @ ETA)
    if np.max(np.abs(lam.imag)) > 1e-8 * max(float(np.max(np.abs(lam.real))), 1.0):
        raise ValueError("complex spectrum")
    lam, V = lam.real, V.real
    n2 = np.einsum("...ak,a,...ak->...k", V, ETA_D, V)
    V = V / np.sqrt(np.abs(n2))[..., None, :]
    return lam, V, np.sign(n2)


def own_kph_cells(A_list, M, roots, weights="h"):
    """per-cell (1/2) sum_ab w_ab f_a^2 f_b^2 B_ab^2 summed over the matrices in A_list;
    weights 'h' = Frobenius (K_P^h), 'plain' = sig_a sig_b (the plain trace K_P)."""
    lam, V, sig = own_eig(M)
    p1, p2 = roots
    f = (lam - p1) * (lam - p2)
    eV = ETA @ V
    w = (f[..., :, None] * f[..., None, :]) ** 2
    if weights == "plain":
        # plain trace form: (1/2) sum_ab f_a^2 f_b^2 sig_a sig_b B_ab^2
        w = w * (sig[..., :, None] * sig[..., None, :])
    out = np.zeros(M.shape[:-2])
    for A in A_list:
        # B_ab = sig_a (eta v_a)^T A (eta v_b)   (S V^T eta A eta V)
        X = np.einsum("...ca,...cd,...db->...ab", eV, A, eV) * sig[..., :, None]
        out += 0.5 * np.sum(w * X * X, axis=(-1, -2))
    return out


def own_kph_static_density(M, cfg, roots=ROOTS, weights="h"):
    h = cfg["h"]
    d = np.zeros(M.shape[:-2])
    for br in ("fwd", "bwd"):
        d += 0.5 * own_kph_cells(own_jets(M, h, br), M, roots, weights)
    return h ** 3 * d


def own_kph_kin_density(M, a0, cfg, roots=ROOTS):
    return cfg["h"] ** 3 * own_kph_cells([a0], M, roots)


def own_comm_eta(A, B):
    return A @ ETA @ B - B @ ETA @ A


def own_inner_eta(F, Gm):
    return np.einsum("...ab,...cd,ac,bd->...", F, Gm, ETA, ETA, optimize=True)


def own_kin_I1_density(M, a0, cfg):
    h = cfg["h"]
    d = np.zeros(M.shape[:-2])
    for br in ("fwd", "bwd"):
        for A in own_jets(M, h, br):
            F = own_comm_eta(a0, A)
            d += 0.5 * 4.0 * own_inner_eta(F, F)
    return h ** 3 * d


def own_e_u_density(M, cfg):
    h = cfg["h"]
    d = np.zeros(M.shape[:-2])
    for br in ("fwd", "bwd"):
        A = own_jets(M, h, br)
        for i in range(3):
            for j in range(i + 1, 3):
                F = own_comm_eta(A[i], A[j])
                d += 0.5 * 4.0 * own_inner_eta(F, F)
    return h ** 3 * d


def own_gap23(M):
    w = np.linalg.eigvalsh(M[..., 1:, 1:])
    return w[..., 1] - w[..., 0]


def radius(n, h):
    x = (np.arange(n) - (n - 1) / 2.0) * h
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    return np.sqrt(X * X + Y * Y + Z * Z)


# ============================================================ the fixed-J functional
def energy_J(M, cfg, c, J, a0=None):
    """own energies; a0 = a0_local(M) unless given (frozen generator)."""
    if a0 is None:
        a0 = COMMON.a0_local(M)
    e_u_own = float(np.sum(own_e_u_density(M, cfg)))
    e_u, v4 = STACK.e_parts(M, cfg)
    ks = float(np.sum(own_kph_static_density(M, cfg)))
    kI = float(np.sum(own_kin_I1_density(M, a0, cfg)))
    kK = float(np.sum(own_kph_kin_density(M, a0, cfg)))
    kin_tot = kI + c * kK
    fixj = J * J / (4.0 * kin_tot)
    return {"E_u": float(e_u), "E_u_own": e_u_own, "V4": float(v4), "E_KP": ks, "kin_I1": kI, "kin_KP": kK,
            "kin_tot": kin_tot, "fixedJ_term": fixj, "omega": J / (2.0 * kin_tot),
            "E_J": float(e_u) + float(v4) + c * ks + fixj}


def grad_J(M, cfg, c, J, a0=None):
    """the protocol's gradient: certified grad + c G_stat - J^2/(4 kin_tot^2) (kin_grad_I1 + c G_kin),
    a0 frozen inside the gradient."""
    if a0 is None:
        a0 = COMMON.a0_local(M)
    Gt = STACK.grad(M, cfg)
    Es, Gs, kK, Gk = TERMS.kp_h_energy_grad(M, cfg, ROOTS, a0)
    kI = float(STACK.kin_of(M, a0, cfg))
    kin_tot = kI + c * kK
    Gt = Gt + c * Gs - (J * J / (4.0 * kin_tot * kin_tot)) * (STACK.kin_grad(M, a0, cfg) + c * Gk)
    return Gt


def fd_gate(M, cfg, c, J, free, rng, eps=1e-4, a0_frozen=True):
    """directional derivative of the OWN energy along a random free symmetric direction vs the
    assembled gradient (a0 frozen at M, as the protocol prescribes; and refreshed, as a probe)."""
    D = rng.normal(size=M.shape)
    D = 0.5 * (D + D.swapaxes(-1, -2)) * free[..., None, None]
    D /= np.sqrt(np.sum(D * D))
    a0 = COMMON.a0_local(M)
    Gt = grad_J(M, cfg, c, J, a0)
    an = float(np.sum(Gt * D))
    ep = energy_J(M + eps * D, cfg, c, J, a0)["E_J"]
    em = energy_J(M - eps * D, cfg, c, J, a0)["E_J"]
    fd_frozen = (ep - em) / (2 * eps)
    ep2 = energy_J(M + eps * D, cfg, c, J)["E_J"]
    em2 = energy_J(M - eps * D, cfg, c, J)["E_J"]
    fd_refr = (ep2 - em2) / (2 * eps)
    return {"analytic": an, "fd_frozen_a0": fd_frozen, "fd_refreshed_a0": fd_refr,
            "rel_err_frozen": abs(fd_frozen - an) / max(abs(an), 1e-300),
            "rel_err_refreshed": abs(fd_refr - an) / max(abs(an), 1e-300)}


def own_fire(M0, cfg, c, J, free, steps, a0_fixed=None, dt0=0.01, dt_max=0.1, tag=""):
    """own FIRE (Bitzek 2006 rules, the protocol's dt0 / dt_max / alpha 0.1 / N_min 5)."""
    M = M0.copy()
    fm = free[..., None, None].astype(float)
    v = np.zeros_like(M)
    dt, alpha, n_up = dt0, 0.1, 0
    F = -grad_J(M, cfg, c, J, a0_fixed) * fm
    e0 = energy_J(M, cfg, c, J, a0_fixed)
    rows = [{"it": 0, "E_J": e0["E_J"], "E_KP": e0["E_KP"], "fmax": float(np.max(np.abs(F))), "dt": dt}]
    for it in range(1, steps + 1):
        P = float(np.sum(F * v))
        if P > 0.0:
            n_up += 1
            vn, fn = np.sqrt(np.sum(v * v)), np.sqrt(np.sum(F * F))
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
        F = -grad_J(M, cfg, c, J, a0_fixed) * fm
        if it % 10 == 0 or it == steps:
            e = energy_J(M, cfg, c, J, a0_fixed)
            rows.append({"it": it, "E_J": e["E_J"], "E_KP": e["E_KP"], "fmax": float(np.max(np.abs(F))), "dt": dt})
            log(f"  {tag} own FIRE it {it:3d} E_J {e['E_J']:.6f} dE {e['E_J'] - e0['E_J']:+.4e} fmax {rows[-1]['fmax']:.3e} dt {dt:.4f}")
    return M, rows


# ============================================================ helpers
def shell_profile(dens, r, edges):
    out = []
    for a, b in zip(edges[:-1], edges[1:]):
        m = (r >= a) & (r < b)
        out.append({"r_lo": float(a), "r_hi": float(b), "n": int(m.sum()),
                    "mean": float(dens[m].mean()) if m.any() else None})
    return out


def rel_diff(A, B):
    return float(np.sqrt(np.sum((A - B) ** 2)) / np.sqrt(np.sum(B * B)))


def load_seed(n, L):
    if n == 32:
        return np.load(SEED32)
    return np.load(os.path.join(CKW, f"seed_n{n}_L{L:g}_it3000.npy"))


# ============================================================ main
def main():
    steps32 = int(ARGV[1]) if len(ARGV) > 1 else 50
    steps48 = int(ARGV[2]) if len(ARGV) > 2 else 50
    assert steps32 <= 100 and steps48 <= 100, "the brief caps own descents at 100 steps"
    rng = np.random.default_rng(1414)
    res = {"argv": ARGV, "steps": {"n32": steps32, "n48": steps48}, "fields": {}, "seeds": {}, "traces": {}}

    # -------- the vacuum reference for kin_KP (own analytic + own numeric)
    cfg32 = COMMON.cfg_of(32, 48.0)
    Mv = np.zeros((2, 2, 2, 4, 4)) + STACK.vac4(cfg32)
    a0v = COMMON.a0_local(Mv)
    dv = own_kph_cells([a0v], Mv, ROOTS)
    f2, f3 = (0.3 + G) * (0.3 - 1.0), (0.0 + G) * (0.0 - 1.0)
    pred = 0.5 * 2.0 * (f2 * f3) ** 2 * 0.3 ** 2
    res["vacuum_kinKP_cell"] = {"own_numeric_no_h3": float(dv[0, 0, 0]), "analytic_no_h3": pred,
                                "a0_local_equals_G1_clock": float(np.max(np.abs(a0v - (G1 @ Mv - Mv @ G1))))}
    log(f"vacuum kin_KP per cell (no h^3): own {dv[0,0,0]:.4f} analytic {pred:.4f}; x h^3 (1.5) = {pred*3.375:.2f}, (1.0) = {pred:.2f}")

    seeds_loaded = {}
    for tag, n, L, c, J, frz in FIELDS:
        cfg = COMMON.cfg_of(n, L)
        h = cfg["h"]
        M = np.load(os.path.join(CKB, tag + ".npy"))
        rec = json.load(open(os.path.join(CKB, tag + ".json")))
        res["traces"][tag] = {"stop": rec["stop"], "iters": rec["iters"], "trace": rec["trace"],
                              "maxit": rec["maxit"], "frozen_generator": rec["frozen_generator"]}
        free = ~STACK.pin_shell(n, h)
        r = radius(n, h)
        key = (n, L)
        if key not in seeds_loaded:
            seeds_loaded[key] = load_seed(n, L)
        S0 = seeds_loaded[key]
        a0_seed = COMMON.a0_local(S0)
        a0_fixed = a0_seed if frz else None
        out = {"n": n, "L": L, "h": h, "c": c, "J": J, "frozen_generator": frz}
        log(f"=== {tag} (n {n}, L {L:g}, h {h:g}, c {c}, J {J}, frozen {frz})")

        # sanity: pinned shell at the vacuum, no boosts
        vac = STACK.vac4(cfg)
        out["sanity"] = {"max_abs_M0i": float(np.abs(M[..., 0, 1:]).max()),
                         "pinned_shell_dev_from_vac": float(np.abs(M[~free] - vac).max()),
                         "seed_pinned_shell_dev": float(np.abs(S0[~free] - vac).max()),
                         "n_free": int(free.sum()), "n_pinned": int((~free).sum())}

        # -------- B1 energies (own) + cross-check against the allowed stack
        e = energy_J(M, cfg, c, J, a0_fixed)
        a0 = a0_fixed if frz else COMMON.a0_local(M)
        Es_mod, _, kK_mod, _ = TERMS.kp_h_energy_grad(M, cfg, ROOTS, a0)
        kI_mod = float(STACK.kin_of(M, a0, cfg))
        e_local = energy_J(M, cfg, c, J)          # local a0 refreshed on the end field (for frz)
        out["B1"] = dict(e)
        out["B1"]["cross_check"] = {"E_KP_module": float(Es_mod), "kin_KP_module": float(kK_mod), "kin_I1_stack": kI_mod,
                                    "rel_E_KP": abs(e["E_KP"] - Es_mod) / Es_mod, "rel_kin_KP": abs(e["kin_KP"] - kK_mod) / kK_mod,
                                    "rel_kin_I1": abs(e["kin_I1"] - kI_mod) / kI_mod, "rel_E_u": abs(e["E_u_own"] - e["E_u"]) / e["E_u"]}
        out["B1"]["with_local_a0_on_end_field"] = {"kin_I1": e_local["kin_I1"], "kin_KP": e_local["kin_KP"], "omega": e_local["omega"]}
        log(f"  B1 E_u {e['E_u']:.4f} (own {e['E_u_own']:.4f}) V4 {e['V4']:.4f} E_KP {e['E_KP']:.2f} (mod {Es_mod:.2f}) "
            f"kin_I1 {e['kin_I1']:.2f} kin_KP {e['kin_KP']:.4e} omega {e['omega']:.4e} fixedJ {e['fixedJ_term']:.4e}")

        # -------- B2 force + FD gate + own continuation
        Gt = grad_J(M, cfg, c, J, a0_fixed)
        F = -Gt * free[..., None, None]
        fmax = float(np.max(np.abs(F)))
        gate = fd_gate(M, cfg, c, J, free, rng)
        fm = free[..., None, None]
        G_cert = STACK.grad(M, cfg) * fm
        G_kp = c * TERMS.kp_h_energy_grad(M, cfg, ROOTS)[1] * fm
        G_fix = (F + G_cert + G_kp) * (-1.0)                      # what remains: the fixed-J part
        cellmax = np.max(np.abs(F), axis=(-1, -2))
        ijk = np.unravel_index(int(np.argmax(cellmax)), cellmax.shape)
        out["B2"] = {"fmax_full": fmax, "fmax_static_only": float(np.max(np.abs(G_cert + G_kp))),
                     "fmax_certified_part_only": float(np.max(np.abs(G_cert))),
                     "fmax_KP_part_only": float(np.max(np.abs(G_kp))), "fmax_fixedJ_part_only": float(np.max(np.abs(G_fix))),
                     "argmax_cell": [int(q) for q in ijk], "argmax_r": float(r[ijk]),
                     "at_argmax_abs_cert": float(np.max(np.abs(G_cert[ijk]))), "at_argmax_abs_KP": float(np.max(np.abs(G_kp[ijk]))),
                     "force_norm_full": float(np.sqrt(np.sum(F * F))), "force_norm_cert": float(np.sqrt(np.sum(G_cert ** 2))),
                     "force_norm_KP": float(np.sqrt(np.sum(G_kp ** 2))),
                     "frac_free_cells_fmax_above_1e-6": float(np.mean(cellmax[free] > 1e-6)),
                     "frac_free_cells_fmax_above_1e-2": float(np.mean(cellmax[free] > 1e-2)),
                     "fd_gate": gate, "trace_last_fmax": rec["trace"][-1]["fmax"]}
        log(f"  B2 fmax {fmax:.4e} (L_cert part {out['B2']['fmax_certified_part_only']:.4e}, K_P part {out['B2']['fmax_KP_part_only']:.4e}, "
            f"fixed-J part {out['B2']['fmax_fixedJ_part_only']:.2e}) at r {out['B2']['argmax_r']:.2f}; norms cert {out['B2']['force_norm_cert']:.3e} "
            f"KP {out['B2']['force_norm_KP']:.3e}; FD gate frozen-a0 rel {gate['rel_err_frozen']:.2e}, refreshed-a0 rel {gate['rel_err_refreshed']:.2e}")
        steps = steps32 if n == 32 else steps48
        M2, rows = own_fire(M, cfg, c, J, free, steps, a0_fixed, tag=tag)
        out["B2"]["own_fire"] = {"steps": steps, "rows": rows, "dE": rows[-1]["E_J"] - rows[0]["E_J"],
                                 "dE_KP": rows[-1]["E_KP"] - rows[0]["E_KP"],
                                 "fmax_end": rows[-1]["fmax"], "monotone": all(rows[i + 1]["E_J"] <= rows[i]["E_J"] for i in range(len(rows) - 1)),
                                 "field_move_rel": rel_diff(M2, M), "field_move_max": float(np.max(np.abs(M2 - M)))}
        log(f"  B2 own FIRE {steps} steps: dE {out['B2']['own_fire']['dE']:+.4e}, fmax end {rows[-1]['fmax']:.3e}")

        # -------- B4 kin_KP density map and shells (own), B7 kin_I1 radial
        dK = own_kph_kin_density(M, a0, cfg)
        dI = own_kin_I1_density(M, a0, cfg)
        vac_cell = pred * h ** 3
        edges = np.arange(0.0, r.max() + 1.5, 1.5)
        prof = shell_profile(dK / vac_cell, r, edges)
        gap = own_gap23(M)
        out["B4"] = {"vac_cell_value": vac_cell, "kinKP_over_vac_shells": prof,
                     "frac_cells_within_5pct_of_vac": float(np.mean(np.abs(dK / vac_cell - 1) < 0.05)),
                     "frac_free_cells_within_5pct_of_vac": float(np.mean(np.abs(dK[free] / vac_cell - 1) < 0.05)),
                     "mean_ratio_free_r_gt_12": float(np.mean(dK[free & (r > 12)] / vac_cell)),
                     "mean_ratio_free_r_gt_18": float(np.mean(dK[free & (r > 18)] / vac_cell)),
                     "mean_ratio_r_lt_12": float(np.mean(dK[r < 12] / vac_cell)),
                     "kinKP_total_over_n3_vac": e["kin_KP"] / (vac_cell * n ** 3),
                     "kinKP_frac_r_gt_12": float(dK[r > 12].sum() / dK.sum()),
                     "kinKP_frac_pinned": float(dK[~free].sum() / dK.sum()),
                     "gap_ratio_sq_mean_free_r_gt_12": float(np.mean((gap[free & (r > 12)] / 0.3) ** 2)),
                     "corr_kinKP_gap2_free": float(np.corrcoef(dK[free], gap[free] ** 2)[0, 1])}
        log(f"  B4 kin_KP/vac: free r>12 mean {out['B4']['mean_ratio_free_r_gt_12']:.3f}, r>18 {out['B4']['mean_ratio_free_r_gt_18']:.3f}, "
            f"r<12 {out['B4']['mean_ratio_r_lt_12']:.3f}; total/(n^3 vac) {out['B4']['kinKP_total_over_n3_vac']:.3f}; "
            f"free cells within 5% {out['B4']['frac_free_cells_within_5pct_of_vac']:.3f}")
        out["B7"] = {"kinI1_frac_r_gt_12": float(dI[r > 12].sum() / dI.sum()),
                     "kinI1_frac_r_gt_6": float(dI[r > 6].sum() / dI.sum()),
                     "kinI1_frac_pinned": float(dI[~free].sum() / dI.sum()),
                     "kinI1_shells": shell_profile(dI, r, edges)}
        log(f"  B7 kin_I1 frac r>12: {out['B7']['kinI1_frac_r_gt_12']:.3f} (r>6 {out['B7']['kinI1_frac_r_gt_6']:.3f})")

        # -------- B5 gap statistics (own)
        out["B5"] = {"gap_mean_free": float(gap[free].mean()), "gap_min_free": float(gap[free].min()),
                     "gap_frac_free_below_0.05": float(np.mean(gap[free] < 0.05)),
                     "gap_frac_free_below_0.1": float(np.mean(gap[free] < 0.1)),
                     "gap_frac_free_above_0.29": float(np.mean(gap[free] > 0.29)),
                     "gap_hist_free": np.histogram(gap[free], bins=[0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.29, 0.31, 0.5, 2.0])[0].tolist(),
                     "gap_mean_free_r_gt_12": float(gap[free & (r > 12)].mean()),
                     "gap_mean_r_lt_6": float(gap[r < 6].mean())}
        log(f"  B5 gap free: mean {out['B5']['gap_mean_free']:.4f} min {out['B5']['gap_min_free']:.4f} frac<0.05 {out['B5']['gap_frac_free_below_0.05']:.4f} "
            f"frac>0.29 {out['B5']['gap_frac_free_above_0.29']:.3f}")

        # -------- mutations
        plain = float(np.sum(own_kph_static_density(M, cfg, ROOTS, "plain")))
        plain_mod = float(TERMS.static_energy("K_P", M, cfg))
        kK_lit = float(np.sum(own_kph_kin_density(M, a0, cfg, ROOTS_LITERAL)))
        Es_lit = float(np.sum(own_kph_static_density(M, cfg, ROOTS_LITERAL)))
        a0g = G1 @ M - M @ G1
        kI_G1 = float(np.sum(own_kin_I1_density(M, a0g, cfg)))
        kK_G1 = float(np.sum(own_kph_kin_density(M, a0g, cfg)))
        out["mutations"] = {"plain_KP_static_own": plain, "plain_KP_static_module": plain_mod,
                            "plain_minus_h_rel": (plain - e["E_KP"]) / e["E_KP"],
                            "literal_roots_kin_KP": kK_lit, "literal_over_intended_kin_KP": kK_lit / e["kin_KP"],
                            "literal_roots_E_KP": Es_lit,
                            "G1_kin_I1": kI_G1, "G1_over_local_kin_I1": kI_G1 / e["kin_I1"],
                            "G1_kin_KP": kK_G1, "G1_over_local_kin_KP": kK_G1 / e["kin_KP"]}
        log(f"  MUT plain K_P {plain:.2f} (module {plain_mod:.2f}) vs K_P^h {e['E_KP']:.2f}; literal kin_KP x{kK_lit/e['kin_KP']:.3e}; "
            f"G1 kin_I1 {kI_G1:.2f} (x{kI_G1/e['kin_I1']:.3f}), G1 kin_KP x{kK_G1/e['kin_KP']:.3f}")
        res["fields"][tag] = out

        # -------- seed statistics (once per (n, L))
        skey = f"seed_n{n}_L{L:g}"
        if skey not in res["seeds"]:
            es = energy_J(S0, cfg, 1.0, 200.0)
            gs = own_gap23(S0)
            dKs = own_kph_kin_density(S0, a0_seed, cfg)
            dIs = own_kin_I1_density(S0, a0_seed, cfg)
            Fs = -grad_J(S0, cfg, 1.0, 200.0) * free[..., None, None]
            res["seeds"][skey] = {"E_u": es["E_u"], "V4": es["V4"], "E_KP": es["E_KP"], "kin_I1": es["kin_I1"], "kin_KP": es["kin_KP"],
                                  "omega_c1_J200": es["omega"], "fmax_c1_J200": float(np.max(np.abs(Fs))),
                                  "fmax_Lcert_only": float(np.max(np.abs(STACK.grad(S0, cfg) * free[..., None, None]))),
                                  "gap_mean_free": float(gs[free].mean()), "gap_min_free": float(gs[free].min()),
                                  "gap_frac_free_below_0.05": float(np.mean(gs[free] < 0.05)),
                                  "gap_hist_free": np.histogram(gs[free], bins=[0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.29, 0.31, 0.5, 2.0])[0].tolist(),
                                  "kinKP_mean_ratio_free_r_gt_12": float(np.mean(dKs[free & (r > 12)] / vac_cell)),
                                  "kinI1_frac_r_gt_12": float(dIs[r > 12].sum() / dIs.sum())}
            s = res["seeds"][skey]
            log(f"  SEED {skey}: E_u {s['E_u']:.4f} E_KP {s['E_KP']:.2f} kin_I1 {s['kin_I1']:.2f} kin_KP {s['kin_KP']:.4e} "
                f"gap mean {s['gap_mean_free']:.4f} frac<0.05 {s['gap_frac_free_below_0.05']:.4f} kinI1 frac r>12 {s['kinI1_frac_r_gt_12']:.3f} "
                f"fmax(c1,J200) {s['fmax_c1_J200']:.3e}")

    # -------- B3: field differences across J and B6: frozen vs refreshed
    F32 = {t: np.load(os.path.join(CKB, t + ".npy")) for t in ("n32_c1_J50", "n32_c1_J200", "n32_c1_J800", "n32_c1_J200_frz",
                                                                "n32_c0.3_J200", "n32_c3_J200")}
    seed32 = seeds_loaded[(32, 48.0)]
    base = F32["n32_c1_J200"]
    res["B3"] = {"rel_J50_vs_J200": rel_diff(F32["n32_c1_J50"], base), "rel_J800_vs_J200": rel_diff(F32["n32_c1_J800"], base),
                 "max_abs_J50_vs_J200": float(np.max(np.abs(F32["n32_c1_J50"] - base))),
                 "max_abs_J800_vs_J200": float(np.max(np.abs(F32["n32_c1_J800"] - base))),
                 "rel_J200_vs_seed": rel_diff(base, seed32), "max_abs_J200_vs_seed": float(np.max(np.abs(base - seed32))),
                 "rel_J50_vs_J200_over_move": rel_diff(F32["n32_c1_J50"], base) / rel_diff(base, seed32),
                 "rel_J800_vs_J200_over_move": rel_diff(F32["n32_c1_J800"], base) / rel_diff(base, seed32),
                 "omega": {t: res["fields"][t]["B1"]["omega"] for t in ("n32_c1_J50", "n32_c1_J200", "n32_c1_J800")},
                 "fixedJ_term": {t: res["fields"][t]["B1"]["fixedJ_term"] for t in ("n32_c1_J50", "n32_c1_J200", "n32_c1_J800")},
                 "rel_c0.3_vs_c1": rel_diff(F32["n32_c0.3_J200"], base), "rel_c3_vs_c1": rel_diff(F32["n32_c3_J200"], base)}
    om = res["B3"]["omega"]
    res["B3"]["omega_over_J"] = {t: om[t] / J for t, J in (("n32_c1_J50", 50.0), ("n32_c1_J200", 200.0), ("n32_c1_J800", 800.0))}
    log(f"B3 rel diff J50/J200 {res['B3']['rel_J50_vs_J200']:.3e}, J800/J200 {res['B3']['rel_J800_vs_J200']:.3e} "
        f"(J200 vs seed {res['B3']['rel_J200_vs_seed']:.3e}); omega/J {res['B3']['omega_over_J']}")
    frz, ref = F32["n32_c1_J200_frz"], base
    cfg = COMMON.cfg_of(32, 48.0)
    a0_seed = COMMON.a0_local(seed32)
    res["B6"] = {"rel_frz_vs_refreshed": rel_diff(frz, ref), "max_abs_frz_vs_refreshed": float(np.max(np.abs(frz - ref))),
                 "rel_frz_vs_refreshed_over_move": rel_diff(frz, ref) / rel_diff(ref, seed32),
                 "kin_I1_frozen_a0_on_frz_field": float(np.sum(own_kin_I1_density(frz, a0_seed, cfg))),
                 "kin_I1_local_a0_on_frz_field": float(np.sum(own_kin_I1_density(frz, COMMON.a0_local(frz), cfg))),
                 "kin_I1_frozen_a0_on_refreshed_field": float(np.sum(own_kin_I1_density(ref, a0_seed, cfg))),
                 "kin_I1_local_a0_on_refreshed_field": float(np.sum(own_kin_I1_density(ref, COMMON.a0_local(ref), cfg))),
                 "kin_I1_seed_a0_on_seed": float(np.sum(own_kin_I1_density(seed32, a0_seed, cfg))),
                 "a0_seed_vs_a0_local_end_rel": rel_diff(a0_seed, COMMON.a0_local(ref)),
                 "kin_KP_frozen_a0_on_frz_field": float(np.sum(own_kph_kin_density(frz, a0_seed, cfg))),
                 "kin_KP_local_a0_on_frz_field": float(np.sum(own_kph_kin_density(frz, COMMON.a0_local(frz), cfg)))}
    log(f"B6 rel frz vs refreshed {res['B6']['rel_frz_vs_refreshed']:.3e}; kin_I1 frozen {res['B6']['kin_I1_frozen_a0_on_frz_field']:.2f} "
        f"local {res['B6']['kin_I1_local_a0_on_frz_field']:.2f}; seed a0 on seed {res['B6']['kin_I1_seed_a0_on_seed']:.2f}")

    # -------- boost-sector probe (not claimed): M0i = 0 is an invariant subspace of the descent
    free32 = ~STACK.pin_shell(32, cfg["h"])
    r32 = radius(32, cfg["h"])
    E0 = energy_J(base, cfg, 1.0, 200.0)
    G0 = grad_J(base, cfg, 1.0, 200.0)
    probes = {}
    dirs = {}
    for k in range(3):
        D = np.zeros_like(base)
        blk = rng.normal(size=base.shape[:3] + (3,)) * free32[..., None]
        D[..., 0, 1:] = blk; D[..., 1:, 0] = blk
        dirs[f"random_{k}"] = D / np.sqrt(np.sum(D * D))
    for nm, rc in (("core_bump", 0.0), ("bulk_bump", 15.0)):
        D = np.zeros_like(base)
        w = np.exp(-((r32 - rc) / 3.0) ** 2) * free32
        if nm == "bulk_bump":
            X = np.stack(np.meshgrid(*[np.arange(32) - 15.5] * 3, indexing="ij"), -1)
            w = w * (np.abs(X[..., 1]) < 3) * (np.abs(X[..., 2]) < 3) * (X[..., 0] > 0)   # one patch on the +x axis
        D[..., 0, 1] = w; D[..., 1, 0] = w
        dirs[nm] = D / np.sqrt(np.sum(D * D))
    for nm, D in dirs.items():
        g1 = float(np.sum(G0 * D))
        row = {"grad_component": g1}
        for eps in (1e-3, 1e-2, 3e-2):
            ep = energy_J(base + eps * D, cfg, 1.0, 200.0)
            em = energy_J(base - eps * D, cfg, 1.0, 200.0)
            row[f"curv_eps_{eps:g}"] = (ep["E_J"] + em["E_J"] - 2 * E0["E_J"]) / eps ** 2
            row[f"curv_KP_eps_{eps:g}"] = (ep["E_KP"] + em["E_KP"] - 2 * E0["E_KP"]) / eps ** 2
            row[f"curv_cert_eps_{eps:g}"] = ((ep["E_u"] + ep["V4"]) + (em["E_u"] + em["V4"]) - 2 * (E0["E_u"] + E0["V4"])) / eps ** 2
        probes[nm] = row
        log(f"BOOST probe {nm}: grad.D {g1:.2e}, curvature eps 1e-3 {row['curv_eps_0.001']:.4e} (K_P part {row['curv_KP_eps_0.001']:.4e}, "
            f"L_cert part {row['curv_cert_eps_0.001']:.4e}), eps 1e-2 {row['curv_eps_0.01']:.4e}")
    # 20 own FIRE steps from a boosted field: does |M0i| grow or decay?
    Mb = base + 1e-2 * dirs["random_0"]
    Mb2, rows_b = own_fire(Mb, cfg, 1.0, 200.0, free32, 20, tag="boost_probe")
    probes["fire_from_boosted"] = {"M0i_norm_start": float(np.sqrt(np.sum(Mb[..., 0, 1:] ** 2))),
                                   "M0i_norm_end": float(np.sqrt(np.sum(Mb2[..., 0, 1:] ** 2))),
                                   "M0i_max_start": float(np.abs(Mb[..., 0, 1:]).max()), "M0i_max_end": float(np.abs(Mb2[..., 0, 1:]).max()),
                                   "E_start": rows_b[0]["E_J"], "E_end": rows_b[-1]["E_J"], "E_unboosted": E0["E_J"], "rows": rows_b}
    res["boost_probe"] = probes
    log(f"BOOST FIRE 20 steps: |M0i| {probes['fire_from_boosted']['M0i_norm_start']:.4e} -> {probes['fire_from_boosted']['M0i_norm_end']:.4e}, "
        f"E {rows_b[0]['E_J']:.4f} -> {rows_b[-1]['E_J']:.4f} (unboosted {E0['E_J']:.4f})")

    # -------- B4 scaling ratios
    f = res["fields"]
    res["B4_scaling"] = {"kinKP_L72_over_L48_h1.5": f["n48L72_c1_J200"]["B1"]["kin_KP"] / f["n32_c1_J200"]["B1"]["kin_KP"],
                         "L3_ratio": (72 / 48) ** 3,
                         "kinKP_h1_over_h1.5_L48": f["n48L48_c1_J200"]["B1"]["kin_KP"] / f["n32_c1_J200"]["B1"]["kin_KP"],
                         "omega_L72_over_L48": f["n48L72_c1_J200"]["B1"]["omega"] / f["n32_c1_J200"]["B1"]["omega"],
                         "omega_h1_over_h1.5": f["n48L48_c1_J200"]["B1"]["omega"] / f["n32_c1_J200"]["B1"]["omega"],
                         "kinKP_over_n3vac": {t: f[t]["B4"]["kinKP_total_over_n3_vac"] for t in f}}
    log(f"B4 scaling: kinKP L72/L48 {res['B4_scaling']['kinKP_L72_over_L48_h1.5']:.3f} (L^3 {res['B4_scaling']['L3_ratio']:.3f}); "
        f"h1/h1.5 {res['B4_scaling']['kinKP_h1_over_h1.5_L48']:.3f}; omega L72/L48 {res['B4_scaling']['omega_L72_over_L48']:.3f}")

    # -------- trace analysis (not claimed but found)
    ta = {}
    for tag, rec in res["traces"].items():
        tr = rec["trace"]
        its = np.array([q["it"] for q in tr], float)
        ekp = np.array([q["E_KP"] for q in tr])
        gm = np.array([q["gap_mean"] for q in tr])
        gf = np.array([q["gap_frac_below_0.05"] for q in tr])
        kk = np.array([q["kin_KP"] for q in tr])
        ki = np.array([q["kin_I1"] for q in tr])
        fr = np.array([q["kinI1_frac_r_gt_12"] for q in tr])
        fx = np.array([q["fmax"] for q in tr])
        dts = np.array([q["dt"] for q in tr])
        # rates over the last third and the first third (per iteration)
        k = max(2, len(tr) // 3)
        r_late = float(np.polyfit(its[-k:], ekp[-k:], 1)[0])
        r_early = float(np.polyfit(its[:k], ekp[:k], 1)[0])
        # power-law fit E_KP = a it^b on the second half (b from log-log)
        hh = len(tr) // 2
        b_pow = float(np.polyfit(np.log(its[hh:]), np.log(ekp[hh:]), 1)[0])
        ta[tag] = {"E_KP_first": float(ekp[0]), "E_KP_last": float(ekp[-1]), "E_KP_drop_rel": float((ekp[0] - ekp[-1]) / ekp[0]),
                   "rate_early_per_it": r_early, "rate_late_per_it": r_late, "rate_ratio_late_over_early": r_late / r_early if r_early else None,
                   "its_to_zero_at_late_rate": float(-ekp[-1] / r_late) if r_late < 0 else None,
                   "loglog_slope_second_half": b_pow,
                   "gap_mean_first": float(gm[0]), "gap_mean_last": float(gm[-1]), "gap_frac_first": float(gf[0]), "gap_frac_last": float(gf[-1]),
                   "kin_KP_first": float(kk[0]), "kin_KP_last": float(kk[-1]), "kin_KP_rel_change": float(kk[-1] / kk[0] - 1),
                   "kin_I1_first": float(ki[0]), "kin_I1_last": float(ki[-1]),
                   "kinI1_frac_first": float(fr[0]), "kinI1_frac_last": float(fr[-1]),
                   "fmax_min": float(fx.min()), "fmax_max": float(fx.max()), "fmax_last": float(fx[-1]),
                   "dt_median": float(np.median(dts)), "dt_max_seen": float(dts.max()),
                   "E_monotone_in_trace": bool(np.all(np.diff([q["E"] for q in tr]) <= 0))}
    res["trace_analysis"] = ta
    for tag, q in ta.items():
        log(f"TRACE {tag}: E_KP {q['E_KP_first']:.0f} -> {q['E_KP_last']:.0f} ({100*q['E_KP_drop_rel']:.2f}%), rate early {q['rate_early_per_it']:.3f}/it "
            f"late {q['rate_late_per_it']:.3f}/it, loglog slope {q['loglog_slope_second_half']:.3f}, gap mean {q['gap_mean_first']:.4f}->{q['gap_mean_last']:.4f}, "
            f"gap frac {q['gap_frac_first']:.4f}->{q['gap_frac_last']:.4f}, kin_KP {q['kin_KP_rel_change']:+.3f}, kinI1 frac {q['kinI1_frac_first']:.3f}->{q['kinI1_frac_last']:.3f}, "
            f"dt median {q['dt_median']:.4f}")

    res["wall_s"] = round(time.time() - T0, 1)
    os.makedirs(DATA, exist_ok=True)
    with open(os.path.join(DATA, "m5_32_r14_b_audit.json"), "w") as fh:
        json.dump(res, fh, indent=1)
    log(f"wrote {os.path.join(DATA, 'm5_32_r14_b_audit.json')}")


if __name__ == "__main__":
    main()
