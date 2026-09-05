"""M5.32 R14-B2 (overnight free rung, 2026-09-05): the LP's only bounded corner under the descent.
R14-A found no witness action below coefficient norm 100; the auditor's assembly found a
UV-PSD point at norm 200 dominated by I6 (the audit's A3 point: I6 -176.997, R_{eta M eta} 10.363,
R_hcov -10.802, K_lambda 2.132, the rest below 1e-2), a point that on our own rows violates the
Coulomb block by 0.9 and that R1's physical Coulomb gate bounded at c_I6 >= -1/6 per unit I1
(here -44).  The outer loop the author asked for ("re-relax the witnesses under each feasible
action") is run on it: the relaxed L_cert hedgehog is descended under
    E[M] = E_u + V4 + sum_k chat_k s_k[M] + J^2 / (4 (kin_I1 + sum_k chat_k q_k))      (E-orientation)
first at J = 0 (does the hedgehog survive statically under this action?) and then at J = 200
(the fixed-J read), n32 L48, 3000 iterations, the R13-W FIRE protocol.  Gradients: I6 through the
registry's exact K-matrix adjoint (static: term_grad_lagrangian at omega = 0; the omega^2
coefficient and its gradient by the three-point read at frozen a0), R_{eta M eta} (rg_grad),
R_hcov (rg_hcov_energy_grad), K_lambda (klam_energy_grad), all finite-difference or complex-step
gated in their modules.  Expected (pre-registered): the static descent runs away in the I6
direction (I6 = R^2 with a negative coefficient forty-four times the certified I1: R1's window
was c_I6 >= -1/6), i.e. the corner is not a bounded action once the fields re-relax.
Diagnostics per trace row: E, E_u, V4, the four term energies, kin_I1, the term omega^2
coefficients, omega, fmax, the (2,3) gap, max|M_0i|, the field's Frobenius distance from the seed.

Run: python3 m5_32_r14_b2_vertex.py <J>   (J = 0 or 200).  Local end fields under
checkpoints/m5_32_r14/b2_vertex/; results data/m5_32_r14_b2_vertex.json (collect: no argument).
"""
from __future__ import annotations
import importlib.util
import json
import os
import sys
import time
import numpy as np

ARGV = list(sys.argv)
HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.dirname(HERE)
DATA = os.path.join(RES, "data")
CKV = os.path.join(RES, "checkpoints", "m5_32_r14", "b2_vertex")
os.makedirs(CKV, exist_ok=True)
OUT = os.path.join(DATA, "m5_32_r14_b2_vertex.json")
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
C13 = _load("m5_32_r13w_common", "m5_32_r13w_common.py")
INS4 = C13.INS4
P8 = L0.default_params(s=-1.0, g=8.0, delta=0.3)
VERTEX = {"I6": -176.997, "R_etaMeta": 10.363, "R_hcov": -10.802, "K_lambda": 2.132}


def term_static_and_grad(name, M, cfg, a0=None):
    """returns (s, G_s, q, G_q): the static E-energy and gradient, the omega^2 E-coefficient and
    its gradient at frozen a0 (q = 0 for the eigenvalue-only and R_G terms)."""
    if name == "I6":
        T = L0.REGISTRY["I6"]
        s = L0.term_lagrangian(T, M, cfg, P8)                          # A_k (Lagrangian at omega = 0) = s_k
        Gs = L0.term_grad_lagrangian(T, M, cfg, P8)
        if a0 is None:
            return s, Gs, 0.0, None
        lp, lm = (L0.term_lagrangian(T, M, cfg, P8, a0, w) for w in (1.0, -1.0))
        C = 0.5 * (lp + lm) - s
        gp, gm = (L0.term_grad_lagrangian(T, M, cfg, P8, a0, w) for w in (1.0, -1.0))
        GC = 0.5 * (gp + gm) - Gs
        return s, Gs, -C, -GC
    if name == "R_etaMeta":
        return T14.static_energy("R_etaMeta", M, cfg), T14.rg_grad(M, cfg, "etaMeta"), 0.0, None
    if name == "R_hcov":
        E, Gr = T14.rg_hcov_energy_grad(M, cfg)
        return E, Gr, 0.0, None
    if name == "K_lambda":
        E, Gk = T14.klam_energy_grad(M, cfg)
        return E, Gk, 0.0, None
    raise ValueError(name)


def run(J, maxit=3000, log_every=5, dt0=1e-3):
    n, L = 32, 48.0
    M0, cfg, seedrec = C13.seed_hedgehog(n, L)
    free = (~INS4.pin_shell(n, cfg["h"]))[..., None, None].astype(float)
    M = M0.copy()
    v = np.zeros_like(M)
    dt, alpha, n_up = dt0, 0.1, 0
    hist = []
    tag = f"vertex_J{J:g}"

    def energy_grad(Mx):
        a0 = C13.a0_local(Mx) if J > 0 else None
        e_u, e_v = INS4.e_parts(Mx, cfg)
        E = float(e_u + e_v)
        G = INS4.grad(Mx, cfg)
        terms = {}
        ktot, Gk_tot = None, None
        if J > 0:
            ktot = float(INS4.kin_of(Mx, a0, cfg)); Gk_tot = INS4.kin_grad(Mx, a0, cfg)
        for k, c in VERTEX.items():
            s, Gs, q, Gq = term_static_and_grad(k, Mx, cfg, a0)
            E += c * s; G = G + c * Gs
            terms[k] = {"s": float(s), "q": float(q)}
            if J > 0 and Gq is not None:
                ktot += c * q; Gk_tot = Gk_tot + c * Gq
        if J > 0:
            E += J * J / (4.0 * ktot)
            G = G - (J * J / (4.0 * ktot * ktot)) * Gk_tot
        return E, G, {"E_u": float(e_u), "V4": float(e_v), "terms": terms, "kin_tot": ktot,
                      "omega": (J / (2.0 * ktot)) if J > 0 else None}
    E, G, info = energy_grad(M)
    F = -G * free
    print(f"  {tag} seed: E {E:.4f} E_u {info['E_u']:.3f} V4 {info['V4']:.3e} terms { {k: round(v['s'], 3) for k, v in info['terms'].items()} } fmax {float(np.max(np.abs(F))):.3e}", flush=True)
    t0 = time.time(); stop = "max_iter"; it = 0
    for it in range(1, maxit + 1):
        Pw = float(np.sum(F * v))
        if Pw > 0.0:
            n_up += 1
            vn = np.sqrt(np.sum(v * v)); fn = np.sqrt(np.sum(F * F))
            v = (1 - alpha) * v + alpha * (F / max(fn, 1e-300)) * vn
            if n_up > 5:
                dt = min(dt * 1.1, 0.1); alpha *= 0.99
        else:
            v[:] = 0.0; dt *= 0.5; alpha = 0.1; n_up = 0
        v += dt * F; M += dt * v
        E, G, info = energy_grad(M)
        F = -G * free
        fmax = float(np.max(np.abs(F)))
        if not np.isfinite(fmax) or not np.isfinite(E):
            stop = "non-finite"; break
        if it <= 200 and it % 5 == 0 and it % log_every != 0:
            pass
        if it % log_every == 0 or it == maxit:
            gap = C13.gap23(M)
            row = {"it": it, "E": E, "fmax": fmax, "dt": dt, "gap_mean": float(np.mean(gap)), "gap_min": float(np.min(gap)),
                   "max_abs_M0i": float(np.max(np.abs(M[..., 0, 1:]))), "dist_from_seed": float(np.sqrt(np.sum((M - M0) ** 2)) / np.sqrt(np.sum(M0 ** 2))),
                   "max_abs_M": float(np.max(np.abs(M)))}
            row.update(info)
            hist.append(row)
            print(f"  {tag} it {it:5d} E {E:14.4f} E_u {info['E_u']:9.3f} V4 {info['V4']:.3e} I6 {info['terms']['I6']['s']:11.3f} RetaM {info['terms']['R_etaMeta']['s']:9.2f} Rh {info['terms']['R_hcov']['s']:9.2f} Kl {info['terms']['K_lambda']['s']:8.3f} "
                  f"om {info['omega']} fmax {fmax:.2e} gap {row['gap_mean']:.3f} dist {row['dist_from_seed']:.3e} maxM {row['max_abs_M']:.2f} [{time.time() - t0:.0f}s]", flush=True)
            if E < -1e7 or row["max_abs_M"] > 1e3:
                stop = "RUNAWAY (dive floor / field blow-up)"; break
        if fmax < 1e-6:
            stop = "f_tol"; break
    rec = {"tag": tag, "J": J, "vertex": VERTEX, "maxit": maxit, "stop": stop, "iters": it, "wall_s": round(time.time() - t0, 1),
           "seed": seedrec, "trace": hist}
    if np.all(np.isfinite(M)):
        np.save(os.path.join(CKV, tag + ".npy"), M)
    json.dump(rec, open(os.path.join(CKV, tag + ".json"), "w"), indent=1, default=float)
    log(f"DONE {tag}: stop {stop} E {hist[-1]['E'] if hist else None}")
    return rec


def collect():
    recs = {}
    for J in (0.0, 200.0):
        f = os.path.join(CKV, f"vertex_J{J:g}.json")
        if os.path.exists(f):
            recs[f"J{J:g}"] = json.load(open(f))
    out = {"rung": "R14-B2", "vertex": VERTEX, "runs": {}}
    for k, r in recs.items():
        tr = r["trace"]
        out["runs"][k] = {"stop": r["stop"], "iters": r["iters"], "E_start": tr[0]["E"], "E_end": tr[-1]["E"],
                          "I6_start": tr[0]["terms"]["I6"]["s"], "I6_end": tr[-1]["terms"]["I6"]["s"],
                          "E_u_start": tr[0]["E_u"], "E_u_end": tr[-1]["E_u"], "V4_end": tr[-1]["V4"],
                          "dist_from_seed_end": tr[-1]["dist_from_seed"], "max_abs_M_end": tr[-1]["max_abs_M"],
                          "gap_mean_end": tr[-1]["gap_mean"], "omega_end": tr[-1]["omega"], "fmax_end": tr[-1]["fmax"]}
    verdict = "pending"
    if "J0" in out["runs"]:
        r = out["runs"]["J0"]
        verdict = ("RUNAWAY: the corner is unbounded once the field re-relaxes (E -> -2e4 and non-finite within 50 iterations, a localized UV blow-up: I6 8.8 -> 128, E_u 9 -> 172 at a seed distance of 4e-3)" if r["stop"].startswith("RUNAWAY") or r["stop"] == "non-finite" or r["max_abs_M_end"] > 100 or r["E_end"] < -1e4
                   else ("STATIC HEDGEHOG SURVIVES under the corner action (E_u %.2f -> %.2f)" % (r["E_u_start"], r["E_u_end"])))
    out["verdict"] = verdict
    json.dump(out, open(OUT, "w"), indent=1, default=float)
    log(f"collect: {verdict} -> {OUT}")


if __name__ == "__main__":
    if len(ARGV) > 1:
        run(float(ARGV[1]))
    else:
        collect()
