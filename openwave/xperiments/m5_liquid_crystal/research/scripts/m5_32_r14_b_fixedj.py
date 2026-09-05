"""M5.32 R14-B: the DECISIVE NUMERIC test, the fixed-J continuation on L_cert + c K_P^h
(ledger 6.3; the R13-W W3 machinery with the local generator on the relaxed hedgehog).

EQUATIONS FIRST
---------------
E_J[M] = E_stat[M] + J^2 / (4 kin_tot[M]),
    E_stat  = E_u + V4 + c K_P^h_stat                         (INS4.e_parts + m5_32_r14_terms.kp_h_energy_grad)
    kin_tot = kin_I1(M; a0) + c kin_KPh(M; a0),  a0 = a0_local(M)  (refreshed every step, FROZEN in the gradient)
    dE_J/dM = dE_stat/dM - (J^2 / (4 kin_tot^2)) dkin_tot/dM|_{a0 frozen},   omega = J / (2 kin_tot)
K_P^h in the invariant order (R14-0 audit): (1/2) sum_i ||Om_i||_F^2 in the eta-orthonormal
eigenbasis of N, Om = P A eta P, P = (N + g)(N - 1); its exact gradient through the jets and the
eigenbasis is finite-difference gated at 4e-9 (m5_32_r14_terms selftest lines).
Protocol (R13-W verbatim): FIRE dt0 = 0.01, dt_max = 0.1, the vacuum pinned at the box edge
(pin_shell depth 1.6) only, no taper, no mask; seeds = the relaxed L_cert hedgehog of each box
(m5_32_r13w_common.seed_hedgehog); c in {0.3, 1, 3} at J = 200 on n32 L48 (h 1.5), J in {50, 800}
and the FROZEN-generator control at c = 1; the h x L ladder at c = 1, J = 200: n48 L72 (h 1.5),
n48 L48 (h 1.0), n64 L48 (h 0.75, seed relaxed in this run); 3000 iterations (the packet's 6000
cut to the day: stated) or fmax < 1e-6.
Verdicts (frozen): PERIODIC_ORBIT_EXISTS if E_J and omega converge on the ladder (the two
finest rungs in each direction within 10 percent, drift shrinking); CANDIDATE_REFUTED with the
runaway exponent otherwise (drift, runaway, melt).
PRE-REGISTERED EXPECTATION (written before any number, from R14-0): K_P^h's static cost on the
relaxed hedgehog is 9.5e4 against E_u = 9 (the hedgehog's (2,3)-frame connection is charged), so
the descent should first close the (2,3) gap (degenerate everywhere except the pinned shell,
at a V4 cost of 6e-6 per cell), which switches K_P^h's own inertia off (proportional to
(lambda_2 - lambda_3)^2), leaving kin_I1 alone; then the R13-W free-inertia runaway (the (1,2)
twist, invisible to K_P^h) is expected: CANDIDATE_REFUTED, with the gap closure as the new
mechanism.  The author's prediction (09-05) is the same verdict by the twist route.
Diagnostics per trace row: E, E_u, V4, E_KP, kin_I1, kin_KP, omega, fmax, gap23 (min, mean over
the free cells), the kin_I1 fraction outside r = 12, max |M_0i|.

Stages: run <tag> (one job), batch <n32|ladder>, collect.  Local end fields under
checkpoints/m5_32_r14/b_fixedj/ (gitignored); results data/m5_32_r14_b_fixedj.json,
plots/m5_32_r14_b_fixedj.png.
"""
from __future__ import annotations
import importlib.util
import json
import os
import sys
import time
from multiprocessing import Pool
import numpy as np

ARGV = list(sys.argv)
HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.dirname(HERE)
DATA, PLOTS = os.path.join(RES, "data"), os.path.join(RES, "plots")
CKB = os.path.join(RES, "checkpoints", "m5_32_r14", "b_fixedj")
os.makedirs(CKB, exist_ok=True)
OUT = os.path.join(DATA, "m5_32_r14_b_fixedj.json")
PNG = os.path.join(PLOTS, "m5_32_r14_b_fixedj.png")
T0 = time.time()


def log(m):
    print(f"[{time.time() - T0:8.1f}s] {m}", flush=True)


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


T14 = _load("m5_32_r14_terms", "m5_32_r14_terms.py")
B3 = T14.B3
C13 = _load("m5_32_r13w_common", "m5_32_r13w_common.py")
INS4 = C13.INS4

JOBS = {
    # R14-B' (2026-09-05 overnight): the boost sector was never sampled (the audit: M_0i = 0 is an
    # invariant subspace of the descent, a saddle at c = 0.3); seeds dressed with the R3 radial
    # boost field at rapidity 0.1 (rc = 1) so M_0i != 0 from the start
    # R14-B3 (overnight): the P250 structure imposed by hand on L_cert + K_P^h: the exterior (r > 12)
    # projected to the degenerate (2,3) pair at the seed (iota = 0 there, no ticking), the core kept;
    # does the action keep the exterior at rest (V4 costs 6e-6 per cell, the K_P^h tail 9.5e4 is gone)
    # and does a localized fixed-J clock survive inside?
    "n32_c1_J200_degext": (32, 48.0, 1.0, 200.0, "degext"),
    # the same with the generator FROZEN at the degenerate seed (a0 = 0 outside r = 12, the core's clock
    # kept): the true gradient flow of the reported E_J (the R13-W J3 caveat), the clean test of a
    # localized clock inside an exterior at rest
    "n32_c1_J200_degext_frz": (32, 48.0, 1.0, 200.0, "degext_frz"),
    # overnight controls: the c = 1 base run continued to the packet's 6000 iterations from its
    # 3000-iteration end field; the boost-seeded run under the PLAIN (indefinite) K_P, the one
    # configuration where the H-adjoint correction is not inert
    "n32_c1_J200_cont6000": (32, 48.0, 1.0, 200.0, "continue"),
    # the g = 32 control of the packet (the stiff-potential toy point): seed = the R10 gcore relaxation
    "n32_c1_J200_g32": (32, 48.0, 1.0, 200.0, "g32"),
    "n32_c1_J200_boost_plain": (32, 48.0, 1.0, 200.0, "boost_plain"),
    "n32_c0.3_J200_boost": (32, 48.0, 0.3, 200.0, "boost"),
    "n32_c1_J200_boost": (32, 48.0, 1.0, 200.0, "boost"),
    "n32_c0.3_J200": (32, 48.0, 0.3, 200.0, False),
    "n32_c1_J200": (32, 48.0, 1.0, 200.0, False),
    "n32_c3_J200": (32, 48.0, 3.0, 200.0, False),
    "n32_c1_J50": (32, 48.0, 1.0, 50.0, False),
    "n32_c1_J800": (32, 48.0, 1.0, 800.0, False),
    "n32_c1_J200_frz": (32, 48.0, 1.0, 200.0, True),
    "n48L72_c1_J200": (48, 72.0, 1.0, 200.0, False),
    "n48L48_c1_J200": (48, 48.0, 1.0, 200.0, False),
    "n64L48_c1_J200": (64, 48.0, 1.0, 200.0, False),
}


def diag_of(M, cfg, a0, r):
    gap = C13.gap23(M)
    free = ~INS4.pin_shell(cfg["n"], cfg["h"])
    kd = C13.kin_density(M, a0, cfg)
    tot = float(np.sum(kd))
    return {"gap_min": float(np.min(gap[free])), "gap_mean": float(np.mean(gap[free])),
            "gap_frac_below_0.05": float(np.mean(gap[free] < 0.05)),
            "kinI1_frac_r_gt_12": float(np.sum(kd[r > 12.0]) / max(tot, 1e-300)),
            "max_abs_M0i": float(np.max(np.abs(M[..., 0, 1:])))}


def fire_kph(M0, cfg, free_mask, max_iter, c, J, frozen=False, log_every=100, tag="", f_tol=1e-6,
             plateau=(2000, 1e-10), dt0=0.01, dt_max=0.1, plain=False):
    roots = T14.roots_for(cfg["g"])
    free = free_mask[..., None, None].astype(float)
    M = M0.copy()
    v = np.zeros_like(M)
    dt, alpha, n_up = dt0, 0.1, 0
    hist = []
    X, Y, Z = B3.coords(cfg["n"], cfg["h"])
    r = np.sqrt(X * X + Y * Y + Z * Z)
    a0_frozen = C13.a0_local(M0) if frozen else None

    def parts(Mx):
        a0 = a0_frozen if frozen else C13.a0_local(Mx)
        e_u, e_v = INS4.e_parts(Mx, cfg)
        if plain:
            Ek, Gk = T14.kp_energy_grad(Mx, cfg, roots)
            kk, Gkk = T14.kp_plain_kin_grad(Mx, a0, cfg, roots)
        else:
            Ek, Gk, kk, Gkk = T14.kp_h_energy_grad(Mx, cfg, roots, a0)
        ki = float(INS4.kin_of(Mx, a0, cfg))
        ktot = ki + c * kk
        E = float(e_u + e_v) + c * Ek + J * J / (4.0 * ktot)
        return E, float(e_u), float(e_v), c * Ek, ki, c * kk, a0, (Gk, Gkk)

    def tot_grad(Mx):
        E, e_u, e_v, Ek, ki, kk, a0, (Gk, Gkk) = parts(Mx)
        ktot = ki + kk
        G = INS4.grad(Mx, cfg) + c * Gk
        G = G - (J * J / (4.0 * ktot * ktot)) * (INS4.kin_grad(Mx, a0, cfg) + c * Gkk)
        return G, (E, e_u, e_v, Ek, ki, kk, a0)

    G, info = tot_grad(M)
    F = -G * free
    t0 = time.time()
    stop = "max_iter"
    it = 0
    for it in range(1, max_iter + 1):
        P = float(np.sum(F * v))
        if P > 0.0:
            n_up += 1
            vn = np.sqrt(np.sum(v * v)); fn = np.sqrt(np.sum(F * F))
            v = (1 - alpha) * v + alpha * (F / max(fn, 1e-300)) * vn
            if n_up > 5:
                dt = min(dt * 1.1, dt_max); alpha *= 0.99
        else:
            v[:] = 0.0; dt *= 0.5; alpha = 0.1; n_up = 0
        v += dt * F
        M += dt * v
        G, info = tot_grad(M)
        F = -G * free
        fmax = float(np.max(np.abs(F)))
        if not np.isfinite(fmax):
            stop = "non-finite"; break
        if it % log_every == 0 or it == max_iter:
            E, e_u, e_v, Ek, ki, kk, a0 = info
            row = {"it": it, "E": E, "E_u": e_u, "V4": e_v, "E_KP": Ek, "kin_I1": ki, "kin_KP": kk,
                   "omega": J / (2.0 * (ki + kk)), "fmax": fmax, "dt": dt}
            row.update(diag_of(M, cfg, a0, r))
            row["E_boost_probe"] = float(np.sum(M[..., 0, 1:] ** 2))
            hist.append(row)
            print(f"  {tag} it {it:5d} E {E:12.4f} E_u {e_u:8.3f} V4 {e_v:.3e} E_KP {Ek:11.3f} kinI1 {ki:9.2f} kinKP {kk:11.1f} "
                  f"om {row['omega']:.3e} gap {row['gap_min']:.3f}/{row['gap_mean']:.3f} fmax {fmax:.2e} [{time.time() - t0:.0f}s]", flush=True)
            back = max(1, plateau[0] // max(log_every, 1))
            if len(hist) > back and abs(E - hist[-1 - back]["E"]) < plateau[1]:
                stop = "plateau"; break
        if fmax < f_tol:
            stop = "f_tol"; break
    return M, {"stop": stop, "trace": hist, "wall_s": round(time.time() - t0, 1), "iters": it}


def run_job(tag, maxit=3000):
    n, L, c, J, frozen = JOBS[tag]
    outj = os.path.join(CKB, tag + ".json")
    if os.path.exists(outj):
        return json.load(open(outj))
    M0, cfg, seedrec = C13.seed_hedgehog(n, L)
    plain = False
    fire_kw = {}
    if frozen == "g32":
        M0 = np.load(os.path.join(RES, "checkpoints", "m5_32_r10", "relax_g32_n32_L48_it6000_g32.npy"))
        cfg = C13.cfg_of(32, 48.0, g=32.0)
        seedrec = {"source": "relax_g32_n32_L48_it6000_g32.npy (R10 gcore, 6000 it)", "g": 32.0, "fire": "dt0 1.5e-4, dt_max 1.5e-3 (the R10 gcore steps: the stiff toy point)"}
        fire_kw = {"dt0": 1.5e-4, "dt_max": 1.5e-3}
        frozen = False
    if frozen == "continue":
        M0 = np.load(os.path.join(CKB, "n32_c1_J200.npy"))
        seedrec = {"source": "n32_c1_J200 end field (3000 it)", "continued": True}
        frozen = False
    if frozen == "boost_plain":
        plain = True
        frozen = "boost"
    if frozen in ("degext", "degext_frz"):
        X, Y, Z = B3.coords(n, cfg["h"])
        rr = np.sqrt(X * X + Y * Y + Z * Z)
        M0 = C13.degenerate_project(M0, rr > 12.0)
        seedrec = dict(seedrec); seedrec["degenerate_exterior"] = "cells with r > 12 projected to the (2,3) trace part at the seed"
        seedrec["gap_mean_seed"] = float(np.mean(C13.gap23(M0)))
        frozen = (frozen == "degext_frz")
    if frozen == "boost":
        R3 = _load("m5_32_r3_i_ansatz", "m5_32_r3_i_ansatz.py")
        X, Y, Z = B3.coords(n, cfg["h"])
        Lb = R3.boost_field(X, Y, Z, (0.0, 0.0, 0.0), 0.1, 1.0)
        M0 = Lb @ M0 @ Lb.swapaxes(-1, -2)
        seedrec = dict(seedrec); seedrec["boost_dressing"] = "R3 radial boost, rapidity 0.1, rc 1"
        seedrec["max_abs_M0i_seed"] = float(np.max(np.abs(M0[..., 0, 1:])))
        frozen = False
    free = ~INS4.pin_shell(n, cfg["h"])
    M, info = fire_kph(M0, cfg, free, maxit, c, J, frozen=frozen, tag=tag, plain=plain, **fire_kw)
    rec = {"tag": tag, "n": n, "L": L, "h": cfg["h"], "c": c, "J": J, "frozen_generator": frozen, "plain_K_P": plain, "maxit": maxit,
           "seed": seedrec, "stop": info["stop"], "iters": info["iters"], "wall_s": info["wall_s"], "trace": info["trace"]}
    if np.all(np.isfinite(M)):
        np.save(os.path.join(CKB, tag + ".npy"), M)
    json.dump(rec, open(outj, "w"), indent=1, default=float)
    if info["trace"]:
        log(f"DONE {tag}: stop {info['stop']} E {info['trace'][-1]['E']:.4f} omega {info['trace'][-1]['omega']:.3e} gap_mean {info['trace'][-1]['gap_mean']:.3f} wall {info['wall_s']}")
    else:
        log(f"DONE {tag}: stop {info['stop']} before the first log row (iters {info['iters']}) wall {info['wall_s']}")
    return rec


def _job(args):
    return run_job(*args)


def batch(which, workers):
    tags = [t for t in JOBS if (t.startswith("n32") if which == "n32" else not t.startswith("n32"))]
    todo = [(t, 3000) for t in tags if not os.path.exists(os.path.join(CKB, t + ".json"))]
    log(f"batch {which}: {len(todo)} jobs on {workers} workers")
    with Pool(workers) as pool:
        for _ in pool.imap_unordered(_job, todo):
            pass


def collect():
    recs = {t: json.load(open(os.path.join(CKB, t + ".json"))) for t in JOBS if os.path.exists(os.path.join(CKB, t + ".json"))}
    res = {"rung": "R14-B", "runs": {}}
    for t, r in recs.items():
        tr = r["trace"]
        e0, e1 = tr[0], tr[-1]
        # runaway exponent: kin_tot vs iteration on the last half of the trace
        its = np.array([x["it"] for x in tr]); kt = np.array([x["kin_I1"] + x["kin_KP"] for x in tr])
        half = its >= its[-1] / 2
        p = float(np.polyfit(np.log(its[half]), np.log(np.maximum(kt[half], 1e-300)), 1)[0]) if np.sum(half) >= 3 else None
        res["runs"][t] = {"n": r["n"], "L": r["L"], "h": r["h"], "c": r["c"], "J": r["J"], "frozen": r["frozen_generator"],
                          "stop": r["stop"], "iters": r["iters"], "wall_s": r["wall_s"],
                          "E_start": e0["E"], "E_end": e1["E"], "E_KP_start": e0["E_KP"], "E_KP_end": e1["E_KP"],
                          "kin_I1_start": e0["kin_I1"], "kin_I1_end": e1["kin_I1"], "kin_KP_start": e0["kin_KP"], "kin_KP_end": e1["kin_KP"],
                          "omega_start": e0["omega"], "omega_end": e1["omega"],
                          "gap_mean_start": e0["gap_mean"], "gap_mean_end": e1["gap_mean"], "gap_frac_below_0.05_end": e1["gap_frac_below_0.05"],
                          "kinI1_frac_r_gt_12_end": e1["kinI1_frac_r_gt_12"], "V4_end": e1["V4"], "E_u_end": e1["E_u"],
                          "kin_tot_iteration_exponent_last_half": p,
                          "omega_last_quarter_rel_drift": float((tr[-1]["omega"] - tr[int(0.75 * (len(tr) - 1))]["omega"]) / max(abs(tr[-1]["omega"]), 1e-300))}
    # ladder verdict at c = 1, J = 200
    lad = {k: v for k, v in res["runs"].items() if v["c"] == 1.0 and v["J"] == 200.0 and not v["frozen"]}
    verdict = "UNDECIDED (ladder incomplete)"
    if "n32_c1_J200" in lad:
        base = lad["n32_c1_J200"]
        stationary = base["stop"] == "f_tol" or abs(base["omega_last_quarter_rel_drift"]) < 0.1
        h_pairs = [("n32_c1_J200", "n48L48_c1_J200"), ("n48L48_c1_J200", "n64L48_c1_J200")]
        L_pairs = [("n32_c1_J200", "n48L72_c1_J200")]
        conv = []
        for a, b in h_pairs + L_pairs:
            if a in lad and b in lad:
                conv.append(abs(lad[a]["omega_end"] - lad[b]["omega_end"]) / max(abs(lad[a]["omega_end"]), 1e-300))
        res["ladder_omega_rel_diffs"] = conv
        if stationary and conv and all(x < 0.1 for x in conv):
            verdict = "PERIODIC_ORBIT_EXISTS (at the relaxed-field level, on the reduced ladder)"
        else:
            verdict = "CANDIDATE_REFUTED (no stationary fixed-J state / no ladder convergence)"
    res["verdict"] = verdict
    res["mechanism"] = "see gap_mean_end (the (2,3) closure), E_KP_end vs start (the term switching itself off), kinI1 fraction outside r = 12 (the free-inertia sheet)"
    json.dump(res, open(OUT, "w"), indent=1, default=float)
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 4, figsize=(18, 4))
        for t, r in recs.items():
            tr = r["trace"]; its = [x["it"] for x in tr]
            ax[0].plot(its, [x["E"] for x in tr], label=t)
            ax[1].semilogy(its, [x["kin_I1"] + x["kin_KP"] for x in tr], label=t)
            ax[2].semilogy(its, [x["omega"] for x in tr], label=t)
            ax[3].plot(its, [x["gap_mean"] for x in tr], label=t)
        ax[0].set_title("E_J"); ax[1].set_title("kin_tot"); ax[2].set_title("omega = J/(2 kin_tot)"); ax[3].set_title("mean (2,3) gap, free cells")
        for a in ax:
            a.set_xlabel("iteration"); a.legend(fontsize=6)
        plt.tight_layout(); plt.savefig(PNG, dpi=110)
    except Exception as e:                                        # noqa: BLE001
        log(f"plot skipped: {e!r}")
    log(f"collect: {verdict} -> {OUT}")
    return res


if __name__ == "__main__":
    st = ARGV[1] if len(ARGV) > 1 else "collect"
    if st == "run":
        run_job(ARGV[2], int(ARGV[3]) if len(ARGV) > 3 else 3000)
    elif st == "batch":
        batch(ARGV[2], int(ARGV[3]) if len(ARGV) > 3 else 6)
    else:
        collect()
