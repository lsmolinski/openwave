"""M5.32 R15-P-iv: fixed-J descent on L_P (the author's prediction (iv): the descent reaches the (1,2) sheet;
both sides expect no minimizer).

Equations first (m5_32_r15_common.py for L_P, m5_32_r13w_common.py for the clock):
    E_J[M] = E_stat[M] + J^2 / (4 kin_tot[M]),   kin_tot = kin_I1 + c_P kin_KP23,   omega = J / (2 kin_tot)
    kin_I1  = 4 h^3 sum_cells sum_i <[a0, A_i]_eta, [a0, A_i]_eta>_eta        (the certified I1 inertia)
    kin_KP23 = (1/2) h^3 sum_cells tr(Om_0^T eta Om_0 eta),  Om_0 = P23 a0 eta P23   (the K_P^23 inertia)
    a0 = a0_local(M) = J_loc M - M J_loc, J_loc the rotation generator about the local leading spatial eigenvector
         (refreshed every step, frozen in the gradient: the R13-W / R14-B protocol)
    E_stat = E_u + V4^dd + mu SPLIT + c_P K_P^23 (static),  c_P = 1, mu = 1e-2, J = 200
FIRE (dt0 0.01, dt_max 0.1), the vacuum shell pinned (depth 1.6), seed = the R15-M relaxed hedgehog of the same
(n, L) at (mu 1e-2, c_P 1), 3000 iterations on n32 L48, then n48 L72.
Diagnostics per 100 iterations: E_J, E_stat, kin_I1, kin_KP23, omega, the max and core-mean split (lambda_3 - lambda_2),
the fraction of the kinetic density (I1 + c_P K_P^23) inside r < L/4, the (1,2)-channel share of the I1 kinetic
density (the kinetic density of the (1,2)-plane rotation a0_12 = G3 M - M G3 relative to a0_local's), and the field
checkpoint (resume-complete).
Pre-registered verdicts: PERIODIC_ORBIT_EXISTS (a stationary state: fmax < 1e-3 or an E_J plateau over 2000 iterations
with the kinetic density localized, fraction inside L/4 >= 0.8); CANDIDATE_REFUTED (no stationary state: E_J still
descending at 3000 iterations, or the fixed-J term invisible, E_J - E_stat < 1e-3 E_stat, or the inertia extensive,
the kinetic fraction inside L/4 < 0.5 and growing); BLIND_BY_THEOREM (the descent moves the inertia into the (1,2)
channel where K_P^23 has no static cost, R15-H H4: the (1,2) share rising above 0.5 while E_stat does not rise).

usage: python3 m5_32_r15_p4_fixedj.py <n> <L> <J> <maxit>
       python3 m5_32_r15_p4_fixedj.py stationarity <n> <L> <J>     (the end field: directional derivatives of the true E_J)
"""
import sys
ARGS = list(sys.argv[1:])
import os, json
import numpy as np
import m5_32_r15_common as C15

INS4, C13, B8 = C15.INS4, C15.C13, C15.B8
log = C15.log
OUT = os.path.join(C15.CK, "p4_fixedj")
os.makedirs(OUT, exist_ok=True)


def kp_kin_density(M, a0, cfg):
    E, _, _ = C15.kp23_cells([a0], M, need_grad=False)
    return cfg["h"] ** 3 * E


def make_diag(cfg):
    n, h, L = cfg["n"], cfg["h"], cfg["L"]
    X, Y, Z = INS4.coords(n, h)
    r = np.sqrt(X * X + Y * Y + Z * Z)
    core = r < L / 4

    def diag(M):
        a0 = C13.a0_local(M)
        ki = C13.kin_density(M, a0, cfg)
        kk = cfg["cP"] * kp_kin_density(M, a0, cfg)
        tot = ki + kk
        a12 = B8.G3 @ M - M @ B8.G3
        k12 = C13.kin_density(M, a12, cfg)
        lam = C15.sorted_spectrum(M)
        split = lam[..., 2] - lam[..., 1]
        T = float(np.sum(tot))
        return {"kin_frac_core": float(np.sum(tot[core]) / T) if T > 0 else None,
                "kin_12_over_kin_local_I1": float(np.sum(k12) / max(float(np.sum(ki)), 1e-300)),
                "max_split": float(np.max(split)), "core_mean_split": float(np.mean(split[core]))}
    return diag


def verdict(hist, info, seed_diag, seed_parts):
    """the pre-registered rules, every comparison against the SEED (not the first logged row)."""
    last = hist[-1]
    if info["stop"] == "non-finite":
        return "RUNAWAY"
    stationary = info["stop"] in ("f_tol", "plateau") or last["fmax"] < 1e-3
    if stationary and last["kin_frac_core"] is not None and last["kin_frac_core"] >= 0.8:
        return "PERIODIC_ORBIT_EXISTS"
    invisible = (last["E_J"] - last["E_stat"]) < 1e-3 * abs(last["E_stat"])
    extensive = last["kin_frac_core"] is not None and last["kin_frac_core"] < 0.5 and last["kin_frac_core"] < seed_diag["kin_frac_core"]
    blind = last["kin_12_over_kin_local_I1"] > 0.5 and last["E_stat"] <= seed_parts["E_stat"] * 1.01 \
        and last["kin_12_over_kin_local_I1"] > seed_diag["kin_12_over_kin_local_I1"]
    if blind:
        return "BLIND_BY_THEOREM"
    if (not stationary) or invisible or extensive:
        return "CANDIDATE_REFUTED"
    return "PERIODIC_ORBIT_EXISTS (delocalized)"


def stationarity(n, L, J, nd=6, eps=1e-4, seed=0):
    """is the end field stationary for the TRUE E_J (a0 refreshed)?  Directional derivatives of E_J along
    random free symmetric directions and along the frozen-a0 gradient, by central differences."""
    mu, cP = 1e-2, 1.0
    cfg = C15.cfg_dd(n, L, mu=mu, cP=cP)
    tag = f"fixedJ_n{n}_L{L:g}_J{J:g}"
    M = np.load(os.path.join(OUT, tag + ".npy"))
    free = (~INS4.pin_shell(n, cfg["h"]))[..., None, None].astype(float)

    def EJ(Mx):
        pp = C15.lp_parts(Mx, cfg, C13.a0_local(Mx))
        return pp["E_stat"] + J * J / (4.0 * pp["kin_tot"]), pp
    e0, p0 = EJ(M)
    rng = np.random.default_rng(seed)
    rows = []
    a0 = C13.a0_local(M)
    Gf = C15.lp_grad(M, cfg) - (J * J / (4.0 * p0["kin_tot"] ** 2)) * C15.lp_kin_grad(M, a0, cfg)
    Gf = Gf * free
    dirs = [("frozen-a0 gradient (normalized)", Gf / np.sqrt(np.sum(Gf * Gf)))]
    for i in range(nd):
        D = rng.normal(size=M.shape); D = 0.5 * (D + np.swapaxes(D, -1, -2)) * free
        dirs.append((f"random {i}", D / np.sqrt(np.sum(D * D))))
    for name, D in dirs:
        ep, _ = EJ(M + eps * D); em, _ = EJ(M - eps * D)
        d1 = (ep - em) / (2 * eps); d2 = (ep - 2 * e0 + em) / (eps * eps)
        rows.append({"direction": name, "dE_J/ds": float(d1), "d2E_J/ds2": float(d2)})
        log(f"stationarity {tag} {name}: dE_J/ds {d1:+.4e}  d2E_J/ds2 {d2:+.4e}")
    rec = {"E_J": e0, "parts": p0, "frozen_gradient_norm": float(np.sqrt(np.sum(Gf * Gf))), "frozen_gradient_max": float(np.max(np.abs(Gf))), "rows": rows, "eps": eps}
    js = os.path.join(OUT, tag + ".json")
    if os.path.exists(js):
        r = json.load(open(js)); r["stationarity"] = rec
        r["verdict"] = verdict(r["trace"], {"stop": r["stop"]}, r["seed_diag"], r["seed_parts"])
        json.dump(r, open(js, "w"), indent=1)
        log(f"{tag}: verdict re-evaluated against the seed: {r['verdict']}")
    return rec


def main():
    if ARGS[0] == "stationarity":
        stationarity(int(ARGS[1]), float(ARGS[2]), float(ARGS[3]))
        return
    n, L, J, maxit = int(ARGS[0]), float(ARGS[1]), float(ARGS[2]), int(ARGS[3])
    mu, cP = 1e-2, 1.0
    cfg = C15.cfg_dd(n, L, mu=mu, cP=cP)
    seed = os.path.join(C15.CK, "m_hedgehog", f"relax_n{n}_L{L:g}_mu{mu:g}_cP{cP:g}.npy")
    M0 = np.load(seed)
    free = ~INS4.pin_shell(n, cfg["h"])
    tag = f"fixedJ_n{n}_L{L:g}_J{J:g}"
    ck = os.path.join(OUT, tag + ".npy")
    diag = make_diag(cfg)
    d0 = diag(M0)
    p0 = C15.lp_parts(M0, cfg, C13.a0_local(M0))
    log(f"{tag} seed {seed}: E_stat {p0['E_stat']:.4f} kin_I1 {p0['kin_I1']:.4e} kin_KP {p0['kin_KP_raw']:.4e} diag {d0}")
    M, info = C15.fire_lp(M0, cfg, free, maxit, J=J, a0_of=C13.a0_local, log_every=100, tag=tag, diag=diag, ck_path=ck, ck_every=500)
    np.save(ck, M)
    v = verdict(info["trace"], info, d0, p0)
    rec = {"n": n, "L": L, "h": cfg["h"], "J": J, "mu": mu, "cP": cP, "maxit": maxit, "seed": seed, "stop": info["stop"], "iters": info["iters"],
           "wall_s": info["wall_s"], "seed_parts": p0, "seed_diag": d0, "trace": info["trace"], "verdict": v, "field": ck,
           "verdict_rule": __doc__.split("Pre-registered verdicts:")[1].split("usage:")[0].strip()}
    json.dump(rec, open(os.path.join(OUT, tag + ".json"), "w"), indent=1)
    last = info["trace"][-1]
    log(f"{tag} END {v}: stop {info['stop']} it {info['iters']} E_J {last['E_J']:.4f} E_stat {last['E_stat']:.4f} kin {last['kin_tot']:.4e} "
        f"omega {last['omega']:.5f} fmax {last['fmax']:.3e} kin_frac_core {last['kin_frac_core']} 12-share {last['kin_12_over_kin_local_I1']:.3f}")


if __name__ == "__main__":
    main()
