"""M5.32 R14-C: the scalar-channel Newton route (ledger 6.3; the author's 09-03 rungs 1 and 2,
with the 09-05 amendments: run on both lambda = 0 (L_cert) and lambda = 1 (I1^h)).

EQUATIONS FIRST
---------------
Arm R_G (relaxation, the R3 arm-ii instrument):  E_c[M] = E_lambda[M] + c_R h^3 sum_cells R_G(A_i, M),
  R_G = sum_{ij} G_cd [(A_i)^{jc} (A_j)^{id} - (A_i)^{ic} (A_j)^{jd}] (static, spatial slots),
  G = eta M eta (the one M-dependent G with an exact gradient here; h_cov: linear response only,
  in the R14-A rows), E_lambda = 4 h^3 sum [(1 - lambda) I1 + lambda I1_h] + V4 (m5_32_r2_b_bounded.energy_grad).
  Seeds: the R3 arm-ii RELAXED same-sign dressed pairs at d = 10, 18 and the relaxed single at
  the same lambda (data/m5_32_r3_ii/), so the R_G term's effect is the shift from the certified
  relaxed state; heals of 1500 accepted FIRE steps (it_cap 3000), the R3 protocol verbatim
  (m5_32_r3_ii_pair.descend, with RB.energy_grad wrapped to add c_R R_G and its gradient).
  Reads: E_int(d) = E(pair, d) - 2 E(single); the SIGN of E_int(18) - E_int(10) against sign(c_R)
  (attraction = E_int increases with d); the dressing amplitude trend; the R_G share of E.
  The author's prediction: the sign of dE/dd follows sign(c_R); the pair law d^-2 (charge-dipole).
Arm K_lambda (linear response, the light-scale exchange):  on the relaxed hedgehog the eigenvalue
  deficit dl_a(x) = lambda_a(x) - lambda_a^vac (sorted spectrum of N = M eta) is the scalar
  field of the eigenvalue channel.  Multipoles: Q0_a = sum dl_a h^3 (monopole), D_a = sum dl_a x h^3
  (dipole), the far-field exponent of the shell average dl_a(r) ~ r^-p.  Linearized exchange for
  L = (1/2) c (d phi)^2 - (1/2) m_s^2 phi^2 + J phi with phi := dl_a: the source J = (-c nabla^2 +
  m_s^2) dl_a (spectral), and the interaction of two cores at separation d along z is
  E_int(d) = - sum_x J(x) dl_a(x - d e_z) h^3 (like scalar sources ATTRACT: E_int < 0, rising with d,
  range 1/m_s).  m_s in {0.1, 0.3, 1}, c = 1, d in {6, 9, 12, 18, 24}: the sign and range read.
  Also the trace read the author asked for: sum_a dl_a(r) through the core (is tr N constant?).

Stages: relax (24 heals: lambda in {0, 1} x c_R in {-1, -0.3, 0.3, 1} x {d = 10, d = 18, single}; 6 workers, ~2 h), klambda (~1 min), collect.  Local end fields under
checkpoints/m5_32_r14/c_relax/ (gitignored); results data/m5_32_r14_c_newton.json.
"""
from __future__ import annotations
import importlib.util
import json
import os
import sys
import time
from multiprocessing import Pool
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.dirname(HERE)
DATA = os.path.join(RES, "data")
CKC = os.path.join(RES, "checkpoints", "m5_32_r14", "c_relax")
os.makedirs(CKC, exist_ok=True)
OUT = os.path.join(DATA, "m5_32_r14_c_newton.json")
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
R3II = _load("m5_32_r3_ii_pair", "m5_32_r3_ii_pair.py")
RB = R3II.RB
RB_BASE = RB.energy_grad          # the UNWRAPPED lambda-family energy/gradient, captured once at import
ETA = T14.ETA
JOBS = [(lam, cR, kind, d) for lam in (0.0, 1.0) for cR in (-1.0, 1.0)
        for (kind, d) in (("same", 10), ("same", 18), ("single", 0))]


def rg_static_energy(M, cfg):
    return T14.static_energy("R_etaMeta", M, cfg)


def wrapped_energy_grad(cR):
    # AUDIT 2026-09-05 (R14-C, F2): wrapping RB.energy_grad in place let a worker that runs several
    # jobs stack the wrappers (c_eff = c1 + c2 + ...); 18 of 24 first-run heals were relaxed under
    # accumulated coefficients. Always wrap the pristine base.
    base = RB_BASE

    def eg(M, cfg, lam):
        E, G, info = base(M, cfg, lam)
        if G is None:
            return E, G, info
        Er = rg_static_energy(M, cfg)
        Gr = T14.rg_grad(M, cfg, "etaMeta")
        info = dict(info); info["E_RG"] = float(Er)
        return E + cR * Er, G + cR * Gr, info
    return eg


def run_job(args):
    lam, cR, kind, d = args
    t0 = time.time()
    cfg = R3II.cfg_of(32, 48.0)
    tag = f"lam{lam:g}_cR{cR:+g}_{kind}_d{d}"
    seed = os.path.join(DATA, "m5_32_r3_ii", f"lam{lam:g}_dr1_{kind}_d{d}_n32.npz")
    M0 = np.load(seed)["M"]
    RB.energy_grad = wrapped_energy_grad(cR)          # process-local patch (fork)
    R3II.RB = RB
    E0, _, info0 = RB.energy_grad(M0, cfg, lam)
    M, des = R3II.descend(M0, cfg, lam, R3II.STEPS_ACC, R3II.IT_CAP, tag, log_every=100)
    E1, _, info1 = RB.energy_grad(M, cfg, lam)
    # the coefficient actually applied, recomputed from the pristine base (the audit's x_eff check)
    Eb, _, _ = RB_BASE(M0, cfg, lam)
    x_eff = (E0 - Eb) / info0["E_RG"] if info0.get("E_RG") else None
    row = {"lam": lam, "cR": cR, "x_eff_check": x_eff, "kind": kind, "d": d, "tag": tag, "seed": os.path.basename(seed),
           "E_seed_total": float(E0), "E_RG_seed": info0.get("E_RG"), "E_end_total": float(E1) if np.isfinite(E1) else None,
           "E_RG_end": info1.get("E_RG") if isinstance(info1, dict) else None,
           "descent": {k: v for k, v in des.items() if k != "trace"}, "trace_tail": des["trace"][-3:],
           "wall_s": round(time.time() - t0, 1)}
    try:
        row["end_reads"] = R3II.block_reads(M, cfg, lam)
        row["seed_amp"] = R3II.dressing_amplitude(M0, cfg, kind, d)
        row["end_amp"] = R3II.dressing_amplitude(M, cfg, kind, d)
    except Exception as e:                                        # noqa: BLE001
        row["reads_error"] = repr(e)
    if np.all(np.isfinite(M)):
        np.save(os.path.join(CKC, tag + ".npy"), M)
    json.dump(row, open(os.path.join(CKC, tag + ".json"), "w"), indent=1, default=float)
    log(f"DONE {tag} stop {des['stop']} E {E0:.3f} -> {row['E_end_total']} wall {row['wall_s']}")
    return row


def stage_relax(workers=6):
    todo = [j for j in JOBS if not os.path.exists(os.path.join(CKC, f"lam{j[0]:g}_cR{j[1]:+g}_{j[2]}_d{j[3]}.json"))]
    log(f"relax: {len(todo)} jobs on {workers} workers")
    with Pool(workers) as pool:
        for row in pool.imap_unordered(run_job, todo):
            pass


def stage_klambda():
    out = {}
    for (n, L, tag) in ((32, 48.0, "n32_L48_it12000"), (48, 72.0, "n48_L72_it3000")):
        if n == 32:
            M = np.load(C13.R10_SEED); cfg = C13.cfg_of(32, 48.0)
        else:
            M, cfg, _ = C13.seed_hedgehog(48, 72.0)
        h = cfg["h"]
        lam = T14.sorted_spectrum(M)                                  # (1, delta, 0, -g) descending
        vac = np.array([1.0, 0.3, 0.0, -8.0])
        dl = lam - vac
        X, Y, Z = B3.coords(n, h)
        r = np.sqrt(X * X + Y * Y + Z * Z)
        rec = {"n": n, "L": L, "h": h, "channels": {}}
        edges = np.arange(0.0, L / 2 + h, h)
        idx = np.digitize(r, edges) - 1
        cnt = np.bincount(idx.ravel(), minlength=len(edges))[: len(edges) - 1]
        for a, nm in enumerate(("lambda_1(=1)", "lambda_2(=delta)", "lambda_3(=0)", "lambda_t(=-g)")):
            f = dl[..., a]
            Q0 = float(np.sum(f) * h ** 3)
            D = [float(np.sum(f * c) * h ** 3) for c in (X, Y, Z)]
            Qxx = float(np.sum(f * (3 * Z * Z - r * r)) * h ** 3)
            prof = np.bincount(idx.ravel(), weights=f.ravel(), minlength=len(edges))[: len(edges) - 1] / np.maximum(cnt, 1)
            rr = (np.arange(len(prof)) + 0.5) * h
            sel = (rr >= 6.0) & (rr <= L / 2 - 2 * h) & (np.abs(prof) > 1e-14)
            p = float(np.polyfit(np.log(rr[sel]), np.log(np.abs(prof[sel])), 1)[0]) if np.sum(sel) >= 3 else None
            ch = {"monopole_Q0": Q0, "dipole": D, "quadrupole_zz": Qxx, "shell_profile": [float(v) for v in prof],
                  "shell_r": [float(v) for v in rr], "far_exponent_r6_to_edge": p,
                  "max_abs_deficit": float(np.max(np.abs(f))), "core_value": float(f[n // 2, n // 2, n // 2])}
            # linearized Yukawa exchange (spectral), c = 1
            k = 2 * np.pi * np.fft.fftfreq(n, d=h)
            KX, KY, KZ = np.meshgrid(k, k, k, indexing="ij")
            k2 = KX * KX + KY * KY + KZ * KZ
            fh = np.fft.fftn(f)
            ex = {}
            for ms in (0.1, 0.3, 1.0):
                Jh = (k2 + ms * ms) * fh                              # J = (-c nabla^2 + m^2) phi
                J = np.real(np.fft.ifftn(Jh))
                vals = {}
                for d in (6.0, 9.0, 12.0, 18.0, 24.0):
                    sh = int(round(d / h))
                    phi2 = np.roll(f, sh, axis=2)                     # the partner's field displaced along z (periodic wrap)
                    vals[f"d{d:g}"] = float(-np.sum(J * phi2) * h ** 3)
                ex[f"m_s={ms:g}"] = vals
            ch["yukawa_exchange_c1"] = ex
            rec["channels"][nm] = ch
        # the trace read: sum_a dl_a(r)
        tr = np.sum(dl, axis=-1)
        prof = np.bincount(idx.ravel(), weights=tr.ravel(), minlength=len(edges))[: len(edges) - 1] / np.maximum(cnt, 1)
        rec["trace_deficit_shell_profile"] = [float(v) for v in prof]
        rec["trace_deficit_core"] = float(tr[n // 2, n // 2, n // 2])
        rec["trace_deficit_max_abs"] = float(np.max(np.abs(tr)))
        out[tag] = rec
        log(f"  {tag}: monopoles {[round(rec['channels'][c]['monopole_Q0'], 4) for c in rec['channels']]} "
            f"trace core {rec['trace_deficit_core']:.4f} max {rec['trace_deficit_max_abs']:.4f}")
    return out


def stage_collect():
    res = {"rung": "R14-C", "arms": {}}
    rows = []
    for j in JOBS:
        f = os.path.join(CKC, f"lam{j[0]:g}_cR{j[1]:+g}_{j[2]}_d{j[3]}.json")
        if os.path.exists(f):
            rows.append(json.load(open(f)))
    byk = {(r["lam"], r["cR"], r["kind"], r["d"]): r for r in rows}
    reads = {}
    for lam in (0.0, 1.0):
        for cR in (-1.0, 1.0):
            s = byk.get((lam, cR, "single", 0))
            p10, p18 = byk.get((lam, cR, "same", 10)), byk.get((lam, cR, "same", 18))
            if not (s and p10 and p18) or None in (s["E_end_total"], p10["E_end_total"], p18["E_end_total"]):
                continue
            Ei10 = p10["E_end_total"] - 2 * s["E_end_total"]
            Ei18 = p18["E_end_total"] - 2 * s["E_end_total"]
            # the seed (certified relaxed, c_R = 0) reference: the R_G energy evaluated on the seeds
            seedi10 = p10["E_seed_total"] - 2 * s["E_seed_total"]
            seedi18 = p18["E_seed_total"] - 2 * s["E_seed_total"]
            reads[f"lam{lam:g}_cR{cR:+g}"] = {
                "E_int_d10": Ei10, "E_int_d18": Ei18, "dE_int_18_minus_10": Ei18 - Ei10,
                "sign_follows_cR": bool(np.sign(Ei18 - Ei10) == np.sign(cR)),
                "attractive": bool(Ei18 - Ei10 > 0),
                "seed_linear_response_dE": seedi18 - seedi10,
                "E_RG_share_end": {f"{kd}_d{dd}": (byk[(lam, cR, kd, dd)]["E_RG_end"] / byk[(lam, cR, kd, dd)]["E_end_total"]
                                       if byk[(lam, cR, kd, dd)]["E_RG_end"] is not None else None)
                                   for (kd, dd) in (("single", 0), ("same", 10), ("same", 18))},
                "stops": {f"{kd}_d{dd}": byk[(lam, cR, kd, dd)]["descent"]["stop"] for (kd, dd) in (("single", 0), ("same", 10), ("same", 18))},
                "amp_end_over_seed": {f"{kd}_d{dd}": (byk[(lam, cR, kd, dd)].get("end_amp", {}).get("grid_max_norm_M0i", None),
                                                      byk[(lam, cR, kd, dd)].get("seed_amp", {}).get("grid_max_norm_M0i", None))
                                      for (kd, dd) in (("single", 0), ("same", 10), ("same", 18))}}
    res["arms"]["R_G_relaxed_pairs"] = {"jobs_done": len(rows), "jobs_total": len(JOBS), "reads": reads, "rows": rows}
    res["arms"]["K_lambda_exchange"] = stage_klambda()
    res["wall_s"] = round(time.time() - T0, 1)
    json.dump(res, open(OUT, "w"), indent=1, default=float)
    log(f"collect -> {OUT}")
    return res


if __name__ == "__main__":
    st = sys.argv[1] if len(sys.argv) > 1 else "collect"
    if st == "relax":
        stage_relax(int(sys.argv[2]) if len(sys.argv) > 2 else 6)
    elif st == "klambda":
        json.dump(stage_klambda(), open(os.path.join(CKC, "klambda.json"), "w"), indent=1, default=float)
    else:
        stage_collect()
