"""M5.32 R26-3: the spin gate with the physical rotation generator, on every stored charge,
on a box ladder, and along a delta arm; the g-gate ingredients; the held pair's charge.

EQUATIONS FIRST
---------------
A rigid rotation of the whole field about z, M(x) -> R M(R^-1 x) R^T, has the velocity
    a_rigid = [J_z, M] - (x d_y - y d_x) M   (the internal rotation minus the transport),
the physical generator of the author's 2026-09-24 reply (it cancels on a radial uniaxial
exterior: R26-0 check g, 0.18 percent on the uniaxial seed). The three inertias
    C_int = kin(M; [J_z, M]),  C_orb = kin(M; (x d_y - y d_x) M),  C_rigid = kin(M; a_rigid),
kin(M; a0) = 4 h^3 sum_i <[a0, A_i]_eta, [a0, A_i]_eta>_eta (the certified stencil), with the
within-r profiles (r 3, 6, 9, 12, 18) that split the core from the halo. The clock numbers of
R25-2 on each: omega_* = sqrt(E / (3 C)), J_* = sqrt(4 C E / 3), E_rot / E = 1 / 4 by
construction (R25-2's gate on E_J = E + J^2 / (4 C)).
THE BOX LADDER: the stored R25-2 charges at h 1.5 (L 48, 72, 96) and h 1 (L 48, 64), each at
its END field and its GATE field (the stage-gate sensitivity = the relative change of C_rigid
between the two).
THE DELTA ARM: the S1 charge at delta 0.1 and 0.03 in the R25-2 descent (n 32, L 48, cap
8000, kick), read the same way; the law of C_rigid against delta beside the author's
"residual growth from biaxial strands scales with delta".
THE g INGREDIENTS (the author's 2026-09-24 16:08 UTC gate g = (4/3)(Omega/omega_e)
<r^2>_charge / lambda_C^2): <r^2> on the energy density and on the director's topological
charge density (the hedgehog Jacobian), omega_* on C_int and on C_rigid; omega_e and lambda_C
are not continuum quantities of this stack, so g is NOT evaluated, the ingredients are.
THE HELD PAIR (report 018 item 3): the hedgehog degree of the director on a sphere of radius
4 around each core and on the outer sphere r 20 of the stored R22-2 unlike pairs (d 6, 8,
12), a statement about where the held pair's charges sit at the end of our descent.

PRE-REGISTERED LABELS
---------------------
    SPIN_PHYS_BOX_CONVERGED   C_rigid(L 96) / C_rigid(L 48) at h 1.5 within 10 percent of 1
                              AND C_rigid(L 64) / C_rigid(L 48) at h 1 within 10 percent
    SPIN_PHYS_BOX_DEPENDENT   either ratio off by more than 10 percent
    SPIN_PHYS_INSUFFICIENT    the stage-gate sensitivity above 10 percent on any ladder row, or a
                              ladder row that never reached the gate (the R25-2 rows beyond n 32)
    the delta law, the g ingredients and the pair degrees are reported, never labeled.

Modes: reads | run_arm [workers] | collect | smoke. Output: data/m5_32_r26_3_spin.json,
the delta-arm arrays in data/m5_32_r26_3/ (local, kept). Regenerate: reads about 5 min,
run_arm 2 about 2 to 3 h.
"""

import importlib.util
import json
import multiprocessing as mp
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
OUT_JSON = os.path.join(DATA, "m5_32_r26_3_spin.json")
OUT_NPZ = os.path.join(DATA, "m5_32_r26_3")
R25_2_DIR = os.path.join(DATA, "m5_32_r25_2")
R25_2_JSON = os.path.join(DATA, "m5_32_r25_2_charge.json")
R22_2_DIR = os.path.join(DATA, "m5_32_r22_2")
T0 = time.time()
G = 8.0
W1S = 25.0
LADDER = {"1.5": [(32, 48.0), (48, 72.0), (64, 96.0)], "1": [(48, 48.0), (64, 64.0)]}
TOL_BOX = 0.10
TOL_STAGE = 0.10


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


F0 = _load("m5_32_r26_0_form", "m5_32_r26_0_form.py")
R25 = _load("m5_32_r25_2_charge", "m5_32_r25_2_charge.py")
CS, R21, R20, B3, R0, W1 = R25.CS, R25.R21, R25.R20, R25.B3, R25.R0, R25.W1
R25.OUT_NPZ = (
    OUT_NPZ  # the delta arm's arrays land here (R25's run_job writes to its module's OUT_NPZ)
)
R25.OUT_JSON = OUT_JSON.replace(".json", "_arm_rows.json")


def log(msg):
    print(f"[{time.time() - T0:8.1f}s] {msg}", flush=True)


def cfg_pot(n, L, delta, w1s=W1S):
    cfg = R21.cfg_of(n, L, G, delta)
    p = R21.params_of(G, delta)
    pot = ("v4", R0.roots_of(cfg), W1 * w1s)
    return cfg, p, pot


# ================= the reads =================
def charge_density(M, cfg):
    """the director's topological charge density (the hedgehog Jacobian, oriented outward), h^3-weighted."""
    n, h = cfg["n"], cfg["h"]
    lam, V = np.linalg.eigh(M[..., 1:, 1:])
    d = V[..., :, 2]
    X, Y, Z = B3.coords(n, h)
    r = np.sqrt(X * X + Y * Y + Z * Z)
    rhat = np.stack([X, Y, Z], -1) / np.maximum(r, 1e-12)[..., None]
    sgn = np.sign(np.einsum("...a,...a->...", d, rhat))
    d = d * np.where(sgn == 0, 1.0, sgn)[..., None]
    dx = F0.dsym(d, 0, h)
    dy = F0.dsym(d, 1, h)
    dz = F0.dsym(d, 2, h)

    # the topological current J_k = (1 / 4 pi) (1/2) eps_kij d . (d_i d x d_j d); rho = div J (for the
    # hedgehog d = r-hat, J = r-hat / (4 pi r^2), the Coulomb field of a unit charge)
    def J_comp(du, dv):
        return np.einsum("...a,...a->...", d, np.cross(du, dv)) / (4.0 * np.pi)

    Jx_, Jy_, Jz_ = J_comp(dy, dz), J_comp(dz, dx), J_comp(dx, dy)
    rho = (
        F0.dsym(Jx_[..., None, None], 0, h)[..., 0, 0]
        + F0.dsym(Jy_[..., None, None], 1, h)[..., 0, 0]
        + F0.dsym(Jz_[..., None, None], 2, h)[..., 0, 0]
    )
    return rho * h**3, r


def r2_moments(M, cfg, p, pot):
    n, h = cfg["n"], cfg["h"]
    X, Y, Z = B3.coords(n, h)
    r2 = X * X + Y * Y + Z * Z
    e = R20.density(M, cfg, pot)
    e = np.maximum(e, 0.0)
    rho, r = charge_density(M, cfg)
    pin = B3.pin_shell(n, h)
    out = {
        "r2_energy_weighted": float(np.sum(e * r2) / max(np.sum(e), 1e-300)),
        "energy_total_density_sum": float(np.sum(e)),
        "r2_charge_weighted_abs": float(
            np.sum(np.abs(rho) * r2) / max(np.sum(np.abs(rho)), 1e-300)
        ),
        "charge_total": float(np.sum(rho[~pin])),
        "charge_abs_total": float(np.sum(np.abs(rho[~pin]))),
        "rho_r2_energy_weighted_xy": float(np.sum(e * (X * X + Y * Y)) / max(np.sum(e), 1e-300)),
    }
    return out


def spin_reads(M, cfg, p, pot, delta, label):
    parts = R20.energy_parts(M, cfg, p, pot)
    E = float(parts["E_total"])
    pg = F0.physical_generator(M, cfg)
    out = {
        "label": label,
        "E": E,
        "E_curv": float(parts["E_curv"]),
        "V": float(parts["V"]),
        "generator": pg,
    }
    for key in ("internal", "rigid", "orbital"):
        C = pg[key]
        if C > 0 and E > 0:
            out[f"omega_star_{key}"] = float(np.sqrt(E / (3.0 * C)))
            out[f"J_star_{key}"] = float(np.sqrt(4.0 * C * E / 3.0))
    out["moments"] = r2_moments(M, cfg, p, pot)
    out["g_note"] = (
        "g = (4/3)(Omega/omega_e) <r^2>/lambda_C^2 NOT evaluated: omega_e and lambda_C are not continuum quantities of this stack"
    )
    return out


def ladder_reads():
    out = {}
    for hk, rows in LADDER.items():
        for n, L in rows:
            tag = f"S1_d0.3_w25_n{n}_L{L:g}"
            rec = {}
            cfg, p, pot = cfg_pot(n, L, 0.3)
            for suffix in ("", "_gate"):
                f = os.path.join(R25_2_DIR, tag + suffix + ".npz")
                if not os.path.exists(f):
                    rec[suffix or "_end"] = "absent"
                    continue
                M = np.load(f)["M"]
                rec[suffix or "_end"] = spin_reads(M, cfg, p, pot, 0.3, tag + (suffix or "_end"))
                log(
                    f"{tag}{suffix}: C_int {rec[suffix or '_end']['generator']['internal']:.1f} C_rigid {rec[suffix or '_end']['generator']['rigid']:.1f}"
                )
            fk = os.path.join(R25_2_DIR, tag + "_kick_stage.npz")
            if os.path.exists(fk):
                # RUN-TIME DEVIATION (2026-09-25): a STABLE row keeps its gate field as its end field, so
                # the end-against-gate sensitivity is zero by construction; the kicked-and-relaxed field
                # (the kick of 2 percent of E, two chunks) is the stored perturbation that prices it
                Zk = np.load(fk, allow_pickle=True)
                Mk = Zk["M"]
                rec["_kick_stage"] = spin_reads(Mk, cfg, p, pot, 0.3, tag + "_kick_stage")
                ce, ck = (
                    rec["_end"]["generator"]["rigid"],
                    rec["_kick_stage"]["generator"]["rigid"],
                )
                rec["kick_stage_sensitivity_rigid"] = abs(ce - ck) / max(abs(ce), 1e-300)
            if isinstance(rec.get("_end"), dict) and isinstance(rec.get("_gate"), dict):
                ce, cg = rec["_end"]["generator"]["rigid"], rec["_gate"]["generator"]["rigid"]
                rec["stage_gate_sensitivity_rigid"] = abs(ce - cg) / max(abs(ce), 1e-300)
                ci, cgi = (
                    rec["_end"]["generator"]["internal"],
                    rec["_gate"]["generator"]["internal"],
                )
                rec["stage_gate_sensitivity_internal"] = abs(ci - cgi) / max(abs(ci), 1e-300)
            rec["h"] = float(hk)
            rec["L"] = L
            rec["n"] = n
            out[tag] = rec
    return out


def pair_reads():
    out = {}
    if not os.path.isdir(R22_2_DIR):
        return {"absent": True}
    for f in sorted(os.listdir(R22_2_DIR)):
        if not f.startswith("unlike_") or not f.endswith(".npz") or "ckpt" in f:
            continue
        tag = f[:-4]
        try:
            parts = tag.split("_")
            d = float(parts[1][1:])
            n = int(parts[2][1:])
            L = float(parts[3][1:])
        except Exception:  # noqa: BLE001
            continue
        nf = np.load(os.path.join(R22_2_DIR, f))["n"]
        cfg, _, _ = cfg_pot(n, L, 0.3)
        Rc = min(4.0, 0.5 * d - 1.0)
        rec = {"d": d, "n": n, "L": L, "core_sphere_R": Rc}
        for xc in (+d / 2, -d / 2):
            rec[f"core_x{xc:+g}"] = F0.degree_reader(
                None, cfg, Rc, center=(xc, 0.0, 0.0), nfield=nf
            )["degree"]
        Ro = 0.5 * L - 4.0
        rec["outer_R"] = Ro
        rec["outer_degree"] = F0.degree_reader(None, cfg, Ro, nfield=nf)["degree"]
        rec["mid_degree_R_half_d"] = F0.degree_reader(None, cfg, 0.5 * d, nfield=nf)["degree"]
        out[tag] = rec
        log(
            f"pair {tag}: cores {rec[f'core_x{+d/2:+g}']:.3f} {rec[f'core_x{-d/2:+g}']:.3f} outer {rec['outer_degree']:.3f}"
        )
    return out


# ================= the delta arm =================
def jobs_arm():
    return [dict(kind="S1", delta=d, w1s=W1S, n=32, L=48.0, cap=8000) for d in (0.1, 0.03)]


def run_arm_job(j):
    row = R25.run_job(j)
    if row.get("status") == "OK":
        cfg, p, pot = cfg_pot(j["n"], j["L"], j["delta"])
        M = np.load(os.path.join(OUT_NPZ, row["tag"] + ".npz"))["M"]
        row["spin_reads"] = spin_reads(M, cfg, p, pot, j["delta"], row["tag"])
        row["own_reads"] = {
            "biaxiality": F0.biaxiality_reads(M, cfg, j["delta"]),
            "gap_tail": F0.gap_tail(M, cfg, j["delta"]),
            "partition_r9": F0.partition_reader(M, cfg, 9.0, j["delta"])["partition"],
            "partition_r18": F0.partition_reader(M, cfg, 18.0, j["delta"])["partition"],
        }
    return row


def load_json():
    if os.path.exists(OUT_JSON):
        with open(OUT_JSON) as f:
            return json.load(f)
    return {"task": "M5.32 R26-3", "ladder": {}, "arm_rows": {}, "pair": {}, "collect": {}}


def save_json(J):
    tmp = OUT_JSON + ".tmp"
    with open(tmp, "w") as f:
        json.dump(J, f, indent=1, default=str)
    os.replace(tmp, OUT_JSON)


def run_arm(workers):
    os.makedirs(OUT_NPZ, exist_ok=True)
    J = load_json()
    pending = [
        j for j in jobs_arm() if J["arm_rows"].get(R25.job_tag(j), {}).get("status") != "OK"
    ]
    log(f"arm pool: {len(pending)} jobs, {workers} workers")
    with ProcessPoolExecutor(max_workers=int(workers), mp_context=mp.get_context("spawn")) as ex:
        futs = [ex.submit(run_arm_job, j) for j in pending]
        for fut in as_completed(futs):
            row = fut.result()
            Jn = load_json()
            Jn["arm_rows"][row["tag"]] = row
            save_json(Jn)
    log("pool done")


def collect():
    J = load_json()
    lad = J.get("ladder", {})
    out = {"box": {}, "delta_law": {}, "stage_gate": {}, "gate_labels": {}}
    worst_stage = 0.0
    ratios = {}
    # RUN-TIME DEVIATION (2026-09-25): the R25-2 rows beyond n 32 never reached the gate (capped,
    # FALLING), so a box ladder on them is INSUFFICIENT by the plan's own gate rule; the ratios are
    # reported beside the label
    try:
        with open(R25_2_JSON) as f:
            r25rows = json.load(f)["rows"]
    except Exception:  # noqa: BLE001
        r25rows = {}
    all_at_gate = True
    for hk, rows in LADDER.items():
        cs = {}
        for n, L in rows:
            r = lad.get(f"S1_d0.3_w25_n{n}_L{L:g}", {})
            e = r.get("_end")
            if isinstance(e, dict):
                cs[f"{L:g}"] = {
                    "C_rigid": e["generator"]["rigid"],
                    "C_int": e["generator"]["internal"],
                    "C_orb": e["generator"]["orbital"],
                    "rigid_within_r": e["generator"]["within_r"],
                    "omega_star_rigid": e.get("omega_star_rigid"),
                    "omega_star_internal": e.get("omega_star_internal"),
                }
            sg = r.get("stage_gate_sensitivity_rigid")
            if sg is not None:
                out["stage_gate"][f"h{hk}_L{L:g}"] = sg
                worst_stage = max(worst_stage, sg)
            ks = r.get("kick_stage_sensitivity_rigid")
            if ks is not None:
                out["stage_gate"][f"h{hk}_L{L:g}_kick_stage"] = ks
                worst_stage = max(worst_stage, ks)
            gl = (r25rows.get(f"S1_d0.3_w25_n{n}_L{L:g}") or {}).get("gate_label")
            out["gate_labels"][f"h{hk}_L{L:g}"] = gl
            all_at_gate = all_at_gate and gl == "AT_GATE"
        out["box"][f"h{hk}"] = cs
        Ls = sorted(cs, key=float)
        if len(Ls) >= 2:
            ratios[hk] = cs[Ls[-1]]["C_rigid"] / cs[Ls[0]]["C_rigid"]
            out["box"][f"h{hk}_ratio_Lmax_over_L48_rigid"] = ratios[hk]
            out["box"][f"h{hk}_ratio_Lmax_over_L48_internal"] = (
                cs[Ls[-1]]["C_int"] / cs[Ls[0]]["C_int"]
            )
    if worst_stage > TOL_STAGE or not all_at_gate:
        label = "SPIN_PHYS_INSUFFICIENT"
    elif len(ratios) == 2 and all(abs(v - 1.0) <= TOL_BOX for v in ratios.values()):
        label = "SPIN_PHYS_BOX_CONVERGED"
    elif len(ratios) == 2:
        label = "SPIN_PHYS_BOX_DEPENDENT"
    else:
        label = "SPIN_PHYS_INSUFFICIENT"
    out["label"] = label
    out["worst_stage_gate_sensitivity"] = worst_stage
    # the delta law: n32 L48 at delta 0.3 (stored) plus the arm
    pts = []
    e = (lad.get("S1_d0.3_w25_n32_L48") or {}).get("_end")
    if isinstance(e, dict):
        pts.append(
            (0.3, e["generator"]["rigid"], e["generator"]["internal"], "R25-2 stored", "STABLE?")
        )
    for tag, r in J.get("arm_rows", {}).items():
        sr = r.get("spin_reads")
        if isinstance(sr, dict):
            pts.append(
                (
                    r["delta"],
                    sr["generator"]["rigid"],
                    sr["generator"]["internal"],
                    tag,
                    r.get("kick_label"),
                )
            )
    pts.sort()
    out["delta_law"]["points_delta_Crigid_Cint_tag_kick"] = pts
    if len(pts) >= 2:
        pf = np.polyfit(np.log([q[0] for q in pts]), np.log([max(q[1], 1e-300) for q in pts]), 1)
        out["delta_law"]["rigid_exponent"] = float(pf[0])
    J["collect"] = out
    save_json(J)
    return out


def smoke():
    out = {}
    cfg, p, pot = cfg_pot(32, 48.0, 0.3)
    M = np.load(os.path.join(R25_2_DIR, "S1_d0.3_w25_n32_L48.npz"))["M"]
    t = time.time()
    sr = spin_reads(M, cfg, p, pot, 0.3, "smoke")
    out["spin_reads_wall_s"] = round(time.time() - t, 1)
    out["C"] = {k: sr["generator"][k] for k in ("internal", "orbital", "rigid")}
    out["moments"] = sr["moments"]
    # the charge density integrates to the degree on the S1 seed
    Ms = R20.seed_axes(cfg, (1.0, 0.3, 0.0))
    rho, r = charge_density(Ms, cfg)
    out["seed_charge_inside_r12"] = float(np.sum(rho[r < 12.0]))
    out["seed_charge_ok"] = abs(out["seed_charge_inside_r12"] - 1.0) < 0.05

    # labels
    def syn(ratio_15, ratio_1, stage=0.01):
        lad = {}
        for hk, rows in LADDER.items():
            for k, (n, L) in enumerate(rows):
                fac = 1.0 if k == 0 else (ratio_15 if hk == "1.5" else ratio_1)
                g = {"rigid": 100.0 * fac, "internal": 1000.0, "orbital": 10.0, "within_r": {}}
                lad[f"S1_d0.3_w25_n{n}_L{L:g}"] = {
                    "_end": {"generator": g},
                    "stage_gate_sensitivity_rigid": stage,
                }
        return lad

    keep = OUT_JSON
    globals()["OUT_JSON"] = keep.replace(".json", "_smoke_rows.json")
    labs = {}
    for name, args in (("conv", (1.05, 1.02)), ("dep", (1.5, 1.0)), ("insuf", (1.0, 1.0, 0.5))):
        save_json(
            {"task": "smoke", "ladder": syn(*args), "arm_rows": {}, "pair": {}, "collect": {}}
        )
        labs[name] = collect()["label"]
    if os.path.exists(OUT_JSON):
        os.remove(OUT_JSON)
    globals()["OUT_JSON"] = keep
    out["labels"] = labs
    out["labels_ok"] = labs == {
        "conv": "SPIN_PHYS_BOX_CONVERGED",
        "dep": "SPIN_PHYS_BOX_DEPENDENT",
        "insuf": "SPIN_PHYS_INSUFFICIENT",
    }
    out["PASS"] = bool(out["seed_charge_ok"] and out["labels_ok"] and sr["generator"]["rigid"] > 0)
    with open(OUT_JSON.replace(".json", "_smoke.json"), "w") as f:
        json.dump(out, f, indent=1, default=str)
    print(json.dumps(out, indent=1, default=str))


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "smoke"
    if mode == "smoke":
        smoke()
    elif mode == "reads":
        J = load_json()
        J["ladder"] = ladder_reads()
        J["pair"] = pair_reads()
        save_json(J)
        collect()
    elif mode == "run_arm":
        run_arm(sys.argv[2] if len(sys.argv) > 2 else 2)
        collect()
    elif mode == "collect":
        print(json.dumps(collect(), indent=1, default=str))
    else:
        raise SystemExit(f"unknown mode {mode}")


if __name__ == "__main__":
    main()
