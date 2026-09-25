"""M5.32 R26-4: the index-partition census of the charge (the author's lepton-hierarchy
proposal of 2026-09-24): one branch per partition of the transverse index 2 (4 half-units)
into carriers on the boundary sphere, relaxed in the static sector, kicked, read.

EQUATIONS FIRST
---------------
The stack and the descent are R25-2's, consumed read-only (m5_32_r25_2_charge: the R23-1
reduced L-BFGS-B at c = 0 with M_00 slaved per cell, the pinned shell of depth 1.6, the gate
fmax_spatial < 1e-4 with the last-chunk drop under 1e-5 abs(E), the post-gate kick sized to
KICK_EXCESS of E, the STABLE / SADDLE / UNRESOLVED / FALLING labels). The energy
    E = 4 h^3 sum_{i<j} <F_ij, F_ij>_eta + V4,  the certified biaxial vacuum diag(8, 1, delta, 0),
    delta 0.3, w = W1 x 25, n 32, L 48 (h 1.5): the R25-2 box of the S1 row.
The seeds (R26-0.seed_partition): the radial hedgehog with the transverse pair (delta, 0) in
the frame (phi-hat, theta-hat) rotated by the pair phase
    psi = sum_k m_k arg(w - w_k),   w = (x + i y) / (r + z),
which moves the frame's index (1 at each pole) onto the chosen carriers; the pair melted along
each carrier's ray with the BPS core of its index (R25-1), the core melted as R20.seed_axes.
The five partitions of 4 half-units: {4} (south pole), {3,1} (south 3, north 1), {2,2} (the S1
seed of R25-2 with the strand melt), {2,1,1} (south 2, two halves on the equator at +-x),
{1,1,1,1} (four halves on the equator at +-x, +-y). The pinned shell HOLDS the seed's
partition on the boundary, so the census compares boundary-conditioned objects; the interior
can still reorganize (R25-2's {2,2} did) and the interior reads say how.

THE READS per branch (R26-0's readers, validated there)
------------------------------------------------------
    the seed's partition on r 9 and r 18 BEFORE relaxation (the reader gate; a seed that
        does not read its intended partition is SEED_GATE_FAIL and is not run);
    the energy at the gate with the last-chunk drop as its error, the kick label;
    the partition on the spheres r 3 to 21 at the end (readable spheres only: the reader
        needs the director outward and both gaps open on a ring near each pole), the
        carriers' distance from the polar axis (report 018 item 2, qualitative);
    beta^2 (vacuum, the block mean and the cell mean inside r 3, 6, 9);
    the gap tail (l = 0, 1, 2 coefficients per shell and their slopes on r 6 to 18);
    the physical rotation generator (C_int, C_orb, C_rigid and the within-r profiles);
    the virial ratio E_u / (3 V4) under the pinned shell;
    R25-2's own reads (windings on the axis loops, tube T(z), the spin gate on [J_z, M]).
Collect: the branch table; the barrier between every pair of branches as the maximum of the
energy along the straight line between their end fields (11 points), an UPPER BOUND on the
barrier of any path; the distinctness of branches by energy (beyond the summed errors) and by
the partition on the largest readable sphere; the labels.

PRE-REGISTERED LABELS
---------------------
    BRANCH_COUNT k       k = the number of rows at the gate (AT_GATE), STABLE or SADDLE-then-
                         reconverged under the kick, pairwise distinct (energy beyond the summed
                         errors, or a different partition on the largest readable sphere)
    CENSUS_INSUFFICIENT  fewer than 4 of the 5 rows at the gate
    MOMENT_CHANNEL_CLOSED / OPEN / UNREAD   the l = 1 gap-tail slope on every certified
                         branch (the relaxed charge keeps no z symmetry, R26-0 audit): CLOSED if every readable slope is below -3,
                         OPEN if any is -2.5 or shallower with the l = 1 coefficient above 2e-2
                         of l = 0 in the window, UNREAD otherwise; the z-symmetric branches
                         are the null controls
    beta^2, the virial and the barrier are reported, never labeled.

Modes: smoke | run [workers] | collect | jobs | reread.
Output: data/m5_32_r26_4_census.json, arrays in data/m5_32_r26_4/ (local, kept).
Regenerate: run 5 about 3 h (the R25-2 n 32 row took 2.5 h); collect seconds; smoke a minute.
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
OUT_JSON = os.path.join(DATA, "m5_32_r26_4_census.json")
OUT_NPZ = os.path.join(DATA, "m5_32_r26_4")
T0 = time.time()
G = 8.0
W1S = 25.0
DELTA = 0.3
CAP = {16: 4, 32: 8000, 48: 8000, 64: 6000}
SPHERES = (3.0, 4.5, 6.0, 9.0, 12.0, 15.0, 18.0, 21.0)
BARRIER_POINTS = 11
Z_ASYM = ("3_1", "2_1_1")


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


F0 = _load("m5_32_r26_0_form", "m5_32_r26_0_form.py")
R25 = _load("m5_32_r25_2_charge", "m5_32_r25_2_charge.py")
CS, R21, R20, B3, R0, W1 = R25.CS, R25.R21, R25.R20, R25.B3, R25.R0, R25.W1


def log(msg):
    print(f"[{time.time() - T0:8.1f}s] {msg}", flush=True)


def jobs_main(n=32, L=48.0):
    return [dict(part=k, delta=DELTA, w1s=W1S, n=n, L=L, cap=CAP[n]) for k in F0.PARTITIONS]


def job_tag(j):
    return f"P{j['part']}_d{j['delta']:g}_w{j['w1s']:g}_n{j['n']}_L{j['L']:g}"


def cfg_pot(j):
    cfg = R21.cfg_of(j["n"], j["L"], G, j["delta"])
    p = R21.params_of(G, j["delta"])
    pot = ("v4", R0.roots_of(cfg), W1 * j["w1s"])
    return cfg, p, pot


def seed_gate(M, cfg, delta, intended, radii=(9.0, 18.0)):
    out = {"intended": intended, "reads": {}}
    ok = True
    for R in radii:
        if R > 0.5 * cfg["L"] - 1.6:
            continue
        rd = F0.partition_reader(M, cfg, R, delta)
        out["reads"][f"R{R:g}"] = {"partition": rd["partition"], "total": rd["total_half_units"]}
        ok = ok and rd["partition"] == intended
    out["PASS"] = bool(ok)
    return out


def own_reads(M, cfg, p, pot, delta, part):
    """the R26-0 readers on an end field."""
    out = {"spheres": {}}
    readable = []
    for R in SPHERES:
        if R > 0.5 * cfg["L"] - 1.6:
            continue
        rd = F0.partition_reader(M, cfg, R, delta)
        rec = {
            "partition": rd["partition"],
            "total": rd["total_half_units"],
            "carriers_hu_rho_z_phi": [
                (c["half_units"], c["rho_from_axis"], c["z"], round(c["phi"], 3))
                for c in rd["carriers"]
            ],
            "unreadable_plaquettes": rd.get("unreadable_plaquettes"),
            "unreadable_by_reason": rd.get("unreadable_samples_by_reason"),
            "cap_rho_N_S": rd.get("cap_rho_N_S"),
            "min_align_director_outward": rd.get("min_align_director_outward"),
            "in_pin": rd["in_pin"],
        }
        out["spheres"][f"R{R:g}"] = rec
        if rd["partition"] is not None and not rd["in_pin"]:
            readable.append((R, rd["partition"]))
    out["largest_readable_sphere"] = readable[-1][0] if readable else None
    out["interior_partition"] = readable[-1][1] if readable else None
    out["r9_partition"] = out["spheres"].get("R9", {}).get("partition")
    out["biaxiality"] = F0.biaxiality_reads(M, cfg, delta)
    out["gap_tail"] = F0.gap_tail(M, cfg, delta)
    out["physical_generator"] = F0.physical_generator(M, cfg)
    out["virial"] = F0.virial_reads(M, cfg, p, pot)
    out["z_asymmetric_seed"] = part in Z_ASYM
    return out


def run_job(j):
    t0 = time.time()
    tag = job_tag(j)
    ug = CS.SL.ups_guard()
    cfg, p, pot = cfg_pot(j)
    mask = R21.free_mask(cfg, True)
    os.makedirs(OUT_NPZ, exist_ok=True)
    stage = os.path.join(OUT_NPZ, tag + "_stage.npz")
    row = dict(j, tag=tag, h=cfg["h"], roots=list(pot[1]), W1_eff=pot[2], gate=R25.GATE)
    try:
        if os.path.exists(stage):
            Zs = np.load(stage, allow_pickle=True)
            M, done, chunks = Zs["M"], int(Zs["done"]), json.loads(str(Zs["chunks"]))
            row["resumed_at"] = done
            _, intended = F0.partition_of(j["part"])
            row["seed_gate"] = {"intended": intended, "resumed": True, "PASS": True}
        else:
            M, intended = F0.seed_partition(cfg, j["delta"], j["part"], w=pot[2])
            if np.abs(M[..., 0, 1:]).max() > 0:
                raise RuntimeError("the start field is not block-diagonal")
            row["seed_gate"] = seed_gate(M, cfg, j["delta"], intended)
            if not row["seed_gate"]["PASS"]:
                row.update(status="SEED_GATE_FAIL", E=None, kick_label="NOT_RUN")
                row["wall_s"] = round(time.time() - t0, 1)
                return row
            np.savez_compressed(os.path.join(OUT_NPZ, tag + "_seed.npz"), M=M.astype(np.float64))
            done = 0
            chunks = [dict(R25.chunk_light(M, cfg, p, pot, mask, 0), note="the start field")]
            np.savez_compressed(stage, M=M, done=0, chunks=json.dumps(chunks))
        row["E_seed"] = chunks[0]["E"]
        M, done, chunks, verdict = R25.descend_chunks(
            M, cfg, p, pot, mask, tag, j["cap"], chunks, done, stage
        )
        row.update(chunks=chunks, iters=done, gate_label=verdict)
        E_gate = chunks[-1]["E"]
        row["E_gate"] = E_gate
        row["E_err"] = abs(chunks[-1]["drop"]) if len(chunks) > 1 else None
        np.savez_compressed(os.path.join(OUT_NPZ, tag + "_gate.npz"), M=M.astype(np.float64))
        if ug is not None and ug.wrap_up():
            row["stop"] = "UPS wrap-up before the kick (resumable)"
        if verdict == "AT_GATE":
            Mk, seed, amp_eff, E_ref = R25.kick_field(
                M, mask, pot, tag, energy=lambda X: R25.spatial_fmax(X, cfg, p, pot, mask)[0]
            )
            E_k0 = R25.spatial_fmax(Mk, cfg, p, pot, mask)[0]
            kchunks = [dict(R25.chunk_light(Mk, cfg, p, pot, mask, 0), note="the kicked field")]
            kstage = os.path.join(OUT_NPZ, tag + "_kick_stage.npz")
            Mk, kdone, kchunks, kverdict = R25.descend_chunks(
                Mk, cfg, p, pot, mask, tag + "_kick", 10**9, kchunks, 0, kstage, R25.KICK_CHUNKS
            )
            E_k = kchunks[-1]["E"]
            row["kick"] = {
                "amp_effective": amp_eff,
                "E_ref_reslaved": E_ref,
                "seed": seed,
                "E_kicked_start": E_k0,
                "E_after": E_k,
                "iters": kdone,
                "chunks": kchunks,
                "verdict": kverdict,
                "fmax_after": kchunks[-1]["fmax_spatial"],
            }
            if E_k < E_gate - R25.KICK_LOWER:
                row["kick_label"] = "SADDLE"
                row["E"] = E_k
                M = Mk
                row["reconverged"] = bool(kchunks[-1]["fmax_spatial"] < R25.GATE)
                row["E_err"] = abs(kchunks[-1]["drop"])
            elif E_k - E_gate < R25.KICK_RETURN * max(1.0, abs(E_gate)):
                row["kick_label"] = "STABLE"
                row["E"] = E_gate
                row["reconverged"] = True
            else:
                row["kick_label"] = "UNRESOLVED"
                row["E"] = E_gate
                row["reconverged"] = False
        else:
            row["kick_label"] = "FALLING"
            row["E"] = E_gate
            row["reconverged"] = False
        row["end_reads"] = R25.reads_row(M, cfg, p, pot)
        row["own_reads"] = own_reads(M, cfg, p, pot, j["delta"], j["part"])
        np.savez_compressed(os.path.join(OUT_NPZ, tag + ".npz"), M=M.astype(np.float64))
        row["status"] = "OK"
    except Exception as e:  # noqa: BLE001
        import traceback

        row.update(status="FAILED", stop=repr(e), traceback=traceback.format_exc(), E=None)
    row["wall_s"] = round(time.time() - t0, 1)
    log(
        f"DONE {tag} status {row['status']} gate {row.get('gate_label')} kick {row.get('kick_label')}"
        f" E {row.get('E')} interior {(row.get('own_reads') or {}).get('interior_partition')}"
        f" wall {row['wall_s']}"
    )
    return row


def load_json():
    if os.path.exists(OUT_JSON):
        with open(OUT_JSON) as f:
            return json.load(f)
    return {"task": "M5.32 R26-4", "rows": {}, "collect": {}}


def save_json(J):
    tmp = OUT_JSON + ".tmp"
    with open(tmp, "w") as f:
        json.dump(J, f, indent=1, default=str)
    os.replace(tmp, OUT_JSON)


def run_pool(jobs, workers):
    workers = min(int(workers), 12)
    rows = load_json()["rows"]
    pending = [j for j in jobs if rows.get(job_tag(j), {}).get("status") != "OK"]
    log(f"pool: {len(pending)} jobs, {workers} workers")
    with ProcessPoolExecutor(max_workers=workers, mp_context=mp.get_context("spawn")) as ex:
        futs = [ex.submit(run_job, j) for j in pending]
        for fut in as_completed(futs):
            row = fut.result()
            Jn = load_json()
            Jn["rows"][row["tag"]] = row
            save_json(Jn)
    log("pool done")


# ================= collect =================
def barrier(Ma, Mb, cfg, p, pot, npts=BARRIER_POINTS):
    lam = np.linspace(0.0, 1.0, npts)
    Es = []
    for t in lam:
        Mt = (1.0 - t) * Ma + t * Mb
        Es.append(float(R20.energy_parts(Mt, cfg, p, pot)["E_total"]))
    Es = np.array(Es)
    return {
        "E_path": Es.tolist(),
        "barrier_upper_bound_from_a": float(Es.max() - Es[0]),
        "barrier_upper_bound_from_b": float(Es.max() - Es[-1]),
        "argmax_lambda": float(lam[int(np.argmax(Es))]),
    }


def certified(r):
    return (
        r.get("status") == "OK"
        and r.get("gate_label") == "AT_GATE"
        and (
            r.get("kick_label") == "STABLE"
            or (r.get("kick_label") == "SADDLE" and r.get("reconverged"))
        )
    )


def collect(rows=None, with_barriers=True):
    J = load_json()
    rows = J["rows"] if rows is None else rows
    table = {}
    for tag, r in rows.items():
        own = r.get("own_reads") or {}
        gt = own.get("gap_tail") or {}
        pg = own.get("physical_generator") or {}
        bx = own.get("biaxiality") or {}
        table[tag] = {
            "part": r.get("part"),
            "status": r.get("status"),
            "gate_label": r.get("gate_label"),
            "kick_label": r.get("kick_label"),
            "reconverged": r.get("reconverged"),
            "iters": r.get("iters"),
            "E_seed": r.get("E_seed"),
            "E": r.get("E"),
            "E_err": r.get("E_err"),
            "certified": certified(r),
            "seed_partition_r18": ((r.get("seed_gate") or {}).get("reads") or {})
            .get("R18", {})
            .get("partition"),
            "interior_partition": own.get("interior_partition"),
            "largest_readable_sphere": own.get("largest_readable_sphere"),
            "r9_partition": own.get("r9_partition"),
            "beta2_vacuum": bx.get("vacuum"),
            "beta2_block_mean_r6": bx.get("block_mean_r6"),
            "beta2_cell_mean_r6": bx.get("cell_mean_r6"),
            "virial_Eu_over_3V": (own.get("virial") or {}).get("E_u_over_3V"),
            "C_int": pg.get("internal"),
            "C_orb": pg.get("orbital"),
            "C_rigid": pg.get("rigid"),
            "slope_l0": gt.get("slope_l0"),
            "slope_l1": gt.get("slope_l1"),
            "slope_l2": gt.get("slope_l2"),
            "l1_over_l0_max": gt.get("l1_over_l0_max_in_window"),
            "z_asymmetric_seed": own.get("z_asymmetric_seed"),
            "wall_s": r.get("wall_s"),
        }
    cert = {t: v for t, v in table.items() if v["certified"]}
    at_gate = sum(1 for v in table.values() if v["gate_label"] == "AT_GATE")
    # distinctness
    tags = sorted(cert)
    distinct_pairs, same_pairs = [], []
    for i in range(len(tags)):
        for k in range(i + 1, len(tags)):
            a, b = table[tags[i]], table[tags[k]]
            err = (a["E_err"] or 0.0) + (b["E_err"] or 0.0)
            by_E = abs((a["E"] or 0.0) - (b["E"] or 0.0)) > err
            by_part = (a["interior_partition"] != b["interior_partition"]) and a[
                "interior_partition"
            ] is not None
            (distinct_pairs if (by_E or by_part) else same_pairs).append(
                {
                    "a": tags[i],
                    "b": tags[k],
                    "dE": abs((a["E"] or 0.0) - (b["E"] or 0.0)),
                    "err_sum": err,
                    "by_energy": by_E,
                    "by_partition": by_part,
                }
            )
    # the count: the largest set of pairwise-distinct certified rows (greedy by energy)
    classes = []
    for t in sorted(tags, key=lambda t: table[t]["E"]):
        placed = False
        for c in classes:
            if any(sp["a"] in (t, u) and sp["b"] in (t, u) for u in c for sp in same_pairs):
                c.append(t)
                placed = True
                break
        if not placed:
            classes.append([t])
    k = len(classes)
    if at_gate < 4:
        label = "CENSUS_INSUFFICIENT"
    else:
        label = f"BRANCH_COUNT {k}"
    # the moment channel on the z-asymmetric certified branches
    # AUDIT CORRECTION (R26-0 audit, 2026-09-25): the relaxed charge is not z-reflection symmetric
    # (the descent does not keep the seed's symmetry), so EVERY certified branch carries the
    # moment read; the seed's z-symmetry is reported as data
    slopes = [
        (t, v["slope_l1"], v["l1_over_l0_max"])
        for t, v in table.items()
        if v["certified"] and v["slope_l1"] is not None
    ]
    if not slopes:
        moment = "MOMENT_CHANNEL_UNREAD"
    elif any(s >= -2.5 and (l1 or 0) > 2e-2 for _, s, l1 in slopes):
        moment = "MOMENT_CHANNEL_OPEN"
    elif all(s < -3.0 for _, s, _ in slopes):
        moment = "MOMENT_CHANNEL_CLOSED"
    else:
        moment = "MOMENT_CHANNEL_UNREAD"
    out = {
        "table": table,
        "n_rows": len(table),
        "n_at_gate": at_gate,
        "n_certified": len(cert),
        "distinct_pairs": distinct_pairs,
        "same_pairs": same_pairs,
        "classes": classes,
        "label": label,
        "moment_label": moment,
        "moment_slopes": slopes,
        "l1_over_l0_by_seed_symmetry": {
            t: (
                v["l1_over_l0_max"],
                "z-asymmetric seed" if v["z_asymmetric_seed"] else "z-symmetric seed",
            )
            for t, v in table.items()
        },
    }
    if with_barriers and len(tags) >= 2:
        bars = {}
        fields = {}
        for t in tags:
            f = os.path.join(OUT_NPZ, t + ".npz")
            if os.path.exists(f):
                fields[t] = np.load(f)["M"]
        for i in range(len(tags)):
            for k2 in range(i + 1, len(tags)):
                a, b = tags[i], tags[k2]
                if a in fields and b in fields:
                    r = rows[a]
                    cfg, p, pot = cfg_pot(dict(n=r["n"], L=r["L"], delta=r["delta"], w1s=r["w1s"]))
                    bars[f"{a}|{b}"] = barrier(fields[a], fields[b], cfg, p, pot)
        out["barriers_upper_bound"] = bars
    J["collect"] = out
    save_json(J)
    return out


def smoke():
    """the pipeline end to end at n 16 (cap 4 iterations), the labels on synthetic rows."""
    global OUT_JSON, OUT_NPZ
    keep_json, keep_npz = OUT_JSON, OUT_NPZ
    OUT_JSON = keep_json.replace(".json", "_smoke_rows.json")
    OUT_NPZ = keep_npz + "_smoke"
    out = {}
    t = time.time()
    rows = {}
    for j in jobs_main(16, 24.0)[:2]:
        r = run_job(j)
        rows[r["tag"]] = r
    out["n16_rows"] = {
        k: {
            "status": r["status"],
            "gate": r.get("gate_label"),
            "kick": r.get("kick_label"),
            "seed_gate": r.get("seed_gate"),
            "E": r.get("E"),
            "interior": (r.get("own_reads") or {}).get("interior_partition"),
            "stop": r.get("stop"),
        }
        for k, r in rows.items()
    }
    out["n16_ok"] = all(r["status"] in ("OK", "SEED_GATE_FAIL") for r in rows.values())
    out["n16_wall_s"] = round(time.time() - t, 1)

    # synthetic label tests
    def syn(E, part, kick="STABLE", gate="AT_GATE", err=1e-6, asym=False, s1=None, l1=None):
        return {
            "part": part,
            "status": "OK",
            "gate_label": gate,
            "kick_label": kick,
            "reconverged": True,
            "E": E,
            "E_err": err,
            "own_reads": {
                "interior_partition": [2, 2],
                "largest_readable_sphere": 18.0,
                "gap_tail": {"slope_l1": s1, "l1_over_l0_max_in_window": l1},
                "z_asymmetric_seed": asym,
                "biaxiality": {},
                "physical_generator": {},
                "virial": {},
            },
            "seed_gate": {"reads": {}},
        }

    five = {
        f"P{k}": syn(10.0 + 0.1 * i, k, asym=(k in Z_ASYM), s1=-4.0, l1=0.01)
        for i, k in enumerate(F0.PARTITIONS)
    }
    c5 = collect(five, with_barriers=False)
    two_same = dict(five)
    two_same["P4"]["E"] = two_same["P3_1"]["E"]
    c4 = collect(two_same, with_barriers=False)
    three = {k: v for k, v in five.items() if k not in ("P4", "P3_1")}
    c3 = collect(three, with_barriers=False)
    opn = {k: dict(v) for k, v in five.items()}
    for k in opn:
        opn[k]["own_reads"] = dict(opn[k]["own_reads"])
        opn[k]["own_reads"]["gap_tail"] = {"slope_l1": -2.0, "l1_over_l0_max_in_window": 0.05}
    co = collect(opn, with_barriers=False)
    out["labels"] = {
        "five_distinct": c5["label"],
        "two_same": c4["label"],
        "three_rows": c3["label"],
        "closed": c5["moment_label"],
        "open": co["moment_label"],
    }
    out["labels_ok"] = (
        c5["label"] == "BRANCH_COUNT 5"
        and c4["label"] == "BRANCH_COUNT 4"
        and c3["label"] == "CENSUS_INSUFFICIENT"
        and c5["moment_label"] == "MOMENT_CHANNEL_CLOSED"
        and co["moment_label"] == "MOMENT_CHANNEL_OPEN"
    )
    tags = [job_tag(j) for j in jobs_main()]
    out["wiring"] = {"jobs": len(tags), "tags_unique": len(set(tags)) == len(tags)}
    out["PASS"] = bool(out["n16_ok"] and out["labels_ok"] and out["wiring"]["tags_unique"])
    if os.path.exists(OUT_JSON):
        os.remove(OUT_JSON)
    OUT_JSON, OUT_NPZ = keep_json, keep_npz
    with open(OUT_JSON.replace(".json", "_smoke.json"), "w") as f:
        json.dump(out, f, indent=1, default=str)
    print(json.dumps({k: v for k, v in out.items() if k != "n16_rows"}, indent=1, default=str))
    print(json.dumps(out["n16_rows"], indent=1, default=str)[:3000])
    return out


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "smoke"
    if mode == "smoke":
        smoke()
    elif mode == "run":
        run_pool(jobs_main(), sys.argv[2] if len(sys.argv) > 2 else 5)
    elif mode == "collect":
        c = collect()
        print(
            json.dumps(
                {k: v for k, v in c.items() if k not in ("table", "barriers_upper_bound")},
                indent=1,
                default=str,
            )
        )
        for t, v in c["table"].items():
            print(
                t,
                v["gate_label"],
                v["kick_label"],
                "E",
                v["E"],
                "+-",
                v["E_err"],
                "interior",
                v["interior_partition"],
                "r9",
                v["r9_partition"],
            )
    elif mode == "reread":
        J = load_json()
        for tag, r in J["rows"].items():
            f = os.path.join(OUT_NPZ, tag + ".npz")
            if r.get("status") != "OK" or not os.path.exists(f):
                continue
            cfg, p, pot = cfg_pot(r)
            M = np.load(f)["M"]
            r["own_reads"] = own_reads(M, cfg, p, pot, r["delta"], r["part"])
            log(f"reread {tag}: interior {r['own_reads']['interior_partition']}")
        save_json(J)
    elif mode == "jobs":
        for j in jobs_main():
            print(job_tag(j), j)
    else:
        raise SystemExit(f"unknown mode {mode}")


if __name__ == "__main__":
    main()
