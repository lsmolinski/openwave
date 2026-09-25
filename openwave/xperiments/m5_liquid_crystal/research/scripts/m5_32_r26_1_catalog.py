"""M5.32 R26-1: the generator-catalog audit. The M5.21.3 catalog (m5_21_3_a_4d.gen_catalog)
builds every velocity field as w (G M - M G^T), which is the ANTICOMMUTATOR for a rotation
and the COMMUTATOR for a boost, an antisymmetric matrix in both cases and not a tangent of the
symmetric field (R26-0 check a, in sympy; R25 found it on rot_z). This rung (1) states the
corrected catalog, a0 = w (G M + M G^T), normalized as the original (unit Frobenius over the
lattice, the null test kept), in THIS script, the platform file untouched; (2) shows the
(anti)symmetry of both on every stored field; (3) re-reads the six inertias kin(M; a0) old
against new on every locally stored 4x4 field (the R19-1 to R25-2 arrays and the M5.21.x end
fields), so every published inertia has its corrected value beside it; (4) censuses the
readers: every script that calls gen_catalog or builds a rotation tangent of its own, with
the construction line quoted and a verdict (AFFECTED: feeds the catalog's field to a read;
OWN: builds its own commutator; FD: differentiates a dressed family; NONE: no kinetic read);
the reads that were not stored are priced, not run.

EQUATIONS. Under M -> Lambda M Lambda^T, Lambda = exp(t G), G eta + eta G^T = 0:
    dM/dt = G M + M G^T   (rotation: [J, M]; boost: {K, M}),  kin(M; a0) = 4 h^3 sum_i <[a0, A_i]_eta, .>_eta.

PRE-REGISTERED LABELS
    CATALOG_FAULT_CONFIRMED   every catalog field antisymmetric to round-off on every stored
                              field, every corrected field symmetric, and at least one
                              stored inertia moves by over 1 percent under the correction
    CATALOG_FAULT_REFUTED     a catalog field is a tangent (symmetric) on some stored field

Modes: run [workers] | census | fields | relabel. Output: data/m5_32_r26_1_catalog.json.
Regenerate: run 4 about 15 min.
"""

import glob
import importlib.util
import json
import multiprocessing as mp
import os
import re
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
OUT_JSON = os.path.join(DATA, "m5_32_r26_1_catalog.json")
T0 = time.time()


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


B3 = _load("m5_21_3_a_4d", "m5_21_3_a_4d.py")


def log(msg):
    print(f"[{time.time() - T0:8.1f}s] {msg}", flush=True)


# ================= the corrected catalog =================
def gen_catalog_tangent(cfg, M):
    """the M5.21.3 catalog with the tangent formula: a0 = w (G M + M G^T) per generator, the
    same envelope, null test and unit-Frobenius normalization as gen_catalog."""
    w = B3.envelope(cfg)[..., None, None]
    lam, V = np.linalg.eigh(M[..., 1:4, 1:4])
    out = {}

    def local_rot(vhat):
        W = np.zeros(vhat.shape[:-1] + (4, 4))
        n1, n2, n3 = vhat[..., 0], vhat[..., 1], vhat[..., 2]
        W[..., 1, 2], W[..., 1, 3] = -n3, n2
        W[..., 2, 1], W[..., 2, 3] = n3, -n1
        W[..., 3, 1], W[..., 3, 2] = -n2, n1
        return W

    out["clock_local"] = local_rot(V[..., :, 2])
    out["plane_1d"] = local_rot(V[..., :, 0])
    Jz = np.zeros((4, 4))
    Jz[1, 2], Jz[2, 1] = -1.0, 1.0
    Jx = np.zeros((4, 4))
    Jx[2, 3], Jx[3, 2] = -1.0, 1.0
    Kz = np.zeros((4, 4))
    Kz[0, 3] = Kz[3, 0] = 1.0
    Kx = np.zeros((4, 4))
    Kx[0, 1] = Kx[1, 0] = 1.0
    for nm, Gm in (("rot_z", Jz), ("rot_x", Jx), ("boost_z", Kz), ("boost_x", Kx)):
        out[nm] = np.broadcast_to(Gm, M.shape)
    a0s = {}
    ref = np.sqrt(np.sum((w * M) ** 2))
    for nm, Gm in out.items():
        a0 = w * (Gm @ M + M @ Gm.swapaxes(-1, -2))
        nrm = np.sqrt(np.sum(a0 * a0))
        if nrm <= B3.GEN_NULL_TOL * ref:
            a0s[nm] = np.zeros_like(a0)
        else:
            a0s[nm] = a0 / nrm
    return a0s


# ================= the stored fields =================
def stored_fields():
    """every local npz holding a 4x4 field under key M (the M5.32 rungs and the M5.21.x ends)."""
    files = sorted(glob.glob(os.path.join(DATA, "**", "*.npz"), recursive=True))
    out = []
    for f in files:
        b = os.path.basename(f)
        if "r26_" in f or "_stage" in b or "_ckpt" in b or "smoke" in f or "fullkick" in f:
            continue
        try:
            with np.load(f, allow_pickle=True) as Z:
                if "M" not in Z.files:
                    continue
                shp = Z["M"].shape
        except Exception:  # noqa: BLE001
            continue
        if len(shp) == 5 and shp[3] == 4 and shp[4] == 4 and shp[0] == shp[1] == shp[2]:
            out.append(f)
    return out


def cfg_of_field(f, n):
    """the box of a stored field from its name (L from the tag, else the R21 default L 48)."""
    m = re.search(r"_L(\d+(?:\.\d+)?)", os.path.basename(f))
    L = float(m.group(1)) if m else 48.0
    m = re.search(r"_d(\d+(?:\.\d+)?)", os.path.basename(f))
    delta = float(m.group(1)) if m else 0.3
    return B3.base_cfg(s=-1.0, g=8.0, n=n, L=L, delta=delta)


def read_one(f):
    t0 = time.time()
    rel = os.path.relpath(f, DATA)
    try:
        M = np.load(f)["M"].astype(np.float64)
        n = M.shape[0]
        cfg = cfg_of_field(f, n)
        old = B3.gen_catalog(cfg, M)
        new = gen_catalog_tangent(cfg, M)
        rec = {"n": n, "L": cfg["L"], "delta": cfg["delta"], "gens": {}}
        for nm in old:
            ao, an = old[nm], new[nm]
            mo, mn = float(np.abs(ao).max()), float(np.abs(an).max())
            rec["gens"][nm] = {
                "old_symmetric_part_over_max": (
                    float(np.abs(ao + ao.swapaxes(-1, -2)).max() / mo) if mo > 0 else None
                ),
                "new_antisymmetric_part_over_max": (
                    float(np.abs(an - an.swapaxes(-1, -2)).max() / mn) if mn > 0 else None
                ),
                "kin_old": float(B3.kin_of(M, ao, cfg)),
                "kin_new": float(B3.kin_of(M, an, cfg)),
                "old_null": bool(mo == 0.0),
                "new_null": bool(mn == 0.0),
            }
            ko, kn = rec["gens"][nm]["kin_old"], rec["gens"][nm]["kin_new"]
            rec["gens"][nm]["ratio_new_over_old"] = kn / ko if ko > 0 else None
        Jz = np.zeros((4, 4))
        Jz[1, 2], Jz[2, 1] = -1.0, 1.0
        rec["kin_raw_commutator_Jz"] = float(B3.kin_of(M, Jz @ M - M @ Jz, cfg))
        rec["kin_raw_anticommutator_Jz"] = float(B3.kin_of(M, Jz @ M + M @ Jz, cfg))
        rec["status"] = "OK"
    except Exception as e:  # noqa: BLE001
        rec = {"status": "FAILED", "stop": repr(e)}
    rec["wall_s"] = round(time.time() - t0, 1)
    rec["file"] = rel
    log(f"{rel}: {rec['status']} {rec['wall_s']} s")
    return rec


# ================= the reader census =================
PATTERNS = {
    "calls_gen_catalog": r"gen_catalog\(|\[\"gen_catalog\"\]\(|\['gen_catalog'\]\(",
    "catalog_seeds_dynamics": r"gen_catalog\(.*\)\[|a0\s*=\s*\w+\.gen_catalog|Mt\s*=\s*om\s*\*\s*a0",
    "own_commutator": r"(J\w*|G\w*|W\w*|Jz|Jx|Jy)\s*@\s*M\w*\s*-\s*M\w*\s*@\s*(J\w*|G\w*|W\w*|Jz|Jx|Jy)",
    "own_anticommutator": r"(J\w*|G\w*|W\w*|Jz|Jx|Jy)\s*@\s*M\w*\s*\+\s*M\w*\s*@\s*(J\w*|G\w*|W\w*|Jz|Jx|Jy)",
    "catalog_form": r"@\s*M\w*\s*-\s*M\w*\s*@\s*\w+\.swapaxes",
    "finite_difference_family": r"a0_unit|dressed|\(M_?p\w*\s*-\s*M_?m\w*\)\s*/",
    "kin_read": r"\bkin_of\(|\bkin_grad\(|\btwist_read\(|\bkin\(|\ba0=[A-Za-z_]",
}


def census():
    rows = {}
    for f in sorted(glob.glob(os.path.join(HERE, "*.py"))):
        b = os.path.basename(f)
        if b.startswith("m5_32_r26_") or b == "m5_21_3_a_4d.py":
            continue
        try:
            src = open(f).read()
        except Exception:  # noqa: BLE001
            continue
        hits = {}
        for k, pat in PATTERNS.items():
            found = []
            for i, line in enumerate(src.split("\n"), 1):
                if re.search(pat, line) and not line.strip().startswith("#"):
                    found.append([i, line.strip()[:140]])
            if found:
                hits[k] = found
        if not hits or ("kin_read" not in hits and "calls_gen_catalog" not in hits):
            continue
        uses_catalog = "calls_gen_catalog" in hits
        own = ("own_commutator" in hits) or ("finite_difference_family" in hits)
        # AUDIT CORRECTION (R26-1 audit, 2026-09-25): a catalog field that seeds a dynamics run is a
        # kinetic read too (m5_21_9_e_larmor.py); the dict-call readers m5_21_3_f_confirm.py and
        # m5_21_3_c_films.py were missed by the first regex
        if uses_catalog and "catalog_seeds_dynamics" in hits and "kin_read" not in hits:
            hits["kin_read"] = hits["catalog_seeds_dynamics"]
        if uses_catalog and "kin_read" in hits:
            verdict = (
                "AFFECTED (the catalog's field reaches a kinetic read)"
                if not own
                else "MIXED (catalog and own tangents both read)"
            )
        elif own:
            verdict = "OWN (builds its own tangent)"
        elif uses_catalog:
            verdict = "CATALOG CALLED, no kinetic read found by pattern"
        else:
            verdict = "kinetic read with a tangent from elsewhere: inspect"
        rows[b] = {"verdict": verdict, "hits": hits}
    return rows


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "run"
    J = {"task": "M5.32 R26-1"}
    if mode in ("census", "run"):
        J["reader_census"] = census()
        log(f"census: {len(J['reader_census'])} scripts")
    if mode == "fields":
        for f in stored_fields():
            print(os.path.relpath(f, DATA))
        return
    if mode == "run":
        files = stored_fields()
        workers = min(int(sys.argv[2]) if len(sys.argv) > 2 else 4, 12)
        log(f"{len(files)} stored fields, {workers} workers")
        recs = []
        with ProcessPoolExecutor(max_workers=workers, mp_context=mp.get_context("spawn")) as ex:
            futs = [ex.submit(read_one, f) for f in files]
            for fut in as_completed(futs):
                recs.append(fut.result())
        recs.sort(key=lambda r: r["file"])
        J["fields"] = recs
        J.update(summarize(recs))
    if mode == "relabel":
        with open(OUT_JSON) as f:
            J = json.load(f)
        J.update(summarize(J["fields"]))
    J["wall_s"] = round(time.time() - T0, 1)
    with open(OUT_JSON, "w") as f:
        json.dump(J, f, indent=1, default=str)
    if "summary" in J:
        print(json.dumps(J["summary"], indent=1), J["label"])
    for b, r in J.get("reader_census", {}).items():
        print(b, "->", r["verdict"])


SYM_TOL = 1e-6  # relative (anti)symmetric part: the stored fields are float64 relaxations with
# round-off symmetric parts up to 1e-9 of the field's max entry (the first run used 1e-12 and
# labeled REFUTED on round-off; RUN-TIME CORRECTION 2026-09-25)


def summarize(recs):
    J = {}
    ok = [r for r in recs if r["status"] == "OK"]
    sym_old = max(
        (v["old_symmetric_part_over_max"] or 0.0) for r in ok for v in r["gens"].values()
    )
    asym_new = max(
        (v["new_antisymmetric_part_over_max"] or 0.0) for r in ok for v in r["gens"].values()
    )
    moved = [
        (r["file"], nm, v["ratio_new_over_old"])
        for r in ok
        for nm, v in r["gens"].items()
        if v["ratio_new_over_old"] is not None and abs(v["ratio_new_over_old"] - 1.0) > 0.01
    ]
    per_gen = {}
    for nm in ("rot_z", "rot_x", "boost_z", "boost_x", "clock_local", "plane_1d"):
        rat = [
            r["gens"][nm]["ratio_new_over_old"]
            for r in ok
            if r["gens"][nm]["ratio_new_over_old"] is not None
        ]
        per_gen[nm] = {
            "n": len(rat),
            "ratio_min": float(min(rat)) if rat else None,
            "ratio_max": float(max(rat)) if rat else None,
            "ratio_median": float(np.median(rat)) if rat else None,
        }
    J["summary"] = {
        "n_fields": len(recs),
        "n_ok": len(ok),
        "n_failed": len(recs) - len(ok),
        "failed": [r["file"] + ": " + str(r.get("stop")) for r in recs if r["status"] != "OK"],
        "sym_tol": SYM_TOL,
        "max_old_symmetric_part": sym_old,
        "max_new_antisymmetric_part": asym_new,
        "n_inertias_moved_over_1pct": len(moved),
        "n_channels_with_positive_kin_old": sum(
            1 for r in ok for v in r["gens"].values() if v["ratio_new_over_old"] is not None
        ),
        "n_channels_negative_kin_old_ratio_null": sum(
            1
            for r in ok
            for v in r["gens"].values()
            if v["ratio_new_over_old"] is None and not v["old_null"]
        ),
        "n_sign_flips": sum(
            1 for r in ok for v in r["gens"].values() if v["kin_old"] * v["kin_new"] < 0
        ),
        "rule_note": "the symmetry half of the rule is an algebraic identity for a symmetric M (it measures storage symmetry); the load-bearing content is the algebra (R26-0 check a) plus the inertia shift (audit)",
        "per_generator": per_gen,
    }
    J["label"] = (
        "CATALOG_FAULT_CONFIRMED"
        if (sym_old < SYM_TOL and asym_new < SYM_TOL and moved)
        else "CATALOG_FAULT_REFUTED"
    )
    J["provenance"] = (
        "AUDIT CORRECTION (R26-1 audit, 2026-09-25): the antisymmetric form was documented at the"
        " M5.21.3 audit itself (m5_21_3_f_confirm.py: `a0_conj = w (G M + M G^T)` the symmetric"
        " physical velocity, the instrument's `w (G M - M G^T)` the antisymmetric variant, audit"
        " section C7.4, the corrected reads stored as `kinconj_*`), restated at M5.21.9 (`a0_conj`)"
        " and flagged at M5.32 R1 and R14-0; R25 re-found a known fault, this rung re-confirms it"
        " on every stored field"
    )
    J["not_run_priced"] = {
        "M5.21.5 / 21.9 / 21.14 / 21.15 / 21.16 reads": "the end fields that exist locally are in `fields`; the reads whose arrays were not kept "
        "would need their scripts re-run (hours each, the M5.21.x pools); priced, not run: the corrected inertia on every KEPT end field is the substitute"
    }
    return J


if __name__ == "__main__":
    main()
