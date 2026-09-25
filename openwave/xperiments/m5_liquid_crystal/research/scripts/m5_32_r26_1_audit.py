"""M5.32 R26-1 adversarial audit: the generator-catalog audit refuted with own code.

The audited claims live in ../data/m5_32_r26_1_catalog.json (the audited script was NOT read):
  A. every catalog field (gen_catalog: a0 = w (G M - M G^T), unit Frobenius) is antisymmetric on
     every stored field; the corrected field (w (G M + M G^T)) is symmetric.
  B. the inertias kin_of(M, a0, cfg) move under the correction (per field, per generator), and
     the raw rigid reads kin(M; [J_z, M]) = 14071 against kin(M; J_z M + M J_z) = 37180 on the
     S1 n 32 L 48 field.
  C. the reader census (a regex pass over scripts/) with a verdict per script.
  D. the label CATALOG_FAULT_CONFIRMED and its rule.

Own instruments here: own eta contraction (the sign-matrix form of tr(eta F eta G^T)), own
stencil, own kinetic sum, own catalog construction (old and corrected), own envelope, own grep
census. The platform (m5_21_3_a_4d.py) is imported only for base_cfg / gen_catalog / kin_of /
envelope, the objects the claims are about, so that the platform can be checked against the own
implementation on one field.

Run: /opt/anaconda3/envs/master312/bin/python3 m5_32_r26_1_audit.py
Out: ../data/m5_32_r26_1_audit.json
"""

from __future__ import annotations

import importlib.util
import json
import os
import re
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
AUDITED = os.path.join(DATA, "m5_32_r26_1_catalog.json")
OUT = os.path.join(DATA, "m5_32_r26_1_audit.json")

T0 = time.time()


def _load(name, fn):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fn))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


B3 = _load("m5_21_3_a_4d", "m5_21_3_a_4d.py")

# ----------------------------------------------------------------------------- own algebra
ETA_D = np.array([-1.0, 1.0, 1.0, 1.0])
ETA = np.diag(ETA_D)
SIGN = np.outer(ETA_D, ETA_D)  # SIGN[b, c] = eta_b eta_c: tr(eta F eta G^T) = sum SIGN F G


def own_inner(F, G):
    """<F, G>_eta per cell = tr(eta F eta G^T) = sum_{bc} eta_b eta_c F_bc G_bc (eta diagonal)."""
    return np.sum(SIGN * F * G, axis=(-2, -1))


def own_comm(a, b):
    return a @ ETA @ b - b @ ETA @ a


def own_d1(f, ax, h, side):
    """one-sided difference with a zero row at the unreachable edge (the platform's fwd/bwd)."""
    out = np.zeros_like(f)
    df = np.diff(f, axis=ax) / h
    sl = [slice(None)] * f.ndim
    sl[ax] = slice(0, -1) if side == "fwd" else slice(1, None)
    out[tuple(sl)] = df
    return out


def own_kin(M, a0, h):
    """4 h^3 sum_i <[a0, A_i]_eta, [a0, A_i]_eta>_eta, symmetric stencil = (fwd + bwd) / 2."""
    k = 0.0
    for side in ("fwd", "bwd"):
        for ax in range(3):
            F = own_comm(a0, own_d1(M, ax, h, side))
            k += 0.5 * 4.0 * float(np.sum(own_inner(F, F)))
    return h**3 * k


def own_envelope(n, h, renv=10.0):
    x = (np.arange(n) - (n - 1) / 2.0) * h
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    r = np.sqrt(X * X + Y * Y + Z * Z)
    return np.exp(-((r / renv) ** 4))


def own_generators(M):
    """the six generator matrices G(x) of the catalog, all antisymmetric (rotations, local
    rotations) or symmetric (boosts); the Lorentz condition G eta + eta G^T = 0 is checked."""
    lam, V = np.linalg.eigh(M[..., 1:4, 1:4])

    def local_rot(v):
        W = np.zeros(v.shape[:-1] + (4, 4))
        n1, n2, n3 = v[..., 0], v[..., 1], v[..., 2]
        W[..., 1, 2], W[..., 2, 1] = -n3, n3
        W[..., 1, 3], W[..., 3, 1] = n2, -n2
        W[..., 2, 3], W[..., 3, 2] = -n1, n1
        return W

    Jz = np.zeros((4, 4))
    Jz[1, 2], Jz[2, 1] = -1.0, 1.0
    Jx = np.zeros((4, 4))
    Jx[2, 3], Jx[3, 2] = -1.0, 1.0
    Kz = np.zeros((4, 4))
    Kz[0, 3] = Kz[3, 0] = 1.0
    Kx = np.zeros((4, 4))
    Kx[0, 1] = Kx[1, 0] = 1.0
    gens = {
        "clock_local": local_rot(V[..., :, 2]),
        "plane_1d": local_rot(V[..., :, 0]),
        "rot_z": np.broadcast_to(Jz, M.shape),
        "rot_x": np.broadcast_to(Jx, M.shape),
        "boost_z": np.broadcast_to(Kz, M.shape),
        "boost_x": np.broadcast_to(Kx, M.shape),
    }
    lor = {}
    for nm, G in gens.items():
        lor[nm] = float(np.max(np.abs(G @ ETA + ETA @ G.swapaxes(-1, -2))))
    return gens, lor


def own_catalog(M, w, sign):
    """sign -1: the catalog's a0 = w (G M - M G^T); sign +1: the corrected w (G M + M G^T).
    Unit Frobenius norm over the lattice, null channel below 1e-10 of |w M|."""
    gens, lor = own_generators(M)
    ref = np.sqrt(np.sum((w * M) ** 2))
    out = {}
    for nm, G in gens.items():
        a = w * (G @ M + sign * M @ G.swapaxes(-1, -2))
        nrm = np.sqrt(np.sum(a * a))
        out[nm] = np.zeros_like(a) if nrm <= 1e-10 * ref else a / nrm
    return out, lor


def sym_part_over_max(a):
    """max |(a + a^T) / 2| over max |a| (0 for a null channel)."""
    mx = float(np.max(np.abs(a)))
    if mx == 0.0:
        return 0.0
    return float(np.max(np.abs(0.5 * (a + a.swapaxes(-1, -2)))) / mx)


def antisym_part_over_max(a):
    mx = float(np.max(np.abs(a)))
    if mx == 0.0:
        return 0.0
    return float(np.max(np.abs(0.5 * (a - a.swapaxes(-1, -2)))) / mx)


def rel(a, b):
    return abs(a - b) / max(abs(a), abs(b), 1e-300)


GENS = ("clock_local", "plane_1d", "rot_z", "rot_x", "boost_z", "boost_x")

# ----------------------------------------------------------------------------- the fields
FIELDS = [
    # (subdir, file, n, L, delta, note)
    ("m5_32_r25_2", "S1_d0.3_w25_n32_L48.npz", 32, 48.0, 0.3, "named in the task"),
    ("m5_32_r25_2", "S1_d0.3_w25_n64_L96.npz", 64, 96.0, 0.3, "named in the task"),
    ("m5_32_r23_1", "rad_pin_d0.3_w25_c0.0001_n32_L48.npz", 32, 48.0, 0.3, "named in the task"),
    ("m5_32_r21_1", "S1_Bseed_g8_d0.3_w25_n32.npz", 32, 48.0, 0.3, "auditor's choice (R21-1)"),
    (
        "m5_32_r19_1",
        "Gam_dr_same_d24_n48_g32.npz",
        48,
        48.0,
        0.3,
        "auditor's choice (R19-1); the audited JSON parsed its `_d24` tag as delta 24.0",
    ),
    (
        "",
        "m5_22_4_p1_deut.npz",
        32,
        48.0,
        0.3,
        "the field carrying the audited max symmetric part 1.2e-9 (stored float32)",
    ),
]

audited = json.load(open(AUDITED))
aud_by_file = {x["file"]: x for x in audited["fields"]}

field_rows = []
for sub, fn, n, L, delta, note in FIELDS:
    t1 = time.time()
    path = os.path.join(DATA, sub, fn) if sub else os.path.join(DATA, fn)
    Z = np.load(path)
    M_raw = Z["M"]
    M = M_raw.astype(np.float64)
    cfg = B3.base_cfg(s=-1.0, g=8.0, n=n, L=L, delta=delta)
    h = cfg["h"]
    w = own_envelope(n, h)[..., None, None]
    env_diff = float(np.max(np.abs(w[..., 0, 0] - B3.envelope(cfg))))
    M_asym_abs = float(np.max(np.abs(M - M.swapaxes(-1, -2))))
    M_max = float(np.max(np.abs(M)))
    old_own, lor = own_catalog(M, w, -1.0)
    new_own, _ = own_catalog(M, w, +1.0)
    old_plat = B3.gen_catalog(cfg, M)
    key = f"{sub}/{fn}" if sub else fn
    key = key if key in aud_by_file else None
    aud = aud_by_file.get(key, {})
    gens_out = {}
    worst_kin_rel = 0.0
    for nm in GENS:
        own_vs_plat = float(np.max(np.abs(old_own[nm] - old_plat[nm])))
        s_old = sym_part_over_max(old_plat[nm])
        a_new = antisym_part_over_max(new_own[nm])
        k_old = float(B3.kin_of(M, old_plat[nm], cfg))
        k_new = float(B3.kin_of(M, new_own[nm], cfg))
        row = {
            "own_catalog_minus_platform_max_abs": own_vs_plat,
            "old_symmetric_part_over_max": s_old,
            "new_antisymmetric_part_over_max": a_new,
            "kin_old": k_old,
            "kin_new": k_new,
            "ratio_new_over_old": k_new / k_old if k_old != 0.0 else None,
            "lorentz_condition_max": lor[nm],
        }
        ag = aud.get("gens", {}).get(nm)
        if ag:
            row["audited_kin_old"] = ag["kin_old"]
            row["audited_kin_new"] = ag["kin_new"]
            row["rel_kin_old"] = rel(k_old, ag["kin_old"])
            row["rel_kin_new"] = rel(k_new, ag["kin_new"])
            row["audited_old_symmetric_part_over_max"] = ag["old_symmetric_part_over_max"]
            row["audited_new_antisymmetric_part_over_max"] = ag["new_antisymmetric_part_over_max"]
            worst_kin_rel = max(worst_kin_rel, row["rel_kin_old"], row["rel_kin_new"])
        gens_out[nm] = row
    Jz = np.zeros((4, 4))
    Jz[1, 2], Jz[2, 1] = -1.0, 1.0
    raw_comm = float(B3.kin_of(M, Jz @ M - M @ Jz, cfg))
    raw_anti = float(B3.kin_of(M, Jz @ M + M @ Jz, cfg))
    fr = {
        "file": fn,
        "subdir": sub,
        "note": note,
        "n": n,
        "L": L,
        "h": h,
        "stored_dtype": str(M_raw.dtype),
        "M_max_abs": M_max,
        "M_asymmetry_max_abs": M_asym_abs,
        "M_asymmetry_over_max": M_asym_abs / M_max,
        "own_envelope_minus_platform": env_diff,
        "in_audited_json": key is not None,
        "audited_n_L_delta": [aud.get("n"), aud.get("L"), aud.get("delta")] if aud else None,
        "gens": gens_out,
        "max_old_symmetric_part": max(v["old_symmetric_part_over_max"] for v in gens_out.values()),
        "max_new_antisymmetric_part": max(
            v["new_antisymmetric_part_over_max"] for v in gens_out.values()
        ),
        "kin_raw_commutator_Jz": raw_comm,
        "kin_raw_anticommutator_Jz": raw_anti,
        "worst_rel_vs_audited_kin": worst_kin_rel if aud else None,
        "wall_s": round(time.time() - t1, 2),
    }
    if aud:
        fr["audited_kin_raw_commutator_Jz"] = aud["kin_raw_commutator_Jz"]
        fr["audited_kin_raw_anticommutator_Jz"] = aud["kin_raw_anticommutator_Jz"]
        fr["rel_raw_commutator"] = rel(raw_comm, aud["kin_raw_commutator_Jz"])
        fr["rel_raw_anticommutator"] = rel(raw_anti, aud["kin_raw_anticommutator_Jz"])
    # the own kinetic contraction against the platform, on the S1 n 32 field (claim B)
    if fn == "S1_d0.3_w25_n32_L48.npz":
        own = {}
        worst = 0.0
        for nm in GENS:
            ko = own_kin(M, old_plat[nm], h)
            kn = own_kin(M, new_own[nm], h)
            own[nm] = {
                "own_kin_old": ko,
                "own_kin_new": kn,
                "rel_old_vs_kin_of": rel(ko, gens_out[nm]["kin_old"]),
                "rel_new_vs_kin_of": rel(kn, gens_out[nm]["kin_new"]),
            }
            worst = max(worst, own[nm]["rel_old_vs_kin_of"], own[nm]["rel_new_vs_kin_of"])
        rc = own_kin(M, Jz @ M - M @ Jz, h)
        ra = own_kin(M, Jz @ M + M @ Jz, h)
        own["raw"] = {
            "own_kin_raw_commutator_Jz": rc,
            "own_kin_raw_anticommutator_Jz": ra,
            "rel_commutator_vs_kin_of": rel(rc, raw_comm),
            "rel_anticommutator_vs_kin_of": rel(ra, raw_anti),
        }
        worst = max(worst, own["raw"]["rel_commutator_vs_kin_of"])
        worst = max(worst, own["raw"]["rel_anticommutator_vs_kin_of"])
        own["worst_rel_own_vs_platform"] = worst
        # delta is not consumed by the catalog or by kin_of: same numbers at delta 24
        cfg24 = B3.base_cfg(s=-1.0, g=8.0, n=n, L=L, delta=24.0)
        k24 = float(B3.kin_of(M, B3.gen_catalog(cfg24, M)["rot_z"], cfg24))
        own["delta_independence_rot_z_kin_old_at_delta_24_rel"] = rel(
            k24, gens_out["rot_z"]["kin_old"]
        )
        fr["own_contraction"] = own
    field_rows.append(fr)
    print(
        f"[field] {fn} n{n} L{L:g} h{h:g} dtype {M_raw.dtype} M_asym/max {fr['M_asymmetry_over_max']:.2e} "
        f"max_sym_old {fr['max_old_symmetric_part']:.2e} max_asym_new {fr['max_new_antisymmetric_part']:.2e} "
        f"raw [J,M] {raw_comm:.1f} raw {{J,M}} {raw_anti:.1f} worst_rel_vs_audited "
        f"{fr['worst_rel_vs_audited_kin'] if aud else 'n/a'} ({fr['wall_s']} s)",
        flush=True,
    )

# ------------------------------------------------ the identity: what the symmetry test can fail on
rng = np.random.default_rng(2611)
Ms = B3.sym4(rng.normal(size=(6, 6, 6, 4, 4)))
Mn = rng.normal(size=(6, 6, 6, 4, 4))  # NOT symmetric
w6 = own_envelope(6, 1.0)[..., None, None]
ident = {}
for tag, Mt in (("symmetric_M", Ms), ("nonsymmetric_M", Mn)):
    o, _ = own_catalog(Mt, w6, -1.0)
    c, _ = own_catalog(Mt, w6, +1.0)
    ident[tag] = {
        "max_old_symmetric_part": max(sym_part_over_max(o[g]) for g in GENS),
        "max_new_antisymmetric_part": max(antisym_part_over_max(c[g]) for g in GENS),
    }

# ---------------------------------------------------------------- re-derivations from the JSON
summ = audited["summary"]
ok_fields = [x for x in audited["fields"] if x["status"] == "OK"]
n_chan = sum(len(x["gens"]) for x in ok_fields)
n_null_ratio = sum(
    1 for x in ok_fields for v in x["gens"].values() if v["ratio_new_over_old"] is None
)
n_null_ratio_neg_old = sum(
    1
    for x in ok_fields
    for v in x["gens"].values()
    if v["ratio_new_over_old"] is None and v["kin_old"] < 0
)
n_null_boost = sum(
    1
    for x in ok_fields
    for g, v in x["gens"].items()
    if v["ratio_new_over_old"] is None and g.startswith("boost")
)
n_boost = sum(1 for x in ok_fields for g in x["gens"] if g.startswith("boost"))
moved_nonnull = sum(
    1
    for x in ok_fields
    for v in x["gens"].values()
    if v["ratio_new_over_old"] is not None and abs(v["ratio_new_over_old"] - 1.0) > 0.01
)
moved_all = sum(
    1
    for x in ok_fields
    for v in x["gens"].values()
    if v["kin_old"] != 0.0 and abs(v["kin_new"] / v["kin_old"] - 1.0) > 0.01
)
sign_flips = sum(
    1 for x in ok_fields for v in x["gens"].values() if v["kin_old"] * v["kin_new"] < 0.0
)
max_sym_json = max(v["old_symmetric_part_over_max"] for x in ok_fields for v in x["gens"].values())
max_asym_json = max(
    v["new_antisymmetric_part_over_max"] for x in ok_fields for v in x["gens"].values()
)
delta_tags = sorted({x["delta"] for x in ok_fields})
n_delta_over_2 = sum(1 for x in ok_fields if x["delta"] > 2.0)
per_gen_n = {
    g: sum(1 for x in ok_fields if x["gens"][g]["ratio_new_over_old"] is not None) for g in GENS
}
json_checks = {
    "n_ok_fields": len(ok_fields),
    "n_channels": n_chan,
    "n_ratio_null": n_null_ratio,
    "n_ratio_null_with_negative_kin_old": n_null_ratio_neg_old,
    "n_boost_channels": n_boost,
    "n_boost_channels_with_null_ratio": n_null_boost,
    "moved_over_1pct_among_nonnull_ratio": moved_nonnull,
    "moved_over_1pct_among_all_nonzero_kin_old": moved_all,
    "sign_flips_old_vs_new": sign_flips,
    "audited_n_inertias_moved_over_1pct": summ["n_inertias_moved_over_1pct"],
    "per_generator_n_recount": per_gen_n,
    "audited_per_generator_n": {g: summ["per_generator"][g]["n"] for g in GENS},
    "max_old_symmetric_part_recount": max_sym_json,
    "max_new_antisymmetric_part_recount": max_asym_json,
    "audited_max_old_symmetric_part": summ["max_old_symmetric_part"],
    "audited_max_new_antisymmetric_part": summ["max_new_antisymmetric_part"],
    "distinct_delta_values_parsed": delta_tags,
    "n_fields_with_delta_over_2": n_delta_over_2,
}

# ------------------------------------------------------------------------ the census, own grep
census = audited["reader_census"]
call_re = re.compile(r"""(?:\.|\[["']|^|\s)gen_catalog(?:["']\])?\s*\(""")
own_callers = {}
for fn in sorted(os.listdir(HERE)):
    if not fn.endswith(".py") or fn.startswith("m5_32_r26_") or fn == "m5_21_3_a_4d.py":
        continue
    with open(os.path.join(HERE, fn), errors="replace") as f:
        lines = f.read().splitlines()
    hits = [(i + 1, ln.strip()[:100]) for i, ln in enumerate(lines) if call_re.search(ln)]
    if hits:
        own_callers[fn] = hits
census_callers = {k for k, v in census.items() if v["hits"].get("calls_gen_catalog")}
missing_from_census = sorted(set(own_callers) - set(census))
call_not_flagged = sorted(set(own_callers) & set(census) - census_callers)
kin_hits = [(k, ln, txt) for k, v in census.items() for ln, txt in v["hits"].get("kin_read", [])]
kin_hits_no_kin = [(k, ln, txt) for k, ln, txt in kin_hits if "kin" not in txt.lower()]
inspect_on_noise = sorted(
    k
    for k, v in census.items()
    if "inspect" in v["verdict"]
    and v["hits"].get("kin_read")
    and all("kin" not in t.lower() for _, t in v["hits"]["kin_read"])
)

# the spot checks: read by the auditor (the sites are quoted by line for the record)
SPOT = {
    "m5_32_lagrangian.py": {
        "audited_verdict": census["m5_32_lagrangian.py"]["verdict"],
        "site": "lines 631-640: a0s = B3.gen_catalog(cfgc, Mc); k_cert = B3.kin_of(Mc, a0s[nm]); "
        "k_reg = term_energy(I1, ..., a0s[nm], 1.0) - e_static; compared to the stored "
        "M5.21.16 CHAN rows",
        "tangent": "the catalog's antisymmetric field (boost_z, boost_x, clock_local)",
        "reaches_kinetic_read": True,
        "auditor_verdict": "RIGHT as a label; the read is a registry-vs-certified-vs-stored "
        "consistency test (the same a0 on both sides), so the fault leaves its pass/fail "
        "untouched; the gradient check at 758 uses a random SYMMETRIC a0 (B3.sym4) and is clean",
    },
    "m5_22_4_b_omega.py": {
        "audited_verdict": census["m5_22_4_b_omega.py"]["verdict"],
        "site": "lines 138-142 (p2: kin_of + twist_read per catalog channel on the stored p1 "
        "fields) and 167-173 (p3: FIRE relaxation with a0 = catalog clock_local on an omega "
        "ladder, E = E_u + E_v + omega^2 kin_end)",
        "tangent": "the catalog's antisymmetric field, as the velocity of a relaxation",
        "reaches_kinetic_read": True,
        "auditor_verdict": "RIGHT; the p2 inertias, the argmin channel and the p3 ladder are "
        "reads on the anticommutator velocity",
    },
    "m5_21_16_b_field.py": {
        "audited_verdict": census["m5_21_16_b_field.py"]["verdict"],
        "site": "stage_chan lines 156-172: a0s = B3.gen_catalog(cfg, M) into the script's own "
        "kin_of (eta and flip contractions), the pass rule boost_kin_eta_all_negative; "
        "stage_dress lines 176-205 uses C14's ec.a0_base (a tangent from m5_21_14, not the "
        "catalog)",
        "tangent": "CHAN stage: the catalog's antisymmetric field; DRESS stage: C14 a0_base",
        "reaches_kinetic_read": True,
        "auditor_verdict": "RIGHT (MIXED); the CHAN pass rule and the stored CHAN rows that "
        "m5_32_lagrangian.py checks against are anticommutator reads",
    },
    "m5_21_9_d_fixedj.py": {
        "audited_verdict": census["m5_21_9_d_fixedj.py"]["verdict"],
        "site": "a0_conj lines 57-75: a0 = w [W, M] (the COMMUTATOR, unit norm), the default "
        "conv='conj'; the conv != 'conj' branch (lines 80-81, 101-102, 132-133) uses the "
        "catalog's clock_local; kin_of at 82, 103, 134",
        "tangent": "conj rows: the commutator (correct); catalog rows: the anticommutator",
        "reaches_kinetic_read": True,
        "auditor_verdict": "RIGHT (MIXED); note the docstring at 57-63 already calls the "
        "catalog's G M - M G^T an antisymmetric probe flow that exits the symmetric "
        "configuration space (M5.21.9)",
    },
    "m5_32_r12_a_ring.py": {
        "audited_verdict": census["m5_32_r12_a_ring.py"]["verdict"],
        "site": "a0_local lines 159-168: a0 = J M - M J, J the local leading-eigenvector "
        "rotation (antisymmetric), so the COMMUTATOR; kin_of at 180, own density at 185-193; "
        "selftest at 270-281 matches B8.a0_unit to 1e-8",
        "tangent": "the commutator, own construction; no gen_catalog call",
        "reaches_kinetic_read": False,
        "auditor_verdict": "RIGHT (OWN)",
    },
    "m5_32_r25_2_charge.py": {
        "audited_verdict": census["m5_32_r25_2_charge.py"]["verdict"],
        "site": "spin_gate_reads lines 466-490: a0 = Jz @ M - M @ Jz (the commutator, raw) into "
        "kin_of for C; C_cat = kin_of on the catalog's rot_z kept as a reference",
        "tangent": "the commutator for the verdict, the catalog's field as a labeled reference",
        "reaches_kinetic_read": True,
        "auditor_verdict": "RIGHT (MIXED)",
    },
    "m5_21_9_e_larmor.py": {
        "audited_verdict": census["m5_21_9_e_larmor.py"]["verdict"],
        "site": "setup lines 148-165: a0 = INS4.gen_catalog(cfg, M)['clock_local']; "
        "Mt = om * a0 * free is the INITIAL VELOCITY of the leapfrog; e_tot at 129-132 reads "
        "the kinetic energy 0.5 h^3 sum(Mt^2)",
        "tangent": "the catalog's antisymmetric field seeds the dynamics",
        "reaches_kinetic_read": True,
        "auditor_verdict": "WRONG: the census says OWN; the catalog field is the dynamics' "
        "initial velocity and its Frobenius kinetic energy is read (the regex looks for "
        "kin_of-shaped reads and the own_commutator hit overrode the catalog call)",
    },
    "m5_21_3_f_confirm.py": {
        "audited_verdict": "NOT LISTED (absent from reader_census)",
        "site": "line 70: a0s = g4['gen_catalog'](cfg, M) (a dict call through runpy, missed by "
        "the census regex); kin_of on the catalog's fields at 73 and 90 (AFFECTED); and "
        "conj_catalog at 32-60 CODES THE CORRECTED FIELD w (G M + M G^T) with kin_of on it at "
        "65-66 (the kinconj_* rows), labeled the audit-adopted variant, audit § C7.4",
        "tangent": "both: the catalog's anticommutator and the corrected conjugation tangent",
        "reaches_kinetic_read": True,
        "auditor_verdict": "MISSING from the census; the M5.21.3 record already held the "
        "corrected reads and named the catalog's form the antisymmetric variant",
    },
    "m5_21_3_c_films.py": {
        "audited_verdict": "NOT LISTED (absent from reader_census)",
        "site": "line 69: a0s = g4['gen_catalog'](cfg, M) (dict call); the orbit itself is built "
        "by exact exponentiation L M L^T at 76-90",
        "tangent": "the exact orbit; whether a0s reaches a read was not traced here",
        "reaches_kinetic_read": None,
        "auditor_verdict": "MISSING from the census (a caller the regex did not see)",
    },
}
# provenance sites found while reading (line-quoted; the auditor read these files)
PROVENANCE = [
    "m5_21_3_f_confirm.py docstring (3): 'THE AUDIT-ADOPTED VARIANT: the conjugation-orbit "
    "tangent a0_conj = w (G M + M G^T) (the symmetric, physical velocity of M -> L M L^T; the "
    "instrument's a0 = w (G M - M G^T) is the antisymmetric variant, audit § C7.4)'",
    "m5_21_9_d_fixedj.py a0_conj docstring (57-63): 'gen_catalog's GM - MG^T probe flow is "
    "antisymmetric and exits the symmetric configuration space'",
    "m5_32_r1_a_symbolic.py 362 and 410: an a0_symmetry column and 'J1 @ d - d @ J1.T  # "
    "antisymmetric (gen_catalog)'",
    "m5_32_r14_0_verify.py 461-464: 'gen_catalog's channels are the antisymmetric probes "
    "X M - M X^T flagged at R1; here every channel is a symmetric tangent'",
]

# ----------------------------------------------------------------------------- the verdicts
named = {fr["file"]: fr for fr in field_rows}
s1 = named["S1_d0.3_w25_n32_L48.npz"]
in_json = [fr for fr in field_rows if fr["in_audited_json"]]
max_sym_mine = max(fr["max_old_symmetric_part"] for fr in field_rows)
max_asym_mine = max(fr["max_new_antisymmetric_part"] for fr in field_rows)
worst_kin = max(fr["worst_rel_vs_audited_kin"] for fr in in_json)
worst_raw = max(max(fr["rel_raw_commutator"], fr["rel_raw_anticommutator"]) for fr in in_json)
worst_own = s1["own_contraction"]["worst_rel_own_vs_platform"]
worst_construct = max(
    v["own_catalog_minus_platform_max_abs"] for fr in field_rows for v in fr["gens"].values()
)
deut = named["m5_22_4_p1_deut.npz"]

B_expected = {
    "clock_local": (1.56, 1.24),
    "plane_1d": (1.55, 1.15),
    "rot_z": (1.72, 1.84),
    "rot_x": (1.35, 1.27),
    "boost_z": (-0.34, -0.37),
    "boost_x": (-0.57, -0.55),
}
B_ok = all(
    abs(s1["gens"][g]["kin_old"] - o) < 0.006 and abs(s1["gens"][g]["kin_new"] - n_) < 0.006
    for g, (o, n_) in B_expected.items()
)
raw_ok = (
    abs(s1["kin_raw_commutator_Jz"] - 14071) < 1
    and abs(s1["kin_raw_anticommutator_Jz"] - 37180) < 1
)

claims = {}
claims["A"] = {
    "method": "own catalog (own generators, own envelope, own normalization) for both signs on "
    "6 stored fields; symmetric/antisymmetric parts over the max entry; the platform's "
    "gen_catalog compared entrywise to the own old field; the identity probed on a random "
    "symmetric and a random nonsymmetric M",
    "mine": {
        "max_old_symmetric_part_6_fields": max_sym_mine,
        "max_new_antisymmetric_part_6_fields": max_asym_mine,
        "own_old_minus_platform_max_abs": worst_construct,
        "lorentz_condition_max": max(
            v["lorentz_condition_max"] for fr in field_rows for v in fr["gens"].values()
        ),
        "deut_field": {
            "stored_dtype": deut["stored_dtype"],
            "M_asymmetry_over_max": deut["M_asymmetry_over_max"],
            "max_old_symmetric_part": deut["max_old_symmetric_part"],
        },
        "identity_probe": ident,
    },
    "audited": {
        "max_old_symmetric_part": summ["max_old_symmetric_part"],
        "max_new_antisymmetric_part": summ["max_new_antisymmetric_part"],
        "n_fields": summ["n_ok"],
    },
    "verdict": (
        "CONFIRMED"
        if max_sym_mine < 1e-6 and max_asym_mine < 1e-6 and worst_construct < 1e-12
        else "REFUTED"
    ),
    "note": "the antisymmetry of w (G M - M G^T) for antisymmetric G (and symmetric G) on a "
    "symmetric M is an algebraic identity, so the numerical pass measures the stored M's own "
    "symmetry: the audited 1.2e-9 ceiling is the float32-stored m5_22_4_p1_deut field "
    f"(M asymmetry {deut['M_asymmetry_over_max']:.1e} of max); a nonsymmetric M fails the test "
    f"at {ident['nonsymmetric_M']['max_old_symmetric_part']:.2f}, so the test can fail only there",
}
claims["B"] = {
    "method": "own construction of the 12 fields, the platform's kin_of on each, on 5 fields the "
    "audited JSON holds; on S1 n 32 the own kinetic contraction (sign-matrix eta inner, own "
    "one-sided stencils) against kin_of; the raw [J_z, M] and {J_z, M} reads with no envelope",
    "mine": {
        "S1_n32_L48": {g: (s1["gens"][g]["kin_old"], s1["gens"][g]["kin_new"]) for g in GENS},
        "S1_n32_L48_raw_commutator_Jz": s1["kin_raw_commutator_Jz"],
        "S1_n32_L48_raw_anticommutator_Jz": s1["kin_raw_anticommutator_Jz"],
        "worst_rel_vs_audited_kin_5_fields": worst_kin,
        "worst_rel_vs_audited_raw_5_fields": worst_raw,
        "own_contraction_vs_kin_of_worst_rel": worst_own,
        "delta_independence_rel": s1["own_contraction"][
            "delta_independence_rot_z_kin_old_at_delta_24_rel"
        ],
    },
    "audited": {
        "S1_n32_L48_stated": B_expected,
        "raw_stated": (14071, 37180),
    },
    "verdict": (
        "CONFIRMED"
        if (B_ok and raw_ok and worst_kin < 1e-6 and worst_raw < 1e-6 and worst_own < 1e-9)
        else "REFUTED"
    ),
    "note": "the platform's kin_of agrees with the own contraction to round-off; the JSON's "
    "ratio_new_over_old is null whenever kin_old < 0 (623 of 2136 channels, 542 of them boosts), "
    "so the boost rows carry kin_old / kin_new only; delta does not enter the catalog or kin_of "
    "(same rot_z inertia at delta 24), so the R19-1 tags parsed as delta 10 to 30 are a "
    "metadata mislabel with no effect on the numbers",
}
spot_right = sum(1 for k in SPOT if SPOT[k]["auditor_verdict"].startswith("RIGHT"))
claims["C"] = {
    "method": "own regex over scripts/*.py for gen_catalog CALLS in both spellings "
    "(.gen_catalog( and ['gen_catalog']( ), compared to the census; the census's kin_read "
    "hit lines screened for the substring kin; the named scripts read at the quoted lines",
    "mine": {
        "own_gen_catalog_callers": sorted(own_callers),
        "n_own_callers": len(own_callers),
        "callers_missing_from_census": missing_from_census,
        "callers_in_census_without_calls_gen_catalog_hit": call_not_flagged,
        "census_size": len(census),
        "kin_read_hits": len(kin_hits),
        "kin_read_hits_without_kin_substring": len(kin_hits_no_kin),
        "inspect_verdicts_resting_only_on_such_hits": inspect_on_noise,
        "spot_checks": SPOT,
        "spot_checks_right_of_named": spot_right,
    },
    "audited": {
        "n_census_calls_gen_catalog": len(census_callers),
        "census_callers": sorted(census_callers),
    },
    "verdict": "QUALIFIED",
    "note": "the six named verdicts hold (lagrangian AFFECTED is a same-a0 consistency test, "
    "so its pass/fail is untouched); the census misses the two callers that reach the catalog "
    "through a runpy dict (m5_21_3_f_confirm.py: AFFECTED reads AND the corrected kinconj_* "
    "reads; m5_21_3_c_films.py), verdicts m5_21_9_e_larmor.py OWN although the catalog field "
    "seeds the dynamics' velocity, and 59 of 246 kin_read hits are regex noise (alpha0, omega0, "
    "a0= signatures); six inspect verdicts rest on that noise alone",
}
rule_holds = (
    max_sym_mine < 1e-6
    and max_asym_mine < 1e-6
    and any(
        v["ratio_new_over_old"] is not None and abs(v["ratio_new_over_old"] - 1.0) > 0.01
        for fr in field_rows
        for v in fr["gens"].values()
    )
)
claims["D"] = {
    "method": "the rule re-applied to the own numbers on 6 fields; the count rules recomputed "
    "from the JSON; the provenance of the antisymmetric form traced in the record's scripts",
    "mine": {
        "rule_holds_on_own_numbers": bool(rule_holds),
        "moved_over_1pct_among_nonnull_ratio": moved_nonnull,
        "moved_over_1pct_among_all_nonzero_kin_old": moved_all,
        "sign_flips": sign_flips,
        "n_ratio_null_negative_kin_old": n_null_ratio_neg_old,
        "provenance": PROVENANCE,
    },
    "audited": {
        "label": audited["label"],
        "sym_tol": summ["sym_tol"],
        "n_inertias_moved_over_1pct": summ["n_inertias_moved_over_1pct"],
    },
    "verdict": "QUALIFIED",
    "note": "the label follows from the numbers and 1e-6 is a safe tolerance (the measured "
    "parts sit at 1e-9 or below; float32 storage is the only source above round-off); but the "
    "symmetry half of the rule is an identity for any symmetric M, so it cannot fail on the "
    "stored fields (it tests storage, not the construction), and the over-1-percent bar is met "
    f"by {moved_nonnull} of {n_chan - n_null_ratio} positive-kin_old channels ({moved_all} of "
    f"{n_chan} counting the negative ones), so the moved count needs its denominator stated; "
    "the wording CONFIRMED is right as a re-confirmation: the antisymmetric form was named "
    "the antisymmetric variant at the M5.21.3 audit § C7.4 with the corrected conjugation "
    "catalog coded in m5_21_3_f_confirm.py, restated at M5.21.9 (a0_conj), flagged at M5.32 R1 "
    "and R14-0 verify; the R25/R26 record should carry that provenance instead of a "
    "2026-09-24 discovery",
}

wall = round(time.time() - T0, 1)
out = {
    "task": "M5.32 R26-1 adversarial audit",
    "audited_json": os.path.relpath(AUDITED, HERE),
    "audited_script_read": False,
    "fields": field_rows,
    "identity_probe": ident,
    "json_rederivations": json_checks,
    "claims": claims,
    "summary_table": [[k, v["verdict"]] for k, v in claims.items()],
    "wall_s": wall,
}
with open(OUT, "w") as f:
    json.dump(out, f, indent=1)

print()
print("| claim | mine | audited | verdict |")
print("| --- | --- | --- | --- |")
print(
    f"| A antisymmetric catalog / symmetric corrected | max sym {max_sym_mine:.1e}, "
    f"max antisym {max_asym_mine:.1e} (6 fields), own-vs-platform {worst_construct:.1e} | "
    f"{summ['max_old_symmetric_part']:.1e} / {summ['max_new_antisymmetric_part']:.1e} "
    f"({summ['n_ok']} fields) | {claims['A']['verdict']} |"
)
s1g = s1["gens"]
print(
    f"| B inertias move | S1 n32: "
    + ", ".join(f"{g} {s1g[g]['kin_old']:.2f} -> {s1g[g]['kin_new']:.2f}" for g in GENS)
    + f"; raw {s1['kin_raw_commutator_Jz']:.0f} vs {s1['kin_raw_anticommutator_Jz']:.0f}; "
    f"worst rel vs JSON {worst_kin:.1e}, own contraction vs kin_of {worst_own:.1e} | "
    f"1.56->1.24, 1.55->1.15, 1.72->1.84, 1.35->1.27, -0.34->-0.37, -0.57->-0.55; 14071 vs 37180 | "
    f"{claims['B']['verdict']} |"
)
print(
    f"| C reader census | {len(own_callers)} callers by own grep, {len(missing_from_census)} missing "
    f"from the census ({', '.join(missing_from_census)}), larmor misverdicted, "
    f"{len(kin_hits_no_kin)}/{len(kin_hits)} kin_read hits are noise; {spot_right} named verdicts right | "
    f"{len(census)} scripts, {len(census_callers)} flagged callers | {claims['C']['verdict']} |"
)
print(
    f"| D label + rule | rule holds on own numbers: {rule_holds}; moved {moved_nonnull}/"
    f"{n_chan - n_null_ratio} (nonnull) or {moved_all}/{n_chan} (all); identity: nonsymmetric M "
    f"fails at {ident['nonsymmetric_M']['max_old_symmetric_part']:.2f} | "
    f"{audited['label']}, tol {summ['sym_tol']}, moved {summ['n_inertias_moved_over_1pct']} | "
    f"{claims['D']['verdict']} |"
)
print(f"\nwall {wall} s; written {os.path.relpath(OUT, HERE)}")
