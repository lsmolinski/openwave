"""
M5.8.2m — N-6a: THE UNIT-FREE ZBW LAW — `ω ∝ H_rest` across a seed-structure
mass family. (v2 — the R_W family; the RC attempt was a measured no-op.)

THE TEST (the strongest pre-calibration ZBW statement, promoted 2026-06-07):
`ω = 2mc²/ℏ` is unit-bound as an absolute claim (the 2j insight), but the
LAW behind it — clock frequency PROPORTIONAL to energy — is unit-FREE
across a family of states with different rest energies: constant
ω₁/(2H) ⇒ the law holds pre-calibration; a clean exponent ≠ 1 in ω ∝ H^α
is the honest alternative finding.

★ THE v1 MISFIRE — RECORDED HONESTLY (2026-06-07, the audit that caught it):
  v1 used the core radius RC as the family knob. The fixed-domain audit
  found H_fixed-domain IDENTICAL to 5 digits (16.740) across "members" —
  **RC never enters the seed field**: in build_grid it only regularizes the
  r = 0 division (no voxel sits at r = 0) and normalize() erases it; RC was
  only ever load-bearing through the MASKS. The three v1 runs were the SAME
  physical field with different clamp/probe domains. Two lessons + a bonus:
  · LESSON 1 (institutionalized below as the KNOB GATE): validate that a
    family knob actually changes the seed (H_fixed-domain spread > 5%)
    BEFORE burning GPU time;
  · LESSON 2: ω extracted from different PROBE REGIONS of the same state
    differs strongly (act(0.8) → 1.20; r>2.2 → 1.70; r>2.8 → 2.03-stable
    — the outer shell reads the 2ω₁ harmonic as dominant): spectral
    content is spatially structured; probe-region must be FIXED across any
    comparison (it is, below);
  · the v1 npz files (_m5_8_2m_RC*.npz) are kept as the record of that
    probe-region finding.

THE v2 FAMILY: R_W ∈ {2.5, 3.5, 4.5} — the boost-dressing width, which
DOES enter the seed (`seed_M: w = exp(−(r/R_W)²)`; the dressed defect's
rest energy varies with its dressing cloud). RC stays 0.8 ⇒ the act mask
and probe region are IDENTICAL for all members (no v1 confound). Per
member: β anchored (2d rule), S₁/S₂ saved (any-β H reconstructable), the
free spontaneous dense run (Ṁ₀ = P₀ = 0 exact, dt = DT/2, t = 48), ω₁ via
the N-1 pipeline, AND H(t)/T(t) probes every 200 steps — the law is tested
against BOTH H_static (the seed rest energy) and H_dyn (the breathing
state's own energy — closer to the ZBW object; the excitation confound
stated, not hidden). R_W is a module constant upstream ⇒ subprocess per
member via the M58_RW env.

RESULTS (2026-06-07 — N-6a v2 COMPLETE: ω IS RIGID — the naive law FAILS,
informatively):
  | R_W | H_static | H_dyn | ω₁ |        exponents: ω ∝ H_static^0.033,
  | 2.5 | 11.01 | 23.8 | 1.152 |                    ω ∝ H_dyn^0.024
  | 3.5 | 16.74 | 43.9 | 1.188 |        (the ZBW law ⇒ 1; measured ≈ 0)
  | 4.5 | 29.13 | 93.7 | 1.191 |
  · ω₁ is CONSTANT within one FFT bin across a 2.6× rest-energy / 3.9×
    state-energy family ⇒ the breathing frequency does NOT track the
    dressing energy — the naive ω ∝ H law FAILS on this family, decisively.
  · THE REFRAME (both honest): (i) frequency RIGIDITY against a large
    energy variation is the DTC robustness property — ω₁ is a property of
    the defect CORE (held fixed across this family: the bare hedgehog
    frame is scale-free and R_W-independent), extending N-1's attractor
    statement from "start-independent" to "dressing-energy-independent";
    (ii) the TRUE ZBW mass-family test must vary the CORE — unavailable
    in the V = 0 stack (no intrinsic core length; the frame is scale-free)
    ⇒ it requires the V-on/Faber-r₀ stack (the resumed M5.9/NG track).
    Scope recorded, not hidden.
  · β anchoring responded per member (1.94/1.56/1.55) — the seeds differ
    structurally (knob gate 177.8% H spread); identical probe regions
    (RC = 0.8 fixed) — the v1 confound excluded by design.

USAGE:  python m5_8_2m_zbw_law.py run            (knob gate + the family)
        python m5_8_2m_zbw_law.py one [steps]    (current env R_W; internal)
        python m5_8_2m_zbw_law.py analyze        (the law table)
"""
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
V8 = HERE.parent / "sandbox_v8"
REPO_ROOT = HERE.parents[4]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from openwave.xperiments.m5_liquid_crystal.research.sandbox_v8.m5_8_2c1_full_evolution import (  # noqa: E402
    DT, R_W, build_grid, central, tw,
)
from openwave.xperiments.m5_liquid_crystal.research.sandbox_vn.m5_8_2h_omega_attractor import (  # noqa: E402
    np_commf, spec, peaks,
)

FAMILY = (2.5, 3.5, 4.5)         # R_W members (3.5 = the validated stack)
DENSE = 20
HPROBE = 200
KNOB_GATE_MIN_SPREAD = 0.05      # H_fixed must vary ≥5% across members


def out_npz(rw):
    return HERE / "data" / f"_m5_8_2m_RW{rw:g}.npz"


def ref_npz(rw):
    return V8 / ("_m5_8_2cb_ref.npz" if rw == 3.5
                 else f"_m5_8_2cb_ref_RW{rw:g}.npz")


def u_density(M, h):
    Mi = [central(M, ax, h) for ax in range(3)]
    u = 0.0
    for i in range(3):
        for j in range(i + 1, 3):
            F = np_commf(Mi[i], Mi[j])
            u = u + 2.0 * np.einsum("...ab,...ab->...", F, tw(F))
    return u


def omega1_from(p2, umin, dsamp):
    umin = np.asarray(umin, float)
    bad = np.where(~np.isfinite(umin) | (umin < -20.0))[0]
    ne = int(bad[0]) if len(bad) else len(umin)
    y = np.asarray(p2)[:ne]
    tops = []
    for T in (2.5, 5.0, 10.0):
        n = max(int(T / dsamp) | 1, 3)
        pf = peaks(*spec(y, dsamp, detrend_n=n), 1.0, np.inf, n=1)
        if pf:
            tops.append(pf[0][0])
    return float(np.median(tops)), tops, ne


def one_main():
    from openwave.xperiments.m5_liquid_crystal.research.sandbox_v8 import (
        m5_8_2d_quartic_saturation as d2,
    )
    steps = int(sys.argv[2]) if len(sys.argv) > 2 else 48000
    dt = DT / 2
    M0, Mth, act, core, _, h = d2.load_seed()
    actb = act > 0.5
    n_act = int(act.sum())
    g = build_grid()
    u0 = d2.np_u_density([central(M0, ax, h) for ax in range(3)])
    inter3 = np.zeros(g["r"].shape, bool)
    inter3[3:-3, 3:-3, 3:-3] = True
    umin0 = float(u0[inter3].min())
    beta = 3.0 / (2.0 * abs(umin0))
    S1 = float(u0[actb].sum()) * h ** 3
    S2 = float((u0 * u0)[actb].sum()) * h ** 3
    H_static = S1 + beta * S2
    print(f"[R_W={R_W:g}] seed u_min={umin0:+.4f} → β={beta:.4f};"
          f" H_static={H_static:.4f} (S1={S1:.4f}, S2={S2:.4f});"
          f" n_act={n_act}; {steps} steps dt={dt}")
    rng = np.random.default_rng(7)
    R = rng.standard_normal(M0.shape)
    R = 0.5 * (R + np.swapaxes(R, -1, -2)) * act[..., None, None]
    R /= np.sqrt(np.einsum("...ab,...ab->...", R, R)[actb].mean())
    d2.Mf.from_numpy(M0.astype(np.float32))
    d2.Mdf.fill(0.0)
    d2.Pf.fill(0.0)                                      # FREE: P = 0 exact
    d2.actf.from_numpy(act.astype(np.int32))
    d2.Bas.from_numpy(d2.SYM_BASIS.astype(np.float32))
    inv2h = 1.0 / (2.0 * h)
    p2s, p3s, dMs, umins = [], [], [], []
    hts, hTs, hHs = [], [], []
    t0 = time.time()
    for n_ in range(steps):
        d2.k_flux(inv2h, beta)
        d2.k_force(inv2h, dt)
        d2.red.fill(0.0)
        d2.k_clamp_sum()
        r = d2.red.to_numpy()
        d2.k_clamp_apply(float(r[0]) / n_act, float(r[1]) / n_act,
                         float(r[2]) / n_act)
        d2.k_solve(inv2h)
        d2.k_update(dt)
        if n_ % DENSE == DENSE - 1:
            M = d2.Mf.to_numpy().astype(np.float64)
            u = d2.uf.to_numpy().astype(np.float64)
            dMv = M - M0
            p2s.append(float(np.einsum("...ab,...ab->...", dMv, M0)[actb].mean()))
            p3s.append(float(np.einsum("...ab,...ab->...", dMv, R)[actb].mean()))
            dMs.append(float(np.sqrt(np.einsum(
                "...ab,...ab->...", dMv, dMv)[actb].mean())))
            umins.append(float(u[actb].min()))
            if not np.isfinite(dMs[-1]):
                print(f"   NON-FINITE at step {n_ + 1}")
                break
        if n_ % HPROBE == HPROBE - 1:
            M = d2.Mf.to_numpy().astype(np.float64)
            Md = d2.Mdf.to_numpy().astype(np.float64)
            u = d2.uf.to_numpy().astype(np.float64)
            Mi = [central(M, ax, h) for ax in range(3)]
            T = 0.0
            for i_ in range(3):
                F0 = np_commf(Md, Mi[i_])
                T = T + 2.0 * np.einsum("...ab,...ab->...", F0, tw(F0))
            hts.append(n_ + 1)
            hTs.append(float(T[actb].sum()) * h ** 3)
            hHs.append(float((T + u + beta * u * u)[actb].sum()) * h ** 3)
        if n_ % 8000 == 7999:
            print(f"   [R_W={R_W:g}] step {n_ + 1:6d} [{time.time() - t0:4.0f}s]"
                  f"  dM={dMs[-1]:.3e}  T={hTs[-1]:.3e}  H={hHs[-1]:.4f}"
                  f"  u_min={umins[-1]:+.3f}")
    np.savez(out_npz(R_W), rw=R_W, beta=beta, umin_seed=umin0, S1=S1, S2=S2,
             H_static=H_static, omega1=np.nan, omega_dial=np.array([]),
             p2=np.array(p2s), p3=np.array(p3s), dM=np.array(dMs),
             umin=np.array(umins), ht=np.array(hts), hT=np.array(hTs),
             hH=np.array(hHs), dt=dt, dense=DENSE, ne=-1)
    print(f"   raw traces saved → {out_npz(R_W).name} (analysis next)")
    om, tops, ne = omega1_from(p2s, umins, DENSE * dt)
    print(f"   [R_W={R_W:g}] ω₁ = {om:.3f} (dial: "
          + ", ".join(f"{v:.3f}" for v in tops)
          + f"; healthy {ne}/{len(umins)})  →  ω₁/(2H_static) ="
          f" {om / (2 * H_static):.5f}")
    np.savez(out_npz(R_W), rw=R_W, beta=beta, umin_seed=umin0, S1=S1, S2=S2,
             H_static=H_static, omega1=om, omega_dial=np.array(tops),
             p2=np.array(p2s), p3=np.array(p3s), dM=np.array(dMs),
             umin=np.array(umins), ht=np.array(hts), hT=np.array(hTs),
             hH=np.array(hHs), dt=dt, dense=DENSE, ne=ne)
    print(f"   saved (with ω) → {out_npz(R_W).name}")
    return 0


def run_main():
    env0 = dict(os.environ)
    # ── THE KNOB GATE (lesson 1): seeds first, verify H_fixed actually varies
    print("[knob gate] generating/validating the family seeds")
    Hf = {}
    g = build_grid()
    inter = np.zeros(g["r"].shape, bool)
    inter[2:-2, 2:-2, 2:-2] = True
    fixed = inter & (g["r"] > 1.6) & (g["rho"] > 0.8)
    for rw in FAMILY:
        env = dict(env0, M58_RW=f"{rw}")
        if not ref_npz(rw).exists():
            print(f"  [seed] R_W={rw:g} → {ref_npz(rw).name}")
            subprocess.run([sys.executable, "-u",
                            str(V8 / "m5_8_2cb_taichi_constrained.py"), "ref"],
                           env=dict(env, CB_STEPS="2"), check=True)
        d = np.load(ref_npz(rw))
        M0, h = d["M_seed"], float(d["h"])
        u0 = u_density(M0, h)
        Hf[rw] = float((u0 + 1.5584 * u0 * u0)[fixed].sum()) * h ** 3
        print(f"  R_W={rw:g}: H_fixed-domain = {Hf[rw]:.4f}")
    spread = max(Hf.values()) / min(Hf.values()) - 1
    print(f"  knob-gate spread: {100 * spread:.1f}%"
          f"  (require ≥ {100 * KNOB_GATE_MIN_SPREAD:.0f}%)")
    if spread < KNOB_GATE_MIN_SPREAD:
        print("  ❌ KNOB GATE FAILED — the knob does not reach the seed;"
              " NOT burning GPU time (the v1 lesson). Aborting.")
        return 1
    print("  ✅ knob gate PASSED — the family is physical\n")
    for rw in FAMILY:
        if out_npz(rw).exists():
            print(f"[skip] R_W={rw:g} — {out_npz(rw).name} exists")
            continue
        print(f"[run] R_W={rw:g}")
        subprocess.run([sys.executable, "-u", __file__, "one"],
                       env=dict(env0, M58_RW=f"{rw}"), check=True)
    return analyze_main()


def analyze_main():
    print("\n" + "=" * 78)
    print("[N-6a v2] THE UNIT-FREE ZBW LAW — ω vs H across the R_W dressing"
          " family (24³, fixed probe region)")
    print("=" * 78)
    rows = []
    for rw in FAMILY:
        if not out_npz(rw).exists():
            print(f"  R_W={rw:g}: no data yet")
            continue
        z = np.load(out_npz(rw))
        om, dial = float(z["omega1"]), z["omega_dial"]
        if not np.isfinite(om):
            om, dial, _ = omega1_from(z["p2"], z["umin"],
                                      float(z["dense"]) * float(z["dt"]))
        umin = np.asarray(z["umin"], float)
        bad = np.where(~np.isfinite(umin) | (umin < -20.0))[0]
        ne = int(bad[0]) if len(bad) else len(umin)
        t_end = (np.arange(len(umin))[ne - 1] + 1) * 20
        Hwin = z["hH"][z["ht"] <= t_end * 1]
        H_dyn = float(np.median(Hwin[:max(len(Hwin) // 3, 1)]))
        rows.append((rw, float(z["beta"]), float(z["H_static"]), H_dyn,
                     om, np.asarray(dial)))
    print("| R_W | β | H_static | H_dyn (early) | ω₁ (dial) |"
          " ω₁/(2H_static) | ω₁/(2H_dyn) |")
    print("| --- | --- | --- | --- | --- | --- | --- |")
    for rw, b, Hs, Hd, om, dial in rows:
        print(f"| {rw:g} | {b:.4f} | {Hs:.3f} | {Hd:+.3f} |"
              f" {om:.3f} ({', '.join(f'{v:.2f}' for v in dial)}) |"
              f" {om / (2 * Hs):.5f} | {om / (2 * abs(Hd)):.5f} |")
    if len(rows) >= 2:
        for label, Hcol in (("H_static", 2), ("H_dyn", 3)):
            lH = np.log([abs(r[Hcol]) for r in rows])
            lo = np.log([r[4] for r in rows])
            if np.ptp(lH) > 1e-6:
                a = np.polyfit(lH, lo, 1)[0]
                print(f"  exponent vs {label}:  ω ∝ {label}^{a:.3f}"
                      f"   (ZBW LAW ⇒ 1)")
    print("\n  read: exponent ≈ 1 on either H ⇒ the law holds; a clean ≠ 1"
          " exponent is the honest")
    print("  alternative; the excitation confound (uncontrolled per member)"
          " is stated, not hidden.")
    print("=" * 78)
    return 0


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "analyze"
    if mode == "run":
        return run_main()
    if mode == "one":
        return one_main()
    if mode == "analyze":
        return analyze_main()
    print(__doc__.split("USAGE:")[1])
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
