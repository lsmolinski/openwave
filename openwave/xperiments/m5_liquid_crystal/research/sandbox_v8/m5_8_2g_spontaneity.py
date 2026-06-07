"""
M5.8.2g — THE FIELD HANDOFF (Track C → 2d): the spontaneity test.

C3 (m5_8_2f3) DERIVED at the reduced-CC level: the dressed static minimum is
UN-SITTABLE — K_bb(a*) = −67.6 < 0 (direct quadrature) with U″ = +1723 > 0:
energetically minimal, kinetically ghost ⇒ compulsory spontaneous motion in
the AMPLITUDE/breathing channel — to which the s(t) clock-overlap probe is
BLIND (diagnosing the 2c-1 "spontaneous start machine-zero" null as probe
blindness, G-2c-3).

PREDICTION UNDER TEST: an UNKICKED dressed seed (Ṁ₀ = P₀ = 0 EXACTLY) under
the 2d quartic dynamics (β = 1.558, the validated stack) spontaneously grows
breathing motion — visible in dM_rms(t) (field departure), T(t) (kinetic
energy from zero), u_min(t) (the fuel dive) — while |s(t)| (the clock
overlap) stays comparatively small at first.

HONEST CAVEATS (both outcomes informative):
  · the constrained integrator's positive-inertia eigenprojection is exactly
    the machinery that may suppress or reshape the ghost-amplitude channel —
    the test probes whether the un-sittability survives (a) the per-mode
    projection and (b) full nonlinearity;
  · C3's K_bb < 0 is a NET (summed) inertia — the per-voxel A(Ṁ) spectrum the
    integrator sees is finer-grained; a null here would localize the
    reduction artifact, a growth confirms the mechanism at field level.

GATES (the 3-outcome map, judge the traces — no single-number cuts):
  S0  rig check: the standard kick (0.05) reproduces the 2d clock (|s| alive)
  S1  THE TEST: kick = 0, jitter = 0 — Mdf/Pf EXACTLY zero; growth must come
      from the seed's own truncation-level asymmetry (f32 ~1e-7)
  S2  jitter ladder (1e-6, 1e-4 symmetric random Ṁ₀): exponential instability
      ⇒ comparable growth RATES, onset shifted ~ln(ratio); stable ⇒ response
      amplitude ∝ jitter (ratio ~100)
  S3  blindness check: |s|/dM in the early window — the clock probe should
      LAG the amplitude departure if the C3 channel assignment is right

PREREQUISITE: `_m5_8_2cb_ref.npz` (CB_STEPS=900 python
m5_8_2cb_taichi_constrained.py ref). Imports the 2d module (single ti.init —
the 2e pattern; do NOT import 2e alongside).

RESULTS (2026-06-06 — SPONTANEITY CONFIRMED AT FIELD LEVEL, dt-converged):
  Batch 1 (S0-S3, 12000 steps each):
    S1 (Ṁ₀ = P₀ = 0 EXACT) reaches the SAME breathing state as the kicked
    control — dM 1.58 vs 1.46, T 80 vs 84, the same fuel dive — while the
    clock probe stays at |s| ~ 1e-18 MACHINE ZERO throughout: the BLINDNESS
    half of the C3 diagnosis is unconditional (the old G-2c-3 "spontaneous
    start machine-zero" null measured a symmetry-protected channel while the
    amplitude channel carried O(1) dynamics). Caveat found honestly:
    dM = 0.031 already at step 100 ⇒ the raw seed slides immediately (it was
    never a discrete equilibrium) — the jitter ladder cannot discriminate
    (ratio jit4/jit6 = 0.9 either way) ⇒ the settle protocol is required.
  S4 (damped settle, 24000 steps, γ = 1e-3/step): T drained 45.8 → 0.72
    (64×) with residual creep 1e-5/step — a quasi-settled discrete config
    (saved → _m5_8_2g_settled.npz). Note: the potential descent regenerated
    T at ~90% of the damping drain for 22k+ steps — the downhill never ended.
  S5 (restart from the settled config, P = 0 EXACT) — THE RESULT:
    T regrows 0 → 2.59 (τ=2000) → 5.76 (τ=4000) with H conserved to 0.5%
    and dM 6× above the residual-slide line, super-linear.
    ★ dt/2 CONFIRM AT MATCHED τ (the discipline that killed every artifact
    this session): dM to 4 SIGNIFICANT DIGITS (0.2500 vs 0.2501), T to 3
    (5.761 vs 5.762), H to 5 — the early-window growth is REAL DYNAMICS.
    ⇒ a settled configuration with exactly zero momentum spontaneously
    moves: THE UN-SITTABLE MINIMUM IS CONFIRMED AT FIELD LEVEL — C3's
    derivation validated; G-2c-3's spontaneous start DEMONSTRATED (in the
    amplitude/breathing channel, where it always lived).
  Honest residuals: (i) past τ ≈ 5000 the deep-floor stiffness cascade
    dominates (u_min → −212, H → 4e5 — the known explicit-stepper issue;
    late S5 numbers are numerical, the claim rests on the dt-converged
    window); (ii) "settled" = quasi-settled (T 0.72 residual) — but the
    dt-invariant exponential regrowth from P = 0 exact is not explicable by
    the 1e-5/step residual slope (linear dM ≈ 0.04 at τ=4000 vs measured
    0.25 with T from zero); (iii) the ω-attractor half of G-2c-3 (the
    breathing settling onto the 2e ω₀-comb) remains — the S1-reaches-the-
    same-state-as-S0 match (dM/T/u_min) is attractor-consistent evidence,
    not yet a frequency measurement.

USAGE:  python m5_8_2g_spontaneity.py [steps]                (default 12000)
        python m5_8_2g_spontaneity.py settle [s4] [s5]       (24000 16000)
        python m5_8_2g_spontaneity.py restart <dts> <steps>  (0.5 8000)
"""
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[4]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from openwave.xperiments.m5_liquid_crystal.research.sandbox_v8 import (  # noqa: E402
    m5_8_2d_quartic_saturation as d2,
)
from openwave.xperiments.m5_liquid_crystal.research.sandbox_v8.m5_8_2c1_full_evolution import (  # noqa: E402
    central, tw,
)

BETA = 1.558                       # the 2d winner (2β|u_min(seed)| ≈ 3)
PROBE = 100
OUT_NPZ = HERE / "_m5_8_2g_traces.npz"

ti = d2.ti


@ti.kernel
def k_damp(g: ti.f32):
    for i, j, k in d2.Pf:
        d2.Pf[i, j, k] *= g


def run_spont(tag, seed, dt, steps, kick=0.0, jitter=0.0, jseed=1,
              damp=0.0, M_start=None, M_ref=None):
    """One evolution; kick=jitter=0 ⇒ Ṁ₀ = P₀ EXACTLY zero (the test arm).

    damp > 0: momentum damping P ← (1−damp-rate)·P per step — the SETTLE arm
    (the S4/S5 protocol: damped flow → a discrete V-minimum, then restart).
    M_start overrides the initial field; M_ref is the dM baseline (defaults
    to M_start — for S5 the departure is measured vs the SETTLED config).
    Probes (every PROBE steps): dM_rms (field departure, the breathing-
    amplitude proxy), T (kinetic), u_min (fuel), H, |s| (clock overlap).
    Returns (traces, M_final)."""
    M0, Mth, act, core, rhat, h = seed
    if M_start is not None:
        M0 = M_start
    if M_ref is None:
        M_ref = M0
    n_act = int(act.sum())
    actb = act > 0.5
    d2.Mf.from_numpy(M0.astype(np.float32))
    d2.Mdf.fill(0.0)
    d2.Pf.fill(0.0)
    d2.actf.from_numpy(act.astype(np.int32))
    d2.Bas.from_numpy(d2.SYM_BASIS.astype(np.float32))
    inv2h = 1.0 / (2.0 * h)
    Md0 = np.zeros_like(M0)
    if kick > 0.0:
        Md0 += kick * Mth
    if jitter > 0.0:
        rng = np.random.default_rng(jseed)
        R = rng.standard_normal(M0.shape)
        Md0 += jitter * 0.5 * (R + np.swapaxes(R, -1, -2))
    if kick > 0.0 or jitter > 0.0:
        Md0 = Md0 * act[..., None, None]
        d2.Mdf.from_numpy(Md0.astype(np.float32))
        Mi = [central(M0, ax, h) for ax in range(3)]
        P0 = np.zeros_like(M0)
        for i_ in range(3):
            P0 += d2.np_commf(tw(d2.np_commf(Md0, Mi[i_])), Mi[i_])
        d2.Pf.from_numpy((4.0 * P0).astype(np.float32))
    tr = dict(t=[], dM=[], T=[], umin=[], H=[], s=[])
    t0 = time.time()
    for n_ in range(steps):
        d2.k_flux(inv2h, BETA)
        d2.k_force(inv2h, dt)
        d2.red.fill(0.0)
        d2.k_clamp_sum()
        r = d2.red.to_numpy()
        d2.k_clamp_apply(float(r[0]) / n_act, float(r[1]) / n_act,
                         float(r[2]) / n_act)
        if damp > 0.0:
            k_damp(1.0 - damp)
        d2.k_solve(inv2h)
        d2.k_update(dt)
        if n_ % PROBE == PROBE - 1:
            M = d2.Mf.to_numpy().astype(np.float64)
            Md = d2.Mdf.to_numpy().astype(np.float64)
            u = d2.uf.to_numpy().astype(np.float64)
            dM = float(np.sqrt(np.einsum("...ab,...ab->...",
                                         M - M_ref, M - M_ref)[actb].mean()))
            Mi = [central(M, ax, h) for ax in range(3)]
            T = 0.0
            for i_ in range(3):
                F0 = d2.np_commf(Md, Mi[i_])
                T = T + 2.0 * np.einsum("...ab,...ab->...", F0, tw(F0))
            Tt = float(T[actb].sum()) * h**3
            H = float((T + u + BETA * u * u)[actb].sum()) * h**3
            s = float(np.einsum("...ab,...ab->...", M - M_ref, Mth)[core].mean())
            tr["t"].append(n_ + 1)
            tr["dM"].append(dM)
            tr["T"].append(Tt)
            tr["umin"].append(float(u[actb].min()))
            tr["H"].append(H)
            tr["s"].append(abs(s))
            if n_ % 2000 == 1999:
                print(f"   [{tag}] step {n_ + 1:6d} [{time.time() - t0:4.0f}s]"
                      f"  dM={dM:.3e}  T={Tt:.3e}  u_min={tr['umin'][-1]:+.4f}"
                      f"  |s|={abs(s):.2e}  H={H:.4f}")
            if not np.isfinite(H):
                print(f"   [{tag}] NON-FINITE at step {n_ + 1}")
                break
    return ({k: np.array(v) for k, v in tr.items()},
            d2.Mf.to_numpy().astype(np.float64))


def grow_rate(t, y):
    """Log-slope of y(t) per 1000 steps over the run's thirds (trend style)."""
    y = np.maximum(np.asarray(y, float), 1e-300)
    n = len(y)
    out = []
    for a, b in ((0, n // 3), (n // 3, 2 * n // 3), (2 * n // 3, n)):
        if b - a > 2:
            sl = np.polyfit(t[a:b], np.log(y[a:b]), 1)[0] * 1000.0
            out.append(sl)
    return out


def settle_main():
    """S4 + S5 — THE DECISIVE DISCRIMINATOR: transient slosh vs un-sittable.

    The S1 seed was never an exact static solution of the DISCRETE equations,
    so its growth could be ordinary off-equilibrium relaxation. S4 damps the
    momentum (P ← γP per step) until the field settles at a genuine discrete
    V-minimum; S5 restarts from THAT config with P = 0 exactly. Growth from a
    settled configuration cannot be a transient: it is the un-sittable
    minimum at field level — or, if S5 sits still, the C3 net-K_bb<0 is a
    reduction artifact (decisive either way). S4's own behavior is also
    diagnostic: if the damped flow never kills T, no static attractor exists."""
    steps4 = int(sys.argv[2]) if len(sys.argv) > 2 else 24000
    steps5 = int(sys.argv[3]) if len(sys.argv) > 3 else 16000
    print("=" * 78)
    print("M5.8.2g settle — S4 damped settle + S5 restart (the transient-vs-"
          "un-sittable discriminator)")
    print(f"  β = {BETA}, dt = {d2.DT}, S4 {steps4} damped steps (γ = 0.999)"
          f" + S5 {steps5} undamped from the settled config, P = 0 exact")
    print("=" * 78)
    seed = d2.load_seed()
    print("\n[S4 — damped settle]")
    tr4, M_set = run_spont("S4-damped", seed, d2.DT, steps4, damp=1e-3)
    print(f"    settle quality: T {tr4['T'][0]:.3e} → {tr4['T'][-1]:.3e}"
          f"   dM-rate(end) = "
          f"{(tr4['dM'][-1] - tr4['dM'][-10]) / (9 * PROBE):.3e}/step"
          f"   u_min(end) = {tr4['umin'][-1]:+.4f}")
    np.savez(HERE / "_m5_8_2g_settled.npz", M_settled=M_set,
             T_end=tr4["T"][-1], steps=steps4)
    print("\n[S5 — restart from the settled config, P = 0 exact]")
    tr5, _ = run_spont("S5-restart", seed, d2.DT, steps5,
                       M_start=M_set, M_ref=M_set)
    print("\n" + "=" * 78)
    print("S5 VERDICT (judge the traces):")
    sl = grow_rate(tr5["t"], np.maximum(tr5["dM"], 1e-12))
    print(f"  dM: {tr5['dM'][0]:.3e} → {tr5['dM'][-1]:.3e}   "
          f"T: {tr5['T'][0]:.3e} → {tr5['T'][-1]:.3e}   "
          f"log-slope dM /1k by thirds: "
          + "  ".join(f"{v:+.3f}" for v in sl))
    grew = tr5["dM"][-1] > 100 * max(tr5["dM"][0], 1e-12) and tr5["dM"][-1] > 1e-3
    if grew:
        print("  → GROWTH FROM A SETTLED CONFIG: the un-sittable minimum is"
              " CONFIRMED AT FIELD LEVEL — no transient explanation survives"
              " (G-2c-3 spontaneity demonstrated; C3's derivation validated).")
    else:
        print("  → the settled config SITS: S1's growth was off-equilibrium"
              " transient; C3's net-K_bb<0 does not bind the discrete field"
              " minimum (a reduction artifact at field level — document).")
    z = dict(np.load(OUT_NPZ)) if OUT_NPZ.exists() else {}
    z.update({f"S4_{q}": v for q, v in tr4.items()})
    z.update({f"S5_{q}": v for q, v in tr5.items()})
    np.savez(OUT_NPZ, **z)
    print(f"  traces appended → {OUT_NPZ.name};  settled config → "
          "_m5_8_2g_settled.npz")
    print("=" * 78)
    return 0


def restart_main():
    """S5 dt-variant from the SAVED settled config — the dt-invariance
    discriminator (the 2c discipline): real dynamics reproduces at matched τ;
    stepper-driven growth shifts with dt.
    USAGE: ... restart <dt_scale> <steps>     (e.g. restart 0.5 8000)"""
    dts = float(sys.argv[2]) if len(sys.argv) > 2 else 0.5
    steps = int(sys.argv[3]) if len(sys.argv) > 3 else 8000
    z = np.load(HERE / "_m5_8_2g_settled.npz")
    M_set = z["M_settled"]
    print("=" * 78)
    print(f"M5.8.2g restart — S5 dt-variant: dt × {dts}, {steps} steps "
          f"(τ ≡ {int(steps * dts)} full-dt steps), P = 0 exact, from the"
          " saved settled config")
    print("=" * 78)
    seed = d2.load_seed()
    tr, _ = run_spont(f"S5-dt{dts}", seed, d2.DT * dts, steps,
                      M_start=M_set, M_ref=M_set)
    print("\n  matched-τ comparison (read against the full-dt S5 trace at"
          " τ/2 of these step counts):")
    for frac in (0.25, 0.5, 0.75, 1.0):
        i = min(int(len(tr["t"]) * frac) - 1, len(tr["t"]) - 1)
        print(f"    τ ≡ {int(tr['t'][i] * dts):6d} full-dt steps:  dM ="
              f" {tr['dM'][i]:.3e}   T = {tr['T'][i]:.3e}   H = {tr['H'][i]:.4f}")
    return 0


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "settle":
        return settle_main()
    if len(sys.argv) > 1 and sys.argv[1] == "restart":
        return restart_main()
    steps = int(sys.argv[1]) if len(sys.argv) > 1 else 12000
    print("=" * 78)
    print("M5.8.2g — THE FIELD HANDOFF: the spontaneity test (unkicked dressed"
          " seed)")
    print(f"  β = {BETA} (the 2d winner), dt = {d2.DT}, {steps} steps/run,"
          f" grid {d2.N}³ Metal f32")
    print("  C3 prediction: spontaneous AMPLITUDE growth (dM, T, u_min) with"
          " the s(t) clock probe lagging")
    print("=" * 78)
    seed = d2.load_seed()
    runs = {}
    for tag, kick, jit in (("S0-control k=0.05", 0.05, 0.0),
                           ("S1-TEST k=0 j=0", 0.0, 0.0),
                           ("S2-jit 1e-6", 0.0, 1e-6),
                           ("S2-jit 1e-4", 0.0, 1e-4)):
        print(f"\n[{tag}]")
        runs[tag], _ = run_spont(tag, seed, d2.DT, steps, kick=kick, jitter=jit)

    print("\n" + "=" * 78)
    print("VERDICT TABLE (judge the data — the trend_report discipline):")
    print(f"{'run':22s} {'dM(first)':>10s} {'dM(last)':>10s} {'T(last)':>10s}"
          f" {'u_min(min)':>11s} {'|s|(last)':>10s}  log-slope dM /1k (thirds)")
    for tag, tr in runs.items():
        sl = grow_rate(tr["t"], tr["dM"])
        print(f"{tag:22s} {tr['dM'][0]:10.3e} {tr['dM'][-1]:10.3e}"
              f" {tr['T'][-1]:10.3e} {tr['umin'].min():11.4f}"
              f" {tr['s'][-1]:10.2e}  " + "  ".join(f"{v:+.3f}" for v in sl))
    s1 = runs["S1-TEST k=0 j=0"]
    j6 = runs["S2-jit 1e-6"]
    j4 = runs["S2-jit 1e-4"]
    print("\n  S2 scaling discriminator: dM(last) ratio jit4/jit6 = "
          f"{j4['dM'][-1] / max(j6['dM'][-1], 1e-300):.1f}"
          "   (≈100 ⇒ LINEAR/stable response; ≈1 ⇒ saturated exponential"
          " instability; in-between ⇒ growing, not yet saturated)")
    third = max(len(s1["t"]) // 3, 1)
    print(f"  S3 blindness (early third of the TEST run): |s|/dM = "
          f"{s1['s'][:third].mean() / max(s1['dM'][:third].mean(), 1e-300):.3e}"
          "   (≪1 supports the amplitude-channel assignment)")
    np.savez(OUT_NPZ, **{f"{k}_{q}": v for k, tr in runs.items()
                         for q, v in tr.items()},
             beta=BETA, steps=steps)
    print(f"  traces saved → {OUT_NPZ.name}")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
