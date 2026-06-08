"""
M5.8.2o — N-6c: ω(E) AND MOLTEN-NESS(E) — walking the clock down the
excitation axis toward the ground state.

TWO LOAD-BEARING QUESTIONS RIDE ON THIS RUN (made sharp by N-6d + N-6a):
  1. ZBW: the ratio test belongs to the GROUND clock — the saturated
     breather carries ~1.8× its rest energy as excitation. Does ω₁(E)
     extrapolate cleanly as E → H_rest, and to WHAT? (The linear 2b-2
     class was 5.86; the saturated state runs 1.1–1.2; N-6a showed ω is
     RIGID against dressing energy — is it also rigid against EXCITATION,
     or does it soften the way 5.86 → 1.1 suggests?)
  2. THERMAL (the molten-clock bridge): N-6d found molten-ness GROWS with
     excitation (D₂ 1.70 low-E vs 2.7–3.0 high-E). Does the clock
     REGULARIZE toward the ground state (λ → bias floor, D₂ → ~1)?
     YES ⇒ a clean rest clock + a thermal sector when excited (both the
     ZBW identification and the SABER heat picture close); NO ⇒ honest
     trouble for de Broglie-phase coherence — reported either way.

DESIGN: the standard seed (R_W = 3.5, β = 1.558), kicked at
kick ∈ {0.0125, 0.025} (excitation ∝ kick² ⇒ ~1/16 and ~1/4 of the
validated 0.05 arm) — dt = DT/2, t = 48, dense p2/umin + H/T probes (the
2m pattern). The HIGH-E points come free from the saved 2h arms
(R0-kicked 0.05, R2-seed, R3-jit) and the 2m RW3.5 member; every point
goes through the IDENTICAL ω pipeline + the IDENTICAL 2n chaos battery.

RESULTS (2026-06-07 — the weak-kick phase + the redesign it forced):
  ★ THE KICK KNOB DOES NOT CONTROL EXCITATION — the excitation FLOOR is
    set by the seed's own spontaneous release (~1.7× H_rest): the weak-kick
    arms landed at the SAME-or-higher H_dyn as the 0.05 control (49.2/52.4
    vs 47.6). The un-sittable minimum seen from a new angle, and a real
    physics statement: a cold defect CANNOT be prepared by gentle kicking —
    cooling requires DAMPING (the settle protocol). Itself report-worthy.
  Within the accessible high-E band (excitation 1.7–2.1× H_rest):
    ω₁ scatters 1.09–1.60 with NO monotone E-trend (the dial-tight 1.597
    at kick 0.0125 reads as MODE COMPOSITION, not energy — consistent with
    N-6a's rigidity); molten metrics uniformly high (D₂ 2.3–3.0,
    λ +0.5–0.75) — no regularization visible because NO arm is cold.
  ⇒ the ground-clock measurement ran through the SETTLE route — RESULTS:
  ★ THE CLOCK REGULARIZES TOWARD THE GROUND STATE (the molten bridge
    closes in the right direction). Deep S4 (48000 damped steps): T
    drained 45.8 → 0.21 (217×, H → −165 still slowly descending — the
    downhill never ends, per 2g). Restart, P = 0 exact, thirds of the
    regrowth (E grows along the run — the un-sittable axis sampled in
    time):
    | window  | T (median) | ω₁ (dial)            | D₂    | λ_max  |
    | third 1 | 2.36 COLD  | 1.583 (3.20*,1.51,1.58) | 1.68 | +0.110 |
    | third 2 | ~1100 heated | 2.731 (dial-tight)  | 15.6† | −0.059 |
    | third 3 | ~1.5e5 CASCADE | 13.60 (dial-tight) | 7.4† | +0.335 |
    (* the 2ω rank-flip; † estimator breakdown on nonstationary/cascade
    data — thirds 2-3 are regime-labeled, not load-bearing)
  · THE COLD ROW IS THE FINDING: λ_max +0.110 (vs +0.5–0.75 hot; the
    800-sample bias floor is ABOVE the 2400-sample +0.06 floor — i.e.
    consistent-with-REGULAR) and D₂ 1.68 (vs 2.7–3.0 hot; two-tone
    control 1.15; R1's transit state 1.70) at T = 2.36 — the coldest
    state this program has measured. Cold defect ⇒ near-regular clock;
    hot defect ⇒ molten. The melting curve now has its cold anchor.
  · FULL coldness is FORBIDDEN by the mechanism itself (un-sittability:
    T regrows immediately) — "the ground clock" is the T → 0 limit of a
    trend that is now measured to bend the right way. A feature, not a
    bug: the clock self-starts BECAUSE it cannot be fully cold.
  · ω along the regrowth: 1.58 (cold) → 2.73 (heated transient) →
    [13.6 numerical]; the saturated attractor sits at 1.09–1.15 (2h) —
    ω is STATE-dependent (mode composition), not simply E-monotone;
    the N-6a rigidity statement applies within the saturated class.
  · ZBW absolute with the COLD ω (1.58): 7.3e19 rad/s — still ~21×
    short; the structural-gap conclusion (N-6b) unchanged.

USAGE:  python m5_8_2o_omega_of_E.py run [steps]     (the 2 weak-kick arms)
        python m5_8_2o_omega_of_E.py settle [steps]  (deep S4 damping)
        python m5_8_2o_omega_of_E.py restart [steps] (the ground-clock readout)
        python m5_8_2o_omega_of_E.py analyze         (the combined table)
"""
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
    DT, central, tw,
)
from openwave.xperiments.m5_liquid_crystal.research.sandbox_vn.m5_8_2h_omega_attractor import (  # noqa: E402
    np_commf, spec, peaks, moving_mean,
)
from openwave.xperiments.m5_liquid_crystal.research.sandbox_vn.m5_8_2n_chaos_battery import (  # noqa: E402
    corr_dim, lyap_rosenstein,
)

KICKS = (0.0125, 0.025)
BETA = 1.558
DENSE = 20
HPROBE = 200
H_STATIC = 16.74                 # the validated member's rest energy (2j)


def out_npz(k):
    return HERE / f"_m5_8_2o_k{k:g}.npz"


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


def molten_metrics(p2, umin, dsamp):
    """The 2n battery's load-bearing pair (detrended D₂ + λ) on the healthy window."""
    umin = np.asarray(umin, float)
    bad = np.where(~np.isfinite(umin) | (umin < -20.0))[0]
    ne = int(bad[0]) if len(bad) else len(umin)
    y = np.asarray(p2)[:ne]
    det = y - moving_mean(y, max(int(5.0 / dsamp) | 1, 3))
    D2, _ = corr_dim(det)
    lam = lyap_rosenstein(det)
    return D2, lam


def run_main():
    from openwave.xperiments.m5_liquid_crystal.research.sandbox_v8 import (
        m5_8_2d_quartic_saturation as d2,
    )
    steps = int(sys.argv[2]) if len(sys.argv) > 2 else 48000
    dt = DT / 2
    M0, Mth, act, core, _, h = d2.load_seed()
    actb = act > 0.5
    n_act = int(act.sum())
    inv2h = 1.0 / (2.0 * h)
    for kick in KICKS:
        if out_npz(kick).exists():
            print(f"[skip] kick={kick:g} exists")
            continue
        print(f"\n[kick = {kick:g}]  ({steps} steps, dt = {dt})")
        d2.Mf.from_numpy(M0.astype(np.float32))
        d2.actf.from_numpy(act.astype(np.int32))
        d2.Bas.from_numpy(d2.SYM_BASIS.astype(np.float32))
        Md0 = (kick * Mth) * act[..., None, None]
        d2.Mdf.from_numpy(Md0.astype(np.float32))
        Mi = [central(M0, ax, h) for ax in range(3)]
        P0 = np.zeros_like(M0)
        for i_ in range(3):
            P0 += np_commf(tw(np_commf(Md0, Mi[i_])), Mi[i_])
        d2.Pf.from_numpy((4.0 * P0).astype(np.float32))
        p2s, umins, hts, hHs = [], [], [], []
        t0 = time.time()
        for n_ in range(steps):
            d2.k_flux(inv2h, BETA)
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
                p2s.append(float(np.einsum("...ab,...ab->...",
                                           M - M0, M0)[actb].mean()))
                umins.append(float(u[actb].min()))
                if not np.isfinite(p2s[-1]):
                    print(f"   NON-FINITE at {n_ + 1}")
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
                hHs.append(float((T + u + BETA * u * u)[actb].sum()) * h ** 3)
            if n_ % 8000 == 7999:
                print(f"   step {n_ + 1:6d} [{time.time() - t0:4.0f}s]"
                      f"  H={hHs[-1]:.4f}  u_min={umins[-1]:+.3f}")
        np.savez(out_npz(kick), kick=kick, p2=np.array(p2s),
                 umin=np.array(umins), ht=np.array(hts), hH=np.array(hHs),
                 dt=dt, dense=DENSE)
        print(f"   saved → {out_npz(kick).name}")
    return analyze_main()


def analyze_main():
    dsamp = DENSE * DT / 2
    print("\n" + "=" * 78)
    print("[N-6c] ω(E) + molten(E) — the excitation axis (all points,"
          " identical pipelines)")
    print("=" * 78)
    pts = []
    for kick in KICKS:                                   # the new weak arms
        if not out_npz(kick).exists():
            continue
        z = np.load(out_npz(kick))
        om, dial, ne = omega1_from(z["p2"], z["umin"], dsamp)
        t_end = ne * 20
        Hwin = z["hH"][z["ht"] <= t_end]
        Hd = float(np.median(Hwin[:max(len(Hwin) // 3, 1)]))
        D2, lam = molten_metrics(z["p2"], z["umin"], dsamp)
        pts.append((f"kick {kick:g}", Hd, om, dial, D2, lam))
    zh = dict(np.load(HERE / "_m5_8_2h_dense.npz"))      # the high-E arms
    for tag in ("R0-kicked", "R2-seed", "R3-jit"):
        om, dial, ne = omega1_from(zh[f"{tag}_p2"], zh[f"{tag}_umin"], dsamp)
        t_end = zh[f"{tag}_t"][ne - 1]
        Hwin = zh[f"{tag}_hH"][zh[f"{tag}_ht"] <= t_end]
        Hd = float(np.median(Hwin[:max(len(Hwin) // 3, 1)]))
        D2, lam = molten_metrics(zh[f"{tag}_p2"], zh[f"{tag}_umin"], dsamp)
        pts.append((tag, Hd, om, dial, D2, lam))
    pts.sort(key=lambda p: p[1])
    print(f"  H_rest (static) = {H_STATIC}; excitation = H_dyn − H_rest")
    print("| arm | H_dyn (early) | excitation/H_rest | ω₁ (dial) | D₂ |"
          " λ_max |")
    print("| --- | --- | --- | --- | --- | --- |")
    for name, Hd, om, dial, D2, lam in pts:
        exc = (Hd - H_STATIC) / H_STATIC
        print(f"| {name} | {Hd:+.2f} | {exc:+.2f} |"
              f" {om:.3f} ({', '.join(f'{v:.2f}' for v in dial)}) |"
              f" {D2:.2f} | {lam:+.3f} |")
    print("\n  read: ω₁ flat in E ⇒ rigidity extends to excitation (with"
          " N-6a: a fully rigid clock);")
    print("  D₂/λ falling toward the two-tone floor (D₂~1.2, λ~+0.06) as"
          " E drops ⇒ the ground clock")
    print("  REGULARIZES (the molten-clock bridge closes); flat/rising ⇒"
          " molten all the way down — honest trouble.")
    print("=" * 78)
    return 0


def settle_main():
    """Deep settle (the 2g S4 protocol, 2× longer) → _m5_8_2o_settled.npz,
    then a dense restart measuring ω/D₂/λ in the EARLY low-E window — the
    route to the GROUND clock (the kick knob cannot lower E: the seed's own
    spontaneous release floors excitation at ~1.7× H_rest — measured)."""
    from openwave.xperiments.m5_liquid_crystal.research.sandbox_v8 import (
        m5_8_2d_quartic_saturation as d2,
    )
    from openwave.xperiments.m5_liquid_crystal.research.sandbox_v8.m5_8_2g_spontaneity import (
        run_spont,
    )
    steps4 = int(sys.argv[2]) if len(sys.argv) > 2 else 48000
    seed = d2.load_seed()
    print(f"[deep settle] S4 damped flow, {steps4} steps (γ = 0.999/step)")
    tr4, M_set = run_spont("S4-deep", seed, d2.DT, steps4, damp=1e-3)
    np.savez(HERE / "_m5_8_2o_settled.npz", M_settled=M_set,
             T_end=tr4["T"][-1], steps=steps4)
    print(f"  T drained {tr4['T'][0]:.3e} → {tr4['T'][-1]:.3e};"
          f" settled config → _m5_8_2o_settled.npz")
    return 0


def restart_main():
    """Dense restart from the DEEP-settled config, P = 0 exact — the low-E
    clock measurement (ω + molten metrics on the early window)."""
    from openwave.xperiments.m5_liquid_crystal.research.sandbox_v8 import (
        m5_8_2d_quartic_saturation as d2,
    )
    steps = int(sys.argv[2]) if len(sys.argv) > 2 else 48000
    dt = DT / 2
    M0s, Mth, act, core, _, h = d2.load_seed()
    zs = np.load(HERE / "_m5_8_2o_settled.npz")
    M_set = zs["M_settled"]
    actb = act > 0.5
    n_act = int(act.sum())
    inv2h = 1.0 / (2.0 * h)
    d2.Mf.from_numpy(M_set.astype(np.float32))
    d2.Mdf.fill(0.0)
    d2.Pf.fill(0.0)
    d2.actf.from_numpy(act.astype(np.int32))
    d2.Bas.from_numpy(d2.SYM_BASIS.astype(np.float32))
    p2s, umins, hts, hHs, hTs = [], [], [], [], []
    t0 = time.time()
    for n_ in range(steps):
        d2.k_flux(inv2h, BETA)
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
            p2s.append(float(np.einsum("...ab,...ab->...",
                                       M - M_set, M_set)[actb].mean()))
            umins.append(float(u[actb].min()))
            if not np.isfinite(p2s[-1]):
                print(f"   NON-FINITE at {n_ + 1}")
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
            hHs.append(float((T + u + BETA * u * u)[actb].sum()) * h ** 3)
        if n_ % 8000 == 7999:
            print(f"   [restart] step {n_ + 1:6d} [{time.time() - t0:4.0f}s]"
                  f"  T={hTs[-1]:.3e}  H={hHs[-1]:.4f}  u_min={umins[-1]:+.3f}")
    np.savez(HERE / "_m5_8_2o_ground.npz", p2=np.array(p2s),
             umin=np.array(umins), ht=np.array(hts), hH=np.array(hHs),
             hT=np.array(hTs), dt=dt, dense=DENSE)
    # the low-E clock readout: thirds of the run (E grows along it — the
    # un-sittable regrowth IS the excitation axis, sampled in time)
    print("\n  [ground-clock readout] thirds of the restart (E grows along"
          " the run):")
    print("  | window | T (median) | ω₁ (dial) | D₂ | λ_max |")
    print("  | --- | --- | --- | --- | --- |")
    n3 = len(p2s) // 3
    for wi, (a, b) in enumerate(((0, n3), (n3, 2 * n3), (2 * n3, 3 * n3))):
        y = np.array(p2s[a:b])
        um = np.array(umins[a:b])
        om, dial, _ = omega1_from(y, um, DENSE * dt)
        D2, lam = molten_metrics(y, um, DENSE * dt)
        Tm = float(np.median(np.array(hTs)[(np.array(hts) >= a * 20)
                                           & (np.array(hts) < b * 20)]))
        print(f"  | third {wi + 1} | {Tm:+.3f} | {om:.3f}"
              f" ({', '.join(f'{v:.2f}' for v in dial)}) | {D2:.2f} |"
              f" {lam:+.3f} |")
    print("  saved → _m5_8_2o_ground.npz")
    return 0


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "analyze"
    if mode == "run":
        return run_main()
    if mode == "settle":
        return settle_main()
    if mode == "restart":
        return restart_main()
    if mode == "analyze":
        return analyze_main()
    print(__doc__.split("USAGE:")[1])
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
