"""
M5.8.2h — N-1: the ω-ATTRACTOR measurement + the rod-localization readout.

THE QUESTION (G-2c-3's closing half): 2g demonstrated SPONTANEITY — a settled
config with P = 0 exact regrows breathing motion, dt-converged. What remains
is the FREQUENCY claim: does the spontaneous state oscillate at the SAME
carrier as the kicked saturated state, independent of how it was started?
Same comb across independent starts ⇒ ω is an ATTRACTOR (a property of the
state, not of the kick) ⇒ G-2c-3 closes.

PRE-REGISTERED SCRUTINY (W1, the `ref` mode — runs on SAVED data, free):
the 2e "ω₀ = 0.262 + exact 2–5× comb" headline is suspect BY CONSTRUCTION:
  · 2e FFT'd the late HALF of a t = 48 run — a 24-τ window. ω₀ = 0.262 ⇒
    period 2π/0.262 ≈ 24 τ = EXACTLY the window: ω₀ sits in FFT BIN 1.
  · every FFT bin is an integer multiple of bin 1 — an "exact 2–5× comb" is
    AUTOMATIC for any non-sinusoidal signal with power at the window scale.
  · meanwhile 2d's saturation evidence ("breathing 30–142, ≥45 periods in
    t ≈ 48") implies a carrier near 2π·45/48 ≈ 5.9 — the 2b-2 linear clock
    class (5.86), NOT 0.262.
W1 re-derives the kicked reference comb with a Hann window + parabolic peak
interpolation, splits fast (ω ≥ 1, well-resolved) from slow (ω < 1, window-
limited — flagged with the bin width Δω) bands, and reports honestly which
peaks are RESOLVED LINES vs window artifacts. This is simultaneously N-3
step 1 (the fast-clock-carrier hunt in the saved dense s(t)).

CHANNEL DESIGN (the rectification trap): dM_rms ≥ 0 is an ENVELOPE — a field
breathing as a(t) = A sin(ωt) about M_ref gives dM ∝ |a(t)| ⇒ FFT power at
2ω, not ω. So the dense probes are:
  dM   — rms departure (envelope; expect the 2ω rectified line)
  s    — clock overlap ⟨(M−M_ref)·Mth⟩_core (SIGNED; the 2e-comparable probe;
         blind to the amplitude channel per 2g — its silence is itself data)
  p2   — ⟨(M−M_ref)·M_ref⟩_act (SIGNED amplitude-channel projection)
  p3   — ⟨(M−M_ref)·R⟩_act, R a FIXED random symmetric field (rng seed 7;
         breaks every symmetry — catches any coherent carrier)
  umin — fuel floor (health marker; deep-floor cascade ⇒ u_min → −O(100))

THE ARMS (W1 OUTCOME FORCED R0: the saved kicked reference's slow "peak"
MOVES with the window length — 0.25/0.32/0.57/0.62 across full/half/quarter
windows — unresolved drift, NOT a line; and its fast band has no stable
carrier in s(t). The kicked state must be re-run with the amplitude-channel
probes for an apples-to-apples comb):
  R0-kicked    the 2e/2g control arm, kick = 0.05 (the comb REFERENCE)
  R1-settled   M_start = the 2g settled config, P = 0 EXACT (the clean arm)
  R2-seed      the raw dressed seed, kick = jitter = 0 (the 2g S1 arm)
  R3-jit       raw seed + 1e-6 symmetric jitter, jseed = 2 (≠ 2g's jseed 1)
All arms: β = 1.558, dt = DT/2 = 0.001, dense sampling every 20 steps
(spacing 0.02 τ, Nyquist ω ≈ 157 ≫ 5.86), default 48000 steps (t = 48).
FFT windows: the HEALTHY window (u_min ≥ −20 and finite H — the cascade cut,
REPORTED not hidden) + thirds-stability check. Verdict = the comb table,
3-outcome (ATTRACTOR / SEED-DEPENDENT / INCONCLUSIVE-window).

ROD-LOCALIZATION READOUT (the M5.8.8 piece, same trajectories, near-free):
the act mask (r > 2RC) & (ρ > RHOC) EXCLUDES the z-axis disclination rod, so
these metrics use the WIDE mask (interior & r > 2RC ONLY — rod included):
  A_z       = ⟨z²⟩_w / (½⟨x²+y²⟩_w)   (1 = spherical, ≫1 = rod-like)
  frac_shell = Σw[r ≤ 4]/Σw            (weight in the core neighborhood)
  frac_rod   = Σw[ρ ≤ 2RHOC, r > 4]/Σw (weight on the far rod segment)
  r_mean     = ⟨r⟩_w
on two weights: w_dep = ‖M−M_ref‖² (where the motion lives) and w_T = |T(x)|
(where the kinetic energy lives), every 4000 steps, against the STATIC
baseline (the seed's |u(x)| — the disclination-rod signature).
LOCALIZED ⇒ frac_shell high, A_z ~ 1 vs baseline; ROD-RIDING ⇒ frac_rod
grows toward the baseline.

PREREQUISITES (in ../sandbox_v8/): _m5_8_2cb_ref.npz (the seed),
_m5_8_2g_settled.npz (the settled config), _m5_8_2e_omega.npz (kicked ref).

RESULTS (2026-06-07 — N-1 COMPLETE; 4 arms × 48000 steps dt/2, ~22 min):
  W1 (ref scrutiny) — THE 2e COMB HEADLINE RETIRED: the saved kicked s(t)
    slow "ω₀" MOVES with the FFT window (0.249/0.323/0.568/0.623 across
    full/half/quarter windows, all ⚠️ ≤2 bins from DC) = unresolved slow
    drift, NOT a line; "ω₀ = 0.262 + exact 2–5× comb" was bin-1-of-window
    + automatic bin arithmetic. The s-only reference was also insufficient
    (no amplitude-channel probes) ⇒ R0-kicked re-run added to the arms.
  W2 (the attractor test) — ATTRACTOR-CONSISTENT, quasi-periodic:
    the saturated state is NOT strictly periodic — it carries a RESOLVED
    fundamental + ~2× harmonic over broadband drift + a strong unresolved
    slow envelope. Line positions are DETREND-STABLE (2.5/5/10 τ filter
    windows move them < 2% — the dial discipline applied to the filter):
      R0-kicked  ω₁ = 1.06–1.10  + 2ω₁ = 2.04–2.10
      R2-seed    ω₁ = 1.15       + 2ω₁ = 2.32–2.34   (unkicked!)
      R3-jit     ω₁ = 1.05–1.09  (+1.33, 1.51)
    Three independent starts (kicked / exactly-unkicked / jittered) share
    the fundamental within one FFT bin (Δω = 0.13) ⇒ ω₁ ≈ 1.1 is a
    property of the STATE, not the start — the attractor signature.
    R1-settled differs (stable 1.79–1.80, +3.90 ≈ 2.2×; thirds chirp
    0.9→1.8) — BUT its spatial readout (W4) shows it never reaches the
    breather within the window (departure weight still flowing INWARD,
    r_mean 6.8→5.2, frac_shell 0.002→0.19 rising; H drifting −157→+562):
    INCONCLUSIVE for that arm (convergence window too short), not a
    refutation. |s| medians: 2.3e-2 (kicked) / 6e-20, 2e-18 (unkicked) —
    the 2g blindness reproduced; rig check: R0 at τ-equiv 12000 gives
    dM 1.478 vs 2g S0's 1.46 (dt-consistent); R1 at τ-equiv 4000 gives
    dM 0.2501 = the 2g dt/2 restart to 4 digits.
  N-3 step 1 (the fast-carrier hunt) — NEGATIVE: no detrend-stable
    5.86-class line in ANY arm/channel (R1's 5.82 appears at one filter
    setting only). The saturated breather's fundamental is ω₁ ≈ 1.1 —
    ~5× softer than the 2b-2 LINEAR clock class (5.86): floor softening
    at saturation. The first ZBW ratio (N-3) must use the measured ω₁.
  W4 (rod localization) — THE BREATHER IS NOT THE ROD: motion weight
    (w_dep, w_T) is near-isotropic, A_z ≈ 1.0–1.3 vs the static seed-|u|
    baseline A_z = 5.58; frac_rod ≈ 0.014–0.038 ≈ the rod VOLUME fraction
    (0.021) — zero excess on the disclination line. Early concentration in
    the core shell (frac_shell 0.69–0.71 at t=4 vs volume fraction 0.225 =
    3× excess) with slow outward spread over t = 48 (→ ~0.27, honest:
    delocalization/heating over long horizons). R1-settled's weight flows
    TOWARD the core all window — the spontaneous motion seeks the defect.
  Honest residuals: (i) ω₁'s absolute value carries one-bin uncertainty
    (±0.13; longer windows blocked by the slow envelope + late heating);
    (ii) R1's attractor convergence unmeasured (needs a longer settled run
    or a pre-converged restart); (iii) the slow band (ω < 1) remains
    window-limited in all arms — no claim made there; (iv) late-window
    heating (H drift, w_T redistribution at t ≳ 36) is the known stepper
    stiffness item, excluded from the line measurements by detrend
    stability, not by cuts.

USAGE:  python m5_8_2h_omega_attractor.py ref            (saved-data scrutiny, free)
        python m5_8_2h_omega_attractor.py run [steps]    (the 3 arms; default 48000)
        python m5_8_2h_omega_attractor.py analyze        (re-analyze saved traces)
"""
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent          # sandbox_vn
V8 = HERE.parent / "sandbox_v8"                 # the v8 reference data
REPO_ROOT = HERE.parents[4]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

# ti-free import chain (geometry + metric — `ref`/`analyze` never init Taichi)
from openwave.xperiments.m5_liquid_crystal.research.sandbox_v8.m5_8_2c1_full_evolution import (  # noqa: E402
    N, L, DT, RC, RHOC, build_grid, central, tw,
)

BETA = 1.558                    # the 2d winner
DENSE = 20                      # dense-probe stride (spacing 0.02 τ at dt/2)
HPROBE = 200                    # H/T probe stride
SPAT = 4000                     # rod-metric snapshot stride
UMIN_CUT = -20.0                # deep-floor cascade marker (healthy: ~−3)
OUT_NPZ = HERE / "_m5_8_2h_dense.npz"


def np_commf(A, B):
    return (np.einsum("...ac,...cb->...ab", A, B)
            - np.einsum("...ac,...cb->...ab", B, A))


# ── geometry: the WIDE mask (rod INCLUDED — the act mask excludes ρ ≤ RHOC) ──

def make_geo():
    g = build_grid()
    xs = np.linspace(-L, L, N)
    X, Y, Z = np.meshgrid(xs, xs, xs, indexing="ij")
    inter = np.zeros(g["r"].shape, bool)
    inter[2:-2, 2:-2, 2:-2] = True
    wide = inter & (g["r"] > 2 * RC)                      # rod kept
    shell = wide & (g["r"] <= 4.0)                        # core neighborhood
    rod = wide & (g["rho"] <= 2 * RHOC) & (g["r"] > 4.0)  # far rod segment
    return dict(h=g["h"], r=g["r"], rho=g["rho"], X=X, Y=Y, Z=Z,
                wide=wide, shell=shell, rod=rod)


def rod_metrics(w, geo):
    """A_z, frac_shell, frac_rod, r_mean of a non-negative weight w(x)."""
    w = np.abs(w)
    tot = w[geo["wide"]].sum()
    if tot <= 0:
        return dict(Az=np.nan, fsh=np.nan, frod=np.nan, rmean=np.nan)
    z2 = (w * geo["Z"] ** 2)[geo["wide"]].sum()
    xy2 = (w * (geo["X"] ** 2 + geo["Y"] ** 2))[geo["wide"]].sum()
    return dict(Az=z2 / (0.5 * xy2),
                fsh=w[geo["shell"]].sum() / tot,
                frod=w[geo["rod"]].sum() / tot,
                rmean=(w * geo["r"])[geo["wide"]].sum() / tot)


# ── spectra: Hann window + parabolic sub-bin peak interpolation ─────────────

def moving_mean(y, n):
    n = min(n, max(len(y) // 4 * 2 - 1, 1))    # clamp to short windows
    k = np.ones(n) / n
    norm = np.convolve(np.ones_like(y), k, "same")
    return np.convolve(y, k, "same") / norm


def spec(y, dsamp, detrend_n=0):
    """(ω, amplitude) of y; detrend_n>0 subtracts a moving mean (envelope kill)."""
    y = np.asarray(y, float)
    y = y - (moving_mean(y, detrend_n) if detrend_n else y.mean())
    w = np.hanning(len(y))
    A = np.abs(np.fft.rfft(y * w))
    om = 2.0 * np.pi * np.fft.rfftfreq(len(y), d=dsamp)
    return om, A


def peaks(om, A, lo, hi, n=5):
    """Top-n local maxima in ω ∈ [lo, hi), parabolic-refined; [(ω, A), ...]."""
    out = []
    for k in range(1, len(A) - 1):
        if not (lo <= om[k] < hi and A[k] > A[k - 1] and A[k] >= A[k + 1]):
            continue
        la, lb, lc = np.log(A[k - 1] + 1e-300), np.log(A[k] + 1e-300), \
            np.log(A[k + 1] + 1e-300)
        den = la - 2 * lb + lc
        d = 0.5 * (la - lc) / den if abs(den) > 1e-12 else 0.0
        out.append((om[k] + d * (om[1] - om[0]), A[k]))
    out.sort(key=lambda p: -p[1])
    return out[:n]


def spec_report(name, y, dsamp, fast_detrend=None):
    """Fast (ω ≥ 1, resolved) + slow (ω < 1, window-limited) peak tables."""
    n = len(y)
    T = n * dsamp
    dbin = 2.0 * np.pi / T
    if fast_detrend is None:
        fast_detrend = max(int(5.0 / dsamp) | 1, 3)      # ~5 τ envelope kill
    om_f, A_f = spec(y, dsamp, detrend_n=fast_detrend)
    om_s, A_s = spec(y, dsamp)
    pf = peaks(om_f, A_f, 1.0, np.inf)
    ps = peaks(om_s, A_s, dbin * 0.5, 1.0)
    print(f"    {name}: window T = {T:.1f} τ ({n} samples), Δω(bin) = {dbin:.4f}")
    print("      fast band (ω ≥ 1, resolved):  "
          + (",  ".join(f"ω={w:.3f}(A={a:.2e})" for w, a in pf) or "—"))
    if pf:
        w0 = pf[0][0]
        print("      ratios to top fast peak:      "
              + ",  ".join(f"{w / w0:.3f}" for w, _ in pf))
    print("      slow band (ω < 1):            "
          + (",  ".join(f"ω={w:.3f}(A={a:.2e})"
                        + ("⚠️bin≤2" if w < 2 * dbin else "")
                        for w, a in ps) or "—"))
    return pf, ps, dbin


# ── W1: scrutiny of the SAVED kicked reference (free; N-3 step 1) ───────────

def ref_main():
    z = np.load(V8 / "_m5_8_2e_omega.npz")
    s = z["s_dense"].astype(float)
    dsamp = DENSE * DT / 2                               # 2e: every 20 @ dt/2
    print("=" * 78)
    print("[W1] scrutiny of the SAVED kicked reference (_m5_8_2e_omega.npz)")
    print(f"     {len(s)} dense s(t) samples, spacing {dsamp} τ, total"
          f" t = {len(s) * dsamp:.0f}")
    print("     2e headline under test: 'ω₀ = 0.262 + exact 2–5× comb' —")
    print("     0.262 = bin 1 of 2e's 24-τ late-half window; bins 2–5 are")
    print("     automatic integer multiples (comb carries no extra info).")
    print("=" * 78)
    half = len(s) // 2
    for name, seg in (("full window  ", s),
                      ("late half (2e's window)", s[half:]),
                      ("late quarter A", s[half:half + (len(s) - half) // 2]),
                      ("late quarter B", s[half + (len(s) - half) // 2:])):
        spec_report(name, seg, dsamp)
    print("\n  [N-3 step 1] the fast-clock-carrier hunt (the 2b-2 linear class"
          " 5.86 / apolar-doubled 11.7): read the fast-band tables above —")
    print("  a resolved line there IS the first ZBW-ratio input (ω_fast).")
    return 0


# ── the dense evolution (run mode only — Taichi inits here) ─────────────────

def run_dense(d2, tag, seed, dt, steps, kick=0.0, jitter=0.0, jseed=1,
              M_start=None, M_ref=None, geo=None, R=None):
    """2g's run_spont + dense signed probes + rod-metric snapshots."""
    M0, Mth, act, core, _, h = seed
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
        Rj = rng.standard_normal(M0.shape)
        Md0 += jitter * 0.5 * (Rj + np.swapaxes(Rj, -1, -2))
    if kick > 0.0 or jitter > 0.0:
        Md0 = Md0 * act[..., None, None]
        d2.Mdf.from_numpy(Md0.astype(np.float32))
        Mi = [central(M0, ax, h) for ax in range(3)]
        P0 = np.zeros_like(M0)
        for i_ in range(3):
            P0 += np_commf(tw(np_commf(Md0, Mi[i_])), Mi[i_])
        d2.Pf.from_numpy((4.0 * P0).astype(np.float32))
    dn = dict(t=[], dM=[], s=[], p2=[], p3=[], umin=[])
    hp = dict(t=[], T=[], H=[])
    sp = []                                              # rod-metric rows
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
            dMv = M - M_ref
            dn["t"].append(n_ + 1)
            dn["dM"].append(float(np.sqrt(np.einsum(
                "...ab,...ab->...", dMv, dMv)[actb].mean())))
            dn["s"].append(float(np.einsum(
                "...ab,...ab->...", dMv, Mth)[core].mean()))
            dn["p2"].append(float(np.einsum(
                "...ab,...ab->...", dMv, M_ref)[actb].mean()))
            dn["p3"].append(float(np.einsum(
                "...ab,...ab->...", dMv, R)[actb].mean()))
            dn["umin"].append(float(u[actb].min()))
            if not np.isfinite(dn["dM"][-1]):
                print(f"   [{tag}] NON-FINITE at step {n_ + 1}")
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
            hp["t"].append(n_ + 1)
            hp["T"].append(float(T[actb].sum()) * h ** 3)
            hp["H"].append(float(
                (T + u + BETA * u * u)[actb].sum()) * h ** 3)
            if n_ % SPAT == SPAT - 1:
                dMv = M - M_ref
                w_dep = np.einsum("...ab,...ab->...", dMv, dMv)
                md, mt = rod_metrics(w_dep, geo), rod_metrics(T, geo)
                sp.append((n_ + 1, md["Az"], md["fsh"], md["frod"],
                           md["rmean"], mt["Az"], mt["fsh"], mt["frod"]))
            if n_ % 4000 == 3999:
                print(f"   [{tag}] step {n_ + 1:6d} [{time.time() - t0:4.0f}s]"
                      f"  dM={dn['dM'][-1]:.3e}  T={hp['T'][-1]:.3e}"
                      f"  u_min={dn['umin'][-1]:+.3f}  H={hp['H'][-1]:.4f}")
    return ({k: np.array(v) for k, v in dn.items()},
            {k: np.array(v) for k, v in hp.items()}, np.array(sp))


def run_main():
    from openwave.xperiments.m5_liquid_crystal.research.sandbox_v8 import (
        m5_8_2d_quartic_saturation as d2,
    )
    steps = int(sys.argv[2]) if len(sys.argv) > 2 else 48000
    dt = DT / 2
    print("=" * 78)
    print("M5.8.2h run — the ω-attractor arms (β = 1.558, dt = DT/2 ="
          f" {dt}, {steps} steps = t {steps * dt:.0f}, dense every {DENSE})")
    print("=" * 78)
    seed = d2.load_seed()
    M0, _, act, _, _, _ = seed
    zs = np.load(V8 / "_m5_8_2g_settled.npz")
    M_set = zs["M_settled"]
    geo = make_geo()
    rng = np.random.default_rng(7)                       # FIXED probe field
    R = rng.standard_normal(M0.shape)
    R = 0.5 * (R + np.swapaxes(R, -1, -2)) * act[..., None, None]
    R /= np.sqrt(np.einsum("...ab,...ab->...", R, R)[act > 0.5].mean())
    # static baselines for the rod metrics
    h = seed[5]
    Mi0 = [central(M0, ax, h) for ax in range(3)]
    u0 = d2.np_u_density(Mi0)
    base_seed = rod_metrics(u0, geo)
    MiS = [central(M_set, ax, h) for ax in range(3)]
    base_set = rod_metrics(d2.np_u_density(MiS), geo)
    print("  static |u| baselines (the disclination-rod signature):")
    print(f"    seed:    A_z={base_seed['Az']:.2f}  frac_shell="
          f"{base_seed['fsh']:.3f}  frac_rod={base_seed['frod']:.3f}"
          f"  r_mean={base_seed['rmean']:.2f}")
    print(f"    settled: A_z={base_set['Az']:.2f}  frac_shell="
          f"{base_set['fsh']:.3f}  frac_rod={base_set['frod']:.3f}"
          f"  r_mean={base_set['rmean']:.2f}")
    arms = (("R0-kicked", dict(kick=0.05)),
            ("R1-settled", dict(M_start=M_set, M_ref=M_set)),
            ("R2-seed", dict()),
            ("R3-jit", dict(jitter=1e-6, jseed=2)))
    out = dict(beta=BETA, dt=dt, steps=steps, dense=DENSE,
               base_seed=np.array([base_seed[k] for k in
                                   ("Az", "fsh", "frod", "rmean")]),
               base_set=np.array([base_set[k] for k in
                                  ("Az", "fsh", "frod", "rmean")]))
    for tag, kw in arms:
        print(f"\n[{tag}]")
        dn, hp, sp = run_dense(d2, tag, seed, dt, steps, geo=geo, R=R, **kw)
        for k, v in dn.items():
            out[f"{tag}_{k}"] = v
        for k, v in hp.items():
            out[f"{tag}_h{k}"] = v
        out[f"{tag}_spat"] = sp
        np.savez(OUT_NPZ, **out)                         # save as we go
        print(f"   saved through {tag} → {OUT_NPZ.name}")
    analyze_main()
    return 0


# ── the analysis (re-runnable from the npz; ti-free) ────────────────────────

def healthy_end(umin):
    """First dense index past the deep-floor cascade cut (else len)."""
    bad = np.where(~np.isfinite(umin) | (umin < UMIN_CUT))[0]
    return int(bad[0]) if len(bad) else len(umin)


def analyze_main():
    z = dict(np.load(OUT_NPZ))
    dsamp = float(z["dense"]) * float(z["dt"])
    print("\n" + "=" * 78)
    print("[W2] the ω-attractor comb tables (healthy windows; judge the data)")
    print("=" * 78)
    combs = {}
    for tag in ("R0-kicked", "R1-settled", "R2-seed", "R3-jit"):
        umin = z[f"{tag}_umin"]
        ne = healthy_end(umin)
        tend = z[f"{tag}_t"][ne - 1] * float(z["dt"])
        print(f"\n  [{tag}] healthy window: {ne}/{len(umin)} samples"
              f" (t ≤ {tend:.1f}; u_min cut {UMIN_CUT});"
              f"  u_min(end) = {umin[ne - 1]:+.3f}")
        print(f"    |s| median = {np.median(np.abs(z[f'{tag}_s'][:ne])):.2e}"
              "  (clock-channel amplitude; machine-zero on unkicked arms ="
              " the 2g blindness)")
        for ch in ("p2", "p3", "dM"):
            y = z[f"{tag}_{ch}"][:ne]
            if len(y) < 64:
                print(f"    {ch}: window too short ({len(y)} samples) — skip")
                continue
            pf, ps, dbin = spec_report(ch, y, dsamp)
            combs[(tag, ch)] = pf
            if len(y) >= 192:                            # thirds stability
                th = len(y) // 3
                tops = []
                for a in range(3):
                    seg = y[a * th:(a + 1) * th]
                    p, _, _ = (peaks(*spec(seg, dsamp,
                               detrend_n=max(int(5.0 / dsamp) | 1, 3)),
                               1.0, np.inf, n=1), None, None)
                    tops.append(p[0][0] if p else np.nan)
                print("      thirds top-fast ω:            "
                      + ",  ".join(f"{v:.3f}" for v in tops))
    print("\n  [W3] rectification sanity: expect ω(dM) ≈ 2 × ω(p2/p3) —")
    print("  the envelope doubles the signed carrier; a match validates the"
          " channel reading.")
    print("\n  cross-arm dominant fast carrier (the attractor test):")
    print("  | arm | ω(p2) | ω(p3) | ω(dM) |")
    print("  | --- | --- | --- | --- |")
    for tag in ("R0-kicked", "R1-settled", "R2-seed", "R3-jit"):
        row = []
        for ch in ("p2", "p3", "dM"):
            pf = combs.get((tag, ch), [])
            row.append(f"{pf[0][0]:.3f}" if pf else "—")
        print(f"  | {tag} | " + " | ".join(row) + " |")
    print("  (kicked reference fast band: see `ref` mode — same pipeline.)")
    print("\n" + "=" * 78)
    print("[W4] rod-localization snapshots (wide mask — rod INCLUDED)")
    print("=" * 78)
    bs, bt = z["base_seed"], z["base_set"]
    print(f"  static |u| baseline (seed):    A_z={bs[0]:.2f} frac_shell="
          f"{bs[1]:.3f} frac_rod={bs[2]:.3f} r_mean={bs[3]:.2f}")
    print(f"  static |u| baseline (settled): A_z={bt[0]:.2f} frac_shell="
          f"{bt[1]:.3f} frac_rod={bt[2]:.3f} r_mean={bt[3]:.2f}")
    for tag in ("R0-kicked", "R1-settled", "R2-seed", "R3-jit"):
        sp = z[f"{tag}_spat"]
        if sp.size == 0:
            continue
        print(f"\n  [{tag}]  (w_dep = ‖M−M_ref‖²;  w_T = |T(x)|)")
        print("  | t | A_z(dep) | f_shell(dep) | f_rod(dep) | r_mean(dep)"
              " | A_z(T) | f_shell(T) | f_rod(T) |")
        print("  | --- | --- | --- | --- | --- | --- | --- | --- |")
        for row in sp:
            print(f"  | {row[0] * float(z['dt']):.0f} | {row[1]:.2f} |"
                  f" {row[2]:.3f} | {row[3]:.3f} | {row[4]:.2f} |"
                  f" {row[5]:.2f} | {row[6]:.3f} | {row[7]:.3f} |")
    print("\n  read: LOCALIZED ⇒ frac_shell ≳ baseline & A_z → ~1;"
          " ROD-RIDING ⇒ frac_rod → baseline.")
    print("=" * 78)
    return 0


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "ref"
    if mode == "ref":
        return ref_main()
    if mode == "run":
        return run_main()
    if mode == "analyze":
        return analyze_main()
    print(__doc__.split("USAGE:")[1])
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
