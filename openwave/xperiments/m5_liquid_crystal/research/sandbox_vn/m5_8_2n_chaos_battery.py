"""
M5.8.2n — N-6d: THE DTC CLASSIFICATION BATTERY — periodic clock / "molten"
chaotic clock / dispersed, on the saved dense traces.

WHY NOW (promoted with N-6): N-1 found the saturated state QUASI-periodic —
a resolved ω₁ + 2ω₁ over broadband drift. The DTC/CTC literature's
classification (the 2D-DTC paper's 3-phase map; the CTC "melting" lesson)
is therefore a LIVE question: does the clock TICK (regular, comb-carrying)
or WANDER (molten/chaotic with a residual spectral ridge)?

THE BATTERY (the M5.8.3 diagnostics-suite spec, on SAVED data — CPU only):
  K    the 0-1 chaos test (Gottwald–Melbourne): K → 0 regular, K → 1
       chaotic; median over 64 random phases c
  D₂   correlation dimension (Grassberger–Procaccia) in a delay embedding
       (τ from the ACF 1/e time, m = 6): low integer ⇒ torus/limit cycle,
       fractional/larger ⇒ chaotic sea
  λ_max  largest Lyapunov exponent (Rosenstein): ≤ 0 regular, > 0 chaotic
       (units: per τ-time, from the dense sampling)

HONEST METHOD NOTES (pre-registered):
  · the traces are NON-STATIONARY (growth + slow envelope) — chaos tests on
    raw nonstationary data give false K → 1. Every statistic runs on the
    DETRENDED series (5 τ moving-mean removed, the N-1 dial's middle
    setting) over the healthy window only; the first-difference variant is
    reported alongside as the robustness check.
  · 2400 samples is SHORT for D₂/λ — the numbers are indicative, anchored
    by CONTROLS run through the IDENTICAL pipeline: a pure two-tone signal
    (ω₁ + 2ω₁, the "ideal clock"), the same + white noise at the measured
    broadband level, and a Lorenz-x chaotic reference. Classification is
    relative to the controls, not absolute.

RESULTS (2026-06-07 — N-6d COMPLETE: the "MOLTEN CLOCK" classification):
  CONTROLS first (the battery calibrates itself):
    two-tone ideal clock:  K ≈ 0 ✓, D₂ ≈ 1.2 (limit cycle ✓ — 2ω₁ is
      commensurate), λ ≈ +0.06 = the short-series BIAS FLOOR for λ;
    two-tone + noise:      K → 1 (the 0-1 test false-flags NOISE as chaos)
      but detrended λ stays ≈ +0.01 — additive noise does NOT inflate λ;
    Lorenz-x:              K ≈ 0 (FALSE-0: oversampled deterministic chaos
      defeats the K-test) but D₂ ≈ 2.2–2.5 ✓ and λ ≈ +0.4–0.7 ✓.
  ⇒ METHOD VERDICT: K is unreliable in BOTH directions on this data class
    (false-1 on noise, false-0 on oversampling) — D₂ + λ carry the
    classification; K reported but not load-bearing.
  THE ARMS (detrended, healthy windows):
    R0-kicked   D₂ 2.85  λ +0.74   |  R2-seed  D₂ 2.97  λ +0.50
    R3-jit      D₂ 2.66  λ +0.54   |  R1-settled (less excited)  D₂ 1.70  λ +0.40
  ⇒ CLASSIFICATION: the saturated breather is a **MOLTEN CLOCK** — a
    persistent ω₁ comb (N-1's attractor) riding on LOW-DIMENSIONAL CHAOTIC
    dressing (λ_max ≈ +0.4–0.7 per τ-time — well above the two-tone bias
    floor AND the noisy-clock control, so the divergence is DYNAMICAL, not
    additive noise; D₂ ≈ 2.7–3.0). Not an ideal periodic crystal; not
    dispersed; the DTC literature's MIDDLE phase.
  ⇒ HYPOTHESIS handed to N-6c: R1 (least excited) is markedly LESS
    dimensional (D₂ 1.70) — molten-ness may grow with excitation; the
    GROUND clock could be regular. The weak-kick ω(E) runs test this.
  Caveats: 2400-sample series (D₂/λ indicative, anchored by controls);
  fixed embedding m = 6; λ slope region not per-arm validated.

USAGE:  python m5_8_2n_chaos_battery.py
"""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[4]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from openwave.xperiments.m5_liquid_crystal.research.sandbox_vn.m5_8_2h_omega_attractor import (  # noqa: E402
    moving_mean, healthy_end,
)

DSAMP = 0.02                     # the 2h dense spacing (20 × dt/2)
RNG = np.random.default_rng(11)


def k01_test(y, n_c=64):
    """Gottwald–Melbourne 0-1 test (regression form). K→0 regular, →1 chaotic."""
    y = np.asarray(y, float)
    N = len(y)
    ncut = N // 10
    Ks = []
    for c in RNG.uniform(0.1 * np.pi, 0.9 * np.pi, n_c):
        ph = np.cumsum(y * np.cos(np.arange(N) * c)), \
             np.cumsum(y * np.sin(np.arange(N) * c))
        M = []
        for n in range(1, ncut):
            dp = ph[0][n:] - ph[0][:-n]
            dq = ph[1][n:] - ph[1][:-n]
            M.append((dp * dp + dq * dq).mean())
        M = np.array(M)
        n_ = np.arange(1, ncut)
        K = np.corrcoef(n_, M)[0, 1]
        Ks.append(K)
    return float(np.median(Ks))


def embed(y, m, tau):
    N = len(y) - (m - 1) * tau
    return np.stack([y[i * tau:i * tau + N] for i in range(m)], axis=1)


def acf_tau(y):
    """Delay = the ACF 1/e crossing (bounded to [2, 50] samples)."""
    y = y - y.mean()
    ac = np.correlate(y, y, "full")[len(y) - 1:]
    ac /= ac[0] + 1e-300
    idx = np.where(ac < 1.0 / np.e)[0]
    return int(np.clip(idx[0] if len(idx) else 10, 2, 50))


def corr_dim(y, m=6):
    """Grassberger–Procaccia D₂ (indicative; short-series caveat applies)."""
    tau = acf_tau(y)
    X = embed((y - y.mean()) / (y.std() + 1e-300), m, tau)
    n = len(X)
    if n > 1200:                                   # cap pair count
        X = X[RNG.choice(n, 1200, replace=False)]
        n = 1200
    d = np.sqrt(((X[:, None, :] - X[None, :, :]) ** 2).sum(-1))
    iu = np.triu_indices(n, k=tau)                 # Theiler window ~ tau
    dv = d[iu]
    rs = np.percentile(dv, [5, 10, 20, 30, 40, 50])
    rs = rs[rs > 0]
    C = np.array([(dv < r).mean() for r in rs])
    ok = C > 0
    if ok.sum() < 3:
        return np.nan, tau
    sl = np.polyfit(np.log(rs[ok]), np.log(C[ok]), 1)[0]
    return float(sl), tau


def lyap_rosenstein(y, m=6, horizon=60):
    """Largest Lyapunov (Rosenstein): slope of mean log nearest-neighbor
    divergence over the linear region; per τ-TIME (uses DSAMP)."""
    tau = acf_tau(y)
    X = embed((y - y.mean()) / (y.std() + 1e-300), m, tau)
    n = len(X)
    d = np.sqrt(((X[:, None, :] - X[None, :, :]) ** 2).sum(-1))
    np.fill_diagonal(d, np.inf)
    for off in range(1, min(tau, n)):              # Theiler exclusion
        idx = np.arange(n - off)
        d[idx, idx + off] = np.inf
        d[idx + off, idx] = np.inf
    nn = np.argmin(d, axis=1)
    div = []
    for k in range(1, horizon):
        v = []
        for i in range(n - k):
            j = nn[i]
            if j + k < n:
                dd = np.linalg.norm(X[i + k] - X[j + k])
                if dd > 0:
                    v.append(np.log(dd))
        div.append(np.mean(v) if v else np.nan)
    div = np.array(div)
    ks = np.arange(1, horizon)
    ok = np.isfinite(div[: horizon // 2])
    if ok.sum() < 5:
        return np.nan
    sl = np.polyfit(ks[: horizon // 2][ok], div[: horizon // 2][ok], 1)[0]
    return float(sl / DSAMP)


def battery(name, y):
    det = y - moving_mean(y, max(int(5.0 / DSAMP) | 1, 3))
    diff = np.diff(y)
    rows = []
    for variant, s in (("detrended", det), ("1st-diff", diff)):
        K = k01_test(s)
        D2, tau = corr_dim(s)
        lam = lyap_rosenstein(s)
        rows.append((variant, K, D2, tau, lam))
    print(f"\n  [{name}]  ({len(y)} samples)")
    print("  | variant | K (0-1 test) | D₂ | τ_delay | λ_max (/τ-time) |")
    print("  | --- | --- | --- | --- | --- |")
    for variant, K, D2, tau, lam in rows:
        print(f"  | {variant} | {K:+.3f} | {D2:.2f} | {tau} | {lam:+.3f} |")
    return rows


def main():
    print("=" * 78)
    print("M5.8.2n — N-6d: the DTC classification battery (controls first —"
          " judge arms vs controls)")
    print("=" * 78)
    n = 2400
    t = np.arange(n) * DSAMP
    # controls through the IDENTICAL pipeline
    two_tone = np.sin(1.1 * t) + 0.4 * np.sin(2.2 * t + 0.7)
    noisy = two_tone + 0.25 * RNG.standard_normal(n)
    x = np.empty(n)                                  # Lorenz-x reference
    s_, v_, w_ = 1.0, 1.0, 1.0
    dt_l = 0.01
    burn = 2000
    for i in range(burn + n):
        ds = 10.0 * (v_ - s_)
        dv = s_ * (28.0 - w_) - v_
        dw = s_ * v_ - 8.0 / 3.0 * w_
        s_, v_, w_ = s_ + dt_l * ds, v_ + dt_l * dv, w_ + dt_l * dw
        if i >= burn:
            x[i - burn] = s_
    battery("CONTROL: two-tone (ideal clock ω₁+2ω₁)", two_tone)
    battery("CONTROL: two-tone + noise (measured broadband level)", noisy)
    battery("CONTROL: Lorenz-x (chaotic reference)", x)
    z = dict(np.load(HERE / "_m5_8_2h_dense.npz"))
    for tag in ("R0-kicked", "R1-settled", "R2-seed", "R3-jit"):
        ne = healthy_end(z[f"{tag}_umin"])
        battery(f"{tag} (p2, healthy {ne})", z[f"{tag}_p2"][:ne])
    print("\n  read: an arm matching the NOISY-clock control row (K"
          " intermediate, D₂ low-ish, λ≈0)")
    print("  ⇒ a TICKING clock with broadband dressing; matching Lorenz"
          " (K→1, λ>0, D₂≈2+) ⇒ MOLTEN.")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
