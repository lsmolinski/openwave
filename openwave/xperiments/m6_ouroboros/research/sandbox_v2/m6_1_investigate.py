"""M6.1 — Investigation sweep for the Q=0 neutral chaoiton ODE.

Goal: characterize the locked-ansatz behavior empirically across the open
structural questions from §12.4 of background. Focus on what CAN be computed
without integrator stalls (locked-sign +/-, asymmetric seed, ω-insertion).

Key prior finding (smoke test):
  - locked-ansatz with m_J²=0, sign=-1 (Werbos's convention) → Bessel J_0
    oscillation, energy DIVERGES linearly with r_max (104% on doubling)
  - locked-ansatz with m_J²≠0, sign=-1 → integrator fails (consistency
    constraint V²-A² = -m_J²/λ violated by seed V(0)=A(0))

This sweep explores:
  1. Sign convention (+1 gives K_0 modified-Bessel decay structure)
  2. Asymmetric seed with m_J²≠0 (resolves the consistency violation)
  3. ω-insertion variants
  4. 3D Laplacian (sandbox-v1 convention) comparison

Each run is wrapped in a timeout-style guard via solve_ivp's max_step.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from m6_1_neutral_chaoiton import (
    G_COUPLING, LAMBDA, OMEGA, R_MIN,
    run_neutral_chaoiton, LOCALIZATION_THRESHOLD,
)


def analyze_decay(r, V, A, Q, J):
    """Tail characterization: oscillating (Bessel) vs decaying (K_0)."""
    n = len(r)
    tail_slice = slice(int(0.8 * n), n)
    r_tail = r[tail_slice]
    V_tail = V[tail_slice]

    sign_changes = int(np.sum(np.diff(np.sign(V_tail)) != 0))

    abs_V = np.abs(V_tail) + 1e-300
    try:
        slope, _ = np.polyfit(r_tail, np.log(abs_V), 1)
    except (ValueError, np.linalg.LinAlgError):
        slope = float("nan")

    return sign_changes, float(slope)


def safe_run(label, **kwargs):
    """Wrap run_neutral_chaoiton with overflow + exception protection."""
    try:
        with np.errstate(over="ignore", invalid="ignore"):
            res = run_neutral_chaoiton(**kwargs)
    except Exception as e:
        return {"label": label, "ok": False, "err": str(e)[:50]}

    if not res["success"]:
        return {"label": label, "ok": False, "err": res.get("message", "?")[:50]}

    # check for NaN / inf in fields
    if not (np.all(np.isfinite(res["V"])) and np.all(np.isfinite(res["A"]))):
        return {"label": label, "ok": False, "err": "NaN/inf in fields"}

    sc, slope = analyze_decay(res["r"], res["V"], res["A"], res["Q"], res["J"])
    return {
        "label": label, "ok": True,
        "H": res["H_code"], "mMeV": res["m_chi_MeV"],
        "tail": res["tail"], "loc": res["localized"],
        "Q_CS": res["Q_CS"],
        "dQV": res["locked_drift_QV"], "dJA": res["locked_drift_JA"],
        "sc": sc, "slope": slope,
        "V_at_rmax": res["V"][-1], "A_at_rmax": res["A"][-1],
    }


def print_table(rows):
    print()
    print(f"{'Label':<54} {'OK':>3} {'loc':>3} {'H':>9} {'mMeV':>9} "
          f"{'tail':>8} {'sc':>3} {'slope':>7} {'V(rmax)':>10}")
    print("-" * 116)
    for r in rows:
        if not r["ok"]:
            print(f"{r['label']:<54} {'NO':>3} (FAIL: {r['err']})")
            continue
        loc = "Y" if r["loc"] else "n"
        slope_str = f"{r['slope']:+.3f}" if math.isfinite(r["slope"]) else "  nan"
        print(f"{r['label']:<54} {'Y':>3} {loc:>3} {r['H']:9.2f} {r['mMeV']:9.3f} "
              f"{r['tail']:8.4f} {r['sc']:>3} {slope_str:>7} {r['V_at_rmax']:+10.4f}")


def main():
    rows = []

    print(f"\n{'='*100}")
    print("M6.1 — Q=0 chaoiton investigation sweep")
    print(f"{'='*100}")

    # ── Sweep 1: locked_sign at m_J²=0, multiple r_max ──────────────────
    print(">>> SWEEP 1: locked_sign (Werbos's J=-A vs J=+A) at m_J²=0")
    for sign in (-1, +1):
        for r_max in (30.0, 60.0, 120.0):
            label = f"S1 sign={sign:+d}, m_J²=0, r_max={r_max:>5.0f}"
            rows.append(safe_run(label,
                g=G_COUPLING, lam=LAMBDA, omega=OMEGA, m_J_sq=0.0,
                A0=0.1, B0=0.1, r_min=R_MIN, r_max=r_max, n_grid=2000,
                laplacian_dim=2, omega_insertion="none",
                locked_sign=sign, method="RK45"))

    # ── Sweep 2: asymmetric seed satisfying V²-A² = -m_J²/λ at r_min ─────
    print(">>> SWEEP 2: asymmetric seed (V(0)² - A(0)² = -m_J²/λ consistency)")
    for sign in (-1, +1):
        for m_val in (0.01, 0.1, 0.5, 1.0):
            # A0 chosen so V0² - A0² = -m_val/λ → A0² = V0² + m_val/λ
            V0 = 0.1
            A0_sq = V0 * V0 + m_val / LAMBDA
            A0 = math.sqrt(A0_sq)
            label = f"S2 sign={sign:+d}, m_J²={m_val:.2f}, A0={A0:.3f}"
            rows.append(safe_run(label,
                g=G_COUPLING, lam=LAMBDA, omega=OMEGA, m_J_sq=m_val,
                A0=A0, B0=V0, r_min=R_MIN, r_max=30.0, n_grid=2000,
                laplacian_dim=2, omega_insertion="none",
                locked_sign=sign, method="RK45"))

    # ── Sweep 3: ω-insertion variants at m_J²=0 ────────────────────────
    print(">>> SWEEP 3: ω-insertion variants at m_J²=0")
    for sign in (-1, +1):
        for omega_ins in ("explicit", "in_m_J"):
            label = f"S3 sign={sign:+d}, m_J²=0, ω_ins={omega_ins}"
            rows.append(safe_run(label,
                g=G_COUPLING, lam=LAMBDA, omega=OMEGA, m_J_sq=0.0,
                A0=0.1, B0=0.1, r_min=R_MIN, r_max=30.0, n_grid=2000,
                laplacian_dim=2, omega_insertion=omega_ins,
                locked_sign=sign, method="RK45"))

    # ── Sweep 4: 3D Laplacian (sandbox-v1 convention) ──────────────────
    print(">>> SWEEP 4: 3D Laplacian at m_J²=0")
    for sign in (-1, +1):
        label = f"S4 sign={sign:+d}, m_J²=0, dim=3"
        rows.append(safe_run(label,
            g=G_COUPLING, lam=LAMBDA, omega=OMEGA, m_J_sq=0.0,
            A0=0.1, B0=0.1, r_min=R_MIN, r_max=30.0, n_grid=2000,
            laplacian_dim=3, omega_insertion="none",
            locked_sign=sign, method="RK45"))

    print_table(rows)

    out = Path(__file__).resolve().parent.parent / "plots" / "m6_1_investigation_sweep.json"
    with open(out, "w") as f:
        json.dump(rows, f, indent=2, default=str)
    print(f"\n[json] {out}")

    # ── Verdict ────────────────────────────────────────────────────────
    localized = [r for r in rows if r["ok"] and r["loc"] and abs(r["slope"]) > 0.05]
    print(f"\n{'='*100}")
    if not localized:
        print("VERDICT: NO localized + decaying solution found in any variant.")
        print("        Locked ansatz appears incompatible with the benchmark ODE form.")
    else:
        print(f"VERDICT: {len(localized)} localized + decaying variants found:")
        for r in sorted(localized, key=lambda x: abs(x["H"]))[:5]:
            print(f"   {r['label']}  →  H={r['H']:.2f}, m_χ={r['mMeV']:.3f} MeV, slope={r['slope']:+.3f}")
    print(f"{'='*100}\n")


if __name__ == "__main__":
    main()
