"""M6.1 — Q=0 neutral chaoiton mass (sandbox v2).

Werbos's 2026-05-17 #1 priority ask. Computes the ground-state mass m_χ of the
lightest neutral chaoiton in the Ouroboros system, using the 4-function radial
ODE from his 2026-05-18 benchmark document.

Reads:
    ../theory/Numerical Benchmark for the Ouroboros Lagrangian.docx (§3, §9)
    ../theory/DarkMatterv1.pdf (§5 — locked-ansatz semi-analytic estimate)
    ../0a_background.md § 12.3 (the ask)
    ../0a_background.md § 12.4 (open structural questions)

Writes:
    ../plots/m6_1_neutral_chaoiton_profile.png
    ../plots/m6_1_neutral_chaoiton_energy.png
    ../plots/m6_1_neutral_chaoiton_summary.json

Approach:
    1. Initialize at r=r_min using locked ansatz: J = -A, Q = -V
       (so Q_CS = 0 by construction at the seed)
    2. Integrate the FULL 4-function ODE outward (don't enforce locked ansatz
       as identity — let it evolve and observe drift)
    3. Sweep initial amplitude to find lowest-energy localized solution
    4. Compute m_χ = H_code, then m_χ_MeV = H_code × ℏc / R^phys
    5. Convergence check by doubling N_r and r_max

The four open structural questions from §12.4 of background are exposed as
CLI flags so we can sweep them once the default run is in hand.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

# ---------------------------------------------------------------------------
# Benchmark parameters (from §6 of Werbos's Numerical Benchmark doc)
# ---------------------------------------------------------------------------

G_COUPLING = 1.0625        # electron calibration coupling
LAMBDA = 1.0               # quartic coupling
OMEGA = 1.0                # electron-calibration oscillation frequency (code units)
A0_DEFAULT = 0.1           # A-field amplitude seed
B0_DEFAULT = 0.1           # J-field amplitude seed (equals A0 by locked ansatz)

# Numerical protocol (benchmark §10)
R_MIN = 0.02               # regularity offset to avoid 1/r singularity
R_MAX_DEFAULT = 30.0       # benchmark recommends 20-50
N_GRID_DEFAULT = 2000      # benchmark recommends 2000-5000
RTOL = 1.0e-10
ATOL = 1.0e-12

# Physical scale (from Electron Calibration §5.1)
R_PHYS_FM = 191.0
HBAR_C_MEVFM = 197.32698   # CODATA
M_E_MEV = 0.51100          # electron rest mass

# Localization threshold (benchmark §2.1 for spectrum, §4.2 for calibration)
LOCALIZATION_THRESHOLD = 0.15   # |V| + |A| + |Q| + |J| at r_max

# Convergence target (benchmark §10)
CONVERGENCE_TOL = 0.001    # H must change <0.1% on doubling N_r and r_max


# ---------------------------------------------------------------------------
# ODE — the full 4-function system from benchmark §3
# ---------------------------------------------------------------------------

def ode_full_4fn(r, y, g, lam, m_J_sq, omega, laplacian_dim=2,
                 omega_insertion="none"):
    """Full 4-function Ouroboros radial ODE per benchmark §3.

    Variables:
        y = [V, dV/dr, A, dA/dr, Q, dQ/dr, J, dJ/dr]

    Benchmark equations:
        Δ_r V = Q
        Δ_r A = J
        Δ_r Q = V + m_J² Q + λ Q (Q² − J²)
        Δ_r J = A − m_J² J − λ J (Q² − J²)

    where Δ_r f = f''(r) + (D-1)/r · f'(r), D=2 (benchmark) or D=3 (sandbox v1).

    Open question §12.4 — where ω enters. Three modes:
        "none"     — equations as written in benchmark (no ω term)
        "explicit" — add −ω² × field on LHS of each equation
        "in_m_J"   — m_J² is replaced by m_J² + ω²
    """
    V, dV, A, dA, Q, dQ, J, dJ = y

    inv_r_factor = (laplacian_dim - 1) / r

    s = Q * Q - J * J  # (Q² − J²) used in the cubic terms

    # Base RHS per benchmark
    rhs_V = Q
    rhs_A = J
    rhs_Q = V + m_J_sq * Q + lam * Q * s
    rhs_J = A - m_J_sq * J - lam * J * s

    # ω insertion choice (open question §12.4)
    if omega_insertion == "explicit":
        # Add −ω² × field, so Δ_r V + ω² V = Q etc.
        rhs_V = rhs_V - omega * omega * V
        rhs_A = rhs_A - omega * omega * A
        rhs_Q = rhs_Q - omega * omega * Q
        rhs_J = rhs_J - omega * omega * J
    elif omega_insertion == "in_m_J":
        # Replace m_J² → m_J² + ω² in mass terms only
        delta = omega * omega
        rhs_Q = V + (m_J_sq + delta) * Q + lam * Q * s
        rhs_J = A - (m_J_sq + delta) * J - lam * J * s

    # Recover f''(r) from f''(r) + ((D-1)/r) f'(r) = rhs  →  f'' = rhs − ((D-1)/r) f'
    d2V = rhs_V - inv_r_factor * dV
    d2A = rhs_A - inv_r_factor * dA
    d2Q = rhs_Q - inv_r_factor * dQ
    d2J = rhs_J - inv_r_factor * dJ

    return [dV, d2V, dA, d2A, dQ, d2Q, dJ, d2J]


# ---------------------------------------------------------------------------
# m_J² choice (open question §12.4)
# ---------------------------------------------------------------------------

def m_J_squared(formula, g, omega):
    """Map open-question formula name → m_J² value.

    The benchmark says 'm_J is determined by g' without giving the formula.
    Candidates to sweep:
        0          → no mass term (matches Calibration paper's f(s)=gs² with g→λ)
        g          → m_J² = g (linear)
        g_sq       → m_J² = g² (quadratic)
        omega2     → m_J² = ω² (oscillation-tied)
        omega_g    → m_J² = ω·g (product)
    """
    formulas = {
        "0": 0.0,
        "g": g,
        "g_sq": g * g,
        "omega2": omega * omega,
        "omega_g": omega * g,
    }
    if formula not in formulas:
        raise ValueError(f"Unknown m_J² formula: {formula!r}. "
                         f"Choices: {list(formulas)}")
    return formulas[formula]


# ---------------------------------------------------------------------------
# Energy functional (benchmark §5)
# ---------------------------------------------------------------------------

def energy_density(V, dV, A, dA, Q, dQ, J, dJ, m_J_sq, lam):
    """Hamiltonian density from benchmark §5.

    H = (2π)² R ∫ r dr [
            (V')² + (A')² + (Q')² + (J')²
            − V·Q + A·J
            + (1/2) m_J² (Q² − J²)
            + (λ/4) (Q² − J²)²
        ]

    Returns the integrand WITHOUT the (2π)² R prefactor; caller multiplies r dr
    and integrates. The (2π)² R prefactor is dimensional — for code-unit mass
    we use a single chaoiton "ring radius" R = 1 (will be absorbed in R^phys
    calibration later if needed).
    """
    s = Q * Q - J * J
    kinetic = dV * dV + dA * dA + dQ * dQ + dJ * dJ
    cross = -V * Q + A * J
    mass = 0.5 * m_J_sq * s
    quartic = 0.25 * lam * s * s
    return kinetic + cross + mass + quartic


def integrate_energy(r, V, dV, A, dA, Q, dQ, J, dJ, m_J_sq, lam, ring_R=1.0):
    """Integrate H = (2π)² R ∫ r dr × density."""
    density = energy_density(V, dV, A, dA, Q, dQ, J, dJ, m_J_sq, lam)
    integrand = r * density
    # (2π)² R prefactor
    prefactor = (2.0 * np.pi) ** 2 * ring_R
    return prefactor * np.trapezoid(integrand, r)


# ---------------------------------------------------------------------------
# Localization + Q_CS = 0 checks
# ---------------------------------------------------------------------------

def is_localized(V, A, Q, J, threshold=LOCALIZATION_THRESHOLD):
    """All four field magnitudes at r_max sum below threshold."""
    tail = abs(V[-1]) + abs(A[-1]) + abs(Q[-1]) + abs(J[-1])
    return tail < threshold, tail


def chern_simons_charge(r, A, dJ, J, dA, ring_R=1.0):
    """Q_CS = (2πR) ∫ r dr [A·J' − J·A']  (benchmark §5, toroidal form).

    For the locked ansatz J=−A: A·J' = A·(−A') = −AA', and J·A' = −A·A',
    so A·J' − J·A' = −AA' − (−AA') = 0 — exactly zero. Any departure from 0
    measures how much the locked ansatz is breaking during integration.
    """
    integrand = r * (A * dJ - J * dA)
    return 2.0 * np.pi * ring_R * np.trapezoid(integrand, r)


# ---------------------------------------------------------------------------
# Single integration run
# ---------------------------------------------------------------------------

def run_neutral_chaoiton(g=G_COUPLING, lam=LAMBDA, omega=OMEGA,
                          m_J_sq=0.0,
                          A0=A0_DEFAULT, B0=B0_DEFAULT,
                          r_min=R_MIN, r_max=R_MAX_DEFAULT, n_grid=N_GRID_DEFAULT,
                          laplacian_dim=2, omega_insertion="none",
                          locked_sign=-1, method="RK45"):
    """Integrate full 4-fn ODE with locked-ansatz seed; return diagnostics.

    locked_sign: -1 → Werbos's J=-A, Q=-V (top eq gives Bessel J_0 oscillation)
                 +1 → J=+A, Q=+V       (top eq gives K_0 modified-Bessel decay)
    """

    # Locked ansatz at r=r_min: Q = s·V, J = s·A with s=locked_sign
    s = locked_sign
    y0 = [
        B0, 0.0,        # V, V'
        A0, 0.0,        # A, A'
        s * B0, 0.0,    # Q = s·V
        s * A0, 0.0,    # J = s·A
    ]

    r_eval = np.linspace(r_min, r_max, n_grid)

    sol = solve_ivp(
        fun=lambda r, y: ode_full_4fn(
            r, y, g, lam, m_J_sq, omega,
            laplacian_dim=laplacian_dim,
            omega_insertion=omega_insertion,
        ),
        t_span=(r_min, r_max),
        y0=y0,
        t_eval=r_eval,
        method=method,
        rtol=RTOL,
        atol=ATOL,
        dense_output=False,
    )

    if not sol.success:
        return {"success": False, "message": sol.message}

    r = sol.t
    V, dV, A, dA, Q, dQ, J, dJ = sol.y

    localized, tail = is_localized(V, A, Q, J)
    H_code = integrate_energy(r, V, dV, A, dA, Q, dQ, J, dJ, m_J_sq, lam)
    Q_CS = chern_simons_charge(r, A, dJ, J, dA)

    # locked-ansatz drift: max|Q − s·V| and max|J − s·A| over r
    locked_drift_QV = float(np.max(np.abs(Q - locked_sign * V)))
    locked_drift_JA = float(np.max(np.abs(J - locked_sign * A)))

    # Convert m_χ_code → MeV via electron calibration scale
    # m_e = H_code_electron × ℏc / R^phys → R^phys / ℏc = H_code_electron / m_e
    # m_χ_MeV = H_code_χ × ℏc / R^phys = H_code_χ × m_e / H_code_electron
    # Calibration paper §5.1: H_code for electron is 0.494. We use that.
    H_CODE_ELECTRON = 0.494
    m_chi_mev = H_code * M_E_MEV / H_CODE_ELECTRON

    return {
        "success": True,
        "r": r, "V": V, "A": A, "Q": Q, "J": J,
        "dV": dV, "dA": dA, "dQ": dQ, "dJ": dJ,
        "H_code": float(H_code),
        "m_chi_MeV": float(m_chi_mev),
        "Q_CS": float(Q_CS),
        "tail": float(tail),
        "localized": bool(localized),
        "locked_drift_QV": locked_drift_QV,
        "locked_drift_JA": locked_drift_JA,
        "params": {
            "g": g, "lam": lam, "omega": omega, "m_J_sq": m_J_sq,
            "A0": A0, "B0": B0,
            "r_min": r_min, "r_max": r_max, "n_grid": n_grid,
            "laplacian_dim": laplacian_dim,
            "omega_insertion": omega_insertion,
        },
    }


# ---------------------------------------------------------------------------
# Convergence check
# ---------------------------------------------------------------------------

def convergence_check(base_result, **run_kwargs):
    """Verify H_code stable to <0.1% on doubling N_r and r_max separately."""
    base_H = base_result["H_code"]
    base_params = base_result["params"]

    kw = {k: v for k, v in base_params.items()}
    # Double N_r
    kw_double_N = {**kw, "n_grid": kw["n_grid"] * 2}
    res_N = run_neutral_chaoiton(**kw_double_N)
    delta_N = abs(res_N["H_code"] - base_H) / abs(base_H)

    # Double r_max
    kw_double_R = {**kw, "r_max": kw["r_max"] * 2}
    res_R = run_neutral_chaoiton(**kw_double_R)
    delta_R = abs(res_R["H_code"] - base_H) / abs(base_H)

    return {
        "base_H": base_H,
        "double_N_H": res_N["H_code"], "delta_N_pct": 100 * delta_N,
        "double_R_H": res_R["H_code"], "delta_R_pct": 100 * delta_R,
        "converged": delta_N < CONVERGENCE_TOL and delta_R < CONVERGENCE_TOL,
    }


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def save_plots(result, out_dir):
    """Field profiles + energy density plots."""
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("[skip] matplotlib not available — no plots")
        return

    r = result["r"]
    fig, ax = plt.subplots(2, 2, figsize=(10, 7))
    ax[0, 0].plot(r, result["V"], label="V(r)")
    ax[0, 0].plot(r, result["A"], label="A(r)")
    ax[0, 0].set_xlabel("r"); ax[0, 0].set_ylabel("V, A")
    ax[0, 0].legend(); ax[0, 0].set_title("Poloidal (V, A)")

    ax[0, 1].plot(r, result["Q"], label="Q(r)")
    ax[0, 1].plot(r, result["J"], label="J(r)")
    ax[0, 1].set_xlabel("r"); ax[0, 1].set_ylabel("Q, J")
    ax[0, 1].legend(); ax[0, 1].set_title("Toroidal (Q, J)")

    ax[1, 0].plot(r, result["Q"] + result["V"], label="Q + V")
    ax[1, 0].plot(r, result["J"] + result["A"], label="J + A")
    ax[1, 0].axhline(0, color="k", lw=0.5)
    ax[1, 0].set_xlabel("r"); ax[1, 0].set_ylabel("locked-ansatz drift")
    ax[1, 0].legend(); ax[1, 0].set_title("Locked-ansatz drift (should ≈ 0)")

    # Energy density
    s = result["Q"] ** 2 - result["J"] ** 2
    m_J_sq = result["params"]["m_J_sq"]
    lam = result["params"]["lam"]
    density = (result["dV"] ** 2 + result["dA"] ** 2
               + result["dQ"] ** 2 + result["dJ"] ** 2
               - result["V"] * result["Q"] + result["A"] * result["J"]
               + 0.5 * m_J_sq * s
               + 0.25 * lam * s * s)
    ax[1, 1].plot(r, r * density)
    ax[1, 1].set_xlabel("r"); ax[1, 1].set_ylabel("r × energy density")
    ax[1, 1].set_title(f"Energy density (H_code = {result['H_code']:.4f})")

    plt.tight_layout()
    out_path = out_dir / "m6_1_neutral_chaoiton_profile.png"
    plt.savefig(out_path, dpi=120)
    plt.close()
    print(f"[plot] {out_path}")


def save_summary(result, convergence, out_dir):
    """Serializable JSON summary."""
    summary = {
        "H_code": result["H_code"],
        "m_chi_MeV": result["m_chi_MeV"],
        "Q_CS": result["Q_CS"],
        "tail": result["tail"],
        "localized": result["localized"],
        "locked_drift_QV": result["locked_drift_QV"],
        "locked_drift_JA": result["locked_drift_JA"],
        "params": result["params"],
        "convergence": convergence,
    }
    out_path = out_dir / "m6_1_neutral_chaoiton_summary.json"
    with open(out_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"[json] {out_path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--m-j-sq", default="0",
                        help="m_J² formula: 0 | g | g_sq | omega2 | omega_g")
    parser.add_argument("--omega-insertion", default="none",
                        choices=["none", "explicit", "in_m_J"])
    parser.add_argument("--laplacian", type=int, default=2, choices=[2, 3])
    parser.add_argument("--A0", type=float, default=A0_DEFAULT)
    parser.add_argument("--r-max", type=float, default=R_MAX_DEFAULT)
    parser.add_argument("--n-grid", type=int, default=N_GRID_DEFAULT)
    parser.add_argument("--skip-convergence", action="store_true")
    args = parser.parse_args()

    m_J_sq_value = m_J_squared(args.m_j_sq, G_COUPLING, OMEGA)

    print(f"\n{'='*72}")
    print(f"M6.1 — Q=0 neutral chaoiton (sandbox v2)")
    print(f"{'='*72}")
    print(f"g            = {G_COUPLING}")
    print(f"λ            = {LAMBDA}")
    print(f"ω            = {OMEGA}")
    print(f"m_J² formula = {args.m_j_sq!r}  →  m_J² = {m_J_sq_value}")
    print(f"ω insertion  = {args.omega_insertion!r}")
    print(f"Laplacian D  = {args.laplacian}")
    print(f"A₀ = B₀      = {args.A0}")
    print(f"r ∈ [{R_MIN}, {args.r_max}],  N_grid = {args.n_grid}")
    print(f"Locked ansatz initial:  J = −A,  Q = −V  →  Q_CS = 0 at seed")
    print()

    result = run_neutral_chaoiton(
        g=G_COUPLING, lam=LAMBDA, omega=OMEGA,
        m_J_sq=m_J_sq_value,
        A0=args.A0, B0=args.A0,    # locked ansatz: B0 = A0
        r_min=R_MIN, r_max=args.r_max, n_grid=args.n_grid,
        laplacian_dim=args.laplacian,
        omega_insertion=args.omega_insertion,
    )

    if not result["success"]:
        print(f"[FAIL] integration failed: {result['message']}")
        return 1

    print(f"H_code              = {result['H_code']:.6f}")
    print(f"m_χ (MeV)           = {result['m_chi_MeV']:.6f}")
    print(f"Q_CS                = {result['Q_CS']:.6e}  (should be 0)")
    print(f"Tail (|V|+|A|+|Q|+|J| @ r_max)")
    print(f"                    = {result['tail']:.6e}  "
          f"(threshold {LOCALIZATION_THRESHOLD})")
    print(f"Localized?          = {result['localized']}")
    print(f"Locked drift |Q+V|  = {result['locked_drift_QV']:.6e}")
    print(f"Locked drift |J+A|  = {result['locked_drift_JA']:.6e}")
    print()

    # Convergence check
    if not args.skip_convergence:
        print("Convergence check (doubling N_r and r_max)…")
        convergence = convergence_check(result)
        print(f"  base H              = {convergence['base_H']:.6f}")
        print(f"  double N_r:  H      = {convergence['double_N_H']:.6f}  "
              f"(Δ = {convergence['delta_N_pct']:.3f}%)")
        print(f"  double r_max: H     = {convergence['double_R_H']:.6f}  "
              f"(Δ = {convergence['delta_R_pct']:.3f}%)")
        print(f"  converged           = {convergence['converged']}")
    else:
        convergence = {"skipped": True}

    # Decision branch
    print()
    print(f"{'='*72}")
    m_chi = result["m_chi_MeV"]
    if not result["localized"]:
        verdict = ("NO-LOCALIZED — locked ansatz may be incompatible with "
                   "current ODE form. Try different m_J² formula or ω insertion.")
    elif 0.5 <= m_chi <= 10.0:
        verdict = "COLD-DM CANDIDATE — m_χ in MeV range; compute σv → Ω_χ h²"
    elif 0.001 <= m_chi <= 0.1:
        verdict = "DARK-RADIATION CANDIDATE — m_χ sub-MeV; different cosmology"
    elif m_chi < 0.001:
        verdict = "OUT-OF-RANGE LOW — check parameter sweep, ODE form"
    else:
        verdict = "OUT-OF-RANGE HIGH — check parameter sweep, ODE form"
    print(f"VERDICT: {verdict}")
    print(f"{'='*72}\n")

    # Outputs
    out_dir = Path(__file__).resolve().parent.parent / "plots"
    out_dir.mkdir(exist_ok=True)
    save_plots(result, out_dir)
    save_summary(result, convergence, out_dir)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
