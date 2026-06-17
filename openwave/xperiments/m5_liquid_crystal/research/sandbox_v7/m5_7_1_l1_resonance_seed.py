"""
M5.7.1 — resonance-hunt pipeline: l=1 director-rotation seed on the biaxial hedgehog

Close's protocol (2026-04): seed an l=1 harmonic perturbation on the matrix-field defect,
sweep its amplitude, and look for a regime where the energy stays localized substantially
longer than it disperses — "an unstable particle or resonance." This script builds the
PIPELINE and validates it at two amplitudes (linear reference vs A/λ=1 test); the full
{0.5, 1, 2} sweep is M5.7.2.

WHAT WE INHERIT FROM M5.6 (do not rebuild — this is the numpy mirror of production):
  • biaxial hedgehog  M = O·D(s(r))·Oᵀ , D=diag(1,δ,0)  (seed_biaxial_hedgehog_M, frame
    {r̂, e_Θ, e_Φ}, radial Faber melt s(r)=r/√(r²+r₀²), z-disclination smoothstep melt).
  • Eq.18 leapfrog   M_tt = c²·div(G) − dV_M ,  G_α = 8Σ_ν[[M_α,M_ν],M_ν]  (= evolve_M +
    compute_curvature_flux). Simple ½‖Ṁ‖² kinetic — OK for LIFETIME (M5.6.5d: its only null
    mode is the trace, never sourced by the traceless curvature force; the faithful O(x)∈SO(3)
    kinetic is needed only for the clock FREQUENCY at M5.8).
  • V-on b=0 amplitude well  V = a·Tr(M²) + c·(Tr(M²))² , min at s₂*=1+δ²  (M5.6.5c). A
    SIMILARITY perturbation M→R M Rᵀ preserves Tr(M²), so V is exactly flat to the seed —
    the perturbation lives purely in the kinetic + curvature (twist) sector. Correct: V
    confines amplitude, the l=1 twist is the QM/δ sector V leaves alone.

THE PERTURBATION (design locked 2026-05-28).  An l=1 director rotation
about ŷ:  M_pert = R_y(α)·M_bg·R_y(α)ᵀ ,  α(x) = δθ_peak · g(x) ,  g = f(r)·(z/r) (the
Y_10 dipole pattern), f a shell localized to the defect's textured region. Rotating the
director IS what an EM-wave drive does (M5.6.4 Maxwell-from-tilts), so the M5.7 seed = the
resonant drive — one code path. A/λ=1 ⇒ δθ_peak=2π (a full turn at the shell peak);
the per-defect length scale (λ=4r₀) matches the per-defect picture (the excitation lives in the
defect's own (A,ω) excess, not box-scale modes).

KEY PHYSICS FINDING (informs the baseline).  The curvature force is CUBIC in ∂M, so around
constant vacuum (∂D=0) it vanishes to all linear orders — there is NO linear wave propagation
in vacuum (only the local V_M mass term). The Skyrme kinetic activates only where the
background already has gradient (the hedgehog texture, C_μν≠0, M5.6.2). So a "Gaussian in
vacuum" is NOT a free disperser here. The correct dispersion baseline is the LINEAR-AMPLITUDE
limit of the SAME seed on the SAME background — exactly Close's amplitude knob.

METRIC (control-subtracted intensity localization — sign-safe, apples-to-apples).
  δM(x,t) = M_perturbed(x,t) − M_control(x,t)   (control = unperturbed biaxial evolving
            alongside; isolates the seed from the hedgehog's own intrinsic slosh, M5.6.2b)
  I(x,t)  = ‖δM‖²_F  (≥0; the matrix analog of |ψ|² — the textbook wave-packet tracer)
  L(t)    = I_local(r<R_core) / I_total          localization fraction, R_core = 3·r₀
  τ       = first t with L(t) ≤ 0.5·L(0)         the localization lifetime
We report τ_res/τ_lin: ≫1 = the nonlinear amplitude self-traps = a resonance signal.
(Secondary, for continuity with the roadmap's energy intent: control-subtracted excess
energy fraction, clipped to the regions the seed ADDED energy.)

USAGE:
    python -m openwave.xperiments.m5_liquid_crystal.research.sandbox_v7.m5_7_1_l1_resonance_seed
"""
import os
import sys

import numpy as np

# ----- substrate (matches _topo_biaxial1_von.py: δ=0.30, K=1.0, b=0 well) ----------
DTYPE = np.float32                       # sandbox: f32 is 2× faster, ample for a lifetime metric
DELTA = 0.30
S2_STAR = 1.0 + DELTA**2                 # Tr(diag(1,δ,0)²) — the V-on amplitude target
N = int(os.environ.get("M57_N", 44))     # grid points/axis (env-overridable for resolution sweeps)
L = 5.0
RC = 0.9                                 # radial core melt = r₀ (Faber scale)
RHOC = 0.9                               # z-axis disclination melt
C_WAVE = 0.3                             # matches the GUI c (5a bounded-not-bug run)
C_COEF = 2.0                             # V-on well stiffness (K=1 analog)
A_COEF = -2.0 * C_COEF * S2_STAR         # a = −2c·s₂* keeps the b=0 well min at s₂*
DT = float(os.environ.get("M57_DT", 0.004))  # < 5c's 0.006 (env-overridable; finer grid → smaller for CFL)
N_STEPS = int(os.environ.get("M57_STEPS", 2000))
SAMPLE_EVERY = 20
R_CORE = 3.0 * RC                        # measurement aperture (= 2.7)

# perturbation shell — localized in the ACTIVE textured region (r>2·RC; the inner core is
# isotropic/regularized and frozen, so the l=1 twist must live in the shell where the frame
# is well-defined). Y_10 dipole pattern (z/r) × Gaussian shell.
SHELL_PEAK = 2.5 * RC                    # 2.25 — just outside the frozen core, inside R_core
SHELL_W = 0.7


# ----- matrix helpers (numpy mirror of the @ti.func ops; identical math to 5c) ------
def matmul(A, B):
    return np.matmul(A, B)               # BLAS batched gemm over (...,3,3) — far faster than einsum


def commf(A, B):
    return np.matmul(A, B) - np.matmul(B, A)


def frob2(A):
    return np.einsum("...ab,...ab->...", A, A)


def tr(A):
    return np.einsum("...aa->...", A)


def central(f, axis, h):
    out = np.zeros_like(f)
    sp, sm, so = [slice(None)] * f.ndim, [slice(None)] * f.ndim, [slice(None)] * f.ndim
    sp[axis], sm[axis], so[axis] = slice(2, None), slice(0, -2), slice(1, -1)
    out[tuple(so)] = (f[tuple(sp)] - f[tuple(sm)]) / (2 * h)
    return out


def V_M(M, a, b, c):
    M2 = matmul(M, M)
    tr2 = tr(M2)
    tr3 = tr(matmul(M2, M))
    return a * tr2 - b * tr3 + c * tr2 * tr2


def dV_M(M, a, b, c):
    M2 = matmul(M, M)
    tr2 = tr(M2)[..., None, None]
    return 2.0 * a * M - 3.0 * b * M2 + 4.0 * c * tr2 * M


def force(M, h, a, b, c):
    Mx, My, Mz = central(M, 0, h), central(M, 1, h), central(M, 2, h)
    cc = lambda A, B: commf(commf(A, B), B)
    Gx = 8.0 * (cc(Mx, My) + cc(Mx, Mz))
    Gy = 8.0 * (cc(My, Mx) + cc(My, Mz))
    Gz = 8.0 * (cc(Mz, Mx) + cc(Mz, My))
    divG = central(Gx, 0, h) + central(Gy, 1, h) + central(Gz, 2, h)
    return C_WAVE**2 * divG - dV_M(M, a, b, c)


def energy_density(M, M_prev, h, a, b, c):
    Mdot = (M - M_prev) / DT
    kinetic = 0.5 * frob2(Mdot)
    Mx, My, Mz = central(M, 0, h), central(M, 1, h), central(M, 2, h)
    curv = 4.0 * (frob2(commf(Mx, My)) + frob2(commf(Mx, Mz)) + frob2(commf(My, Mz)))
    return kinetic + C_WAVE**2 * curv + V_M(M, a, b, c)


# ----- biaxial hedgehog background (production-equivalent; mirror of seed_biaxial_hedgehog_M)
def build_biaxial_M():
    xs = np.linspace(-L, L, N)
    h = xs[1] - xs[0]
    X, Y, Z = np.meshgrid(xs, xs, xs, indexing="ij")
    r = np.sqrt(X**2 + Y**2 + Z**2) + 1e-9
    rho = np.sqrt(X**2 + Y**2) + 1e-9
    rhat = np.stack([X / r, Y / r, Z / r], -1)
    azim = np.stack([-Y / rho, X / rho, np.zeros_like(r)], -1)
    sdisc = np.clip(rho / RHOC, 0, 1)
    shrink = (sdisc**2 * (3 - 2 * sdisc))[..., None]
    ephi = azim * shrink
    ephi = ephi - np.einsum("...a,...a->...", ephi, rhat)[..., None] * rhat
    etheta = np.cross(ephi, rhat)
    srad = (r / np.sqrt(r**2 + RC**2))[..., None, None]
    diso = (1.0 + DELTA) / 3.0
    d0 = diso + srad[..., 0, 0] * (1.0 - diso)
    d1 = diso + srad[..., 0, 0] * (DELTA - diso)
    d2 = diso + srad[..., 0, 0] * (0.0 - diso)
    oprod = lambda v: v[..., :, None] * v[..., None, :]
    M = (d0[..., None, None] * oprod(rhat) + d1[..., None, None] * oprod(etheta)
         + d2[..., None, None] * oprod(ephi))
    interior = np.zeros(r.shape, bool)
    interior[2:-2, 2:-2, 2:-2] = True
    active = (r > 2 * RC) & (rho > RHOC) & interior      # off point-core + disclination (5c mask)
    return dict(M=M.astype(DTYPE), h=DTYPE(h), r=r.astype(DTYPE),
                X=X.astype(DTYPE), Y=Y.astype(DTYPE), Z=Z.astype(DTYPE), active=active)


# ----- l=1 director-rotation perturbation -------------------------------------------
def perturbation_field(bg):
    """Normalized seed geometry g(x) = f(r)·(z/r) (Y_10 dipole × localized shell), peak 1."""
    r, Z = bg["r"], bg["Z"]
    f = np.exp(-((r - SHELL_PEAK) ** 2) / (2 * SHELL_W**2))    # localized active-region shell
    g = f * (Z / r)                                            # l=1 dipole (cos θ)
    g = g / (np.abs(g).max() + 1e-30)                          # peak |g| = 1
    return g


def apply_l1_rotation(M_bg, g, dtheta_peak):
    """M_pert = R_y(α)·M_bg·R_y(α)ᵀ , α = dtheta_peak·g — local l=1 twist about ŷ.
    A similarity transform: preserves Tr(M²) (V-flat) — pure orientation/twist seed."""
    a = (dtheta_peak * g).astype(DTYPE)
    ca, sa = np.cos(a), np.sin(a)
    z = np.zeros_like(a)
    o = np.ones_like(a)
    R = np.stack([np.stack([ca, z, sa], -1),
                  np.stack([z, o, z], -1),
                  np.stack([-sa, z, ca], -1)], -2)            # (...,3,3) R_y(α) per voxel
    return matmul(matmul(R, M_bg), np.swapaxes(R, -1, -2)).astype(DTYPE)


# ----- lockstep leapfrog: control + N seeds stepped together, scalar metrics inline --
def run_lockstep(bg, seeds, n_steps):
    """Evolve [control, *seeds] in lockstep (Dirichlet BC, released from rest), computing
    the control-subtracted intensity metric inline — no full-field snapshots stored.

    seeds: list of (label, a_over_lambda, M_pert0). Returns per-seed scalar series.
    Metric (vs the co-evolved control): I=‖M_seed−M_ctrl‖²; L=I_local(r<R_core)/I_total;
    R_rms of I. The control is the unperturbed biaxial (its own intrinsic slosh, M5.6.2b).

    FLOOR-NORMALIZED lifetime: a fully-dispersed seed spreads ~uniformly over the active
    region, so L → L_floor = (active voxels in core)/(active voxels). The localization
    EXCESS  Lnorm(t) = (L(t)−L_floor)/(L(0)−L_floor)  starts at 1, → 0 on full dispersion —
    this normalizes out both the dispersion floor AND any L(0) difference between seeds, so
    amplitudes are comparable even when the (nonlinear) seed patterns differ. τ = Lnorm↓50%."""
    h, active, r = bg["h"], bg["active"], bg["r"]
    mask = active[..., None, None]
    in_core = active & (r < R_CORE)
    r2_active = (r**2 * active).astype(np.float64)
    L_floor = in_core.sum() / active.sum()                   # uniform-dispersion floor

    inits = [bg["M"]] + [s[2] for s in seeds]                # index 0 = control
    cur = [m.copy() for m in inits]
    prev = [m.copy() for m in inits]
    series = [dict(L=[], Rrms=[], Itot=[], dH=[]) for _ in seeds]
    Hctrl = []                                               # control total H (own intrinsic slosh)
    steps = []

    def total_H(fi):
        return float((energy_density(cur[fi], prev[fi], h, A_COEF, 0.0, C_COEF) * active).sum())

    def sample():
        Mc = cur[0]
        Hc = total_H(0)
        Hctrl.append(Hc)
        for si in range(len(seeds)):
            I = (frob2(cur[si + 1] - Mc) * active).astype(np.float64)
            tot = I.sum() + 1e-30
            series[si]["L"].append(I[in_core].sum() / tot)
            series[si]["Rrms"].append(np.sqrt((I * r2_active).sum() / tot))
            series[si]["Itot"].append(tot)
            series[si]["dH"].append(total_H(si + 1) - Hc)    # EXCESS H from the seed (should be conserved)

    sample(); steps.append(0)
    for step in range(1, n_steps + 1):
        for fi in range(len(cur)):
            f = force(cur[fi], h, A_COEF, 0.0, C_COEF)
            new = np.where(mask, 2 * cur[fi] - prev[fi] + DT**2 * f, inits[fi])
            prev[fi], cur[fi] = cur[fi], new
        if step % SAMPLE_EVERY == 0:
            sample(); steps.append(step)
            if step % 400 == 0:
                Ls = " ".join(f"{seeds[si][0]}:L={series[si]['L'][-1]:.2f}" for si in range(len(seeds)))
                print(f"      step {step:>5}/{n_steps}  {Ls}", flush=True)
    finite = all(np.isfinite(c).all() for c in cur)
    steps = np.array(steps)
    Hctrl = np.array(Hctrl)
    ctrl_drift = float((Hctrl.max() - Hctrl.min()) / (abs(Hctrl[0]) + 1e-30))
    out = []
    for si, (label, a_ol, _) in enumerate(seeds):
        L = np.array(series[si]["L"]); Rrms = np.array(series[si]["Rrms"]); Itot = np.array(series[si]["Itot"])
        dH = np.array(series[si]["dH"])                      # excess seed energy over time
        Lnorm = (L - L_floor) / (L[0] - L_floor + 1e-30)     # localization excess, starts at 1
        tau, censored = lifetime(Lnorm, steps)
        # seed-energy conservation: |ΔH(t) − ΔH(0)| / ΔH(0) — flat ⇒ physical redistribution,
        # growing ⇒ numerical pumping. (ΔH(0) is the energy the seed injected.)
        dH0 = dH[0] if abs(dH[0]) > 1e-12 else 1e-12
        dH_drift = float(np.max(np.abs(dH - dH[0])) / abs(dH0))
        out.append(dict(label=label, a=a_ol, L=L, Lnorm=Lnorm, Rrms=Rrms, Itot=Itot, steps=steps,
                        tau_t=tau * DT, censored=censored, dH0=float(dH[0]), dH_drift=dH_drift))
    return out, finite, L_floor, ctrl_drift


def lifetime(Lnorm, steps):
    """τ = first sampled step where the localization excess Lnorm ≤ 0.5; (tau_steps, censored?)."""
    below = np.where(Lnorm <= 0.5)[0]
    if len(below) == 0:
        return steps[-1], True                               # never decayed within the window
    return steps[below[0]], False


# ====================================================================================
def main():
    print("=" * 78)
    print("M5.7.1 — l=1 resonance-hunt pipeline on the biaxial hedgehog")
    print(f"  grid {N}³  L={L}  c={C_WAVE}  dt={DT}  steps={N_STEPS}  δ={DELTA}  V-on(b=0,s₂*={S2_STAR:.2f})  f32")
    print("=" * 78)
    bg = build_biaxial_M()
    print(f"  active voxels (off core r>2r₀ + disclination): {100*bg['active'].mean():.0f}%   "
          f"R_core(measure)={R_CORE:.2f}  shell peak r={SHELL_PEAK:.2f}", flush=True)

    g = perturbation_field(bg)
    # δθ_peak = π·(A/λ) (recalibrated 2026-05-28): A/λ=1 ⇒ π (director antipodal at peak =
    # genuine max displacement, no 2π wrap-to-self). Linear ref + the {0.5,1,2} sweep.
    amps = [0.05, 0.5, 1.0, 2.0]
    seeds = [(f"A/λ={a}", a, apply_l1_rotation(bg["M"], g, np.pi * a)) for a in amps]
    print(f"\n  lockstep evolve: control + {len(seeds)} seeds  δθ_peak = "
          f"{', '.join(f'{np.pi*a:.2f}' for a in amps)} rad …", flush=True)
    res_list, finite, L_floor, ctrl_drift = run_lockstep(bg, seeds, N_STEPS)

    print(f"\n  RESULTS  (dispersion floor L_floor={L_floor:.3f}; τ = localization-excess Lnorm↓50%):")
    print(f"  control total-H drift over the run = {100*ctrl_drift:.2f}%  "
          f"({'conservative' if ctrl_drift < 0.05 else 'DRIFTING — check CFL/dt'})")
    for rr in res_list:
        cen = "  [CENSORED → lower bound]" if rr["censored"] else ""
        print(f"    [{rr['label']:<8}] τ={rr['tau_t']:6.2f}t   "
              f"Lnorm:{rr['Lnorm'][0]:.2f}→{rr['Lnorm'][-1]:.2f}  R_rms:{rr['Rrms'][0]:.2f}→{rr['Rrms'][-1]:.2f}  "
              f"seed-ΔH drift {100*rr['dH_drift']:+.0f}%{cen}")

    lin = res_list[0]                                        # A/λ=0.05 linear dispersion baseline
    sweep = res_list[1:]
    print("\n" + "=" * 78)
    pipeline_ok = finite and (len(lin["L"]) > 2) and (lin["Lnorm"][0] > 0)
    print(f"M5.7.1 PIPELINE:  finite-runs={finite}  metric-computed={pipeline_ok}  (N={N}, {N_STEPS} steps, dt={DT})")
    print(f"  baseline τ_lin(A/λ=0.05) = {lin['tau_t']:.2f}t   final Lnorm={lin['Lnorm'][-1]:.2f}")
    # rank by τ, but break censored ties by FINAL Lnorm (higher = more localized)
    best = max(sweep, key=lambda d: (d["tau_t"], d["Lnorm"][-1]))
    for rr in sweep:
        ratio = rr["tau_t"] / lin["tau_t"] if lin["tau_t"] > 0 else float("inf")
        flag = "  ← most localized" if rr is best else ""
        cons = "" if rr["dH_drift"] < 0.10 else "  ⚠️ΔH-drift (not energy-conserving — suspect numerical)"
        print(f"    {rr['label']}:  τ/τ_lin={ratio:.1f}  final-Lnorm={rr['Lnorm'][-1]:.2f}"
              f"{'  (censored)' if rr['censored'] else ''}{flag}{cons}")
    best_ratio = best["tau_t"] / lin["tau_t"] if lin["tau_t"] > 0 else float("inf")
    energy_ok = best["dH_drift"] < 0.10 and ctrl_drift < 0.05
    nonmono = best["Lnorm"][-1] > lin["Lnorm"][-1] + 0.15    # nonlinear retains clearly more than linear
    if nonmono and (best["censored"] or best_ratio >= 2.0):
        verdict = ("RESONANCE CANDIDATE (energy-conserving)" if energy_ok
                   else "localization signal but ⚠️ENERGY DRIFT — likely a numerical artifact, not physical")
        print(f"  → {verdict}: {best['label']} retains Lnorm={best['Lnorm'][-1]:.2f} "
              f"vs linear {lin['Lnorm'][-1]:.2f}; τ ≥ {best_ratio:.0f}× baseline.")
    else:
        print("  → No clear self-trapping; amplitudes disperse comparably (Close's 'mostly dispersion').")
    print("PASS" if pipeline_ok else "FAIL — pipeline did not produce a clean metric")
    print("=" * 78)
    return 0 if pipeline_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
