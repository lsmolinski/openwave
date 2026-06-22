"""
M5.9.4 - lepton rest-energy ledger on the PRODUCTION 4x4 engine, WITH the negative
energy contributions (boost-GEM + clock activation). The full production engine replacement
for the 3x3 numpy toy (m5_9_1/2/3) that Dr. Duda correctly flagged as too simplified.

WHY THIS SCRIPT EXISTS (Duda's review, 2026-06-20):
  - "3x3 diag(1,delta,0) ... they definitely require full 4x4 field, allowing negative
     mass/energy contributions from gravity and oscillations."  -> we now run the
     production 4x4 Taichi field with the boost (gravity) axis live.
  CONVENTION (Duda flagged the eta placement, 2026-06-20/21; engine FLIPPED to index-0 on
     2026-06-21, see research/_convention_refactor/): the ENGINE now orders the eigenvalues
     D = diag(g, 1, delta, 0) with g at INDEX 0 (the time/boost axis, Duda's convention), and
     the Minkowski metric (signed_dot4) puts its minus at index 0 too -- so g and the metric
     minus are on the SAME axis (physics correct). The flip was proven physics-neutral by the
     golden master (66 invariants, worst rel 4e-7). The neutrino build #236 uses the same index-0.
  - "Faber's ansatz is only an approximation ... serious calculations use it only as the
     starting point of energy minimization."  -> we MINIMIZE (gradient flow), not evaluate
     a fixed profile.
  - "the Python files are much too simple ... I have no idea what parameters are used."
     -> every parameter is set from the production TOPOLOGY_SEED block and dumped to JSON.

THE LEDGER (the load-bearing new physics):
  The engine's DISPLAY energy compute_energyH_density_M is EUCLIDEAN (Frobenius norm,
  always >= 0) so it cannot show a negative gravity term. The engine's DYNAMICS are
  Minkowski-signed (signed_dot4, eta = diag(-1,1,1,1); the minus sits at index 0 = the g/time
  axis, see CONVENTION above): the (alpha,0)/(0,alpha) "time-axis" components enter with a MINUS sign. We add a SIGNED energy kernel mirroring the display
  one but using signed_dot4 for the curvature. Then per configuration:

    H_euclid  = kinetic + c^2*curv_euclid + (V - v0)        (all positive; the old view)
    H_signed  = kinetic + c^2*curv_signed + (V - v0)        (the physical, conserved energy)
    boost_GEM = H_euclid - H_signed = 2*c^2*sum_(alpha,0) [M_mu,M_nu]^2   (>= 0, the dip)

  Rest energy = H_signed of the MINIMIZED, clock-active configuration. The boost lowers it
  (gravity), and activating the clock lowers it further (oscillation). Both are the negative
  contributions the 3x3 could not see.

Self-contained, headless (CPU Taichi), matplotlib Agg. Run from anywhere:
    python m5_9_4_engine_lepton.py
"""
import os
import sys
import json
import time

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../"))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

import numpy as np
import taichi as ti

ti.init(arch=ti.cpu, default_fp=ti.f32)

from openwave.xperiments.m5_liquid_crystal import medium
from openwave.xperiments.m5_liquid_crystal import engine1_seeds as seeds
from openwave.xperiments.m5_liquid_crystal import engine2_pde as pde
from openwave.xperiments.m5_liquid_crystal import engine3_observables as obs

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))

# ----------------------------------------------------------------------------
# PARAMETERS - mirror the production TOPOLOGY_SEED (_topo_biaxial1_von.py),
# documented to Duda's standard. Every knob is here and dumped to JSON.
# ----------------------------------------------------------------------------
P = dict(
    EDGE=1e-15,            # universe edge (m); internal "am" scale, absolute scale via e_scale
    TARGET_VOXELS=48 ** 3,  # grid (production runs 64^3; 48^3 for the foundation pass)
    R0_FRACTION=0.06,      # radial melt scale r0_vox = R0_FRACTION * N (Faber core)
    RHOC_VOXELS=3.0,       # z-disclination melt scale (voxels)
    BIAXIAL_DELTA=0.30,    # middle eigenvalue delta; spatial D = diag(1, delta, 0)
    B_STAR=0.13,           # boost dressing strength (GEM dip ~0.13 in production)
    RW_FRACTION=0.29,      # boost dressing radial width rw_vox = RW_FRACTION * N
    CLOCK_KICK=0.05,       # clock phase kick (radians) opening the leapfrog clock
    LDG_STIFFNESS_K=1.0,   # LdG potential coefficient (units c^2/dx^4)
    LC_G=8.0,              # the boost/time axis eigenvalue g (medium.LC_G)
    DT_SCALE_4D=0.5,       # dt derate for the 4D leapfrog path
    C_AMRS=0.2998,         # wave speed am/rs (WAVE_SPEED*RONTOSECOND/ATTOMETER, SIM_SPEED=1)
    KM_INERTIA_4D=30.0,    # diagonal faithful-lite inertia (leapfrog path)
    FLOW_STEPS=400,        # gradient-flow (minimization) steps
    FLOW_DAMP=0.6,         # velocity-retention per step in the constrained flow (1=undamped)
    CLOCK_STEPS=600,       # clock-active evolution steps for the time-average
    M_E_MEV=0.510999,      # electron rest energy (calibration anchor)
)

# optional env overrides (reproducibility: document the exact run, vary resolution/steps)
P["TARGET_VOXELS"] = int(os.environ.get("LEPTON_TARGET_VOXELS", P["TARGET_VOXELS"]))
P["FLOW_STEPS"] = int(os.environ.get("LEPTON_FLOW_STEPS", P["FLOW_STEPS"]))
P["CLOCK_STEPS"] = int(os.environ.get("LEPTON_CLOCK_STEPS", P["CLOCK_STEPS"]))
_SUFFIX = os.environ.get("LEPTON_OUT_SUFFIX", "")   # e.g. "_64" -> m5_9_4_results_64.json

# physical lepton targets (ratios are the real test)
M_E, M_MU, M_TAU = 0.510999, 105.658, 1776.86
R_MU, R_TAU = M_MU / M_E, M_TAU / M_E   # 206.77, 3477.2


# ----------------------------------------------------------------------------
# Build the production field + scratch energy fields (placed before any kernel)
# ----------------------------------------------------------------------------
EDGE = P["EDGE"]
tf = medium.TensorField([EDGE, EDGE, EDGE], P["TARGET_VOXELS"], [0.5, 0.5, 0.5], viz_stride=2)
N = tf.nx
C = N // 2
GRID = tf.grid_size

# scratch per-voxel energy densities (interior-only written; boundary stays 0)
e_kin = ti.field(ti.f32, shape=GRID)   # 1/2 ||Mdot||^2  (Euclidean kinetic)
e_kins = ti.field(ti.f32, shape=GRID)  # 1/2 signed_dot4(Mdot,Mdot)  (Minkowski kinetic)
e_ce = ti.field(ti.f32, shape=GRID)    # c^2 * curvature_euclid (Frobenius)
e_cs = ti.field(ti.f32, shape=GRID)    # c^2 * curvature_signed (Minkowski)
e_pot = ti.field(ti.f32, shape=GRID)   # V_M - v0
e_boost = ti.field(ti.f32, shape=GRID)  # per-voxel (alpha,0) curvature weight (the GEM density)

trackers = medium.Trackers(GRID)
observables = medium.FieldObservables(GRID)

# derived numbers
dx_am = float(tf.dx_am)
c_amrs = P["C_AMRS"]
dt_rs = P["DT_SCALE_4D"] * dx_am * 0.95 / (c_amrs * 3 ** 0.5)
delta = P["BIAXIAL_DELTA"]
ldg_c = P["LDG_STIFFNESS_K"] * c_amrs ** 2 / dx_am ** 4
ldg_a = -2.0 * ldg_c * (1.0 + delta ** 2)
ldg_b = 0.0
tr2_vac = 1.0 + delta ** 2
v0 = ldg_a * tr2_vac + ldg_c * tr2_vac ** 2
r0_vox = P["R0_FRACTION"] * N
rhoc_vox = P["RHOC_VOXELS"]
rw_vox = P["RW_FRACTION"] * N

print("=" * 84)
print("M5.9.4 - lepton rest-energy ledger on the production 4x4 engine")
print("=" * 84)
print(f"grid {N}x{N}x{N} (={N**3} voxels)  dx_am={dx_am:.4e}  dt_rs={dt_rs:.4e}")
print(f"LdG: a={ldg_a:.3e} b={ldg_b:.1f} c={ldg_c:.3e}  v0={v0:.3e}")
print(f"r0_vox={r0_vox:.2f}  rhoc_vox={rhoc_vox:.2f}  rw_vox={rw_vox:.2f}  delta={delta}  g={P['LC_G']}")


# ----------------------------------------------------------------------------
# SIGNED energy kernel - the GEM instrumentation (the new physics)
# ----------------------------------------------------------------------------
@ti.kernel
def energy_decompose(tfx: ti.template(),  # type: ignore
                     kin: ti.template(), kins: ti.template(),  # type: ignore
                     ce: ti.template(), cs: ti.template(),  # type: ignore
                     pot: ti.template(), boost: ti.template(),  # type: ignore
                     c_a: ti.f32, dt: ti.f32, a: ti.f32, b: ti.f32,  # type: ignore
                     cc: ti.f32, vv0: ti.f32):  # type: ignore
    nx, ny, nz = tfx.nx, tfx.ny, tfx.nz
    inv_dt = 1.0 / dt
    c2 = c_a * c_a
    inv_2dx = 1.0 / (2.0 * tfx.dx_am)
    for i, j, k in ti.ndrange((1, nx - 1), (1, ny - 1), (1, nz - 1)):
        m = tfx.M_am[i, j, k]
        m_dot = (m - tfx.M_prev_am[i, j, k]) * inv_dt
        kin[i, j, k] = 0.5 * m_dot.norm_sqr()                  # Euclidean kinetic (>= 0)
        kins[i, j, k] = 0.5 * pde.signed_dot4(m_dot, m_dot)    # Minkowski kinetic (time axis < 0)
        mx = (tfx.M_am[i + 1, j, k] - tfx.M_am[i - 1, j, k]) * inv_2dx
        my = (tfx.M_am[i, j + 1, k] - tfx.M_am[i, j - 1, k]) * inv_2dx
        mz = (tfx.M_am[i, j, k + 1] - tfx.M_am[i, j, k - 1]) * inv_2dx
        cxy = pde.commutator(mx, my)
        cxz = pde.commutator(mx, mz)
        cyz = pde.commutator(my, mz)
        curv_e = cxy.norm_sqr() + cxz.norm_sqr() + cyz.norm_sqr()
        curv_s = (pde.signed_dot4(cxy, cxy) + pde.signed_dot4(cxz, cxz)
                  + pde.signed_dot4(cyz, cyz))
        ce[i, j, k] = c2 * 4.0 * curv_e
        cs[i, j, k] = c2 * 4.0 * curv_s
        pot[i, j, k] = pde.V_M(m, a, b, cc) - vv0
        # GEM density = (alpha,0) part of the curvature: curv_euclid - curv_signed = the dip
        boost[i, j, k] = c2 * 4.0 * (curv_e - curv_s)  # = 2*c^2*4*sum (alpha,0) comps^2


@ti.kernel
def scale_field(fld: ti.template(), s: ti.f32):  # type: ignore
    for I in ti.grouped(fld):
        fld[I] = fld[I] * s


def ledger(tag):
    """Compute the full signed/euclid energy ledger of the CURRENT field state.

    H_euclid = kin_euclid + c^2*curv_euclid + (V-v0)   (the old all-positive view)
    H_signed = kin_signed + c^2*curv_signed + (V-v0)   (the physical Minkowski energy)
      boost_GEM      = curv_euclid - curv_signed  (>=0; the negative GRAVITY contribution)
      clock_negative = kin_euclid  - kin_signed   (>=0; the negative OSCILLATION contribution)
    """
    energy_decompose(tf, e_kin, e_kins, e_ce, e_cs, e_pot, e_boost,
                     c_amrs, dt_rs, ldg_a, ldg_b, ldg_c, v0)
    kin = float(e_kin.to_numpy().sum())
    kins = float(e_kins.to_numpy().sum())
    ce = float(e_ce.to_numpy().sum())
    cs = float(e_cs.to_numpy().sum())
    pot = float(e_pot.to_numpy().sum())
    gem = float(e_boost.to_numpy().sum())
    H_euclid = kin + ce + pot
    H_signed = kins + cs + pot
    return dict(tag=tag, kin_euclid=kin, kin_signed=kins, curv_euclid=ce, curv_signed=cs,
                potential=pot, boost_GEM=gem, clock_negative=(kin - kins),
                H_euclid=H_euclid, H_signed=H_signed)


def print_ledger(L):
    print(f"  [{L['tag']}]")
    print(f"    kin_euclid     = {L['kin_euclid']:+.6e}   (Frobenius, positive)")
    print(f"    kin_signed     = {L['kin_signed']:+.6e}   (Minkowski; time-axis clock < 0)")
    print(f"    curv_euclid    = {L['curv_euclid']:+.6e}   (Frobenius, positive)")
    print(f"    curv_signed    = {L['curv_signed']:+.6e}   (Minkowski; lower = boost dip)")
    print(f"    potential V-v0 = {L['potential']:+.6e}")
    print(f"    boost_GEM      = {L['boost_GEM']:+.6e}   (negative gravity: curv_e - curv_s)")
    print(f"    clock_negative = {L['clock_negative']:+.6e}   (negative oscill.: kin_e - kin_s)")
    print(f"    H_euclid       = {L['H_euclid']:+.6e}")
    print(f"    H_signed       = {L['H_signed']:+.6e}   <- the physical rest energy")


# ----------------------------------------------------------------------------
# Seeding + minimization (gradient flow on the full 4x4 field)
# ----------------------------------------------------------------------------
def seed(b_star, kick, r0=None, dlt=None):
    seeds.seed_dressed_hedgehog_M(tf, C, C, C, r0 if r0 else r0_vox, rhoc_vox,
                                  dlt if dlt is not None else delta, b_star, rw_vox, kick)


def setup_constrained():
    """The launcher's Option-B init: act mask (interior, off-core, off-axis) + sym basis."""
    idx = np.indices((N, N, N)).astype(np.float32)
    r_vox = np.sqrt((idx[0] - C) ** 2 + (idx[1] - C) ** 2 + (idx[2] - C) ** 2)
    rho_vox = np.sqrt((idx[0] - C) ** 2 + (idx[1] - C) ** 2)
    act = np.zeros((N, N, N), np.float32)
    act[2:-2, 2:-2, 2:-2] = 1.0
    act *= (r_vox > 6.0) * (rho_vox > 3.0)
    tf.act4d.from_numpy(act)
    basis = np.zeros((10, 4, 4), np.float32)
    for a_ in range(4):
        basis[a_, a_, a_] = 1.0
    bi = 4
    for a_ in range(4):
        for c_ in range(a_ + 1, 4):
            basis[bi, a_, c_] = basis[bi, c_, a_] = 1.0 / np.sqrt(2.0)
            bi += 1
    tf.sym_basis.from_numpy(basis)
    n_act_pl = float(act[N // 2].sum() + act[:, N // 2].sum() + act[:, :, N // 2].sum())
    return n_act_pl


def gradient_flow_euclid(n_steps, hist=None):
    """Overdamped gradient flow: zero velocity each step => pure steepest-descent on the
    EUCLIDEAN energy (V on). Relaxes the Faber ansatz to its true minimum (Duda's point).
    stable_mask=0 => Euclidean flux (robust). Returns H_signed history."""
    tf.stable_mask.fill(0.0)
    pde.compute_tstar(tf)
    n_pl = float(N * N + N * N + N * N)
    for step in range(n_steps):
        tf.M_prev_am.copy_from(tf.M_am)      # zero velocity -> GD step
        pde.compute_curvature_flux_4d(tf)
        pde.sample_v03_drift(tf, dt_rs)
        vm = [tf.v03_sums[a_] / n_pl for a_ in range(3)]
        pde.evolve_M_4d(tf, c_amrs, dt_rs, ldg_c, vm[0], vm[1], vm[2], P["KM_INERTIA_4D"])
        tf.swap_matrix_buffers()
        if hist is not None and step % 10 == 0:
            tf.M_prev_am.copy_from(tf.M_am)  # static read (kinetic=0) for clean H
            hist.append(ledger(f"flow{step}")["H_signed"])
    return hist


def clock_run(n_steps, n_act_pl, record_every=20):
    """Activate the de Broglie clock with the validated CONSTRAINED integrator (the
    B-1-validated dt scale 0.007) and record the full signed ledger over time. With the
    clock running, the velocity Mdot has (alpha,0) components, so kin_signed picks up the
    negative oscillation contribution. Returns a list of per-sample ledgers + a blow-up flag."""
    dt_rs_b = 0.007 * dx_am * 0.95 / (c_amrs * 3 ** 0.5)   # DT_SCALE_4D=0.007 (validated)
    dt_eff_b = c_amrs * dt_rs_b
    cc4d = 0.5 * ldg_c / (c_amrs * c_amrs)
    pde.compute_tstar(tf)
    pde.init_P_4d(tf, P["CLOCK_KICK"])
    series = []
    blew_up = False
    for n in range(n_steps):
        pde.flux_4d_constrained(tf)
        pde.update_P_4d(tf, dt_eff_b, cc4d)
        pde.sample_p03_drift(tf)
        npl = max(n_act_pl, 1.0)
        pde.apply_p03_clamp(tf, tf.v03_sums[0] / npl, tf.v03_sums[1] / npl,
                            tf.v03_sums[2] / npl)
        pde.solve_constrained_4d(tf)
        pde.update_M_4d_constrained(tf, dt_eff_b)
        tf.swap_matrix_buffers()
        if n % record_every == 0:
            mmax = float(np.abs(tf.M_am.to_numpy()).max())
            if not np.isfinite(mmax) or mmax > 50.0:
                print(f"    [guard] blow-up at step {n} (max|M|={mmax:.2e}); bailing")
                blew_up = True
                break
            L = ledger(f"clk{n}")
            L["step"] = n
            L["max_M"] = mmax
            series.append(L)
    return series, blew_up


# ----------------------------------------------------------------------------
# MAIN: the electron rest-energy ledger
# ----------------------------------------------------------------------------
def main():
    t0 = time.time()
    results = dict(params=P, derived=dict(N=N, dx_am=dx_am, dt_rs=dt_rs,
                                          ldg_a=ldg_a, ldg_b=ldg_b, ldg_c=ldg_c, v0=v0,
                                          r0_vox=r0_vox, rhoc_vox=rhoc_vox, rw_vox=rw_vox))

    # --- (1) static dressed seed: the Faber + boost STARTING point ----------------
    print("\n[1] static dressed seed (Faber + boost), kick=0 -> the energy-min STARTING point")
    seed(P["B_STAR"], 0.0)
    L_seed = ledger("seed_dressed")
    print_ledger(L_seed)
    results["seed_dressed"] = L_seed

    # control: UNDRESSED (b*=0) seed -> shows the boost-GEM dip is real (signed drops)
    print("\n[1b] control: UNDRESSED seed (b*=0) -> boost_GEM should ~vanish")
    seed(0.0, 0.0)
    L_undressed = ledger("seed_undressed")
    print_ledger(L_undressed)
    results["seed_undressed"] = L_undressed
    print(f"    => boost_GEM(dressed)={L_seed['boost_GEM']:.3e} vs "
          f"(undressed)={L_undressed['boost_GEM']:.3e}  "
          f"[dressed should be >> undressed = the gravity dip is real]")

    # --- (2) MINIMIZE: gradient flow relaxes the Faber ansatz (Duda's core point) -
    print(f"\n[2] gradient-flow minimization ({P['FLOW_STEPS']} steps) - Faber is only the START")
    seed(P["B_STAR"], 0.0)
    H_hist = []
    gradient_flow_euclid(P["FLOW_STEPS"], hist=H_hist)
    tf.M_prev_am.copy_from(tf.M_am)
    L_min = ledger("minimized")
    print_ledger(L_min)
    results["minimized"] = L_min
    results["flow_history"] = H_hist
    drop = (L_seed["H_signed"] - L_min["H_signed"]) / abs(L_seed["H_signed"]) * 100 if L_seed["H_signed"] else 0
    print(f"    => H_signed dropped {drop:+.1f}% from the fixed-ansatz seed to the minimum")

    # --- (3) clock activation: the constrained integrator, signed kinetic ----------
    print(f"\n[3] clock activation ({P['CLOCK_STEPS']} constrained steps) - the oscillation term")
    n_act_pl = setup_constrained()
    seed(P["B_STAR"], P["CLOCK_KICK"])
    series, blew_up = clock_run(P["CLOCK_STEPS"], n_act_pl)
    H_static = L_min["H_signed"]
    results["clock_series"] = series
    results["clock_blew_up"] = blew_up
    if series:
        # average over the second half (after the kick transient)
        half = series[len(series) // 2:]
        H_clock_avg = float(np.mean([s["H_signed"] for s in half]))
        kin_s_avg = float(np.mean([s["kin_signed"] for s in half]))
        clk_neg_avg = float(np.mean([s["clock_negative"] for s in half]))
        print(f"    sampled {len(series)} ledgers; clock-active (2nd-half mean):")
        print(f"      H_signed       = {H_clock_avg:+.6e}  (static minimized = {H_static:+.6e})")
        print(f"      kin_signed      = {kin_s_avg:+.6e}  (clock oscillation, signed)")
        print(f"      clock_negative  = {clk_neg_avg:+.6e}  (Euclid - signed kinetic = the dip)")
        if blew_up:
            print("      [!] run flagged a blow-up before completing; treat the clock term as")
            print("          regularization-limited (Duda's open Higgs problem), reported honestly.")
    else:
        H_clock_avg = float("nan")
        print("    no stable clock samples (immediate blow-up) - clock term regularization-limited")
    results["H_clock_avg"] = H_clock_avg

    # --- (4) electron calibration: pin e_scale (on the stable static minimum) ------
    # Calibrate on the MINIMIZED static signed energy (robust); the clock term is a
    # measured correction, not the calibration anchor (it is regularization-sensitive).
    H_rest = H_static
    e_scale = P["M_E_MEV"] / H_rest if H_rest else float("nan")
    results["H_rest_electron"] = H_rest
    results["e_scale_MeV_per_unit"] = e_scale
    print(f"\n[4] electron calibration on the static minimum: H_rest={H_rest:.6e}"
          f" -> e_scale={e_scale:.6e} MeV/unit")

    # --- plots --------------------------------------------------------------------
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
    if H_hist:
        ax[0].plot(np.arange(len(H_hist)) * 10, H_hist, "-o", ms=3)
        ax[0].set_title("gradient-flow minimization (H_signed)")
        ax[0].set_xlabel("flow step")
        ax[0].set_ylabel("H_signed")
        ax[0].grid(alpha=0.3)
    if series:
        steps = [s["step"] for s in series]
        ax[1].plot(steps, [s["H_signed"] for s in series], "-o", ms=3, label="H_signed")
        ax[1].plot(steps, [s["kin_signed"] for s in series], "-s", ms=3, label="kin_signed")
        ax[1].axhline(H_static, color="k", ls="--", lw=1, label="static minimized")
        ax[1].set_title("clock activation (constrained)")
        ax[1].set_xlabel("constrained step")
        ax[1].legend()
        ax[1].grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, f"m5_9_4_ledger{_SUFFIX}.png"), dpi=110)
    print(f"\n  plot -> m5_9_4_ledger{_SUFFIX}.png")

    results["wall_seconds"] = time.time() - t0
    with open(os.path.join(HERE, f"m5_9_4_results{_SUFFIX}.json"), "w") as f:
        json.dump(results, f, indent=2, default=float)
    print(f"  json -> m5_9_4_results{_SUFFIX}.json   ({results['wall_seconds']:.1f}s)")

    # --- summary table ------------------------------------------------------------
    print("\n" + "=" * 84)
    print("LEDGER SUMMARY (signed = physical rest energy):")
    print(f"  {'config':>18} {'H_euclid':>14} {'H_signed':>14} {'boost_GEM':>12} {'clock_neg':>12}")
    for key in ("seed_undressed", "seed_dressed", "minimized"):
        L = results[key]
        print(f"  {key:>18} {L['H_euclid']:>14.4e} {L['H_signed']:>14.4e}"
              f" {L['boost_GEM']:>12.4e} {L['clock_negative']:>12.4e}")
    print(f"  {'clock_time_avg':>18} {'-':>14} {H_clock_avg:>14.4e} {'-':>12} {'-':>12}")
    print("=" * 84)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
