"""
M5.9.5 - on the production 4x4 engine: (A) the Higgs / core-volume confiner scale-selection
search, (B) the #220 BONUS - the DYNAMICAL clock-rate scaling omega(r0), (C) the #200 lepton
spectrum attempt + the honest hierarchy-origin verdict.

This is sub-tasks 4 + 5 of the corrected #200 plan, built on the signed-energy ledger from
m5_9_4 (which established the boost-GEM negative gravity term + the signed kinetic clock term).

(A) HIGGS / CONFINER (sub-task 4). The toy m5_9_3 showed the 3x3 Faber functional is
    scale-balanced (E ~ 1/r0, r0 a free modulus) and that a fixed-density confiner
    E_conf = B*integral (1-s^2)^3 ~ B*r0^3 selects a scale (interior minimum). Here we test
    the SAME idea on the real 4x4 signed energy: scan r0, measure the signed curvature E(r0)
    and the confiner integral, and check whether curv_signed + B*conf has an interior minimum
    (scale-selection), now WITH the boost (gravity) term present.

(B) #220 CLOCK SCALING (the bonus). #220 proved E*r0 = const (scale-covariance) ANALYTICALLY
    for the V-on Faber soliton, so omega ~ 1/r0 ~ m as a PRINCIPLE; the residual it left open
    was a DIRECT dynamical omega(r0) readout (the V-on dynamical stack was unbuilt). We now have
    that stack (the constrained integrator). We run the clock at several r0, FFT a field probe
    to get the natural clock frequency omega, and check omega*r0 ~ const.

(C) #200 SPECTRUM (sub-task 5). Test whether the three leptons emerge as eigen-configurations
    (scan the biaxial delta + the scale). Honest expectation from the toy: the mass-vs-scale law
    is clean but the discrete hierarchy ORIGIN (why three, why 206.8 / 3477) stays Yukawa input.
    Report how far off, not a forced match.

Headless CPU Taichi + matplotlib Agg. Run:  python m5_9_5_higgs_clock_spectrum.py
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

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))

# --- parameters (mirror m5_9_4 / the production seed) -----------------------------
EDGE = 1e-15
TARGET_VOXELS = 48 ** 3
RHOC_VOXELS = 3.0
DELTA = 0.30
B_STAR = 0.13
RW_FRACTION = 0.29
C_AMRS = 0.2998
KM = 30.0
LDG_K = 1.0
M_E, M_MU, M_TAU = 0.510999, 105.658, 1776.86
R_MU, R_TAU = M_MU / M_E, M_TAU / M_E

tf = medium.TensorField([EDGE, EDGE, EDGE], TARGET_VOXELS, [0.5, 0.5, 0.5], viz_stride=2)
N = tf.nx
C = N // 2
GRID = tf.grid_size
dx_am = float(tf.dx_am)
dt_rs = 0.5 * dx_am * 0.95 / (C_AMRS * 3 ** 0.5)
ldg_c = LDG_K * C_AMRS ** 2 / dx_am ** 4
ldg_a = -2.0 * ldg_c * (1.0 + DELTA ** 2)
tr2_vac = 1.0 + DELTA ** 2
v0 = ldg_a * tr2_vac + ldg_c * tr2_vac ** 2
rw_vox = RW_FRACTION * N

e_kin = ti.field(ti.f32, shape=GRID)
e_kins = ti.field(ti.f32, shape=GRID)
e_ce = ti.field(ti.f32, shape=GRID)
e_cs = ti.field(ti.f32, shape=GRID)
e_pot = ti.field(ti.f32, shape=GRID)

# grid geometry (numpy) for the confiner integral
_idx = np.indices((N, N, N)).astype(np.float64)
_r = np.sqrt((_idx[0] - C) ** 2 + (_idx[1] - C) ** 2 + (_idx[2] - C) ** 2)
_interior = np.zeros((N, N, N), bool)
_interior[1:-1, 1:-1, 1:-1] = True

print("=" * 84)
print("M5.9.5 - Higgs/confiner + clock scaling omega(r0) + lepton spectrum")
print("=" * 84)
print(f"grid {N}^3  dx_am={dx_am:.4e}  ldg_c={ldg_c:.3e}  delta={DELTA}")


@ti.kernel
def energy_decompose(tfx: ti.template(),  # type: ignore
                     kin: ti.template(), kins: ti.template(),  # type: ignore
                     ce: ti.template(), cs: ti.template(),  # type: ignore
                     pot: ti.template(),  # type: ignore
                     c_a: ti.f32, dt: ti.f32, a: ti.f32, cc: ti.f32, vv0: ti.f32):  # type: ignore
    nx, ny, nz = tfx.nx, tfx.ny, tfx.nz
    inv_dt = 1.0 / dt
    c2 = c_a * c_a
    inv_2dx = 1.0 / (2.0 * tfx.dx_am)
    for i, j, k in ti.ndrange((1, nx - 1), (1, ny - 1), (1, nz - 1)):
        m = tfx.M_am[i, j, k]
        m_dot = (m - tfx.M_prev_am[i, j, k]) * inv_dt
        kin[i, j, k] = 0.5 * m_dot.norm_sqr()
        kins[i, j, k] = 0.5 * pde.signed_dot4(m_dot, m_dot)
        mx = (tfx.M_am[i + 1, j, k] - tfx.M_am[i - 1, j, k]) * inv_2dx
        my = (tfx.M_am[i, j + 1, k] - tfx.M_am[i, j - 1, k]) * inv_2dx
        mz = (tfx.M_am[i, j, k + 1] - tfx.M_am[i, j, k - 1]) * inv_2dx
        cxy = pde.commutator(mx, my)
        cxz = pde.commutator(mx, mz)
        cyz = pde.commutator(my, mz)
        ce[i, j, k] = c2 * 4.0 * (cxy.norm_sqr() + cxz.norm_sqr() + cyz.norm_sqr())
        cs[i, j, k] = c2 * 4.0 * (pde.signed_dot4(cxy, cxy) + pde.signed_dot4(cxz, cxz)
                                  + pde.signed_dot4(cyz, cyz))
        pot[i, j, k] = pde.V_M(m, a, 0.0, cc) - vv0


def seed(r0_vox, b_star=B_STAR, kick=0.0, dlt=DELTA):
    seeds.seed_dressed_hedgehog_M(tf, C, C, C, r0_vox, RHOC_VOXELS, dlt, b_star, rw_vox, kick)


def static_ledger(r0_vox, b_star=B_STAR, dlt=DELTA):
    """Signed-energy ledger of the STATIC dressed hedgehog at scale r0_vox."""
    seed(r0_vox, b_star=b_star, kick=0.0, dlt=dlt)
    tf.M_prev_am.copy_from(tf.M_am)  # static => kinetic 0
    energy_decompose(tf, e_kin, e_kins, e_ce, e_cs, e_pot, C_AMRS, dt_rs, ldg_a, ldg_c, v0)
    cs = float(e_cs.to_numpy().sum())
    ce = float(e_ce.to_numpy().sum())
    pot = float(e_pot.to_numpy().sum())
    return dict(curv_signed=cs, curv_euclid=ce, potential=pot, H_signed=cs + pot)


def confiner_integral(r0_vox):
    """B-coefficient integrand: integral (1 - s^2)^3 over the interior, s = r/sqrt(r0^2+r^2).
    A FIXED coefficient on this term ~ r0^3 (the m5_9_3 toy result), breaking Faber scale-
    covariance and selecting r0. Computed analytically from the grid geometry."""
    s2 = _r ** 2 / (r0_vox ** 2 + _r ** 2)
    integ = ((1.0 - s2) ** 3)[_interior].sum()
    return float(integ)


def fit_pow(x, y):
    lx, ly = np.log(np.asarray(x, float)), np.log(np.abs(np.asarray(y, float)) + 1e-300)
    p, b = np.polyfit(lx, ly, 1)
    yhat = p * lx + b
    ss = np.sum((ly - yhat) ** 2)
    st = np.sum((ly - ly.mean()) ** 2)
    return float(p), float(1.0 - ss / st if st > 0 else np.nan)


def setup_constrained():
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


def clock_frequency(r0_vox, n_act_pl, n_steps=400, kick=0.05):
    """Run the constrained clock and FFT field probes -> the NATURAL clock frequency omega.
    Probe must be INSIDE the act mask (r_vox>6), else Mdot=0 there. We record the 6 clock-
    active components M[0,1],M[0,2],M[1,2],M[0,3],M[1,3],M[2,3] at an active equatorial voxel
    and FFT each; the peak with the most spectral power gives omega (the (alpha,3) components
    carry the time-axis de Broglie clock; the (i,j) carry the biaxial frame rotation).
    Returns (omega_peak, probe_amp, blew_up)."""
    seed(r0_vox, b_star=B_STAR, kick=kick)
    dt_rs_b = 0.007 * dx_am * 0.95 / (C_AMRS * 3 ** 0.5)
    dt_eff_b = C_AMRS * dt_rs_b
    cc4d = 0.5 * ldg_c / (C_AMRS * C_AMRS)
    pde.compute_tstar(tf)
    pde.init_P_4d(tf, kick)
    px, py, pz = C + 8, C, C   # equatorial probe INSIDE the act region (r_vox=8 > 6)
    comps = [(0, 1), (0, 2), (1, 2), (0, 3), (1, 3), (2, 3)]
    series = [[] for _ in comps]
    blew_up = False
    for n in range(n_steps):
        pde.flux_4d_constrained(tf)
        pde.update_P_4d(tf, dt_eff_b, cc4d)
        pde.sample_p03_drift(tf)
        npl = max(n_act_pl, 1.0)
        pde.apply_p03_clamp(tf, tf.v03_sums[0] / npl, tf.v03_sums[1] / npl, tf.v03_sums[2] / npl)
        pde.solve_constrained_4d(tf)
        pde.update_M_4d_constrained(tf, dt_eff_b)
        tf.swap_matrix_buffers()
        mval = tf.M_am[px, py, pz]
        for ci, (a_, b_) in enumerate(comps):
            series[ci].append(float(mval[a_, b_]))
        if n % 100 == 99:
            mmax = float(np.abs(tf.M_am.to_numpy()).max())
            if not np.isfinite(mmax) or mmax > 50.0:
                blew_up = True
                break
    nrec = len(series[0])
    if nrec < 32:
        return float("nan"), float("nan"), blew_up
    win = np.hanning(nrec)
    freqs = np.fft.rfftfreq(nrec, d=dt_eff_b)
    best_pow, best_om, best_amp = -1.0, float("nan"), float("nan")
    for ci in range(len(comps)):
        s = np.asarray(series[ci])
        s = s - s.mean()
        amp = float(np.abs(s).max())
        if amp < 1e-9:
            continue
        spec = np.abs(np.fft.rfft(s * win))
        spec[0] = 0.0
        pk = int(np.argmax(spec))
        if spec[pk] > best_pow:
            best_pow = float(spec[pk])
            best_om = float(freqs[pk]) * 2.0 * np.pi
            best_amp = amp
    return best_om, best_amp, blew_up


def main():
    t0 = time.time()
    out = dict(grid=N, dx_am=dx_am, ldg_c=ldg_c, delta=DELTA)

    # ============ (A) Higgs / confiner scale-selection ============================
    print("\n[A] Higgs / confiner scale-selection (signed energy vs r0)")
    r0_fracs = np.array([0.04, 0.05, 0.06, 0.08, 0.10, 0.13])
    r0_voxs = r0_fracs * N
    curv_s, pots, confs = [], [], []
    for r0v in r0_voxs:
        L = static_ledger(r0v)
        cf = confiner_integral(r0v)
        curv_s.append(L["curv_signed"])
        pots.append(L["potential"])
        confs.append(cf)
        print(f"    r0_vox={r0v:5.2f}: curv_signed={L['curv_signed']:.4e} "
              f"V={L['potential']:.4e} conf_int={cf:.4e}")
    curv_s = np.array(curv_s); pots = np.array(pots); confs = np.array(confs)
    p_curv, r2_curv = fit_pow(r0_voxs, curv_s)
    p_conf, r2_conf = fit_pow(r0_voxs, confs)
    print(f"    -> curv_signed ~ r0^{p_curv:.2f} (R2={r2_curv:.3f}); "
          f"conf_int ~ r0^{p_conf:.2f} (R2={r2_conf:.3f})")
    # choose B so curv + B*conf has a minimum inside the window; report selected r0*
    sel = {}
    for Bfac in (0.5, 1.0, 2.0):
        # set B so the confiner ~ matches curv at the mid scale, times Bfac
        mid = len(r0_voxs) // 2
        B = Bfac * curv_s[mid] / confs[mid]
        Etot = curv_s + B * confs
        imin = int(np.argmin(Etot))
        interior = 0 < imin < len(r0_voxs) - 1
        sel[f"B={Bfac}x"] = dict(B=B, r0_star=float(r0_voxs[imin]), interior=interior,
                                 E_star=float(Etot[imin]))
        print(f"    confiner B={Bfac}x: r0*={r0_voxs[imin]:.2f} interior_min={interior}")
    out["confiner"] = dict(r0_voxs=r0_voxs.tolist(), curv_signed=curv_s.tolist(),
                           confiner=confs.tolist(), potential=pots.tolist(),
                           p_curv=p_curv, p_conf=p_conf, selections=sel)

    # ============ (B) #220 BONUS: dynamical clock frequency omega(r0) =============
    print("\n[B] #220 BONUS: dynamical clock frequency omega(r0) - does omega*r0 ~ const?")
    n_act_pl = setup_constrained()
    r0_clock = np.array([0.05, 0.07, 0.10, 0.13]) * N   # wider span so omega can separate
    omegas, blews = [], []
    for r0v in r0_clock:
        om, amp, blew = clock_frequency(r0v, n_act_pl, n_steps=1200)  # longer => resolve period
        omegas.append(om); blews.append(blew)
        print(f"    r0_vox={r0v:5.2f}: omega={om:.4e}  omega*r0={om * r0v:.4e}  "
              f"blew_up={blew}")
    omegas = np.array(omegas)
    valid = np.isfinite(omegas) & (~np.array(blews)) & (omegas > 0)
    if valid.sum() >= 2:
        p_om, r2_om = fit_pow(r0_clock[valid], omegas[valid])
        prod = omegas[valid] * r0_clock[valid]
        cv = float(np.std(prod) / np.mean(prod)) if np.mean(prod) else float("nan")
        print(f"    -> omega ~ r0^{p_om:.2f} (R2={r2_om:.3f}); omega*r0 CV={cv * 100:.1f}% "
              f"(scale-covariance => p~-1, CV~0)")
    else:
        p_om, r2_om, cv = float("nan"), float("nan"), float("nan")
        print("    -> insufficient valid clock measurements")
    out["clock_scaling"] = dict(r0_voxs=r0_clock.tolist(), omega=omegas.tolist(),
                                blew_up=[bool(b) for b in blews], p_omega=p_om, cv_omega_r0=cv)

    # ============ (C) #200 spectrum: delta scan + hierarchy verdict ===============
    print("\n[C] #200 spectrum attempt: does varying the biaxial structure give 206.8 / 3477?")
    deltas = np.array([0.15, 0.30, 0.50, 0.70])
    Hs = []
    for dlt in deltas:
        L = static_ledger(0.06 * N, dlt=dlt)
        Hs.append(L["H_signed"])
        print(f"    delta={dlt:.2f}: H_signed={L['H_signed']:.4e}")
    Hs = np.array(Hs)
    ratios = Hs / Hs[1]   # normalize to delta=0.30
    print(f"    H_signed ratios (norm @ delta=0.30): {np.round(ratios, 3)}")
    print(f"    target lepton ratios: m_mu/m_e={R_MU:.1f}, m_tau/m_e={R_TAU:.1f}")
    span = float(Hs.max() / Hs.min())
    print(f"    => delta-knob energy span = {span:.2f}x  (need {R_TAU:.0f}x for e->tau)")
    out["spectrum"] = dict(deltas=deltas.tolist(), H_signed=Hs.tolist(),
                           span=span, R_mu=R_MU, R_tau=R_TAU)

    # --- plots --------------------------------------------------------------------
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.5))
    ax[0].loglog(r0_voxs, curv_s, "-o", label="curv_signed")
    ax[0].loglog(r0_voxs, confs / confs[0] * curv_s[0], "-s", label="confiner (scaled)")
    ax[0].set_title(f"(A) confiner: curv~r0^{p_curv:.2f}, conf~r0^{p_conf:.2f}")
    ax[0].set_xlabel("r0_vox"); ax[0].legend(); ax[0].grid(alpha=0.3, which="both")
    if valid.sum() >= 2:
        ax[1].loglog(r0_clock[valid], omegas[valid], "-o")
        ax[1].set_title(f"(B) clock omega~r0^{p_om:.2f} (CV {cv*100:.1f}%)")
    ax[1].set_xlabel("r0_vox"); ax[1].set_ylabel("omega"); ax[1].grid(alpha=0.3, which="both")
    ax[2].semilogy(deltas, Hs, "-o")
    ax[2].set_title(f"(C) delta scan: span {span:.2f}x (need {R_TAU:.0f}x)")
    ax[2].set_xlabel("biaxial delta"); ax[2].set_ylabel("H_signed"); ax[2].grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "m5_9_5_results.png"), dpi=110)

    out["wall_seconds"] = time.time() - t0
    with open(os.path.join(HERE, "m5_9_5_results.json"), "w") as f:
        json.dump(out, f, indent=2, default=float)
    print(f"\n  plot -> m5_9_5_results.png ; json -> m5_9_5_results.json "
          f"({out['wall_seconds']:.1f}s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
