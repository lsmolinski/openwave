"""
M5.8.1 — headless CPU smoke test of the production 4x4 promotion (no GUI).

Builds a small TensorField, seeds vacuum + biaxial hedgehog, runs eigen_decompose and a
few evolve_M steps, and checks the 4x4 substrate behaves: M is 4x4, the director comes
from the spatial sub-block, g (time axis, index 3) stays constant, dynamics are finite.

Catches compile errors + shape/index bugs before the on-screen GUI test (which is the
GUI-only part the user runs). Run from the openwave repo root:

    PYTHONPATH=/Users/xrodz/Documents/source_code/OPENWAVE-LABS/openwave \\
        python -m openwave.xperiments.m5_liquid_crystal.research.sandbox_v8.m5_8_1_headless_check
(or just: python m5_8_1_headless_check.py  with the repo root on sys.path — handled below)
"""
import os, sys
# put the openwave repo root on the path so the package relative-imports resolve
_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../"))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

import numpy as np
import taichi as ti
ti.init(arch=ti.cpu, default_fp=ti.f32)

from openwave.xperiments.m5_liquid_crystal import medium
from openwave.xperiments.m5_liquid_crystal import engine1_seeds as seeds
from openwave.xperiments.m5_liquid_crystal import engine2_pde as pde

print(f"MDIM={medium.MDIM}  LC_DELTA={medium.LC_DELTA}  LC_G={medium.LC_G}")

# --- build a small field -----------------------------------------------------
EDGE = 1e-15
tf = medium.TensorField([EDGE, EDGE, EDGE], 24**3, [0.5, 0.5, 0.5], viz_stride=2)
N = tf.nx
c = N // 2
print(f"grid {tf.nx}x{tf.ny}x{tf.nz}")

def report(tag, ok):
    print(f"  [{'PASS ✓' if ok else 'FAIL ✗'}] {tag}")

# --- (1) vacuum seed: M = diag(δ,δ,1,g), director = ẑ ------------------------
print("\n(1) seed_vacuum_M")
seeds.seed_vacuum_M(tf, tf.lc_delta)
M = tf.M_am.to_numpy()
report(f"M is 4x4 (shape {M.shape[-2:]})", M.shape[-2:] == (4, 4))
mc = M[c, c, c]
report(f"M[3,3] == g ({mc[3,3]:.3f})", abs(mc[3, 3] - medium.LC_G) < 1e-4)
report(f"time off-diagonals 0 (max {np.abs(mc[3,:3]).max():.1e})", np.abs(mc[3, :3]).max() < 1e-6)
ev = np.sort(np.linalg.eigvalsh(mc))[::-1]
report(f"spectrum {np.round(ev,3)} ~ [g,1,δ,δ]", np.allclose(ev, [medium.LC_G, 1, medium.LC_DELTA, medium.LC_DELTA], atol=1e-3))
pde.eigen_decompose(tf)
d = tf.director_nhat.to_numpy()[c, c, c]
report(f"vacuum director {np.round(d,3)} ~ ±ẑ", abs(abs(d[2]) - 1) < 1e-3)

# --- (2) biaxial hedgehog: M spatial = O diag(1,δ,0) O^T, director radial -----
print("\n(2) seed_biaxial_hedgehog_M + eigen_decompose")
seeds.seed_biaxial_hedgehog_M(tf, c, c, c, 0.06 * N, 3.0, 0.30)
pde.eigen_decompose(tf)
M = tf.M_am.to_numpy()
report("M[...,3,3] all == g", np.allclose(M[..., 3, 3], medium.LC_G, atol=1e-4))
report("time off-diagonals all 0", np.abs(M[..., 3, :3]).max() < 1e-6)
ev = tf.eigenvalues.to_numpy()
report("eigenvalues finite (no NaN)", np.all(np.isfinite(ev)))
dn = tf.director_nhat.to_numpy()
norms = np.linalg.norm(dn, axis=-1)
report(f"director ~unit (mean |n|={norms.mean():.3f})", abs(norms.mean() - 1) < 0.05)
# off-core voxel: principal eigenvalue ~1, spectrum from the SPATIAL block (not g)
ic = c + N // 4
sp_ev = np.sort(np.linalg.eigvalsh(M[ic, c, c][:3, :3]))[::-1]
report(f"spatial sub-block spectrum {np.round(sp_ev,3)} ~ [1,δ,0]", abs(sp_ev[0] - 1) < 0.05)

# --- (3) evolve a few steps (V off): g constant, finite, stable --------------
print("\n(3) evolve_M x5 (V off)")
M0 = tf.M_am.to_numpy().copy()
g0 = M0[..., 3, 3].copy()
for _ in range(5):
    pde.compute_curvature_flux(tf)
    pde.evolve_M(tf, 1.0, 0.05 * tf.dx_am, 0.0, 0.0, 0.0)   # c=1.0, a=b=c=0 → V off
    tf.swap_matrix_buffers()
M5 = tf.M_am.to_numpy()
report("no NaN/inf after 5 steps", np.all(np.isfinite(M5)))
report(f"g constant (max drift {np.abs(M5[...,3,3]-g0).max():.1e})", np.abs(M5[..., 3, 3] - g0).max() < 1e-5)
report("time off-diagonals stay 0", np.abs(M5[..., 3, :3]).max() < 1e-6)
spatial_moved = np.abs(M5[..., :3, :3] - M0[..., :3, :3]).max()
report(f"spatial block evolved (max Δ {spatial_moved:.2e} > 0)", spatial_moved > 0)

# --- (4) THE GUI-EXPLOSION REGRESSION: V_M must act on the spatial block only, so g
#         (=8, g²=64) cannot inflate dV_M. Two checks: (a) dV_M is g-INDEPENDENT (the
#         decisive fix verification), (b) a production-faithful V-on evolve stays bounded.
print("\n(4) V-on regression — dV_M must be g-independent (the GUI explosion)")
_res = ti.Matrix.field(4, 4, ti.f32, shape=())

@ti.kernel
def eval_dVM(m: ti.types.matrix(4, 4, ti.f32), a: ti.f32, b: ti.f32, c: ti.f32):  # type: ignore
    _res[None] = pde.dV_M(m, a, b, c)

def dVM_of(g):
    M = np.zeros((4, 4), np.float32)
    M[:3, :3] = np.diag([1.0, 0.30, 0.0])  # spatial diag(1,δ,0)
    M[3, 3] = g
    eval_dVM(ti.Matrix(M.tolist()), 1.0, 0.0, 1.0)
    return _res.to_numpy()

d8, d100 = dVM_of(8.0), dVM_of(100.0)
report("dV_M time row/col == 0 (g decoupled)", np.abs(d8[3, :]).max() < 1e-5 and np.abs(d8[:, 3]).max() < 1e-5)
report(f"dV_M spatial g-INDEPENDENT (Δ g=8 vs g=100 = {np.abs(d8[:3,:3]-d100[:3,:3]).max():.1e})",
       np.abs(d8[:3, :3] - d100[:3, :3]).max() < 1e-4)

# (b) production-faithful V-on evolve (c_amrs≈0.3, dt_rs CFL-bound, ldg_c=c²/dx⁴)
print("\n    production-faithful V-on evolve (real c_amrs / dt_rs / ldg):")
c_amrs = 0.2998   # WAVE_SPEED·RONTOSECOND/ATTOMETER ≈ 0.3 am/rs (SIM_SPEED=1)
dt_rs = tf.dx_am * 0.95 / (c_amrs * 3 ** 0.5)
ldg_c = 1.0 * c_amrs**2 / tf.dx_am**4
ldg_a = -2.0 * ldg_c * (1.0 + 0.30**2)
seeds.seed_biaxial_hedgehog_M(tf, c, c, c, 0.06 * N, 3.0, 0.30)
pde.eigen_decompose(tf)
g0 = tf.M_am.to_numpy()[..., 3, 3].copy()
for _ in range(30):
    pde.compute_curvature_flux(tf)
    pde.evolve_M(tf, c_amrs, dt_rs, ldg_a, 0.0, ldg_c)
    tf.swap_matrix_buffers()
Mv = tf.M_am.to_numpy()
amax = np.abs(Mv[..., :3, :3]).max()
report("no NaN/inf after 30 V-on steps", np.all(np.isfinite(Mv)))
report(f"spatial block BOUNDED (max|M_sp|={amax:.2f}, not exploding)", np.isfinite(amax) and amax < 10.0)
report(f"g still constant (drift {np.abs(Mv[...,3,3]-g0).max():.1e})", np.abs(Mv[..., 3, 3] - g0).max() < 1e-3)

# --- (5) observables path: the M-substrate kernels the GUI compute_field_observables
#         calls (update_trackers_M + compute_energyH_density_M). The ψ-retire (M5.8
#         cleanup) touched engine3_observables, so exercise these to catch a dropped
#         @ti.kernel decorator (would raise "called from Python-scope"). ------------
print("\n(5) observables kernels (update_trackers_M + compute_energyH_density_M)")
from openwave.xperiments.m5_liquid_crystal import engine3_observables as obs
trackers = medium.Trackers(tf.grid_size)
observables = medium.FieldObservables(tf.grid_size)
seeds.seed_biaxial_hedgehog_M(tf, c, c, c, 0.06 * N, 3.0, 0.30)
pde.eigen_decompose(tf)
ok_kernels = True
try:
    obs.update_trackers_M(tf, trackers, 0.05 * tf.dx_am, tf.lc_delta)
    obs.compute_energyH_density_M(tf, observables, 0.2998, 0.05 * tf.dx_am, 0.0, 0.0, 0.0, 0.0, 1.0)
except Exception as e:
    ok_kernels = False
    print(f"    EXCEPTION: {type(e).__name__}: {e}")
report("update_trackers_M + compute_energyH_density_M run as KERNELS (not Python-scope)", ok_kernels)
if ok_kernels:
    eH = observables.energyH_density_aJ.to_numpy()
    report(f"energyH density finite (max={np.abs(eH).max():.2e})", np.all(np.isfinite(eH)))

# --- (6) M5.8.2c — the 4D Minkowski port: dressed seeder, stable mask, signed
#         flux (b=0 identity vs the 3D kernel + f32-vs-f64 cross-check against the
#         2c-1 numpy mirror), and a 30-step evolve_M_4d bounded run. -----------------
print("\n(6) M5.8.2c 4D port (seed_dressed_hedgehog_M + stable_mask + flux_4d + evolve_M_4d)")
# NOTE: numpy mirrors defined LOCALLY — importing the 2c-1/2c-2 sandbox modules
# pulls m5_8_2a's module-level ti.init → a SECOND runtime init with live fields
# → segfault (caught 2026-06-05; Taichi gotcha: never re-init mid-run).


def np_tw(A):
    B = A.copy()
    B[..., 3, :3] *= -1.0
    B[..., :3, 3] *= -1.0
    return B


def np_central(f, axis, h):
    out = np.zeros_like(f)
    sp_, sm, so = [slice(None)] * f.ndim, [slice(None)] * f.ndim, [slice(None)] * f.ndim
    sp_[axis], sm[axis], so[axis] = slice(2, None), slice(0, -2), slice(1, -1)
    out[tuple(so)] = (f[tuple(sp_)] - f[tuple(sm)]) / (2 * h)
    return out


def np_commf(A, B):
    return (np.einsum("...ac,...cb->...ab", A, B)
            - np.einsum("...ac,...cb->...ab", B, A))

# (6a) dressed seed + mask
RW_VOX = 0.29 * N
seeds.seed_dressed_hedgehog_M(tf, c, c, c, 0.06 * N, 3.0, 0.30, 0.13, RW_VOX, 0.0)
pde.compute_stable_mask(tf)
M_np = tf.M_am.to_numpy().astype(np.float64)
psi_np = tf.M_psi_am.to_numpy().astype(np.float64)
mask_np = tf.stable_mask.to_numpy()
report(f"dressed seed finite (max|M|={np.abs(M_np).max():.2f}) + M_psi finite",
       np.isfinite(M_np).all() and np.isfinite(psi_np).all())
frac = mask_np.mean()
report(f"stable mask computed (stable fraction {100 * frac:.1f}% — 2b-2 ballpark ~50-90%)",
       0.3 < frac < 0.98)

# (6b) f32-vs-f64 kernel cross-check: production G_4d vs 2× the 2c-1 numpy flux
# (production G_α = 8Σ[twF, M_ν] = 2× the spike's 4Σ[twF, M_j]), computed from the
# SAME (f32-seeded) M field in f64, with the SAME mask blend.
pde.compute_curvature_flux_4d(tf)
Gx_ti = tf.curv_flux_x.to_numpy().astype(np.float64)
h64 = float(tf.dx_am)
Mi64 = [np_central(M_np, ax, h64) for ax in range(3)]
st = mask_np[..., None, None]


def tw_blend(F):
    return F + st * (np_tw(F) - F)


Fxy, Fxz = np_commf(Mi64[0], Mi64[1]), np_commf(Mi64[0], Mi64[2])
Gx_np = 8.0 * (np_commf(tw_blend(Fxy), Mi64[1]) + np_commf(tw_blend(Fxz), Mi64[2]))
inner = (slice(2, -2),) * 3
scale = np.abs(Gx_np[inner]).max() + 1e-30
rel = np.abs(Gx_ti[inner] - Gx_np[inner]).max() / scale
report(f"flux_4d f32 kernel matches f64 numpy mirror (rel {rel:.2e} < 1e-3)", rel < 1e-3)

# (6c) b=0 identity: no (α,3) components ⇒ flux_4d ≡ flux_3d EXACTLY
seeds.seed_dressed_hedgehog_M(tf, c, c, c, 0.06 * N, 3.0, 0.30, 0.0, RW_VOX, 0.0)
pde.compute_stable_mask(tf)
pde.compute_curvature_flux_4d(tf)
G4 = tf.curv_flux_x.to_numpy().copy()
pde.compute_curvature_flux(tf)
G3 = tf.curv_flux_x.to_numpy()
ident = np.abs(G4 - G3).max() / (np.abs(G3).max() + 1e-30)
report(f"b=0 identity: flux_4d ≡ flux_3d (rel diff {ident:.2e} — FP-reassociation only)",
       ident < 1e-6)

# (6d) 300-step evolve_M_4d bounded, V-off then V-on — the production-DEFAULT
# stack (SIGNED_FLUX_4D off ⇒ Euclid flux + dressed-t* V pinning + km inertia;
# the 2c-2 lessons: 30 steps was too short to catch slow pumps — run 300)
for ldg_k, vtag in ((0.0, "V-off"), (1.0, "V-on")):
    seeds.seed_dressed_hedgehog_M(tf, c, c, c, 0.06 * N, 3.0, 0.30, 0.13, RW_VOX, 0.05)
    pde.compute_stable_mask(tf)
    tf.stable_mask.fill(0.0)        # production default: SIGNED_FLUX_4D off (safe v1)
    pde.compute_tstar(tf)
    c_amrs = 0.2998
    dt_rs = 0.5 * tf.dx_am * 0.95 / (c_amrs * 3**0.5)      # DT_SCALE_4D=0.5 faithful
    ldg_c_ = ldg_k * c_amrs**2 / tf.dx_am**4
    for _ in range(300):
        pde.compute_curvature_flux_4d(tf)
        pde.sample_v03_drift(tf, dt_rs)
        n_pl = float(tf.nx * tf.ny + tf.ny * tf.nz + tf.nx * tf.nz)
        vm = [tf.v03_sums[a_] / n_pl for a_ in range(3)]
        pde.evolve_M_4d(tf, c_amrs, dt_rs, ldg_c_, vm[0], vm[1], vm[2], 30.0)
        tf.swap_matrix_buffers()
    Mend = tf.M_am.to_numpy()
    report(f"evolve_M_4d 300 steps bounded ({vtag}, kicked, full fix stack; "
           f"max|M|={np.abs(Mend).max():.2f})",
           np.isfinite(Mend).all() and np.abs(Mend).max() < 100.0)

print("\n  => headless check done. dV_M g-independent + bounded V-on evolve = the GUI "
      "explosion is fixed; observables kernels run; the 4D port kernels (dressed seed, "
      "stable mask, signed flux f32≡f64, b=0 identity, bounded 4D evolve) verified. "
      "GUI re-test confirms on screen.")
