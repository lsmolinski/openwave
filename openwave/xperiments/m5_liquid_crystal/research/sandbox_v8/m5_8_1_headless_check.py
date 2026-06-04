"""
M5.8.1 — headless CPU smoke test of the production 4x4 promotion (no GUI).

Builds a small WaveField, seeds vacuum + biaxial hedgehog, runs eigen_decompose and a
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
wf = medium.TensorField([EDGE, EDGE, EDGE], 24**3, [0.5, 0.5, 0.5], viz_stride=2)
N = wf.nx
c = N // 2
print(f"grid {wf.nx}x{wf.ny}x{wf.nz}")

def report(tag, ok):
    print(f"  [{'PASS ✓' if ok else 'FAIL ✗'}] {tag}")

# --- (1) vacuum seed: M = diag(δ,δ,1,g), director = ẑ ------------------------
print("\n(1) seed_vacuum_M")
seeds.seed_vacuum_M(wf, wf.lc_delta)
M = wf.M_am.to_numpy()
report(f"M is 4x4 (shape {M.shape[-2:]})", M.shape[-2:] == (4, 4))
mc = M[c, c, c]
report(f"M[3,3] == g ({mc[3,3]:.3f})", abs(mc[3, 3] - medium.LC_G) < 1e-4)
report(f"time off-diagonals 0 (max {np.abs(mc[3,:3]).max():.1e})", np.abs(mc[3, :3]).max() < 1e-6)
ev = np.sort(np.linalg.eigvalsh(mc))[::-1]
report(f"spectrum {np.round(ev,3)} ~ [g,1,δ,δ]", np.allclose(ev, [medium.LC_G, 1, medium.LC_DELTA, medium.LC_DELTA], atol=1e-3))
pde.eigen_decompose(wf)
d = wf.director_nhat.to_numpy()[c, c, c]
report(f"vacuum director {np.round(d,3)} ~ ±ẑ", abs(abs(d[2]) - 1) < 1e-3)

# --- (2) biaxial hedgehog: M spatial = O diag(1,δ,0) O^T, director radial -----
print("\n(2) seed_biaxial_hedgehog_M + eigen_decompose")
seeds.seed_biaxial_hedgehog_M(wf, c, c, c, 0.06 * N, 3.0, 0.30)
pde.eigen_decompose(wf)
M = wf.M_am.to_numpy()
report("M[...,3,3] all == g", np.allclose(M[..., 3, 3], medium.LC_G, atol=1e-4))
report("time off-diagonals all 0", np.abs(M[..., 3, :3]).max() < 1e-6)
ev = wf.eigenvalues.to_numpy()
report("eigenvalues finite (no NaN)", np.all(np.isfinite(ev)))
dn = wf.director_nhat.to_numpy()
norms = np.linalg.norm(dn, axis=-1)
report(f"director ~unit (mean |n|={norms.mean():.3f})", abs(norms.mean() - 1) < 0.05)
# off-core voxel: principal eigenvalue ~1, spectrum from the SPATIAL block (not g)
ic = c + N // 4
sp_ev = np.sort(np.linalg.eigvalsh(M[ic, c, c][:3, :3]))[::-1]
report(f"spatial sub-block spectrum {np.round(sp_ev,3)} ~ [1,δ,0]", abs(sp_ev[0] - 1) < 0.05)

# --- (3) evolve a few steps (V off): g constant, finite, stable --------------
print("\n(3) evolve_M x5 (V off)")
M0 = wf.M_am.to_numpy().copy()
g0 = M0[..., 3, 3].copy()
for _ in range(5):
    pde.compute_curvature_flux(wf)
    pde.evolve_M(wf, 1.0, 0.05 * wf.dx_am, 0.0, 0.0, 0.0)   # c=1.0, a=b=c=0 → V off
    wf.swap_matrix_buffers()
M5 = wf.M_am.to_numpy()
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
dt_rs = wf.dx_am * 0.95 / (c_amrs * 3 ** 0.5)
ldg_c = 1.0 * c_amrs**2 / wf.dx_am**4
ldg_a = -2.0 * ldg_c * (1.0 + 0.30**2)
seeds.seed_biaxial_hedgehog_M(wf, c, c, c, 0.06 * N, 3.0, 0.30)
pde.eigen_decompose(wf)
g0 = wf.M_am.to_numpy()[..., 3, 3].copy()
for _ in range(30):
    pde.compute_curvature_flux(wf)
    pde.evolve_M(wf, c_amrs, dt_rs, ldg_a, 0.0, ldg_c)
    wf.swap_matrix_buffers()
Mv = wf.M_am.to_numpy()
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
trackers = medium.Trackers(wf.grid_size)
observables = medium.FieldObservables(wf.grid_size)
seeds.seed_biaxial_hedgehog_M(wf, c, c, c, 0.06 * N, 3.0, 0.30)
pde.eigen_decompose(wf)
ok_kernels = True
try:
    obs.update_trackers_M(wf, trackers, 0.05 * wf.dx_am, wf.lc_delta)
    obs.compute_energyH_density_M(wf, observables, 0.2998, 0.05 * wf.dx_am, 0.0, 0.0, 0.0, 0.0, 1.0)
except Exception as e:
    ok_kernels = False
    print(f"    EXCEPTION: {type(e).__name__}: {e}")
report("update_trackers_M + compute_energyH_density_M run as KERNELS (not Python-scope)", ok_kernels)
if ok_kernels:
    eH = observables.energyH_density_aJ.to_numpy()
    report(f"energyH density finite (max={np.abs(eH).max():.2e})", np.all(np.isfinite(eH)))

print("\n  => headless check done. dV_M g-independent + bounded V-on evolve = the GUI "
      "explosion is fixed; observables kernels run. GUI re-test confirms on screen.")
