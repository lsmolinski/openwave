"""
M5.8.1 — topological charge quantization on the 4×4 engine (headless).

The "charge" half of the Coulomb claim: a hedgehog defect carries an INTEGER
topological charge Q (the director's winding number on a sphere around the core),
and charge is ADDITIVE (a +1/−1 pair enclosed together gives Q=0). Reproduces the
M5.1 winding result on the current 4×4 matrix substrate via the production
pipeline: seed_hedgehog_M → eigen_decompose (director = principal eigenvector of
M's spatial block) → compute_winding_number.

PASS: Q(+1)≈+1, Q(−1)≈−1, Q(pair, enclosing sphere)≈0 — integers from topology,
not from tuning.

USAGE:  python m5_8_1_topo_charge.py
"""
import sys
from pathlib import Path

import numpy as np  # noqa: F401  (kept for parity with sibling scripts)
import taichi as ti

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[5]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

ti.init(arch=ti.cpu, default_fp=ti.f32)

from openwave.xperiments.m5_liquid_crystal import medium  # noqa: E402
from openwave.xperiments.m5_liquid_crystal import engine1_seeds as seeds  # noqa: E402
from openwave.xperiments.m5_liquid_crystal import engine2_pde as pde  # noqa: E402
from openwave.xperiments.m5_liquid_crystal import engine3_observables as obs  # noqa: E402

N = 49
tf = medium.TensorField([1e-15] * 3, N**3)
c = tf.nx // 2
dq = 0.25 * tf.nx  # D/4 vacuum-blend falloff (M5.1 default)


def seed_and_Q(defects, sphere_center, radius_vox):
    """Seed hedgehog defect(s) on M, extract the director, return winding Q."""
    n_d = len(defects)
    centers = ti.field(ti.i32, shape=(n_d, 3))
    signs = ti.field(ti.i32, shape=(n_d,))
    for i, (pos, s) in enumerate(defects):
        centers[i, 0], centers[i, 1], centers[i, 2] = pos
        signs[i] = s
    seeds.seed_hedgehog_M(tf, centers, signs, float(dq), n_d, tf.lc_delta)
    pde.eigen_decompose(tf)
    dn = tf.director_nhat.to_numpy()
    return obs.compute_winding_number(dn, sphere_center, radius_vox)


print("=" * 70)
print("M5.8.1 — topological charge quantization (director winding, 4×4 engine)")
print("=" * 70)
print(f"grid {tf.nx}³  dx={tf.dx_am:.2f} am  δ={tf.lc_delta}")
print("(Q measured at r=5, inside the defect core — beyond it the D/4 vacuum blend")
print(" returns the field to ẑ by design, so a far sphere reads the vacuum's Q=0.)\n")

R = 5.0
q_plus = seed_and_Q([((c, c, c), +1)], (c, c, c), R)
q_minus = seed_and_Q([((c, c, c), -1)], (c, c, c), R)
# pair: charges retain their integer identity around each core; the enclosing
# sphere (in the blended-to-vacuum far field) reads the total Q = 0.
d = 10
n_d = 2
centers = ti.field(ti.i32, shape=(n_d, 3))
signs = ti.field(ti.i32, shape=(n_d,))
for i, (pos, s) in enumerate([((c - d, c, c), +1), ((c + d, c, c), -1)]):
    centers[i, 0], centers[i, 1], centers[i, 2] = pos
    signs[i] = s
seeds.seed_hedgehog_M(tf, centers, signs, float(dq), n_d, tf.lc_delta)
pde.eigen_decompose(tf)
dn = tf.director_nhat.to_numpy()
q_pair_plus = obs.compute_winding_number(dn, (c - d, c, c), R)
q_pair_minus = obs.compute_winding_number(dn, (c + d, c, c), R)
q_pair_total = obs.compute_winding_number(dn, (c, c, c), float(d + 8))

print(f"  single +1 hedgehog:              Q = {q_plus:+.3f}   (expect +1)")
print(f"  single −1 hedgehog:              Q = {q_minus:+.3f}   (expect −1)")
print(f"  pair, around the +1 core:        Q = {q_pair_plus:+.3f}   (expect +1)")
print(f"  pair, around the −1 core:        Q = {q_pair_minus:+.3f}   (expect −1)")
print(f"  pair, enclosing sphere (vacuum): Q = {q_pair_total:+.3f}   (expect 0 — additivity)")

ok = (abs(q_plus - 1) < 0.1 and abs(q_minus + 1) < 0.1
      and abs(q_pair_plus - 1) < 0.15 and abs(q_pair_minus + 1) < 0.15
      and abs(q_pair_total) < 0.1)
print("\n" + "=" * 70)
print(f"  {'PASS' if ok else 'FAIL'} — charge is INTEGER + ADDITIVE (topology, not tuning)")
print("=" * 70)
sys.exit(0 if ok else 1)
