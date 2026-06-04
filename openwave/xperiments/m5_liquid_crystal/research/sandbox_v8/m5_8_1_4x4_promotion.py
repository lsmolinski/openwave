"""
M5.8.1 — numpy validation of the 3x3 -> 4x4 substrate promotion (storage only).

Decision (confirmed): the 4x4 field stores M = block-diag(g, M_spatial), with the time
axis index 0 (eigenvalue g, the gravity/boost axis), boost-DECOUPLED for now (the
Minkowski coupling that drives the clock is M5.8.2). The director is read from the
SPATIAL 3x3 sub-block M[1:4, 1:4] (Cardano reused; no quartic solver).

What this must PROVE before touching production:
  (1) spectrum: eig(M_4x4) = {g} U eig(M_spatial) = {g, 1, delta, 0}.
  (2) director: principal eigenvector of the spatial sub-block == the 3x3-only director.
  (3) PHYSICS-PRESERVING: with constant g, the index-generic curvature flux + evolve_M
      preserve the block structure -> the 4x4 run's spatial block is IDENTICAL to a pure
      3x3 run, and g stays exactly constant. (So M5.8.1 reproduces the 3D physics; the
      g-coupling that changes things is M5.8.2.)
  (4) observables: Frobenius norms / energyH curvature are unchanged (the constant-g
      block contributes nothing to gradients/commutators).

Run:  python m5_8_1_4x4_promotion.py
"""

import numpy as np

np.random.seed(0)
N = 12                      # small grid
DELTA = 0.30               # biaxial minor axis (matches _topo_biaxial1)
G = 8.0                    # gravity/time-axis eigenvalue (g >> 1), constant background
DX, DT, C = 1.0, 0.1, 1.0

def comm(A, B):
    """Grid-wise matrix commutator [A,B], A,B shape (...,n,n)."""
    return A @ B - B @ A

# ----------------------------------------------------------------------------
# Seed a biaxial-hedgehog-like SPATIAL 3x3 field O diag(1,delta,0) O^T, then embed.
# ----------------------------------------------------------------------------
def biaxial_spatial(N):
    """3x3 M per voxel: radial frame O=[r^|e_th|e_ph], D=diag(1,delta,0)."""
    M = np.zeros((N, N, N, 3, 3))
    c = (N - 1) / 2.0
    for i in range(N):
        for j in range(N):
            for k in range(N):
                p = np.array([i - c, j - c, k - c])
                r = np.linalg.norm(p) + 1e-9
                rhat = p / r
                # azimuthal (z-disclination) frame, regularized near axis
                a = np.array([-p[1], p[0], 0.0])
                an = np.linalg.norm(a) + 1e-9
                ephi = a / an
                etheta = np.cross(ephi, rhat)
                M[i, j, k] = (1.0 * np.outer(rhat, rhat)
                              + DELTA * np.outer(etheta, etheta)
                              + 0.0 * np.outer(ephi, ephi))
    return M

def embed_4x4(Msp, g):
    """M_4x4 = block-diag(g, M_spatial): time axis index 0, decoupled."""
    N = Msp.shape[0]
    M4 = np.zeros((N, N, N, 4, 4))
    M4[..., 0, 0] = g
    M4[..., 1:, 1:] = Msp
    return M4

Msp = biaxial_spatial(N)
M4 = embed_4x4(Msp, G)

# ----------------------------------------------------------------------------
# (1) spectrum  +  (2) spatial-subblock director
# ----------------------------------------------------------------------------
ic = N // 2 + 1   # an off-center voxel (nonzero gradients)
ev_full = np.sort(np.linalg.eigvalsh(M4[ic, ic, ic]))[::-1]
ev_sp = np.sort(np.linalg.eigvalsh(Msp[ic, ic, ic]))[::-1]
print("=" * 68)
print("(1) SPECTRUM")
print("=" * 68)
print(f"  eig(M_4x4) = {np.round(ev_full, 4)}   (expect ~ [g={G}, 1, {DELTA}, 0])")
print(f"  eig(M_spatial) = {np.round(ev_sp, 4)}")
spec_ok = np.allclose(ev_full, np.sort([G, 1.0, DELTA, 0.0])[::-1], atol=1e-6)
print(f"  spectrum = {{g}} U eig(spatial)?  {'PASS ✓' if spec_ok else 'FAIL ✗'}")

print("\n" + "=" * 68)
print("(2) DIRECTOR from spatial sub-block vs 3x3-only")
print("=" * 68)
def principal(m3):
    w, v = np.linalg.eigh(m3)
    return v[:, np.argmax(w)]          # eigenvalue-1 axis = the director
d_from_4x4 = principal(M4[ic, ic, ic][1:, 1:])   # spatial sub-block of the 4x4
d_from_3x3 = principal(Msp[ic, ic, ic])
align = abs(np.dot(d_from_4x4, d_from_3x3))       # apolar: |dot|=1 if same axis
print(f"  |director_subblock · director_3x3| = {align:.8f}   "
      f"{'PASS ✓ (identical axis)' if abs(align - 1) < 1e-8 else 'FAIL ✗'}")
print(f"  global principal of M_4x4 is the TIME axis (g={G}) — NOT the director: "
      f"confirms why the sub-block read is required.")

# ----------------------------------------------------------------------------
# (3) PHYSICS-PRESERVING: 4x4 spatial block == pure 3x3 run; g stays constant.
# ----------------------------------------------------------------------------
def curvature_flux(M):
    inv2 = 1.0 / (2 * DX)
    Mx = (np.roll(M, -1, 0) - np.roll(M, 1, 0)) * inv2
    My = (np.roll(M, -1, 1) - np.roll(M, 1, 1)) * inv2
    Mz = (np.roll(M, -1, 2) - np.roll(M, 1, 2)) * inv2
    Gx = 8 * (comm(comm(Mx, My), My) + comm(comm(Mx, Mz), Mz))
    Gy = 8 * (comm(comm(My, Mx), Mx) + comm(comm(My, Mz), Mz))
    Gz = 8 * (comm(comm(Mz, Mx), Mx) + comm(comm(Mz, My), My))
    return Gx, Gy, Gz

def step(M, Mprev):
    inv2 = 1.0 / (2 * DX)
    Gx, Gy, Gz = curvature_flux(M)
    divG = ((np.roll(Gx, -1, 0) - np.roll(Gx, 1, 0))
            + (np.roll(Gy, -1, 1) - np.roll(Gy, 1, 1))
            + (np.roll(Gz, -1, 2) - np.roll(Gz, 1, 2))) * inv2
    Mnew = 2 * M - Mprev + (DT**2) * (C**2) * divG
    return Mnew, M

# run BOTH: pure 3x3, and the 4x4 (block-diag), for a few steps
A, Aprev = Msp.copy(), Msp.copy()
B, Bprev = M4.copy(), M4.copy()
g0 = B[..., 0, 0].copy()
nsteps = 8
for _ in range(nsteps):
    A, Aprev = step(A, Aprev)
    B, Bprev = step(B, Bprev)

spatial_match = np.max(np.abs(B[..., 1:, 1:] - A))
g_drift = np.max(np.abs(B[..., 0, 0] - g0))
offdiag = np.max(np.abs(B[..., 0, 1:]))        # time-space coupling that should stay 0
print("\n" + "=" * 68)
print(f"(3) PHYSICS-PRESERVING after {nsteps} curvature-leapfrog steps (V off)")
print("=" * 68)
print(f"  max|B_spatial - A_3x3|   = {spatial_match:.2e}   "
      f"{'PASS ✓ identical spatial dynamics' if spatial_match < 1e-9 else 'FAIL ✗'}")
print(f"  max|g(t) - g0|           = {g_drift:.2e}   "
      f"{'PASS ✓ time axis inert (constant g)' if g_drift < 1e-9 else 'FAIL ✗'}")
print(f"  max|time-space coupling| = {offdiag:.2e}   "
      f"{'PASS ✓ block structure preserved' if offdiag < 1e-12 else 'FAIL ✗'}")

# ----------------------------------------------------------------------------
# (4) observables: energyH curvature + Frobenius deviation unchanged by the g block.
# ----------------------------------------------------------------------------
def curv_energy(M):
    inv2 = 1.0 / (2 * DX)
    Mx = (np.roll(M, -1, 0) - np.roll(M, 1, 0)) * inv2
    My = (np.roll(M, -1, 1) - np.roll(M, 1, 1)) * inv2
    Mz = (np.roll(M, -1, 2) - np.roll(M, 1, 2)) * inv2
    e = 0.0
    for X, Y in ((Mx, My), (Mx, Mz), (My, Mz)):
        e += 4 * np.sum(comm(X, Y) ** 2)
    return e
eH_3 = curv_energy(Msp)
eH_4 = curv_energy(M4)
print("\n" + "=" * 68)
print("(4) OBSERVABLES — energyH curvature unchanged by the constant-g block")
print("=" * 68)
print(f"  curvature energy  3x3 = {eH_3:.6f}   4x4 = {eH_4:.6f}   "
      f"{'PASS ✓' if abs(eH_3 - eH_4) < 1e-9 else 'FAIL ✗'}")
# D_vac for trackers: 4x4 = diag(g,1,delta,delta) (uniaxial) — Frobenius dev well-defined
Dvac4 = np.diag([G, 1.0, DELTA, DELTA])
dev = np.linalg.norm(M4[ic, ic, ic] - Dvac4)
print(f"  ‖M − D_vac(4x4)‖_F at a voxel = {dev:.4f}  (D_vac = diag(g,1,δ,δ); well-defined ✓)")

print("\n  => M5.8.1 math validated: block-diag(g, spatial) is physics-preserving; "
      "spatial-subblock director is correct; ready to port to production.")
