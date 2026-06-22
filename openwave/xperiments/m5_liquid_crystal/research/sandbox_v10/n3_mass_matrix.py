#!/usr/bin/env python3
"""
N3 / S1 , the closed-loop FIELD THEORY -> flavour mass matrix (the machine the search drives).

The bridge (checkpoints/04_n3_design.md, 05_n3_derisk.md): the three neutrino FLAVOUR states are the
SAME closed disclination loop at three SO(3) ORIENTATIONS (issue #199: oscillation = an SO(3) rotation
among flavours). Co-located (same center, three orientations) is the cleanest realization , the coupling
then depends ONLY on the relative rotation R_a^T R_b, a pure orientation effect.

The flavour mass matrix is the ENERGY HESSIAN of the LdG functional projected onto the three flavour
field displacements (the standard "mass = second variation of the energy"):

    dM_a(x) = M_a(x) - M_vac(x)                         (flavour-a displacement from the common vacuum)
    K_ab = INT < grad dM_a , grad dM_b >_s  dx          (kinetic / curvature hopping, the structure carrier)
    P_ab = INT < dM_a , dM_b >_s  dx                    (potential / overlap, the on-site Gram)
    M_mass = K + kappa P     (real-symmetric 3x3; kappa weights the on-site potential curvature)
    U = eigenvectors of M_mass   ->   pmns_angles(U)    (n3_derisk canonical assignment)

mu-tau (2<->3) symmetry is GEOMETRIC: e-loop = reference (R_e=I), mu/tau = a MIRROR PAIR
(R_mu=Rx(+alpha), R_tau=Rx(-alpha)). Reflection is an LdG symmetry mapping mu<->tau and fixing e, so by
construction M_e,mu=M_e,tau and M_mu,mu=M_tau,tau -> the matrix lands in the de-risk [[x,y,y],[y,z,w],
[y,w,z]] form. The MAGIC condition (x+y=z+w) is one scalar the geometry must hit -> the S2 scan crosses it
(the TBM gate). delta (index-2 SO(3)-breaking eigenvalue) sources theta13 (S3).

Convention: index-0 (Duda), D=diag(g,1,delta,0), eta=diag(-1,1,1,1). Far-field blends to a COMMON zhat
so dM_a is localized at the core (the overlap is core-core, orientation-dependent, not a trivial far-field
offset). Undressed (g decoupled in [0,0], cancels in dM) -> no 1e20 range here; N1 precision enters only
with boost dressing (S3). Headless numpy f64. LOCAL (#236 N-program, HELD). Run: python3 n3_mass_matrix.py
"""

import json
import os
import numpy as np

from n1_precision_method import SIGN_MAT, _grads, _comm, _sdot
from n2_closed_loop import G_SCALE
from n3_derisk import angles_from_M, tbm_err, TH12_TBM

HERE = os.path.dirname(os.path.abspath(__file__))
DELTA_PHYS = 1.0e-10


# ----------------------------------------------------------------------------------
# rotations
# ----------------------------------------------------------------------------------
def rot_axis(axis, ang):
    """3x3 rotation by `ang` about a unit `axis` (Rodrigues)."""
    a = np.asarray(axis, float); a = a / (np.linalg.norm(a) + 1e-15)
    c, s = np.cos(ang), np.sin(ang)
    K = np.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])
    return np.eye(3) * c + s * K + (1 - c) * np.outer(a, a)


# ----------------------------------------------------------------------------------
# oriented closed-loop seeder (common far-field zhat; only the core winding rotates)
# ----------------------------------------------------------------------------------
def seed_loop_oriented(n, R3, R_loop, delta, q=0.5, core_vox=2.0, blend_pow=4):
    """4x4 LdG field M(x) for a circular disclination loop whose plane + core winding are rotated by
    the 3x3 rotation R3, with the far field blended to the GLOBAL zhat (common to all flavours).

    Local construction (N2 seed_loop_M, in the R3^T-rotated frame): the loop lies in the local z=0
    plane radius R_loop; the director winds by q*psi (psi = meridian angle) around the core. We rotate
    the winding director by R3 and blend it toward the global zhat with a core-localized weight, so the
    displacement dM = M - M_vac is supported near the (rotated) core."""
    c0 = (n - 1) / 2.0
    idx = np.arange(n) - c0
    X, Y, Z = np.meshgrid(idx, idx, idx, indexing="ij")
    P = np.stack([X, Y, Z], axis=-1)                       # (n,n,n,3) global coords
    Ploc = np.einsum("ab,...b->...a", R3.T, P)             # local coords = R3^T x
    xl, yl, zl = Ploc[..., 0], Ploc[..., 1], Ploc[..., 2]

    s = np.arctan2(yl, xl)                                  # azimuth around loop axis (local)
    er = np.stack([np.cos(s), np.sin(s), np.zeros_like(s)], axis=-1)   # local radial
    ez = np.zeros_like(er); ez[..., 2] = 1.0                          # local loop-axis
    rho_r = np.sqrt(xl * xl + yl * yl) - R_loop
    rho_z = zl
    psi = np.arctan2(rho_z, rho_r)
    d_core = np.sqrt(rho_r * rho_r + rho_z * rho_z)
    ang = q * psi
    d_wind_loc = np.cos(ang)[..., None] * er + np.sin(ang)[..., None] * ez   # winding director (local)
    d_wind = np.einsum("ab,...b->...a", R3, d_wind_loc)    # rotate to global frame

    w = 1.0 / (1.0 + (d_core / (3.0 * core_vox))[..., None] ** blend_pow)    # core weight
    zhat = np.zeros_like(d_wind); zhat[..., 2] = 1.0        # COMMON global far-field
    n_blend = w * d_wind + (1.0 - w) * zhat
    n_unit = n_blend / (np.linalg.norm(n_blend, axis=-1, keepdims=True) + 1e-12)

    eye3 = np.eye(3)
    nn = np.einsum("...a,...b->...ab", n_unit, n_unit)
    M_sp = delta * eye3 + (1.0 - delta) * nn               # uniaxial, eigenvalues (1, delta, delta)
    M = np.zeros((n, n, n, 4, 4))
    M[..., 0, 0] = G_SCALE
    M[..., 1:, 1:] = M_sp
    return M


def vacuum_field(n, delta):
    """Common vacuum: g at [0,0], spatial uniaxial along the global zhat (the far field of every loop)."""
    M = np.zeros((n, n, n, 4, 4))
    M[..., 0, 0] = G_SCALE
    zz = np.zeros(3); zz[2] = 1.0
    M_sp = delta * np.eye(3) + (1.0 - delta) * np.outer(zz, zz)
    M[..., 1:, 1:] = M_sp
    return M


# ----------------------------------------------------------------------------------
# the energy-Hessian coupling matrix (kinetic hopping + potential overlap)
# ----------------------------------------------------------------------------------
def _grad_overlap(dA, dB):
    """INT <grad dA, grad dB>_s over the interior (the kinetic/curvature hopping integral)."""
    Ax, Ay, Az = _grads(dA)
    Bx, By, Bz = _grads(dB)
    return float(np.sum(_sdot(Ax, Bx) + _sdot(Ay, By) + _sdot(Az, Bz)))


def _field_overlap(dA, dB):
    """INT <dA, dB>_s over the interior (the on-site potential / Gram overlap)."""
    return float(np.sum(_sdot(dA[1:-1, 1:-1, 1:-1], dB[1:-1, 1:-1, 1:-1])))


def coupling_matrix(dfields, kappa=0.0):
    """Real-symmetric 3x3 mass matrix M = K + kappa P from the three flavour displacements dfields."""
    K = np.zeros((3, 3)); P = np.zeros((3, 3))
    for a in range(3):
        for b in range(a, 3):
            K[a, b] = K[b, a] = _grad_overlap(dfields[a], dfields[b])
            P[a, b] = P[b, a] = _field_overlap(dfields[a], dfields[b])
    return K + kappa * P, K, P


# ----------------------------------------------------------------------------------
# build the three flavour displacements for a given geometry, return the mass matrix + angles
# ----------------------------------------------------------------------------------
def flavour_mass_matrix(n=40, alpha=0.6, R_loop=9.0, delta=DELTA_PHYS, q=0.5,
                        core_vox=2.0, kappa=0.0, tilt_axis=(1.0, 0.0, 0.0)):
    """e-loop = reference (R=I); mu/tau = mirror pair Rx(+alpha)/Rx(-alpha) about tilt_axis. Returns
    (M_mass 3x3, angles, eigenvalues, structure-diagnostics)."""
    Re = np.eye(3)
    Rmu = rot_axis(tilt_axis, +alpha)
    Rtau = rot_axis(tilt_axis, -alpha)
    Mvac = vacuum_field(n, delta)
    fields = [seed_loop_oriented(n, R, R_loop, delta, q=q, core_vox=core_vox)
              for R in (Re, Rmu, Rtau)]
    dfields = [f - Mvac for f in fields]
    M_mass, K, P = coupling_matrix(dfields, kappa=kappa)
    U, ang, w = angles_from_M(M_mass)
    # structure diagnostics vs the de-risk [[x,y,y],[y,z,w],[y,w,z]] form
    mutau_viol = float(abs(M_mass[0, 1] - M_mass[0, 2]) + abs(M_mass[1, 1] - M_mass[2, 2]))
    x, y, z, ww = M_mass[0, 0], M_mass[0, 1], M_mass[1, 1], M_mass[1, 2]
    magic_resid = float((x + y) - (z + ww))                # = 0 at the TBM gate
    return {
        "M_mass": M_mass.tolist(), "angles": ang, "eigenvalues": w.tolist(),
        "mutau_violation": mutau_viol, "magic_residual": magic_resid,
        "tbm_err": float(tbm_err(ang)),
        "params": {"n": n, "alpha": alpha, "R_loop": R_loop, "delta": delta, "q": q,
                   "core_vox": core_vox, "kappa": kappa},
    }


def main():
    print("=" * 80)
    print("N3 / S1 , closed-loop field theory -> flavour mass matrix (machine check)")
    print("=" * 80)

    # a representative geometry; S2 scans alpha (+ R, core, kappa) for the magic crossing
    for alpha in (0.3, 0.6, 0.9, 1.2):
        r = flavour_mass_matrix(n=40, alpha=alpha, R_loop=9.0, core_vox=2.0, kappa=0.0)
        M = np.array(r["M_mass"])
        a = r["angles"]
        print(f"\n[alpha={alpha:.2f}]  mu-tau viol={r['mutau_violation']:.2e}  "
              f"magic resid (x+y)-(z+w)={r['magic_residual']:+.4e}")
        print(f"            M_mass = [[{M[0,0]:+.3e},{M[0,1]:+.3e},{M[0,2]:+.3e}],")
        print(f"                      [{M[1,0]:+.3e},{M[1,1]:+.3e},{M[1,2]:+.3e}],")
        print(f"                      [{M[2,0]:+.3e},{M[2,1]:+.3e},{M[2,2]:+.3e}]]")
        print(f"            angles: th12={a['theta12']:.3f} th23={a['theta23']:.3f} "
              f"th13={a['theta13']:.3f}  (TBM err {r['tbm_err']:.3f} deg)")

    # quick check: does the mu-tau form hold by construction (it must, by the mirror symmetry)?
    r0 = flavour_mass_matrix(n=40, alpha=0.6)
    mutau_ok = r0["mutau_violation"] < 1e-8 * abs(np.array(r0["M_mass"])).max()
    print("\n" + "=" * 80)
    print(f"S1 machine: mu-tau form holds by construction = {mutau_ok}  "
          f"(viol {r0['mutau_violation']:.2e}); magic residual is the knob S2 drives to 0.")
    print("  next: n3_search.py scans the geometry for magic_residual -> 0 = the TBM gate.")
    print("=" * 80)

    summary = {
        "construction": "3 co-located loops, e=ref, mu/tau=Rx(+-alpha) mirror pair; "
                        "M_mass = grad-overlap + kappa field-overlap (energy Hessian projection)",
        "convention": "index-0 Duda; far field common zhat; undressed (N1 precision enters with dressing)",
        "mutau_form_by_construction": bool(mutau_ok),
        "sample_alpha_0p6": r0,
    }
    with open(os.path.join(HERE, "n3_mass_matrix_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    print(f"summary -> {os.path.join(HERE, 'n3_mass_matrix_summary.json')}")
    return mutau_ok


if __name__ == "__main__":
    raise SystemExit(0 if main() else 1)
