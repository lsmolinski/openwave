"""
M5.6.5b — biaxial-ellipsoid glyph geometry (headless de-risk before the GGUI kernel)

The "biaxial top at each voxel" showcase (4b Option c, deferred from M5.4 step 5 because
uniaxial M's two minor axes were degenerate). Now that M5.6 is genuinely biaxial
(D = diag(1, δ, 0), three distinct eigenvalues) the ellipsoid surface is meaningful.

KEY SIMPLIFICATION: a symmetric M maps the unit sphere to EXACTLY the ellipsoid whose
principal axes are M's eigenvectors with semi-axis lengths = M's eigenvalues:

    x = M · u ,  ‖u‖ = 1   ⇒   xᵀ M⁻² x = 1   (ellipsoid, semi-axes λ_i along e_i)

So the glyph needs NO eigendecomposition — per template sphere-vertex u, the ellipsoid
vertex is  p + scale · (M · u).  This script validates:
  1. icosphere template (vertices on the unit sphere + triangle faces, closed surface)
  2. M·u reproduces the correct ellipsoid (semi-axes = eigenvalues, oriented by eigenvectors)
  3. multi-glyph mesh assembly (per-voxel vertex block + offset face indices) — the kernel logic
  4. the λ_min visibility floor (D=(1,δ,0) ⇒ flat disk; a small floor keeps the glyph 3D-visible)

USAGE:
    python -m openwave.xperiments.m5_liquid_crystal.research.sandbox_v6.m5_6_5b_ellipsoid_geometry
"""
import numpy as np


def icosphere(subdiv=1):
    """Unit-sphere mesh by subdividing an icosahedron. Returns (verts[V,3], faces[F,3])."""
    t = (1.0 + np.sqrt(5.0)) / 2.0
    verts = np.array([
        [-1, t, 0], [1, t, 0], [-1, -t, 0], [1, -t, 0],
        [0, -1, t], [0, 1, t], [0, -1, -t], [0, 1, -t],
        [t, 0, -1], [t, 0, 1], [-t, 0, -1], [-t, 0, 1],
    ], float)
    faces = np.array([
        [0, 11, 5], [0, 5, 1], [0, 1, 7], [0, 7, 10], [0, 10, 11],
        [1, 5, 9], [5, 11, 4], [11, 10, 2], [10, 7, 6], [7, 1, 8],
        [3, 9, 4], [3, 4, 2], [3, 2, 6], [3, 6, 8], [3, 8, 9],
        [4, 9, 5], [2, 4, 11], [6, 2, 10], [8, 6, 7], [9, 8, 1],
    ], int)
    for _ in range(subdiv):
        mid = {}
        new_faces = []
        vlist = list(verts)

        def midpoint(a, b):
            key = (min(a, b), max(a, b))
            if key not in mid:
                vlist.append((vlist[a] + vlist[b]) / 2.0)
                mid[key] = len(vlist) - 1
            return mid[key]

        for (a, b, c) in faces:
            ab, bc, ca = midpoint(a, b), midpoint(b, c), midpoint(c, a)
            new_faces += [[a, ab, ca], [b, bc, ab], [c, ca, bc], [ab, bc, ca]]
        verts = np.array(vlist)
        faces = np.array(new_faces)
    verts = verts / np.linalg.norm(verts, axis=1, keepdims=True)   # project to unit sphere
    return verts, faces


def main():
    print("=" * 72)
    print("M5.6.5b — biaxial-ellipsoid glyph geometry")
    print("=" * 72)
    u, faces = icosphere(subdiv=1)
    V, F = len(u), len(faces)
    closed = np.isclose(np.linalg.norm(u, axis=1), 1.0).all() and (np.bincount(faces.ravel()).min() >= 5)
    print(f"[1] icosphere(subdiv=1): {V} verts, {F} faces; on unit sphere + closed = {closed}")

    # --- 2. M·u reproduces the ellipsoid: semi-axes = eigenvalues, along eigenvectors -----
    th = 0.7
    Oa = np.array([[np.cos(th), -np.sin(th), 0], [np.sin(th), np.cos(th), 0], [0, 0, 1.0]])
    Ob = np.array([[1.0, 0, 0], [0, np.cos(0.5), -np.sin(0.5)], [0, np.sin(0.5), np.cos(0.5)]])
    O = Oa @ Ob
    DELTA = 0.3
    D = np.diag([1.0, DELTA, 0.0])
    M = O @ D @ O.T
    evals, evecs = np.linalg.eigh(M)                               # ascending
    print("\n[2] ellipsoid x = M·u for M = O·diag(1,δ,0)·Oᵀ:")
    # exact: max_u|M u · e_i| = ‖M e_i‖ = λ_i. The discrete max undershoots (no vertex on the
    # pole) but converges to λ_i as the icosphere is refined ⇒ confirms the geometry is exact.
    print(f"    {'eigenvalue':>12} {'subdiv1':>10} {'subdiv3':>10} {'subdiv3 rel.err':>16}")
    u3, _ = icosphere(subdiv=3)
    conv = True
    for i in range(3):
        e = evecs[:, i]
        ext1 = np.abs((u @ M.T) @ e).max()
        ext3 = np.abs((u3 @ M.T) @ e).max()
        rel = abs(ext3 - evals[i]) / max(evals[i], 1e-9)
        conv &= (rel < 0.02) or (evals[i] < 1e-9)
        print(f"    {evals[i]:>12.4f} {ext1:>10.4f} {ext3:>10.4f} {rel:>16.1e}")
    axis_ok = conv
    print(f"    → discrete extent → λ_i as the sphere refines (semi-axes = eigenvalues): {axis_ok}")

    # --- 3. multi-glyph mesh assembly (the kernel's per-voxel block + index offset) -------
    centers = np.array([[0.2, 0.5, 0.5], [0.5, 0.5, 0.5], [0.8, 0.5, 0.5]])
    Ms = [O @ np.diag([1.0, 0.5, 0.5]) @ O.T, M, O @ np.diag([1.0, 0.6, 0.2]) @ O.T]  # uni, biax, biax
    scale = 0.04
    all_v, all_f = [], []
    for g, (p, Mg) in enumerate(zip(centers, Ms)):
        all_v.append(p[None, :] + scale * (u @ Mg.T))             # V verts for this glyph
        all_f.append(faces + g * V)                                # faces offset into the block
    all_v = np.concatenate(all_v)
    all_f = np.concatenate(all_f)
    idx_ok = all_f.max() == len(all_v) - 1 and all_f.min() == 0
    print(f"\n[3] multi-glyph assembly: {len(centers)} glyphs → {len(all_v)} verts, {len(all_f)} faces")
    print(f"    face-index range [{all_f.min()}, {all_f.max()}] matches vert count {len(all_v)}: {idx_ok}")

    # --- 4. λ_min visibility floor: D=(1,δ,0) is a flat disk; floor keeps it 3D ----------
    flat_thick = np.abs((u @ M.T) @ evecs[:, 0]).max()             # extent along the λ=0 axis
    FLOOR = 0.12
    Df = np.diag([1.0, max(DELTA, FLOOR), max(0.0, FLOOR)])
    Mf = O @ Df @ O.T
    floor_thick = np.abs((u @ Mf.T) @ evecs[:, 0]).max()
    print(f"\n[4] λ_min floor: raw min-axis extent = {flat_thick:.3f} (≈0, flat disk); "
          f"with floor {FLOOR} → {floor_thick:.3f} (3D-visible)")
    print(f"    (viz-only floor on the smallest eigenvalue; physics M unchanged)")

    ok = closed and axis_ok and idx_ok and floor_thick > 0.05
    print("\n" + "=" * 72)
    print("M5.6.5b geometry: M·u = ellipsoid (no eigendecomp needed); multi-glyph mesh assembly")
    print("validated. Production kernel writes, per sampled voxel v and template vertex k:")
    print("   vert[v·V+k] = voxel_center_norm + scale·(M_floored · u[k]) ;  face[v·F+f] = tmpl+v·V")
    print("PASS" if ok else "PARTIAL — inspect above")
    print("=" * 72)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
