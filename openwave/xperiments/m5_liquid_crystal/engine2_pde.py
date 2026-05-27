"""
M5 ENGINE — PDE (engine2_pde)

The complete PDE layer: the reusable differential operators AND the equations
integrated with them — merged from the former operators + PDE modules
(2026-05-25 SoC consolidation). Two internal sections below:
  DIFFERENTIAL OPERATORS — @ti.func ∇², ∇·, ∇×, ∇×∇× primitives
  PDE DYNAMICS           — evolve_psi (leapfrog), V_psi (potential hook),
                           relax_director_step (dissipative gradient descent)

Concern order: engine1_seeds → engine2_pde → engine3_observables → engine4_render.
M5.4 watch-item: if eigen_decompose ends up shared by render + observables + the
PDE, re-extract a small operators module at that point.
"""

import taichi as ti

# ================================================================
# DIFFERENTIAL OPERATORS
# ================================================================
# Vector calculus on the Vector(3) field ψ stored in psi_am. All operators
# are 2nd-order central-difference @ti.func kernels — caller is responsible
# for skipping boundary voxels (each operator's halo requirement is noted in
# its docstring).
#
#   compute_laplacian   — ∇²ψ          → Vector(3); 1-cell halo
#   compute_divergence  — ∇·ψ          → scalar;    1-cell halo
#   compute_curl        — ∇×ψ          → Vector(3); 1-cell halo
#   compute_curl_curl   — ∇×(∇×ψ)      → Vector(3); 2-cell halo (uses identity)
#
# The curl-curl is implemented via the vector identity
#       ∇×(∇×ψ) = ∇(∇·ψ) − ∇²ψ
# rather than applying compute_curl twice. Two reasons:
#   1. Reuses validated Laplacian — only the gradient-of-divergence part is
#      new code, easier to verify
#   2. 2-cell halo (vs. 4-cell for nested curl) keeps the loop range tight
#      and matches Exp 7 v2's implementation (the reference)
#
# Why @ti.func not @ti.kernel: these are inner-loop primitives meant to be
# called from kernels that iterate over voxels (evolve_psi, future force
# kernels, divergence-cleaning projections, etc.), not standalone passes.


@ti.func
def compute_laplacian(
    wave_field: ti.template(),  # type: ignore
    i: ti.i32,  # type: ignore
    j: ti.i32,  # type: ignore
    k: ti.i32,  # type: ignore
):
    """
    Compute the vector Laplacian ∇²ψ at voxel [i, j, k] using a 6-point stencil.

    Reads from `wave_field.psi_am` (the current-time buffer, ψ at t).
    Returns a Vector(3) — the discrete Laplacian acts component-wise:
        ∇²ψ = (∇²ψ_x, ∇²ψ_y, ∇²ψ_z)

    Standard 2nd-order finite-difference stencil:
        ∇²ψ ≈ (face_sum − 6·center) / dx²
    where face_sum = ψ[i±1] + ψ[j±1] + ψ[k±1] (six face neighbors).

    Note: the L/T (longitudinal vs transverse) decomposition of the field is a
    property of the field configuration (radial vs tangential displacement),
    NOT of the operator. M5's vector field handles both together; the L/T view
    is computed post hoc by projecting ψ onto radial/tangential directions.

    Args:
        wave_field: WaveField instance (reads psi_am)
        i, j, k: Voxel indices (must be interior: 0 < i < nx-1, etc. — caller
            handles boundary skip)

    Returns:
        ∇²ψ as Vector(3) in units [am / am²] = [1/am]
    """
    face_sum = (
        wave_field.psi_am[i + 1, j, k]
        + wave_field.psi_am[i - 1, j, k]
        + wave_field.psi_am[i, j + 1, k]
        + wave_field.psi_am[i, j - 1, k]
        + wave_field.psi_am[i, j, k + 1]
        + wave_field.psi_am[i, j, k - 1]
    )
    center = wave_field.psi_am[i, j, k]
    laplacian_psi_am = (face_sum - 6.0 * center) / (wave_field.dx_am**2)
    return laplacian_psi_am


@ti.func
def compute_divergence(
    wave_field: ti.template(),  # type: ignore
    i: ti.i32,  # type: ignore
    j: ti.i32,  # type: ignore
    k: ti.i32,  # type: ignore
):
    """
    Compute the scalar divergence ∇·ψ at voxel [i, j, k].

    Standard 2nd-order central-difference stencil:
        ∇·ψ ≈ (∂_x ψ_x + ∂_y ψ_y + ∂_z ψ_z)
            = [ ψ_x[i+1] − ψ_x[i−1]
              + ψ_y[j+1] − ψ_y[j−1]
              + ψ_z[k+1] − ψ_z[k−1] ] / (2·dx)

    Used by:
      - M5.1 winding-number tracker (∇·n on the director field is an
        observable near hedgehog defects)
      - M5.2 Close Eq. 23 — the constraint ∇·s = 0 must be enforced;
        divergence-cleaning projection reads this each step
      - M5.0e curl_curl identity ∇×(∇×ψ) = ∇(∇·ψ) − ∇²ψ

    Args:
        wave_field: WaveField instance (reads psi_am)
        i, j, k: Voxel indices (must be interior: 0 < i < nx-1, etc. —
            caller handles boundary skip; 1-cell halo)

    Returns:
        ∇·ψ as ti.f32 scalar in units [am / am] = [dimensionless]
    """
    inv_2dx = 1.0 / (2.0 * wave_field.dx_am)
    dpsi_x_dx = (wave_field.psi_am[i + 1, j, k][0] - wave_field.psi_am[i - 1, j, k][0]) * inv_2dx
    dpsi_y_dy = (wave_field.psi_am[i, j + 1, k][1] - wave_field.psi_am[i, j - 1, k][1]) * inv_2dx
    dpsi_z_dz = (wave_field.psi_am[i, j, k + 1][2] - wave_field.psi_am[i, j, k - 1][2]) * inv_2dx
    return dpsi_x_dx + dpsi_y_dy + dpsi_z_dz


@ti.func
def compute_curl(
    wave_field: ti.template(),  # type: ignore
    i: ti.i32,  # type: ignore
    j: ti.i32,  # type: ignore
    k: ti.i32,  # type: ignore
):
    """
    Compute the vector curl ∇×ψ at voxel [i, j, k].

    Standard 2nd-order central-difference stencil:
        (∇×ψ)_x = ∂_y ψ_z − ∂_z ψ_y
        (∇×ψ)_y = ∂_z ψ_x − ∂_x ψ_z
        (∇×ψ)_z = ∂_x ψ_y − ∂_y ψ_x

    Used by:
      - M5.2 Close Eq. 19 linear limit  ∂²_t Q = −c² · ∇×(∇×Q)
      - Spin-density observables (s = curl of velocity field analog)
      - Magnetic-channel diagnostics in M7+ (transverse component extraction)

    Args:
        wave_field: WaveField instance (reads psi_am)
        i, j, k: Voxel indices (must be interior: 0 < i < nx-1, etc. —
            caller handles boundary skip; 1-cell halo)

    Returns:
        ∇×ψ as Vector(3) in units [am / am] = [dimensionless]
    """
    inv_2dx = 1.0 / (2.0 * wave_field.dx_am)
    # ∂_y ψ_z, ∂_z ψ_y
    dpsi_z_dy = (wave_field.psi_am[i, j + 1, k][2] - wave_field.psi_am[i, j - 1, k][2]) * inv_2dx
    dpsi_y_dz = (wave_field.psi_am[i, j, k + 1][1] - wave_field.psi_am[i, j, k - 1][1]) * inv_2dx
    # ∂_z ψ_x, ∂_x ψ_z
    dpsi_x_dz = (wave_field.psi_am[i, j, k + 1][0] - wave_field.psi_am[i, j, k - 1][0]) * inv_2dx
    dpsi_z_dx = (wave_field.psi_am[i + 1, j, k][2] - wave_field.psi_am[i - 1, j, k][2]) * inv_2dx
    # ∂_x ψ_y, ∂_y ψ_x
    dpsi_y_dx = (wave_field.psi_am[i + 1, j, k][1] - wave_field.psi_am[i - 1, j, k][1]) * inv_2dx
    dpsi_x_dy = (wave_field.psi_am[i, j + 1, k][0] - wave_field.psi_am[i, j - 1, k][0]) * inv_2dx
    return ti.Vector(
        [
            dpsi_z_dy - dpsi_y_dz,
            dpsi_x_dz - dpsi_z_dx,
            dpsi_y_dx - dpsi_x_dy,
        ]
    )


@ti.func
def compute_curl_curl(
    wave_field: ti.template(),  # type: ignore
    i: ti.i32,  # type: ignore
    j: ti.i32,  # type: ignore
    k: ti.i32,  # type: ignore
):
    """
    Compute ∇×(∇×ψ) at voxel [i, j, k] via the vector identity
        ∇×(∇×ψ) = ∇(∇·ψ) − ∇²ψ

    The gradient of divergence is computed by central-differencing
    compute_divergence at the six face neighbors — that's what makes
    this a 2-cell halo operator (each face neighbor itself needs a 1-cell
    halo for its own divergence stencil, so [i±2, j, k] etc. are read).

    The Laplacian piece reuses the validated compute_laplacian.

    Used by M5.2 Close Eq. 19 linear limit  ∂²_t Q = −c² · ∇×(∇×Q).

    Args:
        wave_field: WaveField instance (reads psi_am)
        i, j, k: Voxel indices (must be interior with 2-cell halo:
            1 < i < nx-2, etc. — caller handles boundary skip)

    Returns:
        ∇×(∇×ψ) as Vector(3) in units [am / am²] = [1/am]
    """
    inv_2dx = 1.0 / (2.0 * wave_field.dx_am)

    # ∇(∇·ψ): central-difference the divergence at the six face neighbors
    div_xp = compute_divergence(wave_field, i + 1, j, k)
    div_xm = compute_divergence(wave_field, i - 1, j, k)
    div_yp = compute_divergence(wave_field, i, j + 1, k)
    div_ym = compute_divergence(wave_field, i, j - 1, k)
    div_zp = compute_divergence(wave_field, i, j, k + 1)
    div_zm = compute_divergence(wave_field, i, j, k - 1)

    grad_div = ti.Vector(
        [
            (div_xp - div_xm) * inv_2dx,
            (div_yp - div_ym) * inv_2dx,
            (div_zp - div_zm) * inv_2dx,
        ]
    )

    laplacian = compute_laplacian(wave_field, i, j, k)

    return grad_div - laplacian


# ================================================================
# MATRIX OPERATORS (M5.4 — Landau–de Gennes substrate M = O·D·O^T)
# ================================================================
# The Vector(3) toolkit above operates on the retiring ψ field. These operate on
# the real-symmetric 3×3 order parameter M stored in wave_field.M_am. All three
# primitives were de-risked in research/sandbox_v3/m5_3_matrix_feasibility.py
# (storage round-trip, commutator vs analytic, ti.sym_eig director recovery 0.9995)
# before this production port.
#
#   commutator         — [A, B] = A·B − B·A           pure matrix @ti.func
#   principal_director — M → (principal eigenvector, eigenvalues) via ti.sym_eig
#   matrix_laplacian   — ∇²M component-wise            → 3×3; 1-cell halo
#   antisym_A_mu       — A_μ = [M, ∂_μ M] (paper Eq.19) → (A_x, A_y, A_z); 1-cell halo
#   eigen_decompose    — full-grid kernel: M_am → director_nhat + eigenvalues
#
# INDEX-GENERIC NOTE (M5.8): these are written for 3×3 so the eventual 4×4 SO(1,3)
# promotion (D = diag(g,1,δ,0)) is a type change, not a rewrite. matrix_laplacian
# and antisym_A_mu are dimension-agnostic; only principal_director / eigen_decompose
# assume a 3-vector director (the 4D clock readout gets its own kernel in M5.8).
#
# WHAT M5.4 USES vs DEFERS: the M5.4 Coulomb gate is director-equivalent (relax the
# eigenvector with the existing M5.1 relax, rebuild M), so eigen_decompose is on the
# critical path and commutator / matrix_laplacian / antisym_A_mu are built + unit-
# tested here but NOT yet driving dynamics. The curvature F_μν = ∂_μA_ν − ∂_νA_μ
# (Eq.20) and the full Eq.18 action that consume A_μ land in M5.5.


@ti.func
def commutator(a, b):  # type: ignore
    """Matrix commutator [A, B] = A·B − B·A — building block of Eq.19 A_μ + Eq.20 F_μν."""
    return a @ b - b @ a


@ti.func
def principal_director(m):  # type: ignore
    """M (3×3 symmetric) → (principal eigenvector n̂, eigenvalues vec3 [λ₁≥λ₂≥λ₃]).

    ANALYTIC symmetric-3×3 eigensolver: Cardano closed-form eigenvalues + the
    principal eigenvector as the largest row cross-product of (M − λ₁I). Replaces
    `ti.sym_eig` (M5.6.5a-fix), which is accurate for UNIAXIAL (degenerate) M but
    CATASTROPHICALLY WRONG for BIAXIAL M on Metal/f32 — eigenvalue error ~0.48 vs
    numpy for distinct eigenvalues (and f64 is unavailable: SPIRV codegen fails).
    This analytic solver matches numpy to f32 precision on biaxial M (eigenvalue
    err ~6e-6, director err ~2e-7 over 20k random symmetric matrices); the M5.6
    biaxial substrate (3 distinct eigenvalues) needs it. Returns eigenvalues already
    sorted descending (λ₁≥λ₂≥λ₃).
    """
    a00, a11, a22 = m[0, 0], m[1, 1], m[2, 2]
    a01, a02, a12 = m[0, 1], m[0, 2], m[1, 2]
    q = (a00 + a11 + a22) / 3.0
    p1 = a01 * a01 + a02 * a02 + a12 * a12
    p2 = (a00 - q) ** 2 + (a11 - q) ** 2 + (a22 - q) ** 2 + 2.0 * p1
    p = ti.sqrt(ti.max(p2 / 6.0, 1e-30))                  # guard isotropic (p2→0)
    # B = (M − qI)/p ;  r = det(B)/2, clamped to [−1, 1] (numerical safety for acos)
    b00, b11, b22 = (a00 - q) / p, (a11 - q) / p, (a22 - q) / p
    b01, b02, b12 = a01 / p, a02 / p, a12 / p
    detB = (b00 * (b11 * b22 - b12 * b12)
            - b01 * (b01 * b22 - b12 * b02)
            + b02 * (b01 * b12 - b11 * b02))
    r = ti.max(-1.0, ti.min(1.0, detB / 2.0))
    phi = ti.acos(r) / 3.0
    two_pi_3 = 2.0943951023931953                        # 2π/3
    eig1 = q + 2.0 * p * ti.cos(phi)                     # largest
    eig3 = q + 2.0 * p * ti.cos(phi + two_pi_3)          # smallest
    eig2 = 3.0 * q - eig1 - eig3                         # middle
    # principal eigenvector = null space of (M − eig1·I): largest of the 3 row cross-products
    row0 = ti.Vector([a00 - eig1, a01, a02])
    row1 = ti.Vector([a01, a11 - eig1, a12])
    row2 = ti.Vector([a02, a12, a22 - eig1])
    c01, c02, c12 = row0.cross(row1), row0.cross(row2), row1.cross(row2)
    n01, n02, n12 = c01.norm_sqr(), c02.norm_sqr(), c12.norm_sqr()
    nvec = c01
    best = n01
    if n02 > best:
        nvec, best = c02, n02
    if n12 > best:
        nvec, best = c12, n12
    n = nvec / (ti.sqrt(best) + 1e-20)
    return n, ti.Vector([eig1, eig2, eig3])


@ti.func
def matrix_laplacian(
    wave_field: ti.template(),  # type: ignore
    i: ti.i32,  # type: ignore
    j: ti.i32,  # type: ignore
    k: ti.i32,  # type: ignore
):
    """∇²M at voxel [i, j, k] — the 6-point stencil applied entry-wise to M.

    Taichi matrix arithmetic is element-wise for + and scalar ·, so the same
    stencil that gives the Vector(3) Laplacian gives the matrix Laplacian with no
    per-component loop. Reads wave_field.M_am; returns a 3×3 matrix in [1/am].
    1-cell halo (caller skips boundary). Used by the matrix evolution (M5.5) and
    by operator-verification tests in M5.4.
    """
    face_sum = (
        wave_field.M_am[i + 1, j, k]
        + wave_field.M_am[i - 1, j, k]
        + wave_field.M_am[i, j + 1, k]
        + wave_field.M_am[i, j - 1, k]
        + wave_field.M_am[i, j, k + 1]
        + wave_field.M_am[i, j, k - 1]
    )
    center = wave_field.M_am[i, j, k]
    return (face_sum - 6.0 * center) / (wave_field.dx_am**2)


@ti.func
def antisym_A_mu(
    wave_field: ti.template(),  # type: ignore
    i: ti.i32,  # type: ignore
    j: ti.i32,  # type: ignore
    k: ti.i32,  # type: ignore
):
    """Paper Eq.19 antisymmetric connection A_μ = [M, ∂_μ M] for μ ∈ {x, y, z}.

    Central-differences M along each axis, then commutes with M at the voxel.
    Returns the three antisymmetric 3×3 matrices (A_x, A_y, A_z). The M5.3 spike
    verified the spatial-derivative commutator [∂_xM, ∂_yM] over a field; this is
    the same machinery in the [M, ∂_μM] arrangement. 1-cell halo. Consumed by the
    Eq.20 curvature F_μν and the Eq.18 action — both assembled in M5.5.
    """
    inv_2dx = 1.0 / (2.0 * wave_field.dx_am)
    m = wave_field.M_am[i, j, k]
    dM_x = (wave_field.M_am[i + 1, j, k] - wave_field.M_am[i - 1, j, k]) * inv_2dx
    dM_y = (wave_field.M_am[i, j + 1, k] - wave_field.M_am[i, j - 1, k]) * inv_2dx
    dM_z = (wave_field.M_am[i, j, k + 1] - wave_field.M_am[i, j, k - 1]) * inv_2dx
    return commutator(m, dM_x), commutator(m, dM_y), commutator(m, dM_z)


@ti.func
def compute_laplacian_director(
    wave_field: ti.template(),  # type: ignore
    i: ti.i32,  # type: ignore
    j: ti.i32,  # type: ignore
    k: ti.i32,  # type: ignore
):
    """∇²n̂ at voxel [i, j, k] — 6-point stencil on the DERIVED director field.

    Identical form to compute_laplacian, but reads wave_field.director_nhat (the
    matrix-substrate director, eigenvector of M) rather than the retiring ψ. The
    M5.4 director-equivalent relaxation gradient-descends the Frank energy on this
    field. 1-cell halo (caller skips boundary). (compute_laplacian on ψ stays for
    the retiring wave leapfrog until the M5.4 cleanup deletes it.)
    """
    face_sum = (
        wave_field.director_nhat[i + 1, j, k]
        + wave_field.director_nhat[i - 1, j, k]
        + wave_field.director_nhat[i, j + 1, k]
        + wave_field.director_nhat[i, j - 1, k]
        + wave_field.director_nhat[i, j, k + 1]
        + wave_field.director_nhat[i, j, k - 1]
    )
    center = wave_field.director_nhat[i, j, k]
    return (face_sum - 6.0 * center) / (wave_field.dx_am**2)


@ti.kernel
def rebuild_M_from_director(wave_field: ti.template(), delta: ti.f32):  # type: ignore
    """Rebuild the uniaxial order parameter M = δ·I + (1−δ)·n̂⊗n̂ from director_nhat.

    Called after the director-equivalent relaxation so M_am, the `‖M−D‖_F`/`‖Ṁ‖_F`
    trackers, and the flux-mesh orientation-deviation mode stay consistent with the
    relaxed director. Writes all three M buffers (M_prev = M = M_new) — BC-consistent
    for the M5.5 matrix leapfrog, and `Ṁ = 0` for a static relaxation (clock ω = 0).
    The uniaxial inverse of eigen_decompose.
    """
    eye = ti.Matrix.identity(ti.f32, 3)
    for i, j, k in wave_field.director_nhat:
        n = wave_field.director_nhat[i, j, k]
        m = delta * eye + (1.0 - delta) * n.outer_product(n)
        wave_field.M_am[i, j, k] = m
        wave_field.M_prev_am[i, j, k] = m
        wave_field.M_new_am[i, j, k] = m


@ti.kernel
def eigen_decompose(wave_field: ti.template()):  # type: ignore
    """THE LYNCHPIN (4b_rendering_features.md): refresh director_nhat + eigenvalues from M_am.

    Per voxel: eigen-decompose M_am, store the sorted spectrum (λ₁≥λ₂≥λ₃) in
    `eigenvalues` and the principal eigenvector in `director_nhat`. The whole
    rendering stack (glyphs, vector_warp, granule director) and the redefined
    amplitude tracker (‖M−D‖_F) read these derived fields, so this runs once/frame
    after M changes — never in an inner loop (M5.3 cost: ~1.1× a Vector(3)
    Laplacian at 256³).

    SIGN CONTINUITY: a nematic director is apolar (n ≡ −n), and `ti.sym_eig`
    returns an arbitrary eigenvector sign per voxel. A spurious sign flip between
    neighbours fakes a huge ∇n̂ (wrong Frank energy). We resolve it by aligning the
    fresh eigenvector to the EXISTING director_nhat (flip if n̂·n̂_prev < 0). The
    matrix seeders write director_nhat with the correct physical sign first, so the
    branch stays smooth from seed onward; on a cold field (n̂_prev = 0) the dot is 0
    and no flip occurs.
    """
    for i, j, k in wave_field.M_am:
        n, evals = principal_director(wave_field.M_am[i, j, k])
        # sort eigenvalues descending λ₁≥λ₂≥λ₃ (3-element selection sort)
        e0, e1, e2 = evals[0], evals[1], evals[2]
        lo = ti.min(e0, ti.min(e1, e2))
        hi = ti.max(e0, ti.max(e1, e2))
        mid = e0 + e1 + e2 - lo - hi
        wave_field.eigenvalues[i, j, k] = ti.Vector([hi, mid, lo])
        # sign-continuous principal director
        if n.dot(wave_field.director_nhat[i, j, k]) < 0.0:
            n = -n
        wave_field.director_nhat[i, j, k] = n


# ================================================================
# MATRIX EVOLUTION (M5.5) — Eq.18 action leapfrog
# ================================================================
# The first LIVE matrix-field dynamics: evolve M under Duda's Eq.18 action with a
# simple ½‖Ṁ‖² kinetic + the FAITHFUL Eq.18 potential (decision 2026-05-26):
#     U(M) = 4 Σ_{μ<ν} ‖[M_μ, M_ν]‖²_F  +  V_M(M)          (M_μ = ∂_μ M)
#     EOM:  ∂²_t M = c²·force_curv − dV_M/dM ,   force_curv = Σ_α ∂_α G_α
#     G_α = ∂U_curv/∂M_α = 8 Σ_ν [[M_α,M_ν], M_ν]          (symmetric)
# Validated in sandbox_v5 (m5_5_2/3); see 5a §9. The faithful curvature kinetic
# (F_μ0² = 4‖[M_μ,Ṁ]‖², degenerate metric) is the M5.6 refinement.
#
# Two-pass per step: compute_curvature_flux (G_α everywhere, 1-cell halo) → evolve_M
# (divergence of G + V-force, leapfrog, 2-cell halo) → wave_field.swap_matrix_buffers().


@ti.func
def V_M(m, a: ti.f32, b: ti.f32, c: ti.f32):  # type: ignore
    """Eq.13 Landau–de Gennes Higgs-like potential V = a·Tr(M²) − b·Tr(M³) + c·(Tr(M²))².

    Exact (a,b,c) giving the Λ=(1,δ,0) vacuum is Q7 (Duda open); a=b=c=0 → free-curvature
    dynamics (V off). Rotation-invariant ⇒ acts only on the eigenvalue/regularization
    sector (m5_5_3 finding).
    """
    m2 = m @ m
    tr2 = m2.trace()
    tr3 = (m2 @ m).trace()
    return a * tr2 - b * tr3 + c * tr2 * tr2


@ti.func
def dV_M(m, a: ti.f32, b: ti.f32, c: ti.f32):  # type: ignore
    """∂V_LG/∂M = 2a·M − 3b·M² + 4c·Tr(M²)·M  (clean matrix derivative for symmetric M)."""
    m2 = m @ m
    tr2 = m2.trace()
    return 2.0 * a * m - 3.0 * b * m2 + 4.0 * c * tr2 * m


@ti.kernel
def compute_curvature_flux(wave_field: ti.template()):  # type: ignore
    """Curvature flux G_α = 8 Σ_ν [[M_α,M_ν],M_ν] = ∂U_curv/∂M_α → curv_flux_{x,y,z}.

    M_α = ∂_α M (central diff of M_am). [[A,B],B] is symmetric (antisym⊗sym), so each
    G_α is symmetric. 1-cell halo (caller skips boundary). Run before evolve_M.
    """
    nx, ny, nz = wave_field.nx, wave_field.ny, wave_field.nz
    inv_2dx = 1.0 / (2.0 * wave_field.dx_am)
    for i, j, k in ti.ndrange((1, nx - 1), (1, ny - 1), (1, nz - 1)):
        mx = (wave_field.M_am[i + 1, j, k] - wave_field.M_am[i - 1, j, k]) * inv_2dx
        my = (wave_field.M_am[i, j + 1, k] - wave_field.M_am[i, j - 1, k]) * inv_2dx
        mz = (wave_field.M_am[i, j, k + 1] - wave_field.M_am[i, j, k - 1]) * inv_2dx
        # G_α = 8 Σ_{ν≠α} [[M_α,M_ν],M_ν]  (ν=α term is zero)
        wave_field.curv_flux_x[i, j, k] = 8.0 * (
            commutator(commutator(mx, my), my) + commutator(commutator(mx, mz), mz))
        wave_field.curv_flux_y[i, j, k] = 8.0 * (
            commutator(commutator(my, mx), mx) + commutator(commutator(my, mz), mz))
        wave_field.curv_flux_z[i, j, k] = 8.0 * (
            commutator(commutator(mz, mx), mx) + commutator(commutator(mz, my), my))


@ti.kernel
def evolve_M(
    wave_field: ti.template(),  # type: ignore
    c_amrs: ti.f32,  # type: ignore
    dt_rs: ti.f32,  # type: ignore
    a: ti.f32,  # type: ignore
    b: ti.f32,  # type: ignore
    c: ti.f32,  # type: ignore
):
    """Leapfrog M one step under the Eq.18 action: ∂²_t M = c²·(Σ_α ∂_α G_α) − dV_M/dM.

    M_new = 2M − M_prev + (dt)²·[c²·div(G) − dV_M(M)].  Reads M_am, M_prev_am,
    curv_flux_* (from compute_curvature_flux); writes M_new_am. Caller then calls
    wave_field.swap_matrix_buffers(). 2-cell halo (G is a 1-cell-halo quantity, its
    divergence needs one more). Dirichlet BC: boundary M_new left untouched (seeders
    write all three M buffers, per the triple-buffer BC rule).
    """
    nx, ny, nz = wave_field.nx, wave_field.ny, wave_field.nz
    inv_2dx = 1.0 / (2.0 * wave_field.dx_am)
    c2 = c_amrs * c_amrs
    dt2 = dt_rs * dt_rs
    for i, j, k in ti.ndrange((2, nx - 2), (2, ny - 2), (2, nz - 2)):
        div_G = (
            (wave_field.curv_flux_x[i + 1, j, k] - wave_field.curv_flux_x[i - 1, j, k])
            + (wave_field.curv_flux_y[i, j + 1, k] - wave_field.curv_flux_y[i, j - 1, k])
            + (wave_field.curv_flux_z[i, j, k + 1] - wave_field.curv_flux_z[i, j, k - 1])
        ) * inv_2dx
        force = c2 * div_G - dV_M(wave_field.M_am[i, j, k], a, b, c)
        wave_field.M_new_am[i, j, k] = (
            2.0 * wave_field.M_am[i, j, k] - wave_field.M_prev_am[i, j, k] + dt2 * force
        )


# ================================================================
# ψ EVOLVING ENGINE — LAGRANGIAN FIELD EVOLUTION  (DORMANT LEGACY, M5.4)
# ================================================================
# DORMANT LEGACY (M5.4 vector→matrix migration): evolve_psi (wave leapfrog) + V_psi
# are RETIRED from the live path — the launcher's compute_propagation is now a no-op,
# and the matrix field M is static in M5.4 (seed + relax). The matrix-field leapfrog
# (∂²_t M from the Eq.18 action) replaces this in M5.5. Retained, not deleted, because
# the M5.0–M5.3 research validation scripts (sandbox_v2/v3) drive evolve_psi directly,
# and the shared ψ differential operators (compute_laplacian / divergence / curl /
# curl_curl, above) get repointed onto M for Close's Eq.23 at M5.7.
#
# evolve_psi evolves the field ψ via leapfrog/Verlet integration of the
# wave equation:
#     ∂²ψ/∂t² = c²·∇²ψ − ∂V/∂ψ
#
# In M5.0d this lands with V(ψ) = 0 (free-wave). The nonlinear potential term
# −∂V/∂ψ is added in M5.2 (Close's Eq. 23 + Klein-Gordon mass term + LdG).
#
# Why "evolve_psi" not "propagate_wave": the operation is field evolution,
# not wave-specific. Leapfrog also evolves topology (defect drift), gradient
# descent for relaxation reuses the same buffers, etc. Naming follows what the
# function operates on (the field ψ), not the method (leapfrog) or one of the
# behaviors (wave propagation).
#
# Caller MUST call wave_field.swap_buffers() after this kernel returns, to
# cycle the triple buffer: psi_prev ← psi, psi ← psi_new.


@ti.kernel
def evolve_psi(
    wave_field: ti.template(),  # type: ignore
    c_amrs: ti.f32,  # type: ignore
    dt_rs: ti.f32,  # type: ignore
    m_freq_rs: ti.f32,  # type: ignore
    lambda_phi4: ti.f32,  # type: ignore
):
    """
    Evolve ψ one timestep via leapfrog/Verlet integration.

    Klein-Gordon + Mexican-hat φ⁴ (V = ½·m²·|ψ|² + ¼·λ·(|ψ|² − 1)²):
        ∂²ψ/∂t² = c²·∇²ψ − m²·ψ − λ·(|ψ|² − 1)·ψ

    Discrete (leapfrog/Verlet):
        ψ_new = 2·ψ − ψ_prev + (c·dt)²·∇²ψ − (m·dt)²·ψ − (dt)²·λ·(|ψ|² − 1)·ψ

    Setting both m_freq_rs = 0 and lambda_phi4 = 0 recovers the free-wave.

    Reads:  wave_field.psi_am, wave_field.psi_prev_am
    Writes: wave_field.psi_new_am

    Boundary: fixed-value Dirichlet — boundary voxels are never updated;
    they stay at whatever value the active seeder wrote. The interpretation
    of that value depends on the seeder:
      - M5.0 wave seeders (seed_gaussian, seed_dispersion_modes) leave
        boundaries at ψ ≈ 0 (zero-displacement BC; quiet wave boundary)
      - M5.1+ topology seeders (seed_vacuum, seed_hedgehog) leave boundaries
        at ψ = ẑ (vacuum-director BC; far-field defect reference)
    Both are valid Dirichlet BCs at different fixed values. The propagator
    is agnostic — it just iterates (1, n−1) on each axis (the 6-point
    Laplacian's 1-cell halo requirement) and never touches the boundary.

    BC consistency requirement: ALL THREE BUFFERS (psi_prev_am, psi_am,
    psi_new_am) must hold the same boundary value before the first call.
    psi_new_am's boundary is never written by the propagator, so swap_buffers'
    `psi_am.copy_from(psi_new_am)` would otherwise clobber the seeded
    boundary value with psi_new_am's default zero. Seeders fix this by
    writing all three buffers; verified bug 2026-05-09 with seed_hedgehog.

    CFL stability: dt ≤ dx / (c·√3) for 3D wave equation. Caller is
    responsible for sizing dt below this bound (see _launcher CFL eval).

    Args:
        wave_field: WaveField instance (reads psi_am, psi_prev_am; writes psi_new_am)
        c_amrs: wave speed in scaled units (am/rs); already includes any
            slow-motion factor for visualization
        dt_rs: timestep in rontoseconds, sized below the CFL bound
        m_freq_rs: Klein-Gordon mass-frequency m·c²/ℏ in rad/rs storage units
            (electron: ~7.76e-7 rad/rs at SIM_SPEED=1; 0 disables the mass term)
        lambda_phi4: Mexican-hat coupling in 1/(am²·rs²). Set in launcher to
            (c_amrs/dx_am)² for a natural grid-scale restoring force; 0
            disables the φ⁴ term.
    """
    nx, ny, nz = wave_field.nx, wave_field.ny, wave_field.nz
    c_dt_squared = (c_amrs * dt_rs) ** 2
    m_dt_squared = (m_freq_rs * dt_rs) ** 2
    lambda_dt_squared = lambda_phi4 * dt_rs * dt_rs

    # Interior voxels only (Dirichlet BC; 6-point Laplacian needs 1-cell halo)
    for i, j, k in ti.ndrange((1, nx - 1), (1, ny - 1), (1, nz - 1)):
        laplacian_psi_am = compute_laplacian(wave_field, i, j, k)
        psi_ijk = wave_field.psi_am[i, j, k]
        psi_sq = psi_ijk.norm_sqr()
        phi4_factor = lambda_dt_squared * (psi_sq - 1.0)

        # Leapfrog/Verlet update:
        # ψ_new = 2·ψ − ψ_prev + (c·dt)²·∇²ψ − (m·dt)²·ψ − (dt)²·λ·(|ψ|²−1)·ψ
        wave_field.psi_new_am[i, j, k] = (
            2.0 * psi_ijk
            - wave_field.psi_prev_am[i, j, k]
            + c_dt_squared * laplacian_psi_am
            - m_dt_squared * psi_ijk
            - phi4_factor * psi_ijk
        )


# ================================================================
# POTENTIAL ENERGY HOOK — V(ψ)
# ================================================================
# Default V(ψ) = 0 (free-wave) in M5.0g. M5.2 swaps this implementation in to
# add the nonlinear physics:
#   - Klein-Gordon mass:    V = ½ m² |ψ|²
#   - Close Eq. 23 nonlin:  V = … (from −u·∇s + w×s couplings)
#   - LdG biaxial:          V = polynomial in (Q, axes δ/1/g)
#
# --- NATURAL-UNIT BASIS — applied here in V(ψ) -----------------------------
# Basis convention: ℏ = c = λ̄_C(electron) = 1, where λ̄_C = ℏ/(m_e·c) is the
# reduced Compton wavelength of the electron (constants.COMPTON_WAVELENGTH_
# REDUCED_ELECTRON_AM). This makes m_e = 1 in natural units and aligns
# dimensional couplings with textbook field-theory form.
#
# Single bridge value: the Klein-Gordon mass-frequency m·c²/ℏ in rad/rs.
# Computed once in launcher as
#     m_freq_kg_rs = c_amrs / COMPTON_WAVELENGTH_REDUCED_ELECTRON_AM
# (electron rest-frame angular frequency, ~7.76e-7 rad/rs at SIM_SPEED=1).
# Threaded into V_psi via compute_energyH_density and (Step 2+) evolve_psi.
# For ψ in am, V_storage = ½·(m_freq_rs)²·|ψ|² produces an energy density in
# (am/rs)² — same units as the kinetic and gradient terms.
#
# M5.2 Step 4a (this version): plain KG + Mexican-hat φ⁴ Mexican-hat active —
#     V(ψ) = ½·m²·|ψ|² + ¼·λ·(|ψ|² − 1)²
# EL derivative: ∂V/∂ψ = m²·ψ + λ·(|ψ|² − 1)·ψ
# evolve_psi adds −(dt)²·[m²·ψ + λ·(|ψ|² − 1)·ψ] to the leapfrog.
#
# The φ⁴ Mexican-hat moves the potential minimum from ψ = 0 (plain KG) to the
# unit sphere |ψ| = 1 — the right shape for stabilizing a unit-vector
# director field. Setting both m_freq_rs = 0 AND lambda_phi4 = 0 recovers
# the free-wave (V = 0) used by the M5.0h dispersion tests.
#
# DERRICK'S THEOREM CAVEAT: φ⁴ alone (without a 4th-derivative Skyrme term)
# does NOT prevent the hedgehog core from collapsing to a point in 3D
# (E_kinetic ~ R + E_potential·R³ both monotone in R). M5.2 Step 4b will
# add the Skyrme term if 4a shows the hedgehog still collapses.
#
# Linear kernels (Laplacian, divergence, curl, curl-curl) stay in storage
# units (am, rs, rHz) throughout — dimensionally self-balancing, don't
# benefit from natural units. See M5.0f decision-record in 2a_path_to_m5.md.


@ti.func
def V_psi(
    psi: ti.template(),  # type: ignore
    m_freq_rs: ti.f32,  # type: ignore
    lambda_phi4: ti.f32,  # type: ignore
):
    """
    Scalar potential V(ψ) at one voxel — KG mass + Mexican-hat φ⁴ (M5.2 Step 4a).

    V(ψ) = ½·m²·|ψ|² + ¼·λ·(|ψ|² − 1)²
    EL contribution: ∂V/∂ψ = m²·ψ + λ·(|ψ|² − 1)·ψ
    evolve_psi adds −(dt)²·[m²·ψ + λ·(|ψ|² − 1)·ψ] to the leapfrog.

    Mexican-hat term: minimum is the unit sphere |ψ| = 1 (V = 0); maximum at
    ψ = 0 (V = ¼λ). Right shape to stabilize a director field whose vacuum
    is the unit sphere. Setting λ = 0 reduces to plain KG; setting both = 0
    recovers the free-wave (V = 0). |ψ|² = 1 is the implicit vacuum convention
    of seed_hedgehog / seed_vacuum (boundary at n̂ = ẑ, |n̂| = 1).

    Args:
        psi: Vector(3) field value at the voxel
        m_freq_rs: Klein-Gordon mass-frequency m·c²/ℏ in rad/rs storage units
            (electron: ~7.76e-7 rad/rs at SIM_SPEED=1; 0 disables the KG term)
        lambda_phi4: Mexican-hat coupling in 1/(am²·rs²). Set in launcher to
            (c_amrs/dx_am)² → φ⁴ "restoring force" matches grid scale. 0
            disables the term.

    Returns:
        scalar V(ψ) in the same units as kinetic and gradient terms (am²/rs²)
    """
    psi_sq = psi.norm_sqr()
    v_kg = 0.5 * m_freq_rs * m_freq_rs * psi_sq
    diff = psi_sq - 1.0
    v_phi4 = 0.25 * lambda_phi4 * diff * diff
    return v_kg + v_phi4


# ================================================================
# GRADIENT-DESCENT RELAXATION (M5.1 task 6)
# ================================================================
# Direct port of Exp 2's relax() inner loop to Taichi:
#     dn = ∇²n − (n·∇²n) n         (tangent projection onto unit sphere)
#     n_new = n + τ · dn           (gradient-descent step)
#     n_new ← n_new / |n_new|       (unit-length constraint)
#     n_new[defect_core] = ±ẑ      (soft Dirichlet pin)
#
# Stability: τ < dx²/(2·dim·K) by the standard heat-equation CFL condition.
# For dim=3 in 3D the bound is τ < dx²/(6·K). Use ~50% of the bound for
# headroom; the launcher computes this from wave_field.dx_am at runtime.
#
# Output buffer convention: writes to psi_new_am (same buffer leapfrog uses).
# The caller copies psi_new_am to BOTH psi_am AND psi_prev_am so that
# subsequent leapfrog steps see ψ̇ = 0 (no spurious time-derivative from the
# relaxation update — relaxation is a STATIC operation, not a dynamics step).
# This is different from swap_buffers, which rotates buffers for time evolution.
#
# Boundary handling: psi_new_am boundary voxels are copied from psi_am (no
# update) — preserves the fixed-value Dirichlet BC (e.g. ẑ at the universe
# edge from seed_vacuum).
#
# Core pin: soft Dirichlet at the SINGLE closest voxel to each defect center.
# Without this, gradient descent can numerically dissolve the defect on a
# discrete grid (topology is not strictly preserved by discretization).
#
# M7 FORWARD-LINK: this kernel is mathematically the γ → ∞ limit of a damped
# wave equation `∂²_t ψ = c²∇²ψ − γ·∂_t ψ`. Although the relaxation use is
# purely numerical (no physics interpretation in M5.1), the same primitives
# — Laplacian stencil, tangent projection, soft-core pin — will be reused
# in M7 thermal-modulation kernels where γ becomes a PHYSICAL damping
# coefficient (radiation loss, phonon coupling, EM-load impedance). See
# `research/2a_path_to_m5.md § Beyond M6 — thermal mechanics pathway`
# for the infrastructure-foundation discussion.


@ti.kernel
def relax_director_step(
    wave_field: ti.template(),  # type: ignore
    tau: ti.f32,  # type: ignore
    pin_centers: ti.template(),  # type: ignore
    pin_signs: ti.template(),  # type: ignore
    n_defects: ti.i32,  # type: ignore
):
    """
    One gradient-descent step on the Frank elastic energy `H_F = (K/2)·|∇n̂|²`.

    Reads:  wave_field.director_nhat
    Writes: wave_field.director_nhat_new  (caller copies → director_nhat)

    M5.4: repointed from the retiring ψ buffers onto the matrix-substrate director
    (director_nhat = principal eigenvector of M). The math is byte-identical to the
    M5.1 relaxation — tangent-projected gradient descent on (K/2)|∇n̂|² + soft-core
    pin — so the Coulomb result carries over. Caller rebuilds M from the relaxed
    director via rebuild_M_from_director.

    Args:
        wave_field: WaveField with the director field in director_nhat
        tau: step size (must satisfy τ < dx²/(6·K) for stability)
        pin_centers: ti.field shape (n_defects, 3) i32 — defect centers in voxel coords
        pin_signs: ti.field shape (n_defects,) i32 — ±1 per defect (pin direction)
        n_defects: number of active defects in pin_centers/signs

    Note: K_frank is implicit in `tau`'s CFL derivation — the kernel itself
    doesn't multiply by K because gradient descent on `(K/2)·|∇n̂|²` gives
    `∂n/∂τ = K·∇²n` and the K can be absorbed into the step size (τ_eff = K·τ).
    For K=1 (M5.1 default) this is τ_eff = τ.
    """
    nx, ny, nz = wave_field.nx, wave_field.ny, wave_field.nz

    # ---- Interior: tangent-projected gradient step + renormalization ----
    for i, j, k in ti.ndrange((1, nx - 1), (1, ny - 1), (1, nz - 1)):
        n = wave_field.director_nhat[i, j, k]
        lap = compute_laplacian_director(wave_field, i, j, k)
        n_dot_lap = n.dot(lap)
        dn = lap - n_dot_lap * n  # tangent projection (remove component along n)
        n_new = n + tau * dn
        norm = n_new.norm() + 1e-12
        wave_field.director_nhat_new[i, j, k] = n_new / norm

    # ---- Boundary: preserve Dirichlet BC by copying from current director ----
    for j, k in ti.ndrange(ny, nz):
        wave_field.director_nhat_new[0, j, k] = wave_field.director_nhat[0, j, k]
        wave_field.director_nhat_new[nx - 1, j, k] = wave_field.director_nhat[nx - 1, j, k]
    for i, k in ti.ndrange(nx, nz):
        wave_field.director_nhat_new[i, 0, k] = wave_field.director_nhat[i, 0, k]
        wave_field.director_nhat_new[i, ny - 1, k] = wave_field.director_nhat[i, ny - 1, k]
    for i, j in ti.ndrange(nx, ny):
        wave_field.director_nhat_new[i, j, 0] = wave_field.director_nhat[i, j, 0]
        wave_field.director_nhat_new[i, j, nz - 1] = wave_field.director_nhat[i, j, nz - 1]

    # ---- Pin defect cores (soft Dirichlet at single closest voxel) ----
    for d in range(n_defects):
        ci = pin_centers[d, 0]
        cj = pin_centers[d, 1]
        ck = pin_centers[d, 2]
        sgn = ti.cast(pin_signs[d], ti.f32)
        wave_field.director_nhat_new[ci, cj, ck] = ti.Vector([0.0, 0.0, sgn])
