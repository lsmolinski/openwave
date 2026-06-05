"""
M5 ENGINE — PDE (engine2_pde)

The complete PDE layer: the reusable differential operators AND the equations
integrated with them — merged from the former operators + PDE modules
(2026-05-25 SoC consolidation). Two internal sections below:
  MATRIX OPERATORS       — @ti.func commutator, matrix ∇², antisym Â_μ, eigensolve
  PDE DYNAMICS           — evolve_M (matrix leapfrog), V_M / dV_M (LdG potential),
                           compute_curvature_flux, eigen_decompose,
                           relax_director_step (dissipative gradient descent)

Concern order: engine1_seeds → engine2_pde → engine3_observables → engine4_render.
M5.4 watch-item: if eigen_decompose ends up shared by render + observables + the
PDE, re-extract a small operators module at that point.
"""

import taichi as ti

# ================================================================
# MATRIX OPERATORS (M5.4 — Landau–de Gennes substrate M = O·D·O^T)
# ================================================================
# The Vector(3) toolkit above operates on the retiring ψ field. These operate on
# the real-symmetric 3×3 order parameter M stored in tensor_field.M_am. All three
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
    p = ti.sqrt(ti.max(p2 / 6.0, 1e-30))  # guard isotropic (p2→0)
    # B = (M − qI)/p ;  r = det(B)/2, clamped to [−1, 1] (numerical safety for acos)
    b00, b11, b22 = (a00 - q) / p, (a11 - q) / p, (a22 - q) / p
    b01, b02, b12 = a01 / p, a02 / p, a12 / p
    detB = (
        b00 * (b11 * b22 - b12 * b12)
        - b01 * (b01 * b22 - b12 * b02)
        + b02 * (b01 * b12 - b11 * b02)
    )
    r = ti.max(-1.0, ti.min(1.0, detB / 2.0))
    phi = ti.acos(r) / 3.0
    two_pi_3 = 2.0943951023931953  # 2π/3
    eig1 = q + 2.0 * p * ti.cos(phi)  # largest
    eig3 = q + 2.0 * p * ti.cos(phi + two_pi_3)  # smallest
    eig2 = 3.0 * q - eig1 - eig3  # middle
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
def eigvec_for(m, lam):  # type: ignore
    """Unit eigenvector of symmetric 3×3 `m` for a KNOWN eigenvalue `lam`.

    Null space of (m − lam·I), taken as the largest of the 3 row cross-products
    (same construction as the principal eigenvector in `principal_director`, but
    for an arbitrary supplied eigenvalue — used for the MIDDLE eigenvector, the
    δ "clock-hand" axis, VIZ.3). Robust when `lam` is well-separated; degenerate
    (repeated-eigenvalue) cases give an arbitrary-but-unit vector in the
    eigenspace, which is the correct apolar behaviour for a uniaxial M.
    """
    a00, a11, a22 = m[0, 0], m[1, 1], m[2, 2]
    a01, a02, a12 = m[0, 1], m[0, 2], m[1, 2]
    row0 = ti.Vector([a00 - lam, a01, a02])
    row1 = ti.Vector([a01, a11 - lam, a12])
    row2 = ti.Vector([a02, a12, a22 - lam])
    c01, c02, c12 = row0.cross(row1), row0.cross(row2), row1.cross(row2)
    n01, n02, n12 = c01.norm_sqr(), c02.norm_sqr(), c12.norm_sqr()
    vec = c01
    best = n01
    if n02 > best:
        vec, best = c02, n02
    if n12 > best:
        vec, best = c12, n12
    return vec / (ti.sqrt(best) + 1e-20)


@ti.func
def matrix_laplacian(
    tensor_field: ti.template(),  # type: ignore
    i: ti.i32,  # type: ignore
    j: ti.i32,  # type: ignore
    k: ti.i32,  # type: ignore
):
    """∇²M at voxel [i, j, k] — the 6-point stencil applied entry-wise to M.

    Taichi matrix arithmetic is element-wise for + and scalar ·, so the same
    stencil that gives the Vector(3) Laplacian gives the matrix Laplacian with no
    per-component loop. Reads tensor_field.M_am; returns a 3×3 matrix in [1/am].
    1-cell halo (caller skips boundary). Used by the matrix evolution (M5.5) and
    by operator-verification tests in M5.4.
    """
    face_sum = (
        tensor_field.M_am[i + 1, j, k]
        + tensor_field.M_am[i - 1, j, k]
        + tensor_field.M_am[i, j + 1, k]
        + tensor_field.M_am[i, j - 1, k]
        + tensor_field.M_am[i, j, k + 1]
        + tensor_field.M_am[i, j, k - 1]
    )
    center = tensor_field.M_am[i, j, k]
    return (face_sum - 6.0 * center) / (tensor_field.dx_am**2)


@ti.func
def antisym_A_mu(
    tensor_field: ti.template(),  # type: ignore
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
    inv_2dx = 1.0 / (2.0 * tensor_field.dx_am)
    m = tensor_field.M_am[i, j, k]
    dM_x = (tensor_field.M_am[i + 1, j, k] - tensor_field.M_am[i - 1, j, k]) * inv_2dx
    dM_y = (tensor_field.M_am[i, j + 1, k] - tensor_field.M_am[i, j - 1, k]) * inv_2dx
    dM_z = (tensor_field.M_am[i, j, k + 1] - tensor_field.M_am[i, j, k - 1]) * inv_2dx
    return commutator(m, dM_x), commutator(m, dM_y), commutator(m, dM_z)


@ti.func
def compute_laplacian_director(
    tensor_field: ti.template(),  # type: ignore
    i: ti.i32,  # type: ignore
    j: ti.i32,  # type: ignore
    k: ti.i32,  # type: ignore
):
    """∇²n̂ at voxel [i, j, k] — 6-point stencil on the DERIVED director field.

    Identical form to compute_laplacian, but reads tensor_field.director_nhat (the
    matrix-substrate director, eigenvector of M) rather than the retiring ψ. The
    M5.4 director-equivalent relaxation gradient-descends the Frank energy on this
    field. 1-cell halo (caller skips boundary). (compute_laplacian on ψ stays for
    the retiring wave leapfrog until the M5.4 cleanup deletes it.)
    """
    face_sum = (
        tensor_field.director_nhat[i + 1, j, k]
        + tensor_field.director_nhat[i - 1, j, k]
        + tensor_field.director_nhat[i, j + 1, k]
        + tensor_field.director_nhat[i, j - 1, k]
        + tensor_field.director_nhat[i, j, k + 1]
        + tensor_field.director_nhat[i, j, k - 1]
    )
    center = tensor_field.director_nhat[i, j, k]
    return (face_sum - 6.0 * center) / (tensor_field.dx_am**2)


@ti.kernel
def rebuild_M_from_director(tensor_field: ti.template(), delta: ti.f32):  # type: ignore
    """Rebuild the uniaxial order parameter M = δ·I + (1−δ)·n̂⊗n̂ from director_nhat.

    Called after the director-equivalent relaxation so M_am, the `‖M−D‖_F`/`‖Ṁ‖_F`
    trackers, and the flux-mesh orientation-deviation mode stay consistent with the
    relaxed director. Writes all three M buffers (M_prev = M = M_new) — BC-consistent
    for the M5.5 matrix leapfrog, and `Ṁ = 0` for a static relaxation (clock ω = 0).
    The uniaxial inverse of eigen_decompose.
    """
    eye = ti.Matrix.identity(ti.f32, 3)
    g = tensor_field.lc_g  # M5.8.1: time-axis (index 3) eigenvalue
    for i, j, k in tensor_field.director_nhat:
        n = tensor_field.director_nhat[i, j, k]
        msp = delta * eye + (1.0 - delta) * n.outer_product(n)  # 3×3 spatial
        m = ti.Matrix.zero(ti.f32, 4, 4)  # embed block-diag(spatial, g)
        for a in ti.static(range(3)):
            for bb in ti.static(range(3)):
                m[a, bb] = msp[a, bb]
        m[3, 3] = g
        tensor_field.M_am[i, j, k] = m
        tensor_field.M_prev_am[i, j, k] = m
        tensor_field.M_new_am[i, j, k] = m


@ti.kernel
def eigen_decompose(tensor_field: ti.template()):  # type: ignore
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
    for i, j, k in tensor_field.M_am:
        m = tensor_field.M_am[i, j, k]
        n, evals = principal_director(m)
        # sort eigenvalues descending λ₁≥λ₂≥λ₃ (3-element selection sort)
        e0, e1, e2 = evals[0], evals[1], evals[2]
        lo = ti.min(e0, ti.min(e1, e2))
        hi = ti.max(e0, ti.max(e1, e2))
        mid = e0 + e1 + e2 - lo - hi
        tensor_field.eigenvalues[i, j, k] = ti.Vector([hi, mid, lo])
        # sign-continuous principal director
        if n.dot(tensor_field.director_nhat[i, j, k]) < 0.0:
            n = -n
        tensor_field.director_nhat[i, j, k] = n
        # VIZ.3: middle (δ) eigenvector = the "clock-hand" axis that sweeps around
        # the director under the Zitterbewegung twist. Same apolar sign-continuity.
        nm = eigvec_for(m, mid)
        if nm.dot(tensor_field.director_mid[i, j, k]) < 0.0:
            nm = -nm
        tensor_field.director_mid[i, j, k] = nm


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
# (divergence of G + V-force, leapfrog, 2-cell halo) → tensor_field.swap_matrix_buffers().


@ti.func
def V_M(m, a: ti.f32, b: ti.f32, c: ti.f32):  # type: ignore
    """Eq.13 Landau–de Gennes Higgs-like potential V = a·Tr(M²) − b·Tr(M³) + c·(Tr(M²))².

    Exact (a,b,c) giving the Λ=(1,δ,0) vacuum is Q7 (Duda open); a=b=c=0 → free-curvature
    dynamics (V off). Rotation-invariant ⇒ acts only on the eigenvalue/regularization
    sector (m5_5_3 finding).
    """
    # M5.8.1 — act on the SPATIAL 3×3 block ONLY. The time axis (index 3, eigenvalue g)
    # is boost-decoupled and must NOT feed Tr(M²)/Tr(M³): including g²(=64)/g³(=512)
    # would inflate the LdG potential and its force by orders of magnitude → blow-up.
    msp = ti.Matrix(
        [[m[0, 0], m[0, 1], m[0, 2]], [m[1, 0], m[1, 1], m[1, 2]], [m[2, 0], m[2, 1], m[2, 2]]]
    )
    m2 = msp @ msp
    tr2 = m2.trace()
    tr3 = (m2 @ msp).trace()
    return a * tr2 - b * tr3 + c * tr2 * tr2


@ti.func
def dV_M(m, a: ti.f32, b: ti.f32, c: ti.f32):  # type: ignore
    """∂V_LG/∂M = 2a·M − 3b·M² + 4c·Tr(M²)·M, on the SPATIAL 3×3 block; time row/col
    force = 0 (g decoupled — M5.8.1). Mirrors V_M's spatial-only restriction."""
    msp = ti.Matrix(
        [[m[0, 0], m[0, 1], m[0, 2]], [m[1, 0], m[1, 1], m[1, 2]], [m[2, 0], m[2, 1], m[2, 2]]]
    )
    m2 = msp @ msp
    tr2 = m2.trace()
    dsp = 2.0 * a * msp - 3.0 * b * m2 + 4.0 * c * tr2 * msp  # 3×3
    d4 = ti.Matrix.zero(ti.f32, 4, 4)
    for i in ti.static(range(3)):
        for j in ti.static(range(3)):
            d4[i, j] = dsp[i, j]
    return d4


@ti.kernel
def compute_curvature_flux(tensor_field: ti.template()):  # type: ignore
    """Curvature flux G_α = 8 Σ_ν [[M_α,M_ν],M_ν] = ∂U_curv/∂M_α → curv_flux_{x,y,z}.

    M_α = ∂_α M (central diff of M_am). [[A,B],B] is symmetric (antisym⊗sym), so each
    G_α is symmetric. 1-cell halo (caller skips boundary). Run before evolve_M.
    """
    nx, ny, nz = tensor_field.nx, tensor_field.ny, tensor_field.nz
    inv_2dx = 1.0 / (2.0 * tensor_field.dx_am)
    for i, j, k in ti.ndrange((1, nx - 1), (1, ny - 1), (1, nz - 1)):
        mx = (tensor_field.M_am[i + 1, j, k] - tensor_field.M_am[i - 1, j, k]) * inv_2dx
        my = (tensor_field.M_am[i, j + 1, k] - tensor_field.M_am[i, j - 1, k]) * inv_2dx
        mz = (tensor_field.M_am[i, j, k + 1] - tensor_field.M_am[i, j, k - 1]) * inv_2dx
        # G_α = 8 Σ_{ν≠α} [[M_α,M_ν],M_ν]  (ν=α term is zero)
        tensor_field.curv_flux_x[i, j, k] = 8.0 * (
            commutator(commutator(mx, my), my) + commutator(commutator(mx, mz), mz)
        )
        tensor_field.curv_flux_y[i, j, k] = 8.0 * (
            commutator(commutator(my, mx), mx) + commutator(commutator(my, mz), mz)
        )
        tensor_field.curv_flux_z[i, j, k] = 8.0 * (
            commutator(commutator(mz, mx), mx) + commutator(commutator(mz, my), my)
        )


@ti.kernel
def evolve_M(
    tensor_field: ti.template(),  # type: ignore
    c_amrs: ti.f32,  # type: ignore
    dt_rs: ti.f32,  # type: ignore
    a: ti.f32,  # type: ignore
    b: ti.f32,  # type: ignore
    c: ti.f32,  # type: ignore
):
    """Leapfrog M one step under the Eq.18 action: ∂²_t M = c²·(Σ_α ∂_α G_α) − dV_M/dM.

    M_new = 2M − M_prev + (dt)²·[c²·div(G) − dV_M(M)].  Reads M_am, M_prev_am,
    curv_flux_* (from compute_curvature_flux); writes M_new_am. Caller then calls
    tensor_field.swap_matrix_buffers(). 2-cell halo (G is a 1-cell-halo quantity, its
    divergence needs one more). Dirichlet BC: boundary M_new left untouched (seeders
    write all three M buffers, per the triple-buffer BC rule).
    """
    nx, ny, nz = tensor_field.nx, tensor_field.ny, tensor_field.nz
    inv_2dx = 1.0 / (2.0 * tensor_field.dx_am)
    c2 = c_amrs * c_amrs
    dt2 = dt_rs * dt_rs
    for i, j, k in ti.ndrange((2, nx - 2), (2, ny - 2), (2, nz - 2)):
        div_G = (
            (tensor_field.curv_flux_x[i + 1, j, k] - tensor_field.curv_flux_x[i - 1, j, k])
            + (tensor_field.curv_flux_y[i, j + 1, k] - tensor_field.curv_flux_y[i, j - 1, k])
            + (tensor_field.curv_flux_z[i, j, k + 1] - tensor_field.curv_flux_z[i, j, k - 1])
        ) * inv_2dx
        force = c2 * div_G - dV_M(tensor_field.M_am[i, j, k], a, b, c)
        m_new = 2.0 * tensor_field.M_am[i, j, k] - tensor_field.M_prev_am[i, j, k] + dt2 * force
        # M5.8.1 — freeze the time axis (index 3): a constant, boost-decoupled g
        # background. The curvature force already preserves the block; this also pins
        # it under V-on (where dV_M would otherwise nudge M[3,3]). M5.8.2 replaces this
        # with the real Minkowski time dynamics.
        g_here = tensor_field.M_am[i, j, k][3, 3]
        m_new[0, 3] = 0.0
        m_new[3, 0] = 0.0
        m_new[1, 3] = 0.0
        m_new[3, 1] = 0.0
        m_new[2, 3] = 0.0
        m_new[3, 2] = 0.0
        m_new[3, 3] = g_here
        tensor_field.M_new_am[i, j, k] = m_new


# ================================================================
# M5.8.2c — 4D MINKOWSKI EVOLUTION (flag-gated; the time axis goes LIVE)
# ================================================================
# The signed curvature flux is the 3D production flux with F → ηFη
# (η = diag(1,1,1,−1), time = matrix index 3) — the 5a §10d ℋ sign rule:
# spatial matrix pairs positive, (α,3) pairs negative. Sandbox anchors:
# m5_8_2a (fuel) → 2b-1 (CC) → 2b-2 (field clock) → 2c-1 (full nonlinear,
# f64 numpy — THE cross-validation reference for these f32 kernels).
# Ghost guard: the per-voxel stable_mask (K>0 ∧ Q(∇ψ) positive-definite,
# computed ONCE post-seed) selects signed flux on the stable region and the
# always-stable Euclidean flux on the fuel shell (the 2b-2 λ≈15.6/t linear
# runaway sector — its propulsion physics needs the constrained faithful
# kernel, future work; v1 keeps it Euclidean-stable). Global guard: the
# coherent (α,3) velocity drift is sampled on 3 mid-planes (Metal-safe —
# NO full-grid reductions, the atomics lesson) and subtracted in evolve.


@ti.func
def eta_twist_masked(f, stable: ti.f32):  # type: ignore
    """F → ηFη blended by the stable mask: (α,3)/(3,α) comps × (1 − 2·stable).

    stable=1 → the Minkowski-signed twist; stable=0 → identity (Euclidean
    fallback on the fuel shell)."""
    s = 1.0 - 2.0 * stable
    out = f
    for a_ in ti.static(range(3)):
        out[a_, 3] = f[a_, 3] * s
        out[3, a_] = f[3, a_] * s
    return out


@ti.func
def signed_dot4(a4, b4):  # type: ignore
    """⟨A,B⟩_s = Σ A∘(ηBη): spatial-pair comps +, (α,3) comps − (full-matrix sum)."""
    acc = 0.0
    for p_ in ti.static(range(4)):
        for q_ in ti.static(range(4)):
            sgn = 1.0
            if ti.static((p_ == 3) != (q_ == 3)):
                sgn = -1.0
            acc += a4[p_, q_] * b4[p_, q_] * sgn
    return acc


@ti.kernel
def compute_stable_mask(tensor_field: ti.template()):  # type: ignore
    """Per-voxel ghost guard (run ONCE post-seed): stable_mask = 1 where the
    signed clock kernel is well-posed — K(x) > 0 AND the 3×3 gradient-stiffness
    form Q(∇ψ) is positive-definite (the 2b-2 criterion, exact):

        P_i = [M_i, M_ψ] ,  p_ij = ⟨P_i, P_j⟩_s ,  K = 4 Σ_i p_ii
        Q_ii = Σ_(j≠i) p_jj ,  Q_ij = −p_ij  (i≠j) ;  stable ⇔ K>0 ∧ minEig(Q)>0

    minEig via the analytic Cardano solver (NOT ti.sym_eig — the Metal/f32
    lesson). Boundary voxels default to 0 (Euclidean). Seed-time cost only.
    """
    nx, ny, nz = tensor_field.nx, tensor_field.ny, tensor_field.nz
    inv_2dx = 1.0 / (2.0 * tensor_field.dx_am)
    for i, j, k in ti.ndrange(nx, ny, nz):
        tensor_field.stable_mask[i, j, k] = 0.0
    for i, j, k in ti.ndrange((1, nx - 1), (1, ny - 1), (1, nz - 1)):
        m_psi = tensor_field.M_psi_am[i, j, k]
        mx = (tensor_field.M_am[i + 1, j, k] - tensor_field.M_am[i - 1, j, k]) * inv_2dx
        my = (tensor_field.M_am[i, j + 1, k] - tensor_field.M_am[i, j - 1, k]) * inv_2dx
        mz = (tensor_field.M_am[i, j, k + 1] - tensor_field.M_am[i, j, k - 1]) * inv_2dx
        p0 = commutator(mx, m_psi)
        p1 = commutator(my, m_psi)
        p2 = commutator(mz, m_psi)
        p00 = signed_dot4(p0, p0)
        p11 = signed_dot4(p1, p1)
        p22 = signed_dot4(p2, p2)
        p01 = signed_dot4(p0, p1)
        p02 = signed_dot4(p0, p2)
        p12 = signed_dot4(p1, p2)
        kin = 4.0 * (p00 + p11 + p22)
        q = ti.Matrix([[p11 + p22, -p01, -p02], [-p01, p00 + p22, -p12], [-p02, -p12, p00 + p11]])
        _, eigs = principal_director(q)  # sorted λ₁≥λ₂≥λ₃
        qscale = ti.abs(q[0, 0]) + ti.abs(q[1, 1]) + ti.abs(q[2, 2]) + 1e-30
        if kin > 0.0 and eigs[2] > 1e-6 * qscale:
            tensor_field.stable_mask[i, j, k] = 1.0


@ti.kernel
def compute_tstar(tensor_field: ti.template()):  # type: ignore
    """Seed the per-voxel V(M) amplitude target t*(x) = Tr(M_sp²) of the CURRENT
    field (call once, right after seeding): the dressed seed becomes the exact
    V-equilibrium — zero static V-force (the 2c-0 'pin the well to the DRESSED
    amplitude' design input; fixes the 2c-2 slow V-pump on the off-minimum
    dressed shell). The well then confines amplitude EXCURSIONS, not the seed."""
    nx, ny, nz = tensor_field.nx, tensor_field.ny, tensor_field.nz
    for i, j, k in ti.ndrange(nx, ny, nz):
        m = tensor_field.M_am[i, j, k]
        msp = ti.Matrix(
            [[m[0, 0], m[0, 1], m[0, 2]], [m[1, 0], m[1, 1], m[1, 2]], [m[2, 0], m[2, 1], m[2, 2]]]
        )
        tensor_field.ldg_tstar[i, j, k] = (msp @ msp).trace()


@ti.func
def dV_M_dressed(m, cc: ti.f32, tstar: ti.f32):  # type: ignore
    """∂V/∂M for the DRESSED well V = cc·(Tr(M_sp²) − t*(x))²: 4·cc·(t−t*)·M_sp,
    spatial block only (time row/col force = 0 — the M5.8.1 rule)."""
    msp = ti.Matrix(
        [[m[0, 0], m[0, 1], m[0, 2]], [m[1, 0], m[1, 1], m[1, 2]], [m[2, 0], m[2, 1], m[2, 2]]]
    )
    t = (msp @ msp).trace()
    dsp = 4.0 * cc * (t - tstar) * msp
    d4 = ti.Matrix.zero(ti.f32, 4, 4)
    for i in ti.static(range(3)):
        for j in ti.static(range(3)):
            d4[i, j] = dsp[i, j]
    return d4


@ti.kernel
def compute_curvature_flux_4d(tensor_field: ti.template()):  # type: ignore
    """The 4D signed curvature flux G_α = 8 Σ_ν [tw(F_αν), M_ν], tw = ηFη on the
    stable region / identity on the fuel shell (stable_mask blend). At b=0 the
    seed has no (α,3) components ⇒ tw is a no-op ⇒ EXACTLY compute_curvature_flux
    (the headless identity check). Same two-pass contract as the 3D kernel."""
    nx, ny, nz = tensor_field.nx, tensor_field.ny, tensor_field.nz
    inv_2dx = 1.0 / (2.0 * tensor_field.dx_am)
    for i, j, k in ti.ndrange((1, nx - 1), (1, ny - 1), (1, nz - 1)):
        st = tensor_field.stable_mask[i, j, k]
        mx = (tensor_field.M_am[i + 1, j, k] - tensor_field.M_am[i - 1, j, k]) * inv_2dx
        my = (tensor_field.M_am[i, j + 1, k] - tensor_field.M_am[i, j - 1, k]) * inv_2dx
        mz = (tensor_field.M_am[i, j, k + 1] - tensor_field.M_am[i, j, k - 1]) * inv_2dx
        fxy = eta_twist_masked(commutator(mx, my), st)
        fxz = eta_twist_masked(commutator(mx, mz), st)
        fyz = eta_twist_masked(commutator(my, mz), st)
        tensor_field.curv_flux_x[i, j, k] = 8.0 * (commutator(fxy, my) + commutator(fxz, mz))
        tensor_field.curv_flux_y[i, j, k] = 8.0 * (
            commutator(-1.0 * fxy, mx) + commutator(fyz, mz)
        )
        tensor_field.curv_flux_z[i, j, k] = 8.0 * (
            commutator(-1.0 * fxz, mx) + commutator(-1.0 * fyz, my)
        )


@ti.kernel
def sample_v03_drift(tensor_field: ti.template(), dt_rs: ti.f32):  # type: ignore
    """Sample the coherent (α,3) velocity drift on 3 mid-planes (Metal-safe:
    plane-scale atomics only — the full-grid-reduction lesson). Results land in
    tensor_field.v03_sums[0..2]; the caller divides by the plane-voxel count."""
    nx, ny, nz = tensor_field.nx, tensor_field.ny, tensor_field.nz
    for a_ in ti.static(range(3)):
        tensor_field.v03_sums[a_] = 0.0
    inv_dt = 1.0 / dt_rs
    im, jm, km = nx // 2, ny // 2, nz // 2
    for j, k in ti.ndrange(ny, nz):
        d = (tensor_field.M_am[im, j, k] - tensor_field.M_prev_am[im, j, k]) * inv_dt
        for a_ in ti.static(range(3)):
            tensor_field.v03_sums[a_] += d[a_, 3]
    for i, k in ti.ndrange(nx, nz):
        d = (tensor_field.M_am[i, jm, k] - tensor_field.M_prev_am[i, jm, k]) * inv_dt
        for a_ in ti.static(range(3)):
            tensor_field.v03_sums[a_] += d[a_, 3]
    for i, j in ti.ndrange(nx, ny):
        d = (tensor_field.M_am[i, j, km] - tensor_field.M_prev_am[i, j, km]) * inv_dt
        for a_ in ti.static(range(3)):
            tensor_field.v03_sums[a_] += d[a_, 3]


@ti.kernel
def evolve_M_4d(
    tensor_field: ti.template(),  # type: ignore
    c_amrs: ti.f32,  # type: ignore
    dt_rs: ti.f32,  # type: ignore
    ldg_cc: ti.f32,  # type: ignore
    vm0: ti.f32,  # type: ignore
    vm1: ti.f32,  # type: ignore
    vm2: ti.f32,  # type: ignore
    km: ti.f32,  # type: ignore
):
    """The M5.8.2 leapfrog: the M5.8.1 time-freeze clamp is REPLACED by the soft
    global guard — the time axis is LIVE. Identical to evolve_M except:

    (i) no (α,3) zeroing and no M[3,3] pin (boost dressing + clock evolve freely);
    (ii) the sampled coherent (α,3) drift (vm0..2, from sample_v03_drift) is
        subtracted — the 2b-1 ghost channel (a free GLOBAL dressing mode);
    (iii) DIAGONAL FAITHFUL-LITE INERTIA m(x) = 1 + km·dx²·Σ_i‖M_i‖²_F dividing
        the acceleration — the scalar shadow of the faithful kinetic operator
        A(Ṁ) = 4Σ[η[Ṁ,M_i]η,M_i] that stabilized the 2c-1 spike: inertia grows
        with local structure (steepening modes get HEAVY and decelerate — the
        self-regulation the simple kinetic lacks, which let the signed
        stiffness pump the core, the 2c-2 repro). m → 1 in weak-gradient
        regions ⇒ the 3D/vacuum limit is untouched; km=0 disables.

    V(M) uses the DRESSED well dV_M_dressed with the per-voxel target t*(x)
    (compute_tstar, run once post-seed): the seed is the exact V-equilibrium —
    zero static V-force (the 2c-2 fix for the slow V-pump on the off-minimum
    dressed shell; 2c-0 energetics unchanged: confinement of excursions kept)."""
    nx, ny, nz = tensor_field.nx, tensor_field.ny, tensor_field.nz
    inv_2dx = 1.0 / (2.0 * tensor_field.dx_am)
    c2 = c_amrs * c_amrs
    dt2 = dt_rs * dt_rs
    dx2 = tensor_field.dx_am * tensor_field.dx_am
    for i, j, k in ti.ndrange((2, nx - 2), (2, ny - 2), (2, nz - 2)):
        div_G = (
            (tensor_field.curv_flux_x[i + 1, j, k] - tensor_field.curv_flux_x[i - 1, j, k])
            + (tensor_field.curv_flux_y[i, j + 1, k] - tensor_field.curv_flux_y[i, j - 1, k])
            + (tensor_field.curv_flux_z[i, j, k + 1] - tensor_field.curv_flux_z[i, j, k - 1])
        ) * inv_2dx
        force = c2 * div_G - dV_M_dressed(
            tensor_field.M_am[i, j, k], ldg_cc, tensor_field.ldg_tstar[i, j, k]
        )
        mx = (tensor_field.M_am[i + 1, j, k] - tensor_field.M_am[i - 1, j, k]) * inv_2dx
        my = (tensor_field.M_am[i, j + 1, k] - tensor_field.M_am[i, j - 1, k]) * inv_2dx
        mz = (tensor_field.M_am[i, j, k + 1] - tensor_field.M_am[i, j, k - 1]) * inv_2dx
        mloc = 1.0 + km * dx2 * (mx.norm_sqr() + my.norm_sqr() + mz.norm_sqr())
        m_new = (
            2.0 * tensor_field.M_am[i, j, k]
            - tensor_field.M_prev_am[i, j, k]
            + dt2 * force * (1.0 / mloc)
        )
        # guard (a): remove the coherent global (α,3) drift accumulated this step
        m_new[0, 3] -= vm0 * dt_rs
        m_new[3, 0] -= vm0 * dt_rs
        m_new[1, 3] -= vm1 * dt_rs
        m_new[3, 1] -= vm1 * dt_rs
        m_new[2, 3] -= vm2 * dt_rs
        m_new[3, 2] -= vm2 * dt_rs
        tensor_field.M_new_am[i, j, k] = m_new


# ================================================================
# M5.8.2c OPTION B — THE CONSTRAINED SPECTRAL-PROJECTION INTEGRATOR
# ================================================================
# RESTORED 2026-06-05 late evening — this section was erased from the working
# tree by a stale-IDE-buffer save (the file was open with a pre-Option-B
# buffer; the GUI "even worse" run executed WITHOUT these kernels). Content
# identical to the B-2 wiring + the velocity-kick fix.
#
# The Minkowski-SIGNED dynamics, stable (INTEGRATOR_4D: "constrained"). The
# structural finding (2c-2): the signed flux under ANY cheap positive inertia
# has slow growing modes — only the faithful kinetic works. Per voxel, the
# inertia operator A(Ṁ) = 4 Σ_i [η[Ṁ,M_i]η, M_i] is eigendecomposed in the
# 10-dim orthonormal symmetric-matrix basis (OWN cyclic Jacobi on local
# matrices — NOT ti.sym_eig, the Metal/f32 lesson); only positive-inertia
# directions evolve (λ > EPS·max|λ|), and P is PROJECTED onto the kept
# subspace every step (frozen directions must not accumulate momentum — the
# 2c-1 v2 constraint-switching pump). State is (M, P): symplectic-Euler
# P += dt·force → clamp → solve Ṁ = A⁺P → M += dt·Ṁ. Validated machine-exact
# against the f64 numpy reference (sandbox_v8/m5_8_2cb_taichi_constrained.py:
# field diff 5.5e-15 @890 steps; f32 1e-5 benign; Metal clean; 8.3 ms/step
# solve at 64³ — B-1 gates [J1][B1][B2][B2m][B3][B4] all PASS, 2026-06-05).
# ⚠️ KNOWN OPEN (2cb-2 long-horizon): the f32 kernel develops a negative-H
# ghost runaway past ~1100 steps even in the validated 24³ config (H plateau
# → collapse; align 1.0 → 0.55 by 6000) — the f64-vs-f32 6000-step
# scheme-vs-precision discrimination is the active investigation.
#
# UNITS + NORMALIZATION: these kernels run the 2c-1 natural convention
# (c = 1, τ = c·t) — the caller passes dt_eff = c_amrs·dt_rs (am). Fluxes and
# A both use the 2c-1 4× normalization (the production 8× leapfrog flux is 2×
# this; a COMMON factor cancels exactly in Ṁ = A⁺P, so the dynamics here is
# the validated one verbatim). V enters the τ-units force as dV/(2c²) — the
# launcher passes ldg_cc_4d = LDG_STIFFNESS_K·0.5/dx⁴ (the c² cancels). The
# constrained path is FULLY signed (no stable-mask blend — the spectral
# projection IS the exact pointwise ghost treatment) and ignores
# KM_INERTIA_4D (the faithful A replaces the diagonal shadow).

EPS_EIG_4D = 0.05  # positive-inertia keep threshold (× max|λ| per voxel)
VCAP_4D = 5.0  # ‖Ṁ‖_F backstop (never engaged in the validated 900-step runs)
JACOBI_SWEEPS = 20  # max cyclic sweeps (converges in ~6–8)
JTOL2_4D = (2e-7) ** 2  # relative off-diagonal² stop (f32 eps scale)


@ti.func
def eta_twist(f):  # type: ignore
    """F → ηFη, unmasked — the constrained path is fully Minkowski-signed."""
    out = f
    for a_ in ti.static(range(3)):
        out[a_, 3] = -f[a_, 3]
        out[3, a_] = -f[3, a_]
    return out


@ti.kernel
def init_P_4d(tensor_field: ti.template(), kick: ti.f32):  # type: ignore
    """P₀ = A(Ṁ₀) with Ṁ₀ = kick·M_ψ — the EXACT 2c-1 kick semantics (a
    VELOCITY in τ-units, dt-INDEPENDENT; kick = 0 ⇒ P₀ = 0 exactly). NOT the
    buffer encoding (M − M_prev)/dt: that made the kick velocity scale as
    1/dt — 16× hotter than validated at the constrained dt (the step-31
    v-cap engagement in the 2cb-2 repro). Run ONCE post-seed. Md_am gets Ṁ₀
    on the act region (read by the first flux call — the 2c-1 step-0
    semantics)."""
    nx, ny, nz = tensor_field.nx, tensor_field.ny, tensor_field.nz
    inv_2dx = 1.0 / (2.0 * tensor_field.dx_am)
    for i, j, k in ti.ndrange(nx, ny, nz):
        tensor_field.P_am[i, j, k] = ti.Matrix.zero(ti.f32, 4, 4)
        tensor_field.Md_am[i, j, k] = ti.Matrix.zero(ti.f32, 4, 4)
    for i, j, k in ti.ndrange((1, nx - 1), (1, ny - 1), (1, nz - 1)):
        md = kick * tensor_field.M_psi_am[i, j, k]
        mx = (tensor_field.M_am[i + 1, j, k] - tensor_field.M_am[i - 1, j, k]) * inv_2dx
        my = (tensor_field.M_am[i, j + 1, k] - tensor_field.M_am[i, j - 1, k]) * inv_2dx
        mz = (tensor_field.M_am[i, j, k + 1] - tensor_field.M_am[i, j, k - 1]) * inv_2dx
        a4 = 4.0 * (
            commutator(eta_twist(commutator(md, mx)), mx)
            + commutator(eta_twist(commutator(md, my)), my)
            + commutator(eta_twist(commutator(md, mz)), mz)
        )
        tensor_field.P_am[i, j, k] = a4
        tensor_field.Md_am[i, j, k] = md * tensor_field.act4d[i, j, k]


@ti.kernel
def flux_4d_constrained(tensor_field: ti.template()):  # type: ignore
    """The 2c-1 flux (C1-pinned signs): G_i = ∂U/∂M_i − ∂T/∂M_i with
    dU_i = 4 Σ_{j≠i} [η F̃_ij η, M_j] and dT_i = −4 [η F_0i η, Ṁ] — fully
    signed, NO mask. Writes curv_flux_x/y/z (border stays zero)."""
    nx, ny, nz = tensor_field.nx, tensor_field.ny, tensor_field.nz
    inv_2dx = 1.0 / (2.0 * tensor_field.dx_am)
    for i, j, k in ti.ndrange((1, nx - 1), (1, ny - 1), (1, nz - 1)):
        md = tensor_field.Md_am[i, j, k]
        mx = (tensor_field.M_am[i + 1, j, k] - tensor_field.M_am[i - 1, j, k]) * inv_2dx
        my = (tensor_field.M_am[i, j + 1, k] - tensor_field.M_am[i, j - 1, k]) * inv_2dx
        mz = (tensor_field.M_am[i, j, k + 1] - tensor_field.M_am[i, j, k - 1]) * inv_2dx
        fxy = commutator(mx, my)
        fxz = commutator(mx, mz)
        fyz = commutator(my, mz)
        tensor_field.curv_flux_x[i, j, k] = 4.0 * (
            commutator(eta_twist(fxy), my)
            + commutator(eta_twist(fxz), mz)
            + commutator(eta_twist(commutator(md, mx)), md)
        )
        tensor_field.curv_flux_y[i, j, k] = 4.0 * (
            commutator(eta_twist(-1.0 * fxy), mx)
            + commutator(eta_twist(fyz), mz)
            + commutator(eta_twist(commutator(md, my)), md)
        )
        tensor_field.curv_flux_z[i, j, k] = 4.0 * (
            commutator(eta_twist(-1.0 * fxz), mx)
            + commutator(eta_twist(-1.0 * fyz), my)
            + commutator(eta_twist(commutator(md, mz)), md)
        )


@ti.kernel
def update_P_4d(tensor_field: ti.template(), dt_eff: ti.f32, ldg_cc4d: ti.f32):  # type: ignore
    """P += dt_eff·(div G − dV_dressed) — the τ-units force step (ldg_cc4d
    pre-scaled by the launcher to the 4× convention: K·0.5/dx⁴)."""
    nx, ny, nz = tensor_field.nx, tensor_field.ny, tensor_field.nz
    inv_2dx = 1.0 / (2.0 * tensor_field.dx_am)
    for i, j, k in ti.ndrange((1, nx - 1), (1, ny - 1), (1, nz - 1)):
        div_g = (
            (tensor_field.curv_flux_x[i + 1, j, k] - tensor_field.curv_flux_x[i - 1, j, k])
            + (tensor_field.curv_flux_y[i, j + 1, k] - tensor_field.curv_flux_y[i, j - 1, k])
            + (tensor_field.curv_flux_z[i, j, k + 1] - tensor_field.curv_flux_z[i, j, k - 1])
        ) * inv_2dx
        force = div_g - dV_M_dressed(
            tensor_field.M_am[i, j, k], ldg_cc4d, tensor_field.ldg_tstar[i, j, k]
        )
        tensor_field.P_am[i, j, k] += dt_eff * force


@ti.kernel
def sample_p03_drift(tensor_field: ti.template()):  # type: ignore
    """3-mid-plane act-region sums of the (α,3) MOMENTUM components → v03_sums
    (the 2c-1 global clamp, plane-sampled — Metal-safe, NO full-grid atomics).
    The caller divides by the act-plane voxel count (launcher-precomputed)."""
    nx, ny, nz = tensor_field.nx, tensor_field.ny, tensor_field.nz
    for a_ in ti.static(range(3)):
        tensor_field.v03_sums[a_] = 0.0
    im, jm, km = nx // 2, ny // 2, nz // 2
    for j, k in ti.ndrange(ny, nz):
        if tensor_field.act4d[im, j, k] > 0.5:
            p = tensor_field.P_am[im, j, k]
            for a_ in ti.static(range(3)):
                tensor_field.v03_sums[a_] += p[a_, 3]
    for i, k in ti.ndrange(nx, nz):
        if tensor_field.act4d[i, jm, k] > 0.5:
            p = tensor_field.P_am[i, jm, k]
            for a_ in ti.static(range(3)):
                tensor_field.v03_sums[a_] += p[a_, 3]
    for i, j in ti.ndrange(nx, ny):
        if tensor_field.act4d[i, j, km] > 0.5:
            p = tensor_field.P_am[i, j, km]
            for a_ in ti.static(range(3)):
                tensor_field.v03_sums[a_] += p[a_, 3]


@ti.kernel
def apply_p03_clamp(tensor_field: ti.template(), m0: ti.f32, m1: ti.f32, m2: ti.f32):  # type: ignore
    """Guard (a): subtract the act-mean (α,3) momentum (act region only), then
    restore P symmetry everywhere — the exact 2c-1 clamp order."""
    nx, ny, nz = tensor_field.nx, tensor_field.ny, tensor_field.nz
    for i, j, k in ti.ndrange(nx, ny, nz):
        if tensor_field.act4d[i, j, k] > 0.5:
            tensor_field.P_am[i, j, k][0, 3] -= m0
            tensor_field.P_am[i, j, k][1, 3] -= m1
            tensor_field.P_am[i, j, k][2, 3] -= m2
        tensor_field.P_am[i, j, k][3, 0] = tensor_field.P_am[i, j, k][0, 3]
        tensor_field.P_am[i, j, k][3, 1] = tensor_field.P_am[i, j, k][1, 3]
        tensor_field.P_am[i, j, k][3, 2] = tensor_field.P_am[i, j, k][2, 3]


@ti.kernel
def solve_constrained_4d(tensor_field: ti.template()):  # type: ignore
    """THE constrained solve (the B-1-validated kernel), fused per voxel:
    build the 10×10 inertia operator A in the sym basis → cyclic Jacobi on
    LOCAL matrices → positive-inertia keep → project P AND solve Ṁ = A⁺P on
    the kept subspace → act mask + ‖Ṁ‖ backstop. Zero-gradient voxels (grid
    border) give A = 0 ⇒ all directions frozen ⇒ P, Ṁ → 0 (the 2c-1 limit)."""
    nx, ny, nz = tensor_field.nx, tensor_field.ny, tensor_field.nz
    inv_2dx = 1.0 / (2.0 * tensor_field.dx_am)
    for i, j, k in ti.ndrange(nx, ny, nz):
        mx = ti.Matrix.zero(ti.f32, 4, 4)
        my = ti.Matrix.zero(ti.f32, 4, 4)
        mz = ti.Matrix.zero(ti.f32, 4, 4)
        if 0 < i < nx - 1:
            mx = (tensor_field.M_am[i + 1, j, k] - tensor_field.M_am[i - 1, j, k]) * inv_2dx
        if 0 < j < ny - 1:
            my = (tensor_field.M_am[i, j + 1, k] - tensor_field.M_am[i, j - 1, k]) * inv_2dx
        if 0 < k < nz - 1:
            mz = (tensor_field.M_am[i, j, k + 1] - tensor_field.M_am[i, j, k - 1]) * inv_2dx
        # --- build A[kk,li] = ⟨A_op(E_li), E_kk⟩, then symmetrize ------------
        aloc = ti.Matrix.zero(ti.f32, 10, 10)
        for li in range(10):
            el = tensor_field.sym_basis[li]
            ae = 4.0 * (
                commutator(eta_twist(commutator(el, mx)), mx)
                + commutator(eta_twist(commutator(el, my)), my)
                + commutator(eta_twist(commutator(el, mz)), mz)
            )
            for kk in range(10):
                aloc[kk, li] = (ae * tensor_field.sym_basis[kk]).sum()
        for p_ in range(10):
            for q_ in range(p_ + 1, 10):
                v = 0.5 * (aloc[p_, q_] + aloc[q_, p_])
                aloc[p_, q_] = v
                aloc[q_, p_] = v
        # --- cyclic Jacobi: A ← JᵀAJ, Q ← QJ ---------------------------------
        qloc = ti.Matrix.identity(ti.f32, 10)
        normf2 = (aloc * aloc).sum()
        off2 = normf2
        sweep = 0
        while sweep < JACOBI_SWEEPS and off2 > JTOL2_4D * normf2:
            for p_ in range(9):
                for q_ in range(p_ + 1, 10):
                    apq = aloc[p_, q_]
                    if ti.abs(apq) > 0:
                        tau = (aloc[q_, q_] - aloc[p_, p_]) / (2.0 * apq)
                        t = 1.0 / (ti.abs(tau) + ti.sqrt(1.0 + tau * tau))
                        if tau < 0:
                            t = -t
                        c = 1.0 / ti.sqrt(1.0 + t * t)
                        s = t * c
                        for r_ in range(10):  # A ← A·J (cols p,q)
                            arp = aloc[r_, p_]
                            arq = aloc[r_, q_]
                            aloc[r_, p_] = c * arp - s * arq
                            aloc[r_, q_] = s * arp + c * arq
                        for r_ in range(10):  # A ← Jᵀ·A (rows p,q)
                            apr = aloc[p_, r_]
                            aqr = aloc[q_, r_]
                            aloc[p_, r_] = c * apr - s * aqr
                            aloc[q_, r_] = s * apr + c * aqr
                        for r_ in range(10):  # Q ← Q·J
                            qrp = qloc[r_, p_]
                            qrq = qloc[r_, q_]
                            qloc[r_, p_] = c * qrp - s * qrq
                            qloc[r_, q_] = s * qrp + c * qrq
            off2 = 0.0
            for p_ in range(9):
                for q_ in range(p_ + 1, 10):
                    off2 += 2.0 * aloc[p_, q_] * aloc[p_, q_]
            sweep += 1
        # --- keep mask + projection + A⁺ solve --------------------------------
        lmax = 0.0
        for a_ in range(10):
            lmax = ti.max(lmax, ti.abs(aloc[a_, a_]))
        pm = tensor_field.P_am[i, j, k]
        pc = ti.Vector.zero(ti.f32, 10)
        for a_ in range(10):
            pc[a_] = (pm * tensor_field.sym_basis[a_]).sum()
        pproj = ti.Vector.zero(ti.f32, 10)
        cdot = ti.Vector.zero(ti.f32, 10)
        for kk in range(10):
            lamk = aloc[kk, kk]
            qp = 0.0
            for a_ in range(10):
                qp += qloc[a_, kk] * pc[a_]
            if lamk > EPS_EIG_4D * (lmax + 1e-30):
                for a_ in range(10):
                    pproj[a_] += qloc[a_, kk] * qp
                    cdot[a_] += qloc[a_, kk] * qp / lamk
        pnew = ti.Matrix.zero(ti.f32, 4, 4)
        mdnew = ti.Matrix.zero(ti.f32, 4, 4)
        for a_ in range(10):
            pnew += pproj[a_] * tensor_field.sym_basis[a_]
            mdnew += cdot[a_] * tensor_field.sym_basis[a_]
        mdnew = mdnew * tensor_field.act4d[i, j, k]
        vn = ti.sqrt((mdnew * mdnew).sum())
        if vn > VCAP_4D:
            mdnew = mdnew * (VCAP_4D / (vn + 1e-30))
        tensor_field.P_am[i, j, k] = pnew
        tensor_field.Md_am[i, j, k] = mdnew


@ti.kernel
def update_M_4d_constrained(tensor_field: ti.template(), dt_eff: ti.f32):  # type: ignore
    """M_new = M + dt_eff·Ṁ (all voxels; Ṁ = 0 outside act ⇒ border/core copy).
    The caller swaps the triple buffer as usual — render paths read M_am."""
    nx, ny, nz = tensor_field.nx, tensor_field.ny, tensor_field.nz
    for i, j, k in ti.ndrange(nx, ny, nz):
        tensor_field.M_new_am[i, j, k] = (
            tensor_field.M_am[i, j, k] + dt_eff * tensor_field.Md_am[i, j, k]
        )


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
# headroom; the launcher computes this from tensor_field.dx_am at runtime.
#
# Output buffer convention: writes director_nhat_new. The caller swaps it into
# director_nhat, then rebuild_M_from_director writes all three M buffers equal so
# subsequent Evolve PDE sees Ṁ = 0 (no spurious time-derivative from the relaxation
# update — relaxation is a STATIC operation, not a dynamics step). This is different
# from swap_matrix_buffers, which rotates the M buffers for time evolution.
#
# Boundary handling: director_nhat_new boundary voxels are copied from director_nhat
# (no update) — preserves the fixed-value Dirichlet BC (e.g. ẑ at the universe edge).
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
    tensor_field: ti.template(),  # type: ignore
    tau: ti.f32,  # type: ignore
    pin_centers: ti.template(),  # type: ignore
    pin_signs: ti.template(),  # type: ignore
    n_defects: ti.i32,  # type: ignore
):
    """
    One gradient-descent step on the Frank elastic energy `H_F = (K/2)·|∇n̂|²`.

    Reads:  tensor_field.director_nhat
    Writes: tensor_field.director_nhat_new  (caller copies → director_nhat)

    M5.4: repointed from the retiring ψ buffers onto the matrix-substrate director
    (director_nhat = principal eigenvector of M). The math is byte-identical to the
    M5.1 relaxation — tangent-projected gradient descent on (K/2)|∇n̂|² + soft-core
    pin — so the Coulomb result carries over. Caller rebuilds M from the relaxed
    director via rebuild_M_from_director.

    Args:
        tensor_field: TensorField with the director field in director_nhat
        tau: step size (must satisfy τ < dx²/(6·K) for stability)
        pin_centers: ti.field shape (n_defects, 3) i32 — defect centers in voxel coords
        pin_signs: ti.field shape (n_defects,) i32 — ±1 per defect (pin direction)
        n_defects: number of active defects in pin_centers/signs

    Note: K_frank is implicit in `tau`'s CFL derivation — the kernel itself
    doesn't multiply by K because gradient descent on `(K/2)·|∇n̂|²` gives
    `∂n/∂τ = K·∇²n` and the K can be absorbed into the step size (τ_eff = K·τ).
    For K=1 (M5.1 default) this is τ_eff = τ.
    """
    nx, ny, nz = tensor_field.nx, tensor_field.ny, tensor_field.nz

    # ---- Interior: tangent-projected gradient step + renormalization ----
    for i, j, k in ti.ndrange((1, nx - 1), (1, ny - 1), (1, nz - 1)):
        n = tensor_field.director_nhat[i, j, k]
        lap = compute_laplacian_director(tensor_field, i, j, k)
        n_dot_lap = n.dot(lap)
        dn = lap - n_dot_lap * n  # tangent projection (remove component along n)
        n_new = n + tau * dn
        norm = n_new.norm() + 1e-12
        tensor_field.director_nhat_new[i, j, k] = n_new / norm

    # ---- Boundary: preserve Dirichlet BC by copying from current director ----
    for j, k in ti.ndrange(ny, nz):
        tensor_field.director_nhat_new[0, j, k] = tensor_field.director_nhat[0, j, k]
        tensor_field.director_nhat_new[nx - 1, j, k] = tensor_field.director_nhat[nx - 1, j, k]
    for i, k in ti.ndrange(nx, nz):
        tensor_field.director_nhat_new[i, 0, k] = tensor_field.director_nhat[i, 0, k]
        tensor_field.director_nhat_new[i, ny - 1, k] = tensor_field.director_nhat[i, ny - 1, k]
    for i, j in ti.ndrange(nx, ny):
        tensor_field.director_nhat_new[i, j, 0] = tensor_field.director_nhat[i, j, 0]
        tensor_field.director_nhat_new[i, j, nz - 1] = tensor_field.director_nhat[i, j, nz - 1]

    # ---- Pin defect cores (soft Dirichlet at single closest voxel) ----
    for d in range(n_defects):
        ci = pin_centers[d, 0]
        cj = pin_centers[d, 1]
        ck = pin_centers[d, 2]
        sgn = ti.cast(pin_signs[d], ti.f32)
        tensor_field.director_nhat_new[ci, cj, ck] = ti.Vector([0.0, 0.0, sgn])
