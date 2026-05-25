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
# ψ EVOLVING ENGINE — LAGRANGIAN FIELD EVOLUTION
# ================================================================
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

    Reads:  wave_field.psi_am
    Writes: wave_field.psi_new_am  (caller copies → psi_am AND psi_prev_am)

    Args:
        wave_field: WaveField with the director field in psi_am
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
        n = wave_field.psi_am[i, j, k]
        lap = compute_laplacian(wave_field, i, j, k)
        n_dot_lap = n.dot(lap)
        dn = lap - n_dot_lap * n  # tangent projection (remove component along n)
        n_new = n + tau * dn
        norm = n_new.norm() + 1e-12
        wave_field.psi_new_am[i, j, k] = n_new / norm

    # ---- Boundary: preserve Dirichlet BC by copying from current ψ ----
    for j, k in ti.ndrange(ny, nz):
        wave_field.psi_new_am[0, j, k] = wave_field.psi_am[0, j, k]
        wave_field.psi_new_am[nx - 1, j, k] = wave_field.psi_am[nx - 1, j, k]
    for i, k in ti.ndrange(nx, nz):
        wave_field.psi_new_am[i, 0, k] = wave_field.psi_am[i, 0, k]
        wave_field.psi_new_am[i, ny - 1, k] = wave_field.psi_am[i, ny - 1, k]
    for i, j in ti.ndrange(nx, ny):
        wave_field.psi_new_am[i, j, 0] = wave_field.psi_am[i, j, 0]
        wave_field.psi_new_am[i, j, nz - 1] = wave_field.psi_am[i, j, nz - 1]

    # ---- Pin defect cores (soft Dirichlet at single closest voxel) ----
    for d in range(n_defects):
        ci = pin_centers[d, 0]
        cj = pin_centers[d, 1]
        ck = pin_centers[d, 2]
        sgn = ti.cast(pin_signs[d], ti.f32)
        wave_field.psi_new_am[ci, cj, ck] = ti.Vector([0.0, 0.0, sgn])


