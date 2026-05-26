"""
M5 ENGINE — FIELD OBSERVABLES (engine3_observables)

Derived measurements computed from the field each frame:
  - update_trackers          per-voxel amplitude/frequency EMA (Trackers)
  - compute_energyH_density  Hamiltonian energy density (FieldObservables)
  - compute_energyF_density  Frank elastic energy density
  - compute_winding_number   topological charge Q on a sphere (CPU numpy)
  - sample_avg_* + helpers   3-plane global aggregates

Depends on engine2_pde (V_psi, used by compute_energyH_density).
"""

import taichi as ti
import numpy as np

from . import engine2_pde as pde

# ================================================================
# TRACKERS — PER-VOXEL AMPLITUDE & FREQUENCY
# ================================================================
# Trackers run after evolve_psi each step. They derive observable scalars
# from ψ (the field) without touching it — read-only on psi_am, write to the
# tracker fields in-place. EMA / zero-crossing logic preserved from M4 but
# adapted to vector ψ:
#   - amp_local_emarms_am  : RMS of |ψ| over recent history (EMA on |ψ|²)
#   - freq_local_cross_rHz : zero-crossing rate of ψ_z (chosen polarization axis)
#
# PERFORMANCE NOTE (deferred to M5.0i): this is a SECOND full-grid pass over
# ψ — evolve_psi already touched every voxel; we re-stream ψ and ψ_prev
# from memory just to compute trackers. At 100M voxels this costs ~5 GB of
# extra memory bandwidth per step (~25–30% of total per-step traffic).
#
# M2 avoided this by computing trackers INSIDE the same loop as the wave
# update, but M2's wave kernel was analytical (no neighbor dependency), so
# trackers could read the just-written ψ trivially. M5's leapfrog can do the
# same: in a merged kernel, write ψ_new then immediately use it (with
# ψ_prev already in registers) for tracker EMA / zero-crossing — gives true
# central-difference ψ̇ for free, more accurate than the forward diff used
# below (psi_new gets consumed by swap_buffers before update_trackers
# can see it, so we fall back to (ψ − ψ_prev) / dt here).
#
# Splitting was deliberate for M5.0d: cleaner kernel boundaries, AMR-ready
# (octree retrofit in M5.6/M5.8 prefers single-purpose kernels), and trackers
# are skippable in tests (e.g., M5.0h dispersion regression). Profile-first
# in M5.0i — alongside swap_buffers' rotating-pointer optimization (see
# medium.py:swap_buffers docstring), this is the second-largest known
# bandwidth target on the per-step path.


@ti.kernel
def update_trackers(
    wave_field: ti.template(),  # type: ignore
    trackers: ti.template(),  # type: ignore
    dt_rs: ti.f32,  # type: ignore
    elapsed_t_rs: ti.f32,  # type: ignore
):
    """
    Update per-voxel amplitude and frequency trackers from ψ.

    Reads: wave_field.psi_am, wave_field.psi_prev_am
    Writes: trackers.amp_local_emarms_am,
            trackers.freq_local_cross_rHz,
            trackers.last_crossing

    Amplitude:
        |ψ|² is the squared vector magnitude (scalar). EMA:
            rms² ← α·|ψ|² + (1−α)·rms²_old
        with a small unconditional decay so stale regions fade.

    Frequency:
        Zero crossing tracked on ψ_z (polarization axis chosen by convention).
        Detect sign change between psi_prev and psi: positive-going crossing
        marks one period boundary. Period_rs = elapsed − last_crossing,
        smoothed via EMA. If no crossing: untouched, decays naturally.

    Note: vector polarization is unknown a priori — using ψ_z as the canonical
    crossing channel. For radial wave-center fields (annihilation1) this is
    the longitudinal component near the equator; for wave seeds along x
    with ŷ polarization, ψ_z stays at 0, so freq tracker stays at base value.
    M5.0e will replace this with magnitude-zero-crossing or per-component tracking.
    """
    nx, ny, nz = wave_field.nx, wave_field.ny, wave_field.nz
    decay = ti.cast(0.999, ti.f32)  # ~1000 frames to ~37%
    alpha_rms = ti.cast(0.005, ti.f32)
    alpha_freq = ti.cast(0.05, ti.f32)
    min_period_rs = ti.cast(2.0, ti.f32) * dt_rs

    for i, j, k in ti.ndrange(nx, ny, nz):
        # RMS amplitude via EMA on |ψ|²
        psi_curr = wave_field.psi_am[i, j, k]
        psi_prev = wave_field.psi_prev_am[i, j, k]
        psi2 = psi_curr.norm_sqr()
        rms2_old = trackers.amp_local_emarms_am[i, j, k] ** 2
        rms2_new = alpha_rms * psi2 + (1.0 - alpha_rms) * rms2_old
        trackers.amp_local_emarms_am[i, j, k] = ti.sqrt(rms2_new) * decay

        # Zero-crossing frequency on ψ_z (canonical channel)
        prev_z = psi_prev[2]
        curr_z = psi_curr[2]
        if prev_z < 0.0 and curr_z >= 0.0:
            period_rs = elapsed_t_rs - trackers.last_crossing[i, j, k]
            if period_rs > min_period_rs:
                measured_freq = 1.0 / period_rs
                old_freq = trackers.freq_local_cross_rHz[i, j, k]
                trackers.freq_local_cross_rHz[i, j, k] = (
                    alpha_freq * measured_freq + (1.0 - alpha_freq) * old_freq
                )
            trackers.last_crossing[i, j, k] = elapsed_t_rs

        # Unconditional frequency decay (counteracted by zero-crossing updates)
        trackers.freq_local_cross_rHz[i, j, k] = trackers.freq_local_cross_rHz[i, j, k] * decay


# ================================================================
# MATRIX-SUBSTRATE TRACKERS (M5.4) — amplitude ‖M−D‖_F + clock ‖Ṁ‖_F
# ================================================================
# The matrix-world replacement for update_trackers' ψ observables (4b
# tracker-redefinitions table). Reuses the SAME Trackers fields, EMA cadence, and
# 3-plane aggregation (sample_avg_trackers) — only the sampled quantity changes:
#   - amplitude = ‖M − D_vac‖_F  (Frobenius deviation from the ẑ-vacuum order
#       parameter D_vac = diag(δ, δ, 1)) → thermal A of SABER's joint (A, ω)
#   - frequency = ‖Ṁ‖_F = ‖M − M_prev‖_F / dt  (order-parameter rotation rate)
#       → de Broglie clock ω (M5.8 headline) / thermal ω. Zero for a static
#         configuration; becomes the genuine clock once M is leapfrogged (M5.6+).
# The ψ_z zero-crossing it replaces was always a convention hack (the update_trackers
# docstring admits it) — the matrix rotation rate is the real thing. This is the
# "measurement infra lands EARLY" deliverable from the M5.3 thermal-prereq analysis.


@ti.kernel
def update_trackers_M(
    wave_field: ti.template(),  # type: ignore
    trackers: ti.template(),  # type: ignore
    dt_rs: ti.f32,  # type: ignore
    delta: ti.f32,  # type: ignore
):
    """Matrix-substrate amplitude/frequency trackers (see section header).

    Reads:  wave_field.M_am, wave_field.M_prev_am
    Writes: trackers.amp_local_emarms_am  (EMA of ‖M − D_vac‖_F — thermal A)
            trackers.freq_local_cross_rHz (EMA of ‖Ṁ‖_F        — clock ω / thermal ω)

    Args:
        wave_field: WaveField (reads M_am, M_prev_am)
        trackers: Trackers (writes the per-voxel amp/freq EMA fields)
        dt_rs: timestep (rs) — divides ‖M − M_prev‖_F to give the rate ‖Ṁ‖_F
        delta: uniaxial minor-axis eigenvalue (wave_field.lc_delta) → D_vac diagonal
    """
    nx, ny, nz = wave_field.nx, wave_field.ny, wave_field.nz
    decay = ti.cast(0.999, ti.f32)
    alpha_rms = ti.cast(0.005, ti.f32)
    alpha_freq = ti.cast(0.05, ti.f32)
    inv_dt = 1.0 / dt_rs
    # ẑ-vacuum order parameter D_vac = δ·I + (1−δ)·ẑ⊗ẑ = diag(δ, δ, 1)
    d_vac = ti.Matrix([[delta, 0.0, 0.0], [0.0, delta, 0.0], [0.0, 0.0, 1.0]])

    for i, j, k in ti.ndrange(nx, ny, nz):
        m = wave_field.M_am[i, j, k]

        # Amplitude: EMA of ‖M − D_vac‖_F² (Frobenius), then √ — same RMS form as
        # update_trackers but on the order-parameter deviation instead of |ψ|².
        dev2 = (m - d_vac).norm_sqr()
        rms2_old = trackers.amp_local_emarms_am[i, j, k] ** 2
        rms2_new = alpha_rms * dev2 + (1.0 - alpha_rms) * rms2_old
        trackers.amp_local_emarms_am[i, j, k] = ti.sqrt(rms2_new) * decay

        # Frequency: EMA of ‖Ṁ‖_F = ‖M − M_prev‖_F / dt (frame rotation rate).
        mdot_mag = (m - wave_field.M_prev_am[i, j, k]).norm() * inv_dt
        old_freq = trackers.freq_local_cross_rHz[i, j, k]
        trackers.freq_local_cross_rHz[i, j, k] = (
            alpha_freq * mdot_mag + (1.0 - alpha_freq) * old_freq
        ) * decay


# ================================================================
# ENERGY DENSITY (HAMILTONIAN) — PER-VOXEL ENERGY FIELD
# ================================================================
# H(x) = ½|ψ̇|² + ½c²|∇ψ|² + V(ψ)
#
# Naming convention: the quantity stored is ENERGY density (aJ per voxel).
# The "_H" suffix tags which formula was used to compute it (Hamiltonian).
# Future formulas would parallel: `_L` for Lagrangian density, `_K` for
# kinetic-only, etc.
#
# Populated by compute_energyH_density into observables.energyH_density_aJ
# each step; consumed by:
#   - force_motion.compute_force_vector  → F = −∇E per wave-center
#   - sample_avg_trackers (below)         → 3-plane mean → dashboard total
#   - launcher WAVE_MENU=4 flux mesh      → visualization
#
# This is a third full-grid pass per step (alongside evolve_psi and
# update_trackers). Performance optimization (merge into evolve_psi)
# is M5.0i territory; we land it as a separate kernel for clarity in M5.0g.


@ti.kernel
def compute_energyH_density(
    wave_field: ti.template(),  # type: ignore
    observables: ti.template(),  # type: ignore
    c_amrs: ti.f32,  # type: ignore
    dt_rs: ti.f32,  # type: ignore
    m_freq_rs: ti.f32,  # type: ignore
    lambda_phi4: ti.f32,  # type: ignore
):
    """
    Compute per-voxel energy density (Hamiltonian formula) into
    observables.energyH_density_aJ.
    H = ½|ψ̇|² + ½c²|∇ψ|² + V(ψ)

    Reads:  wave_field.psi_am, wave_field.psi_prev_am
    Writes: observables.energyH_density_aJ

    Boundary handling: 1-cell halo (gradient stencil); boundary voxels left
    at zero (Dirichlet BC in ψ ⇒ ψ ≈ 0 at edges ⇒ H ≈ 0 trivially).

    Kinetic ψ̇ uses forward difference (ψ − ψ_prev)/dt because psi_new is
    consumed by swap_buffers before this kernel runs — same compromise as
    update_trackers. Central-difference ψ̇ becomes available once those
    two kernels are merged into evolve_psi (M5.0i optimization).

    PHYSICAL-SCALING TODO — M5.2: the field is currently named `_aJ` but
    actually stores `(am/rs)²` per voxel — the kernel writes raw kinematic
    Hamiltonian density without the physical-energy conversion factor
    (ρ_medium × voxel_volume_am³ × INTERNAL_ENERGY_TO_AJ ≈ matches M4's
    `E = ρ·V·(f·A)²` formula). Sufficient for M5.0g–M5.0i because:
      (a) F = −∇E only depends on the gradient (ratios survive scaling)
      (b) Tests against Exp 4 KG dispersion check ω(k), not absolute E
    But once M5.2 introduces V(ψ) terms with explicit dimensional
    couplings (Klein-Gordon mass, Close Eq. 23, LdG potential), we MUST
    apply the physical-scaling factor here so the V term and the
    kinetic+gradient terms add in the same units. See dashboard "(rel.)"
    labels in _launcher.py for the corresponding display caveat.
    """
    nx, ny, nz = wave_field.nx, wave_field.ny, wave_field.nz
    inv_dt = 1.0 / dt_rs
    half_c2 = 0.5 * c_amrs**2
    inv_2dx = 1.0 / (2.0 * wave_field.dx_am)

    for i, j, k in ti.ndrange((1, nx - 1), (1, ny - 1), (1, nz - 1)):
        psi = wave_field.psi_am[i, j, k]
        psi_prev = wave_field.psi_prev_am[i, j, k]

        # Kinetic: ½ |ψ̇|²
        psi_dot = (psi - psi_prev) * inv_dt
        kinetic = 0.5 * psi_dot.norm_sqr()

        # Gradient: ½ c² |∇ψ|²  (central-difference on each axis, component-wise)
        d_dx = (wave_field.psi_am[i + 1, j, k] - wave_field.psi_am[i - 1, j, k]) * inv_2dx
        d_dy = (wave_field.psi_am[i, j + 1, k] - wave_field.psi_am[i, j - 1, k]) * inv_2dx
        d_dz = (wave_field.psi_am[i, j, k + 1] - wave_field.psi_am[i, j, k - 1]) * inv_2dx
        gradient = half_c2 * (d_dx.norm_sqr() + d_dy.norm_sqr() + d_dz.norm_sqr())

        # Potential: V(ψ) — KG mass term (Step 2) + Mexican-hat φ⁴ (Step 4a)
        potential = pde.V_psi(psi, m_freq_rs, lambda_phi4)

        observables.energyH_density_aJ[i, j, k] = kinetic + gradient + potential


@ti.kernel
def compute_energyH_density_M(
    wave_field: ti.template(),  # type: ignore
    observables: ti.template(),  # type: ignore
    c_amrs: ti.f32,  # type: ignore
    dt_rs: ti.f32,  # type: ignore
    a: ti.f32,  # type: ignore
    b: ti.f32,  # type: ignore
    c: ti.f32,  # type: ignore
    e_scale: ti.f32,  # type: ignore
):
    """Matrix Hamiltonian density (M5.5.4) — the Eq.18 ℋ for the matrix substrate.

        H = ½‖Ṁ‖²_F  +  c²·4 Σ_{μ<ν}‖[M_μ,M_ν]‖²_F  +  V_M(M)
            kinetic        curvature (= ¼Σ‖F_μν‖², the page-18 energy)   potential

    Replaces the dormant-ψ placeholder (uniform ¼λ) for matrix xperiments — resolves
    the M5.4 WAVE_MENU=4 carry-over. Matches the evolve_M force so energy is conserved
    (simple ½‖Ṁ‖² kinetic per the 2026-05-26 decision). `e_scale` = physical-energy
    factor (ρ_medium × voxel_volume_am³ × INTERNAL_ENERGY_TO_AJ); pass 1.0 for bare units.

    Reads:  wave_field.M_am, wave_field.M_prev_am
    Writes: observables.energyH_density_aJ
    """
    nx, ny, nz = wave_field.nx, wave_field.ny, wave_field.nz
    inv_dt = 1.0 / dt_rs
    c2 = c_amrs * c_amrs
    inv_2dx = 1.0 / (2.0 * wave_field.dx_am)

    for i, j, k in ti.ndrange((1, nx - 1), (1, ny - 1), (1, nz - 1)):
        m = wave_field.M_am[i, j, k]
        m_dot = (m - wave_field.M_prev_am[i, j, k]) * inv_dt
        kinetic = 0.5 * m_dot.norm_sqr()                      # ½‖Ṁ‖²_F

        mx = (wave_field.M_am[i + 1, j, k] - wave_field.M_am[i - 1, j, k]) * inv_2dx
        my = (wave_field.M_am[i, j + 1, k] - wave_field.M_am[i, j - 1, k]) * inv_2dx
        mz = (wave_field.M_am[i, j, k + 1] - wave_field.M_am[i, j, k - 1]) * inv_2dx
        cxy = pde.commutator(mx, my)
        cxz = pde.commutator(mx, mz)
        cyz = pde.commutator(my, mz)
        curvature = 4.0 * (cxy.norm_sqr() + cxz.norm_sqr() + cyz.norm_sqr())

        potential = pde.V_M(m, a, b, c)
        observables.energyH_density_aJ[i, j, k] = e_scale * (kinetic + c2 * curvature + potential)


# Note: the global energy aggregate (mean per voxel + grid total) is computed
# by the 3-plane sample_avg_trackers below alongside amp / freq — single pass
# over the slice planes, single GPU↔CPU sync per dashboard cadence. The launcher
# multiplies the returned mean by voxel_count to get the grid-total scalar.


# ================================================================
# FRANK ELASTIC ENERGY DENSITY (M5.1)
# ================================================================
# H_F = (K/2) · |∇n̂|²   per voxel scalar
# |∇n̂|² = Σᵢ Σ_α (∂ᵢ n̂_α)²   (9 central-difference terms: 3 axes × 3 components)
#
# This is the same gradient term that's already implicit in the wave equation's
# c²·∇²ψ Laplacian. The kernel doesn't add new physics — it provides the SCALAR
# E we need for the M5.1 downstream tests:
#   - task 6 (gradient descent):     monitor E falling monotonically
#   - task 7 (Coulomb 1/d fit):       E(d) data points across separations
#   - WAVE_MENU=5 visualization:      color the flux mesh by elastic stress
#
# Naming convention: parallels energyH_density_aJ — quantity is energy density,
# _F suffix tags the Frank formula (vs _H Hamiltonian, future _L Lagrangian).
# K_frank coupling currently hard-coded to 1.0 (Exp 2 baseline); M5.6 plumbs in
# the physical elastic constants alongside the LdG axis hierarchy.


K_FRANK = 1.0  # Frank elastic coupling — Exp 2 baseline; physical scaling in M5.6


@ti.kernel
def compute_energyF_density(
    wave_field: ti.template(),  # type: ignore
    observables: ti.template(),  # type: ignore
    K_frank: ti.f32,  # type: ignore
):
    """
    Compute per-voxel Frank elastic energy density H_F = (K/2)·|∇n̂|² into
    observables.energyF_density_aJ.

    Reads:  wave_field.director_nhat  (the matrix-substrate director n̂ = principal
            eigenvector of M; |n̂|=1 by construction. M5.4: repointed from the
            retiring ψ — Frank energy is now the elastic energy of the eigenvector.)
    Writes: observables.energyF_density_aJ

    Same central-difference gradient stencil as compute_energyH_density's
    gradient term, sans the c² factor and sans the kinetic/potential terms.

    Boundary handling: 1-cell halo (gradient stencil); boundary voxels left
    at zero. This is consistent with energyH_density_aJ — both kernels skip
    the outermost voxels.

    PHYSICAL-SCALING TODO — M5.6: the field is named `_aJ` aspirationally but
    actually stores `(K_frank · (am/am)²)` per voxel (dimensionless until the
    LdG elastic constants land in physical units). Sufficient for M5.1 because:
      (a) gradient descent on F depends only on |∇n̂|² gradient (relative)
      (b) Coulomb 1/d fit cares about the FUNCTIONAL form, not absolute scale
    Once M5.6 introduces the biaxial LdG potential with physical elastic
    constants (k1, k2, k3 — splay/twist/bend), apply the appropriate factor
    here so F adds to the H total in the same physical units.
    """
    nx, ny, nz = wave_field.nx, wave_field.ny, wave_field.nz
    half_K = 0.5 * K_frank
    inv_2dx = 1.0 / (2.0 * wave_field.dx_am)

    for i, j, k in ti.ndrange((1, nx - 1), (1, ny - 1), (1, nz - 1)):
        # Gradient: 9 central-difference terms — ∂_x/y/z applied to each director component
        d_dx = (wave_field.director_nhat[i + 1, j, k] - wave_field.director_nhat[i - 1, j, k]) * inv_2dx
        d_dy = (wave_field.director_nhat[i, j + 1, k] - wave_field.director_nhat[i, j - 1, k]) * inv_2dx
        d_dz = (wave_field.director_nhat[i, j, k + 1] - wave_field.director_nhat[i, j, k - 1]) * inv_2dx
        grad_n_sqr = d_dx.norm_sqr() + d_dy.norm_sqr() + d_dz.norm_sqr()
        observables.energyF_density_aJ[i, j, k] = half_K * grad_n_sqr


# ================================================================
# WINDING-NUMBER DIAGNOSTIC (M5.1 task 8)
# ================================================================
# Compute the topological charge Q of a director field on a sphere around
# a given center. For a +1 hedgehog Q ≈ +1, anti-hedgehog Q ≈ −1, vacuum
# Q ≈ 0. Pure CPU numpy; called occasionally (once per dashboard frame) as
# a diagnostic, not per-step. Port of Exp 3's `winding_number()`.
#
# Math:
#     Q = (1/4π) ∮_S n̂ · (∂_θ n̂ × ∂_φ n̂) dθ dφ
# Where S is a sphere of given radius around the defect center.
# Method: trilinear-sample ψ on a regular (θ, φ) grid; finite-difference
# the angular derivatives; integrate.


def compute_winding_number(
    psi_np,
    center_vox,
    radius_vox,
    n_theta=32,
    n_phi=64,
):
    """Compute topological winding number Q on a sphere of radius `radius_vox`
    around `center_vox` in voxel coordinates.

    Args:
        psi_np: (nx, ny, nz, 3) numpy array of director-field values
                (typically `wave_field.psi_am.to_numpy()`).
        center_vox: (cx, cy, cz) tuple of voxel-coord defect center.
        radius_vox: sphere radius in voxel units. Should be small enough that
                    the sphere fits inside the grid but large enough to be
                    away from the defect core (e.g., 4-10 voxels).
        n_theta: number of polar samples (default 32).
        n_phi: number of azimuthal samples (default 64).

    Returns:
        Q (float): the winding number. ≈ +1 for hedgehog, −1 for anti-hedgehog,
                   0 for vacuum-only enclosure. Discrete grids give Q within
                   ~5% of the integer for a typical relaxed defect.
    """
    nx, ny, nz, _ = psi_np.shape
    cx, cy, cz = center_vox

    # Spherical sample grid
    theta = np.linspace(1e-3, np.pi - 1e-3, n_theta)
    phi = np.linspace(0.0, 2.0 * np.pi, n_phi, endpoint=False)
    dth = float(theta[1] - theta[0])
    dph = float(phi[1] - phi[0])
    TH, PH = np.meshgrid(theta, phi, indexing="ij")
    sx = cx + radius_vox * np.sin(TH) * np.cos(PH)
    sy = cy + radius_vox * np.sin(TH) * np.sin(PH)
    sz = cz + radius_vox * np.cos(TH)

    # Trilinear sample of ψ at the sphere points
    i0 = np.clip(np.floor(sx).astype(int), 0, nx - 2)
    j0 = np.clip(np.floor(sy).astype(int), 0, ny - 2)
    k0 = np.clip(np.floor(sz).astype(int), 0, nz - 2)
    fx = np.clip(sx - i0, 0.0, 1.0)
    fy = np.clip(sy - j0, 0.0, 1.0)
    fz = np.clip(sz - k0, 0.0, 1.0)

    n_s = np.zeros(sx.shape + (3,), dtype=np.float64)
    for dxi in (0, 1):
        for dyi in (0, 1):
            for dzi in (0, 1):
                w = (fx if dxi else 1 - fx) * (fy if dyi else 1 - fy) * (fz if dzi else 1 - fz)
                n_s += w[..., None] * psi_np[i0 + dxi, j0 + dyi, k0 + dzi, :]

    # Re-normalize the sampled directors (interpolation slightly shrinks |n̂|)
    norm = np.linalg.norm(n_s, axis=-1, keepdims=True)
    n_s = n_s / np.maximum(norm, 1e-12)

    # Angular derivatives (central differences; φ wraps periodic)
    dn_dtheta = np.zeros_like(n_s)
    dn_dphi = np.zeros_like(n_s)
    dn_dtheta[1:-1] = (n_s[2:] - n_s[:-2]) / (2 * dth)
    dn_dtheta[0] = (n_s[1] - n_s[0]) / dth
    dn_dtheta[-1] = (n_s[-1] - n_s[-2]) / dth
    dn_dphi[:, :] = (np.roll(n_s, -1, axis=1) - np.roll(n_s, 1, axis=1)) / (2 * dph)

    # Integrand: n̂ · (∂_θ n̂ × ∂_φ n̂). The sinθ factor is implicit in the
    # cross product for a unit-vector field on a parameterized sphere.
    cross = np.cross(dn_dtheta, dn_dphi, axis=-1)
    integrand = (n_s * cross).sum(axis=-1)
    return float(integrand.sum() * dth * dph / (4.0 * np.pi))


# ================================================================
# 3-PLANE SAMPLING FOR AVERAGE TRACKERS
# ================================================================
# PERFORMANCE NOTE: Full GPU reduction (atomic_add over all voxels) causes
# severe performance issues due to atomic contention with millions of voxels.
# The 3-plane sampling is a deliberate compromise:
# - Samples ~3N² voxels instead of N³ (e.g., 3% for 100³ grid)
# - Assumes isotropic field distribution (valid for most wave scenarios)
# - Acceptable accuracy vs massive performance gain
#
# AUTO-REDUCE CAVEAT (M5.0h, 2026-05-08): Taichi's "implicit sum reduction"
# pattern (`s = 0; for ... in ti.ndrange: s += ...`) lowers to atomic_add on
# the Metal backend and hits the same contention wall as explicit atomic_add.
# At 63³ × ~5 reductions/step the test pinned the GPU at 100% and stalled
# at ~25 steps/min. Workaround used in m5_0h_dispersion: replace the global
# inner product with sparse point sampling (one voxel read per mode of
# interest, FFT recovers ω from the mixed time series). Don't assume "auto-
# reduction" sidesteps this on Metal — it doesn't.
# ================================================================

# Cached slice buffers — separate caches per domain (Trackers vs FieldObservables)
# so each domain owns its own 3-plane sampling pass. Per `feedback_visual_rendering_priority`,
# each class is self-contained: Trackers samples amp/freq, FieldObservables samples H/F.
_slice_xy_amp = None
_slice_xy_freq = None
_slice_xz_amp = None
_slice_xz_freq = None
_slice_yz_amp = None
_slice_yz_freq = None
_slice_xy_energyH = None
_slice_xy_energyF = None
_slice_xz_energyH = None
_slice_xz_energyF = None
_slice_yz_energyH = None
_slice_yz_energyF = None


# --- TRACKERS slice-copy kernels (amp + freq) -----------------------------


@ti.kernel
def _copy_slice_xy_trackers(
    trackers: ti.template(),  # type: ignore
    slice_amp: ti.template(),  # type: ignore
    slice_freq: ti.template(),  # type: ignore
    mid_z: ti.i32,  # type: ignore
):
    """Copy amp / freq at the XY center plane (z=mid_z) to 2D buffers."""
    for i, j in slice_amp:
        slice_amp[i, j] = trackers.amp_local_emarms_am[i, j, mid_z]
        slice_freq[i, j] = trackers.freq_local_cross_rHz[i, j, mid_z]


@ti.kernel
def _copy_slice_xz_trackers(
    trackers: ti.template(),  # type: ignore
    slice_amp: ti.template(),  # type: ignore
    slice_freq: ti.template(),  # type: ignore
    mid_y: ti.i32,  # type: ignore
):
    """Copy amp / freq at the XZ center plane (y=mid_y) to 2D buffers."""
    for i, k in slice_amp:
        slice_amp[i, k] = trackers.amp_local_emarms_am[i, mid_y, k]
        slice_freq[i, k] = trackers.freq_local_cross_rHz[i, mid_y, k]


@ti.kernel
def _copy_slice_yz_trackers(
    trackers: ti.template(),  # type: ignore
    slice_amp: ti.template(),  # type: ignore
    slice_freq: ti.template(),  # type: ignore
    mid_x: ti.i32,  # type: ignore
):
    """Copy amp / freq at the YZ center plane (x=mid_x) to 2D buffers."""
    for j, k in slice_amp:
        slice_amp[j, k] = trackers.amp_local_emarms_am[mid_x, j, k]
        slice_freq[j, k] = trackers.freq_local_cross_rHz[mid_x, j, k]


# --- FIELD-OBSERVABLES slice-copy kernels (energyH + energyF) -------------


@ti.kernel
def _copy_slice_xy_observables(
    observables: ti.template(),  # type: ignore
    slice_energyH: ti.template(),  # type: ignore
    slice_energyF: ti.template(),  # type: ignore
    mid_z: ti.i32,  # type: ignore
):
    """Copy energyH / energyF at the XY center plane (z=mid_z) to 2D buffers."""
    for i, j in slice_energyH:
        slice_energyH[i, j] = observables.energyH_density_aJ[i, j, mid_z]
        slice_energyF[i, j] = observables.energyF_density_aJ[i, j, mid_z]


@ti.kernel
def _copy_slice_xz_observables(
    observables: ti.template(),  # type: ignore
    slice_energyH: ti.template(),  # type: ignore
    slice_energyF: ti.template(),  # type: ignore
    mid_y: ti.i32,  # type: ignore
):
    """Copy energyH / energyF at the XZ center plane (y=mid_y) to 2D buffers."""
    for i, k in slice_energyH:
        slice_energyH[i, k] = observables.energyH_density_aJ[i, mid_y, k]
        slice_energyF[i, k] = observables.energyF_density_aJ[i, mid_y, k]


@ti.kernel
def _copy_slice_yz_observables(
    observables: ti.template(),  # type: ignore
    slice_energyH: ti.template(),  # type: ignore
    slice_energyF: ti.template(),  # type: ignore
    mid_x: ti.i32,  # type: ignore
):
    """Copy energyH / energyF at the YZ center plane (x=mid_x) to 2D buffers."""
    for j, k in slice_energyH:
        slice_energyH[j, k] = observables.energyH_density_aJ[mid_x, j, k]
        slice_energyF[j, k] = observables.energyF_density_aJ[mid_x, j, k]


# --- Sampler functions ----------------------------------------------------


def sample_avg_trackers(
    wave_field,
    trackers,
):
    """
    Estimate global Trackers aggregates by sampling 3 orthogonal center planes.

    Computes:
        trackers.amp_global_emarms_am   ← √⟨amp²⟩  (RMS over the 3 planes)
        trackers.freq_global_avg_rHz    ← ⟨freq⟩   (mean over the 3 planes)

    Why 3-plane sampling: full 3D reductions with GPU atomic_add cause severe
    atomic contention at million-voxel scale and stalled the GUI for tens of
    seconds in early M5.0d.2 testing. 3-plane sampling reads ~3·N² voxels
    (3% of a 100³ grid), is exact for spatially-uniform fields, and is a sound
    estimator for smooth wave packets.

    Companion sampler: `sample_avg_observables` (handles energyH / energyF
    via the FieldObservables class). Each domain owns its own pass for clean
    separation of concerns.

    Args:
        wave_field: WaveField (grid dimensions)
        trackers: Trackers (reads amp/freq per-voxel fields; writes the two
            global aggregates)
    """
    global _slice_xy_amp, _slice_xy_freq
    global _slice_xz_amp, _slice_xz_freq
    global _slice_yz_amp, _slice_yz_freq

    nx, ny, nz = wave_field.nx, wave_field.ny, wave_field.nz

    # Initialize slice buffers once
    if _slice_xy_amp is None:
        _slice_xy_amp = ti.field(dtype=ti.f32, shape=(nx, ny))
        _slice_xy_freq = ti.field(dtype=ti.f32, shape=(nx, ny))
        _slice_xz_amp = ti.field(dtype=ti.f32, shape=(nx, nz))
        _slice_xz_freq = ti.field(dtype=ti.f32, shape=(nx, nz))
        _slice_yz_amp = ti.field(dtype=ti.f32, shape=(ny, nz))
        _slice_yz_freq = ti.field(dtype=ti.f32, shape=(ny, nz))

    # Copy 3 center slices to 2D buffers (parallel kernels)
    mid_x, mid_y, mid_z = nx // 2, ny // 2, nz // 2
    _copy_slice_xy_trackers(trackers, _slice_xy_amp, _slice_xy_freq, mid_z)
    _copy_slice_xz_trackers(trackers, _slice_xz_amp, _slice_xz_freq, mid_y)
    _copy_slice_yz_trackers(trackers, _slice_yz_amp, _slice_yz_freq, mid_x)

    # Transfer 2D slices to CPU for numpy operations (exclude boundary voxels)
    xy_amp = _slice_xy_amp.to_numpy()[1:-1, 1:-1]
    xy_freq = _slice_xy_freq.to_numpy()[1:-1, 1:-1]
    xz_amp = _slice_xz_amp.to_numpy()[1:-1, 1:-1]
    xz_freq = _slice_xz_freq.to_numpy()[1:-1, 1:-1]
    yz_amp = _slice_yz_amp.to_numpy()[1:-1, 1:-1]
    yz_freq = _slice_yz_freq.to_numpy()[1:-1, 1:-1]

    n_samples = xy_amp.size + xz_amp.size + yz_amp.size

    # RMS amplitude: √(⟨A²⟩) — energy-weighting convention
    total_amp_squared = (xy_amp**2).sum() + (xz_amp**2).sum() + (yz_amp**2).sum()
    trackers.amp_global_emarms_am[None] = float(np.sqrt(total_amp_squared / n_samples))

    # Plain mean for frequency
    total_freq = xy_freq.sum() + xz_freq.sum() + yz_freq.sum()
    trackers.freq_global_avg_rHz[None] = float(total_freq / n_samples)


def sample_avg_observables(
    wave_field,
    observables,
):
    """
    Estimate global FieldObservables aggregates by sampling 3 orthogonal planes.

    Computes:
        observables.energyH_global_avg_aJ ← ⟨E_H⟩  (mean Hamiltonian density)
        observables.energyF_global_avg_aJ ← ⟨E_F⟩  (mean Frank elastic density)

    Caller (launcher) derives the grid total trivially:
        energyH_total_aJ = energyH_global_avg_aJ × voxel_count

    Same 3-plane sampling rationale as `sample_avg_trackers`. Separate pass so
    FieldObservables is self-contained.

    Args:
        wave_field: WaveField (grid dimensions)
        observables: FieldObservables (reads energyH/energyF per-voxel fields;
            writes the two global aggregates)
    """
    global _slice_xy_energyH, _slice_xy_energyF
    global _slice_xz_energyH, _slice_xz_energyF
    global _slice_yz_energyH, _slice_yz_energyF

    nx, ny, nz = wave_field.nx, wave_field.ny, wave_field.nz

    # Initialize slice buffers once
    if _slice_xy_energyH is None:
        _slice_xy_energyH = ti.field(dtype=ti.f32, shape=(nx, ny))
        _slice_xy_energyF = ti.field(dtype=ti.f32, shape=(nx, ny))
        _slice_xz_energyH = ti.field(dtype=ti.f32, shape=(nx, nz))
        _slice_xz_energyF = ti.field(dtype=ti.f32, shape=(nx, nz))
        _slice_yz_energyH = ti.field(dtype=ti.f32, shape=(ny, nz))
        _slice_yz_energyF = ti.field(dtype=ti.f32, shape=(ny, nz))

    # Copy 3 center slices to 2D buffers (parallel kernels)
    mid_x, mid_y, mid_z = nx // 2, ny // 2, nz // 2
    _copy_slice_xy_observables(observables, _slice_xy_energyH, _slice_xy_energyF, mid_z)
    _copy_slice_xz_observables(observables, _slice_xz_energyH, _slice_xz_energyF, mid_y)
    _copy_slice_yz_observables(observables, _slice_yz_energyH, _slice_yz_energyF, mid_x)

    # Transfer 2D slices to CPU for numpy operations (exclude boundary voxels)
    xy_energyH = _slice_xy_energyH.to_numpy()[1:-1, 1:-1]
    xy_energyF = _slice_xy_energyF.to_numpy()[1:-1, 1:-1]
    xz_energyH = _slice_xz_energyH.to_numpy()[1:-1, 1:-1]
    xz_energyF = _slice_xz_energyF.to_numpy()[1:-1, 1:-1]
    yz_energyH = _slice_yz_energyH.to_numpy()[1:-1, 1:-1]
    yz_energyF = _slice_yz_energyF.to_numpy()[1:-1, 1:-1]

    n_samples = xy_energyH.size + xz_energyH.size + yz_energyH.size

    # Plain mean for energy density H (per voxel); launcher multiplies by
    # voxel_count to get the grid total.
    total_energyH = xy_energyH.sum() + xz_energyH.sum() + yz_energyH.sum()
    observables.energyH_global_avg_aJ[None] = float(total_energyH / n_samples)

    # Plain mean for Frank elastic energy density F (per voxel).
    total_energyF = xy_energyF.sum() + xz_energyF.sum() + yz_energyF.sum()
    observables.energyF_global_avg_aJ[None] = float(total_energyF / n_samples)
