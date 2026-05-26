"""
LIQUID-CRYSTAL Model Medium Data-Grid

Object Classes @spacetime module.

LIQUID-CRYSTAL evolves the field ψ via a Lagrangian-derived PDE:
    ∂²_t ψ = c²·∇²ψ − ∂V/∂ψ

The field stores granule displacement (M3/M4 conceptual continuity) but is named ψ
to reflect its M5 role as the configuration variable of the Lagrangian. Time
evolution uses a triple-buffer leapfrog scheme (psi_prev_am, psi_am, psi_new_am).

────────────────────────────────────────────────────────────────────────────
M5.4 MATRIX-FIELD SUBSTRATE MIGRATION (in progress, 2026-05-26)
────────────────────────────────────────────────────────────────────────────
M5.2 closed as an informative negative: the Vector(3) director ψ cannot host
the paper's physics. The substrate becomes the Landau–de Gennes order parameter
    M(x) = O(x)·D·O^T(x)            (real-symmetric 3×3)
where D = diag(g, 1, δ) is the frozen global eigenvalue spectrum and O(x) is the
per-voxel dynamical rotation. M5.4 carries the M5.1 static-topology results onto
this substrate (feasibility proven in research/sandbox_v3/m5_3_matrix_feasibility.py).

Migration is ADDITIVE-then-cleanup to keep the engine runnable at every step:
  - M_am / M_prev_am / M_new_am  — the matrix triple buffer (the new substrate)
  - director_nhat / eigenvalues  — DERIVED per-frame via eigen_decompose (the
    principal eigenvector of M is the director the M5.1 machinery consumes)
  - psi_am / psi_prev_am / psi_new_am  — RETIRING. The wave-displacement role
    dies (no displacement vector in an orientation field); the director role
    moves to director_nhat. These Vector(3) buffers are removed in the final
    M5.4 cleanup once every consumer reads director_nhat. See 4b_rendering_features.md.
"""

import taichi as ti

from openwave.common import colormap, constants, utils

# Uniaxial minor-axis eigenvalue for M5.4 (placeholder). The order parameter is
# M = δ·I + (1−δ)·n̂⊗n̂ → eigenvalues (1, δ, δ), principal eigenvector = n̂.
# The physical biaxial δ ~ ℏ (and the gravity boost g) land in M5.6 / M5.8; M5.4
# only needs a uniaxial spectrum to reproduce the M5.1 Coulomb topology on M.
LC_DELTA = 0.5


@ti.data_oriented
class WaveField:
    """
    Lagrangian wave field on a cell-centered grid with attometer scaling.

    Triple-buffer leapfrog convention:
        psi_prev_am — ψ at t−dt (history buffer; read by stencil)
        psi_am      — ψ at t    (current buffer; read by stencil)
        psi_new_am  — ψ at t+dt (output buffer; written by leapfrog kernel)

    After each leapfrog step, swap_buffers() cycles: prev ← curr, curr ← new.
    psi_new_am is overwritten by the next step.

    AMR-readiness convention:
        Kernels MUST read grid dimensions via wave_field.nx / .ny / .nz attributes
        (or the .grid_size tuple). Do NOT bake fixed (nx, ny, nz) constants into
        @ti.kernel signatures — that would prevent the M5.6 / M5.8 AMR retrofit
        from swapping in an octree-based field-storage layer without rewriting
        kernels. Field access via wave_field.psi_am[i, j, k] is the canonical pattern.

    This class:
    - Cell-centered cubic grid
    - Attometer scaling for numerical precision (f32 fields)
    - Computed positions from indices (memory efficient)
    - Wave properties stored at each voxel
    - Asymmetric universe support (nx ≠ ny ≠ nz allowed)

    Initialization Strategy:
    1. User specifies init_universe_size [x, y, z] in meters (can be asymmetric)
    2. Compute universe volume and target voxel count
    3. Calculate cubic voxel size: dx = (volume / target_voxels)^(1/3)
    4. Compute grid dimensions: nx = int(x_size / dx), ny = int(y_size / dx), nz = int(z_size / dx)
    5. Recalculate actual universe size to fit integer voxel counts
    6. Initialize scalar and vector fields with attometer scaling for f32 precision
    """

    def __init__(
        self,
        init_universe_size,
        target_voxels,
        flux_mesh_planes=[0.5, 0.5, 0.5],
        viz_stride=4,
    ):
        """
        Initialize WaveField from universe size with automatic voxel sizing.

        Args:
            init_universe_size: Simulation domain size [x, y, z] in meters.
                Can be asymmetric. Will be adjusted to fit integer voxel counts.
            target_voxels: Desired total voxel count (impacts memory and performance).

        Note:
            Voxel size (dx) is cubic (same for all axes) to preserve wave physics.
            Grid counts (nx, ny, nz) can differ for asymmetric domain shapes.
        """
        # Compute initial grid properties (before rounding and grid symmetry)
        init_universe_volume = (
            init_universe_size[0] * init_universe_size[1] * init_universe_size[2]
        )

        # Calculate cubic voxel size from target voxel count
        # CRITICAL: voxels must remain cubic (same edge length on all axes)
        # This preserves wave physics isotropy. Only the NUMBER of voxels varies per axis.
        self.voxel_volume = init_universe_volume / target_voxels  # cubic voxels
        self.voxel_edge = self.voxel_volume ** (1 / 3)  # same as dx, dx³ = voxel volume
        self.voxel_edge_am = self.voxel_edge / constants.ATTOMETER  # in attometers
        self.dx = self.voxel_edge  # additional alias for simplicity
        self.dx_am = self.voxel_edge_am  # additional alias for simplicity

        # Calculate grid dimensions (number of complete voxels per dimension) - asymmetric
        # Uses nearest odd integer to ensure grid symmetry with unique central voxel:
        # 1. User-specified universe size is arbitrary (any float value)
        # 2. voxel_edge comes from cube root, rarely divides evenly into universe size
        # 3. Ensures integer count needed for array indexing and loop bounds
        # 4. Rounds to nearest odd integer for symmetric grid with central voxel
        # 5. Actual universe size recalculated below to fit integer voxel count
        self.grid_size = [
            utils.round_to_nearest_odd(init_universe_size[0] / self.dx),
            utils.round_to_nearest_odd(init_universe_size[1] / self.dx),
            utils.round_to_nearest_odd(init_universe_size[2] / self.dx),
        ]  # same as (nx, ny, nz)
        self.nx = self.grid_size[0]  # additional alias for simplicity
        self.ny = self.grid_size[1]  # additional alias for simplicity
        self.nz = self.grid_size[2]  # additional alias for simplicity
        self.max_grid_size = max(self.nx, self.ny, self.nz)
        self.min_grid_size = min(self.nx, self.ny, self.nz)

        # Compute total voxels (asymmetric grid)
        self.voxel_count = self.nx * self.ny * self.nz

        # Recompute actual universe dimensions to fit integer number of cubic voxels
        self.universe_size = [self.nx * self.dx, self.ny * self.dx, self.nz * self.dx]
        self.universe_size_am = [size / constants.ATTOMETER for size in self.universe_size]
        self.max_universe_edge = max(self.nx * self.dx, self.ny * self.dx, self.nz * self.dx)
        self.max_universe_edge_am = self.max_universe_edge / constants.ATTOMETER
        self.universe_volume = self.voxel_count * self.voxel_volume

        # ================================================================
        # DATA STRUCTURE & INITIALIZATION
        # ================================================================
        # PROPAGATED VECTOR FIELDS (values in attometers for f32 precision)
        # This avoids catastrophic cancellation in difference calculations
        # Scales 1e-17 m values to ~10 am, well within f32 range
        #
        # Triple-buffer leapfrog scheme:
        #     psi_prev_am — ψ at t−dt (history; read-only during step)
        #     psi_am      — ψ at t    (current; read-only during step + 6-point Laplacian)
        #     psi_new_am  — ψ at t+dt (output; written by leapfrog kernel)
        # After each step, swap_buffers() cycles prev ← curr, curr ← new.
        self.psi_am = ti.Vector.field(3, dtype=ti.f32, shape=self.grid_size)  # am, ψ at t
        self.psi_prev_am = ti.Vector.field(3, dtype=ti.f32, shape=self.grid_size)  # am, ψ at t−dt
        self.psi_new_am = ti.Vector.field(3, dtype=ti.f32, shape=self.grid_size)  # am, ψ at t+dt
        self.position_render = ti.Vector.field(3, dtype=ti.f32, shape=(self.nx * self.ny))  # flat

        # ================================================================
        # MATRIX-FIELD SUBSTRATE (M5.4) — the Landau–de Gennes order parameter
        # ================================================================
        # M(x) = O(x)·D·O^T(x), real-symmetric 3×3, stored as the full 9-component
        # ti.Matrix.field(3,3) (decided 2026-05-26: matches the proven feasibility
        # spike, ti.sym_eig-ready, no reassembly; the 6-component symmetric packing
        # can swap in behind this accessor later if memory binds). Triple buffer
        # mirrors the psi leapfrog convention — M_prev/M/M_new — for the matrix
        # evolution that lands in M5.5. In M5.4 only M_am is populated (by the
        # matrix seeders) and read (by eigen_decompose); M_prev/M_new are allocated
        # now so the buffer layout is final.
        self.M_am = ti.Matrix.field(3, 3, dtype=ti.f32, shape=self.grid_size)  # M at t
        self.M_prev_am = ti.Matrix.field(3, 3, dtype=ti.f32, shape=self.grid_size)  # M at t−dt
        self.M_new_am = ti.Matrix.field(3, 3, dtype=ti.f32, shape=self.grid_size)  # M at t+dt

        # DERIVED per-frame from M_am by engine2_pde.eigen_decompose (the lynchpin
        # kernel — see 4b_rendering_features.md). director_nhat is the principal
        # eigenvector (the apolar nematic director the M5.1 Frank/Coulomb/glyph
        # machinery consumes); eigenvalues are the (λ₁≥λ₂≥λ₃) spectrum, used by the
        # ‖M−D‖_F amplitude tracker and the biaxial-ellipsoid render option.
        # These are stateless caches: recomputed every frame, valid even when paused.
        self.director_nhat = ti.Vector.field(3, dtype=ti.f32, shape=self.grid_size)  # eigvec
        self.eigenvalues = ti.Vector.field(3, dtype=ti.f32, shape=self.grid_size)  # λ₁≥λ₂≥λ₃/voxel
        # Working buffer for the director-equivalent Coulomb relaxation (M5.4 gate):
        # relax_director_step writes the next director here, then M is rebuilt from it.
        self.director_nhat_new = ti.Vector.field(3, dtype=ti.f32, shape=self.grid_size)
        self.lc_delta = LC_DELTA  # uniaxial minor-axis eigenvalue (M = δI + (1−δ)n̂⊗n̂)

        # TODO: check need for velocity field = pressure / density
        # Wave velocity vector field (v = dψ/dt)
        # self.velocity_am = ti.Vector.field(3, dtype=ti.f32, shape=self.grid_size)  # am/s

        # TODO: Implement DERIVED SCALAR FIELDS
        # WAVELENGTH, PERIOD, phase, momentum

        # TODO: Implement DERIVED VECTOR FIELDS (directions normalized to unit vectors)
        # energy_flux, wave_direction, displacement_direction, wave_mode, wave_type

        # ================================================================
        # Grid Lines: data structures & initialization
        # ================================================================
        # Grid: optimized grid lines for rendering
        # Each line spans the entire grid dimension (e.g., from x=0 to x=1 in normalized coords)

        # Line count per direction:
        # - X-direction (parallel to X): (ny+1) × (nz+1) lines
        # - Y-direction (parallel to Y): (nx+1) × (nz+1) lines
        # - Z-direction (parallel to Z): (nx+1) × (ny+1) lines
        # Total vertices = 2 × (sum of lines)
        self.line_count = (
            (self.ny + 1) * (self.nz + 1)  # X-parallel lines
            + (self.nx + 1) * (self.nz + 1)  # Y-parallel lines
            + (self.nx + 1) * (self.ny + 1)  # Z-parallel lines
        )
        self.grid_lines = ti.Vector.field(3, dtype=ti.f32, shape=self.line_count * 2)
        self.edge_lines = ti.Vector.field(3, dtype=ti.f32, shape=24)  # 12 edges, 2 vertices each

        # Populate the grid and edge lines with normalized positions (ready for rendering)
        self.populate_grid_lines()  # initialize grid lines (already normalized)
        self.populate_edge_lines()  # initialize edge lines (already normalized)

        # ================================================================
        # Flux Mesh: data structures & initialization
        # ================================================================
        # Three orthogonal meshes intersecting at universe center
        # Each flux mesh has resolution matching simulation voxel grid
        # Vertices positioned at voxel centers for direct property sampling

        # XY Flux Mesh (at z = center): spans x and y dimensions
        # - Vertices: nx × ny grid points
        # - Indices: (nx-1) × (ny-1) quads × 6 indices (2 triangles each)
        self.fluxmesh_xy_vertices = ti.Vector.field(n=3, dtype=ti.f32, shape=(self.nx, self.ny))
        self.fluxmesh_xy_colors = ti.Vector.field(n=3, dtype=ti.f32, shape=(self.nx, self.ny))
        self.fluxmesh_xy_indices = ti.field(dtype=ti.i32, shape=(self.nx - 1, self.ny - 1, 6))

        # XZ Flux Mesh (at y = center): spans x and z dimensions
        # - Vertices: nx × nz grid points
        # - Indices: (nx-1) × (nz-1) quads × 6 indices (2 triangles each)
        self.fluxmesh_xz_vertices = ti.Vector.field(n=3, dtype=ti.f32, shape=(self.nx, self.nz))
        self.fluxmesh_xz_colors = ti.Vector.field(n=3, dtype=ti.f32, shape=(self.nx, self.nz))
        self.fluxmesh_xz_indices = ti.field(dtype=ti.i32, shape=(self.nx - 1, self.nz - 1, 6))

        # YZ Flux Mesh (at x = center): spans y and z dimensions
        # - Vertices: ny × nz grid points
        # - Indices: (ny-1) × (nz-1) quads × 6 indices (2 triangles each)
        self.fluxmesh_yz_vertices = ti.Vector.field(n=3, dtype=ti.f32, shape=(self.ny, self.nz))
        self.fluxmesh_yz_colors = ti.Vector.field(n=3, dtype=ti.f32, shape=(self.ny, self.nz))
        self.fluxmesh_yz_indices = ti.field(dtype=ti.i32, shape=(self.ny - 1, self.nz - 1, 6))

        # Initialize flux mesh (vertices, indices, colors)
        self.flux_mesh_planes = flux_mesh_planes  # positions in normalized coords
        self.fm_plane_x_idx = int(flux_mesh_planes[0] * self.nx)  # x index of flux mesh plane
        self.fm_plane_y_idx = int(flux_mesh_planes[1] * self.ny)  # y index of flux mesh plane
        self.fm_plane_z_idx = int(flux_mesh_planes[2] * self.nz)  # z index of flux mesh plane

        self.create_flux_mesh()

        # ================================================================
        # Director Glyphs (M5.1): line-segment visualization of n̂ = ψ/|ψ|
        # ================================================================
        # The flux mesh renders one scalar per voxel (magnitude or signed
        # component); for a director field, the *direction* is the structure
        # we need to see. Director glyphs draw a line segment from each
        # sampled voxel to voxel + L·n̂, with color = signed-component RGB
        # so opposite directions look opposite (red↔cyan, green↔magenta,
        # blue↔yellow). Design doc: research/2b_director_glyph_rendering.md.
        #
        # Same 3-plane sampling pattern as flux_mesh, sampled every VIZ_STRIDE
        # voxels. VIZ_STRIDE is the consolidated xparameter that ALSO drives
        # the SHOW_GRANULES sampling stride — same density for both overlays
        # so granules and glyphs visually align (granules show medium motion
        # at each sample point; glyphs show director orientation at the same
        # point). Round-up `(n + stride - 1) // stride` mirrors the granule
        # indexing in _launcher.py so last-row coverage matches.
        self.GLYPH_STRIDE = viz_stride
        nx_s = max(1, (self.nx + viz_stride - 1) // viz_stride)
        ny_s = max(1, (self.ny + viz_stride - 1) // viz_stride)
        nz_s = max(1, (self.nz + viz_stride - 1) // viz_stride)
        self.glyph_nx_s = nx_s
        self.glyph_ny_s = ny_s
        self.glyph_nz_s = nz_s
        self.n_glyphs = nx_s * ny_s + nx_s * nz_s + ny_s * nz_s
        # Two vertices per line segment (start + end); GGUI renders consecutive
        # vertex pairs as a line.
        self.director_glyph_vertices = ti.Vector.field(3, ti.f32, (2 * self.n_glyphs))
        # Per-vertex color matches the line endpoints. Both endpoints get the
        # same color so each glyph reads as a uniform-colored line.
        self.director_glyph_colors = ti.Vector.field(3, ti.f32, (2 * self.n_glyphs))
        # Flat-array offsets into the (2·n_glyphs) buffer per plane:
        self.glyph_offset_xy = 0
        self.glyph_offset_xz = nx_s * ny_s
        self.glyph_offset_yz = nx_s * ny_s + nx_s * nz_s

        # Half-arrowhead barb segments — one extra line per glyph from the
        # shaft tip to a perpendicular offset, so direction (head vs tail) is
        # visually unambiguous. Same indexing scheme + plane offsets as the
        # main shaft buffer; rendered separately by the launcher (gated on
        # engine4_render.SHOW_DIRECTOR_ARROWHEAD) so it can be toggled
        # without re-running the kernel.
        self.director_glyph_arrow_vertices = ti.Vector.field(3, ti.f32, (2 * self.n_glyphs))
        self.director_glyph_arrow_colors = ti.Vector.field(3, ti.f32, (2 * self.n_glyphs))

    def swap_buffers(self):
        """
        Cyclic shift of the triple-buffer leapfrog state.

        Called once per timestep, AFTER the leapfrog kernel has written psi_new_am.
        Shift:
            psi_prev_am ← psi_am      (the just-completed step's "current" becomes "previous")
            psi_am      ← psi_new_am  (the just-completed step's "next" becomes "current")
            psi_new_am  is left as-is — it will be overwritten by the next leapfrog step

        Implementation note: uses Taichi's field.copy_from() — two full-grid GPU copies
        per step. M5.0i profile (2026-05-08): copies are ~34 % of step time at 384³
        (≈6.9 ms of 19.4 ms total).

        WHY NOT ROTATING-POINTER (investigated and rejected in M5.0i):
            The naive optimization — rotate Python attribute names so the same three
            ti.field allocations cycle through prev/curr/new roles — silently breaks
            correctness. Cause: Taichi's `ti.template()` caches attribute lookups at
            kernel-compilation time, NOT at each call. So when @ti.kernel reads
            `wave_field.psi_am` during tracing, Taichi binds that to the SPECIFIC
            ti.field instance that wave_field.psi_am pointed to AT FIRST CALL. Later
            calls with the same `wave_field` Python object reuse the cached binding —
            attribute rotation on `wave_field` itself is invisible to the cached
            kernel. M5.0h dispersion test reproduced this: after rotation, c² recovery
            regressed from ±0.5 % to −19 % systematically, because the leapfrog kept
            reading the originally-bound prev/curr fields whose data was now in the
            wrong roles.
            The proper fix (deferred): pass fields explicitly to kernels —
            `evolve_psi(psi_prev, psi_curr, psi_new, ...)` — so each rotation
            creates a distinct template-arg tuple that Taichi compiles separately
            (3 cyclic permutations → 3 cached compilations). That's a 2–3 hr refactor
            of every kernel that touches the triple buffer; not justified by current
            budget (51 fps at 384³ is well over the 20 fps floor). Re-evaluate when
            M5.2's V(ψ) makes per-step work heavier and the 35 % win matters.
        """
        self.psi_prev_am.copy_from(self.psi_am)
        self.psi_am.copy_from(self.psi_new_am)

    @ti.kernel
    def populate_grid_lines(self):
        """
        Create optimized grid lines for GGUI rendering and visualization.
        Draws continuous lines across the entire grid face-to-face.

        For a 2×2×2 voxel grid (3×3×3 grid points):
        - 9 lines parallel to X-axis (one per (j,k) pair on YZ plane)
        - 9 lines parallel to Y-axis (one per (i,k) pair on XZ plane)
        - 9 lines parallel to Z-axis (one per (i,j) pair on XY plane)
        - Total: 27 lines × 2 vertices = 54 vertices (vs 108 in non-optimized version)

        Line calculation:
        - X-direction: (ny+1) × (nz+1) lines
        - Y-direction: (nx+1) × (nz+1) lines
        - Z-direction: (nx+1) × (ny+1) lines

        Positions stored directly in normalized coordinates (0-1 range) ready for rendering.
        Uses max_grid_size for uniform normalization across asymmetric grids.
        """
        # Grid dimensions and normalization factor
        max_dim = ti.cast(self.max_grid_size, ti.f32)

        # Calculate line counts per direction
        x_lines = (self.ny + 1) * (self.nz + 1)
        y_lines = (self.nx + 1) * (self.nz + 1)
        # z_lines = (self.nx + 1) * (self.ny + 1)  # implicit, noted as reminder

        # Parallelize over all lines using single outermost loop
        for line_idx in range(self.line_count):
            vertex_idx = line_idx * 2  # Each line has 2 vertices

            if line_idx < x_lines:
                # X-parallel lines: decode (j, k) position
                temp = line_idx
                j = temp // (self.nz + 1)
                k = temp % (self.nz + 1)

                # Line from x=0 to x=nx (normalized to 0 to nx/max)
                y_norm = ti.cast(j, ti.f32) / max_dim
                z_norm = ti.cast(k, ti.f32) / max_dim

                self.grid_lines[vertex_idx] = ti.Vector([0.0, y_norm, z_norm])
                self.grid_lines[vertex_idx + 1] = ti.Vector(
                    [ti.cast(self.nx, ti.f32) / max_dim, y_norm, z_norm]
                )

            elif line_idx < x_lines + y_lines:
                # Y-parallel lines: decode (i, k) position
                temp = line_idx - x_lines
                i = temp // (self.nz + 1)
                k = temp % (self.nz + 1)

                # Line from y=0 to y=ny (normalized to 0 to ny/max)
                x_norm = ti.cast(i, ti.f32) / max_dim
                z_norm = ti.cast(k, ti.f32) / max_dim

                self.grid_lines[vertex_idx] = ti.Vector([x_norm, 0.0, z_norm])
                self.grid_lines[vertex_idx + 1] = ti.Vector(
                    [x_norm, ti.cast(self.ny, ti.f32) / max_dim, z_norm]
                )

            else:
                # Z-parallel lines: decode (i, j) position
                temp = line_idx - x_lines - y_lines
                i = temp // (self.ny + 1)
                j = temp % (self.ny + 1)

                # Line from z=0 to z=nz (normalized to 0 to nz/max)
                x_norm = ti.cast(i, ti.f32) / max_dim
                y_norm = ti.cast(j, ti.f32) / max_dim

                self.grid_lines[vertex_idx] = ti.Vector([x_norm, y_norm, 0.0])
                self.grid_lines[vertex_idx + 1] = ti.Vector(
                    [x_norm, y_norm, ti.cast(self.nz, ti.f32) / max_dim]
                )

    def populate_edge_lines(self):
        """
        Populate bounding box edge lines for GGUI wireframe rendering.

        Creates the 12 edges of the grid bounding box, each defined by 2 vertices.
        Total: 12 edges × 2 vertices = 24 vertices stored in self.edge_lines.

        Edge layout:
        - Edges 0-3: X-parallel (4 edges at YZ face corners)
        - Edges 4-7: Y-parallel (4 edges at XZ face corners)
        - Edges 8-11: Z-parallel (4 edges at XY face corners)

        Coordinates normalized by max_grid_size for uniform rendering of
        asymmetric grids (nx ≠ ny ≠ nz supported).

        Note: CPU assignment used instead of GPU kernel - for 24 fixed values,
        kernel launch overhead exceeds computation time.
        """
        # Normalized extents for asymmetric grid
        max_dim = float(self.max_grid_size)
        x_max = self.nx / max_dim
        y_max = self.ny / max_dim
        z_max = self.nz / max_dim

        # Define 12 edges as (start, end) vertex pairs
        edges = [
            # X-parallel edges (at 4 YZ corners)
            ([0, 0, 0], [x_max, 0, 0]),
            ([0, y_max, 0], [x_max, y_max, 0]),
            ([0, 0, z_max], [x_max, 0, z_max]),
            ([0, y_max, z_max], [x_max, y_max, z_max]),
            # Y-parallel edges (at 4 XZ corners)
            ([0, 0, 0], [0, y_max, 0]),
            ([x_max, 0, 0], [x_max, y_max, 0]),
            ([0, 0, z_max], [0, y_max, z_max]),
            ([x_max, 0, z_max], [x_max, y_max, z_max]),
            # Z-parallel edges (at 4 XY corners)
            ([0, 0, 0], [0, 0, z_max]),
            ([x_max, 0, 0], [x_max, 0, z_max]),
            ([0, y_max, 0], [0, y_max, z_max]),
            ([x_max, y_max, 0], [x_max, y_max, z_max]),
        ]

        # Populate field directly (no kernel needed for 24 values)
        for i, (start, end) in enumerate(edges):
            self.edge_lines[i * 2] = start
            self.edge_lines[i * 2 + 1] = end

    @ti.kernel
    def create_flux_mesh(self):
        """
        Initialize normalized flux mesh for all three orthogonal planes.

        Creates vertex positions and triangle indices for XY, XZ, and YZ flux mesh
        positioned at the xperiment parameter. Each plane is a 2D mesh that samples wave
        properties from the voxel grid.

        Coordinate system:
        - Positions stored in normalized coordinates [0, 1] for rendering
        - Each vertex corresponds to a voxel center for direct property sampling
        - Meshes intersect at xperiment parameter

        Mesh structure:
        - Vertices: Grid of 3D positions matching voxel resolution
        - Indices: Triangle pairs forming quads (2 triangles × 3 vertices = 6 indices)
        - Colors: Initialized to colormap.COLOR_FLUXMESH, updated by update_flux_mesh_values()
        """

        # Planes position in normalized coordinates
        max_dim = ti.cast(self.max_grid_size, ti.f32)
        fm_plane_x = self.flux_mesh_planes[0] * (self.nx / max_dim)  # x position normalized
        fm_plane_y = self.flux_mesh_planes[1] * (self.ny / max_dim)  # y position normalized
        fm_plane_z = self.flux_mesh_planes[2] * (self.nz / max_dim)  # z position normalized

        # ================================================================
        # XY Plane (at z = center): spans (0→1, 0→1, 0.5)
        # ================================================================
        for i, j in ti.ndrange(self.nx, self.ny):
            # Normalized coordinates for rendering
            x_norm = (ti.cast(i, ti.f32) + 0.5) / max_dim
            y_norm = (ti.cast(j, ti.f32) + 0.5) / max_dim

            # Vertex position
            self.fluxmesh_xy_vertices[i, j] = ti.Vector([x_norm, y_norm, fm_plane_z])

            # Initialize color (will be changed by update_flux_mesh_values)
            self.fluxmesh_xy_colors[i, j] = ti.Vector(colormap.COLOR_FLUXMESH[1])

        # Triangle indices for XY plane
        for i, j in ti.ndrange(self.nx - 1, self.ny - 1):
            # Each quad = 2 triangles
            # Triangle 1: (i,j) → (i+1,j) → (i,j+1)
            self.fluxmesh_xy_indices[i, j, 0] = i * self.ny + j
            self.fluxmesh_xy_indices[i, j, 1] = (i + 1) * self.ny + j
            self.fluxmesh_xy_indices[i, j, 2] = i * self.ny + (j + 1)

            # Triangle 2: (i+1,j) → (i+1,j+1) → (i,j+1)
            self.fluxmesh_xy_indices[i, j, 3] = (i + 1) * self.ny + j
            self.fluxmesh_xy_indices[i, j, 4] = (i + 1) * self.ny + (j + 1)
            self.fluxmesh_xy_indices[i, j, 5] = i * self.ny + (j + 1)

        # ================================================================
        # XZ Plane (at y = center): spans (0→1, 0.5, 0→1)
        # ================================================================
        for i, k in ti.ndrange(self.nx, self.nz):
            # Normalized coordinates for rendering
            x_norm = (ti.cast(i, ti.f32) + 0.5) / max_dim
            z_norm = (ti.cast(k, ti.f32) + 0.5) / max_dim

            # Vertex position
            self.fluxmesh_xz_vertices[i, k] = ti.Vector([x_norm, fm_plane_y, z_norm])

            # Initialize color (will be changed by update_flux_mesh_values)
            self.fluxmesh_xz_colors[i, k] = ti.Vector(colormap.COLOR_FLUXMESH[1])

        # Triangle indices for XZ plane
        for i, k in ti.ndrange(self.nx - 1, self.nz - 1):
            # Each quad = 2 triangles
            # Triangle 1: (i,k) → (i+1,k) → (i,k+1)
            self.fluxmesh_xz_indices[i, k, 0] = i * self.nz + k
            self.fluxmesh_xz_indices[i, k, 1] = (i + 1) * self.nz + k
            self.fluxmesh_xz_indices[i, k, 2] = i * self.nz + (k + 1)

            # Triangle 2: (i+1,k) → (i+1,k+1) → (i,k+1)
            self.fluxmesh_xz_indices[i, k, 3] = (i + 1) * self.nz + k
            self.fluxmesh_xz_indices[i, k, 4] = (i + 1) * self.nz + (k + 1)
            self.fluxmesh_xz_indices[i, k, 5] = i * self.nz + (k + 1)

        # ================================================================
        # YZ Plane (at x = center): spans (0.5, 0→1, 0→1)
        # ================================================================
        for j, k in ti.ndrange(self.ny, self.nz):
            # Normalized coordinates for rendering
            y_norm = (ti.cast(j, ti.f32) + 0.5) / max_dim
            z_norm = (ti.cast(k, ti.f32) + 0.5) / max_dim

            # Vertex position
            self.fluxmesh_yz_vertices[j, k] = ti.Vector([fm_plane_x, y_norm, z_norm])

            # Initialize color (will be changed by update_flux_mesh_values)
            self.fluxmesh_yz_colors[j, k] = ti.Vector(colormap.COLOR_FLUXMESH[1])

        # Triangle indices for YZ plane
        for j, k in ti.ndrange(self.ny - 1, self.nz - 1):
            # Each quad = 2 triangles
            # Triangle 1: (j,k) → (j+1,k) → (j,k+1)
            self.fluxmesh_yz_indices[j, k, 0] = j * self.nz + k
            self.fluxmesh_yz_indices[j, k, 1] = (j + 1) * self.nz + k
            self.fluxmesh_yz_indices[j, k, 2] = j * self.nz + (k + 1)

            # Triangle 2: (j+1,k) → (j+1,k+1) → (j,k+1)
            self.fluxmesh_yz_indices[j, k, 3] = (j + 1) * self.nz + k
            self.fluxmesh_yz_indices[j, k, 4] = (j + 1) * self.nz + (k + 1)
            self.fluxmesh_yz_indices[j, k, 5] = j * self.nz + (k + 1)


@ti.data_oriented
class Trackers:
    """
    Wave-property trackers for each voxel — STATEFUL, time-integrated quantities.

    Tracks amplitude envelope (EMA RMS of |ψ|) and frequency (zero-crossing
    detection) at each grid point. These observables have temporal state
    (EMA accumulator, last-crossing timestamps) so they only update meaningfully
    while dynamics are running.

    Separation of concerns (2026-05-11 refactor): per-voxel ENERGY DENSITY
    fields (`energyH_density_aJ`, `energyF_density_aJ`) and their global
    aggregates moved to `FieldObservables` below — those are *instantaneous*
    derived scalars from the current ψ (no temporal state), so they
    conceptually belong in a different class.
    """

    def __init__(self, grid_size):
        """
        Initialize tracker fields for wave-property monitoring.

        Args:
            grid_size: Grid dimensions [nx, ny, nz] matching WaveField.
        """
        # LOCAL FIELDS per voxel (stateful)
        # Amplitude tracks A via EMA of |ψ| and RMS calculation
        # Frequency tracks local oscillation rate via zero-crossing detection
        self.amp_local_emarms_am = ti.field(dtype=ti.f32, shape=grid_size)  # am, rms amp
        self.last_crossing = ti.field(dtype=ti.f32, shape=grid_size)  # rs, last zero crossing
        self.freq_local_cross_rHz = ti.field(dtype=ti.f32, shape=grid_size)  # rHz, local frequency

        # GLOBAL AVERAGES for visualization scaling — initialized to zero; the
        # update_trackers EMA + sample_avg_trackers fill these in during sim.
        # M5 has no universal reference scale, so we let the simulation discover them.
        self.amp_global_emarms_am = ti.field(dtype=ti.f32, shape=())  # RMS all voxels
        self.freq_global_avg_rHz = ti.field(dtype=ti.f32, shape=())  # avg frequency all voxels


@ti.data_oriented
class FieldObservables:
    """
    Per-voxel derived scalar fields — STATELESS, instantaneous from current ψ.

    Each field here is a pure function of the current ψ configuration. No EMA,
    no temporal accumulator, no zero-crossing history — values are recomputed
    fresh every frame by their respective kernels. This means they are valid
    *even when dynamics are paused*: pause the sim, read the observable, get
    the current ψ's value of that observable.

    Companion to `Trackers`; the split was introduced 2026-05-11 to separate
    M3/M4-style wave statistics (stateful) from M5 field-theoretic observables
    (stateless). Each class owns its own 3-plane sampling pass for global
    aggregates (`sample_avg_trackers` vs `sample_avg_observables`).

    Future hooks (post-M5.1) — each lands as: new ti.field + dedicated kernel
    + extend `sample_avg_observables` for the global mean:
      - winding_number_local (M5.1 task 8)
      - divergence_psi (M5.2 ∇·s = 0 diagnostic)
      - curl_psi_magnitude / spin_density_s (M5.2 Close Eq. 23)
      - energy_flux (M5.3+ Poynting analog)
    """

    def __init__(self, grid_size):
        """
        Initialize derived-scalar fields.

        Args:
            grid_size: Grid dimensions [nx, ny, nz] matching WaveField.
        """
        # Per-voxel ENERGY density (HAMILTONIAN formula):
        #     H = ½|ψ̇|² + ½c²|∇ψ|² + V(ψ)
        # Naming: quantity is *energy*; the "_H" prefix denotes the formula
        # (Hamiltonian). Future formulas would parallel: `_L` for Lagrangian
        # density, `_K` for kinetic-only.
        # Populated each frame by engine3_observables.compute_energyH_density.
        # Consumed by:
        #   - force_motion.compute_force_vector  →  F = −∇E (replaces M4's
        #     postulated F = −∇E with E = ρV(fA)² formula). E is computed
        #     via the Hamiltonian formula (hence _H prefix); the physics
        #     statement F = −∇E is canonical regardless of how E is derived.
        #   - sample_avg_observables  →  3-plane-sampled global mean
        #   - launcher WAVE_MENU=4 flux-mesh visualization
        # In M5.0g–M5.1 the V(ψ) term is zero; M5.2 plugs in Klein-Gordon
        # mass + Close Eq. 23 + LdG via the V_psi hook in engine2_pde.
        self.energyH_density_aJ = ti.field(dtype=ti.f32, shape=grid_size)  # aJ-per-voxel

        # Per-voxel FRANK ELASTIC energy density (M5.1):
        #     H_F = (K/2) · |∇n̂|²
        # Computed by engine3_observables.compute_energyF_density. Same _aJ
        # convention as energyH_density_aJ (the "_F" denotes the Frank
        # formula). The K_frank coupling is hard-coded to 1.0 (Exp 2 baseline);
        # physical value lands with M5.6 LdG elastic constants.
        # Consumed by:
        #   - launcher WAVE_MENU=5 flux-mesh visualization
        #   - M5.1 task 6 gradient-descent diagnostic (monotone decrease)
        #   - M5.1 task 7 Coulomb 1/d fit (sum over volume → E(d))
        self.energyF_density_aJ = ti.field(dtype=ti.f32, shape=grid_size)  # aJ-per-voxel

        # GLOBAL AVERAGES for visualization scaling — filled by
        # sample_avg_observables (M5.1) using its own 3-plane sampling pass.
        self.energyH_global_avg_aJ = ti.field(dtype=ti.f32, shape=())  # mean energy density (H)
        self.energyF_global_avg_aJ = ti.field(dtype=ti.f32, shape=())  # mean energy density (F)


if __name__ == "__main__":
    print("\n================================================================")
    print("SMOKE TEST: DATA-GRID MODULE")
    print("================================================================")

    ti.init(arch=ti.gpu)

    # ================================================================
    # Parameters & Subatomic Objects Instantiation
    # ================================================================

    UNIVERSE_SIZE = [
        2e-15,
        2e-15,
        2e-15,
    ]  # m, simulation domain [x, y, z] dimensions (can be asymmetric)

    wave_field = WaveField(
        UNIVERSE_SIZE, target_voxels=3.5e8
    )  # 350M voxels (~14GB), 1B voxels (~40GB)

    print(f"\nGrid Statistics:")
    print(
        f"  Requested universe: [{UNIVERSE_SIZE[0]:.1e}, {UNIVERSE_SIZE[1]:.1e}, {UNIVERSE_SIZE[2]:.1e}] m"
    )
    print(
        f"  Actual universe: [{wave_field.universe_size[0]:.1e}, {wave_field.universe_size[1]:.1e}, {wave_field.universe_size[2]:.1e}] m"
    )
    print(f"  Grid size: {wave_field.nx} x {wave_field.ny} x {wave_field.nz} voxels")
    print(f"  Voxel edge: {wave_field.dx:.2e} m (cubic - same for all axes)")
    print(f"  Voxel count: {wave_field.voxel_count:,}")
    print(f"  Voxel edge (am): {wave_field.dx_am:.2f} am")
    print(f"  Universe volume: {wave_field.universe_volume:.2e} m³")
    print(f"  Note: voxels-per-wavelength resolution is now xperiment-driven")
    print(f"        (declared via TEST_SEED or defect Compton wavelength)")

    print("\n================================================================")
    print("END SMOKE TEST: DATA-GRID MODULE")
    print("================================================================")

    # Properly exit
    import sys

    sys.exit(0)
