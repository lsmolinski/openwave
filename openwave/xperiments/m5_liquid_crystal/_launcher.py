"""
Xperiment Launcher

Unified launcher for wave-field xperiments featuring:
- UI-based xperiment selection and switching
- Single source of truth for rendering and UI code
- Xperiment-specific parameters in /xparameters directory
"""

import webbrowser
import importlib
import sys
import os
import time
from pathlib import Path

import taichi as ti
import numpy as np

from openwave.common import colormap, constants
from openwave.i_o import flux_mesh, render, video

import openwave.xperiments.m5_liquid_crystal.medium as medium
import openwave.xperiments.m5_liquid_crystal.particle as particle
import openwave.xperiments.m5_liquid_crystal.engine1_seeds as seeds
import openwave.xperiments.m5_liquid_crystal.engine2_pde as pde
import openwave.xperiments.m5_liquid_crystal.engine3_observables as observables
import openwave.xperiments.m5_liquid_crystal.engine4_render as viz
import openwave.xperiments.m5_liquid_crystal.force_motion as force_motion
import openwave.xperiments.m5_liquid_crystal.instrumentation as instrument

# ================================================================
# XPERIMENT PARAMETERS MANAGEMENT
# ================================================================


class XperimentManager:
    """Manages loading and switching between xperiment parameters."""

    def __init__(self):
        self.available_xperiments = self._discover_xperiments()
        self.xperiment_display_names = {}  # Cache display names from meta
        self.current_xperiment = None
        self.current_xparameters = None

    def _discover_xperiments(self):
        """Discover all available xperiment parameters in /xparameters directory."""
        parameters_dir = Path(__file__).parent / "xparameters"

        if not parameters_dir.exists():
            return []

        xperiment_files = [
            file.stem for file in parameters_dir.glob("*.py") if file.name != "__init__.py"
        ]

        return sorted(xperiment_files)

    def load_xperiment(self, xperiment_name):
        """Load xperiment parameters by name.

        Args:
            xperiment_name: Parameter file name without .py extension

        Returns:
            dict: Parameters dictionary or None if loading fails
        """
        try:
            module_path = f"openwave.xperiments.m5_liquid_crystal.xparameters.{xperiment_name}"
            parameters_module = importlib.import_module(module_path)
            importlib.reload(parameters_module)  # Reload for fresh parameters

            self.current_xperiment = xperiment_name
            self.current_xparameters = parameters_module.XPARAMETERS

            # Cache display name from meta
            self.xperiment_display_names[xperiment_name] = parameters_module.XPARAMETERS["meta"][
                "X_NAME"
            ]

            return self.current_xparameters

        except Exception as e:
            print(f"Error loading xperiment '{xperiment_name}': {e}")
            return None

    def get_xperiment_display_name(self, xperiment_name):
        """Get display name from parameter meta or fallback to filename conversion."""
        if xperiment_name in self.xperiment_display_names:
            return self.xperiment_display_names[xperiment_name]

        # Fallback: try to load just for the name
        try:
            module_path = f"openwave.xperiments.m5_liquid_crystal.xparameters.{xperiment_name}"
            parameters_module = importlib.import_module(module_path)
            display_name = parameters_module.XPARAMETERS["meta"]["X_NAME"]
            self.xperiment_display_names[xperiment_name] = display_name
            return display_name
        except:
            # Final fallback: convert filename
            return " ".join(word.capitalize() for word in xperiment_name.split("_"))


# ================================================================
# SIMULATION STATE
# ================================================================


class SimulationState:
    """Manages the state of the simulation."""

    def __init__(self):
        self.wave_field = None
        self.trackers = None  # Trackers (M3/M4 wave statistics: amp, freq, EMA-RMS)
        self.observables = None  # FieldObservables (M5 derived scalars: energyH, energyF)
        self.c_amrs = 0.0
        self.dt_rs = 0.0
        self.cfl_factor = 0.0
        # V(ψ) coupling values — populated in compute_timestep, threaded through
        # to V_psi via evolve_psi + compute_energyH_density. Set to 0 here so
        # any pre-init kernel call gets free-wave behavior.
        #   m_freq_kg_rs : Klein-Gordon mass-frequency m·c²/ℏ for the electron
        #                  (rad/rs storage units; ~7.76e-7 at SIM_SPEED=1)
        #                  M5.2 Step 2 activated this in V_psi + evolve_psi.
        #   lambda_phi4  : Mexican-hat φ⁴ coupling (1/(am²·rs²)). Set to
        #                  (c_amrs/dx_am)² for grid-scale restoring force.
        #                  M5.2 Step 4a activated this — V = ¼λ(|ψ|²−1)²
        #                  moves the potential minimum to the unit sphere.
        self.m_freq_kg_rs = 0.0
        self.lambda_phi4 = 0.0
        # M5.5.4 — Eq.13 LdG V(M)=a·Tr(M²)−b·Tr(M³)+c·(Tr(M²))² couplings for evolve_M /
        # compute_energyH_density_M. Default 0 → free curvature dynamics (V off); the exact
        # Λ=(1,δ,0)-producing values are Q7 (Duda open). Expose as xparameters later.
        self.ldg_a = 0.0
        self.ldg_b = 0.0
        self.ldg_c = 0.0
        self.ldg_v0 = 0.0  # vacuum potential V_M(D), subtracted in the energyH display (M5.6.5c)
        self.elapsed_t_rs = 0.0
        self.clock_start_time = time.time()
        self.frame = 1
        # Field aggregates — populated by sample_avg_trackers; physical units
        # except for the energy attrs (see notes below).
        self.amp_global_rms = 0.0  # m
        self.freq_global_avg = 0.0  # Hz
        self.wavelength_global_avg = 0.0  # m, derived from freq_global_avg
        # Energy attrs — currently in SCALED units (per-voxel (am/rs)² ×
        # ATTOJOULE), not true J. The kernel that populates them
        # (compute_energyH_density) doesn't yet apply the physical scaling
        # factor (ρ_medium × voxel_volume × INTERNAL_ENERGY_TO_AJ). Useful as
        # *relative* observables for trends and force-gradient sourcing
        # (F = −∇E only cares about gradients). REVIEW IN M5.2 when nonlinear
        # V(ψ) couplings land alongside the physical-scaling factor.
        self.energyH_global_avg = 0.0  # per-voxel mean in scaled "(rel.)" units
        self.energyH_total = 0.0  # mean × voxel_count, scaled "(rel.)" units
        # Frank elastic density (M5.1 task 5) — populated by sample_avg_trackers.
        self.energyF_global_avg = 0.0  # per-voxel mean in scaled "(rel.)" units
        # Resolution metric: voxels-per-wavelength (xperiment-driven).
        # 0.0 means "no reference λ declared" → dashboard shows n/a
        self.wave_res = 0.0

        # Current xperiment parameters
        self.X_NAME = ""
        self.CAM_INIT = [2.00, 1.50, 1.75]
        self.UNIVERSE_SIZE = []
        self.TARGET_VOXELS = 1e8
        self.NUM_SOURCES = 1
        self.SOURCES_POSITION = []
        self.SOURCES_OFFSET_DEG = []
        self.INIT_VELOCITY = None
        self.APPLY_MOTION = True

        # UI control variables
        self.SHOW_AXIS = False
        self.TICK_SPACING = 0.25
        self.SHOW_GRID = False
        self.SHOW_EDGES = False
        self.VIZ_STRIDE = 4  # sampling stride for directors AND granules (every Nth voxel)
        self.SHOW_GLYPHS = 0  # director-glyph planes (0=off, 1=XY, 2=+XZ, 3=all)
        self.GLYPH_VECTOR = (
            0  # glyph state — 0=director only, 1=director+delta, 2=E field, 3=B field (∇×n̂)
        )
        self.GLYPH_SIZE = (
            0  # glyph shaft — 0=unit (director), 1=field magnitude (E:|∇·n̂|, B:‖∇×n̂‖)
        )
        self.GLYPH_COLOR = (
            1  # glyph color — 0=single (COLOR_MEDIUM, see far field), 1=field gradient
        )
        self.FLUX_MESH_PLANES = [0.5, 0.5, 0.5]
        self.SHOW_FLUX_MESH = 0
        self.WARP_MESH = 300
        self.SHOW_GRANULES = False
        self.SIM_SPEED = 1.0
        self.PAUSED = False

        # Color control variables
        self.COLOR_THEME = "OCEAN"
        self.WAVE_MENU = 1
        # VIZ.2 — WM7 EM-curl color: 0=orange ‖∇×n̂‖ magnitude (honest static default),
        # 1=bluered signed (∇×n̂)·axis → N=red/S=blue poles. Axis default ẑ (auto-axis = M5.8).
        self.CURL_COLOR = 0
        self.CURL_AXIS = ti.Vector([0.0, 0.0, 1.0], dt=ti.f32)
        # VIZ.4 — magnetic-dipole placeholder sample (M5.6.5f stage 1). When True,
        # the curl observables are overwritten each frame with an analytic dipole B
        # about DIPOLE_AXIS at DIPOLE_CENTER so the N/S coloring + B glyphs + moment
        # glyph can be validated before the real circulating B exists (M5.8). Set by
        # the _viz_sample_dipole xparameter; off for all physics runs.
        self.DIPOLE_SAMPLE = False
        self.DIPOLE_AXIS = ti.Vector([0.0, 0.0, 1.0], dt=ti.f32)
        self.DIPOLE_CENTER = [0.5, 0.5, 0.5]  # fractional grid coords
        self.DIPOLE_R0_VOX = 3.0  # core regularization radius (voxels)

        # Data Analytics & video export toggles
        self.INSTRUMENTATION = False
        self.EXPORT_VIDEO = False
        self.VIDEO_FRAMES = 24

        # Optional wave seed (test xperiment only)
        self.TEST_SEED = None
        # Optional topology seed (M5.1+ vacuum / hedgehog xperiments)
        self.TOPOLOGY_SEED = None

        # M5.1 task 6 — gradient-descent relaxation
        # pin_centers / pin_signs / n_defects are captured at seed time so the
        # relaxation kernel can soft-pin the defect cores at ±ẑ each step
        # (prevents numerical dissolution of topology on the discrete grid).
        self.pin_centers = None  # ti.field (n_defects, 3) i32 — voxel coords
        self.pin_signs = None  # ti.field (n_defects,) i32 — ±1
        self.n_defects = 0
        # Auto-relax on seed — set by xperiment via AUTO_RELAX_STEPS in topology_seed.
        # Relaxation is a numerical ground-state-finder, not physics, so there's
        # no demo / interactive use case — every seed gets relaxed once,
        # silently, before the user sees the field. M7 (thermal modulation)
        # will introduce proper interactive heat-modulation controls.
        self.AUTO_RELAX_STEPS = 0

    def apply_xparameters(self, params):
        """Apply parameters from xperiment parameter dictionary."""
        # Meta
        self.X_NAME = params["meta"]["X_NAME"]

        # Camera
        self.CAM_INIT = params["camera"]["INITIAL_POSITION"]

        # Universe
        universe = params["universe"]
        self.UNIVERSE_SIZE = list(universe["SIZE"])
        self.TARGET_VOXELS = universe["TARGET_VOXELS"]

        # Wave Centers
        sources = params["wave_centers"]
        self.NUM_SOURCES = sources["COUNT"]
        self.SOURCES_POSITION = sources["POSITION"]
        self.SOURCES_OFFSET_DEG = sources["PHASE_OFFSETS_DEG"]
        self.INIT_VELOCITY = sources.get("INIT_VELOCITY", None)
        self.APPLY_MOTION = sources["APPLY_MOTION"]

        # UI defaults
        ui = params["ui_defaults"]
        self.SHOW_AXIS = ui["SHOW_AXIS"]
        self.TICK_SPACING = ui["TICK_SPACING"]
        self.SHOW_GRID = ui["SHOW_GRID"]
        self.SHOW_EDGES = ui["SHOW_EDGES"]
        self.VIZ_STRIDE = ui.get("VIZ_STRIDE", 4)
        self.SHOW_GLYPHS = ui.get("SHOW_GLYPHS", 0)
        self.GLYPH_VECTOR = ui.get("GLYPH_VECTOR", 0)  # 0=director 1=+delta 2=E 3=B
        self.GLYPH_SIZE = ui.get("GLYPH_SIZE", 0)  # 0=unit, 1=field magnitude
        self.GLYPH_COLOR = ui.get("GLYPH_COLOR", 1)  # 0=single, 1=field gradient
        self.FLUX_MESH_PLANES = ui["FLUX_MESH_PLANES"]
        self.SHOW_FLUX_MESH = ui["SHOW_FLUX_MESH"]
        self.WARP_MESH = ui["WARP_MESH"]
        self.SHOW_GRANULES = ui["SHOW_GRANULES"]
        self.SIM_SPEED = ui.get("SIM_SPEED", 1.0)
        self.PAUSED = ui["PAUSED"]

        # Color defaults
        color = params["color_defaults"]
        self.COLOR_THEME = color["COLOR_THEME"]
        self.WAVE_MENU = color["WAVE_MENU"]
        self.CURL_COLOR = color.get("CURL_COLOR", 0)  # VIZ.2 WM7: 0=orange mag, 1=bluered N/S

        # Data Analytics & video export toggles
        diag = params["analytics"]
        self.INSTRUMENTATION = diag["INSTRUMENTATION"]
        self.EXPORT_VIDEO = diag["EXPORT_VIDEO"]
        self.VIDEO_FRAMES = diag["VIDEO_FRAMES"]

        # Optional wave seed (only present in _test_smoke xperiment)
        self.TEST_SEED = params.get("test_seed", None)
        # Optional topology seed (only present in _test_topology xperiment, M5.1+)
        self.TOPOLOGY_SEED = params.get("topology_seed", None)
        # M5.1 task 6 — auto-relaxation (optional; 0 = manual via RELAX button only)
        topo_seed = params.get("topology_seed") or {}
        self.AUTO_RELAX_STEPS = topo_seed.get("AUTO_RELAX_STEPS", 0)

        # VIZ.4 — magnetic-dipole placeholder sample (optional; off unless the
        # _viz_sample_dipole xparameter sets it). DIPOLE_AXIS doubles as CURL_AXIS
        # so the bluered N/S projection is taken along the moment direction.
        self.DIPOLE_SAMPLE = bool(topo_seed.get("DIPOLE_SAMPLE", False))
        if self.DIPOLE_SAMPLE:
            ax = topo_seed.get("DIPOLE_AXIS", [0.0, 0.0, 1.0])
            self.DIPOLE_AXIS = ti.Vector([float(ax[0]), float(ax[1]), float(ax[2])], dt=ti.f32)
            self.CURL_AXIS = ti.Vector([float(ax[0]), float(ax[1]), float(ax[2])], dt=ti.f32)
            self.DIPOLE_CENTER = topo_seed.get("CENTER", [0.5, 0.5, 0.5])
            self.DIPOLE_R0_VOX = float(topo_seed.get("DIPOLE_R0_VOX", 3.0))

    def initialize_grid(self):
        """Initialize or reinitialize the wave-field grid and wave-centers."""
        self.wave_field = medium.WaveField(
            self.UNIVERSE_SIZE,
            self.TARGET_VOXELS,
            self.FLUX_MESH_PLANES,
            viz_stride=self.VIZ_STRIDE,
        )
        self.trackers = medium.Trackers(self.wave_field.grid_size)
        self.observables = medium.FieldObservables(self.wave_field.grid_size)

        # Resolution metric: voxels-per-wavelength, declared by the active
        # xperiment. Sources (in priority order):
        #   1. TEST_SEED["VOXELS_PER_WAVELENGTH"]   — seed-driven xperiments
        #   2. (M5.2+) defect Compton wavelength    — particle xperiments
        #   3. None / 0.0                           — no reference (vacuum tests)
        if self.TEST_SEED is not None:
            self.wave_res = float(self.TEST_SEED.get("VOXELS_PER_WAVELENGTH", 0.0))
        else:
            self.wave_res = 0.0

        # Initialize wave-centers
        self.wave_center = particle.WaveCenter(
            self.wave_field.grid_size,
            self.NUM_SOURCES,
            self.SOURCES_POSITION,
            self.SOURCES_OFFSET_DEG,
            self.INIT_VELOCITY,
        )

    def compute_timestep(self):
        """Compute timestep from the 3D CFL stability bound.

        CFL condition for the leapfrog wave-equation solver:
            dt ≤ dx / (c · √3)

        We size dt at 95% of the bound — sitting exactly at 1/3 puts us on
        the marginal-stability boundary, where f32 rounding alone is enough
        to cross into amplification. cfl_factor = (SIM_SPEED · 0.95)² / 3,
        so at SIM_SPEED = 1 the factor is ≈ 0.301 (well under 1/3 ≈ 0.333).

        SIM_SPEED ∈ (0, 1] is a slow-motion knob that scales the effective
        wave speed for visualization without affecting stability.
        """
        cfl_safety = 0.95
        self.c_amrs = (
            constants.WAVE_SPEED / constants.ATTOMETER * constants.RONTOSECOND * self.SIM_SPEED
        )  # am/rs (slowed)
        # Tight CFL bound against physical c (= c_amrs / SIM_SPEED), with safety
        # factor — so cfl_factor saturates at (cfl_safety² / 3) at SIM_SPEED=1.
        self.dt_rs = (
            self.wave_field.dx_am * cfl_safety / (self.c_amrs / self.SIM_SPEED * (3**0.5))
        )  # rs
        self.cfl_factor = round((self.c_amrs * self.dt_rs / self.wave_field.dx_am) ** 2, 7)
        # V(ψ) couplings — only meaningful for director-class (topology) seeds,
        # where vacuum is |ψ|=1 (unit sphere). For wave-class seeds (TEST_SEED
        # plane wave / Gaussian — small ψ around 0), the φ⁴ Mexican-hat with
        # min at |ψ|=1 creates positive feedback (the (|ψ|²−1) factor ≈ −1
        # everywhere, so the EL term pulls ψ AWAY from 0 toward the unit
        # sphere, blowing up the small-amplitude plane wave). Gate on
        # TOPOLOGY_SEED presence to keep wave-class xperiments (M5.0 smoke
        # tests) on V(ψ)=0 free-wave dynamics.
        if self.TOPOLOGY_SEED is not None:
            # Klein-Gordon mass-frequency for the electron (m·c²/ℏ in rad/rs).
            # Scales with SIM_SPEED via c_amrs so the mass-oscillation timescale
            # stays consistent with the wave timescale under slow-motion playback.
            self.m_freq_kg_rs = self.c_amrs / constants.COMPTON_WAVELENGTH_REDUCED_ELECTRON_AM
            # Mexican-hat φ⁴ coupling. Natural grid scale is (c/dx)², but the
            # nonlinear feedback amplifies |ψ| excursions caused by the Laplacian
            # near the defect core — at λ = (c/dx)² the system NaNs around step 40
            # of a hedgehog seed. Empirically (research/scripts/m5_2_phi4_defect_survival),
            # 0.1·(c/dx)² is comfortably stable (|ψ| stays in [0.83, 1.18] over
            # 400 steps vs [0.73, 1.21] for free-wave). 0.3·(c/dx)² is the
            # marginal stability ceiling at the production grid.
            self.lambda_phi4 = 0.1 * (self.c_amrs / self.wave_field.dx_am) ** 2
            # M5.6.5c — Eq.13 LdG V(M) coupling for the matrix substrate (evolve_M /
            # compute_energyH_density_M). The b=0 amplitude well V = a·Tr(M²) + c·(Tr(M²))²
            # pins Tr(M²)→s₂*=1+δ² — confines the energy seen diluting under Evolve PDE —
            # while leaving the biaxiality δ EXACTLY flat (the canonical 3-term Eq.13 has NO
            # biaxial minimum: for any (a,b,c) the nonzero eigenvalues collapse to one λ*;
            # 5a §5f). Natural cubic-balance unit is c²/dx⁴ (production sweep m5_6_5c_prod_scale:
            # confines 3.3× across k∈[0.5,25] with no blow-up, dt²-stable). OFF unless the
            # xperiment sets LDG_STIFFNESS_K > 0. b=0 is the interim — a fully biaxial-STABLE
            # vacuum needs an extra invariant in V (Q7, flagged to Duda).
            ldg_k = float(self.TOPOLOGY_SEED.get("LDG_STIFFNESS_K", 0.0))
            if ldg_k > 0.0:
                delta = float(self.TOPOLOGY_SEED.get("BIAXIAL_DELTA", 0.3))
                self.ldg_c = ldg_k * self.c_amrs**2 / self.wave_field.dx_am**4
                self.ldg_a = -2.0 * self.ldg_c * (1.0 + delta * delta)  # min at s₂*=1+δ²
                self.ldg_b = 0.0
                # vacuum potential V_M(D=diag(1,δ,0)) — subtracted in the energyH display so
                # the negative well bottom doesn't swamp the structure (kernel docstring).
                tr2_vac = 1.0 + delta * delta
                tr3_vac = 1.0 + delta * delta * delta
                self.ldg_v0 = self.ldg_a * tr2_vac - self.ldg_b * tr3_vac + self.ldg_c * tr2_vac**2
            else:
                self.ldg_a = self.ldg_b = self.ldg_c = self.ldg_v0 = 0.0
        else:
            # Wave-class seed (TEST_SEED plane wave) — disable V(ψ) entirely
            # so evolve_psi runs the free-wave equation `∂²ψ/∂t² = c²·∇²ψ`.
            self.m_freq_kg_rs = 0.0
            self.lambda_phi4 = 0.0
            self.ldg_a = self.ldg_b = self.ldg_c = self.ldg_v0 = 0.0

    def reset_sim(self):
        """Reset simulation state."""
        self.wave_field = None
        self.trackers = None
        self.c_amrs = 0.0
        self.dt_rs = 0.0
        self.cfl_factor = 0.0
        self.m_freq_kg_rs = 0.0
        self.lambda_phi4 = 0.0
        self.ldg_a = self.ldg_b = self.ldg_c = self.ldg_v0 = 0.0
        self.elapsed_t_rs = 0.0
        self.clock_start_time = time.time()
        self.frame = 1
        self.amp_global_rms = 0.0
        self.freq_global_avg = 0.0
        self.wavelength_global_avg = 0.0
        self.energyH_global_avg = 0.0
        self.energyH_total = 0.0
        self.energyF_global_avg = 0.0
        self.wave_res = 0.0
        self.initialize_grid()
        self.compute_timestep()
        initialize_xperiment(self)


# ================================================================
# UI OVERLAY WINDOWS
# ================================================================


def display_xperiment_launcher(xperiment_mgr, state):
    """Display xperiment launcher UI with selectable xperiments.

    Args:
        xperiment_mgr: XperimentManager instance
        state: SimulationState instance

    Returns:
        str or None: Selected xperiment name or None
    """
    selected_xperiment = None

    with render.gui.sub_window("XPERIMENT LAUNCHER", 0.00, 0.00, 0.14, 0.33) as sub:
        sub.text("(needs window reload)", color=colormap.LIGHT_BLUE[1])
        for xp_name in xperiment_mgr.available_xperiments:
            display_name = xperiment_mgr.get_xperiment_display_name(xp_name)
            is_current = xp_name == xperiment_mgr.current_xperiment

            if sub.checkbox(display_name, is_current) and not is_current:
                selected_xperiment = xp_name

        if sub.button("Close Launcher (esc)"):
            render.window.running = False

    # TODO: remove hardcoded WIP notice and implement proper xperiment status handling
    with render.gui.sub_window("WORK-IN-PROGRESS XPERIMENT", 0.40, 0.00, 0.20, 0.08) as sub:
        sub.text("*** MODEL STILL UNDER DEVELOPMENT ***", color=colormap.RED[1])

    return selected_xperiment


def display_controls(state):
    """Display the controls UI overlay."""
    with render.gui.sub_window("CONTROLS", 0.00, 0.33, 0.16, 0.39) as sub:
        state.SHOW_AXIS = sub.checkbox(f"Axis (ticks: {state.TICK_SPACING})", state.SHOW_AXIS)
        state.SHOW_EDGES = sub.checkbox("Sim Universe Edges", state.SHOW_EDGES)
        state.SHOW_GLYPHS = sub.slider_int("Glyphs", state.SHOW_GLYPHS, 0, 3)
        if sub.checkbox("Director Vector", state.GLYPH_VECTOR == 0):
            state.GLYPH_VECTOR = 0
        if sub.checkbox("Director + Delta Vectors", state.GLYPH_VECTOR == 1):
            state.GLYPH_VECTOR = 1
        if sub.checkbox("Electric Field (divergence)", state.GLYPH_VECTOR == 2):
            state.GLYPH_VECTOR = 2
        if sub.checkbox("Magnetic Field (curl)", state.GLYPH_VECTOR == 3):
            state.GLYPH_VECTOR = 3
        # shaft length: unit (director, all glyphs visible) vs field magnitude (E: |∇·n̂| charge
        # density, B: ‖∇×n̂‖); color: single flat COLOR_MEDIUM (see far field) vs field gradient.
        if sub.checkbox("Glyph Size (unit/magnitude)", state.GLYPH_SIZE == 1):
            state.GLYPH_SIZE = 1
        else:
            state.GLYPH_SIZE = 0
        if sub.checkbox("Glyph Color (single/gradient)", state.GLYPH_COLOR == 1):
            state.GLYPH_COLOR = 1
        else:
            state.GLYPH_COLOR = 0
        state.SHOW_FLUX_MESH = sub.slider_int("Flux Mesh", state.SHOW_FLUX_MESH, 0, 3)
        state.WARP_MESH = sub.slider_int("Warp Mesh", state.WARP_MESH, 0, 50)
        # VIZ.3 — 4-state glyph select (mutually exclusive checkboxes):
        #   0=Director Vector (principal axis n̂ only),
        #   1=Director + Delta Vectors (n̂ + delta cross-bar = biaxial ellipsoid frame),
        #   2=Electric Field (+→− barb, charge-colored), 3=Magnetic Field (∇×n̂).
        # Greek δ is spelled "Delta" — GGUI cannot render Greek glyphs.
        state.SHOW_GRANULES = sub.checkbox("Show Granule Motion", state.SHOW_GRANULES)
        state.APPLY_MOTION = sub.checkbox("Apply Motion", state.APPLY_MOTION)
        state.SIM_SPEED = sub.slider_float("Speed", state.SIM_SPEED, 0.5, 1.0)
        if state.PAUSED:
            if sub.button(">> EVOLVE PDE >>"):
                state.PAUSED = False
        else:
            if sub.button("Pause"):
                state.PAUSED = True
        if sub.button("Reset Simulation"):
            state.reset_sim()


def display_wave_menu(state):
    """Display wave properties selection menu."""
    with render.gui.sub_window("WAVE MENU", 0.00, 0.79, 0.15, 0.21) as sub:
        if sub.checkbox("Deviation (Magnitude)", state.WAVE_MENU == 1):
            state.WAVE_MENU = 1
            state.wave_field.create_flux_mesh()
        if sub.checkbox("Thermal Amp (EMA RMS)", state.WAVE_MENU == 2):
            state.WAVE_MENU = 2
            state.wave_field.create_flux_mesh()
        if sub.checkbox("Thermal Clock (omega)", state.WAVE_MENU == 3):
            state.WAVE_MENU = 3
            state.wave_field.create_flux_mesh()
        if sub.checkbox("ENERGY (Hamiltonian)", state.WAVE_MENU == 4):
            state.WAVE_MENU = 4
            state.wave_field.create_flux_mesh()
        if sub.checkbox("ENERGY (Frank Elastic)", state.WAVE_MENU == 5):
            state.WAVE_MENU = 5
            state.wave_field.create_flux_mesh()
        if sub.checkbox("EM div (charge/E)", state.WAVE_MENU == 6):
            state.WAVE_MENU = 6
            state.wave_field.create_flux_mesh()
        if sub.checkbox("EM curl (rotation/B)", state.WAVE_MENU == 7):
            state.WAVE_MENU = 7
            state.wave_field.create_flux_mesh()
        # VIZ.2: WM7 color toggle — orange ‖∇×n̂‖ magnitude (off) vs bluered signed (∇×n̂)·ẑ
        # → N=red/S=blue poles (on). The vector-warp (fabric-twist) is always on for WM7.
        if state.WAVE_MENU == 7:
            if sub.checkbox("  > B color N/S (bluered)", state.CURL_COLOR == 1):
                state.CURL_COLOR = 1
            else:
                state.CURL_COLOR = 0
        # Display gradient palette with 2× average range for headroom (allows peak visualization)
        if state.WAVE_MENU == 1:  # Displacement on orange gradient
            render.canvas.triangles(og_palette_vertices, per_vertex_color=og_palette_colors)
            with render.gui.sub_window("displacement", 0.00, 0.73, 0.08, 0.06) as sub:
                sub.text(f"0       {state.amp_global_rms*2:.0e}m")
        if state.WAVE_MENU == 2:  # Thermal Amplitude (EMA RMS) on ironbow gradient
            render.canvas.triangles(ib_palette_vertices, per_vertex_color=ib_palette_colors)
            with render.gui.sub_window("amplitude", 0.00, 0.73, 0.08, 0.06) as sub:
                sub.text(f"0       {state.amp_global_rms*2:.0e}m")
        if state.WAVE_MENU == 3:  # Thermal Clock on blueprint gradient
            render.canvas.triangles(bp_palette_vertices, per_vertex_color=bp_palette_colors)
            with render.gui.sub_window("omega", 0.00, 0.73, 0.08, 0.06) as sub:
                sub.text(f"0       {state.freq_global_avg*2:.0e}Hz")
        if state.WAVE_MENU == 4:  # Energy density (Hamiltonian) on ironbow gradient
            render.canvas.triangles(ib_palette_vertices, per_vertex_color=ib_palette_colors)
            with render.gui.sub_window("energyH", 0.00, 0.73, 0.08, 0.06) as sub:
                # Unit label "rel." — value is per-voxel mean × 4 (matches the
                # colormap range max in update_flux_mesh_values). Underlying field
                # is in scaled (am/rs)² units, not physical J/m³ — REVIEW IN M5.2
                # when physical-scaling factor (ρ × voxel_volume × correction)
                # is wired into compute_energyH_density.
                sub.text(f"0      {state.energyH_global_avg*4:.0e}rel.")
        if state.WAVE_MENU == 5:  # Frank elastic density on ironbow gradient
            render.canvas.triangles(ib_palette_vertices, per_vertex_color=ib_palette_colors)
            with render.gui.sub_window("energyF", 0.00, 0.73, 0.08, 0.06) as sub:
                # Same "rel." caveat as WAVE_MENU=4 — K_frank=1.0 dimensionless
                # until M5.6 physical elastic constants land.
                sub.text(f"0      {state.energyF_global_avg*4:.0e}rel.")
        if state.WAVE_MENU == 6:  # EM divergence ∇·n̂ on greenyellow diverging gradient (signed)
            render.canvas.triangles(gy_palette_vertices, per_vertex_color=gy_palette_colors)
            with render.gui.sub_window("div (charge/E)", 0.00, 0.73, 0.08, 0.06) as sub:
                sub.text(" -           +")
        if state.WAVE_MENU == 7:  # EM curl: orange magnitude OR bluered signed (∇×n̂)·ẑ N/S
            if state.CURL_COLOR == 1:
                render.canvas.triangles(br_palette_vertices, per_vertex_color=br_palette_colors)
                with render.gui.sub_window("B (S/N)", 0.00, 0.73, 0.08, 0.06) as sub:
                    sub.text(" S           N")
            else:
                render.canvas.triangles(og_palette_vertices, per_vertex_color=og_palette_colors)
                with render.gui.sub_window("curl (rot/B)", 0.00, 0.73, 0.08, 0.06) as sub:
                    sub.text("0           max")


def display_level_specs(state, level_bar_vertices):
    """Display OpenWave level specifications overlay."""
    render.canvas.triangles(level_bar_vertices, color=colormap.ORANGE[1])
    with render.gui.sub_window("LIQUID-CRYSTAL MODEL (M5)", 0.84, 0.01, 0.16, 0.16) as sub:
        sub.text("Medium: Indexed Voxel Grid")
        sub.text("Data-Structure: Vector Field")
        sub.text("Coupling: Non-linear Lagrangian")
        sub.text("Propagation: PDE Solver")
        sub.text("Boundary: Dirichlet Condition")
        if sub.button("Wave Notation Guide"):
            webbrowser.open(
                "https://github.com/openwave-labs/openwave/blob/main/openwave/wave_notation.md"
            )


def display_data_dashboard(state):
    """Display simulation data dashboard."""
    clock_time = time.time() - state.clock_start_time
    sim_time_years = clock_time / (state.elapsed_t_rs * constants.RONTOSECOND or 1) / 31_536_000

    with render.gui.sub_window("DATA-DASHBOARD", 0.84, 0.17, 0.16, 0.60) as sub:
        state.INSTRUMENTATION = sub.checkbox("Instrumentation ON/OFF", state.INSTRUMENTATION)
        sub.text("--- SPACETIME ---", color=colormap.LIGHT_BLUE[1])
        sub.text(f"Medium Density: {constants.MEDIUM_DENSITY:.1e} kg/m³")
        sub.text(f"Wave Speed (c): {constants.WAVE_SPEED:.1e} m/s")

        sub.text("\n--- SIMULATION DOMAIN ---", color=colormap.LIGHT_BLUE[1])
        sub.text(f"Universe: {state.wave_field.max_universe_edge:.1e} m")
        sub.text(f"Voxel Count: {state.wave_field.voxel_count:,}")
        sub.text(
            f"Grid Size: {state.wave_field.nx} x {state.wave_field.ny} x {state.wave_field.nz}"
        )
        sub.text(f"Voxel Edge: {state.wave_field.dx:.2e} m")

        sub.text("\n--- RESOLUTION ---", color=colormap.LIGHT_BLUE[1])
        # wave_res is xperiment-driven; 0.0 means no reference λ declared.
        # Stable leapfrog needs ≥12 voxels/λ; <10 is undersampled.
        if state.wave_res > 0:
            sub.text(f"Wave: {state.wave_res:.1f} voxels/wave (>12)")
            if state.wave_res < 10:
                sub.text(f"*** WARNING: Undersampling! ***", color=(1.0, 0.0, 0.0))
        else:
            sub.text("Wave: n/a (no wave declared)")

        sub.text("\n--- WAVE-FIELD ---", color=colormap.LIGHT_BLUE[1])
        sub.text(f"Amplitude (avg): {state.amp_global_rms:.1e} m")
        sub.text(f"Frequency (avg): {state.freq_global_avg:.1e} Hz")
        sub.text(f"Wavelength (avg): {state.wavelength_global_avg:.1e} m")
        # Unit label "rel." — energyH_total is `mean × voxel_count × ATTOJOULE`
        # but the mean is in scaled (am/rs)² units, not actual aJ — so the value
        # isn't physically J. REVIEW IN M5.2 when physical-scaling factor lands.
        sub.text(f"Total ENERGY (H): {state.energyH_total:.1e} rel.")

        sub.text("\n--- TIME MICROSCOPE ---", color=colormap.LIGHT_BLUE[1])
        sub.text(f"Sim Steps (frames): {state.frame:,}")
        sub.text(f"Sim Time: {state.elapsed_t_rs:,.0f} rs")
        sub.text(f"Clock Time: {clock_time:.2f} s")
        sub.text(f"(1s sim time takes {sim_time_years:.0e}y)")

        sub.text("\n--- CFL TIMESTEP ---", color=colormap.LIGHT_BLUE[1])
        sub.text(f"c_amrs: {state.c_amrs:.3f} am/rs")
        sub.text(
            f"dt_rs: {state.dt_rs:.3f} rs",
            color=(1.0, 1.0, 1.0) if state.cfl_factor <= (1 / 3) else (1.0, 0.0, 0.0),
        )
        sub.text(
            f"CFL Factor: {state.cfl_factor:.3f} (target < 1/3)",
            color=(1.0, 1.0, 1.0) if state.cfl_factor <= (1 / 3) else (1.0, 0.0, 0.0),
        )


# ================================================================
# XPERIMENT RENDERING
# ================================================================


def initialize_xperiment(state):
    """Initialize color palettes, wave charger and instrumentation.

    Args:
        state: SimulationState instance with xperiment parameters
    """
    global og_palette_vertices, og_palette_colors
    global ib_palette_vertices, ib_palette_colors
    global bp_palette_vertices, bp_palette_colors
    global gy_palette_vertices, gy_palette_colors
    global br_palette_vertices, br_palette_colors
    global level_bar_vertices

    # Initialize color palette scales for gradient rendering and level indicator
    og_palette_vertices, og_palette_colors = colormap.get_palette_scale(
        colormap.orange, 0.00, 0.72, 0.079, 0.01
    )
    ib_palette_vertices, ib_palette_colors = colormap.get_palette_scale(
        colormap.ironbow, 0.00, 0.72, 0.079, 0.01
    )
    bp_palette_vertices, bp_palette_colors = colormap.get_palette_scale(
        colormap.blueprint, 0.00, 0.72, 0.079, 0.01
    )
    gy_palette_vertices, gy_palette_colors = colormap.get_palette_scale(
        colormap.greenyellow, 0.00, 0.72, 0.079, 0.01
    )
    br_palette_vertices, br_palette_colors = colormap.get_palette_scale(
        colormap.bluered, 0.00, 0.72, 0.079, 0.01
    )
    level_bar_vertices = colormap.get_level_bar_geometry(0.84, 0.00, 0.159, 0.01)

    # Optional initial-condition seed for test xperiments
    if state.TEST_SEED is not None:
        seed = state.TEST_SEED
        polarization = ti.Vector(seed["POLARIZATION"], dt=ti.f32)
        seeds.seed_gaussian(
            state.wave_field,
            state.c_amrs,
            state.dt_rs,
            seed["AMPLITUDE_AM"],
            seed["VOXELS_PER_WAVELENGTH"],
            polarization,
            seed["DIRECTION_AXIS"],
        )

    # Optional topology seed (M5.1+ vacuum / hedgehog xperiments)
    if state.TOPOLOGY_SEED is not None:
        topo = state.TOPOLOGY_SEED
        seed_mode = topo.get("MODE", "vacuum")  # "vacuum" or "hedgehog"

        if seed_mode == "vacuum":
            seeds.seed_vacuum_M(state.wave_field, state.wave_field.lc_delta)
            print("[M5.4] seeded matrix vacuum (M = δI + (1−δ)ẑ⊗ẑ everywhere)")

        elif seed_mode == "hedgehog":
            # Defects are specified in normalized [0,1] domain coordinates;
            # convert to voxel coords and bundle into ti.fields for the kernel.
            defects = topo["DEFECTS"]  # list of {"CENTER": [x,y,z], "SIGN": ±1}
            n_defects = len(defects)
            wf = state.wave_field

            centers_np = np.zeros((n_defects, 3), dtype=np.int32)
            signs_np = np.zeros(n_defects, dtype=np.int32)
            for d, defect in enumerate(defects):
                cx_norm, cy_norm, cz_norm = defect["CENTER"]
                centers_np[d, 0] = int(round(cx_norm * (wf.nx - 1)))
                centers_np[d, 1] = int(round(cy_norm * (wf.ny - 1)))
                centers_np[d, 2] = int(round(cz_norm * (wf.nz - 1)))
                signs_np[d] = int(defect["SIGN"])

            centers_field = ti.field(dtype=ti.i32, shape=(n_defects, 3))
            signs_field = ti.field(dtype=ti.i32, shape=(n_defects,))
            centers_field.from_numpy(centers_np)
            signs_field.from_numpy(signs_np)

            # D/4 in voxel units — w_vac falloff radius for the vacuum blend.
            # Default from Exp 2 is one quarter of the largest grid extent.
            domain_quarter_fraction = topo.get("DOMAIN_QUARTER_FRACTION", 0.25)
            domain_quarter_voxels = float(domain_quarter_fraction * max(wf.nx, wf.ny, wf.nz))

            seeds.seed_hedgehog_M(
                wf, centers_field, signs_field, domain_quarter_voxels, n_defects, wf.lc_delta
            )
            # Stash pin info for relaxation kernel (M5.1 task 6)
            state.pin_centers = centers_field
            state.pin_signs = signs_field
            state.n_defects = n_defects
            print(
                f"[M5.1] seeded {n_defects} hedgehog defect(s); "
                f"D/4 = {domain_quarter_voxels:.1f} voxels"
            )

        elif seed_mode == "biaxial_hedgehog":
            # M5.6.5a: a single BIAXIAL hedgehog M = O·D(s(r))·Oᵀ, D = diag(1, δ, 0)
            # (frame O=[r̂|e_Θ|e_Φ] + disclination smoothstep + radial eigenvalue melt).
            wf = state.wave_field
            center = topo.get("CENTER", [0.5, 0.5, 0.5])
            cx = int(round(center[0] * (wf.nx - 1)))
            cy = int(round(center[1] * (wf.ny - 1)))
            cz = int(round(center[2] * (wf.nz - 1)))
            r0_vox = float(topo.get("R0_FRACTION", 0.06) * max(wf.nx, wf.ny, wf.nz))
            rhoc_vox = float(topo.get("RHOC_VOXELS", 3.0))
            biaxial_delta = float(topo.get("BIAXIAL_DELTA", 0.3))
            seeds.seed_biaxial_hedgehog_M(wf, cx, cy, cz, r0_vox, rhoc_vox, biaxial_delta)
            # NO auto-relax: the biaxial M is constructed directly; relax_director_step would
            # rebuild M uniaxially from the director and destroy the biaxial structure.
            print(
                f"[M5.6.5a] seeded biaxial hedgehog D=diag(1,{biaxial_delta},0) at "
                f"({cx},{cy},{cz}); r0={r0_vox:.1f}, ρc={rhoc_vox:.1f} voxels (no relax)"
            )

        else:
            print(f"[M5.1] WARNING: unknown TOPOLOGY_SEED mode: {seed_mode!r}")

        # Optional auto-relaxation (M5.1 task 6)
        if state.AUTO_RELAX_STEPS > 0 and state.n_defects > 0:
            print(f"[M5.1] auto-relaxing for {state.AUTO_RELAX_STEPS} gradient-descent steps...")
            relax_field(state, state.AUTO_RELAX_STEPS)
            print("[M5.1] auto-relax complete")

        # VIZ.3: populate the derived eigenframe (director_nhat + director_mid +
        # eigenvalues) from the seeded M so a PAUSED boot renders the δ-clock-hand
        # glyph correctly before the first Evolve-PDE step. Cheap, runs once.
        pde.eigen_decompose(state.wave_field)

    if state.INSTRUMENTATION:
        print("\n" + "=" * 64)
        print("INSTRUMENTATION ENABLED")
        print("=" * 64)


def relax_field(state, n_steps):
    """Run N gradient-descent relaxation steps on the director field (M5.1 task 6).

    Lowers the Frank elastic energy by smoothing the seeded hedgehog blend
    zone, while preserving topology via soft-core pinning. After all N steps,
    psi_am and psi_prev_am hold the same relaxed field so subsequent
    EVOLVE PSI sees ψ̇ = 0 (no spurious time-derivative artifact).

    Args:
        state: SimulationState with wave_field + pin_centers/signs/n_defects
        n_steps: number of relax_director_step iterations to run

    Step size: τ = 0.4 · dx²/6 (40% of the heat-equation CFL bound for headroom).
    Empirically Exp 2 used τ = 0.008 with dx = 0.34 → 41% of CFL; we match.
    """
    if state.pin_centers is None or state.n_defects == 0:
        # No defects → nothing to pin; gradient descent on a uniform field
        # has nothing to do (∇²n = 0 everywhere). Skip.
        return
    wf = state.wave_field
    cfl_bound = (wf.dx_am**2) / 6.0  # τ < dx²/(2·dim·K) with dim=3, K=1
    tau = 0.4 * cfl_bound
    for _ in range(n_steps):
        pde.relax_director_step(wf, tau, state.pin_centers, state.pin_signs, state.n_defects)
        # M5.4: relaxation operates on director_nhat (the principal eigenvector of M).
        # Copy the relaxed director back. NOT swap_buffers — static update.
        wf.director_nhat.copy_from(wf.director_nhat_new)
    # Rebuild M from the relaxed director so M_am, the ‖M−D‖/‖Ṁ‖ trackers, and the
    # flux-mesh orientation mode stay consistent with the relaxed director.
    pde.rebuild_M_from_director(wf, wf.lc_delta)
    # Refresh observables so dashboard + flux mesh reflect the relaxed field
    compute_field_observables(state)


def compute_propagation(state):
    """Per-step field evolution — M5.5.4: the Eq.18 matrix-action leapfrog ("Evolve PDE").

    Evolves the Landau–de Gennes order parameter M under Duda's Eq.18 action
    (simple ½‖Ṁ‖² kinetic + faithful potential, validated symplectic in
    sandbox_v5/m5_5_4_matrix_evolution_check.py):

        ∂²_t M = c²·Σ_α ∂_α G_α − dV_M(M) ,   G_α = 8 Σ_ν [[M_α,M_ν],M_ν]

    Two-pass: compute_curvature_flux (G_α) → evolve_M (leapfrog) → swap_matrix_buffers,
    then eigen_decompose refreshes the derived director_nhat from the evolved M so the
    render + the ‖M−D‖_F/‖Ṁ‖_F trackers reflect the new field. V off (ldg_a/b/c=0) →
    free curvature dynamics; the faithful curvature kinetic (degenerate metric, O-DoF)
    is the M5.6 refinement. (The retired ψ leapfrog stays dormant; see engine2_pde.)
    """
    wf = state.wave_field
    pde.compute_curvature_flux(wf)
    pde.evolve_M(wf, state.c_amrs, state.dt_rs, state.ldg_a, state.ldg_b, state.ldg_c)
    wf.swap_matrix_buffers()
    pde.eigen_decompose(wf)  # refresh director_nhat from the evolved M (for render + trackers)


def compute_field_observables(state):
    """Compute per-voxel energy densities + global aggregates from the current ψ.

    Runs every frame regardless of pause state so the flux-mesh visualization
    and dashboard reflect the current field — including the static seeded
    state when PAUSED=True. The kernels here read ψ but don't modify it, so
    they are safe to run while dynamics are halted.

    When paused at the seeded state, ψ == ψ_prev, so the kinetic term in H
    is zero and H reduces to ½c²|∇ψ|² + V(ψ) — structurally similar to F up
    to a c² factor (with V=0 in M5.1).
    """
    # PER-VOXEL MATRIX HAMILTONIAN (M5.5.4) =============================
    # H = ½‖Ṁ‖² + c²·4Σ‖[M_μ,M_ν]‖² + V_M(M)  → observables.energyH_density_aJ.
    # Replaces the dormant-ψ placeholder (uniform ¼λ) — resolves the M5.4 WAVE_MENU=4
    # carry-over. e_scale=1.0 (bare units; the physical-energy calibration is tied to a
    # reference mass scale → deferred to M5.9). Consumed by flux-mesh WAVE_MENU=4.
    observables.compute_energyH_density_M(
        state.wave_field,
        state.observables,
        state.c_amrs,
        state.dt_rs,
        state.ldg_a,
        state.ldg_b,
        state.ldg_c,
        state.ldg_v0,
        1.0,
    )

    # PER-VOXEL FRANK ELASTIC ENERGY DENSITY (M5.1 task 5) ===============
    # F = (K/2)·|∇n̂|²  → observables.energyF_density_aJ.
    # Consumed by flux-mesh WAVE_MENU=5; M5.1 task 6 (gradient-descent
    # monotone-decrease diagnostic); M5.1 task 7 (Coulomb 1/d fit).
    observables.compute_energyF_density(state.wave_field, state.observables, observables.K_FRANK)

    # EM-FROM-TILTS OBSERVABLES (M5.6.5b "see EM") ======================
    # ∇·n̂ (splay = Coulomb-charge-like) + ‖∇×n̂‖ (twist+bend = B-like circulation).
    # Consumed by flux-mesh WAVE_MENU 6 (∇·n̂, bluered) / 7 (‖∇×n̂‖, ironbow).
    observables.compute_director_em(state.wave_field, state.observables)
    # VIZ.4: overwrite the curl (B) field with an analytic dipole placeholder so the
    # N/S coloring + B glyphs render correctly before the real circulating B exists
    # (M5.8). div (charge) is left as the real seeded splay for context. Same center-
    # voxel convention as the biaxial_hedgehog seed so the dipole sits on the defect.
    if state.DIPOLE_SAMPLE:
        wf = state.wave_field
        cx = state.DIPOLE_CENTER[0] * (wf.nx - 1)
        cy = state.DIPOLE_CENTER[1] * (wf.ny - 1)
        cz = state.DIPOLE_CENTER[2] * (wf.nz - 1)
        observables.fill_dipole_sample_B(
            wf, state.observables, state.DIPOLE_AXIS, cx, cy, cz, state.DIPOLE_R0_VOX, 1.0
        )
    observables.compute_director_em_scale(state.wave_field, state.observables)

    # MATRIX-SUBSTRATE TRACKERS (M5.4) ==================================
    # ‖M−D‖_F amplitude (thermal A) + ‖Ṁ‖_F clock-ω. EMA on the current M; runs
    # every frame (incl. paused) so the flux-mesh "Thermal Amp"/"Thermal Clock"
    # modes + dashboard reflect the static seeded/relaxed state. Replaces the
    # ψ-based update_trackers (which now only fires on the retiring EVOLVE PSI path).
    observables.update_trackers_M(
        state.wave_field, state.trackers, state.dt_rs, state.wave_field.lc_delta
    )

    # IN-FRAME DATA SAMPLING & ANALYTICS ================================
    # Frame skip reduces GPU->CPU transfer overhead during dynamics; when
    # paused we run every frame because the field doesn't change and the
    # cost is dominated by cached slice copies, not the reduction.
    # Two separate samplers per the 2026-05-11 SoC refactor: Trackers and
    # FieldObservables each own their own 3-plane pass.
    if state.frame % 60 == 0 or state.frame == 10 or state.PAUSED:
        observables.sample_avg_trackers(state.wave_field, state.trackers)
        observables.sample_avg_observables(state.wave_field, state.observables)
        state.energyH_total = (
            state.observables.energyH_global_avg_aJ[None] * state.wave_field.voxel_count
        ) * constants.ATTOJOULE  # J
    state.amp_global_rms = state.trackers.amp_global_emarms_am[None] * constants.ATTOMETER  # m
    state.freq_global_avg = state.trackers.freq_global_avg_rHz[None] / constants.RONTOSECOND  # Hz
    state.wavelength_global_avg = constants.WAVE_SPEED / (state.freq_global_avg or 1)  # no 0 div
    state.energyH_global_avg = state.observables.energyH_global_avg_aJ[None] * constants.ATTOJOULE
    state.energyF_global_avg = state.observables.energyF_global_avg_aJ[None] * constants.ATTOJOULE

    if state.INSTRUMENTATION:
        instrument.log_timestep_data(state.frame, state.wave_field, state.trackers)
        if state.frame == 500:
            instrument.plot_probe_wave_profile(state.wave_field)


def compute_force_motion(state):
    """
    Compute forces and update particle motion.

    Physics:
    - Force = -grad(E) where E = rho * V * (f * A)^2
    - Motion: Euler integration of F = m * a
    """

    # Compute force from energy gradient, then integrate motion
    force_motion.compute_force_vector(
        state.wave_field,
        state.observables,
        state.wave_center,
    )
    if state.APPLY_MOTION:
        force_motion.integrate_motion_leapfrog(
            state.wave_field,
            state.wave_center,
            state.dt_rs,
        )
    else:
        # Zero-out velocities if not integrating force to motion
        for wc_idx in range(state.wave_center.num_sources):
            state.wave_center.velocity_amrs[wc_idx] = ti.Vector([0.0, 0.0, 0.0], dt=ti.f32)

    # Annihilation naturally occurs from wave physics, but needs numerical precision check
    # Detect and handle particle annihilation (opposite phase WCs meeting).
    # Threshold: WCs can be at grid diagonal positions and dt_rs may cause larger jumps.
    # M5.0d.3: hardcoded to 6 voxels (was state.wave_field.ewave_res / 2). Defects
    # don't physically have a single universal interaction radius; M5.2 will replace
    # this with a per-defect-type Compton-wavelength-based threshold.
    annihilation_threshold = 6.0  # in voxels
    force_motion.detect_annihilation(state.wave_center, annihilation_threshold)


def _curl_projection(state):
    """`(curl_radial, curl_center)` for the bluered N/S projection of ∇×n̂ (B).

    Shared by the WM7 flux mesh AND the B-field glyphs so their N/S coloring is
    identical. For the dipole sample → radial `(∇×n̂)·r̂` about the defect center
    (true N=red/S=blue poles, matching a bar magnet); general runs → fixed-axis
    projection (radial=0, center unused). center is in voxel coords.
    """
    if state.DIPOLE_SAMPLE:
        wf = state.wave_field
        return 1, ti.Vector(
            [
                state.DIPOLE_CENTER[0] * (wf.nx - 1),
                state.DIPOLE_CENTER[1] * (wf.ny - 1),
                state.DIPOLE_CENTER[2] * (wf.nz - 1),
            ],
            dt=ti.f32,
        )
    return 0, ti.Vector([0.0, 0.0, 0.0], dt=ti.f32)


def render_elements(state):
    """Render grid, edges, flux mesh and test particles."""
    if state.SHOW_GRID:
        render.scene.lines(state.wave_field.grid_lines, width=1, color=colormap.COLOR_MEDIUM[1])

    if state.SHOW_EDGES:
        render.scene.lines(state.wave_field.edge_lines, width=1, color=colormap.COLOR_MEDIUM[1])

    if state.SHOW_FLUX_MESH > 0 and state.SHOW_GRANULES == False:
        # VIZ.4: for the dipole sample, color the bluered B by the RADIAL projection
        # (∇×n̂)·r̂ about the defect center → true N/S poles (red N above / blue S below,
        # matching a bar magnet). General runs use the fixed-axis projection (radial=0).
        curl_radial, curl_center = _curl_projection(state)
        viz.update_flux_mesh_values(
            state.wave_field,
            state.trackers,
            state.observables,
            state.WAVE_MENU,
            state.WARP_MESH,
            state.CURL_COLOR,
            state.CURL_AXIS,  # VIZ.2: fixed-axis projection for bluered N/S (default ẑ)
            curl_radial,  # VIZ.4: 1 = radial (∇×n̂)·r̂ → true N/S poles (dipole sample)
            curl_center,  # defect center in voxel coords (radial mode only)
        )
        flux_mesh.render_flux_mesh(render.scene, state.wave_field, state.SHOW_FLUX_MESH)

    # M5.1 director-glyph overlay — line segments showing n̂ orientation,
    # signed-component RGB so opposite directions are visually opposite.
    # Renders on top of (or instead of) the flux mesh; fast (~768 segments).
    if state.SHOW_GLYPHS > 0:
        glyph_length = 0.02  # ~2% of universe edge in normalized [0,1] coords
        # VIZ.3 — 4-state glyph select (GLYPH_VECTOR):
        #   0 = Director n̂ only (principal axis; apolar, agnostic to size/color)
        #   1 = Director + Delta (n̂ + CYAN delta cross-bar = ellipsoid frame; apolar)
        #   2 = E-field (n̂ + +→− barb, charge-colored; honors size/color)   ← polar
        #   3 = B-field (∇×n̂, em-vector kernel)                             ← polar
        if state.GLYPH_VECTOR == 3:  # B-field vectors (∇×n̂) — em-vector kernel
            em_scale = max(
                state.observables.director_div_absmax[None],
                state.observables.director_curl_max[None],
            )  # shared distortion scale (matches WAVE_MENU 7)
            # gradient color uses the SAME signed N/S projection as the WM7 mesh
            curl_radial, curl_center = _curl_projection(state)
            viz.update_em_vector_glyphs(
                state.wave_field,
                state.observables,
                glyph_length,
                em_scale,
                state.SHOW_GLYPHS,
                state.GLYPH_SIZE,
                state.GLYPH_COLOR,
                state.CURL_AXIS,
                curl_radial,
                curl_center,
            )
            render.scene.lines(
                state.wave_field.director_glyph_vertices,
                per_vertex_color=state.wave_field.director_glyph_colors,
                width=2.0,
            )
            # half-arrow barb pass — B vectors are polar; the tip shows circulation direction
            render.scene.lines(
                state.wave_field.director_glyph_arrow_vertices,
                per_vertex_color=state.wave_field.director_glyph_arrow_colors,
                width=2.0,
            )
        else:  # Director-only (0) / Director+Delta (1) / E-field (2) — centered glyph kernel
            div_scale = state.observables.director_div_absmax[None]
            glyph_mode = 1 if state.GLYPH_VECTOR == 2 else 0  # 0=Director, 1=E-field
            show_delta = 1 if state.GLYPH_VECTOR == 1 else 0  # delta cross-bar only in state 1
            viz.update_director_glyphs(
                state.wave_field,
                state.observables,
                glyph_length,
                state.SHOW_GLYPHS,
                div_scale,
                state.GLYPH_SIZE,
                state.GLYPH_COLOR,
                glyph_mode,
                show_delta,
            )
            render.scene.lines(
                state.wave_field.director_glyph_vertices,
                per_vertex_color=state.wave_field.director_glyph_colors,
                width=2.0,
            )
            # Second-segment pass: delta cross-bar (state 1) OR E +→− barb (mode 1).
            # The arrow buffer is always written (blanked in director-only state 0).
            render.scene.lines(
                state.wave_field.director_glyph_arrow_vertices,
                per_vertex_color=state.wave_field.director_glyph_arrow_colors,
                width=2.0,
            )

    # VIZ.4 — magnetic-moment vector glyph (the placeholder dipole's axis marker).
    # Rendered independently of SHOW_GLYPHS so the moment is always visible alongside
    # the N/S-colored field. YELLOW + thick so it reads as the principal axis.
    if state.DIPOLE_SAMPLE:
        wf = state.wave_field
        viz.update_moment_glyph(
            wf,
            state.DIPOLE_AXIS,
            state.DIPOLE_CENTER[0] * (wf.nx - 1),
            state.DIPOLE_CENTER[1] * (wf.ny - 1),
            state.DIPOLE_CENTER[2] * (wf.nz - 1),
            0.18,  # normalized length — larger than voxel glyphs (~0.02)
            ti.Vector(
                [colormap.YELLOW[1][0], colormap.YELLOW[1][1], colormap.YELLOW[1][2]], dt=ti.f32
            ),
        )
        render.scene.lines(
            wf.moment_glyph_vertices,
            per_vertex_color=wf.moment_glyph_colors,
            width=2.0,
        )

    # Render granule positional displacement (only when flux mesh is active, since position
    # is sampled from full-grid displacement data)
    if state.SHOW_GRANULES and state.SHOW_FLUX_MESH > 0:
        granule_radius = 0.001  # in screen space (relative to max universe edge and scale factor)
        amp_boost = state.WARP_MESH  # Boost granule displacement for better visibility
        nx, ny = state.wave_field.nx, state.wave_field.ny
        # Shared VIZ_STRIDE drives both granules and director glyphs — same density,
        # so granules visually align with the underlying glyphs (M5.1 consolidation).
        stride = max(1, state.VIZ_STRIDE)
        sampled_nx = (nx + stride - 1) // stride
        sampled_ny = (ny + stride - 1) // stride
        num_render = sampled_nx * sampled_ny
        viz.sample_position_to_render(state.wave_field, amp_boost, stride, num_render)
        pos_np = state.wave_field.position_render.to_numpy()[:num_render]
        render.scene.particles(pos_np, granule_radius, color=colormap.COLOR_MEDIUM[1])


# ================================================================
# MAIN LOOP
# ================================================================


def main():
    """Main entry point for xperiment launcher."""
    selected_xperiment_arg = sys.argv[1] if len(sys.argv) > 1 else None

    # Initialize Taichi
    ti.init(arch=ti.gpu, log_level=ti.WARN)  # GPU preferred, suppress info logs

    # Initialize xperiment manager and state
    xperiment_mgr = XperimentManager()
    state = SimulationState()

    # Load xperiment from CLI argument or default
    default_xperiment = selected_xperiment_arg or "_topo_biaxial1_von"
    if default_xperiment not in xperiment_mgr.available_xperiments:
        print(f"Error: Xperiment '{default_xperiment}' not found!")
        return

    params = xperiment_mgr.load_xperiment(default_xperiment)
    if not params:
        return

    state.apply_xparameters(params)
    state.initialize_grid()
    state.compute_timestep()
    initialize_xperiment(state)

    # Initialize GGUI rendering
    render.init_UI(state.UNIVERSE_SIZE, state.TICK_SPACING, state.CAM_INIT)

    # Main rendering loop
    while render.window.running:
        render.init_scene(state.SHOW_AXIS)  # Initialize scene with lighting and camera

        # Handle ESC key for window close
        if render.window.is_pressed(ti.ui.ESCAPE):
            render.window.running = False
            break

        # Display UI overlays
        new_xperiment = display_xperiment_launcher(xperiment_mgr, state)
        display_controls(state)

        # Handle xperiment switching via process replacement
        if new_xperiment:
            print("\n================================================================")
            print("XPERIMENT LAUNCH")
            print(f"Now running: {new_xperiment}\n")

            sys.stdout.flush()
            sys.stderr.flush()
            render.window.running = False

            # os.execv replaces current process (macOS may show harmless warning)
            os.execv(sys.executable, [sys.executable, __file__, new_xperiment])

        # Always recompute timestep so observables have a valid dt_rs even
        # while paused (compute_energyH_density needs it for the kinetic term).
        state.compute_timestep()

        if not state.PAUSED:
            # Run simulation step and update time
            compute_propagation(state)
            compute_force_motion(state)
            state.elapsed_t_rs += state.dt_rs  # Accumulate simulation time
            state.frame += 1

        # Always compute observables — flux-mesh + dashboard reflect the
        # current ψ even when paused (user can inspect the seeded state).
        compute_field_observables(state)

        # Render scene elements
        render_elements(state)

        # Display additional UI elements and scene
        display_wave_menu(state)
        display_data_dashboard(state)
        display_level_specs(state, level_bar_vertices)
        render.show_scene()

        # Capture frame for video export (stops after VIDEO_FRAMES)
        if state.EXPORT_VIDEO:
            video.export_frame(state.frame, state.VIDEO_FRAMES)

    if state.INSTRUMENTATION:
        instrument.generate_plots()


# ================================================================
# ENTRY POINT
# ================================================================
if __name__ == "__main__":
    main()
