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

import openwave.xperiments.m5_lagrangian_field.medium as medium
import openwave.xperiments.m5_lagrangian_field.particle as particle
import openwave.xperiments.m5_lagrangian_field.lagrangian_engine as lagrange
import openwave.xperiments.m5_lagrangian_field.xforce_motion as force_motion
import openwave.xperiments.m5_lagrangian_field.instrumentation as instrument

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
            module_path = f"openwave.xperiments.m5_lagrangian_field.xparameters.{xperiment_name}"
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
            module_path = f"openwave.xperiments.m5_lagrangian_field.xparameters.{xperiment_name}"
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
        self.FLUX_MESH_PLANES = [0.5, 0.5, 0.5]
        self.SHOW_FLUX_MESH = 0
        self.WARP_MESH = 300
        self.SHOW_DIRECTORS = 0  # M5.1: director-glyph overlay (0=off, 1=XY, 2=+XZ, 3=all)
        self.VIZ_STRIDE = (
            4  # M5.1: shared sampling stride for directors AND granules (every Nth voxel)
        )
        self.PARTICLE_SHELL = False
        self.SHOW_GRANULES = False
        self.SIM_SPEED = 1.0
        self.PAUSED = False

        # Color control variables
        self.COLOR_THEME = "OCEAN"
        self.WAVE_MENU = 1

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
        self.FLUX_MESH_PLANES = ui["FLUX_MESH_PLANES"]
        self.SHOW_FLUX_MESH = ui["SHOW_FLUX_MESH"]
        self.WARP_MESH = ui["WARP_MESH"]
        self.SHOW_DIRECTORS = ui.get("SHOW_DIRECTORS", 0)
        self.VIZ_STRIDE = ui.get("VIZ_STRIDE", 4)
        self.PARTICLE_SHELL = ui["PARTICLE_SHELL"]
        self.SHOW_GRANULES = ui["SHOW_GRANULES"]
        self.SIM_SPEED = ui.get("SIM_SPEED", 1.0)
        self.PAUSED = ui["PAUSED"]

        # Color defaults
        color = params["color_defaults"]
        self.COLOR_THEME = color["COLOR_THEME"]
        self.WAVE_MENU = color["WAVE_MENU"]

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

    def reset_sim(self):
        """Reset simulation state."""
        self.wave_field = None
        self.trackers = None
        self.c_amrs = 0.0
        self.dt_rs = 0.0
        self.cfl_factor = 0.0
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

    with render.gui.sub_window("XPERIMENT LAUNCHER", 0.00, 0.00, 0.14, 0.35) as sub:
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
        sub.text("*** METHOD STILL UNDER DEVELOPMENT ***", color=colormap.RED[1])

    return selected_xperiment


def display_controls(state):
    """Display the controls UI overlay."""
    with render.gui.sub_window("CONTROLS", 0.00, 0.35, 0.16, 0.33) as sub:
        state.SHOW_AXIS = sub.checkbox(f"Axis (ticks: {state.TICK_SPACING})", state.SHOW_AXIS)
        state.SHOW_EDGES = sub.checkbox("Sim Universe Edges", state.SHOW_EDGES)
        state.SHOW_FLUX_MESH = sub.slider_int("Flux Mesh", state.SHOW_FLUX_MESH, 0, 3)
        state.WARP_MESH = sub.slider_int("Warp Mesh", state.WARP_MESH, 0, 300)
        state.SHOW_DIRECTORS = sub.slider_int("Directors", state.SHOW_DIRECTORS, 0, 3)
        state.PARTICLE_SHELL = sub.checkbox("Particle Shell", state.PARTICLE_SHELL)
        state.SHOW_GRANULES = sub.checkbox("Show Granule Motion", state.SHOW_GRANULES)
        state.SIM_SPEED = sub.slider_float("Speed", state.SIM_SPEED, 0.5, 1.0)
        state.APPLY_MOTION = sub.checkbox("Apply Motion", state.APPLY_MOTION)
        if state.PAUSED:
            if sub.button(">> EVOLVE PSI >>"):
                state.PAUSED = False
        else:
            if sub.button("Pause"):
                state.PAUSED = True
        if sub.button("Reset Simulation"):
            state.reset_sim()


def display_wave_menu(state):
    """Display wave properties selection menu."""
    with render.gui.sub_window("WAVE MENU", 0.00, 0.80, 0.15, 0.20) as sub:
        if sub.checkbox("Displacement (Magnitude)", state.WAVE_MENU == 1):
            state.WAVE_MENU = 1
            state.wave_field.create_flux_mesh()
        if sub.checkbox("Amplitude (EMA RMS)", state.WAVE_MENU == 2):
            state.WAVE_MENU = 2
            state.wave_field.create_flux_mesh()
        if sub.checkbox("Frequency (L&T)", state.WAVE_MENU == 3):
            state.WAVE_MENU = 3
            state.wave_field.create_flux_mesh()
        if sub.checkbox("ENERGY (Hamiltonian)", state.WAVE_MENU == 4):
            state.WAVE_MENU = 4
            state.wave_field.create_flux_mesh()
        if sub.checkbox("ENERGY (Frank Elastic)", state.WAVE_MENU == 5):
            state.WAVE_MENU = 5
            state.wave_field.create_flux_mesh()
        # Display gradient palette with 2× average range for headroom (allows peak visualization)
        if state.WAVE_MENU == 1:  # Displacement on orange gradient
            render.canvas.triangles(og_palette_vertices, per_vertex_color=og_palette_colors)
            with render.gui.sub_window("displacement", 0.00, 0.74, 0.08, 0.06) as sub:
                sub.text(f"0       {state.amp_global_rms*2:.0e}m")
        if state.WAVE_MENU == 2:  # Amplitude (EMA RMS) on ironbow gradient
            render.canvas.triangles(ib_palette_vertices, per_vertex_color=ib_palette_colors)
            with render.gui.sub_window("amplitude", 0.00, 0.74, 0.08, 0.06) as sub:
                sub.text(f"0       {state.amp_global_rms*2:.0e}m")
        if state.WAVE_MENU == 3:  # Frequency (L&T) on blueprint gradient
            render.canvas.triangles(bp_palette_vertices, per_vertex_color=bp_palette_colors)
            with render.gui.sub_window("frequency", 0.00, 0.74, 0.08, 0.06) as sub:
                sub.text(f"0       {state.freq_global_avg*2:.0e}Hz")
        if state.WAVE_MENU == 4:  # Energy density (Hamiltonian) on ironbow gradient
            render.canvas.triangles(ib_palette_vertices, per_vertex_color=ib_palette_colors)
            with render.gui.sub_window("energyH", 0.00, 0.74, 0.08, 0.06) as sub:
                # Unit label "rel." — value is per-voxel mean × 4 (matches the
                # colormap range max in update_flux_mesh_values). Underlying field
                # is in scaled (am/rs)² units, not physical J/m³ — REVIEW IN M5.2
                # when physical-scaling factor (ρ × voxel_volume × correction)
                # is wired into compute_energyH_density.
                sub.text(f"0      {state.energyH_global_avg*4:.0e}rel.")
        if state.WAVE_MENU == 5:  # Frank elastic density on ironbow gradient
            render.canvas.triangles(ib_palette_vertices, per_vertex_color=ib_palette_colors)
            with render.gui.sub_window("energyF", 0.00, 0.74, 0.08, 0.06) as sub:
                # Same "rel." caveat as WAVE_MENU=4 — K_frank=1.0 dimensionless
                # until M5.6 physical elastic constants land.
                sub.text(f"0      {state.energyF_global_avg*4:.0e}rel.")


def display_level_specs(state, level_bar_vertices):
    """Display OpenWave level specifications overlay."""
    render.canvas.triangles(level_bar_vertices, color=colormap.ORANGE[1])
    with render.gui.sub_window("LAGRANGIAN-FIELD METHOD (M5)", 0.84, 0.01, 0.16, 0.16) as sub:
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
    global level_bar_vertices

    # Initialize color palette scales for gradient rendering and level indicator
    og_palette_vertices, og_palette_colors = colormap.get_palette_scale(
        colormap.orange, 0.00, 0.73, 0.079, 0.01
    )
    ib_palette_vertices, ib_palette_colors = colormap.get_palette_scale(
        colormap.ironbow, 0.00, 0.73, 0.079, 0.01
    )
    bp_palette_vertices, bp_palette_colors = colormap.get_palette_scale(
        colormap.blueprint, 0.00, 0.73, 0.079, 0.01
    )
    level_bar_vertices = colormap.get_level_bar_geometry(0.84, 0.00, 0.159, 0.01)

    # Optional initial-condition seed for test xperiments
    if state.TEST_SEED is not None:
        seed = state.TEST_SEED
        polarization = ti.Vector(seed["POLARIZATION"], dt=ti.f32)
        lagrange.seed_gaussian(
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
            lagrange.seed_vacuum(state.wave_field)
            print("[M5.1] seeded vacuum (n = ẑ everywhere)")

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

            lagrange.seed_hedgehog(
                wf, centers_field, signs_field, domain_quarter_voxels, n_defects
            )
            # Stash pin info for relaxation kernel (M5.1 task 6)
            state.pin_centers = centers_field
            state.pin_signs = signs_field
            state.n_defects = n_defects
            print(
                f"[M5.1] seeded {n_defects} hedgehog defect(s); "
                f"D/4 = {domain_quarter_voxels:.1f} voxels"
            )

        else:
            print(f"[M5.1] WARNING: unknown TOPOLOGY_SEED mode: {seed_mode!r}")

        # Optional auto-relaxation (M5.1 task 6)
        if state.AUTO_RELAX_STEPS > 0 and state.n_defects > 0:
            print(f"[M5.1] auto-relaxing for {state.AUTO_RELAX_STEPS} gradient-descent steps...")
            relax_field(state, state.AUTO_RELAX_STEPS)
            print("[M5.1] auto-relax complete")

    if state.INSTRUMENTATION:
        print("\n" + "=" * 64)
        print("INSTRUMENTATION ENABLED")
        print("=" * 64)


def compute_propagation(state):
    """Step ψ one timestep via leapfrog, then update trackers.

    Dynamics-only: runs the wave propagation + buffer swap + per-voxel
    amp/freq EMA. Per-voxel ENERGY DENSITIES (H, F) and the global aggregates
    are computed by `compute_field_observables`, which runs every frame
    regardless of pause state so visualization reflects the current ψ
    (including the static seeded state when PAUSED=True).
    """

    # ψ PROPAGATION =======================================
    # Leapfrog/Verlet step: ψ_new = 2·ψ − ψ_prev + (c·dt)²·∇²ψ
    lagrange.evolve_psi(state.wave_field, state.c_amrs, state.dt_rs)
    # Cycle the triple buffer: psi_prev ← psi, psi ← psi_new
    state.wave_field.swap_buffers()

    # TRACKERS (per-voxel amp / freq from ψ) ==============
    # Time-dependent EMA — only meaningful during dynamics, so kept here.
    lagrange.update_trackers(state.wave_field, state.trackers, state.dt_rs, state.elapsed_t_rs)


def relax_field(state, n_steps):
    """Run N gradient-descent relaxation steps on the director field (M5.1 task 6).

    Lowers the Frank elastic energy by smoothing the seeded hedgehog blend
    zone, while preserving topology via soft-core pinning. After all N steps,
    psi_am and psi_prev_am hold the same relaxed field so subsequent
    PROPAGATE WAVE sees ψ̇ = 0 (no spurious time-derivative artifact).

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
        lagrange.relax_director_step(wf, tau, state.pin_centers, state.pin_signs, state.n_defects)
        # Copy relaxed field back into psi_am (and psi_prev_am to keep ψ̇=0).
        # NOT swap_buffers — that would cycle prev←curr which we don't want
        # for a static (non-time-evolving) update.
        wf.psi_am.copy_from(wf.psi_new_am)
        wf.psi_prev_am.copy_from(wf.psi_new_am)
    # Refresh observables so dashboard + flux mesh reflect the relaxed field
    compute_field_observables(state)


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
    # PER-VOXEL ENERGY DENSITY (HAMILTONIAN) ============================
    # H = ½|ψ̇|² + ½c²|∇ψ|² + V(ψ)  → observables.energyH_density_aJ.
    # Consumed by xforce_motion (F = −∇E) and flux-mesh WAVE_MENU=4.
    lagrange.compute_energyH_density(
        state.wave_field, state.observables, state.c_amrs, state.dt_rs
    )

    # PER-VOXEL FRANK ELASTIC ENERGY DENSITY (M5.1 task 5) ===============
    # F = (K/2)·|∇n̂|²  → observables.energyF_density_aJ.
    # Consumed by flux-mesh WAVE_MENU=5; M5.1 task 6 (gradient-descent
    # monotone-decrease diagnostic); M5.1 task 7 (Coulomb 1/d fit).
    lagrange.compute_energyF_density(state.wave_field, state.observables, lagrange.K_FRANK)

    # IN-FRAME DATA SAMPLING & ANALYTICS ================================
    # Frame skip reduces GPU->CPU transfer overhead during dynamics; when
    # paused we run every frame because the field doesn't change and the
    # cost is dominated by cached slice copies, not the reduction.
    # Two separate samplers per the 2026-05-11 SoC refactor: Trackers and
    # FieldObservables each own their own 3-plane pass.
    if state.frame % 60 == 0 or state.frame == 10 or state.PAUSED:
        lagrange.sample_avg_trackers(state.wave_field, state.trackers)
        lagrange.sample_avg_observables(state.wave_field, state.observables)
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


def render_elements(state):
    """Render grid, edges, flux mesh and test particles."""
    if state.SHOW_GRID:
        render.scene.lines(state.wave_field.grid_lines, width=1, color=colormap.COLOR_MEDIUM[1])

    if state.SHOW_EDGES:
        render.scene.lines(state.wave_field.edge_lines, width=1, color=colormap.COLOR_MEDIUM[1])

    if state.SHOW_FLUX_MESH > 0 and state.SHOW_GRANULES == False:
        lagrange.update_flux_mesh_values(
            state.wave_field,
            state.trackers,
            state.observables,
            state.WAVE_MENU,
            state.WARP_MESH,
        )
        flux_mesh.render_flux_mesh(render.scene, state.wave_field, state.SHOW_FLUX_MESH)

    # M5.1 director-glyph overlay — line segments showing n̂ orientation,
    # signed-component RGB so opposite directions are visually opposite.
    # Renders on top of (or instead of) the flux mesh; fast (~768 segments).
    if state.SHOW_DIRECTORS > 0:
        glyph_length = 0.02  # ~2% of universe edge in normalized [0,1] coords
        arrow_length = (
            glyph_length * lagrange.ARROWHEAD_LENGTH_FRAC
            if lagrange.SHOW_DIRECTOR_ARROWHEAD
            else 0.0
        )
        lagrange.update_director_glyphs(
            state.wave_field, glyph_length, arrow_length, state.SHOW_DIRECTORS
        )
        render.scene.lines(
            state.wave_field.director_glyph_vertices,
            per_vertex_color=state.wave_field.director_glyph_colors,
            width=2.0,
        )
        # Half-arrowhead barb pass — toggled by lagrangian_engine.SHOW_DIRECTOR_ARROWHEAD.
        if lagrange.SHOW_DIRECTOR_ARROWHEAD:
            render.scene.lines(
                state.wave_field.director_glyph_arrow_vertices,
                per_vertex_color=state.wave_field.director_glyph_arrow_colors,
                width=2.0,
            )

    if state.PARTICLE_SHELL:
        # Convert wave-centers positions from [ijk] to [screen_normalization]
        # Use position_float for smooth rendering (position_grid is integer, causes jumpy motion)
        # Normalize by max_grid_size to respect asymmetric universes (like flux_mesh does)
        max_dim = float(state.wave_field.max_grid_size)
        for wc_idx in range(state.wave_center.num_sources):
            # Skip inactive (annihilated) WCs
            if state.wave_center.active[wc_idx] == 0:
                continue

            wc_pos_screen = ti.Vector(
                [
                    state.wave_center.position_float[wc_idx][0] / max_dim,
                    state.wave_center.position_float[wc_idx][1] / max_dim,
                    state.wave_center.position_float[wc_idx][2] / max_dim,
                ],
                dt=ti.f32,
            )
            position = np.array(
                [[wc_pos_screen[0], wc_pos_screen[1], wc_pos_screen[2]]], dtype=np.float32
            )
            # Particle shell radius — fixed-fraction-of-universe-edge default.
            # M5.2 replaces with per-defect-type Compton wavelength sizing.
            radius = 0.02  # ~2% of normalized universe edge
            color = (
                colormap.COLOR_PARTICLE[1]
                if state.SOURCES_OFFSET_DEG[wc_idx] == 180
                else colormap.COLOR_ANTI[1]
            )
            # Render particle shell at wave-center position
            render.scene.particles(position, radius, color=color)

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
        lagrange.sample_position_to_render(state.wave_field, amp_boost, stride, num_render)
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
    default_xperiment = selected_xperiment_arg or "_test4_topology"
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
