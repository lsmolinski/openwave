"""
XPERIMENT PARAMETERS

This XPERIMENT showcases:
- 10 WC electron in tetrahedral arrangement
- WCs placed at standing wave node separations (lambda/2)
- Tests lock-in stability with leapfrog integrator

STATUS: UNSTABLE — tetrahedral geometry is incompatible with uniform lambda/2
node lattice. 15/45 WC pairs sit at non-node distances (sqrt(3)*lam/2,
sqrt(2)*lam/2) creating destabilizing forces. Two paths to fix:
  1. Non-linear lambda(r) from Yee & Hauger shells (Phase 1d) — variable node
     spacing may accommodate the tetrahedral geometry
  2. M4 vector method (Phase 1c) — L/T spin conversion may create different
     equilibrium conditions than scalar M3
"""

import numpy as np
from openwave.common import constants

# Standing wave node spacing: lambda/2 is the first lock-in well
# WCs should sit at standing wave nodes of each other's fields
EWAVE_LENGTH = constants.EWAVE_LENGTH  # m, ~2.854e-17 m

UNIVERSE_EDGE = 1e-15  # m, universe edge length in meters
TARGET_VOXELS = 100_000_000  # Target voxel count (impacts performance)

# Lock spacing calibrated to lambda/2 (first standing wave node)
# In normalized coords: (lambda/2) / UNIVERSE_EDGE * grid_size / grid_size = lambda/2 / UNIVERSE_EDGE
# But we need it in normalized (0-1) coords based on grid_size
# grid_n ~ 464 for 100M voxels, dx ~ 2.16e-18 m
# lambda/2 = 1.427e-17 m, in voxels = 6.6, in normalized = 6.6/464 ~ 0.014
LOCK_SPACING = (EWAVE_LENGTH / 2) / UNIVERSE_EDGE  # ~0.0143, = lambda/2 in normalized coords


def tetrahedral_10(center=(0.5, 0.5, 0.5)):
    """Generate 10 positions in 1-3-6 tetrahedral arrangement."""
    layer_h = LOCK_SPACING * np.sqrt(2 / 3)  # tetrahedral layer spacing
    cx, cy, cz = center  # Reference center position

    Rb = LOCK_SPACING * 2 / np.sqrt(3)  # base vertex radius
    Rm = LOCK_SPACING / np.sqrt(3)  # midpoint/middle layer radius

    angles_v = np.radians([90, 210, 330])  # vertex angles
    angles_m = np.radians([30, 150, 270])  # midpoint angles

    return [
        # Base: 3 vertices + 3 midpoints
        *[[cx + Rb * np.cos(a), cy + Rb * np.sin(a), cz] for a in angles_v],
        *[[cx + Rm * np.cos(a), cy + Rm * np.sin(a), cz] for a in angles_m],
        # Middle: 3 elements
        *[[cx + Rm * np.cos(a), cy + Rm * np.sin(a), cz + layer_h] for a in angles_v],
        # Apex: 1 element
        [cx, cy, cz + 2 * layer_h],
    ]


XPARAMETERS = {
    "meta": {
        "X_NAME": f"  /Electron Test",
        "DESCRIPTION": "Energy Wave Charging, Propagation and Interaction",
    },
    "camera": {
        "INITIAL_POSITION": [0.27, 1.62, 0.90],  # [x, y, z] in normalized coordinates
    },
    "universe": {
        "SIZE": [UNIVERSE_EDGE, UNIVERSE_EDGE, UNIVERSE_EDGE],  # m, simulation domain [x, y, z]
        "TARGET_VOXELS": TARGET_VOXELS,  # Simulation voxel count (impacts performance)
    },
    "wave_centers": {
        "COUNT": 10,  # Number of wave-centers for this xperiment
        # Wave-Center positions: normalized coordinates (0-1 range, relative to max universe edge)
        "POSITION": tetrahedral_10(center=(0.60, 0.73, 0.50)),
        # Phase offsets for each wave-center (integer degrees, converted to radians internally)
        "PHASE_OFFSETS_DEG": [180, 180, 180, 180, 180, 180, 180, 180, 180, 180],
        "APPLY_MOTION": False,  # Toggle to apply motion at wave-centers, from force at each iteration
    },
    "ui_defaults": {
        "SHOW_AXIS": False,  # Toggle to show/hide axis lines
        "TICK_SPACING": 0.25,  # Axis tick marks spacing for position reference
        "SHOW_GRID": False,  # Toggle to show/hide the voxel data-grid
        "SHOW_EDGES": False,  # Toggle to show/hide universe edges
        "FLUX_MESH_PLANES": [0.5, 0.5, 0.5],  # [x, y, z] positions relative to universe size
        "SHOW_FLUX_MESH": 2,  # Flux Mesh toggle, 0: none, 1: xy, 2: xy+xz, 3: xy+xz+yz
        "WARP_MESH": 300,  # Visual warp mesh effect intensity
        "PARTICLE_SHELL": True,  # Toggle to enable/disable particle shell rendering
        "SHOW_GRANULES": False,  # Toggle to show/hide granule particles (rendered as points)
        "TIMESTEP": 2.0,  # Simulation timestep in rontoseconds (reduced for stability)
        "PAUSED": False,  # Pause/Start simulation toggle
    },
    "color_defaults": {
        "COLOR_THEME": "OCEAN",  # Choose color theme for rendering (OCEAN, DESERT, FOREST)
        "WAVE_MENU": 1,  # Check _launcher.py display_wave_menu() for wave_menu indexing
    },
    "analytics": {
        "INSTRUMENTATION": False,  # Toggle data acquisition and analytics
        "EXPORT_VIDEO": False,  # Toggle frame image export to video directory
        "VIDEO_FRAMES": 24,  # Target frame number to stop recording and finalize video export
    },
}
