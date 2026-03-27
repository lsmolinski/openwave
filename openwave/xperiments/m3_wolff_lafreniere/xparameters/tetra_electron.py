"""
XPERIMENT PARAMETERS

Validation suite for K=10 electron/positron stability.
Tests that the Combined Wolff-LaFreniere tetrahedron stability is real physics,
not a lucky grid alignment.

TESTS (switch by uncommenting the desired POSITION/PHASE block):
  1. ROTATED ELECTRON — same K=10 tetrahedron, rotated 45° on all axes
  2. POSITRON — same geometry, opposite phase (0° instead of 180°)
  3. TRANSLATED — moved to a different position in the universe
"""

import numpy as np
from openwave.common import constants

UNIVERSE_EDGE = 1e-15  # m, universe edge length in meters
TARGET_VOXELS = 100_000_000  # Target voxel count (impacts performance)


def tetrahedron_10(univ_edge, center=(0.5, 0.5, 0.5), rotation=(0, 0, 0)):
    """Generate 10 positions in 1-3-6 tetrahedral arrangement.

    Args:
        center: (x, y, z) center position in normalized coords
        rotation: (rx, ry, rz) rotation in degrees around x, y, z axes
    """
    # Lock spacing calibrated to lambda (first standing wave node)
    # In normalized coords: lambda / UNIVERSE_EDGE
    LOCK_SPACING = constants.EWAVE_LENGTH / univ_edge  # lambda in normalized coords

    layer_h = LOCK_SPACING * np.sqrt(2 / 3)  # tetrahedral layer spacing
    cx, cy, cz = center

    Rb = LOCK_SPACING * 2 / np.sqrt(3)  # base vertex radius
    Rm = LOCK_SPACING / np.sqrt(3)  # midpoint/middle layer radius

    angles_v = np.radians([90, 210, 330])  # vertex angles
    angles_m = np.radians([30, 150, 270])  # midpoint angles

    # Generate positions relative to center (0,0,0)
    local_positions = [
        # Base: 3 vertices + 3 midpoints
        *[[Rb * np.cos(a), Rb * np.sin(a), 0] for a in angles_v],
        *[[Rm * np.cos(a), Rm * np.sin(a), 0] for a in angles_m],
        # Middle: 3 elements
        *[[Rm * np.cos(a), Rm * np.sin(a), layer_h] for a in angles_v],
        # Apex: 1 element
        [0, 0, 2 * layer_h],
    ]

    # Apply rotation if specified
    rx, ry, rz = np.radians(rotation)
    Rx = np.array([[1, 0, 0], [0, np.cos(rx), -np.sin(rx)], [0, np.sin(rx), np.cos(rx)]])
    Ry = np.array([[np.cos(ry), 0, np.sin(ry)], [0, 1, 0], [-np.sin(ry), 0, np.cos(ry)]])
    Rz = np.array([[np.cos(rz), -np.sin(rz), 0], [np.sin(rz), np.cos(rz), 0], [0, 0, 1]])
    R = Rz @ Ry @ Rx

    return [[cx + p[0], cy + p[1], cz + p[2]] for p in (R @ np.array(local_positions).T).T]


# ── TEST 1: Rotated electron (45° on all axes) ──────────────────────────────
POSITIONS = tetrahedron_10(UNIVERSE_EDGE, center=(0.50, 0.50, 0.50), rotation=(0, 0, 0))
# POSITIONS = tetrahedron_10(UNIVERSE_EDGE, center=(0.35, 0.65, 0.45), rotation=(0, 0, 0))
PHASES = [180] * 10  # electron (all same phase = 180°)

XPARAMETERS = {
    "meta": {
        "X_NAME": f"  /Tetrahedron Electron",
        "DESCRIPTION": "K=10 tetrahedron stability validation",
    },
    "camera": {
        "INITIAL_POSITION": [0.36, 1.20, 0.75],  # [x, y, z] in normalized coordinates
    },
    "universe": {
        "SIZE": [UNIVERSE_EDGE, UNIVERSE_EDGE, UNIVERSE_EDGE],  # m, simulation domain [x, y, z]
        "TARGET_VOXELS": TARGET_VOXELS,  # Simulation voxel count (impacts performance)
    },
    "wave_centers": {
        "COUNT": 10,  # Number of wave-centers for this xperiment
        "POSITION": POSITIONS,
        "PHASE_OFFSETS_DEG": PHASES,
        "APPLY_MOTION": True,
    },
    "ui_defaults": {
        "SHOW_AXIS": False,
        "TICK_SPACING": 0.25,
        "SHOW_GRID": False,
        "SHOW_EDGES": False,
        "FLUX_MESH_PLANES": [0.5, 0.5, 0.5],
        "SHOW_FLUX_MESH": 2,
        "WARP_MESH": 150,
        "PARTICLE_SHELL": True,
        "TIMESTEP": 2.0,
        "PAUSED": False,
    },
    "color_defaults": {
        "COLOR_THEME": "OCEAN",
        "WAVE_MENU": 4,
    },
    "analytics": {
        "INSTRUMENTATION": False,
        "EXPORT_VIDEO": False,
        "VIDEO_FRAMES": 24,
    },
}
