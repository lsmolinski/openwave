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
from openwave.xperiments.m3_wolff_lafreniere.xparameters.tetra_electron import tetrahedron_10

UNIVERSE_EDGE = 1e-15  # m, universe edge length in meters
TARGET_VOXELS = 100_000_000  # Target voxel count (impacts performance)

# ── TEST 1: Rotated electron (45° on all axes) ──────────────────────────────
POSITIONS = tetrahedron_10(univ_edge=UNIVERSE_EDGE, center=(0.50, 0.50, 0.50), rotation=(0, 0, 0))
PHASES = [180] * 10  # electron (all same phase = 180°)

# ── TEST 2: Positron (same geometry, phase = 0°) ────────────────────────────
# POSITIONS = tetrahedron_10(univ_edge=UNIVERSE_EDGE, center=(0.5, 0.5, 0.5), rotation=(45, 45, 45))
# PHASES = [0] * 10  # positron (all same phase = 0°)

# ── TEST 3: Translated (off-center, different grid alignment) ────────────────
# POSITIONS = tetrahedron_10(univ_edge=UNIVERSE_EDGE, center=(0.35, 0.65, 0.45), rotation=(0, 0, 0))
# PHASES = [180] * 10

XPARAMETERS = {
    "meta": {
        "X_NAME": f"  /Tetrahedrons",
        "DESCRIPTION": "K=10 tetrahedron stability validation",
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
