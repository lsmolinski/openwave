"""
XPERIMENT PARAMETERS

Two K=10 tetrahedra interaction tests.
Tests particle-particle behavior at distance.

TESTS (switch by uncommenting the desired block):
  1. ELECTRON + ELECTRON — same phase (180°+180°), expect repulsion
  2. ELECTRON + POSITRON — opposite phase (180°+0°), expect attraction
  3. POSITRON + POSITRON — same phase (0°+0°), expect repulsion
"""

from openwave.xperiments.m3_wolff_lafreniere.xparameters.formation02 import tetrahedron_10

UNIVERSE_EDGE = 1e-15  # m, universe edge length in meters
TARGET_VOXELS = 100_000_000  # Target voxel count (impacts performance)

# Separation between the two tetrahedra centers (~10λ apart along x-axis)
# Each tetrahedron has radius ~1λ, so 10λ gives ~8λ clear space between them
SEP_X = 0.30  # normalized coords (~10λ at 1e-15 universe)

# ── TEST 1: Electron + Electron (same phase → expect repulsion) ─────────────
POSITIONS_A = tetrahedron_10(UNIVERSE_EDGE, (0.5 - SEP_X / 2, 0.50, 0.50), (0, 0, 0))
POSITIONS_B = tetrahedron_10(UNIVERSE_EDGE, (0.5 + SEP_X / 2, 0.50, 0.50), (60, 30, 0))
PHASES_A = [180] * 10  # electron
PHASES_B = [180] * 10  # electron

# ── TEST 2: Electron + Positron (opposite phase → expect attraction) ─────────
# POSITIONS_A = tetrahedron_10(UNIVERSE_EDGE, (0.5 - SEP_X / 2, 0.50, 0.50), (0, 0, 0))
# POSITIONS_B = tetrahedron_10(UNIVERSE_EDGE, (0.5 + SEP_X / 2, 0.50, 0.50), (60, 30, 0))
# PHASES_A = [180] * 10  # electron
# PHASES_B = [0] * 10  # positron

# ── TEST 3: Positron + Positron (same phase → expect repulsion) ──────────────
# POSITIONS_A = tetrahedron_10(UNIVERSE_EDGE, (0.5 - SEP_X / 2, 0.50, 0.50), (0, 0, 0))
# POSITIONS_B = tetrahedron_10(UNIVERSE_EDGE, (0.5 + SEP_X / 2, 0.50, 0.50), (60, 30, 0))
# PHASES_A = [0] * 10  # positron
# PHASES_B = [0] * 10  # positron

# Combine both tetrahedra into a single WC list
POSITIONS = POSITIONS_A + POSITIONS_B
PHASES = PHASES_A + PHASES_B

XPARAMETERS = {
    "meta": {
        "X_NAME": f"  /2 Tetrahedra (WIP)",
        "DESCRIPTION": "Two K=10 tetrahedra interaction",
    },
    "camera": {
        "INITIAL_POSITION": [0.27, 1.62, 0.90],
    },
    "universe": {
        "SIZE": [UNIVERSE_EDGE, UNIVERSE_EDGE, UNIVERSE_EDGE],
        "TARGET_VOXELS": TARGET_VOXELS,
    },
    "wave_centers": {
        "COUNT": 20,  # 10 + 10 = two tetrahedra
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
