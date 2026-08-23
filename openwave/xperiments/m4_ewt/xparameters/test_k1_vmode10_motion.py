# -*- coding: utf-8 -*-
"""
TEST: K=1, V_MODE=10 – check if WC moves.
This tests the motion engine with a single WC.
"""

UNIVERSE_EDGE = 5e-15
TARGET_VOXELS = 55_000_000

K = 1
POSITIONS = [[0.5, 0.5, 0.5]]
PHASES = [180]

XPARAMETERS = {
    "meta": {
        "X_NAME": "🧪 TEST: K=1 V_MODE=10 MOTION",
        "DESCRIPTION": "Single WC – check if motion works with V_MODE=10",
    },
    "camera": {"INITIAL_POSITION": [0.94, 0.91, 0.69]},
    "universe": {
        "SIZE": [UNIVERSE_EDGE, UNIVERSE_EDGE, UNIVERSE_EDGE],
        "TARGET_VOXELS": TARGET_VOXELS,
    },
    "wave_centers": {
        "COUNT": K,
        "POSITION": POSITIONS,
        "PHASE_OFFSETS_DEG": PHASES,
        "APPLY_MOTION": True,
    },
    "engine": {
        # ============================================================
        # P1: Base wave seed – VERY WEAK
        # ============================================================
        "SEED_MODE": 2,
        "SEED_BOOST": 0.0001,

        # ============================================================
        # P2: Non-linear potential – V_MODE=10 (Gaussian + quintic)
        # ============================================================
        "V_MODE": 10,
        "V_C1": -0.01,           # Very weak focusing
        "V_C2": 1.0,             # Strong saturation

        # ============================================================
        # P3: Wave-center interaction – SOFT DRIVE (adds asymmetry)
        # ============================================================
        "WC_INTERACT_MODE": 3,
        "WC_BOOST": 0.01,        # Weak drive – just enough to perturb
        "WC_RADIUS": 2,
        "WC_SIGMA": 1.5,

        # ============================================================
        # EMC density profile – PRESSURE ACTIVE
        # ============================================================
        "R_WALL": 100.0,
        "WALL_HEIGHT": 1.2,
        "DEFICIT_DEPTH": 0.9,
        "R_SOLITON": 35.0,
        "SIGMA": 3.0,
        "PRESSURE_STRENGTH": 0.001,   # Pressure force ON

        # ============================================================
        # Numerical stability
        # ============================================================
        "CFL_SAFETY": 0.01,
        "VELOCITY_DAMPING": 0.999,
    },
    "ui_defaults": {
        "SHOW_AXIS": False,
        "TICK_SPACING": 0.25,
        "SHOW_GRID": False,
        "SHOW_EDGES": False,
        "FLUX_MESH_PLANES": [0.5, 0.5, 0.5],
        "SHOW_FLUX_MESH": 1,
        "WARP_MESH": 30,
        "SHOW_GRANULES": False,
        "PARTICLE_SHELL": True,
        "SIM_SPEED": 0.25,
        "PAUSED": False,
    },
    "color_defaults": {"COLOR_THEME": "OCEAN", "WAVE_MENU": 1},
    "analytics": {
        "INSTRUMENTATION": True,
        "EXPORT_VIDEO": False,
        "VIDEO_FRAMES": 24,
    },
}