# -*- coding: utf-8 -*-
"""
TEST: K=1, V_MODE=10 with FIXED SIGN in V_C2.
This should be stable and show WC motion.
"""

UNIVERSE_EDGE = 5e-15
TARGET_VOXELS = 55_000_000

K = 10
POSITIONS = [[0.5, 0.5, 0.5]]
PHASES = [180]

XPARAMETERS = {
    "meta": {
        "X_NAME": "🧪 TEST: K=1 V_MODE=10 FIXED SIGN",
        "DESCRIPTION": "Fixed V_C2 sign – should be stable",
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
        # P1: Base wave seed – moderate
        # ============================================================
        "SEED_MODE": 2,
        "SEED_BOOST": 0.001,

        # ============================================================
        # P2: V_MODE=10 with FIXED SIGN
        # ============================================================
        "V_MODE": 10,
        "V_C1": -0.5,            # Focusing 
        "V_C2": 0.1,            # POSITIVE – saturation (with fixed formula)


        # ============================================================
        # P3: Soft drive – reduced boost (dt is 100× smaller)
        # ============================================================
        "WC_INTERACT_MODE": 0,
        "WC_BOOST": 0,     
        "WC_RADIUS": 2,
        "WC_SIGMA": 1.5,

        # ============================================================
        # EMC density profile
        # ============================================================
        "R_WALL": 100.0,
        "WALL_HEIGHT": 1.2,
        "DEFICIT_DEPTH": 1.0,
        "R_SOLITON": 50.0,
        "SIGMA": 3.0,
        "PRESSURE_STRENGTH": 0.05,

        # ============================================================
        # Numerical stability
        # ============================================================
        "CFL_SAFETY": 0.01,
        "VELOCITY_DAMPING": 0.99,
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
        "PARTICLE_SHELL": False,
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