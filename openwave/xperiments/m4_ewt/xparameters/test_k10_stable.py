# -*- coding: utf-8 -*-
"""
TEST: Stability of K=10 (1-3-6 geometry) with V_MODE=10
Safe starting configuration – amplitude should not explode.
"""

UNIVERSE_EDGE = 2e-15  # m, larger universe to avoid boundary reflections
TARGET_VOXELS = 55_000_000

from openwave.xperiments.m4_ewt.xparameters.utils.geometry import tetrahedron_10

K = 10
POSITIONS = tetrahedron_10(UNIVERSE_EDGE, center=(0.5, 0.5, 0.5))
PHASES = [180] * K

XPARAMETERS = {
    "meta": {
        "X_NAME": "🧪 TEST: K=10 STABLE (V_MODE=10)",
        "DESCRIPTION": "First stability test for K=10 – safe parameters",
    },
    "camera": {
        "INITIAL_POSITION": [0.94, 0.91, 0.69],
    },
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
        # P1: Base wave seed – VERY WEAK, so WC are the main source
        # ============================================================
        "SEED_MODE": 2,
        "SEED_BOOST": 0.0001,      # Extremely low to avoid dominating

        # ============================================================
        # P2: Non-linear potential – SAFE values
        # ============================================================
        "V_MODE": 0,             # Gaussian + quintic saturation
        "V_C1": -0.5,             # Weaker focusing (safer)
        "V_C2": 0.5,              # Strong saturation – prevents blow‑up

        # ============================================================
        # P3: Wave-center interaction – WEAK to avoid overheating
        # ============================================================
        "WC_INTERACT_MODE": 3,    # soft (additive)
        "WC_BOOST": 0.1,          # Weak drive
        "WC_RADIUS": 2,
        "WC_SIGMA": 1.5,

        # ============================================================
        # EMC density profile
        # ============================================================
        "R_WALL": 100.0,
        "WALL_HEIGHT": 1.2,
        "DEFICIT_DEPTH": 0.9,
        "R_SOLITON": 35.0,
        "SIGMA": 3.0,
        "PRESSURE_STRENGTH": 0.001,

        # ============================================================
        # Numerical stability – VERY CONSERVATIVE
        # ============================================================
        "CFL_SAFETY": 0.02,       # Very small → stable
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
        "PARTICLE_SHELL": False,
        "SIM_SPEED": 0.5,         # Slower → more stable
        "PAUSED": False,
    },
    "color_defaults": {
        "COLOR_THEME": "OCEAN",
        "WAVE_MENU": 1,
    },
    "analytics": {
        "INSTRUMENTATION": True,
        "EXPORT_VIDEO": False,
        "VIDEO_FRAMES": 24,
    },
}