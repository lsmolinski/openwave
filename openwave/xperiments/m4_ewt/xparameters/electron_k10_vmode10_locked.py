"""
electron_k10_vmode10_locked.py
Configuration for Electron K=10 in V_MODE=10 (Gaussian profile + quintic saturation)
using locked EWT geometry (1-3-6 layout).
"""

from openwave.xperiments.m4_ewt.xparameters.utils.geometry import (
    generate_positions_by_EWT_geometry_locked,
)

# ================================================================
# V_MODE=10: Gaussian profile + quintic saturation (K=10 Locked)
# ================================================================
UNIVERSE_EDGE = 2e-15
TARGET_VOXELS = 80_000_000
K = 10
PERTURBATION = 0.01

POSITIONS = generate_positions_by_EWT_geometry_locked(
    UNIVERSE_EDGE, K, center=(0.5, 0.5, 0.5), perturbation=PERTURBATION
)
PHASES = [180] * K

XPARAMETERS = {
    "meta": {
        "X_NAME": "Electron K=10 V_MODE=10",
        "DESCRIPTION": "K=10 locked 1-3-6 geometry with V_MODE=10 (Gaussian + quintic saturation)",
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
        "SEED_MODE": 2,
        "SEED_BOOST": 0.00000001,
        "V_MODE": 10,  
        "V_C1": -1.15,  # focusing (cubic)
        "V_C2": 0.05,  # saturation (quintic)
        "WC_INTERACT_MODE": 3,
        "WC_BOOST": 0.0005,
        "WC_RADIUS": 5,
        "WC_SIGMA": 1.44,
        "R_WALL": 100.0,
        "WALL_HEIGHT": 1.2,
        "DEFICIT_DEPTH": 0.8,
        "R_SOLITON": 45.0,
        "SIGMA": 1.1,
        "PRESSURE_STRENGTH": 0.1,  # vacuum pressure (active for V_MODE >= 4)
        "CFL_SAFETY": 0.005,
        "VELOCITY_DAMPING": 0.99999
    },
    "ui_defaults": {
        "SHOW_AXIS": False,
        "TICK_SPACING": 0.25,
        "SHOW_GRID": False,
        "SHOW_EDGES": False,
        "FLUX_MESH_PLANES": [0.5, 0.5, 0.5],
        "SHOW_FLUX_MESH": 1,
        "WARP_MESH": 100,
        "SHOW_GRANULES": False,
        "PARTICLE_SHELL": False,
        "SIM_SPEED": 1.0,
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