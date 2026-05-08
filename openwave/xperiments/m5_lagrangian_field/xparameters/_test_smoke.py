"""
XPERIMENT PARAMETERS — Plane Wave UI Smoke Test

This XPERIMENT seeds a transverse-polarized traveling plane wave along x and
lets propagate_psi (leapfrog) evolve it. The visual goal is to confirm:

- The wavefront propagates in the +x direction (not symmetrically)
- The wavelength matches the seed (12 voxels/λ)
- The phase speed visually tracks c (one λ traveled per period)
- No spurious checkerboard / dispersion until reflections arrive at the BC

Caveat: Dirichlet BCs at the universe edge will reflect the wave back ~halfway
through the simulation domain crossing time. The first ~50% of crossing time
is the clean regime for visual inspection.

Wave-centers are present (annihilation1 baseline) but APPLY_MOTION=False and
the seed overrides any wave-center contribution at t=0. The propagate_psi
kernel does not consult wave_center, so they stay decorative.
"""

UNIVERSE_EDGE = 1e-15  # m
TARGET_VOXELS = 100_000_000

# Wave seed parameters (consumed by _launcher when xperiment name == "_test_smoke")
WAVE_SEED = {
    "AMPLITUDE_AM": 5.0,  # peak displacement, am
    "VOXELS_PER_WAVELENGTH": 60.0,  # spatial period in voxels (≥12 for stable f32)
    "POLARIZATION": [0.0, 1.0, 0.0],  # transverse: ê_y for x-propagating wave
    "DIRECTION_AXIS": 0,  # 0=x, 1=y, 2=z
}

XPARAMETERS = {
    "meta": {
        "X_NAME": "Smoke Test",
        "DESCRIPTION": "Seeded traveling plane wave — visual leapfrog verification",
    },
    "camera": {
        "INITIAL_POSITION": [0.97, 1.44, 1.11],
    },
    "universe": {
        "SIZE": [UNIVERSE_EDGE, UNIVERSE_EDGE, UNIVERSE_EDGE],
        "TARGET_VOXELS": TARGET_VOXELS,
    },
    "wave_centers": {
        "COUNT": 1,
        "POSITION": [[0.50, 0.50, 0.50]],
        "PHASE_OFFSETS_DEG": [0],
        "INIT_VELOCITY": [[0.0, 0.0, 0.0]],
        "APPLY_MOTION": False,
    },
    "ui_defaults": {
        "SHOW_AXIS": False,
        "TICK_SPACING": 0.25,
        "SHOW_GRID": False,
        "SHOW_EDGES": True,
        "FLUX_MESH_PLANES": [0.5, 0.5, 0.5],
        "SHOW_FLUX_MESH": 3,
        "WARP_MESH": 40,
        "PARTICLE_SHELL": False,
        "SHOW_GRANULES": False,
        "SIM_SPEED": 1.0,
        "PAUSED": False,
    },
    "color_defaults": {
        "COLOR_THEME": "OCEAN",
        "WAVE_MENU": 1,
    },
    "analytics": {
        "INSTRUMENTATION": False,
        "EXPORT_VIDEO": False,
        "VIDEO_FRAMES": 24,
    },
    "wave_seed": WAVE_SEED,
}
