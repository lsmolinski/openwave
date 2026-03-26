"""
XPERIMENT PARAMETERS

This XPERIMENT showcases:
- Opposite-phase WC annihilation (head-on, close range)
- WCs placed inside the first sinc barrier (separation < λ/2)
- The deepest well at r=0 pulls them together → wave cancellation

PHYSICS: For opposite-phase WCs, the energy landscape has:
  well at r=0 (deepest), barrier at λ/2, well at λ, barrier at 3λ/2 ...
  Starting inside the first barrier means the gradient naturally pulls to r=0.
  Far-field barriers (from phasor energy ∝ A_self × ∇A_other) are too strong
  to cross — real pair annihilation requires particles created at close range.
"""

UNIVERSE_EDGE = 1e-15  # m, universe edge length in meters
TARGET_VOXELS = 100_000_000  # Target voxel count (impacts performance)

XPARAMETERS = {
    "meta": {
        "X_NAME": f"Particle Annihilation",
        "DESCRIPTION": "Opposite-phase WC annihilation",
    },
    "camera": {
        "INITIAL_POSITION": [0.37, 1.45, 1.24],  # [x, y, z] in normalized coordinates
    },
    "universe": {
        "SIZE": [UNIVERSE_EDGE, UNIVERSE_EDGE, UNIVERSE_EDGE],  # m, simulation domain [x, y, z]
        "TARGET_VOXELS": TARGET_VOXELS,  # Simulation voxel count (impacts performance)
    },
    "wave_centers": {
        "COUNT": 2,  # Number of wave-centers for this xperiment
        # Wave-Center positions: normalized coordinates (0-1 range, relative to max universe edge)
        "POSITION": [
            [0.15, 0.50, 0.50],
            [0.85, 0.50, 0.50],
        ],
        # Phase offsets for each wave-center (integer degrees, converted to radians internally)
        "PHASE_OFFSETS_DEG": [0, 180],
        # Initial velocity [vx, vy, vz] in am/rs (c = 0.3 am/rs)
        "INIT_VELOCITY": [
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
        ],
        "APPLY_MOTION": True,  # Toggle to apply motion at wave-centers, from force at each iteration
    },
    "ui_defaults": {
        "SHOW_AXIS": False,  # Toggle to show/hide axis lines
        "TICK_SPACING": 0.25,  # Axis tick marks spacing for position reference
        "SHOW_GRID": False,  # Toggle to show/hide the voxel data-grid
        "SHOW_EDGES": False,  # Toggle to show/hide universe edges
        "FLUX_MESH_PLANES": [0.5, 0.5, 0.5],  # [x, y, z] positions relative to universe size
        "SHOW_FLUX_MESH": 1,  # Flux Mesh toggle, 0: none, 1: xy, 2: xy+xz, 3: xy+xz+yz
        "WARP_MESH": 500,  # Visual warp mesh effect intensity
        "PARTICLE_SHELL": True,  # Toggle to enable/disable particle shell rendering
        "TIMESTEP": 5.0,  # Simulation timestep in rontoseconds
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
