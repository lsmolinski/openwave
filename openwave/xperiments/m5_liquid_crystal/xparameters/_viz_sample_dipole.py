"""
XPERIMENT PARAMETERS — Magnetic-dipole VIZ SAMPLE (VIZ.4 / M5.6.5f stage 1)

A RENDER UNIT-TEST, not physics. A static hedgehog is a pure ELECTRIC charge
(∇×n̂ ≈ 0) — no circulating B, no poles — so the bluered N/S coloring (VIZ.2
WM7) and the B-field glyphs (VIZ.3 state 3) have nothing to show until a
twisting/spinning defect generates a real circulating B (the Zitterbewegung
clock, M5.8). This xperiment overwrites the curl observable with an ideal
ANALYTIC dipole field

    B(r) = amp·[3(m̂·r̂)r̂ − m̂] / max(r, r0)³        (m̂ = DIPOLE_AXIS = +ẑ)

so the render path can be built + visually approved NOW. At M5.8 the same render
points at the real circulating B with no render change (the DIPOLE_SAMPLE flag
goes away). See 4b §4.5 (sample-first) + §5.3 (placeholder strategy).

WHAT YOU SHOULD SEE (boots paused, WAVE_MENU 7 = EM curl, bluered N/S on):
  • A bar-magnet N/S: RED top hemisphere (N — B flows OUT) + BLUE bottom (S — B
    flows IN), white at the equator. This is the RADIAL projection (∇×n̂)·r̂ ∝ cosθ
    (the dipole sample sets curl_radial=1) — the coloring that actually shows poles,
    matching Duda's slide. (The fixed-axis projection (∇×n̂)·ẑ would instead light
    RED at BOTH ±ẑ ends + a blue equatorial belt — the field's axial component, real
    but not a bar magnet — which is why the dipole sample uses radial.)
  • The B-field glyphs (Magnetic Field state) trace the closing dipole field lines.
  • A YELLOW magnetic-MOMENT vector glyph μ marks the moment axis at the center.

A biaxial hedgehog is still seeded underneath so the ∇·n̂ charge (WM6 / Electric
Field glyphs) shows real structure for context — only the curl/B channel is the
placeholder.

Repo discipline (4b §5.3): the placeholder is a physics-shaped analytic function
(fine in open-source OpenWave); it validates the renderer, which is open infra.
"""

UNIVERSE_EDGE = 1e-15  # m
TARGET_VOXELS = 64**3  # ~262k voxels

TOPOLOGY_SEED = {
    "MODE": "biaxial_hedgehog",
    "CENTER": [0.50, 0.50, 0.50],
    "R0_FRACTION": 0.06,
    "RHOC_VOXELS": 3.0,
    "BIAXIAL_DELTA": 0.30,
    "AUTO_RELAX_STEPS": 0,
    "LDG_STIFFNESS_K": 0.0,  # V off — this is a render test, no dynamics needed
    # --- VIZ.4 dipole placeholder ---
    "DIPOLE_SAMPLE": True,
    "DIPOLE_AXIS": [0.0, 0.0, 1.0],  # moment direction m̂ (also the bluered projection axis)
    "DIPOLE_R0_VOX": 3.0,  # core regularization radius (voxels) — B saturates inside
}


XPARAMETERS = {
    "meta": {
        "X_NAME": "VIZ Magnetic Dipole",
        "DESCRIPTION": "Analytic dipole B placeholder — validate N/S coloring + B glyphs + moment glyph (VIZ.4)",
    },
    "camera": {
        "INITIAL_POSITION": [1.10, 1.46, 0.81],
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
        "SHOW_EDGES": False,
        "FLUX_MESH_PLANES": [0.5, 0.5, 0.5],
        "SHOW_FLUX_MESH": 3,  # all 3 planes — see the full dipole pattern
        "WARP_MESH": 0,  # no vector-warp (raw 1/r³ field warps wildly near core); color-only
        "SHOW_GLYPHS": 3,  # all 3 planes
        "GLYPH_VECTOR": 3,  # Magnetic Field (∇×n̂) glyphs — trace the dipole field lines
        "GLYPH_SIZE": 1,  # magnitude-scaled → declutters the far field where B→0
        "GLYPH_COLOR": 1,  # gradient color
        "VIZ_STRIDE": 2,
        "SHOW_GRANULES": False,
        "SIM_SPEED": 1.0,
        "PAUSED": True,
    },
    "color_defaults": {
        "COLOR_THEME": "OCEAN",
        "WAVE_MENU": 7,  # EM curl (B) — the dipole channel
        "CURL_COLOR": 1,  # bluered signed (∇×n̂)·axis → N=red / S=blue
    },
    "analytics": {
        "INSTRUMENTATION": False,
        "EXPORT_VIDEO": False,
        "VIDEO_FRAMES": 24,
    },
    "topology_seed": TOPOLOGY_SEED,
}
