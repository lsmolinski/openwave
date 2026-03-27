"""
XPERIMENT PARAMETERS

K=2 through K=10 progressive particle formation test.
Verifies EWT prediction: K=2..9 are UNSTABLE (decay/fly apart),
K=10 is the first stable standalone particle (electron tetrahedron).

Switch K by uncommenting the desired configuration below.
Each K places WCs at λ/2 spacing in the simplest possible geometry:
  K=2:  line (2 WCs along x-axis)
  K=3:  equilateral triangle
  K=4:  regular tetrahedron (4 vertices)
  K=5:  trigonal bipyramid
  K=6:  octahedron
  K=7:  pentagonal bipyramid
  K=8:  cube (dual tetrahedra — should be neutral/unstable differently)
  K=9:  tricapped trigonal prism
  K=10: 1-3-6 tetrahedron (STABLE — the electron)
"""

import numpy as np
from openwave.common import constants

UNIVERSE_EDGE = 1e-15  # m, universe edge length in meters
TARGET_VOXELS = 100_000_000  # Target voxel count (impacts performance)

# Lock spacing: λ in normalized coords (first lock-in well for Combined W-L)
# The envelope 2|sin(kr/2)|/r has zeros at r=nλ — these are the energy minima
# where same-phase WCs lock in. (Wolff-original sinc had wells at λ/2.)
EWAVE_LENGTH = constants.EWAVE_LENGTH
LOCK_SPACING = EWAVE_LENGTH / UNIVERSE_EDGE

# ════════════════════════════════════════════════════════════════════════════
# SELECT K VALUE HERE — uncomment one line
# ════════════════════════════════════════════════════════════════════════════
K = 2
# K = 2    # Line — EXPECT: unstable
# K = 3    # Triangle — EXPECT: unstable
# K = 4    # Tetrahedron (4) — EXPECT: unstable
# K = 5    # Trigonal bipyramid — EXPECT: unstable
# K = 6    # Octahedron — EXPECT: unstable
# K = 7    # Pentagonal bipyramid — EXPECT: unstable
# K = 8    # Cube (dual tetra) — EXPECT: unstable (neutral in EWT)
# K = 9    # Tricapped prism — EXPECT: unstable
# K = 10   # 1-3-6 tetrahedron — EXPECT: STABLE (electron)


def generate_K_positions(K, center=(0.5, 0.5, 0.5)):
    """Generate K WC positions in a compact geometry at λ/2 spacing.

    Uses simple symmetric arrangements. Not all are the natural EWT geometry
    (only K=10 has the 1-3-6 tetrahedron), but they test whether arbitrary
    K-body arrangements are stable or decay.
    """
    cx, cy, cz = center
    s = LOCK_SPACING  # shorthand: λ/2 in normalized coords

    if K == 2:
        # Line along x
        return [[cx - s / 2, cy, cz], [cx + s / 2, cy, cz]]

    elif K == 3:
        # Equilateral triangle in xy plane
        angles = np.radians([90, 210, 330])
        r = s / np.sqrt(3)  # circumradius for edge = s
        return [[cx + r * np.cos(a), cy + r * np.sin(a), cz] for a in angles]

    elif K == 4:
        # Regular tetrahedron (4 vertices)
        # Vertices of a regular tetrahedron with edge length s
        h = s * np.sqrt(2 / 3)
        r = s / np.sqrt(3)
        return [
            [cx, cy + r, cz],
            [cx - s / 2, cy - r / 2, cz],
            [cx + s / 2, cy - r / 2, cz],
            [cx, cy, cz + h],
        ]

    elif K == 5:
        # Trigonal bipyramid: triangle + 2 poles
        angles = np.radians([90, 210, 330])
        r = s / np.sqrt(3)
        positions = [[cx + r * np.cos(a), cy + r * np.sin(a), cz] for a in angles]
        positions.append([cx, cy, cz + s / 2])  # top pole
        positions.append([cx, cy, cz - s / 2])  # bottom pole
        return positions

    elif K == 6:
        # Octahedron: 6 vertices along axes
        d = s / np.sqrt(2)
        return [
            [cx + d, cy, cz],
            [cx - d, cy, cz],
            [cx, cy + d, cz],
            [cx, cy - d, cz],
            [cx, cy, cz + d],
            [cx, cy, cz - d],
        ]

    elif K == 7:
        # Pentagonal bipyramid: pentagon + 2 poles
        angles = np.radians([i * 72 for i in range(5)])
        r = s / (2 * np.sin(np.pi / 5))  # circumradius for edge = s
        positions = [[cx + r * np.cos(a), cy + r * np.sin(a), cz] for a in angles]
        positions.append([cx, cy, cz + s / 2])
        positions.append([cx, cy, cz - s / 2])
        return positions

    elif K == 8:
        # Cube: 8 vertices (dual tetrahedra in EWT — K=8 muon neutrino)
        d = s / 2
        return [
            [cx + d, cy + d, cz + d],
            [cx + d, cy + d, cz - d],
            [cx + d, cy - d, cz + d],
            [cx + d, cy - d, cz - d],
            [cx - d, cy + d, cz + d],
            [cx - d, cy + d, cz - d],
            [cx - d, cy - d, cz + d],
            [cx - d, cy - d, cz - d],
        ]

    elif K == 9:
        # Tricapped trigonal prism: 6 (trigonal prism) + 3 caps
        angles = np.radians([90, 210, 330])
        r = s / np.sqrt(3)
        h = s / 2
        # Top and bottom triangles
        positions = []
        for z_off in [-h, h]:
            for a in angles:
                positions.append([cx + r * np.cos(a), cy + r * np.sin(a), cz + z_off])
        # 3 equatorial caps (midpoints of prism edges, pushed outward)
        cap_angles = np.radians([30, 150, 270])
        cap_r = r * 1.5
        for a in cap_angles:
            positions.append([cx + cap_r * np.cos(a), cy + cap_r * np.sin(a), cz])
        return positions

    elif K == 10:
        # 1-3-6 tetrahedron (from tetra_electron.py)
        from openwave.xperiments.m3_wolff_lafreniere.xparameters.tetra_electron import (
            tetrahedron_10,
        )

        return tetrahedron_10(center=center)

    else:
        raise ValueError(f"K={K} not implemented")


POSITIONS = generate_K_positions(K)
PHASES = [180] * K  # all same phase (electron-like)

XPARAMETERS = {
    "meta": {
        "X_NAME": f"Particle Formation K={K}",
        "DESCRIPTION": f"K={K} stability test — {'STABLE' if K == 10 else 'expect UNSTABLE'}",
    },
    "camera": {
        "INITIAL_POSITION": [0.36, 1.20, 0.75],
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
