"""Pytest configuration for FORCE_UNIFICATION tests."""

import sys
from pathlib import Path

# Add this directory to sys.path so wave_engine_1D_v2 is importable
sys.path.insert(0, str(Path(__file__).parent))
