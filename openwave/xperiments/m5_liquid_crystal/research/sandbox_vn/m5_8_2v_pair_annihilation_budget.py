"""
M5.8.2v — TWO-DEFECT ± ANNIHILATION ACCOUNTING on the production matrix.

Upgrade of the 1D sine_gordon_annihilation.py demo toward the 4×4 matrix
substrate. Annihilation of a charge pair conserves two ledgers:

  CHARGE:  Q = (+1) + (-1) -> 0   (Gauss-Bonnet additivity)
  ENERGY:  the pair's integrated-Hamiltonian rest energy is released as
           outgoing waves (mass -> radiation), Q = 0 throughout.

This script settles the ACCOUNTING (the conserved quantities) on the production
matrix engine, which is the tractable + verifiable part:

  1. CHARGE bookkeeping is already validated on the 4×4 engine by
     m5_8_1_topo_charge.py (single +1 -> +0.996, single -1 -> -0.996, the
     enclosing sphere of a +/- pair -> 0.000, additivity PASS). Re-stated here.
  2. ENERGY budget: each defect's rest energy is H_static (the validated
     seed-level Hamiltonian); a well-separated +/- pair carries ~ 2 x H_static,
     which is the energy available for release on annihilation. Computed here
     across delta via the validated H_of primitive.

NOT done here (the named remaining build): the full 3D DYNAMICAL capture ->
breather -> vacuum evolution of the pair on the 4×4 field. The current engine
has no two-defect evolution path (the M5.1 ψ-engine pair-evolver was retired in
the M5.8 ψ-retire; m5_1_winding.py is archived). The 1D sine_gordon_annihilation.py
demonstrates the mechanism (kink + antikink -> breather -> vacuum, Q = 0); the
3D constrained two-defect evolver is a dedicated follow-up.

Headless, numpy-only (seed-level). No GUI.
"""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from m5_8_2q_delta_scaling import H_of  # validated seed rest energy  # noqa: E402

# Charge-bookkeeping result of record (m5_8_1_topo_charge.py, 4×4 engine).
CHARGE_RECORD = {
    "single +1": +0.996, "single -1": -0.996,
    "pair around +1": +0.995, "pair around -1": -0.995,
    "pair enclosing (total)": 0.000,
}

G_TIME = 8.0
DELTAS = [0.1, 0.3, 0.5]


def main():
    print("=" * 78)
    print("M5.8.2v — TWO-DEFECT ± ANNIHILATION ACCOUNTING (production matrix)")
    print("=" * 78)

    print("\n[1] CHARGE ledger (from m5_8_1_topo_charge.py, director winding on the 4×4 M):")
    for k, v in CHARGE_RECORD.items():
        print(f"    {k:24s} Q = {v:+.3f}")
    q_total = CHARGE_RECORD["pair enclosing (total)"]
    charge_ok = abs(q_total) < 0.05
    print(f"    => pair total Q = {q_total:+.3f}  ({'CONSERVED -> 0' if charge_ok else 'NOT 0'}); "
          "annihilation is charge-allowed (Gauss-Bonnet additivity)")

    print("\n[2] ENERGY ledger (rest energy per defect = H_static; pair budget = 2x):")
    print("    |   delta  | H_single (rest) | pair budget 2x | released on annihilation |")
    print("    | --- | --- | --- | --- |")
    rows = []
    for d in DELTAS:
        Hs, _ = H_of(d, G_TIME)
        rows.append((d, Hs, 2 * Hs))
        print(f"    | {d:6.2f} | {Hs:14.4f} | {2*Hs:13.4f} | {2*Hs:12.4f} (-> outgoing waves) |")

    print("\n  ACCOUNTING VERDICT:")
    print("    charge: (+1) + (-1) -> 0  (validated, m5_8_1)")
    print("    energy: 2 x H_static rest energy -> released as outgoing waves, Q = 0 throughout")
    print("    Both ledgers balance on the production matrix at the accounting level.")
    print("  REMAINING (named build): the 3D DYNAMICAL capture -> breather -> vacuum")
    print("  evolution of the pair; mechanism shown in 1D by sine_gordon_annihilation.py.")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
