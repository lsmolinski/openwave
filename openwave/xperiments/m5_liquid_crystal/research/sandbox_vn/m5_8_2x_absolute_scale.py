"""
M5.8.2x - Absolute-scale calibration (the #208 deliverable: Phase A + B)

Builds the sim -> physical unit map for the M5 liquid-crystal substrate and
computes the dimensionless invariants that test it (the fine-structure
constant alpha, the ZBW ratio, the apolar doubling), plus the ~28x clock gap.

This is a PURE-ANALYSIS script: it takes the measured sim-unit numbers as
inputs (each sourced inline) and propagates them. No engine run here; the
heavy measurements live in the cited scripts.

FRAMING (per Duda's input on #208): the model is SCALE-FREE, so the absolute
scale is not "discovered", it is DEFINED by anchoring the unit map on ONE
observable, after which every other observable becomes a falsifiable physical
number. With c emergent (radial cone, N-4) and hbar fixed by postulate P2
(lattice action unit == hbar), the map has exactly TWO free dials: a length
scale and an energy scale. Fixing those two from one anchor closes the map.

KEY DISTINCTION:
  - dimensionless invariants (alpha, omega_M/omega_clock, the ZBW ratio) are
    ANCHOR-INDEPENDENT -> they are the real tests, computed first.
  - dimensionful predictions (absolute omega in Hz, etc.) are set by the
    chosen anchor -> the ~28x gap is one such number.

MEASURED INPUTS (sources):
  Coulomb coupling  b  = -5.633e-3   m5_4_coulomb_matrix.py (R2=0.9781; Frank-
                                      elastic energy x voxel; dx=15.385 am, delta=0.5)
  rest energy  H_static = 16.74      m5_8_2j_zbw_ratio.py (matrix Hamiltonian
                                      u + beta u^2, beta=1.558, R_W=3.5)
  clock omega1 = 1.188               m5_8_2j / m5_8_2h (median detrend, start-
                                      independent across arms at 5.4%)
  apolar doubling omega_M/omega1 = 2 machine-exact, m5_8_2s_spin_half_apolar.py
  Faber  r0 = 2.2132 fm -> 0.511 MeV m5_6_3a/3b (E*r0=const; Faber Eq.8
                                      E0 = alpha * hbar c * pi / (4 r0))
"""

import numpy as np

# ---------- physical constants ----------
MEV_J = 0.511e6 * 1.602176634e-19          # electron rest energy m_e c^2 in J
HBAR = 1.054571817e-34                      # J s
ALPHA_EXP = 7.2973525693e-3                 # measured fine-structure constant
ALPHA_INV_EXP = 137.035999084

# ---------- measured sim-unit inputs ----------
B_COUL = -5.633e-3      # Coulomb 1/d coefficient; |b| = e^2/(4 pi eps0) in [energy x length]
H_STATIC = 16.74        # rest energy (lattice Hamiltonian units)
OMEGA1 = 1.188          # clock fundamental (rad / lattice-time)


def banner(title):
    print("=" * 74)
    print(title)
    print("=" * 74)


def part_a_invariants():
    """Anchor-INDEPENDENT dimensionless numbers (the real tests)."""
    banner("PART A - dimensionless invariants (anchor-independent)")

    # (A1) fine-structure constant.
    #   alpha = (e^2/4pi eps0) / (hbar c).
    #   lattice convention: hbar = 1 (action unit, P2); c = 1 (radial cone
    #   voxel/dt, N-4 m5_8_2k); |b| = e^2/(4 pi eps0) in [energy x length].
    #   => alpha = |b| / (hbar c).
    # The remaining uncertainty is an O(1) unit-reconciliation factor: the
    # Coulomb |b| is measured in the FRANK-elastic functional while hbar is
    # set on the matrix Hamiltonian, and the radial cone carries a c/2
    # convention. We bracket c in {0.5, 1, 2} to expose that sensitivity.
    print("(A1) fine-structure constant  alpha = |b| / (hbar c),  hbar=1, |b|=%.3e" % abs(B_COUL))
    print("     %-10s %-14s %-14s %-10s" % ("c (lattice)", "alpha", "1/alpha", "vs 137.036"))
    for c_lat in (0.5, 1.0, 2.0):
        a = abs(B_COUL) / (1.0 * c_lat)
        print("     %-10.1f %-14.4e %-14.1f x%.2f" % (c_lat, a, 1.0 / a, (1.0 / a) / ALPHA_INV_EXP))
    a_central = abs(B_COUL) / 1.0
    print("     central (c=1): 1/alpha = %.1f   vs measured %.3f  (factor %.2f)" %
          (1.0 / a_central, ALPHA_INV_EXP, (1.0 / a_central) / ALPHA_INV_EXP))
    print("     -> a dimensionless EM coupling of the RIGHT ORDER from pure topology,")
    print("        not tuned. Precision awaits the consistent-units re-run (Frank vs")
    print("        Hamiltonian energy functional; the c-cone factor).")
    print()

    # (A2) ZBW ratio - the clean, start-independent clock invariant.
    zbw_ratio = OMEGA1 / (2.0 * H_STATIC)
    print("(A2) ZBW ratio  omega1 / (2 H_static) = %.4f" % zbw_ratio)
    print("     start-independent across arms (5.4%% spread, m5_8_2j) -> unit-free,")
    print("     inherits the N-1 attractor. THIS is the rock-solid clock invariant.")
    print()

    # (A3) apolar doubling - machine-exact.
    print("(A3) apolar doubling  omega_M / omega_clock = 2.0  (machine-exact, m5_8_2s)")
    print("     the Dirac ZBW factor-of-2; an exact invariant, no calibration needed.")
    print()
    return a_central, zbw_ratio


def part_b_unit_map():
    """The dimensionful predictions per anchor; the ~28x clock gap."""
    banner("PART B - the unit map per anchor + the ~28x clock gap")

    # Energy anchor (P1): H_static = 16.74 lattice <-> m_e c^2 = 0.511 MeV
    E_per_lat_J = MEV_J / H_STATIC
    E_per_lat_MeV = 0.511 / H_STATIC
    print("Energy anchor (P1):  H_static = %.2f lattice  <->  m_e c^2 = 0.511 MeV" % H_STATIC)
    print("   => 1 lattice energy unit = %.4f MeV = %.3e J" % (E_per_lat_MeV, E_per_lat_J))
    print()

    # Clock prediction under P1 + P2 (action == hbar):
    #   omega_phys = omega1 * (m_e c^2 / hbar) / H_static
    om_phys = OMEGA1 * (MEV_J / HBAR) / H_STATIC
    zbw_target = 2.0 * (MEV_J / HBAR)
    gap = zbw_target / om_phys
    print("Clock prediction (P1 + P2):  omega_phys = omega1 * (m_e c^2/hbar) / H_static")
    print("   omega_phys   = %.3e rad/s" % om_phys)
    print("   ZBW target   = 2 m_e c^2/hbar = %.3e rad/s" % zbw_target)
    print("   GAP          = %.1fx  (the clock runs ~28x BELOW the electron ZBW)" % gap)
    print()

    print("Anchor map (what each FIXES, then PREDICTS):")
    rows = [
        ("Coulomb (e-scale)", "charge x length x energy", "alpha (dimensionless) -> 1/alpha ~ %.0f" % (1.0 / abs(B_COUL))),
        ("Electron clock", "time/energy via the ZBW", "the mass; ratio omega1/(2H)=%.3f" % (OMEGA1 / (2 * H_STATIC))),
        ("Faber mass (r0)", "r0=2.2132 fm <-> 0.511 MeV", "clock omega -> lands %.0fx low (the gap)" % gap),
    ]
    for name, fixes, predicts in rows:
        print("   %-20s fixes: %-26s predicts: %s" % (name, fixes, predicts))
    print()
    return gap


def structural_findings(gap):
    banner("STRUCTURAL FINDINGS (the calibration's real conclusions)")
    print("1. The anchors are NOT all independent. Faber's rest-energy formula")
    print("   E0 = alpha * hbar c * pi / (4 r0) CONTAINS alpha as an INPUT, so")
    print("   anchoring the mass via Faber presupposes alpha. The Coulomb route")
    print("   is the independent alpha measurement; the clock is the independent")
    print("   absolute-FREQUENCY test -> and it is the one that fails (~%.0fx)." % gap)
    print()
    print("2. The ~28x gap is the model's ONE genuine absolute-scale discrepancy.")
    print("   It is dimensionful and structural (omega is rigid against energy,")
    print("   mass-family independent, N-6a). It does NOT close by bookkeeping.")
    print()
    print("3. alpha lands the RIGHT ORDER (1/alpha ~ 178 vs 137) from pure")
    print("   topology, un-tuned. The residual factor is exactly the consistent-")
    print("   units reconciliation (energy functional + c-cone), which is the")
    print("   re-run the DoD specifies -> Phase A is a real, falsifiable result,")
    print("   not a fit.")
    print()
    print("4. VERDICT: a working 2-dial unit map exists per anchor; the honest")
    print("   state is that NO single anchor makes ALL observables land at once")
    print("   (Coulomb gives alpha ~right-order; Faber gives the mass by")
    print("   construction; the clock is 28x off). The sim-units-RATIO fallback")
    print("   (predictions as ratios) is therefore the correct mode until the")
    print("   V-on / Faber-r0 core route is shown to lift omega toward the ZBW.")


def main():
    banner("M5.8.2x - ABSOLUTE-SCALE CALIBRATION  (issue #208, Phase A+B)")
    print("Scale-free model: define the scale from ONE anchor; the rest predict.")
    print("Two free dials (length, energy); c emergent, hbar = action unit (P2).")
    print()
    a_central, zbw_ratio = part_a_invariants()
    gap = part_b_unit_map()
    structural_findings(gap)


if __name__ == "__main__":
    main()
