"""
M5.8.2y - The 28x absolute-omega gap: lever budget (issue #217)

#208 quantified the gap (clock runs 28.2x below the electron ZBW under the
energy + action postulates) and named two candidate physics levers: the V-on
LdG core and the faithful F-commutator kinetic. This script settles, from
VALIDATED diagnoses, how much of the 28x each lever can close.

The two levers' signs are already established by the M5.6.5 diagnoses:

  LEVER 1 - V(M) on (the LdG amplitude potential, off-by-default)
    M5.6.5c finding: V(M) is ROTATION-INVARIANT. It pins the amplitude
    (Tr M^2) but NOT the frame ORIENTATION; the directors still slosh with V
    on. The de Broglie clock IS that orientation/twist slosh. So V-on does
    NOT touch the clock frequency => NULL lever on the gap (factor 1.00).
    [refutes the #208 'core-geometry lifts omega' guess for the V lever]
    Source: 0b_M5_roadmap.md M5.6.5c; sandbox_v6/m5_6_5c_potential_confinement.py

  LEVER 2 - the faithful F-commutator kinetic (vs the shipped 1/2||Mdot||^2)
    M5.6.5d finding (re-run 2026-06-17): the faithful metric G's 5 PHYSICAL
    inertia eigenvalues span [0.054, 1.447] vs the simple kinetic's uniform
    0.5, so the clock frequency is mis-set by x[sqrt(0.5/1.447), sqrt(0.5/
    0.054)] = x[0.59, 3.05]. Best case (clock = softest mode) RAISES omega by
    x3.05. So the faithful kinetic closes AT MOST x3.05 of the 28x.
    Source: sandbox_v6/m5_6_5d_faithful_kinetic.py (Stage 2c)
"""

import numpy as np

GAP = 28.2                      # #208: ZBW / omega_phys
INERTIA_SPAN = (0.054, 1.447)   # 5d: faithful physical-mode inertia (5-95%)
INERTIA_MEDIAN_PHYS = 0.265     # 5d median sorted physical eigenvalue (mid mode)
SIMPLE_INERTIA = 0.5            # the shipped 1/2||Mdot||^2 uniform inertia


def freq_factor(inertia):
    """omega ~ 1/sqrt(inertia); faithful vs simple => sqrt(simple/faithful)."""
    return np.sqrt(SIMPLE_INERTIA / inertia)


def main():
    print("=" * 74)
    print("M5.8.2y - the 28x absolute-omega gap: LEVER BUDGET  (issue #217)")
    print("=" * 74)
    print(f"Gap to close: {GAP:.1f}x  (clock runs 28x BELOW the electron ZBW, #208)\n")

    # Lever 1: V-on
    v_on_factor = 1.00
    print("LEVER 1  V(M) on  ->  factor %.2f  (NULL)" % v_on_factor)
    print("   V is rotation-invariant (M5.6.5c): confines amplitude, not the")
    print("   orientation/twist sector the clock lives in. Does not move omega.\n")

    # Lever 2: faithful kinetic
    f_best = freq_factor(INERTIA_SPAN[0])    # softest mode -> largest omega rise
    f_worst = freq_factor(INERTIA_SPAN[1])   # stiffest mode -> omega drop
    f_mid = freq_factor(INERTIA_MEDIAN_PHYS)
    print("LEVER 2  faithful F-commutator kinetic  ->  factor x[%.2f, %.2f]" % (f_worst, f_best))
    print("   physical-mode inertia span [%.3f, %.3f] vs simple %.1f (M5.6.5d)." %
          (INERTIA_SPAN[0], INERTIA_SPAN[1], SIMPLE_INERTIA))
    print("   best case (clock = softest mode): omega rises x%.2f" % f_best)
    print("   mid case  (clock = median mode):  omega x%.2f" % f_mid)
    print("   sign is clock-mode dependent; the BOUND is x%.2f (cannot exceed).\n" % f_best)

    # Residual budget
    print("RESIDUAL after the levers (best case = V-on null x faithful x%.2f):" % f_best)
    res_best = GAP / (v_on_factor * f_best)
    res_mid = GAP / (v_on_factor * f_mid)
    print("   best case:  %.1fx / %.2f = %.1fx STILL UNACCOUNTED" % (GAP, f_best, res_best))
    print("   mid  case:  %.1fx / %.2f = %.1fx STILL UNACCOUNTED" % (GAP, f_mid, res_mid))
    print()

    print("=" * 74)
    print("VERDICT")
    print("=" * 74)
    print("The 28x gap is LARGELY STRUCTURAL. The two named physics levers cannot")
    print("close it:")
    print("  - V-on is a NULL lever on frequency (rotation-invariant).")
    print("  - the faithful kinetic is the ONE real lever, but bounded to x%.2f," % f_best)
    print("    so it closes at most ~3x, leaving >= %.0fx unaccounted." % res_best)
    print("The dominant residual (~%.0f-%.0fx) is NOT a dynamics error: it is the" % (res_best, res_mid))
    print("missing LENGTH ANCHOR at V=0 (no independent length scale fixes the")
    print("absolute frequency; the action<->hbar postulate carries it). Closing it")
    print("needs the Faber r0 length chain (E*r0=const, r0=2.2132 fm), i.e. an")
    print("INDEPENDENT length calibration -- not the V-on or kinetic physics.")
    print()
    print("So #208's 'structural' read is confirmed and quantified: lever budget")
    print("~3x physical (faithful kinetic), ~9-20x structural (length anchor). The")
    print("sim-units-ratio fallback stands; the absolute Hz needs the r0 length fix.")


if __name__ == "__main__":
    main()
