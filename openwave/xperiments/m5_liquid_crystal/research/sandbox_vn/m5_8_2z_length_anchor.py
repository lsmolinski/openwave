"""
M5.8.2z - The Faber r0 LENGTH-anchor calibration (issue #218)

#208 anchored the unit map on ENERGY (H_static <-> 0.511 MeV) + action<->hbar
and got the clock 28x BELOW the electron ZBW. #217 showed neither V-on nor the
faithful kinetic closes it, and pointed at a missing LENGTH anchor. This script
does the length route: anchor the map on the Faber r0 (the validated length
2.2132 fm) + the radial cone c, and recompute the absolute clock omega.

Two INDEPENDENT determinations of omega_phys:
  ENERGY route (#208):  omega_phys = omega1 * (m_e c^2 / hbar) / H_static
  LENGTH route (#218):  omega_phys = omega1 * c_phys * r0_sim / (c_sim * r0_phys)
If the model is scale-consistent they agree; the 28x gap IS their disagreement.

MEASURED INPUTS (sim units, sourced):
  omega1   = 1.188     clock fundamental, rad/sim-time   (m5_8_2j / m5_8_2h)
  H_static = 16.74     rest energy, hbar=1               (m5_8_2j)
  r0_sim   = 0.69      Faber regularization radius (s=0.5 half-melt x sqrt(3)),
                       measured from the seed _m5_8_2cb_ref.npz (h=0.5217)
  c_sim    = 2.0       radial tilt cone c(r^)=2 (m5_8_2k; gradient per-sim-len
                       via inv2h). Bracketed over {1,2} for the cone convention.
  r0_phys  = 2.2132 fm Faber radius -> 0.511 MeV (m5_6_3a/3b, E*r0=const)
"""

import numpy as np

# physical constants
MC2_J = 0.511e6 * 1.602176634e-19          # m_e c^2 (J)
HBAR = 1.054571817e-34                      # J s
C_PHYS = 2.99792458e8                       # m/s
FM = 1e-15                                  # m
ZBW = 2.0 * MC2_J / HBAR                     # electron Zitterbewegung 2 m_e c^2/hbar
HBARC_MEVFM = 197.3269804                    # hbar c (MeV fm)
ALPHA = 1.0 / 137.035999

# measured sim inputs
OMEGA1 = 1.188
H_STATIC = 16.74
R0_SIM = 0.69
R0_PHYS = 2.2132 * FM


def banner(t):
    print("=" * 74); print(t); print("=" * 74)


def main():
    banner("M5.8.2z - the Faber r0 LENGTH-anchor calibration (issue #218)")
    print(f"electron ZBW target: 2 m_e c^2/hbar = {ZBW:.3e} rad/s\n")

    # ENERGY route (#208)
    om_E = OMEGA1 * (MC2_J / HBAR) / H_STATIC
    gap_E = ZBW / om_E
    print("ENERGY route (#208):  omega_phys = omega1*(m_ec^2/hbar)/H_static")
    print(f"   = {om_E:.3e} rad/s   ->  {gap_E:.1f}x BELOW the ZBW\n")

    # LENGTH route (#218), bracketed over the cone convention
    print("LENGTH route (#218):  omega_phys = omega1*c_phys*r0_sim/(c_sim*r0_phys)")
    print(f"   {'c_sim':>6} {'omega_phys':>14} {'vs ZBW':>12}")
    om_L = {}
    for c_sim in (1.0, 2.0):
        o = OMEGA1 * C_PHYS * R0_SIM / (c_sim * R0_PHYS)
        om_L[c_sim] = o
        print(f"   {c_sim:>6.1f} {o:>14.3e} {('x%.1f HIGH' % (o/ZBW)):>12}")
    print()

    # The two anchors bracket the ZBW. Their geometric mean:
    banner("THE KEY RESULT - the two anchors BRACKET the ZBW")
    for c_sim in (2.0, 1.0):
        gm = np.sqrt(om_E * om_L[c_sim])
        split = np.sqrt(om_L[c_sim] / om_E)   # the factor each anchor is off (opposite signs)
        print(f"  c_sim={c_sim:.0f}:  energy {om_E:.2e} (/{gap_E:.0f})  x  length {om_L[c_sim]:.2e} (x{om_L[c_sim]/ZBW:.0f})")
        print(f"          geometric mean = {gm:.3e} rad/s  vs ZBW {ZBW:.3e}  ->  x{gm/ZBW:.2f}")
        print(f"          (each anchor off by ~x{split:.0f} in OPPOSITE directions)\n")

    banner("DIAGNOSTIC - the model's E*r product vs the Compton relation")
    # model E*r0 in physical: 0.511 MeV * 2.2132 fm
    Er_model = 0.511 * 2.2132                  # MeV fm
    print(f"  model:   H_static<->0.511 MeV, r0<->2.2132 fm  =>  E*r0 = {Er_model:.3f} MeV fm")
    print(f"  Faber:   E*r0 = alpha*(pi/4)*hbar c = {ALPHA*np.pi/4*HBARC_MEVFM:.3f} MeV fm  (the CLASSICAL-radius relation)")
    print(f"  Compton: E*lambda_C = hbar c = {HBARC_MEVFM:.1f} MeV fm  (what a ZBW clock needs)")
    print(f"  => the model's particle obeys the classical-radius relation (E*r0 = alpha*pi/4 * hbar c),")
    print(f"     NOT the Compton relation; r0 = alpha*(pi/4)*lambda_C is ~{1/(ALPHA*np.pi/4):.0f}x below lambda_C.")
    print()

    banner("VERDICT")
    gm2 = np.sqrt(om_E * om_L[2.0])
    print("The LENGTH anchor does NOT independently close the gap (it OVERSHOOTS")
    print(f"~36x at the measured cone c=2, mirror-image of the energy route's 28x).")
    print("BUT the two anchors BRACKET the ZBW nearly symmetrically, and their")
    print(f"GEOMETRIC MEAN reproduces it to ~{abs(gm2/ZBW-1)*100:.0f}% at c=2 (within ~2x over c in [1,2]).")
    print()
    print("=> the absolute-scale error is a SINGLE factor ~32 carried in OPPOSITE")
    print("   directions by the energy and length handles, NOT an irreducible")
    print("   structural deficit. The ZBW scale IS recoverable when energy and")
    print("   length are calibrated JOINTLY (the geometric mean), rather than")
    print("   from one anchor. Recommendation for the unit map: anchor on BOTH")
    print("   (the Faber E*r0=const line) and read the scale off the joint fit,")
    print("   not the energy postulate alone. The residual factor ~2 and the cone")
    print("   convention are the next refinements.")


if __name__ == "__main__":
    main()
