# N3 checkpoint 05 , de-risk PASS + what it tells the loop search

`n3_derisk.py` PASS (summary: `n3_derisk_summary.json`). The flavour-space mass-matrix <-> PMNS
scaffold is verified, so any later TBM miss is loop PHYSICS, not bookkeeping.

| Test | Result |
| --- | --- |
| D1 bridge round-trip | `M = U_TBM diag(m) U_TBM^T` -> diagonalize -> recover TBM angles, err **2.84e-14**. M_TBM IS magic + mu-tau (both violations ~1e-16). |
| D2 Z3 democratic | `M = m0 I + t(J-I)` -> **sin^2 th12 = 1/3 exactly** (trimaximal solar), but the other doublet is **degenerate** (th23/th13 undetermined). |
| D3 magic+mu-tau scan | 400 random magic+mu-tau matrices -> **all EXACT TBM** (worst err 8.7e-13). TBM is the SYMMETRY, not a tuning. |
| D4 breaking -> theta13 | both mu-tau breaking and magic breaking turn on theta13 **linearly** (mu-tau slope 46.3 deg/eps); range up to 4.6 deg over the scan. |

## The structure the loop overlaps MUST produce (this is the S2 target, made precise)

The general magic + mu-tau flavour matrix is
```text
  M = [[x, y, y],
       [y, z, w],     with the MAGIC constraint  x + y = z + w
       [y, w, z]]
  eigenvalues: m2 = x+2y (trimaximal/solar col (1,1,1)/sqrt3),
               m1 = x-y  (col (2,-1,-1)/sqrt6),
               m3 = z-w  (reactor col (0,1,-1)/sqrt2)
  degenerate doublet (=> NOT TBM) iff x=z and y=w (the pure-democratic point).
```

So the loop field theory must give a coupling matrix with:
1. **mu-tau (2<->3) symmetry**: `M_e,mu = M_e,tau` and `M_mu,mu = M_tau,tau`. GEOMETRIC realization: the
   e-loop sits on a symmetry axis, the mu/tau loops are a MIRROR PAIR. (Pure 3-fold Z3 of three identical
   loops gives a CIRCULANT = democratic = degenerate -> NOT TBM; the e-loop must be distinguished.)
2. **magic** (`x+y=z+w`, equal row sums): ONE extra scalar condition. Codimension-1 surface in the loop
   geometry -> generically crossed by a 1D scan of any geometric knob that moves `(x+y)-(z+w)` through 0.
   This crossing IS the ★ TBM gate.
3. **away from pure-democratic** (`x != z`): lifts the doublet -> picks out th23=45, th13=0. Automatic once
   e is distinguished from mu/tau.

## theta13 channel (D4 -> S3)
Near TBM, theta13 is LINEAR in both the mu-tau-breaking and the magic-breaking. The LC `delta` (the index-2
SO(3)-breaking eigenvalue) is the physical source S3 will switch on. D4 sets the conversion: ~46 deg per unit
mu-tau asymmetry, so an effective breaking ~0.18 gives 8.5 deg. This reframes the central tension precisely:
does the LC delta=1e-10 produce an EFFECTIVE flavour-space breaking ~0.18 (needs ~1.8e9 amplification, e.g.
a near-degenerate m1<->m3 gap), or does the mixing see a larger effective delta? -> the S3 question.

## Next: S1 `n3_mass_matrix.py`
Seed 3 closed loops (mu-tau mirror arrangement, e on axis, tunable geometry) as LdG 4x4 fields (index-0,
N1 precision-safe when dressed); compute the 3x3 symmetric coupling matrix from the LdG energy overlaps;
diagonalize -> U -> angles. Verify the matrix lands in the [[x,y,y],[y,z,w],[y,w,z]] mu-tau form by
construction, then S2 scans the geometry for the magic crossing (the gate).
