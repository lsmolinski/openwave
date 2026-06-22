# Frozen sandboxes (index-3 era, pre 2026-06-21 flip)

These sandbox directories are HISTORICAL RECORDS of experiments run against the legacy index-3
engine (`D = diag(1, δ, 0, g)`, time axis = array index 3). They were **deliberately NOT edited**
by the 2026-06-21 index-0 flip , editing them would falsify the record of what actually ran
(their results were produced by index-3 code). Perfect provenance: byte-identical to when they ran.

| Frozen dir | Era |
| --- | --- |
| `sandbox_v1` .. `sandbox_v5` | early 3×3 / matrix-feasibility models (pre 4×4 engine) |
| `sandbox_v6` .. `sandbox_v8` | 4×4 promotion + M5.8 clock/GEM, index-3 engine |
| `sandbox_vn` | scratch / PMNS-SO(3) (#199), m5_8_2q/2r , index-3 engine |

## How to re-run a frozen sandbox faithfully

These scripts import the shared engine, which is now index-0. To run them with the index-3 engine
they were written against, check out the repository at the **last index-3-engine commit** (the
commit immediately BEFORE the convention flip; the flip commit message references
`research/_convention_refactor/`). As of this writing the flip is an uncommitted working-tree
change, so `HEAD` (`6a33d8c`) still carries the index-3 engine:

```text
git stash            # set aside the index-0 flip (if uncommitted)
#   or: git checkout <pre-flip-commit>
python3 research/sandbox_v8/<script>.py
git stash pop        # restore the flip
```

The physics is identical either way (the flip is proven physics-neutral, see
[`CONVENTION.md`](CONVENTION.md) + [`golden_master.py`](golden_master.py)); re-running against the
index-3 commit only reproduces the exact stored arithmetic/labels of the original run.

## NOT frozen (flipped to index-0, active)

`sandbox_v9` (#200 lepton mass, active) and `sandbox_v10` (neutrino #236) were flipped to index-0
along with the engine. See [`CONVENTION.md`](CONVENTION.md).
