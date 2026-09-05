"""M5.32 R14 filler (owed since R12, PAUSE RECORD 3): the 12000-iteration ring ladder.
Runs m5_32_r12_a_ring.relax on the R12 ring seed (ring_seed(cfg, a)) at n32 L48, a in {6, 9},
under the R10 protocol, to 12000 iterations (cached under checkpoints/m5_32_r12/<tag>_it12000.*),
and reads the end state with R12's read_state (cord radius, E_u, kin).  Question (R12 owed item):
does the cord park near 5.4 (the audit's geometric extrapolation) or slide.
Run: python3 m5_32_r14_filler_ring12000.py <a>   (one a per process)
"""
import importlib.util, json, os, sys, time
ARGV = list(sys.argv)
HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("r12", os.path.join(HERE, "m5_32_r12_a_ring.py"))
R12 = importlib.util.module_from_spec(spec); spec.loader.exec_module(R12)
a = float(ARGV[1]) if len(ARGV) > 1 else 6.0
cfg = R12.cfg_of(32, 48.0)
tag = f"n32_L48_a{a:g}"
M0 = R12.ring_seed(cfg, a)
t0 = time.time()
M, info = R12.relax(cfg, M0, 12000, tag)
rec = {"tag": tag, "a_seed": a, "maxit": 12000, "info": info, "seed": R12.read_state(M0, cfg, tag + " seed"),
       "end": R12.read_state(M, cfg, f"{tag} it12000"), "wall_s": round(time.time() - t0, 1)}
out = os.path.join(R12.CK, f"{tag}_it12000_read.json")
json.dump(rec, open(out, "w"), indent=1, default=float)
print("DONE", tag, json.dumps(rec["end"], default=float)[:400], flush=True)
