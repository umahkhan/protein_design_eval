"""Clean timing benchmark (run on an otherwise idle machine): 5 timed designs per structure per model at the middle temperature.
Campaign wall-clock is kept separately because background load inflated it."""
import os, time
import pandas as pd
from lib import *

R = os.path.join(ROOT, "results")
d = pd.read_csv(os.path.join(R, "design_results.csv"))
N = 5
rows = []
for cls in (MPNN, ESMIF):
    m = cls()
    T = TEMPS[m.name][1]
    ctxs = {p: m.prepare(p, c) for p, c in PROTEINS.items()}
    m.design(next(iter(ctxs.values())), T)                                          # untimed warm-up
    for pdb, ctx in ctxs.items():
        set_seed(7)
        t0 = time.time()
        for _ in range(N):
            m.design(ctx, T)
        dt = (time.time() - t0) / N
        rows.append({"pdb_id": pdb, "model": m.name, "length": len(ctx["native_seq"]), "temperature": T, "n_timed": N, "sec_per_design": dt})
        print(f"{m.name:12s} {pdb} L={len(ctx['native_seq']):3d}: {dt:.3f} s/design", flush=True)
pd.DataFrame(rows).to_csv(os.path.join(R, "timing_benchmark.csv"), index=False)
