"""Hard-gate check: re-run one setting per model from its recorded seed and compare to results/design_results.csv."""
import json, os
import pandas as pd
from lib import *

d = pd.read_csv(os.path.join(ROOT, "results", "design_results.csv"))
res = {}
for cls, T_idx in [(MPNN, 0), (ESMIF, 1)]:
    m = cls()
    sub = d[(d.model == m.name) & (d.temperature == TEMPS[m.name][T_idx])]
    pdb = sorted(sub.pdb_id.unique(), key=lambda p: sub[sub.pdb_id == p].length.iloc[0])[0]      # shortest structure keeps this fast
    T = TEMPS[m.name][T_idx]
    ctx = m.prepare(pdb, PROTEINS[pdb])
    first = d[d.model == m.name].pdb_id.iloc[0]                                                     # campaign warm-up used the first structure
    m.design(m.prepare(first, PROTEINS[first]), TEMPS[m.name][0])
    ref = d[(d.pdb_id == pdb) & (d.model == m.name) & (d.temperature == T)]
    set_seed(int(ref.seed.iloc[0]))
    again = [m.design(ctx, T)[0] for _ in range(len(ref))]
    res[f"{m.name}|{pdb}|T={T}"] = again == list(ref.sequence)
    print(m.name, pdb, T, "identical to recorded run:", res[f"{m.name}|{pdb}|T={T}"])
json.dump(res, open(os.path.join(ROOT, "results", "repro_check.json"), "w"), indent=2)
