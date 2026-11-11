"""Re-time settings whose wall-clock was inflated by machine load: per-design time > 1.4x the idle-machine benchmark
(results/timing_benchmark.csv) for that structure and model, or > 3x the median of its own settings. Sequences are regenerated from the recorded seed and must match exactly; only timings change.
Originals are kept in *_original columns and the change is logged in run_metadata.json."""
import json, os, time
import pandas as pd
from lib import *

R = os.path.join(ROOT, "results")
d = pd.read_csv(os.path.join(R, "design_results.csv"))
c = pd.read_csv(os.path.join(R, "config_runtime.csv"))
meta = json.load(open(os.path.join(R, "run_metadata.json")))
if "config_seconds_original" not in c.columns:          # keep the first-pass originals if re-run
    c["config_seconds_original"] = c.config_seconds
    d["gen_seconds_original"] = d.gen_seconds
bench = pd.read_csv(os.path.join(R, "timing_benchmark.csv")).set_index(["pdb_id", "model"]).sec_per_design
med = c.groupby(["pdb_id", "model"]).config_seconds.transform("median")
ref_s = [bench[(p, m)] * 20 for p, m in zip(c.pdb_id, c.model)]
flag = c[(c.config_seconds > 1.4 * pd.Series(ref_s, index=c.index)) | (c.config_seconds > 3 * med)]
print("flagged:", flag[["pdb_id", "model", "temperature", "config_seconds"]].to_dict("records"))
retimed = []
for model, cls in (("ProteinMPNN", MPNN), ("ESM-IF1", ESMIF)):
    fm = flag[flag.model == model]
    if fm.empty:
        continue
    m = cls()
    first = d[d.model == model].pdb_id.iloc[0]
    m.design(m.prepare(first, PROTEINS[first]), TEMPS[model][0])                # untimed warm-up, as in the campaign
    for _, r in fm.iterrows():
        ctx = m.prepare(r.pdb_id, PROTEINS[r.pdb_id])
        ref = d[(d.pdb_id == r.pdb_id) & (d.model == model) & (d.temperature == r.temperature)]
        set_seed(int(r.seed)); t0 = time.time()
        out = [m.design(ctx, r.temperature) for _ in range(len(ref))]
        dt = time.time() - t0
        assert [o[0] for o in out] == list(ref.sequence), f"sequence mismatch for {r.pdb_id} {r.temperature}"
        d.loc[ref.index, "gen_seconds"] = [o[2] for o in out]
        c.loc[(c.pdb_id == r.pdb_id) & (c.model == model) & (c.temperature == r.temperature), "config_seconds"] = dt
        retimed.append({"pdb": r.pdb_id, "model": model, "temperature": float(r.temperature), "original_seconds": float(r.config_seconds), "retimed_seconds": dt})
        print(f"re-timed {r.pdb_id} {model} T={r.temperature}: {r.config_seconds:.0f}s -> {dt:.0f}s (sequences identical)", flush=True)
meta["retimed_settings"] = meta.get("retimed_settings", []) + retimed
meta.setdefault("campaign_seconds_original_wallclock", meta["campaign_seconds"])
d.to_csv(os.path.join(R, "design_results.csv"), index=False)
c.to_csv(os.path.join(R, "config_runtime.csv"), index=False)
json.dump(meta, open(os.path.join(R, "run_metadata.json"), "w"), indent=2)
