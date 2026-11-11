"""Compute per-design extras and summary_results.csv (structure, family and overall levels)."""
import json
import os

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from lib import ROOT, pairwise_identity_mean

R = os.path.join(ROOT, "results")
d = pd.read_csv(os.path.join(R, "design_results.csv"))
cfg = pd.read_csv(os.path.join(R, "config_runtime.csv"))
nat = pd.read_csv(os.path.join(R, "native_scores.csv"))
meta = json.load(open(os.path.join(R, "run_metadata.json")))

# Orient scores so HIGHER = more model-preferred within each model (ProteinMPNN reports NLL, lower=better).
d["score_loglik"] = np.where(d.model == "ProteinMPNN", -d.model_score, d.model_score)
# Within-model, within-structure z-score pooled over the 3 temperatures. Never compared across models.
g = d.groupby(["model", "pdb_id"]).score_loglik
d["score_z"] = (d.score_loglik - g.transform("mean")) / g.transform("std")
d.to_csv(os.path.join(R, "design_results.csv"), index=False)

rows = []
for (pdb, model, T), x in d.groupby(["pdb_id", "model", "temperature"]):
    v = x[x.valid]
    pid = pairwise_identity_mean(list(x.sequence))
    c = cfg[(cfg.pdb_id == pdb) & (cfg.model == model) & (cfg.temperature == T)].config_seconds.iloc[0]
    rows.append({"level": "structure", "family": x.family.iloc[0], "pdb_id": pdb, "length": int(x.length.iloc[0]),
                 "model": model, "temperature": T, "n_designs": len(x),
                 "valid_pct": 100 * x.valid.mean(), "n_unique": int(x.sequence.nunique()), "unique_pct": 100 * x.sequence.nunique() / len(x),
                 "recovery_mean": v.native_identity.mean(), "recovery_median": v.native_identity.median(), "recovery_sd": v.native_identity.std(),
                 "score_loglik_mean": v.score_loglik.mean(), "score_z_mean": v.score_z.mean(),
                 "pairwise_identity_mean": pid, "pairwise_distance_mean": 1 - pid,
                 "gen_sec_per_design": x.gen_seconds.mean(), "config_seconds": c})
s = pd.DataFrame(rows)
nm = nat.set_index(["pdb_id", "model"]).native_score
s["native_score_raw"] = [nm[(p, m)] for p, m in zip(s.pdb_id, s.model)]

num = ["valid_pct", "unique_pct", "recovery_mean", "recovery_median", "recovery_sd", "score_loglik_mean", "score_z_mean",
       "pairwise_identity_mean", "pairwise_distance_mean", "gen_sec_per_design", "length"]


def pool(keys, level, **fixed):
    p = s.groupby(keys)[num].mean().reset_index()
    p["n_designs"] = s.groupby(keys).n_designs.sum().values
    p["n_unique"] = s.groupby(keys).n_unique.sum().values
    p["n_structures"] = s.groupby(keys).pdb_id.nunique().values
    p["config_seconds"] = s.groupby(keys).config_seconds.sum().values
    p["level"] = level
    for k, v in fixed.items():
        p[k] = v
    return p


fam = pool(["family", "model", "temperature"], "family", pdb_id="ALL")
over = pool(["model", "temperature"], "overall", family="ALL", pdb_id="ALL")
out = pd.concat([s, fam, over], ignore_index=True)
out.round(5).to_csv(os.path.join(R, "summary_results.csv"), index=False)

# ---- console digest (used to write the interpretation) ----
pd.set_option("display.width", 250, "display.max_columns", 30)
cols = ["model", "temperature", "recovery_mean", "pairwise_distance_mean", "score_loglik_mean", "gen_sec_per_design", "config_seconds"]
print("n structures", s.pdb_id.nunique(), "designs", len(d), "valid", int(d.valid.sum()), "unique/setting total", int(s.n_unique.sum()))
print(over[cols].round(3).to_string(), "\n")
print(fam.pivot_table(index=["family", "model"], columns="temperature", values=["recovery_mean", "pairwise_distance_mean"]).round(3).to_string(), "\n")
for m, x in s.groupby("model"):
    r, p = spearmanr(x.pairwise_distance_mean, x.recovery_mean)
    print(f"{m}: Spearman(diversity, recovery) over {len(x)} structure x T points = {r:.2f} (p={p:.3g})")
    print("  length vs sec/design Spearman", round(spearmanr(x.length, x.gen_sec_per_design)[0], 2))
print()
print("\nlikelihood-vs-recovery design-level rho:", {m: round(spearmanr(x.native_identity, x.score_z)[0], 2) for m, x in d.groupby("model")})
print("native identical designs:", int((d.native_identity == 1).sum()))
print("campaign min", round(meta["campaign_seconds"] / 60, 1), "load", meta["model_load_seconds"])
