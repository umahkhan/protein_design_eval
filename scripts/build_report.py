"""Build the self-contained HTML report (report/report.html) from results/*.csv."""
import json
import os
import re

import numpy as np
import pandas as pd
import plotly
from scipy.stats import spearmanr

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.path.join(ROOT, "results")
d = pd.read_csv(os.path.join(R, "design_results.csv"))
s_all = pd.read_csv(os.path.join(R, "summary_results.csv"))
meta = json.load(open(os.path.join(R, "run_metadata.json")))
repro = json.load(open(os.path.join(R, "repro_check.json")))
S = s_all[s_all.level == "structure"]
F = s_all[s_all.level == "family"]
O = s_all[s_all.level == "overall"]
M, E = "ProteinMPNN", "ESM-IF1"
FAMS = ["Nanobody (VHH)", "Cytokine / checkpoint", "Viral glycoprotein", "Biocatalysis enzyme"]


def ov(model, T, col):
    return float(O[(O.model == model) & (O.temperature == T)][col].iloc[0])


k = {}
Tm = {M: sorted(O[O.model == M].temperature), E: sorted(O[O.model == E].temperature)}
for tag, m in (("m", M), ("e", E)):
    lo, mid, hi = Tm[m]
    k[f"{tag}_lo"], k[f"{tag}_mid"], k[f"{tag}_hi"] = lo, mid, hi
    for col, name in (("recovery_mean", "rec"), ("pairwise_distance_mean", "div"), ("score_loglik_mean", "ll")):
        k[f"{tag}_{name}_lo"], k[f"{tag}_{name}_mid"], k[f"{tag}_{name}_hi"] = ov(m, lo, col), ov(m, mid, col), ov(m, hi, col)
    x = S[S.model == m]
    k[f"{tag}_rho"], k[f"{tag}_p"] = spearmanr(x.pairwise_distance_mean, x.recovery_mean)
    dd = d[d.model == m]
    k[f"{tag}_rho_score"] = spearmanr(dd.native_identity, dd.score_z)[0]
    k[f"{tag}_drop"] = k[f"{tag}_rec_lo"] - k[f"{tag}_rec_hi"]
    k[f"{tag}_len_rho"] = spearmanr(x.length, x.gen_sec_per_design)[0]
    k[f"{tag}_sec100"] = 100 * np.polyfit(x.length, x.gen_sec_per_design, 1)[0]
k["m_mid_drec"] = k["m_rec_mid"] - k["m_rec_lo"]
k["m_mid_ddiv"] = k["m_div_mid"] - k["m_div_lo"]
k["e_mid_drec"] = k["e_rec_mid"] - k["e_rec_lo"]


def matched(sub_e):
    """ESM-IF1 recovery at ProteinMPNN's diversity, interpolated per structure along that structure's ESM-IF1 temperature curve."""
    out = []
    for pdb, e in sub_e.groupby("pdb_id"):
        e = e.sort_values("pairwise_distance_mean")
        for _, r in S[(S.model == M) & (S.pdb_id == pdb)].iterrows():
            out.append(np.interp(r.pairwise_distance_mean, e.pairwise_distance_mean, e.recovery_mean) - r.recovery_mean)
    return out


adv = matched(S[S.model == E])
k["adv_n"], k["adv_pos"], k["adv_mean"], k["adv_min"], k["adv_max"] = len(adv), int(sum(a > 0 for a in adv)), float(np.mean(adv)), min(adv), max(adv)
for i, fam in enumerate(FAMS):
    a_ = matched(S[(S.model == E) & (S.family == fam)])
    k[f"advf{i}"], k[f"advf{i}_pos"], k[f"advf{i}_n"] = float(np.mean(a_)), int(sum(x > 0 for x in a_)), len(a_)
gen = d.groupby("model").gen_seconds.mean()
k["m_sec"], k["e_sec"], k["ratio"] = gen[M], gen[E], gen[E] / gen[M]
k["m_total"] = float(S[S.model == M].config_seconds.sum())
k["e_total"] = float(S[S.model == E].config_seconds.sum())
k["campaign"] = meta["campaign_seconds"]
k["campaign_min"] = meta["campaign_seconds"] / 60
k["m_dph"], k["e_dph"] = 3600 / gen[M], 3600 / gen[E]
k["m_load"], k["e_load"] = meta["model_load_seconds"][M], meta["model_load_seconds"][E]
k["n_designs"], k["n_valid"] = len(d), int(d.valid.sum())
k["n_struct"] = d.pdb_id.nunique()
k["n_fam"] = d.family.nunique()
k["n_excl"] = len(meta.get("excluded", []))
k["n_err"] = int(d.error.fillna("").astype(bool).sum())
k["n_unique_settings"] = int((S.n_unique == S.n_designs).sum())
k["n_settings"] = len(S)
k["n_repro"] = len(repro)
k["repro_ok"] = all(repro.values())
k["campaign_clean_min"] = (float(S.config_seconds.sum()) + sum(meta["model_load_seconds"].values())) / 60
k["n_retimed"] = len(meta.get("retimed_settings", []))
k["len_min"], k["len_max"] = int(d.length.min()), int(d.length.max())
sec_T = d.groupby(["model", "temperature"]).gen_seconds.mean()
k["m_secT_spread"], k["e_secT_spread"] = float(sec_T[M].max() - sec_T[M].min()), float(sec_T[E].max() - sec_T[E].min())

cols_d = ["pdb_id", "family", "model", "temperature", "native_identity", "score_z", "gen_seconds"]
data = {
    "S": S[["family", "pdb_id", "length", "model", "temperature", "recovery_mean", "pairwise_distance_mean",
            "score_z_mean", "score_loglik_mean", "gen_sec_per_design", "valid_pct", "unique_pct"]].round(5).to_dict("records"),
    "proteins": meta["proteins"],
    "k": {a: (bool(b) if isinstance(b, (bool, np.bool_)) else float(b)) for a, b in k.items()},
}
tpl = open(os.path.join(ROOT, "scripts", "report_template.html")).read()
plotly_js = open(os.path.join(os.path.dirname(plotly.__file__), "package_data", "plotly.min.js")).read()
html = (tpl.replace("/*PLOTLY_JS*/", plotly_js.replace("</script>", "<\\/script>"))
           .replace("/*DATA*/", json.dumps(data))
           .replace("/*META*/", json.dumps({"hardware": meta["hardware"], "os": meta["os"], "python": meta["python"], "packages": meta["packages"],
                                              "started": meta["started"], "excluded": meta.get("excluded", []),
                                              "sha": {a: b["sha256"][:16] for a, b in meta["checkpoints"].items()}})))


def fmt(val, f):
    special = {"pct": f"{100 * val:.0f}%", "pp": f"{100 * val:.0f}", "c": f"{val:,.0f}", "abs": f"{abs(val):.2f}", "abspp": f"{abs(100 * val):.0f}"}
    return special[f] if f in special else f"{val:.{int(f)}f}"


for name, val in k.items():                         # {{name:fmt}} tokens in the narrative
    if isinstance(val, (bool, np.bool_)):
        continue
    for f in ("0", "1", "2", "3", "pct", "pp", "c", "abs", "abspp"):
        tok = "{{%s:%s}}" % (name, f)
        if tok in html:
            html = html.replace(tok, fmt(val, f))
os.makedirs(os.path.join(ROOT, "report"), exist_ok=True)
out = os.path.join(ROOT, "report", "report.html")
open(out, "w").write(html)
print("wrote", out, f"{len(html) / 1e6:.1f} MB", "unresolved tokens:", re.findall(r"\{\{[^}]*\}\}", html))
