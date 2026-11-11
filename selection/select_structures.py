"""Candidate query + diversity filter for the structure panel (RCSB read-only queries and PDB file checks).

Per family: filtered RCSB candidates -> farthest-point sampling so picks are sequence-dissimilar -> selection/selected.csv.
NOTE: the committed selection/selected.csv is the frozen panel used for all results. It came from an earlier version of this
script that also balanced picks by similarity to older PDB entries; that analysis was removed, so re-running this script
against live RCSB data will not reproduce the same list.
"""
import csv, json, os, re, sys, urllib.request, urllib.error
import numpy as np
from Bio import Align
from Bio.Align import substitution_matrices

HERE = os.path.dirname(os.path.abspath(__file__))
import json, os, re, sys, urllib.request, urllib.error
import numpy as np
from Bio import Align
from Bio.Align import substitution_matrices

HERE = os.path.dirname(os.path.abspath(__file__))
N_PER_FAMILY, MAX_ID = 8, 0.40

def post(u, b):
    return json.load(urllib.request.urlopen(urllib.request.Request(u, json.dumps(b).encode(), {"Content-Type": "application/json"}), timeout=90))
def T(a, op, v): return {"type": "terminal", "service": "text", "parameters": {"attribute": a, "operator": op, "value": v}}
def AND(*n): return {"type": "group", "logical_operator": "and", "nodes": list(n)}
def OR(*n): return {"type": "group", "logical_operator": "or", "nodes": list(n)}
DESC = "rcsb_polymer_entity.pdbx_description"

# family -> (min len, max len, search keywords, include regex, exclude regex, extra filters, max identity allowed between picks)
FAMILIES = {
 "Nanobody (VHH)": (110, 135, ["nanobody", "VHH", "single domain antibody", "sybody"],
                    r"nanobod|vhh|single.domain|sybody", r"fusion|de novo", [], 1.0),
 "Cytokine / checkpoint": (80, 260, ["interleukin", "interferon", "tumor necrosis factor", "granulocyte", "erythropoietin", "growth hormone",
                    "vascular endothelial growth factor", "transforming growth factor", "programmed cell death protein 1", "programmed cell death 1 ligand 1",
                    "cytotoxic T-lymphocyte protein 4", "lymphocyte activation gene 3", "oncostatin", "leukemia inhibitory factor", "thymic stromal"],
                    r"interleukin|interferon|tumor necrosis factor|granulocyte|erythropoietin|growth hormone|vascular endothelial growth|transforming growth|programmed cell death|cytotoxic t-lymphocyte|lymphocyte activation|oncostatin|leukemia inhibitory|stromal lymphopoietin",
                    r"receptor|regulatory|binding protein|binder|induced|inhibitor|kinase|\bfusion|hepatoma|propeptide|latency|nuclear factor|mimic|,", [], MAX_ID),
 "Viral glycoprotein": (90, 350, ["spike", "glycoprotein", "hemagglutinin", "envelope", "fusion glycoprotein", "receptor-binding domain", "attachment"],
                    r"spike|glycoprotein|hemagglutinin|envelope|fusion|receptor.binding|attachment", r"capsid|nucleocapsid|binder|nanobody|\bfab\b|antibody|chimeri|vp25|,|sortase|coatomer",
                    [T("rcsb_entity_source_organism.taxonomy_lineage.name", "exact_match", "Viruses")], MAX_ID),
 "Biocatalysis enzyme": (150, 350, ["PETase", "cutinase", "transaminase", "aminotransferase", "ketoreductase", "imine reductase", "dehalogenase",
                    "peroxygenase", "ene reductase", "old yellow enzyme", "lipase", "esterase", "nitrilase", "amine dehydrogenase", "carbonyl reductase"],
                    r"petase|cutinase|transaminase|aminotransferase|ketoreductase|carbonyl reductase|imine reductase|ene.reductase|old yellow|dehalogenase|peroxygenase|lipase|nitrilase|amine dehydrogenase",
                    r"receptor|fusion|inhibitor protein|artificial|designed|de novo", [], MAX_ID),
}
GQ = """query($ids:[String!]!){polymer_entities(entity_ids:$ids){rcsb_id
 entity_poly{pdbx_seq_one_letter_code_can rcsb_sample_sequence_length}
 rcsb_polymer_entity_container_identifiers{auth_asym_ids}
 rcsb_polymer_entity{pdbx_description}
 entry{rcsb_accession_info{initial_release_date} rcsb_entry_info{resolution_combined}}}}"""

aligner = Align.PairwiseAligner(mode="global", substitution_matrix=substitution_matrices.load("BLOSUM62"), open_gap_score=-10, extend_gap_score=-0.5)
def identity(a, b):
    al = aligner.align(a, b)[0]
    same = sum(sum(x == y for x, y in zip(a[s1:e1], b[s2:e2])) for (s1, e1), (s2, e2) in zip(*al.aligned))
    return same / min(len(a), len(b))

def candidates(fam):
    lo, hi, kws, inc, exc, extra, _ = FAMILIES[fam]
    q = {"query": AND(T("rcsb_accession_info.initial_release_date", "greater_or_equal", "2022-01-01T00:00:00Z"),
                      T("exptl.method", "exact_match", "X-RAY DIFFRACTION"), T("rcsb_entry_info.resolution_combined", "less_or_equal", 2.5),
                      T("entity_poly.rcsb_sample_sequence_length", "range", {"from": lo, "to": hi, "include_lower": True, "include_upper": True}),
                      T("entity_poly.rcsb_entity_polymer_type", "exact_match", "Protein"), OR(*[T(DESC, "contains_words", w) for w in kws]), *extra),
         "return_type": "polymer_entity", "request_options": {"paginate": {"start": 0, "rows": 1000}}}
    ids = [r["identifier"] for r in post("https://search.rcsb.org/rcsbsearch/v2/query", q)["result_set"]]
    rows = []
    for k in range(0, len(ids), 100):
        rows += [e for e in post("https://data.rcsb.org/graphql", {"query": GQ, "variables": {"ids": ids[k:k + 100]}})["data"]["polymer_entities"] if e]
    out = []
    for e in rows:
        desc = e["rcsb_polymer_entity"]["pdbx_description"]
        chains = [c for c in e["rcsb_polymer_entity_container_identifiers"]["auth_asym_ids"] if len(c) == 1]
        seq = e["entity_poly"]["pdbx_seq_one_letter_code_can"]
        if not chains or not re.search(inc, desc, re.I) or re.search(exc, desc, re.I) or set(seq) - set("ACDEFGHIKLMNPQRSTVWY"): continue
        res = (e["entry"]["rcsb_entry_info"]["resolution_combined"] or [9])[0]
        out.append({"pdb": e["rcsb_id"].split("_")[0], "chain": chains[0], "len": len(seq), "seq": seq, "desc": desc[:60],
                    "release": e["entry"]["rcsb_accession_info"]["initial_release_date"][:10], "res": res})
    return sorted(out, key=lambda r: r["res"])

def pdb_ok(pdb):
    try:
        urllib.request.urlopen(urllib.request.Request(f"https://files.rcsb.org/download/{pdb}.pdb", method="HEAD"), timeout=30); return True
    except urllib.error.URLError: return False


selected = []
for fam in FAMILIES:
    cands, maxid = candidates(fam)[:150], FAMILIES[fam][6]
    picks = []
    while cands and len(picks) < N_PER_FAMILY:
        if not picks:
            best = cands.pop(0)
        else:
            scored = sorted((max(identity(c["seq"], p["seq"]) for p in picks), c["res"], i) for i, c in enumerate(cands))
            if scored[0][0] > maxid:
                break
            best = cands.pop(scored[0][2])
        if best["pdb"] in {p["pdb"] for p in picks} or not pdb_ok(best["pdb"]):
            continue
        best["family"] = fam
        picks.append(best)
    selected += picks
    print(fam, [(p["pdb"], p["chain"], p["len"]) for p in picks], flush=True)
cols = ["family", "pdb", "chain", "len", "desc", "release", "res", "seq"]
with open(os.path.join(HERE, "selected.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, cols, extrasaction="ignore"); w.writeheader(); w.writerows(selected)
print(len(selected), "structures selected")
