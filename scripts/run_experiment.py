"""Full campaign. N structures x 2 models x 3 temperatures x 20 designs. Writes partial results after every setting."""
import hashlib, json, os, platform, subprocess, sys, time, urllib.request
from importlib.metadata import version

import pandas as pd
from lib import *

N = 20
OUT = os.path.join(ROOT, "results")
os.makedirs(os.path.join(OUT, "raw"), exist_ok=True)
FAM = {r["pdb"]: r for r in PANEL}


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def fetch(pdb):
    path = os.path.join(DATA, f"{pdb}.pdb")
    if not os.path.exists(path):
        urllib.request.urlretrieve(f"https://files.rcsb.org/download/{pdb}.pdb", path)


def main():
    t_campaign = time.time()
    rows, configs, native_rows, excluded = [], [], [], []
    meta = {"started": time.strftime("%Y-%m-%d %H:%M:%S"), "designs_per_config": N, "device": str(DEVICE),
            "hardware": subprocess.getoutput("sysctl -n machdep.cpu.brand_string") + ", 32 GB RAM",
            "os": platform.platform(), "python": sys.version.split()[0],
            "packages": {p: version(p) for p in ["torch", "numpy", "pandas", "scipy", "plotly", "biotite", "torch_geometric", "fair-esm", "biopython"]},
            "checkpoints": {"ProteinMPNN": {"file": "vanilla_model_weights/v_48_020.pt", "sha256": sha256(MPNN_WEIGHTS)},
                            "ESM-IF1": {"file": "esm_if1_gvp4_t16_142M_UR50.pt", "sha256": sha256(ESM_WEIGHTS)}},
            "temperatures": TEMPS,
            "mpnn_params": {"model": "v_48_020 (vanilla, 0.20 A training noise)", "backbone_noise": 0.0, "omit_AAs": "X", "batch_size": 1, "decoding": "random order"},
            "esm_params": {"sampling": "multinomial, logits/T, no beam search", "partial_seq": None},
            "model_load_seconds": {}, "proteins": {}}
    panel = []
    for r in PANEL:                                   # pre-check: both models must parse the same, X-free native sequence
        fetch(r["pdb"])
        try:
            a, b = MPNN.prepare(None, r["pdb"], r["chain"])["native_seq"], ESMIF.prepare(None, r["pdb"], r["chain"])["native_seq"]
            if a != b or "X" in a:
                raise ValueError(f"native sequences differ between models or contain X (len {len(a)} vs {len(b)})")
            panel.append(r)
        except Exception as e:
            excluded.append({"pdb": r["pdb"], "reason": repr(e)})
            print("EXCLUDED", r["pdb"], repr(e), flush=True)
    seed = 2024
    for cls in (MPNN, ESMIF):
        m = cls()
        meta["model_load_seconds"][m.name] = m.load_time
        ctxs = {}
        for r in panel:
            try:
                ctxs[r["pdb"]] = m.prepare(r["pdb"], r["chain"])
            except Exception as e:
                excluded.append({"pdb": r["pdb"], "model": m.name, "reason": repr(e)})
        m.design(next(iter(ctxs.values())), TEMPS[m.name][0])           # untimed warm-up
        for r in panel:
            pdb, chain = r["pdb"], r["chain"]
            if pdb not in ctxs:
                continue
            ctx = ctxs[pdb]; native = ctx["native_seq"]
            meta["proteins"].setdefault(pdb, {"chain": chain, "length": len(native), "native_sequence": native, "family": r["family"],
                                              "description": r["desc"], "release": r["release"], "resolution": float(r["res"]),
                                              "source_url": f"https://files.rcsb.org/download/{pdb}.pdb"})
            native_rows.append({"pdb_id": pdb, "chain": chain, "model": m.name, "native_score": m.score_native(ctx)})
            for T in TEMPS[m.name]:
                seed += 1
                set_seed(seed)
                t0, seqs = time.time(), []
                for i in range(N):
                    err = ""
                    try:
                        seq, score, dt = m.design(ctx, T)
                    except Exception as e:                       # recorded, never silently dropped
                        seq, score, dt, err = "", float("nan"), float("nan"), repr(e)
                    seen_before = seq in seqs
                    seqs.append(seq)
                    valid = is_valid(seq, len(native))
                    rows.append({"pdb_id": pdb, "chain": chain, "family": r["family"], "length": len(native),
                                 "model": m.name, "temperature": T, "seed": seed, "design_id": f"{pdb}_{m.name}_T{T}_{i + 1:02d}",
                                 "sequence": seq, "valid": valid, "native_identity": identity(seq, native) if valid else float("nan"),
                                 "model_score": score, "score_definition": m.score_name, "gen_seconds": dt,
                                 "unique": bool(seq) and not seen_before, "error": err})
                configs.append({"pdb_id": pdb, "model": m.name, "temperature": T, "seed": seed, "n": N, "config_seconds": time.time() - t0})
                with open(os.path.join(OUT, "raw", f"{pdb}_{m.name}_T{T}.fasta"), "w") as f:
                    f.write(f">native {pdb}_{chain}\n{native}\n")
                    for i, s in enumerate(seqs):
                        f.write(f">{pdb}_{m.name}_T{T}_{i + 1:02d}\n{s}\n")
                print(f"{m.name:12s} {pdb} L={len(native):3d} T={T}: {time.time() - t0:6.1f}s  rec {pd.DataFrame(rows[-N:]).native_identity.mean():.3f}  [{(time.time() - t_campaign) / 60:.1f} min]", flush=True)
            pd.DataFrame(rows).to_csv(os.path.join(OUT, "design_results.partial.csv"), index=False)   # crash-safe checkpoint
        del m
    meta["campaign_seconds"] = time.time() - t_campaign
    meta["finished"] = time.strftime("%Y-%m-%d %H:%M:%S")
    meta["excluded"] = excluded
    pd.DataFrame(rows).to_csv(os.path.join(OUT, "design_results.csv"), index=False)
    pd.DataFrame(configs).to_csv(os.path.join(OUT, "config_runtime.csv"), index=False)
    pd.DataFrame(native_rows).to_csv(os.path.join(OUT, "native_scores.csv"), index=False)
    json.dump(meta, open(os.path.join(OUT, "run_metadata.json"), "w"), indent=2)
    print(f"done: {len(rows)} designs in {meta['campaign_seconds'] / 60:.1f} min; excluded: {excluded}")


if __name__ == "__main__":
    main()
