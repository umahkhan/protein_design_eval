"""Milestone 2: load each model, design a few sequences for one protein at one temperature."""
import sys
import numpy as np
from lib import *

import os, urllib.request
pdb, chain, n = PANEL[0]["pdb"], PANEL[0]["chain"], 5
if not os.path.exists(os.path.join(DATA, f"{pdb}.pdb")):
    urllib.request.urlretrieve(f"https://files.rcsb.org/download/{pdb}.pdb", os.path.join(DATA, f"{pdb}.pdb"))
for cls, T in [(MPNN, 0.3), (ESMIF, 1.0)]:
    m = cls()
    print(f"\n=== {m.name}: load {m.load_time:.2f}s  score={m.score_name}")
    ctx = m.prepare(pdb, chain)
    native = ctx["native_seq"]
    print(f"native ({len(native)}): {native}")
    print(f"native score: {m.score_native(ctx):.4f}")
    set_seed(42)
    seqs = []
    for i in range(n):
        seq, score, dt = m.design(ctx, T)
        seqs.append(seq)
        print(f"T={T} #{i} valid={is_valid(seq, len(native))} len={len(seq)} rec={identity(seq, native):.3f} "
              f"score={score:.4f} t={dt:.2f}s\n   {seq}")
    print(f"unique {len(set(seqs))}/{n}  mean pairwise identity {pairwise_identity_mean(seqs):.3f}")
    # reproducibility check
    set_seed(42)
    again = m.design(ctx, T)[0]
    print("reproducible (first seq matches on reseed):", again == seqs[0])
