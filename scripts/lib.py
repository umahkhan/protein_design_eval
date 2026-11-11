"""Thin wrappers around the official ProteinMPNN and ESM-IF1 code (single chain, CPU)."""
import copy
import itertools
import os
import random
import sys
import time

import numpy as np
import torch

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "vendor", "ProteinMPNN"))

DATA = os.path.join(ROOT, "data")
AA20 = set("ACDEFGHIKLMNPQRSTVWY")
MPNN_ALPHABET = "ACDEFGHIKLMNPQRSTVWYX"
MPNN_WEIGHTS = os.path.join(ROOT, "vendor", "ProteinMPNN", "vanilla_model_weights", "v_48_020.pt")
ESM_WEIGHTS = os.path.join(DATA, "esm_if1_gvp4_t16_142M_UR50.pt")
DEVICE = torch.device("cpu")

def load_panel():
    """Structure panel written by selection/select_structures.py (family, pdb, chain, len, desc ...)."""
    import csv
    path = os.path.join(ROOT, "selection", "selected.csv")
    return list(csv.DictReader(open(path))) if os.path.exists(path) else []


PANEL = load_panel()
PROTEINS = {r["pdb"]: r["chain"] for r in PANEL}
TEMPS = {"ProteinMPNN": [0.1, 0.3, 0.5], "ESM-IF1": [0.5, 1.0, 1.5]}


def set_seed(seed):
    torch.manual_seed(seed)
    random.seed(seed)
    np.random.seed(seed)


def chain_file(pdb_id, chain):
    """Single-chain, first-model, ATOM-only copy of a PDB entry, so both models receive identical input
    (ProteinMPNN's parser fails on some multi-chain files; ESM-IF1 would otherwise be given the same chain anyway)."""
    out_dir = os.path.join(DATA, "chains")
    os.makedirs(out_dir, exist_ok=True)
    out = os.path.join(out_dir, f"{pdb_id}_{chain}.pdb")
    if not os.path.exists(out):
        keep, prev_key, new_num = [], None, 0
        for line in open(os.path.join(DATA, f"{pdb_id}.pdb")):
            if line.startswith("ENDMDL"):
                break
            if line.startswith("ATOM  ") and line[21] == chain and line[16] in " A":
                key = line[22:27]
                if key != prev_key:
                    # strictly consecutive numbering: separates insertion codes (52, 52A, ...) and closes up missing-residue gaps,
                    # so both models read exactly the same residues (ProteinMPNN pads gaps with '-', ESM-IF1 skips them)
                    new_num += 1
                    prev_key = key
                keep.append(f"{line[:22]}{new_num:4d} {line[27:]}")
        open(out, "w").write("".join(keep) + "END\n")
    return out


class MPNN:
    name = "ProteinMPNN"
    score_name = "mean_nll_designed_residues (lower=better)"

    def __init__(self):
        from protein_mpnn_utils import ProteinMPNN
        t0 = time.time()
        ckpt = torch.load(MPNN_WEIGHTS, map_location=DEVICE)
        self.model = ProteinMPNN(ca_only=False, num_letters=21, node_features=128, edge_features=128,
                                 hidden_dim=128, num_encoder_layers=3, num_decoder_layers=3,
                                 augment_eps=0.0, k_neighbors=ckpt["num_edges"])
        self.model.to(DEVICE)
        self.model.load_state_dict(ckpt["model_state_dict"])
        self.model.eval()
        self.load_time = time.time() - t0
        self.omit = np.array([a in "X" for a in MPNN_ALPHABET]).astype(np.float32)
        self.bias = np.zeros(len(MPNN_ALPHABET))

    def prepare(self, pdb_id, chain):
        from protein_mpnn_utils import parse_PDB, tied_featurize
        path = chain_file(pdb_id, chain)
        d = parse_PDB(path)[0]
        chain_id_dict = {d["name"]: ([chain], [])}
        feats = tied_featurize([copy.deepcopy(d)], DEVICE, chain_id_dict, None, None, None, None, None, ca_only=False)
        return {"feats": feats, "native_seq": d[f"seq_chain_{chain}"]}

    @torch.no_grad()
    def design(self, ctx, temperature):
        """Returns (sequence, score, seconds). Score = mean NLL of designed residues."""
        from protein_mpnn_utils import _S_to_seq, _scores
        (X, S, mask, lengths, chain_M, chain_encoding_all, _cl, _vl, _ml, _mcl, chain_M_pos, omit_AA_mask,
         residue_idx, _dm, _tp, pssm_coef, pssm_bias, pssm_log_odds_all, bias_by_res_all, _tb) = ctx["feats"]
        t0 = time.time()
        randn_2 = torch.randn(chain_M.shape, device=X.device)
        pssm_log_odds_mask = (pssm_log_odds_all > 0.0).float()
        out = self.model.sample(X, randn_2, S, chain_M, chain_encoding_all, residue_idx, mask=mask,
                                temperature=temperature, omit_AAs_np=self.omit, bias_AAs_np=self.bias,
                                chain_M_pos=chain_M_pos, omit_AA_mask=omit_AA_mask, pssm_coef=pssm_coef,
                                pssm_bias=pssm_bias, pssm_multi=0.0, pssm_log_odds_flag=False,
                                pssm_log_odds_mask=pssm_log_odds_mask, pssm_bias_flag=False,
                                bias_by_res=bias_by_res_all)
        S_s = out["S"]
        log_probs = self.model(X, S_s, mask, chain_M * chain_M_pos, residue_idx, chain_encoding_all, randn_2,
                               use_input_decoding_order=True, decoding_order=out["decoding_order"])
        score = _scores(S_s, log_probs, mask * chain_M * chain_M_pos).item()
        seq = _S_to_seq(S_s[0], chain_M[0])
        return seq, score, time.time() - t0

    @torch.no_grad()
    def score_native(self, ctx):
        from protein_mpnn_utils import _scores
        (X, S, mask, lengths, chain_M, chain_encoding_all, _cl, _vl, _ml, _mcl, chain_M_pos, omit_AA_mask,
         residue_idx, *_rest) = ctx["feats"]
        set_seed(0)
        randn_1 = torch.randn(chain_M.shape, device=X.device)
        lp = self.model(X, S, mask, chain_M * chain_M_pos, residue_idx, chain_encoding_all, randn_1)
        return _scores(S, lp, mask * chain_M * chain_M_pos).item()


class ESMIF:
    name = "ESM-IF1"
    score_name = "mean_loglik_per_residue (higher=better)"

    def __init__(self):
        import esm
        t0 = time.time()
        import argparse
        # torch>=2.6 defaults to weights_only=True; the official ckpt stores an argparse.Namespace
        with torch.serialization.safe_globals([argparse.Namespace]):
            self.model, self.alphabet = esm.pretrained.load_model_and_alphabet_local(ESM_WEIGHTS)
        self.model.eval()
        self.load_time = time.time() - t0

    def prepare(self, pdb_id, chain):
        import esm.inverse_folding
        coords, seq = esm.inverse_folding.util.load_coords(chain_file(pdb_id, chain), chain)
        return {"coords": coords, "native_seq": seq}

    @torch.no_grad()
    def _score(self, ctx, seq):
        from esm.inverse_folding.util import score_sequence
        _full, withcoord = score_sequence(self.model, self.alphabet, ctx["coords"], seq)
        return float(withcoord)

    @torch.no_grad()
    def design(self, ctx, temperature):
        t0 = time.time()
        seq = self.model.sample(ctx["coords"], temperature=temperature, device=DEVICE)
        score = self._score(ctx, seq)
        return seq, score, time.time() - t0

    def score_native(self, ctx):
        return self._score(ctx, ctx["native_seq"])


# ---- metrics ----------------------------------------------------------------
def is_valid(seq, native_len):
    return bool(seq) and set(seq) <= AA20 and len(seq) == native_len


def identity(a, b):
    """Fraction of identical positions for equal-length sequences (ungapped, position-wise)."""
    if len(a) != len(b) or not a:
        return float("nan")
    return sum(x == y for x, y in zip(a, b)) / len(a)


def pairwise_identity_mean(seqs):
    seqs = [s for s in seqs if s]
    pairs = list(itertools.combinations(seqs, 2))
    vals = [identity(a, b) for a, b in pairs]
    vals = [v for v in vals if v == v]
    return float(np.mean(vals)) if vals else float("nan")
