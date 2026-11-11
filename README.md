# ProteinMPNN vs. ESM-IF1 across therapeutic and biocatalytic protein families

A reproducible pipeline that compares two inverse-folding models as sequence-design engines on **32 experimental backbones from four protein families** where protein engineering matters: nanobodies, cytokines and immune-checkpoint proteins, viral glycoproteins, and biocatalysis enzymes. Results are presented as a visual report.

**Question:** how should a protein-engineering platform evaluate competing sequence-design models when candidate quality, diversity and compute cost trade off instead of reducing to one pass/fail accuracy number?

**Live report:** **[View the interactive results in your browser](https://umahkhan.github.io/protein_design_eval/report/report.html)** (no install needed; deep links such as [`#enzyme`](https://umahkhan.github.io/protein_design_eval/report/report.html#enzyme) jump to a family).

**Deliverable:** `report/report.html`, a self-contained interactive report (also works offline when downloaded). Supporting data is in `results/`. Interpretation and the acceptance framework are in [`MODEL_EVALUATION.md`](MODEL_EVALUATION.md).

> Computational evaluation only. No wet-lab, folding or stability validation was done, and nothing here shows that a design works.

## Experiment at a glance

| | |
|---|---|
| Inputs | **32 PDB chains**, 8 per family, 106–324 residues, all released 2022 or later (panel table below) |
| Models | ProteinMPNN (`v_48_020`) and ESM-IF1 (`esm_if1_gvp4_t16_142M_UR50`); pretrained weights only, no training or fine-tuning |
| Temperatures | ProteinMPNN 0.1 / 0.3 / 0.5 · ESM-IF1 0.5 / 1.0 / 1.5 (model-specific, not equivalent) |
| Designs | 20 per structure × model × temperature = **3,840** |
| Hardware | MacBook Pro M2 Pro, 32 GB, **CPU only**, batch size 1 |
| Runtime | ≈70 min of generation (99 min measured wall-clock including load-induced slowdowns, see Methodology) |

## What goes in, and what comes out

Both models are **inverse folding** (design) models, not structure predictors. The input is a single chain's **backbone coordinates** from an experimental PDB file (ESM-IF1 uses N, CA, C; ProteinMPNN uses N, CA, C, O and builds a virtual Cβ). The native sequence is **not** given to either model; it is used afterwards only to compute recovery. No structure was predicted, and no ESMFold/AlphaFold model is used. Each PDB entry is cut down to the chosen chain (first model, ATOM records only) and residues are renumbered consecutively, so both models read exactly the same residues. Missing residues are closed up, so chain breaks are not signalled to either model. Ligands, other chains and partner proteins are not part of the input.

## Models

- **ProteinMPNN** ([dauparas/ProteinMPNN](https://github.com/dauparas/ProteinMPNN)): a message-passing network that designs residues autoregressively in random order. Weights `vanilla_model_weights/v_48_020.pt` (6.7 MB), trained with 0.20 Å backbone noise; run with no extra noise and `X` omitted.
- **ESM-IF1** ([facebookresearch/esm](https://github.com/facebookresearch/esm), archived): a GVP-Transformer encoder-decoder that samples sequences residue by residue. Checkpoint `esm_if1_gvp4_t16_142M_UR50.pt` (1.70 GB). This checkpoint has **142M parameters**, not the 124M sometimes quoted.

Temperature: both models divide logits by T before sampling. ProteinMPNN's documentation suggests 0.1–0.3 (0.5 is outside that range); ESM-IF1's README suggests a very low T (e.g. 1e-6) for maximum recovery.

## Structure panel

32 chains, 8 per family, listed in `selection/selected.csv` (frozen; used for all results). Candidates came from RCSB queries by `selection/select_structures.py`:

- X-ray, ≤2.5 Å, released 2022-01 or later, single-letter chain ID, natural proteins (de novo, fusion/chimeric and antibody-complex entries removed by keyword filters).
- Within each family, a farthest-point rule keeps picks sequence-dissimilar (≤40% identity between picks; nanobodies share a conserved framework and are exempt).
- The committed list came from an earlier version of the script that also balanced picks by similarity to older PDB entries. That analysis was removed, so the script now selects on diversity alone and re-running it against live RCSB data will not reproduce the list.
- **Training overlap was not analysed.** Both models were trained on PDB-derived data, and some panel structures are likely re-determinations of proteins they saw, which can inflate recovery. A first trial with ubiquitin gave ESM-IF1 about 95–99% recovery, probably memorisation, so classic textbook proteins were not used.

| Family | PDB | Chain | Residues | Released | Description |
|---|---|---|---|---|---|
| Nanobody (VHH) | 8P88 | C | 121 | 2024-01 | Nanobody C4 |
| Nanobody (VHH) | 9UP9 | A | 106 | 2026-01 | nanobody protein |
| Nanobody (VHH) | 8YVJ | A | 129 | 2024-12 | H5.2 nanobody (VHH) |
| Nanobody (VHH) | 9Y3I | E | 121 | 2026-09 | Nanobody VHH108 |
| Nanobody (VHH) | 10ZO | A | 109 | 2026-03 | Anti-Marburgvirus Nucleoprotein Single Domain Antibody A |
| Nanobody (VHH) | 7TGF | B | 129 | 2022-10 | VHH camelid antibody |
| Nanobody (VHH) | 7UBX | C | 121 | 2022-11 | Nanobody VHH AA6 |
| Nanobody (VHH) | 8QF4 | E | 127 | 2024-03 | Nanobody H11 |
| Cytokine / checkpoint | 7TZH | B | 173 | 2022-05 | Lymphocyte activation gene 3 protein |
| Cytokine / checkpoint | 9XNL | A | 142 | 2026-07 | Interferon 4 |
| Cytokine / checkpoint | 7TZE | A | 203 | 2022-05 | Lymphocyte activation gene 3 protein |
| Cytokine / checkpoint | 8SLR | D | 155 | 2024-04 | Tumor necrosis factor ligand superfamily member 10 |
| Cytokine / checkpoint | 8ALX | A | 115 | 2023-08 | Programmed cell death 1 ligand 1 |
| Cytokine / checkpoint | 8RYS | A | 150 | 2024-03 | Interleukin-1 beta |
| Cytokine / checkpoint | 9EHT | F | 111 | 2025-12 | Programmed cell death protein 1 |
| Cytokine / checkpoint | 7DV4 | A | 115 | 2022-01 | Cytotoxic T-lymphocyte protein 4 |
| Viral glycoprotein | 9CP0 | A | 230 | 2024-10 | Spike protein |
| Viral glycoprotein | 8V5Q | G | 116 | 2024-12 | Envelope glycoprotein E |
| Viral glycoprotein | 9H0B | A | 277 | 2025-08 | Spike glycoprotein |
| Viral glycoprotein | 9GHJ | B | 143 | 2025-06 | Glycoprotein G2 |
| Viral glycoprotein | 7QU1 | C | 152 | 2022-02 | Pre-glycoprotein polyprotein GP complex |
| Viral glycoprotein | 8YSH | B | 209 | 2024-09 | Spike glycoprotein |
| Viral glycoprotein | 7T1V | A | 315 | 2022-01 | Hemagglutinin HA1 chain |
| Viral glycoprotein | 8W0W | C | 221 | 2024-03 | Envelope glycoprotein E2 |
| Biocatalysis enzyme | 9L37 | A | 270 | 2025-10 | Cellulose-binding GDSL lipase/acylhydrolase |
| Biocatalysis enzyme | 7YPD | A | 237 | 2024-02 | Carbonyl reductase |
| Biocatalysis enzyme | 9N8R | A | 278 | 2026-02 | Lipase 1 |
| Biocatalysis enzyme | 9C9E | A | 211 | 2025-08 | Fluoroacetate dehalogenase 4A |
| Biocatalysis enzyme | 9RHW | A | 279 | 2025-08 | Monoacylglycerol lipase |
| Biocatalysis enzyme | 9NFH | A | 118 | 2025-04 | Cis-3-chloroacrylic acid dehalogenase |
| Biocatalysis enzyme | 8AHU | A | 283 | 2022-08 | Aminotransferase class IV |
| Biocatalysis enzyme | 9EVI | A | 324 | 2025-04 | Lipase B |

Panels are small (8 per family), so family-level results are descriptive.

## Setup and run

Requires macOS or Linux, Python 3.11, `curl`, and about 3 GB of free disk (1.7 GB checkpoint plus ~1.2 GB virtual environment). **Run the full pipeline on an otherwise idle machine**: about 2 hours on an M2 Pro.

```bash
./run_all.sh            # venv + pinned deps, downloads, smoke test, 3,840 designs, timing benchmark, re-timing, repro check, analysis, report
./run_all.sh report     # only re-run analysis and rebuild the report from existing results/
open report/report.html
```

Steps, in order:

1. Create `.venv` and install `requirements.txt` (pinned; ESM is installed from the official GitHub repo at a fixed commit because the PyPI 2.0.0 release imports `torch_scatter`, which is not needed and is hard to build on Apple Silicon).
2. Download the model files and verify checkpoint SHA-256 hashes.
3. `scripts/smoke_test.py`: loads both models and designs a few sequences.
4. `scripts/run_experiment.py`: the 3,840-design campaign (downloads the panel's PDB files).
5. `scripts/benchmark_timing.py`: idle-machine timing reference (5 designs per structure per model).
6. `scripts/retime.py`: re-measures settings slowed by machine load; sequences are regenerated from the recorded seed and must match exactly.
7. `scripts/repro_check.py`: re-runs one setting per model from its seed.
8. `scripts/analyze.py`, `scripts/build_report.py`: metrics, summaries and `report/report.html`.

## Methodology

- **Generation:** one fixed seed per structure × model × temperature (seeds 2025–2216); 20 designs sampled one at a time; an untimed warm-up design precedes each model's timed run.
- **Validity:** 20-letter amino-acid alphabet, length equals the native chain, non-empty, no inference error. Failures would be recorded, not dropped (there were none).
- **Native sequence recovery:** position-by-position identity to the experimental sequence. Reported per structure (mean of its 20 designs); **the structure is the unit of replication**, so overall numbers are means over 32 structures (equal weight).
- **Model score:** each model's own per-residue likelihood of the designed sequence given the backbone. ProteinMPNN reports mean negative log-likelihood (lower is better) and ESM-IF1 reports mean log-likelihood (higher is better). `score_loglik` flips ProteinMPNN's sign; `score_z` z-scores it within model and structure. **Raw scores are never compared across models.**
- **Diversity:** mean pairwise distance, 1 − (mean position-wise identity over all 190 pairs of the 20 designs). Unique-sequence counts are also recorded (100% in every setting, so they do not discriminate).
- **Efficiency:** model load is timed separately. Per-design time includes scoring that design with the model's own likelihood. Campaign wall-clock was 99 min, but background load slowed some settings (up to ~2×, one 17×). Settings more than 1.4× slower than the idle-machine benchmark (8 of 192) were re-measured with the same seeds (sequences verified identical); the reported efficiency numbers use those timings, and the original timings are kept in `*_original` columns. `results/timing_benchmark.csv` is the independent idle-machine reference.

## Outputs

| File | Contents |
|---|---|
| `report/report.html` | Interactive dashboard for a product audience: family filter, KPIs, five charts (match vs. variety, by family at default settings, match across proteins, model confidence, time per design), summary table with quality checks, and a goal-fit tool. Deep link to a family with `report.html#enzyme` |
| `results/design_results.csv` | One row per design (3,840): `pdb_id, chain, family, length, model, temperature, seed, design_id, sequence, valid, native_identity, model_score, score_definition, gen_seconds, unique, error, max_id_pre2021, score_loglik, score_z, gen_seconds_original`. `seed` is the setting's seed; `unique` marks the first occurrence of a sequence within its setting. |
| `results/summary_results.csv` | Summaries at three levels (`level` = structure, family, overall) |
| `results/raw/*.fasta` | Native sequence and the 20 designs for each of the 192 settings |
| `results/run_metadata.json` | Hardware, OS, package versions, checkpoint hashes, parameters, load times, native sequences, re-timed settings |
| `results/timing_benchmark.csv`, `config_runtime.csv`, `native_scores.csv`, `repro_check.json` | Timing reference, per-setting runtimes, native-sequence scores, reproducibility check |
| `selection/` | Panel selection code and the frozen `selected.csv` |
| `archive/pilot_3protein_results/` | Superseded 3-protein pilot (8YK7, 8GQQ, 8GC6) kept for transparency |

## Reproducibility

Seeds, package versions, checkpoint hashes (`be4ba36e…` ESM-IF1, `c9cb4a67…` ProteinMPNN) and hardware are in `results/run_metadata.json`. Re-running one setting per model from its recorded seed reproduced the sequences exactly (2 of 192 settings checked, plus the 8 re-timed settings, which also matched). CPU timings vary between runs and machines; sequences should reproduce on the same package versions and platform. The panel selection depends on live RCSB data, so re-running `select_structures.py` later may pick different structures; use the committed `selected.csv`.

## Licences and redistribution

This repository does not include model weights or PDB files; `run_all.sh` downloads them. Check the licence terms of the ProteinMPNN and ESM-IF1 weights before redistributing them. PDB entries are subject to wwPDB/RCSB terms.

## Limitations

Native recovery is not functional success; model likelihood is not measured stability; diversity does not guarantee useful novelty; a backbone-compatible sequence may not fold; no wet-lab validation was done; 32 structures (8 per family) cannot establish general superiority of either model; temperatures are not equivalent between models; training overlap with the models' training data was not analysed; structure-by-temperature points are not independent; and computational metrics cannot on their own establish experimental design success. Nanobody CDR loops are known to be hard for inverse-folding models and were not analysed separately. See `MODEL_EVALUATION.md`.

## Layout

```
run_all.sh  requirements.txt  README.md  MODEL_EVALUATION.md
scripts/    lib.py, smoke_test.py, run_experiment.py, benchmark_timing.py, retime.py, repro_check.py, analyze.py, build_report.py, report_template.html
selection/  select_structures.py, selected.csv
results/    CSVs, raw FASTA, metadata        report/  report.html        archive/  superseded pilot
data/ vendor/   (downloaded; not committed)
```
