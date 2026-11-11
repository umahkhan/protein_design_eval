# Model evaluation: ProteinMPNN vs. ESM-IF1 across protein families

A framework and an evidence set for deciding how to evaluate sequence-design models on an AI protein-engineering platform. Scope: **32 backbones in four families** (nanobodies, cytokines/checkpoint proteins, viral glycoproteins, biocatalysis enzymes), 20 designs per setting, computational metrics only. Labels: **[Measured]** computed from the 3,840 designs, **[Interpretation]** a reading of those numbers, **[Hypothesis]** a product decision the data might support, **[Validate]** evidence still needed.

## 1. Evaluation framework

| Tier | Purpose | Treatment |
|---|---|---|
| **A. Hard gates** | Is the output usable at all? | Explicit pass/fail |
| **B. Performance indicators** | How does the model behave? | Measured and reported; **no universal thresholds** |
| **C. Decision-oriented metrics** | Which setting fits which campaign? | Trade-offs weighed against campaign goals |

### A. Hard gates (all passed)

| Gate | Result |
|---|---|
| Valid sequences (20-letter alphabet, correct length) | 3,840 / 3,840 |
| Successful inference | 3,840 / 3,840, 0 errors |
| No empty or malformed output | 0 found |
| Fully unique designs within each setting | 192 / 192 settings |
| Reproducible from seed | 2 of 192 settings re-run (one per model): identical; the 8 settings re-timed also regenerated identically. The rest were not re-run. |

### B. Performance indicators

Means over the 32 structures (equal weight); ± is the SD across structures of each structure's mean over its 20 designs. Within-model z-scores are not comparable across models.

| Model | T | Native recovery | Diversity (pairwise distance) | Likelihood z (within model) | s / design |
|---|---|---|---|---|---|
| ProteinMPNN | 0.1 | 0.480 ± 0.077 | 0.254 ± 0.037 | +0.99 | 0.28 |
| ProteinMPNN | 0.3 | 0.469 ± 0.076 | 0.344 ± 0.047 | +0.27 | 0.28 |
| ProteinMPNN | 0.5 | 0.447 ± 0.070 | 0.445 ± 0.052 | −1.26 | 0.28 |
| ESM-IF1 | 0.5 | 0.554 ± 0.087 | 0.320 ± 0.069 | +1.14 | 1.92 |
| ESM-IF1 | 1.0 | 0.498 ± 0.083 | 0.493 ± 0.077 | +0.08 | 1.92 |
| ESM-IF1 | 1.5 | 0.429 ± 0.082 | 0.612 ± 0.075 | −1.22 | 1.92 |

Recovery by family (mean of the family's 8 structures), lowest → highest temperature:

| Family | ProteinMPNN (T 0.1 → 0.5) | ESM-IF1 (T 0.5 → 1.5) | ESM-IF1 gain at matched diversity |
|---|---|---|---|
| Nanobody (VHH) | 0.558 → 0.515 | 0.570 → 0.425 | +2 pts (ahead in 16 of 24 comparisons) |
| Cytokine / checkpoint | 0.407 → 0.380 | 0.481 → 0.359 | +7 pts (24 of 24) |
| Viral glycoprotein | 0.438 → 0.412 | 0.513 → 0.410 | +7 pts (24 of 24) |
| Biocatalysis enzyme | 0.514 → 0.480 | 0.652 → 0.521 | +12 pts (24 of 24) |

### C. Decision-oriented metrics

- **Distinct candidates:** every setting produced 20 of 20 distinct sequences, so uniqueness does not separate settings at this sample size.
- **Quality–diversity trade-off:** see §3 and the dashboard's "Diversity vs. recovery" chart.
- **Cost per candidate:** about 0.28 s (ProteinMPNN) vs 1.92 s (ESM-IF1) per design on one CPU process, roughly 12,900 vs 1,900 designs per CPU-hour at batch size 1, growing roughly linearly with length (+0.13 s and +1.04 s per extra 100 residues). Generation for the whole campaign took about 70 minutes. GPU and batching were not measured.
- **Campaign fit:** see §5.

## 2. Temperature sensitivity findings

**[Measured]**
- Raising temperature raised diversity and lowered native recovery for both models, in every family. Pairwise distance rose 0.25 → 0.45 (ProteinMPNN, T 0.1 → 0.5) and 0.32 → 0.61 (ESM-IF1, T 0.5 → 1.5). Mean recovery fell 48.0% → 44.7% and 55.4% → 42.9%.
- The price of diversity differs. ProteinMPNN T 0.1 → 0.3 added 0.09 of pairwise distance for −1.1 points of recovery, and 0.3 → 0.5 cost 2.2 more. ESM-IF1's recovery fell by 5.6 and 6.9 points per step.
- Diversity and recovery are negatively associated across the 96 structure × temperature points for ESM-IF1 (Spearman ρ = −0.86) and ProteinMPNN (ρ = −0.54). The three temperatures of one structure are not independent, so no p-value is claimed for these.
- Both models assign lower likelihood to higher-temperature designs. At the design level, likelihood and native recovery correlate more strongly for ESM-IF1 (ρ = 0.50) than ProteinMPNN (ρ = 0.18).
- Temperature had no effect on speed.

**[Interpretation]**
- Temperature behaves as a diversity control with a model-specific price. Diversity does come at the expense of native recovery and model likelihood in this data, but not uniformly: it was nearly free for ProteinMPNN between T = 0.1 and 0.3.
- Neither high nor low temperature is inherently better. Lower temperature favours reproducing known sequence characteristics; higher temperature widens the candidate pool at the cost of similarity to the native sequence and to the model's own most-likely designs.
- Likelihood is a poor stand-in for recovery or quality, especially for ProteinMPNN.

## 3. Model tradeoffs

**[Measured]** At matched diversity (interpolated per structure along its ESM-IF1 temperature curve), ESM-IF1 designs had higher native recovery than ProteinMPNN designs in 88 of 96 structure-level comparisons, by about 7 points on average. The gap depends on the family (table above): large for enzymes, moderate for cytokines/checkpoints and viral glycoproteins, and small for nanobodies. ESM-IF1 is about 7× slower per design. Recovery varies by family for both models: cytokines/checkpoint proteins are lowest for both, while the highest family is enzymes for ESM-IF1 and nanobodies for ProteinMPNN.

**[Interpretation]** ESM-IF1 sits closer to the native sequence at a given amount of sequence variation, at higher compute cost, except on nanobodies where the two are nearly equivalent. That is a statement about similarity to a known sequence, not about design quality. No overall winner is declared, and a family-specific choice would need family-specific criteria.

## 4. Exploratory vs. exploitation-focused campaigns

Illustrative only: the criteria below are placeholders to show the mechanism, **not validated product requirements** and need scientific and stakeholder validation.

| | Exploitation-style | Exploration-style |
|---|---|---|
| Goal | Stay near known sequence characteristics | Obtain a wide range of alternatives |
| Example criteria (illustrative) | Mean recovery ≥ 0.50 | Diversity ≥ 0.40 and recovery ≥ 0.42 |
| Settings meeting them (all-structure means) | ESM-IF1 T = 0.5 only | ProteinMPNN T = 0.5, ESM-IF1 T = 1.0, ESM-IF1 T = 1.5 |
| Cares most about | Recovery, likelihood | Diversity, distinct candidates, cost per candidate |
| Cares less about | Diversity | Closeness to native |

The same placeholder threshold gives different answers by family, which is why thresholds must be set per campaign with scientists: under "recovery ≥ 0.50" no setting passes for cytokines/checkpoint proteins (best 0.481), while for enzymes ESM-IF1 at all three temperatures and ProteinMPNN at T = 0.1 and 0.3 pass. The dashboard's Campaign fit sliders reproduce this interactively and follow the family filter.

## 5. Product implications

**[Hypothesis]**
- Expose temperature as an explore/exploit control with model-specific presets, labelled by measured diversity and recovery rather than raw T.
- Low-temperature settings as defaults for exploitation campaigns; higher-temperature settings for exploration, especially where diversity is gained cheaply.
- Choose the model per family: ESM-IF1's recovery advantage was large for enzymes and small for nanobodies, so a nanobody campaign might favour the faster ProteinMPNN at similar recovery. This needs confirmation with a larger panel and CDR-aware analysis.
- Pooling candidates from both models might broaden coverage. Cross-model sequence overlap was **not measured**.
- Treat compute as secondary at this scale; revisit for very large libraries or latency-sensitive use.

## 6. Scientific limitations

- Native sequence recovery is not functional success. Model likelihood is not experimentally validated stability. Diversity does not guarantee useful novelty.
- A sequence compatible with a backbone is not necessarily foldable. **No wet-lab validation was performed**, and no structure re-prediction was used.
- 32 structures (8 per family) cannot establish general superiority of either model, and family-level comparisons rest on 8 structures each. Structure × temperature points are not independent.
- Model-specific temperatures are not equivalent.
- **Training overlap was not analysed.** Both models were trained on PDB-derived data, and some panel structures are likely re-determinations of proteins they saw, which can inflate recovery. An early ubiquitin trial gave ESM-IF1 about 95–99% recovery, consistent with memorisation. ESM-IF1 also trained on predicted structures.
- Inputs are single chains stripped of ligands, partners and numbering gaps, which can matter for enzymes, glycoproteins and antibody-like binders. Nanobody CDR loops, known to be difficult, were not analysed separately.
- Computational metrics cannot independently establish experimental design success.
- Timings are single-process CPU, batch size 1, one machine. Campaign wall-clock was inflated by background load; slowed settings were re-measured and the independent benchmark agrees.

## 7. Recommendations for prospective wet-lab validation

**[Validate]**
1. **Add independent computational checks first:** structure re-prediction (designability), stability or energy proxies, and deep-mutational-scanning comparison where available. These were intentionally excluded here.
2. **Strengthen the panel:** more structures per family, a similarity-controlled split against both models' training data (including predicted structures where possible), and CDR vs framework analysis for nanobodies.
3. **Select designs across settings:** sample from each model and temperature, not just the top-scoring designs.
4. **Include controls:** the native sequence, scrambled or random sequences, and a simple baseline.
5. **Measure the outcomes that matter:** expression, solubility, stability and, where relevant, binding or activity, chosen per family (e.g. antigen binding for nanobodies, activity for enzymes).
6. **Agree acceptance criteria with scientists before testing,** separately for exploration and exploitation and per family, then estimate the experimental hit rate per setting and use it to replace the illustrative thresholds above. Set sample sizes with a statistician.

## 8. Reproduction

See `README.md`; `./run_all.sh` regenerates every number here, and `report/report.html` is the interactive dashboard.
