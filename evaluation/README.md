# Evaluation

Experiments E1–E5 (spec §8). Everything reads the corpus or stored run
artifacts; there is no separate analysis pipeline, so the paper's numbers come
from the same code path as the demo's.

```bash
# E1 + E2 (need the corpus only)
python -m aptarank evaluate run --corpus data/corpus/validated.csv \
    --generated data/generated/rnagenesis_ndm1.csv \
    --experiments E1,E2 --folds 5 -o evaluation/results.json

# E4 + E5 (need one run artifact per target)
python -m aptarank evaluate run --corpus data/corpus/validated.csv \
    --experiments E4,E5 --artifacts runs/ndm1.json runs/kpc2.json \
    -o evaluation/results_targets.json

# figures, regenerated from stored results only
python -m aptarank evaluate figures evaluation/results.json -o evaluation/figures
```

## Comparison groups

| group | definition | role |
| --- | --- | --- |
| `validated` | real aptamers from the reference corpus | positive control |
| `random` | IID sampling at corpus letter frequencies, corpus length distribution | negative control |
| `shuffled` | dinucleotide shuffles of the validated group | **hard** negative control |
| `generated` | fine-tuned RNAGenesis output | the actual use case |

The `shuffled` group is an addition to the spec's three. Matching only
single-letter frequencies makes the negative control easy to beat on
composition alone; a dinucleotide-preserving shuffle of a real aptamer has the
same length, the same letters *and* the same adjacent-pair frequencies, so
beating it is a claim about arrangement.

Generated sequences that appear verbatim in the corpus are dropped and counted:
a memorised training sequence would otherwise be scored as if it were novel.

## E1 is out-of-fold, and why

As specified, E1 is circular: the validated aptamers define the corpus
percentiles and are then used as the positive group, so they are scored against
a distribution they themselves produced.

Here, folds are grouped **by target**, not by sequence. Two aptamers selected
against the same protein are often near-relatives, so a random split would
leave a held-out sequence's own family in the reference distribution. Each
validated sequence is scored against a corpus excluding its target's fold, and
that fold's controls are scored against the same reduced reference.

The `generated` group has no fold structure — it is scored against the whole
corpus, exactly as a user would score it.

## E3 uses control percentiles, never raw scores

A raw geometric score is not comparable across targets: each target has a
different cavity and therefore a different control distribution, so pooling raw
scores into one AUROC would be invalid. E3 compares
`tier2_control_percentile` — each score's position in *its own target's* control
distribution.

Confidence intervals bootstrap the **aptamers**, not the individual
aptamer/decoy pairs: each aptamer contributes several decoys, and treating
those as independent observations would produce an interval far too narrow.
Top-1 retrieval accuracy and MRR are reported alongside AUROC, because
"does the tool point at the right protein?" is the question a biologist asks.

Decoy targets are *presumed* non-matched, not confirmed non-binders. A null
result is publishable: it becomes an honest limitation and reinforces the
plausibility-not-prediction framing. Run it early — it affects how Tier 2 is
described throughout the paper.

## E4 caveat

Identical annotations across two targets are not automatically a software
failure: two cavities with similar dimensions *should* produce similar
geometric annotations. Compare `d_pocket_A` before reading a null as a bug, and
prefer demo targets whose pocket descriptors differ measurably.

## Reproducibility

Every experiment records the resolved config, the corpus hashes and the seed.
Figures are regenerated from a results JSON with no folding, so a reviewer can
redraw every figure from the artifacts in this repository.

`results_dev.json` and `figures/` in this directory are generated from the
**synthetic placeholder corpus**. They demonstrate that the machinery runs;
they are not findings.

## Methodological rule

Composite weights stay equal and configurable, and are never fitted to these
groups. Tuning weights to maximise validated-vs-random separation and then
reporting that separation inflates the result by construction.
