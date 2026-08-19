# Descriptor validation

Refinements spec §6 asks for the new Band 2/3 descriptors to be sanity-checked
against the two demo targets **before** they are relied on in evaluation
figures. This is that check, and the numbers below are what it produced.

Regenerate it with:

```bash
python scripts/validate_descriptors.py \
    --targets '~/aptarank-data/cache/targets/*.bundle.json' \
    --candidates data/demo_candidates.csv \
    -c configs/server.yaml \
    --development-corpus data/corpus/dev_placeholder_corpus.csv
```

Run of 2026-08-14, against the prepared IGFBP3 and NDM-1 targets and the
199-sequence synthetic demo batch. The candidate set is synthetic, so treat the
candidate-side numbers as a shape check on the descriptors, not as a result.

---

## 1. Do the shape descriptors distinguish the two targets?

| target | mode | measurement | elongation | planarity |
| --- | --- | --- | --- | --- |
| 3SPU chain A (NDM-1) | pocket | cavity, d_pocket 16.7 Å, d_equiv 10.7 Å | 2.59 | 2.4 Å |
| 7WRQ chain B (IGFBP3) | surface | patch, 2519 Å² over 31 residues | 1.80 | 10.4 Å |
| 7WRQ chain B, top cavity | (cross-reference) | d_pocket 20.5 Å, d_equiv 11.8 Å | 3.32 | 2.3 Å |

Both descriptors respond, and in the expected directions: the interface patch is
rounder (elongation 1.80) than either detected cavity, and the two cavities are
elongated clefts rather than spheres (2.59 and 3.32, against 1.0 for a ball).

**A trap worth naming.** `planarity_A` is computed from whatever point cloud the
mode measures — alpha-sphere centres for a cavity, residue atoms for a patch.
Those are different populations, so the cavity's 2.4 Å and the patch's 10.4 Å are
**not comparable to each other**. A cavity's alpha spheres sit in a thin shell
inside it; a patch's atoms fill a slab of protein surface. Compare planarity
across targets *within* a mode, never across modes.

## 2. Does the radius of gyration add anything over length?

| descriptor | min | median | max |
| --- | --- | --- | --- |
| length (nt) | 20 | 35 | 70 |
| radius of gyration (Å) | 10.3 | 17.1 | 35.4 |
| footprint, Rg model (Å²) | 333 | 916 | 3932 |
| footprint, length model (Å²) | 720 | 1260 | 2520 |

On this candidate set the two are **0.966 correlated**, and candidates of equal
length differ in Rg by a median of only 0.9 Å (max 2.0 Å). So for *these*
sequences the refinement changes little — the batch is structurally homogeneous
(hairpins of varying quality), and in a hairpin the stem length tracks the
sequence length almost exactly.

That is a property of the demo batch, not of the descriptor. On structures that
genuinely differ at equal length the two models diverge sharply: a 100-nt single
helix and a 100-nt four-way junction have Rg 40.9 Å and 26.3 Å, a 2.4× difference
in footprint area, where the length model calls them identical. The real test is
the FM-generated set, which is where structural diversity is expected.

## 3. Does the choice of model change the answer?

Against the IGFBP3 patch (2519 Å²), across 199 candidates:

| model | strong | moderate | weak | "strong" candidates are… |
| --- | --- | --- | --- | --- |
| radius_of_gyration | 10 | 24 | 165 | 54–69 nt |
| length | 19 | 23 | 157 | **60–70 nt — the longest in the set** |

The two agree on 87% of bands. Where they disagree, the reason is structural:

**The length model degenerates into a length ranking here.** Its footprint is
36 Å² × N, so the best possible match to a 2519 Å² patch is N ≈ 70 — which is
exactly the longest candidate in the batch. Every "strong" call it makes is one
of the longest sequences, and its correlation between length and control
percentile is **+0.976**. Under the Rg model the optimum (Rg ≈ 28 Å) sits inside
the range rather than at its edge, the longest candidate is *not* strong, and the
correlation falls to +0.919.

That degeneracy is the substantive reason `radius_of_gyration` is the default.
It is not that the length proxy is inaccurate — both are coarse — but that its
optimum lands on the end of the candidate size range for a patch this large, so
it answers "which is longest?" while appearing to answer "which fits?".

## What this does not establish

* Neither model is calibrated against measured interface areas. They are coarse
  size proxies, and the band is a rank against shuffled controls, not an
  absolute claim about contact area.
* π·Rg² over-states a long thin molecule, which really contacts a surface along a
  strip about one helix wide. Elongated candidates therefore read larger than
  they are.
* The secondary structure fixes no angles, so multiloop arms are assumed to
  splay rather than stack coaxially. Real junctions stack, which makes real
  molecules somewhat more extended than this estimate.
