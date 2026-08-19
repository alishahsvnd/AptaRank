# AptaRank

Two-tier evaluation and ranking of foundation-model-generated RNA aptamer candidates.

You supply a file of candidate RNA sequences, a target protein identifier and
**how you believe an aptamer would engage that target**. AptaRank scores every
candidate, returns a ranked shortlist, and shows the evidence behind each
ranking.

**AptaRank does not predict binding.** Tier 1 measures intrinsic structural
quality against a corpus of experimentally validated aptamers. Tier 2 measures
whether a candidate's geometry is *plausible* against the target's binding site,
using the comparison appropriate to the binding mode you asserted. Tier 2
annotates the ranking and never reorders it.

> Tier 2 provides a control-calibrated, target-specific geometric
> size-agreement annotation — not a binding prediction — and never alters the
> target-independent Tier 1 ranking.

A "strong" band means **strong control-relative geometric agreement**: this
candidate matches this binding site's dimensions better than ~95% of
dinucleotide-shuffled controls do. It does not mean "strong candidate", and it
certainly does not mean "binds". The band colours are steps of one hue rather
than a green/amber/red status palette for exactly that reason.

### Two vocabularies, deliberately

The interface and the paper name the same two things differently, and each
surface uses one of them consistently.

| Concept | User-facing (the dashboard) | Paper / spec / JSON artifact |
| --- | --- | --- |
| Tier 1 | aptamer-likeness | intrinsic structural quality |
| Tier 2 | aptamer-target compatibility | control-relative geometric agreement |

Artifact field names stay in the paper's vocabulary (`tier1_*`, `tier2_*`,
`geometric_agreement_*`): the artifact is the reproducibility record behind the
paper, not a screen.

### Binding modes

The old Tier 2 assumed one mechanism — an unpaired loop plugging a cavity — and
presented it as universal. It is not. Aptamers also lie a helix along a groove,
or cover a flat surface patch through shape and charge complementarity. Baking
in the pocket model overstated what the comparison knew.

**The tool does not infer the mechanism; the expert asserts it, and the
comparison adapts.** That division of labour is the point: the biologist
supplies the biological judgement, the tool supplies fast, reproducible, coarse
geometry.

| Mode | Premise | Compares | Demo target |
| --- | --- | --- | --- |
| `pocket` | a flexible loop engages a concave cavity | loop reach (Å) vs cavity width | NDM-1 (3SPU) |
| `surface` | the molecule covers a flat-ish surface patch | footprint (Å²) vs patch area, plus patch charge | IGFBP3 (7WRQ chain B) |
| `groove` | a helical stem lies along a channel | **future work — described, not built** | — |

Implements `IMPLEMENTATION_SPEC.md` as amended by `REFINEMENTS_SPEC.md`; see
[Deviations from the spec](#deviations-from-the-spec) for the places where the
implementation deliberately differs from both.

---

![dashboard](docs/dashboard_top.png)

## Status

| Component | State |
| --- | --- |
| Ingest & validation (§4.3) | done |
| Tier 1 features: composition, folding, elements, ensemble sampling (§4.4–4.7) | done |
| Shuffled controls (§4.8) | done |
| Corpus precompute + cache, percentile scoring (§3.5.1, §4.9) | done |
| Composite scoring and ranking (§4.9, §6) | done |
| Identifier-based target input, prepared server-side (R§3) | done |
| Shared measurement layer: cavity + patch geometry, freeSASA (R§4) | done |
| Mode A pocket engagement (R§5.1) | done |
| Mode C surface-patch recognition, incl. APBS charge term (R§5.2) | done |
| Mode B groove/helical presentation (R§5.3) | future work, refused with that answer |
| Radius-of-gyration footprint from the element graph (R§6) | done, now the default |
| Descriptor validation against both demo targets (R§6) | done — [docs/descriptor_validation.md](docs/descriptor_validation.md) |
| Calibration bank and banding, per mode (§5.8) | done |
| Spearman tier-independence diagnostic (§6.4) | done, every run |
| Run artifact (§7.2), mode-aware explanations (§7.3), diagrams (§7.4) | done |
| Streamlit dashboard, panels (a)–(f) (§7.1) | done |
| Evaluation E1, E2, E4, E5 (§8) | done |
| Evaluation E3 (matched vs decoy target) | machinery done; needs labelled pairs + ≥2 prepared targets |

**Nothing external is blocking any more.** Both halves of a publishable run exist:

* **Reference library.** Each user uploads their own in step 2 — there is no one
  institutional corpus. The lab server carries a curated example: 336
  experimentally validated RNA aptamers across 146 targets, from the UTexas
  Aptamer Database, with a `reference_train.manifest.json` recording source,
  curator and date. A library that parses but carries no manifest still works
  and is marked `unverified_corpus_provenance` — accepted, not publishable.
* **Targets.** Prepared server-side from a PDB or UniProt identifier (R§3) on
  any machine with fpocket. `scripts/make_dev_target.py` still writes a
  *synthetic* pocket bundle for machines without it, and any run using one is
  ineligible for publication.

A run is publishable only when its reference distributions and its target
evidence are both real; `development_reasons` in the artifact says which are not,
and the dashboard shows that same verdict **before** the run starts.

---

## Install

Python 3.10+. On Windows, use a CPython install (not the MSYS2 one).

```bash
python -m venv .venv
.venv/Scripts/python -m pip install -e ".[dev]"
```

### ushuffle on Windows

`pip install ushuffle` fails on Python ≥ 3.9: the sdist ships pre-generated
Cython C that references `tp_print`, removed from CPython in 3.9. Build it from
the `.pyx` instead — `setup.py` re-cythonizes when Cython is importable, which
build isolation prevents:

```bash
.venv/Scripts/python -m pip install cython setuptools wheel
.venv/Scripts/python -m pip download ushuffle --no-binary :all: --no-deps -d /tmp/ush
tar -xzf /tmp/ush/ushuffle-1.1.2.tar.gz -C /tmp/ush
.venv/Scripts/python -m pip install --no-build-isolation /tmp/ush/ushuffle-1.1.2
```

Requires the MSVC build tools. Verified working with Cython 3.2.9 / MSVC 14.38.

### fpocket and APBS (Tier 2 only)

Both are Linux CLIs with no Windows build, needed only to *prepare* a target —
never to score against one. `deploy/server_install.sh` installs both under
`~/.local` without root. freeSASA, which surface mode needs, is a normal Python
dependency and installs everywhere.

Without APBS, surface mode still scores on size and records the charge term as
`not_computed`. Without fpocket, pocket mode cannot prepare a target at all.

---

## Quickstart

```bash
# 1. synthetic development data (NOT scientific data)
python scripts/make_dev_data.py
python scripts/make_dev_target.py            # synthetic pocket target

# 2. precompute the reference feature matrix once
python -m aptarank corpus build --development-corpus data/corpus/dev_placeholder_corpus.csv

# 3. score a candidate batch, Tier 1 only
python -m aptarank run data/demo_candidates.csv \
    --development-corpus data/corpus/dev_placeholder_corpus.csv

# 4. ...or both tiers, against a target described in a small file
python -m aptarank run data/demo_candidates.csv \
    --development-corpus data/corpus/dev_placeholder_corpus.csv \
    --target-file igfbp3.txt

# 5. inspect a run, then open the dashboard
python -m aptarank artifact summary runs/run_*.json
streamlit run dashboard/streamlit_app.py
```

A target is described by an identifier, a chain and a binding mode (R§3.2). The
structure is fetched and measured **server-side**; nobody builds or uploads a
JSON by hand:

```yaml
# igfbp3.txt
target_name: IGFBP3
source: pdb                 # pdb | alphafold
id: 7WRQ                    # PDB ID, or UniProt accession for alphafold
chain: B                    # which chain is the target protein
binding_mode: surface       # pocket | surface
partner_chain: C            # used to confirm the site, then stripped (R§3.4)
strip_hetatm: true
target_site_residues: [7, 8, 9, 12, 38, 55, 57, 75, 187, 210, 225, 227]
```

Residue numbers are plain integers in the structure's **author numbering** for
the chain named above — write `42`, not `K42`. They are the depositor's labels,
not positions in the sequence, and nothing is re-indexed.

The same thing from flags, or as a standalone step:

```bash
python -m aptarank run candidates.csv --corpus data/corpus/validated.csv \
    --target 3SPU --target-chain A --binding-mode pocket \
    --set 'tier2.target.target_site_residues=[120,122,189]' \
    --set 'tier2.target.retain_hetero_resnames=["ZN"]'

python -m aptarank target build --target-file igfbp3.txt
python scripts/verify_target_bundle.py cache/targets/7WRQ_*.bundle.json --rebuild
```

Every tunable lives in `configs/default.yaml` and can be overridden per run:

```bash
python -m aptarank run candidates.csv --corpus data/corpus/validated.csv \
    --set tier1.shuffle.n_shuffles=99 --set output.n_diagrams=100
```

`--fast` skips shuffled controls and ensemble sampling for exploratory runs.

---

## Repository layout

```
configs/default.yaml            every tunable; copied verbatim into each artifact
src/aptarank/
  config.py                     layered config + the criteria registry
  ingest.py                     §4.3 read, normalise, validate, deduplicate
  provenance.py                 tool versions, content hashes, seed derivation
  pipeline.py                   headless orchestration
  cli.py                        run / corpus build / artifact summary
  tier1/
    folding.py                  the only ViennaRNA entry point (§4.5, §4.7, §7.4)
    elements.py                 forgi structural elements, molecular size (§4.6, R§6)
    features.py                 per-sequence features + parallel batch driver
    shuffles.py                 dinucleotide-preserving controls (§4.8)
    corpus.py                   reference corpus load, precompute, cache (§3.5.1)
    scoring.py                  corpus percentiles + composite (§4.4, §4.9)
    service.py                  Tier 1 orchestration and the ranking
  tier2/
    target.py                   fetch (PDB/AlphaFold), chain + partner order,
                                hetero policy (§5.3, R§3)
    fpocket.py                  strict runner + parser for fpocket output (§5.4)
    geometry.py                 cavity extent, aptamer reach, footprint (§5.4–5.6)
    patch.py                    binding-site patch area and shape, freeSASA (R§4.1)
    selection.py                which cavity is the functional one (§5.4)
    modes.py                    what each binding mode compares (R§5)
    bundle.py                   immutable checksummed target evidence
    build.py                    the one module that needs fpocket (Linux)
    calibration.py              fixed control bank, percentiles, bands (§5.8)
    apbs.py                     electrostatics: cavity or patch potential (§5.7)
    service.py                  annotate the survivors; never reorder them
  artifacts/
    io.py                       run artifact assembly and I/O (§7.2)
    explanation.py              rule-based explanations (§7.3)
    rendering.py                embedded structure diagrams (§7.4)
  evaluation/
    groups.py                   comparison groups + target-grouped folds (§8.1)
    experiments.py              E1–E5 (§8.2)
    stats.py                    effect sizes, bootstrap CIs, retrieval metrics
    figures.py                  paper figures, drawn from stored results only
dashboard/                      Streamlit app; reads artifacts, computes nothing
scripts/                        dev data, dev target, bundle verifier, descriptor check
.github/workflows/              tests (Linux + Windows) and target-bundle builds
tests/                          255 tests, including the §4.5–4.6 spec fixtures
```

Boundaries that matter: workers fold, the parent scores (the corpus is never
shipped to a subprocess); the dashboard and evaluation scripts will read
artifacts only and compute nothing of their own.

---

## Deviations from the spec

Each of these is a deliberate change, not an oversight.

### 1. `tier1_score` is absolute, not a within-batch rank aggregate (§4.9)

The spec's §4.9 snippet computes the composite from `rankdata` *within the
submitted batch*. That makes a candidate's score depend on what else was
submitted alongside it: a batch of one scores 1.0 on every criterion, small
batches get very coarse scores, and two runs cannot be compared.

It also breaks evaluation E1, which compares `tier1_score` distributions across
three separately-scored groups — under batch ranking those distributions are
identical uniform distributions by construction unless all three groups are
pooled into one scoring batch.

Implemented instead: each criterion is scored as a corpus percentile
(`criteria.X.score`, absolute and batch-independent), the composite is their
weighted mean, and the submitted batch determines only the ordinal `rank`.
`batch_rank_fraction` is stored separately for display. The literal spec
behaviour remains available via `tier1.composite.method: batch_rank_aggregation`.

This preserves the spec's stated intent — *"ranking is computed within the
submitted batch, percentile scores are computed against the corpus"* — which
the snippet contradicts.

### 2. Shuffle pass uses a Monte-Carlo p-value, not a 0.95 percentile (§4.8)

Beating 19 of 20 controls is a descriptive 95% win rate, not a 0.05-level test:
with the candidate pooled among M controls, the top two positions already
account for 2/(M+1) ≈ 9.5% of the ranks. `shuffle_pass` is now

```
p = (1 + #{control >= real}) / (M + 1)      pass iff p <= alpha
```

With M=20 and α=0.05 that requires beating all 20 controls (p = 0.048). Both
`percentile` and `p_value` are stored. Config validation rejects an α that the
configured M can never attain. **Publication runs should use `n_shuffles: 99`.**

### 3. Structural sub-score excludes GC and length (§4.8)

Dinucleotide shuffling preserves both exactly, so they contribute an identical
constant to a candidate and all of its controls, compressing the margin without
adding information.

### 4. The artifact embeds SVGs instead of referencing paths (§7.2)

A JSON file containing `svg_path` entries pointing outside itself is not
self-contained. Diagrams are embedded as SVG text; `output.embed_svg: false`
disables them.

### 5. Additional artifact provenance

Beyond §7.2: `artifact_schema_version`, git commit + dirty flag, input file
SHA-256, corpus hashes and `is_placeholder` / `publication_eligible`, the
scoring signature, per-candidate `batch_rank_fraction`, `rank_min` /
`rank_is_tied`, and shuffle `p_value` / `complete` / `n_unique`. Unavailable
values are `null` — NaN is not valid JSON.

### 7. Silent-failure guards

Every one of these turns a class of silent wrong answer into a loud one:

* Each shuffle is checked to preserve the **k-mer multiset**, not merely the
  letter multiset — that is the property §4.8 actually claims.
* A missing criterion weight is an error, never an implicit zero: a typo would
  otherwise change what the composite means while still producing a number.
* Non-finite feature or reference values are rejected before scoring. NaN sorts
  to the top of a percentile and would read as a perfect candidate.
* The corpus cache key includes the full ViennaRNA model settings
  (temperature, dangles, noLP, uniq_ML, …), not just tool versions — those can
  change every cached number while leaving version strings identical.
* Sampled structures are validated for length, balance and alphabet before
  forgi sees them; a dropped sample would bias the loop-size distribution that
  is Tier 2's only input.
* Per-candidate seeds are derived from `(run_seed, candidate_id, sequence,
  purpose)`, so results do not depend on process-pool scheduling and a reused
  input id cannot collide two sequences onto one random stream.

### 6. Missing corpus is a hard failure, and provenance is required to publish

There is no automatic fallback. `--development-corpus` is the only way to run
without real reference data, and it marks the artifact ineligible for
publication and prints a warning on every run.

A run is `publication_eligible` only if its reference library has a
`<name>.manifest.json` recording `source`, `curator` and `curated_date`.
Correct columns are not provenance: the software cannot tell a curated corpus
from somebody's scratch file, and a pipeline that calls both publishable is
lying about the one thing this project is trying to be careful about. The
artifact's `development_reasons` names whichever condition failed.

### 8. `d_pocket` is a robust extent, not twice the RMS spread (§5.4)

The spec's `2 * singular_value / sqrt(n)` is twice the RMS spread of the
alpha-sphere centres — a dispersion statistic, not a cavity width. We compute a
5–95% quantile extent along each principal axis, projecting sphere *surfaces*
(`z ± r`) rather than centres, since what a loop must reach across is the
cavity rather than the set of points at its middle. The RMS spreads, the
centres-only extents and the equivalent-sphere diameter are all stored beside
it, and a cavity whose extent and equivalent diameter disagree by more than 2×
is flagged rather than silently scored.

### 9. The flexible descriptor is primary, contour length is the check (§5.5)

`6 Å × L` is contour length: it assumes a fully extended loop and therefore
over-penalises long loops, which bend. The √L flexible-chain proxy is primary;
contour length is reported as the upper-bound sensitivity analysis. Neither
`flex_c` nor `sigma` is ever tuned on E3.

### 10. Bands come from a fixed control bank, not the submitted batch (§5.8)

The spec derives band thresholds from the shuffles of whatever batch was
submitted, which makes a candidate's band depend on its co-submissions and
leaves the thresholds undefined for small batches. Instead a fixed bank of
corpus-derived dinucleotide shuffles is folded once and reused for every target
and batch; only the pocket arithmetic is redone per target.

Bands are defined on the **control percentile of the absolute mismatch**, not
on the Gaussian score — the Gaussian underflows at large mismatches and depends
on `sigma`, and a display parameter must not be able to move a candidate
between bands. `tier2_control_percentile` is also what evaluation E3 compares,
since raw scores are not comparable across targets.

### 11. Mode C's charge term is target-level, and cannot re-band anything (R§5.2)

The refinements spec combines size/coverage and charge complementarity into one
surface-mode score. Implemented as specified — with one property made explicit
everywhere the number appears, because it is easy to misread.

The patch potential is a property of the **target**. It is identical for every
candidate scored against that target, so it can shift the reported agreement
score but can never reorder two candidates and can never move one between bands.
`disagreement` — the quantity the control percentile and therefore the band is
computed from — is the area mismatch alone.

This is the same reasoning that already keeps `sigma` out of the banding: a term
that cannot discriminate between candidates must not be able to look as if it
did. The record says so in its own text (`charge_is_target_level`,
`charge_note`), the dashboard says so under the metric, and the explanation
sentence attributes the charge to the site rather than to the candidate.

### 12. The surface footprint is read from the fold, not from the length (R§4.2, R§6)

The Band 2 model was `(a_per_nt_ss × footprint_scale)² × length` — each
nucleotide assumed to cover a ~36 Å² square. R§6 asks for it to be replaced by a
radius-of-gyration proxy from the forgi element graph, and it has been.

`radius_of_gyration_A` is the textbook mass-weighted definition evaluated on the
secondary structure: each element carries its nucleotides, distances between
elements are measured along the structure graph, stems are rigid A-form rods and
unpaired regions are Gaussian segments. `footprint_area_A2 = π·(Rg·scale)²`.
It reproduces `L/√12` for a pure helix and `b√(N/6)` for an unstructured strand,
and both anchors are pinned in `tests/test_rna_size.py`.

What it buys: two 100-nt candidates are no longer the same size. A single long
helix (Rg 40.9 Å) and a four-way junction (Rg 26.3 Å) differ 2.4× in footprint
area where the length model called them identical.

Why it is the default rather than an option — measured, not assumed. Against the
2519 Å² IGFBP3 patch the length model's optimum is N ≈ 70 nt, which is exactly
the longest candidate in the demo batch: every "strong" call it makes is one of
the longest sequences, and its correlation between length and control percentile
is +0.976. It answers "which is longest?" while appearing to answer "which
fits?". See [docs/descriptor_validation.md](docs/descriptor_validation.md).

`tier2.surface.footprint_model: length` restores the old behaviour for
comparison. Both are coarse, and neither is calibrated against measured
interface areas; the band is a rank against shuffled controls either way.

### 13. The bundle id excludes fpocket's Monte-Carlo volume (R§4.1)

fpocket estimates cavity volume by Monte-Carlo integration (300 iterations by
default) and exposes no seed, so the same structure measured twice gives volumes
differing by a few percent — we observe up to ~4% across 17 cavities on 7WRQ.

The bundle id previously covered that number, which meant no two builds of the
same evidence could ever agree and the CI `--rebuild` check could not pass on a
real target. The volume and the three quantities derived from it are now
excluded from the id, **declared in the bundle** under `nondeterministic_fields`
rather than left as a silent convention, and the verifier measures the drift
between builds and fails if it exceeds 10%. Nothing scored depends on them:
`d_pocket_A`, which drives every pocket-mode band, comes from alpha-sphere
coordinates and is exactly reproducible.

Two related fixes fall under the same heading — recorded command lines and the
structure's cache path are reduced to filenames for hashing (where a tool ran is
not evidence about a protein), and pocket-overlap evidence is sorted, because
iterating a Python set gives an order that varies between processes.

### 14. E1 is scored out of fold, split by target (§8.2)

As specified, E1 is circular: the validated aptamers define the corpus
percentiles and are then the positive group. Folds are grouped by target, so no
held-out sequence is scored against its own target's family. A dinucleotide-
shuffled hard negative is added alongside the IID random control. See
[evaluation/README.md](evaluation/README.md).

---

## Targets

The user supplies an identifier, a chain and a binding mode. Everything else
happens server-side (R§3): fetch, chain selection, hetero policy, cavity
detection, patch measurement, electrostatics. The heavy, biology-literate step
moved off the biologist and into the pipeline, where §5.3–5.4's logic already
lived.

### Complexes, and the order that matters

For an interface target the binding partner sits *on* the surface being
measured, so the order is not arbitrary (R§3.4):

1. with the complex still loaded, confirm the configured binding-site residues
   and record which of them actually contact the partner;
2. strip the partner chains, exposing the interface;
3. measure the isolated target chain.

The partner says *where* to look; the partner itself must not be there when we
look. Chain removal is a one-pass filter on the chain label — the same operation
as water-stripping, not structural surgery.

Running this against 7WRQ reproduced the supplied interface list exactly: all 31
IGFBP3 residues within 4.5 Å of IGF1, independently computed.

The honest limitation is surfaced, not solved: a surface shaped by a bound
partner can relax once that partner is gone, so the measurement is of the
unbound state. Every multi-chain target carries that sentence into the artifact
and onto the screen.

### Predicted structures

AlphaFold models are accepted, resolved through the AFDB API rather than a
pinned URL (the `model_v4` suffix in the spec now 404s; current is v6). They are
recorded as `structure_source: predicted` everywhere, and surface mode warns
prominently: a predicted model may not show an interface that only forms when a
binding partner is present.

### The bundle

Per §5.2 every target-derived quantity is computed once per run — the
per-candidate part is arithmetic on two numbers — so the heavy tools are
isolated behind a **target bundle**: an immutable, checksummed JSON containing
every detected pocket, the alpha-sphere coordinates inline, the measured
binding-site patch, the full selection evidence, the partner-contact evidence,
the hetero policy actually applied, and tool provenance. Scoring consumes the
bundle and needs no external tools.

A bundle records the mode it was measured for, and refuses to be scored in the
other one: a surface bundle has no cavity and a pocket bundle has no patch, so
the alternative to refusing is comparing against whatever number happens to be
present. Both can exist side by side for one protein — the mode is in the
filename.

Bundles are also produced by `.github/workflows/target-bundle.yml` on a pinned
`ubuntu-24.04` runner. That workflow *is* the reproducibility story: it
recomputes every stored geometric quantity from the inline spheres, re-derives
the bundle id, and rebuilds the whole bundle to confirm the scientific payload
is stable, before uploading it. Anyone can re-run it and check the id. See
deviation 13 for what "stable" can and cannot mean when one of fpocket's numbers
is a Monte-Carlo estimate.

`bundle_id` is a SHA-256 over the reproducible payload **plus the tool identity
that produced it** — the fpocket version, command line and exit code are inside
the hash, so they cannot be edited without detection; timestamps, git state, CI
run ids, working directories and fpocket's Monte-Carlo volume estimate are
outside it, so two honest builds of the same evidence agree. What the id does
not cover is listed in the bundle itself under `nondeterministic_fields`.

### Testing fpocket without fpocket

The parser is pinned against hand-built fixtures covering CRLF, negative and
scientific numbers, inconsistent label spacing, five-digit residue numbers,
missing required fields, sphere-count disagreements and missing files — every
one of which must raise rather than produce a plausible wrong number.

Hand-built fixtures only prove the parser understands the format we *imagined*,
so the CI job also runs the real tool and packages its complete output as a
downloadable fixture archive. **Commit that archive under
`tests/fixtures/fpocket/real/` once per fpocket version bump** — an expiring
Actions artifact is not durable auditability.

Local development uses synthetic bundles from `scripts/make_dev_target.py`,
flagged `synthetic: true`, which taints the run's publication eligibility.

---

## Development

```bash
.venv/Scripts/python -m pytest -q      # 255 tests
```

`tests/test_folding.py` pins the spec §4.5–4.6 fixture values
(mfe −17.9, mfe_norm −0.4475, ensemble defect 0.0092, and the exact forgi
element decomposition). If those move, something upstream changed and every
downstream number is suspect.

CI runs the suite on Linux and Windows across Python 3.10 and 3.11, including
the ushuffle source build, and smoke-tests the headless pipeline end to end.

### One bug worth knowing about

`ushuffle.Shuffler` stores a raw `char*` to the bytes it is constructed with and
does not hold a reference to them. Passing `sequence.encode()` inline lets that
temporary be collected, after which every subsequent `shuffle()` returns reused
memory. On Windows this surfaced loudly as a `UnicodeDecodeError` for 28% of
candidates — but bytes that happened to stay in the ASCII range would have
passed silently and quietly invalidated every §4.8 statistic. The buffer is now
held, and every shuffle is checked to preserve the k-mer multiset.
