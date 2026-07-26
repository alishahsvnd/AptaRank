# Implementation Specification

## RNA Aptamer Candidate Ranking Tool

---

# 1. Overview

**Working title:** AptaRank — a two-tier evaluation and ranking framework for foundation-model–generated RNA aptamers. *(Name is a placeholder; final name TBC.)*

**What we are building.** A batch pipeline plus a dashboard. The user supplies (i) a file of candidate RNA sequences and (ii) an identifier for a target protein. The system scores every candidate, returns a ranked shortlist, and shows the evidence behind each ranking in an interpretable form.

**Why it is needed.** Generative models can now produce tens of thousands of candidate RNA sequences on demand. Each candidate can only be *confirmed* to work by a laboratory experiment, which is slow and expensive. This creates a bottleneck in deciding which handful of the many generated candidates are worth testing. This tool computes interpretable structural and geometric evidence for every candidate so a biologist can prioritise between them.

**Target venue.** IEEE ICDM 2026, Demo Track. 4 pages, IEEE 2-column, single-blind.

**Dates.**


| Date            | Milestone                                                              |
| --------------- | ---------------------------------------------------------------------- |
| 20 Jul 2026     | Kickoff, work allocation                                               |
| 11 Aug 2026     | Implementation complete enough to produce paper figures; writing paper |
| 13 Aug 2026     | Full draft to supervisor                                               |
| 19 Aug 2026     | Submit                                                                 |
| **20 Aug 2026** | **Hard deadline**                                                      |


---



# 2. Background and Related Work



## 2.1 The five things you need to know about the biology

You do not need to understand RNA chemistry. You need these five facts, which fully determine the computation:

1. **An RNA sequence is a string over the 4-letter alphabet** `{A, C, G, U}`**.** Our candidates are 20–100 characters long. That is the entire input data type.
2. **The string folds back on itself.** Certain character pairs stick together (`G`–`C`, `A`–`U`, and more weakly `G`–`U`). A folded RNA is therefore described by a *pairing* of string positions. The standard encoding is **dot-bracket notation**: a string of the same length where `.` means "this position is unpaired" and matched `(` / `)` mean "these two positions are paired with each other."
  ```
   sequence:   GUUCCAUGGGCCUUGACUUGCUGUGUCAUCACCAUGGGAC
   structure:  (((((((((....((((.......))))...)))))))))
  ```
   Predicting this pairing from the sequence is a solved, fast, deterministic computation. We call a library (ViennaRNA) to do it.
3. **An "aptamer" is an RNA whose folded shape sticks to a specific target molecule** — for us, always a protein. That is the whole definition. Useful for diagnostics and therapeutics because binding to a target lets you detect or block it.
4. **A protein target is a 3D structure**, distributed as a `.pdb` text file of atom coordinates from a public database. The relevant feature for us is that its surface has **cavities** (pockets) — concave regions where another molecule can sit. Detecting and measuring these is also a solved computation; we call a tool (fpocket) to do it.
5. **We cannot reliably predict binding.** There is very little labelled binding data, and the tools that claim to predict binding are unreliable. This is the single most important constraint on the design and is why the tool is built the way it is: **we compute interpretable evidence, not a binding prediction.** Anything target-related is presented as *plausibility*, never as proof. This must be reflected in the code (naming, output fields, UI copy), not just in the paper.



## 2.2 Prior work

Existing computational approaches to prioritising aptamer candidates fall into four families:


| Family                         | Approach                                                                                                         | Limitation for us                                                                                         |
| ------------------------------ | ---------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------- |
| Sequence clustering / counting | Group candidates by sequence similarity, rank by how often each appeared in a lab enrichment experiment          | Requires experimental read counts. Generated sequences have none.                                         |
| Structural motif tools         | Add secondary structure prediction (ViennaRNA) on top of clustering, rank by combined sequence + stability score | Still assumes lab-derived pools; ranks on stability alone, no target awareness                            |
| ML/DL binding predictors       | Learn from sequence pairs to predict interaction probability (e.g. RaptScore)                                    | Needs labelled binding data, which does not exist for most targets of interest                            |
| Docking / molecular dynamics   | Physically simulate the aptamer–target interaction                                                               | Very expensive; high and poorly-characterised error on RNA; only viable for a handful of final candidates |




## 2.3 The gap we address

> No existing approach provides a **systematic, modular, target-swappable pipeline for prioritising *generated* aptamer candidates** that reports **biologist-interpretable evidence** and is **explicit about the limits of target-aware scoring** rather than presenting it as a binding prediction.

Three properties follow from this, and every implementation decision below serves at least one:

- **P1 — Target-swappable.** Changing the target must be a config change (one identifier), not a code change. Only Tier 2 may depend on the target.
- **P2 — Interpretable.** Every number shown to a user must be traceable to a named, explainable criterion. No opaque composite that cannot be decomposed.
- **P3 — Honest.** Target-aware output is labelled as plausibility. Tier 2 never overrides Tier 1 (see §6).

---



# 3. Architecture

Architecture

*Figure 1 — the ranking engine (A) and the dashboard (B). Numbered stages correspond to §4–§7 below.*

## 3.1 Flow

```
 (1) INPUTS         candidate sequences (.txt/.csv) + target PDB ID
      |
 (2) TIER 1         intrinsic quality — target-agnostic
      |             scores ALL candidates, produces the ranking
      |             [top ~1000 pass down]
 (3) TIER 2         target-aware plausibility — target-swappable
      |             scores the survivors only, produces badges
 (4) EVIDENCE       per-candidate record: scores, structures, explanation text
      |
 (5) OUTPUT         ranked shortlist (Tier 1 order, Tier 2 badges) + dashboard
```



## 3.2 System inputs


| Input               | Format                                                                               | Constraints                                        |
| ------------------- | ------------------------------------------------------------------------------------ | -------------------------------------------------- |
| Candidate sequences | `.txt` (one sequence per line) or `.csv` (required column `sequence`, optional `id`) | Alphabet `ACGU` after normalisation; length 20–100 |
| Target identifier   | 4-character PDB ID string, e.g. `3SPU`                                               | Must resolve in the RCSB PDB                       |
| Chain selector      | string, optional                                                                     | Defaults to the first protein chain                |
| Config              | YAML/JSON                                                                            | Thresholds, weights, sample counts — see §9        |




## 3.3 System outputs

A single self-describing **run artifact** (see §7.2 for the schema), containing:

- ranked candidate table with per-criterion score breakdown
- per-candidate secondary structure (dot-bracket + rendered diagram)
- Tier 2 band per candidate + target-level evidence
- shuffled-control outcome per candidate
- rule-based natural-language explanation per candidate
- run-level metadata: config used, tool versions, target info, timings



## 3.4 Tools and libraries


| Tool                   | Role                                                                                                                     | Interface               | Install                                    |
| ---------------------- | ------------------------------------------------------------------------------------------------------------------------ | ----------------------- | ------------------------------------------ |
| **ViennaRNA**          | RNA secondary structure prediction: dot-bracket, stability energies, ensemble statistics, stochastic sampling, SVG plots | Python bindings         | `pip install ViennaRNA` (verified: 2.7.2)  |
| **forgi**              | Parses dot-bracket into named structural elements (stems, loops)                                                         | Python                  | `pip install forgi` — see gotcha in §4.5   |
| **ushuffle**           | Dinucleotide-preserving sequence shuffling, for negative controls                                                        | Python                  | `pip install ushuffle` (verified: 1.1.2)   |
| **fpocket**            | Detects and measures cavities on a protein structure                                                                     | CLI, parse output files | `apt install fpocket` or build from source |
| **PDB2PQR + APBS**     | Assigns atomic charges and computes an electrostatic potential grid for the target                                       | CLI                     | `pip install pdb2pqr`; APBS binary         |
| **Biopython**          | PDB file parsing, chain selection, coordinate handling                                                                   | Python                  | `pip install biopython`                    |
| numpy / scipy / pandas | Numerics, ranking, statistics                                                                                            | Python                  | standard                                   |


**Not used, deliberately:** docking software, molecular dynamics, learned binding-affinity predictors, 3D RNA structure prediction. Each is either out of scope for the timeline or introduces error we cannot characterise. 

---



# 3.5 Datasets



## 3.5.1 Reference aptamer corpus

- **What:** ~1,000-5,000 experimentally validated RNA aptamer sequences, curated by Laura (will provide clean data ASAP). Data has known protein targets (sequence ↔ PDB ID pairs) — this is what makes §8 evaluation possible.
- **Format:** CSV — `id, sequence, target_name, target_pdb_id (nullable), source_reference`.
- **Role:** defines the **reference distributions** that Tier 1 scores against. We do not use hard-coded thresholds anywhere; "is this candidate normal?" always means "where does it sit in the corpus distribution?"
- **Delivery:** Laura provides, by end of this week. Tier 1 percentile scoring is blocked without it. 

**Precompute once, cache to disk:** run the full Tier 1 feature extraction (§4) over the corpus and store the resulting feature matrix. Every candidate run then reads this cached matrix instead of recomputing ~5,000 folds. This is the single most important performance decision in Tier 1.

## 3.5.2 Protein Data Bank (PDB)

- **What:** public repository of experimentally determined 3D protein structures.
- **Access:** `https://files.rcsb.org/download/{PDB_ID}.pdb`
- **Role:** supplies the target structure for Tier 2.
- **Caching:** cache downloaded structures and all derived artifacts (fpocket output, APBS grid) keyed by `{pdb_id}_{chain}_{config_hash}`. Target-derived quantities do not depend on candidates, so they are computed **once per run**, not once per candidate.



## 3.5.3 Demo targets


| Target                 | PDB ID      | Role                                           |
| ---------------------- | ----------- | ---------------------------------------------- |
| NDM-1                  | `3SPU`      | Primary demo target                            |
| KPC-2                  | TBC (Laura) | Second target, to evidence target-swappability |
| Corpus-derived targets | various     | Evaluation (§8), labelled aptamer–target pairs |


---



# 4. Tier 1 — Intrinsic Quality (target-agnostic)



## 4.1 What Tier 1 is for

Tier 1 answers: *is this a well-formed, well-behaved aptamer-like sequence at all?* It never looks at the target. It produces **the ranking** — the primary ordering the user sees.

**Inputs:** candidate sequences; cached corpus feature matrix.
**Outputs, per candidate:** feature values, per-criterion percentile scores, composite Tier 1 score, rank, shuffled-control outcome, ensemble loop-size statistics (consumed by Tier 2).

## 4.2 Components


| #   | Component                                 | Library   | Produces                                                         |
| --- | ----------------------------------------- | --------- | ---------------------------------------------------------------- |
| 4.3 | Ingest & validation                       | —         | clean sequence table                                             |
| 4.4 | Composition features                      | —         | length, GC fraction                                              |
| 4.5 | Secondary structure & ensemble statistics | ViennaRNA | dot-bracket, normalised MFE, ensemble defect, positional entropy |
| 4.6 | Structural element features               | forgi     | stem fraction, loop counts and sizes                             |
| 4.7 | Ensemble loop-size distribution           | ViennaRNA | median / p90 loop size → **input to Tier 2**                     |
| 4.8 | Shuffled-control check                    | ushuffle  | pass/fail + margin                                               |
| 4.9 | Composite score                           | —         | Tier 1 score and rank                                            |


---



## 4.3 Step 1 — Ingest and validation

1. Read `.txt` (one sequence per line) or `.csv` (column `sequence`; optional `id`).
2. Normalise: strip whitespace, uppercase, replace `T` → `U`. *(Why: the same molecule is written with* `T` *in DNA-style files and* `U` *in RNA-style files. Same data, different convention.)*
3. Validate: alphabet is exactly `{A,C,G,U}`; length in `[20, 100]`. Reject with a per-row reason; surface rejects in the run artifact, never silently drop.
4. Deduplicate exact matches; retain `duplicate_count`.
5. Assign stable `candidate_id` if not supplied.

**Output:** `DataFrame[candidate_id, sequence, length, duplicate_count]`.

---



## 4.4 Step 2 — Composition features

```python
gc_fraction = (seq.count("G") + seq.count("C")) / len(seq)
```

Pure character counting — no biology needed. Compute `length` and `gc_fraction`.

**Scoring these two is different from everything else.** For GC and length we do **not** want "higher is better" — we want *typical of real aptamers*. Score with a **two-sided percentile**:

```python
def two_sided_score(x, reference_values):
    """1.0 at the centre of the reference distribution, → 0 at either tail."""
    p = empirical_cdf(reference_values, x)     # in [0,1]
    return 1.0 - 2.0 * abs(p - 0.5)
```

For the directional criteria in §4.5 use a **one-sided percentile** instead (higher = better after sign correction).

---



## 4.5 Step 3 — Secondary structure and ensemble statistics (ViennaRNA)

**What this step does.** Given the sequence string, compute the predicted pairing (dot-bracket) and three numbers describing how *confidently* the sequence folds into that one shape. A sequence that folds into a single well-defined shape is more likely to behave predictably; a sequence that flops between many shapes is a poorer candidate.

**Verified API** (ViennaRNA 2.7.2):

```python
import RNA

md = RNA.md()
md.uniq_ML = 1                       # REQUIRED — without this §4.7 sampling returns nothing
fc = RNA.fold_compound(seq, md)

mfe_struct, mfe = fc.mfe()           # str (dot-bracket), float (kcal/mol, negative)
fc.exp_params_rescale(mfe)           # numerical stability before partition function
pf_struct, efe = fc.pf()             # must be called before ensemble_defect / bpp
ed  = fc.ensemble_defect(mfe_struct) # float in [0,1]
pe  = fc.positional_entropy()        # list, length n+1, valid indices 1..n
```

**Features to store:**


| Feature                   | Computation                      | Direction    | Meaning in plain terms                                                                       |
| ------------------------- | -------------------------------- | ------------ | -------------------------------------------------------------------------------------------- |
| `mfe`                     | `fc.mfe()[1]`                    | —            | stability energy, more negative = more stable                                                |
| `mfe_norm`                | `mfe / length`                   | lower better | **use this, not raw MFE** — raw MFE grows with length, so ranking on it just ranks by length |
| `ensemble_defect`         | `fc.ensemble_defect(mfe_struct)` | lower better | 0 = folds into exactly one shape; 1 = no consistent shape                                    |
| `positional_entropy_mean` | `mean(pe[1:n+1])`                | lower better | average per-position folding uncertainty                                                     |


**Sanity values** for the worked example sequence in §2.1: `mfe = -17.9`, `mfe_norm = -0.448`, `ensemble_defect = 0.0092`. Use this as a regression test fixture.

**Ordering constraint:** `fc.pf()` must be called before `ensemble_defect()` or `bpp()`, or you get garbage. Enforce this in a single wrapper function; do not call the raw API from multiple places.

**Performance:** ~1–5 ms per sequence. Fully parallel across candidates — use a process pool, not threads (the bindings are C extensions and release the GIL inconsistently).

---



## 4.6 Step 4 — Structural element features (forgi)

**What this step does.** The dot-bracket string is a flat encoding; forgi parses it into named *elements*. This matters because domain experts told us that the interpretable features are the **loops** — the unpaired regions — since those are the parts that physically stick out and can contact a target. So we need loops as first-class objects, not just a count of `.` characters.

Element type codes:


| Code     | Element       | Plain meaning                                           |
| -------- | ------------- | ------------------------------------------------------- |
| `s`      | stem          | a run of consecutive paired positions (a "rung ladder") |
| `h`      | hairpin loop  | unpaired run closed by a stem at one end                |
| `i`      | interior loop | unpaired region between two stems                       |
| `m`      | multiloop     | junction where 3+ stems meet                            |
| `f`, `t` | 5′ / 3′ tails | unpaired run at either end of the string                |


**Verified API:**

```python
from forgi.graph.bulge_graph import BulgeGraph

bg = BulgeGraph.from_dotbracket(db, seq=seq)
bg.defines                 # {'s0': [1,9,32,40], 'h0': [18,24], 'i0': [10,13,29,31], ...}
bg.element_length('h0')    # 7   -> nucleotides in the element
bg.stem_length('s0')       # 9   -> base PAIRS (stems only; element_length gives 18 = 9*2)
bg.to_element_string()     # 'sssssssssiiiisssshhhhhhhssssiiisssssssss'
```

⚠️ **Install gotcha:** `import forgi` emits a NumPy 2.x compatibility traceback from its optional 3D submodule. It is non-fatal and we do not use the 3D module. Import the specific module (`from forgi.graph.bulge_graph import BulgeGraph`) and suppress the warning, or pin `numpy<2` if it proves noisy.

**Features to store:**


| Feature                                   | Computation                                            |
| ----------------------------------------- | ------------------------------------------------------ |
| `n_hairpins`, `n_interior`, `n_multiloop` | count of `h*`, `i*`, `m*` keys in `bg.defines`         |
| `stem_fraction`                           | (total nucleotides in `s*` elements) / length          |
| `longest_stem_bp`                         | `max(bg.stem_length(e) for e in stems)`                |
| `max_loop_nt`                             | `max(bg.element_length(e) for e in h*/i*/m* elements)` |
| `total_unpaired`                          | count of `.` in dot-bracket                            |


Score `stem_fraction` two-sided against the corpus (§4.4) — both extremes are bad: all-stem means nothing is available to contact a target, all-loop means no stable shape.

---



## 4.7 Step 5 — Ensemble loop-size distribution (feeds Tier 2)

**Why this exists.** The single MFE structure is one prediction. Real molecules move between shapes. Rather than trusting one structure's loop sizes, sample many plausible structures and report the *distribution* of loop sizes. Tier 2's geometric comparison then consumes a distribution statistic instead of a point estimate.

```python
# md.uniq_ML = 1 must have been set at fold_compound construction (§4.5)
fc.pf()
samples = fc.pbacktrack(N_SAMPLES)     # N_SAMPLES = 100 (config)
```

For each sampled dot-bracket: parse with forgi, take `max_loop_nt`. Across all samples compute:


| Feature          | Use                                                                              |
| ---------------- | -------------------------------------------------------------------------------- |
| `loop_nt_median` | primary input to Tier 2 geometric compatibility                                  |
| `loop_nt_p90`    | conservative upper estimate, shown in the UI                                     |
| `loop_nt_iqr`    | dispersion — a wide IQR means the loop size is itself uncertain, worth surfacing |


**Regression test:** with `md.uniq_ML = 0` (the default), `pbacktrack` returns an empty tuple with no error. Assert `len(samples) == N_SAMPLES` and fail loudly.

---



## 4.8 Step 6 — Shuffled-control check

**Why this exists.** A high structural score could simply reflect letter composition — GC-rich strings fold tightly regardless of arrangement. To show the score reflects *arrangement* rather than *composition*, we compare each candidate against shuffled versions of itself that preserve composition.

```python
import ushuffle
shuffled = ushuffle.shuffle(seq.encode(), 2).decode()   # k=2 preserves dinucleotide composition
```

*(k=2 means every adjacent letter-pair frequency is preserved, which is a stricter and more honest control than just preserving single-letter frequencies.)*

Procedure per candidate:

1. Generate `M_SHUFFLES = 20` shuffles (config).
2. Run §4.5 + §4.6 on each — **not** §4.7 (too expensive; sampling is only needed for real candidates).
3. Compute the structural sub-score (the §4.9 composite restricted to the structural criteria) for each shuffle.
4. `shuffle_percentile` = fraction of shuffles the real candidate beats.
5. `shuffle_pass` = `shuffle_percentile >= 0.95`.
6. `shuffle_margin` = real score − median(shuffle scores) — used in the explanation text.

**Cost check:** 20 shuffles × 2 folds each ≈ 40 extra folds per candidate. At ~3 ms/fold that is ~0.12 s per candidate — fine for 12k candidates on the H200 box in parallel, but it dominates Tier 1 runtime. Make `M_SHUFFLES` configurable and expose a `--fast` mode that skips shuffling for exploratory runs.

---



## 4.9 Step 7 — Tier 1 composite score by rank aggregation

**Method: rank aggregation with equal weights.** Not a weighted sum of raw values — the features are on wildly different scales (energies in the tens, defects in [0,1], counts as small integers), and a weighted sum on raw values is dominated by whichever feature happens to have the largest numeric range. Ranking removes scale entirely and is robust to outliers.

**Criteria** (all converted to "higher is better"):


| #   | Criterion              | Source                    | Transform                                                        |
| --- | ---------------------- | ------------------------- | ---------------------------------------------------------------- |
| 1   | Structure stability    | `mfe_norm`                | one-sided percentile, sign-flipped (lower energy → higher score) |
| 2   | Fold definition        | `ensemble_defect`         | one-sided percentile, sign-flipped                               |
| 3   | Fold certainty         | `positional_entropy_mean` | one-sided percentile, sign-flipped                               |
| 4   | Structural composition | `stem_fraction`           | two-sided vs corpus                                              |
| 5   | Sequence composition   | `gc_fraction`             | two-sided vs corpus                                              |


```python
# per criterion: rank within the candidate batch, normalise to [0,1]
r_i = scipy.stats.rankdata(values_i) / len(values_i)
tier1_score = sum(w_i * r_i for i in criteria) / sum(w_i)   # w_i = 1.0 default
```

**Requirements:**

- Weights live in config and default to equal. Do **not** fit weights to the evaluation data — see §8.4.
- Persist every per-criterion score alongside the composite. The composite alone is useless to the UI (P2, interpretability).
- Ranking is computed within the submitted batch. Percentile scores are computed against the **corpus**. Do not confuse the two — they answer different questions ("best of what you gave me" vs "normal compared to known real aptamers").
- `shuffle_pass` is reported but **not** part of the composite — it is a diagnostic about the score, not a component of it.

**Output:** full-batch table sorted descending by `tier1_score`, with rank column. Top `N_TIER2 = 1000` (config) are marked for Tier 2.

---



# 5. Tier 2 — Target-Aware Plausibility (target-swappable)



## 5.1 What Tier 2 is for

Tier 2 answers: *given this specific target, is this candidate's shape geometrically plausible as something that could sit in the target's cavity?*

**This is not a binding prediction and must never be described, named, or displayed as one.** Field-leading domain experts advised explicitly that target-aware computational signals of this kind carry substantial uncertainty and should be treated as one signal among several, not as evidence of binding. That advice is the reason this tier is graded into three bands rather than reported as a number-that-looks-like-a-probability, and the reason it cannot reorder the Tier 1 ranking (§6).

**Inputs:** top ~1000 Tier 1 survivors (with their `loop_nt_median` from §4.7); target PDB ID + chain. (for discussion: should we input all the tier 1 survivors? depends on the compute/also on what biologists actually want)

**Outputs:** per-candidate geometric compatibility score and band; target-level pocket description and electrostatic flag.

## 5.2 The key efficiency insight

Split the work by what it depends on:


| Computed                                  | Depends on         | When                             |
| ----------------------------------------- | ------------------ | -------------------------------- |
| Structure fetch, chain selection, cleanup | target only        | **once per run**                 |
| fpocket cavity detection and measurement  | target only        | **once per run**                 |
| PDB2PQR + APBS electrostatics             | target only        | **once per run**                 |
| Geometric compatibility score             | target × candidate | per candidate (cheap arithmetic) |


Only the last row scales with candidate count, and it is arithmetic on two numbers. Tier 2 is therefore *fast* despite involving heavier tools — provided you structure it this way. Do not put fpocket or APBS inside a per-candidate loop.

---



## 5.3 Step 1 — Fetch and prepare the target structure

1. Download `https://files.rcsb.org/download/{PDB_ID}.pdb`; cache by ID.
2. Parse with Biopython (`Bio.PDB.PDBParser`).
3. Select the chain: config-specified, else the first protein chain. *(A PDB file often contains several copies of the same protein; we only need one.)*
4. Clean: remove water molecules (residue name `HOH`). **Keep other non-protein atoms by default** — for our demo targets these include a metal ion that sits in the functional cavity, and removing it changes the cavity geometry. Make this a config flag `strip_hetatm: false`.
5. Write the cleaned single-chain structure to a working `.pdb`.

**Validation:** assert the output has > 100 atoms and at least one chain; fail with a clear message if the PDB ID resolved to something unexpected (e.g. a nucleic-acid-only structure).

---



## 5.4 Step 2 — Cavity detection and measurement (fpocket)

**What fpocket does.** It fills the protein's surface concavities with spheres ("alpha spheres") and reports each resulting pocket with a score and geometric descriptors. We use it purely as a measuring instrument: *where are the cavities and how big are they?*

```bash
fpocket -f target_clean.pdb
```

Produces `target_clean_out/`:


| File                         | Contents                                                                                                 |
| ---------------------------- | -------------------------------------------------------------------------------------------------------- |
| `target_clean_info.txt`      | Per-pocket summary: Score, Druggability Score, Number of alpha spheres, Volume, mean alpha-sphere radius |
| `pockets/pocket{N}_atm.pdb`  | Protein atoms lining pocket N                                                                            |
| `pockets/pocket{N}_vert.pqr` | Alpha-sphere centre coordinates for pocket N                                                             |


**Parsing:** `_info.txt` is a blank-line-delimited block format, one block per pocket. Write a small parser returning `list[dict]`; unit-test it against a checked-in fixture file, since the format is whitespace-sensitive and will silently give wrong numbers if parsed loosely.

**Computing the characteristic dimension.** fpocket gives Volume, but volume alone does not tell us whether a loop can reach in. Compute an extent from the alpha-sphere centres in `pocket{N}_vert.pqr`:

```python
coords = parse_pqr_coordinates(f"pockets/pocket{N}_vert.pqr")   # (k, 3) array
centred = coords - coords.mean(axis=0)
_, sv, _ = np.linalg.svd(centred, full_matrices=False)
extents = 2 * sv / np.sqrt(len(coords))          # principal-axis extents, Å
d_pocket = float(extents[0])                     # characteristic dimension
d_equiv  = 2 * (3 * volume / (4 * np.pi)) ** (1/3)   # equivalent-sphere diameter, cross-check
```

Store both. If they disagree by more than ~2× the pocket is oddly shaped — flag it in the run artifact rather than silently proceeding.

**Pocket selection — do not fully automate this.** fpocket returns many cavities and ranks them by its own score, which is not guaranteed to identify the functionally relevant one. Selection order:

1. If config supplies `target.active_site_residues` (a list of residue numbers), select the pocket whose lining residues overlap that list most. **This is the default path for our demo targets** — Laura supplies literature-confirmed residue numbers for NDM-1 and KPC-2.
2. Otherwise fall back to the highest fpocket Score, and set `pocket_selection: "automatic"` in the artifact so the UI can caveat it.

**[DISCUSS]** If fpocket's cavity list turns out ambiguous for KPC-2, a fallback is to rescore its pockets with P2Rank. Do not build this now; flag it if it happens.

---



## 5.5 Step 3 — Aptamer geometric descriptor

We need the candidate's loop size in the same physical units as the pocket (Ångström, Å). From §4.7 we have `loop_nt_median` — a count of characters. Convert with an explicit, documented approximation:

```python
A_PER_NT = 6.0      # config; extended single-stranded rise per nucleotide, Å
d_apt_extended = A_PER_NT * loop_nt_median
d_apt_flexible = A_PER_NT * math.sqrt(loop_nt_median) * FLEX_C   # config, conservative variant
```

**Be explicit in the code comments and the docs about what this is.** It is a deliberately coarse, transparent geometric approximation with a known failure mode: it assumes an extended loop, so it *overestimates* reach for long loops, where a flexible chain's real end-to-end distance scales closer to √L. That is why both variants are computed.

**Why not predict the 3D structure instead?** We considered it and rejected it. RNA 3D structure prediction is least accurate precisely on flexible loop regions, is trained mostly on natural sequences (our candidates are model-generated and non-natural), and returns a single static snapshot with confidence scores that are not reliable indicators of accuracy on out-of-distribution input. Stacking that onto an already-approximate pocket measurement would produce a compatibility score with three poorly-characterised error sources instead of one nameable simplification. 

---



## 5.6 Step 4 — Geometric compatibility score

Compare the two dimensions with a smooth fit function so that near-misses degrade gracefully rather than falling off a cliff:

```python
SIGMA = 6.0   # Å, config
geometric_score = math.exp(-((d_apt - d_pocket) ** 2) / (2 * SIGMA ** 2))   # in (0, 1]
```

Compute for both `d_apt_extended` and `d_apt_flexible`; report the extended variant as primary and the flexible one as a sensitivity check in the run artifact.

Store `d_apt`, `d_pocket`, and the raw difference alongside the score — the explanation text (§7.3) quotes these directly, and a bare score is not interpretable (P2).

---



## 5.7 Step 5 — Electrostatic compatibility (target-level)

**What this is.** RNA carries a uniformly negative electrical charge along its backbone. A cavity whose surface is positively charged is therefore more electrostatically hospitable to *any* RNA; a negatively charged one is repulsive. Because this depends only on the target, it is a **single target-level flag, not a per-candidate score.**

```bash
pdb2pqr30 --ff=AMBER --apbs-input=target.in target_clean.pdb target_clean.pqr
apbs target.in          # writes an OpenDX potential grid, e.g. target_pot.dx
```

Then sample the grid at the selected pocket's alpha-sphere centre coordinates (from §5.4) via trilinear interpolation, and compute the mean potential in units of kT/e.


| Output                         | Meaning            |
| ------------------------------ | ------------------ |
| `electrostatic_mean_potential` | continuous value   |
| `electrostatic_compatible`     | `True` if mean > 0 |


**Keep this out of the Tier 2 composite.** It does not vary per candidate, so folding it into a per-candidate score would add a constant to every candidate and change nothing about the ordering while making the score harder to explain. Display it as a target-level badge in the dashboard (§7.1e).

*Stretch:* if APBS integration proves painful, ship the geometric score alone and report electrostatics as future work. It is confirmed in scope but it is the most droppable item here.

---



## 5.8 Step 6 — Bands from shuffled controls

The Tier 2 score is on an arbitrary (0,1] scale whose meaning depends on the target's pocket size. A fixed cutoff like "top 25%" or "score > 0.6" would be arbitrary and, worse, would behave inconsistently when the target changes — which would silently break target-swappability (P1).

Instead, derive thresholds empirically, **per run**:

1. Take the shuffled sequences already generated in §4.8 (reuse them — do not regenerate).
2. Run §4.7 and §5.5–5.6 on a sample of them (config `N_SHUFFLE_TIER2 = 2000`) against the same target.
3. Compute the 75th and 95th percentiles of the resulting shuffled Tier 2 score distribution.

```
score >= p95            -> "strong"
p75 <= score < p95      -> "moderate"
score <  p75            -> "weak"
```

**Requirements:**

- Recompute thresholds for every target. Never persist a threshold across targets.
- Store `p75`, `p95`, and the shuffled score distribution in the run artifact — the dashboard's scatter panel (§7.1f) draws these as reference lines, and the paper needs them.
- Three graded bands, not a binary flag: a hard in/out cutoff implies more precision than this signal has.

---



# 6. Ranking Across Tiers



## 6.1 The rule

> **Tier 1 determines the order. Tier 2 annotates it. Tier 2 never reorders anything.**



## 6.2 Why

If the two scores were merged into one number, a mediocre candidate with a flattering target-plausibility score could outrank a structurally excellent one. That would treat the *less* reliable signal as capable of overriding the *more* reliable one — exactly backwards, and directly contrary to the expert advice that motivated the design (P3). Keeping the tiers separate also means the primary ranking is stable when the target changes, which is what makes the target-swappability claim meaningful.

## 6.3 Implementation

```python
ranked = tier1_results.sort_values("tier1_score", ascending=False)     # THE ordering
ranked["rank"] = range(1, len(ranked) + 1)

survivors = ranked.head(N_TIER2)                                       # default 1000
survivors["tier2_score"] = geometric_scores
survivors["tier2_band"]  = assign_band(survivors["tier2_score"], p75, p95)

# rejoin; candidates below the Tier 2 cut get band = "not_evaluated" (NOT "weak")
final = ranked.merge(survivors[["candidate_id", "tier2_score", "tier2_band"]],
                     on="candidate_id", how="left")
final["tier2_band"] = final["tier2_band"].fillna("not_evaluated")
```

**"not evaluated" ≠ "weak".** A candidate outside the top 1000 has no Tier 2 evidence at all. Displaying that as "weak" would fabricate evidence. Use a visually distinct neutral state in the UI.

## 6.4 Required diagnostic — Spearman correlation

Compute and store `spearman(tier1_rank, tier2_score)` over the survivors, every run.

- **Low or near-zero correlation (expected):** the tiers measure different things, the annotation carries independent information, and the dashboard's scatter panel will show a meaningful spread.
- **High correlation:** Tier 2 is largely restating Tier 1, which would undermine the two-tier design. If this happens, tell Laura before building further on it — it changes what the paper can claim.

This is a five-line computation that directly determines whether a core design claim holds. Run it as soon as both tiers produce output — **do not leave it until evaluation week.**

---



# 7. Output and Dashboard

Dashboard

*Figure 2 — target dashboard layout. Panels (a)–(f) below.*

## 7.1 Panels


| Panel                            | Contents                                                                                                                    | Data source      |
| -------------------------------- | --------------------------------------------------------------------------------------------------------------------------- | ---------------- |
| **(a) Run configuration**        | Input filename, candidate count, target PDB ID, length bounds, seed, "Run pipeline" action                                  | Run metadata     |
| **(b) Ranked candidate list**    | Rank, truncated sequence, Tier 1 score, Tier 2 band, shuffle result. Selectable rows.                                       | §4.9, §5.8, §4.8 |
| **(c) Candidate detail**         | Sequence, dot-bracket, rendered structure diagram, per-criterion percentile bars, composite + badges                        | §4.5–4.9         |
| **(d) Explanation panel**        | Rule-based natural-language paragraph, evidence chips, plausibility caveat line                                             | §7.3             |
| **(e) Target panel**             | PDB ID and name, pocket visual, pocket count and volume, characteristic dimension, electrostatic badge, target swap control | §5.3–5.7         |
| **(f) Tier 1 vs Tier 2 scatter** | One point per survivor, coloured by band, with p75/p95 reference lines                                                      | §5.8, §6.4       |




## 7.2 Run artifact schema

Produce a single JSON artifact. **The pipeline must be fully usable headless** — the dashboard reads this file and holds no computation of its own. This decoupling is what lets the paper's figures be generated from the same runs the demo uses.

```jsonc
{
  "run_id": "...",
  "created_utc": "...",
  "config": { /* full resolved config, including all thresholds and weights */ },
  "versions": { "viennarna": "2.7.2", "forgi": "...", "fpocket": "...", "apbs": "..." },
  "input": { "filename": "...", "n_submitted": 12480, "n_valid": 12401, "n_rejected": 79,
             "rejections": [{ "row": 42, "reason": "invalid character 'N'" }] },
  "target": {
    "pdb_id": "3SPU", "chain": "A", "name": "...",
    "n_pockets": 7, "pocket_selection": "active_site_residues",
    "selected_pocket": { "index": 1, "volume_A3": 412.0,
                         "d_pocket_A": 18.4, "d_equiv_A": 9.2 },
    "electrostatic_mean_potential": 1.83, "electrostatic_compatible": true
  },
  "tier2_thresholds": { "p75": 0.41, "p95": 0.62, "n_shuffled_scored": 2000 },
  "diagnostics": { "spearman_tier1_tier2": -0.04, "runtime_seconds": { "tier1": 812, "tier2": 47 } },
  "candidates": [
    {
      "candidate_id": "c00001", "sequence": "GUUCC...", "length": 40,
      "rank": 1, "tier1_score": 0.87,
      "criteria": { "mfe_norm": {"value": -0.448, "score": 0.91},
                    "ensemble_defect": {"value": 0.0092, "score": 0.96},
                    "positional_entropy_mean": {"value": 0.13, "score": 0.88},
                    "stem_fraction": {"value": 0.65, "score": 0.78},
                    "gc_fraction": {"value": 0.58, "score": 0.64} },
      "structure": { "dot_bracket": "(((((...", "svg_path": "structures/c00001.svg" },
      "elements": { "n_hairpins": 1, "n_interior": 1, "max_loop_nt": 7,
                    "loop_nt_median": 7, "loop_nt_p90": 9, "loop_nt_iqr": 2 },
      "shuffle": { "pass": true, "percentile": 1.0, "margin": 0.31 },
      "tier2": { "score": 0.71, "band": "strong",
                 "d_apt_A": 21.0, "d_pocket_A": 18.4, "difference_A": 2.6 },
      "explanation": "Folds into a single dominant structure — ..."
    }
  ]
}
```



## 7.3 Explanation generation (rule-based)

**Rule-based, not an LLM call.** Deterministic, auditable, no inference dependency, and every sentence is guaranteed to be supported by a stored number.

Build as an ordered list of `(condition, template)` rules over the candidate record. Fire the top ~4 by priority and join into a paragraph:

```python
RULES = [
  (lambda c: c["criteria"]["ensemble_defect"]["score"] > 0.9,
   "Folds into a single dominant structure — ensemble defect {ensemble_defect:.3f} "
   "sits in the top {pct:.0f}% of the reference corpus."),

  (lambda c: c["tier2"]["band"] == "strong",
   "Largest accessible loop spans approximately {d_apt:.0f} Å, compatible with the "
   "{d_pocket:.1f} Å cavity detected on chain {chain} of {pdb_id}."),

  (lambda c: c["shuffle"]["pass"],
   "Outscores its dinucleotide-shuffled counterparts by {margin:.2f} on structural "
   "criteria, so the ranking is not driven by composition alone."),
  # ... including negative/cautionary rules
]
```

**Requirements:**

- Every emitted number must come from the stored record — never recompute or round differently in the text than in the table.
- Include cautionary rules, not only positive ones. A candidate that fails the shuffle check must say so.
- Always append the fixed caveat line: *"Tier 2 reflects geometric plausibility only — it is not evidence of binding."*

Laura: likely need to add to this section to provide enough detail with the rules. 

## 7.4 Structure diagrams

ViennaRNA renders these directly:

```python
RNA.svg_rna_plot(seq, db, f"structures/{candidate_id}.svg")   # verified working
```

Generate for the top `N_DIAGRAMS = 50` only (config). For the paper figure we used custom matplotlib rendering from `RNA.naview_xy_coordinates(db)` for finer styling control — the built-in SVG is fine for the dashboard.

## 7.5 **[DISCUSS]** How to build the dashboard

Discuss in a meeting. The decision is *only* about presentation — §7.2 guarantees the backend is unaffected either way.


| Option                                 | Pros                                                                         | Cons                                                                             |
| -------------------------------------- | ---------------------------------------------------------------------------- | -------------------------------------------------------------------------------- |
| **Streamlit**                          | Fastest to working; Python-only; used in comparable demo papers              | Less control over layout; the reference figure's density is awkward to reproduce |
| **FastAPI + React**                    | Full control; matches the mockup; closer to comparable systems in this venue | Meaningfully more work in a 3-week build                                         |
| Hybrid: FastAPI backend + Streamlit UI | Backend reusable if UI is rebuilt later                                      | Two things to maintain                                                           |


**Laura's steer:** the dashboard screenshot is the figure demo-track reviewers look at first, so visual quality matters more than it would for an internal tool — but not at the cost of the pipeline being finished. Recommend deciding based on your honest estimate of remaining capacity after Tier 1 and Tier 2 are working.

---



# 8. Evaluation

Evaluation runs produce the paper's quantitative claims. Everything here reads the §7.2 artifacts — no separate analysis pipeline.

## 8.1 Comparison groups


| Group                  | Definition                                                                                 | Role                | Labelled?                         |
| ---------------------- | ------------------------------------------------------------------------------------------ | ------------------- | --------------------------------- |
| **Random RNA**         | Generated strings, matched to the corpus length distribution and single-letter frequencies | Negative control    | n/a                               |
| **Validated aptamers** | Real sequences from the reference corpus, with known protein targets                       | Positive control    | **Yes — sequence ↔ target pairs** |
| **FM-generated**       | Fine-tuned RNAGenesis output                                                               | The actual use case | No                                |


Generate random RNA by sampling letters at the corpus's empirical single-letter frequencies, with lengths drawn from the corpus length distribution — this is important, since an unmatched control would let the tool "win" trivially on composition rather than structure.

## 8.2 Experiments

**E1 — Tier 1 discriminates real from random.**
Score all three groups; compare `tier1_score` distributions. Report medians and a Mann-Whitney U test for validated-vs-random. Expect validated > random. If FM-generated lands between them, that is an interesting and reportable result — do not treat it as a bug.
→ *Paper: violin plot.*

**E2 — Structure, not composition.**
Shuffled-control pass rate per group. Expect validated aptamers to pass at a high rate and random RNA near chance.
→ *Paper: bar chart or table.*

**E3 — Matched vs mismatched target (the strongest available test).**
This is the one place we have genuine labels. For each validated aptamer with a known target: score it against (i) its true target and (ii) K decoy targets drawn from other corpus entries. If Tier 2 carries real signal, true-target scores should exceed decoy scores.
Report the paired difference and an AUROC over true-vs-decoy.
→ *Paper: paired comparison plot or single AUROC figure.*
**Note:** this is the highest-value evaluation and also the most likely to return a null result. A null is publishable — it becomes an honest limitation and reinforces the "plausibility, not prediction" framing. **Run it early**, because it affects how Tier 2 is described throughout the paper.

**E4 — Target swappability.**
Score one fixed candidate set against NDM-1 and KPC-2. Compare band assignments and rank-order of `tier2_score`. Expect meaningfully different results — identical results would mean the target input is not actually influencing anything.
→ *Paper: side-by-side comparison.*

**E5 — Tier independence.**
Spearman correlation between Tier 1 rank and Tier 2 score (§6.4), reported per target.
→ *Paper: cited in text and visible in dashboard panel (f).*

## 8.3 Reproducibility

Fixed seeds for shuffling, sampling, and random-RNA generation; seeds recorded in the run artifact. Every figure regenerable from a stored artifact by a single script. Check the analysis notebooks into the repo — the paper cites a public repository, so the repo has to be presentable.

## 8.4 One methodological rule

**Do not fit the composite weights to the evaluation groups.** If weights are tuned to maximise separation between validated and random sequences and we then report that separation as a result, the result is inflated by construction and the reviewer criticism is fatal. Weights stay equal and configurable.

*Optional, only with a held-out split:* fitting weights on a held-out subset and showing they roughly agree with equal weighting is a legitimate robustness check. Treat as stretch, and only if corpus size allows a clean split.

---



# 9. Configuration Reference

All tunables live in one config file. Nothing in this list may be hard-coded.

```yaml
input:
  min_length: 20
  max_length: 100

tier1:
  n_shuffles: 20
  n_ensemble_samples: 100
  weights: { mfe_norm: 1.0, ensemble_defect: 1.0, positional_entropy_mean: 1.0,
             stem_fraction: 1.0, gc_fraction: 1.0 }
  shuffle_pass_percentile: 0.95

tier2:
  n_candidates: 1000
  strip_hetatm: false
  a_per_nt: 6.0
  flex_c: 1.0
  sigma_A: 6.0
  n_shuffled_for_thresholds: 2000
  band_percentiles: { moderate: 0.75, strong: 0.95 }
  target:
    pdb_id: "3SPU"
    chain: "A"
    active_site_residues: []      # supplied by Laura per target

output:
  n_diagrams: 50

seed: 42
```

---



# 10. Build Order and Ownership


| Week    | Deliverable                                                                                                                         | Definition of done                                                                                             |
| ------- | ----------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------- |
| **1**   | Repo skeleton, config loader, ingest/validation (§4.3), Tier 1 features (§4.4–4.6), corpus precompute + cache (§3.5.1)              | Tier 1 features computed for the full corpus and a test candidate batch; §4.5 fixture values reproduce exactly |
| **1–2** | Ensemble sampling (§4.7), shuffled controls (§4.8), composite + ranking (§4.9)                                                      | End-to-end Tier 1 producing a ranked table and a valid run artifact                                            |
| **2**   | Target prep (§5.3), fpocket integration + parser (§5.4), geometric score (§5.5–5.6), banding (§5.8), **Spearman diagnostic (§6.4)** | Tier 2 runs for NDM-1; bands assigned; diagnostic reported to Laura                                            |
| **2–3** | Electrostatics (§5.7, droppable), dashboard (§7), explanation rules (§7.3)                                                          | Dashboard renders a real run artifact                                                                          |
| **3**   | Evaluation runs E1–E5 (§8), figure generation scripts, repo tidy-up                                                                 | All figures regenerable from stored artifacts                                                                  |




## Suggested split

Roughly one engineer per tier through week 2 — Tier 1 is the larger and more time-critical piece (it blocks everything downstream), Tier 2 is smaller but has more external-tool integration risk (fpocket, APBS). Converge on the dashboard and evaluation in week 3. **Confirm in meeting.**

## Flag to Laura immediately if

- The reference corpus is delayed past week 1 (blocks all percentile scoring)
- fpocket returns ambiguous cavities for either demo target
- The Spearman diagnostic (§6.4) comes back high
- E3 (matched vs mismatched target) returns a null result
- APBS integration threatens the week-2 milestone

