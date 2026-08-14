# AptaRank — Refinements Specification

**Companion to:** `IMPLEMENTATION_SPEC.md` (the original build spec). This document  
supersedes the target-input and Tier 2 sections of that spec where they conflict;  
everything else in the original spec still holds.

---

# 0. Priority bands

Band 1 is cheap and fixes things that damage trust in a live demo or review; Band 2 is the core scientific redesign and the bulk of the work.


| Band                                       | Scope                                                                 | Why this order                                                                                                      |
| ------------------------------------------ | --------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------- |
| **1 — quick trust/credibility fixes** (§1) | Provenance/eligibility bug, filename mismatch, copy/labelling changes | Cheap, and each is embarrassing if a reviewer or demo viewer hits it. Bank them before the big redesign.            |
| **2 — core redesign** (§2–§5)              | Tier 2 binding modes + target-input redesign, built together          | The real work. This is what answers the supervisor's review question and makes the paper's contribution defensible. |


---



# BAND 1



# 1. Quick fixes (copy, labelling, provenance)



## 1.1 Terminology — two audiences, never mixed

We now use different terms user-facing vs paper/spec/artifact. **Within each surface,
use one term consistently** — the current bug is that the sidebar and results page
disagree with each other.


| Concept | User-facing (UI copy)            | Paper / spec / JSON artifact fields      |
| ------- | -------------------------------- | ---------------------------------------- |
| Tier 1  | **aptamer-likeness**             | **intrinsic structural quality**         |
| Tier 2  | **aptamer-target compatibility** | **control-relative geometric agreement** |


- Sweep **all** UI copy to the user-facing terms. The sidebar blurb, the results
page, tooltips, panel headers — all must agree.
- Keep JSON artifact field names in paper-side vocabulary (`tier1_*`,
`intrinsic_*`, `tier2_*`, `geometric_agreement_*`). The artifact is the
reproducibility record that backs the paper.
- New sidebar blurb: *"Tier 1 ranks candidates on aptamer-likeness. Tier 2 annotates them with aptamer-target compatibility by performing a geometric agreement check against a target. Tier 2 is not a binding prediction, and it does not change the Tier 1 ranking."*



## 1.2 Reference library label

- Change the description to: **"experimentally validated RNA aptamers, or your own
aptamer reference dataset."**



## 1.3 Protein target section — rename and re-describe

- Rename **"active site"** → **"target binding site"** everywhere. ("Active site" is an enzyme-only term; most aptamer targets are not enzymes — see §3.1.)
- New section description: *"AptaRank annotates each candidate with how well its shape agrees with measurements of the protein target's binding site. Specify a target and binding mode, and optionally target binding-site residues to choose the exact site. AptaRank does not predict binding."*



## 1.4 Analysis-thoroughness copy

- Under "How thorough should the analysis be?", change **"Needed for statistics quoted
in a paper"** → **"Best for statistics quoted in a paper."** (The rigorous setting is
a recommendation, not a hard requirement.)



## 1.5 Provenance / publication-eligibility bug (the important one)

**Symptom:** the New Analysis page displayed `publication_eligible` while the inputs
(reference library with no provenance record) actually make the run **not** eligible.
The Results page then correctly showed `publication_eligible: false` with a warning.
The pre-run promise contradicts the post-run verdict.

**Fix:** compute eligibility from the actual inputs **at or before submission**, and
show the same verdict on the New Analysis page that the Results page will show. The
rule is unchanged — a run is publication-eligible only if every input (candidate set,
reference corpus, target bundle) carries a provenance record and none is marked
synthetic/dev. What changes is *when* the check runs: it must run up front, not only
after completion.

- If any input lacks provenance at submission time, surface the "development run —
not a result" state **on the New Analysis page**, before the user runs anything.
- Never show `publication_eligible` and later retract it.



## 1.6 Provenance filename mismatch

**Symptom:** the Run Provenance tab shows input/corpus filenames that don't match what
the user uploaded.

**Cause:** almost certainly the backend is displaying internal cached/hashed filenames
rather than the user's original filenames.

**Fix:** store the original upload filename alongside the internal one. Show the
original to the user; keep the hash internally for integrity checks. The provenance
record should show both: `uploaded as "my_candidates.csv" (sha256:…)`.

---



# BAND 2



# 2. Why Tier 2 is being redesigned



## 2.1 The problem a reviewer raised

The current Tier 2 assumes a single binding mechanism: that an aptamer binds by
inserting an unpaired **loop** into a protein **pocket**, so compatibility ≈ does the
loop size fit the cavity size. A senior reviewer asked, in effect, *why would an
unpaired loop region have binding properties at all?* — and the honest answer is that
loop-into-pocket is only **one** of several ways aptamers bind proteins, and not the
most common one for non-enzyme targets.

The literature is clear that protein-binding aptamers recognise targets through
several distinct geometries: presenting a flexible loop/single-stranded region;
lying a helical stem along a surface groove; or covering a flat-ish surface patch
via shape and charge complementarity (plus mechanisms we are explicitly not
modelling, like induced-fit reshaping). The current tool bakes in the first and
presents it as universal. That is the weakness we are fixing.

## 2.2 The fix, in one sentence

**Keep the coarse-geometric-compatibility concept, but make the *type* of geometric
comparison depend on a binding mode that the expert user selects.**

The tool does **not** predict the binding mode — that would be a research project. The
expert biologist asserts the mode (they are well-placed to judge what kind of site
they are targeting), and the tool assesses geometric compatibility *appropriate to
that mode*. This division of labour is the whole point: expert supplies the biological
judgment, tool supplies the fast, reproducible, coarse geometry.

This **strengthens** the "coarse compatibility, not docking" framing rather than
weakening it — we are now explicitly modelling different physics coarsely, instead of
pretending one coarse model is universal.

## 2.3 What stays exactly the same

Everything around the comparison is unchanged. The blast radius is contained to "what
geometry gets compared."

- Tier 1 (unchanged)
- Shuffled-control banding into strong / moderate / weak (unchanged)
- "Tier 2 annotates, never reorders" ranking rule (unchanged)
- Dashboard structure, provenance, "not a binding prediction" caveat (unchanged)
- The Gaussian agreement + control-percentile machinery (unchanged — only its inputs change)



## 2.4 Scope for the paper: two modes

Build **Mode A (pocket)** and **Mode C (surface patch)** for the paper. These cover the
two headline demo targets (NDM-1 → A, IGFBP3 → C) and prove the core claim: *the
compatibility check adapts to the expert-specified binding mode.* **Mode B (groove) and
other mechanisms are described as future work, not built.** Two well-executed modes
prove the claim as well as five and fit the deadline.

**Headline demo targets (locked):**

- **Mode A — NDM-1 (PDB** `3SPU`**).** Enzyme with a real catalytic pocket; keeps the
AMR/thesis framing.
- **Mode C — IGFBP3 (PDB** `7WRQ`**, IGFBP3 chain).** Binding/transport protein with an
interface rather than a catalytic pocket; also the target for which fine-tuned
RNAGenesis generation results exist.

**Generation conditions for the demo (three-way):**

- **FM-generated:** RNAGenesis fine-tuned on aptamer data, generating candidates against
IGFBP3 (the target for which fine-tuned generation results currently exist).
- **Validated aptamers:** real aptamer–target pairs, as a positive control (see §9.2).
- **Random RNA:** length- and composition-matched, as the negative control.

> **Note — evaluation targets are not build inputs.** The pipeline only needs to *accept*
> a target config (§3.2) and score it; which validated aptamers are then used as
> positive controls is an evaluation-time decision. The positive-control list finalises
> with the reference corpus (§9.2) and does **not** block implementation.

---



# 3. Target input redesign



## 3.1 Concept

Replace the "upload a prepared `.json` target file" flow (which required the user to run
`aptarank target build` on Linux/GitHub) with a simple identifier-based flow. The user
provides an ID and a chain; **the backend does all structure preparation server-side**.
The heavy biology-literate step moves off the user and into the pipeline, where the
original spec's §5.3–5.4 logic already lives.

The `DEMO / fabricated cavity` option and its provenance guarding (SYNTHETIC warnings,
`publication_eligible: false`) are well-designed and **stay unchanged**.

## 3.2 What the user specifies

Two sources only (UniProt auto-lookup was considered but dropped to save time):


| Source           | ID                         | What it is                           | Key handling                                                                          |
| ---------------- | -------------------------- | ------------------------------------ | ------------------------------------------------------------------------------------- |
| **PDB**          | 4-char PDB ID, e.g. `7WRQ` | Experimental structure               | May be a **complex** — several chains in one file. Handled by chain selection (§3.4). |
| **AlphaFold DB** | UniProt ID, e.g. `P17936`  | **Predicted** single-chain structure | Always chain `A`. Flag as predicted (§3.5).                                           |


User-facing config:

```yaml
target_name: IGFBP3
source: pdb                 # pdb | alphafold
id: 7WRQ                    # PDB ID for source=pdb; UniProt ID for source=alphafold
chain: B                    # which chain is the target protein
binding_mode: surface       # pocket | surface   (see §5)
target_site_residues: [42, 87, 119, 121]   # optional; integers, chain-B numbering
```

- `target_site_residues` are **integers** (residue sequence numbers), keyed to the
chain already named in `chain:`. Show the user an example and explicitly state **not**
to write `K42` — just `42`. (The amino-acid letter is redundant; the number + chain
select the residue. See §3.6.)
- If `target_site_residues` is empty, fall back to automatic pocket selection (§5.2,
step 4 fallback) and mark `pocket_selection: "automatic"` in the artifact so the UI
can caveat it.



## 3.3 Backend flow

1. **Fetch structure.**
  - PDB: `https://files.rcsb.org/download/{ID}.pdb`
  - AlphaFold: AlphaFold DB API, `https://alphafold.ebi.ac.uk/files/AF-{UNIPROT_ID}-F1-model_v4.pdb`
  (confirm current version suffix at build time)
  - Cache by `{source}_{id}`.
2. **Select chain**, keep only that chain (§3.4).
3. **Strip waters** (residue `HOH`); **keep other hetero groups by default**
  (`strip_hetatm: false`) — for demo targets these include a functional metal ion that
   shapes the cavity. (Same as original spec §5.3.)
4. **Run the shared measurement layer** (§4).
5. **Run the mode-specific comparison** (§5).
6. Produce the target bundle **server-side**. The user never sees or uploads a JSON.



## 3.4 Complex handling (chain removal)

**This is straightforward — it is a filter, not structural surgery.** A PDB complex is
multiple chains in one file; every atom line carries a chain ID. "Keep chain B, drop the
rest" is a one-pass filter on that label (Biopython: iterate chains, detach those that
are not the requested one). Mechanically identical to water-stripping, keyed on chain
instead of residue name.

**Critical ordering — remove the partner, but let it point first.** For an interface
target like IGFBP3 in surface mode, the partner chain (IGF1) sits *on* the interface we
want to measure. So:

1. Use the user's `target_site_residues` (which mark the interface) to locate the
  region of interest **while the full complex is still loaded**.
2. **Then** strip the partner chains, exposing the interface surface.
3. Run the measurement layer on the isolated target chain.

The partner's location informs *where* to look; the partner itself is removed so the
surface is exposed for measurement. This is why residue specification and chain
selection work together.

**Honest limitation to display (not a blocker):** if the target chain's functional
surface is shaped *by* the partner and partially reshapes once the partner is removed,
the measured geometry is of the unbound state, which can differ from the bound state.
This is the same induced-fit caveat that applies to the whole coarse-geometry approach —
we surface it, we do not try to solve it. Add a short note in the target panel when
`source=pdb` and the selected chain was part of a multi-chain file.

## 3.5 AlphaFold handling

- **Allowed**, with a **mode-aware warning:**
  - `source=alphafold` **and** `binding_mode=surface`: warn prominently — *"Predicted
  structures may not display an interface that only forms when a binding partner is
  present. Surface-mode results on a predicted structure should be treated with
  extra caution."*
  - `source=alphafold` **and** `binding_mode=pocket`: milder note — catalytic/pocket
  geometry is often well predicted, but it is still a model, not an experiment.
- **Always** set `structure_source: "predicted"` in the run artifact (vs `"experimental"`
for PDB). This is visible downstream and consistent with the existing
provenance/eligibility machinery.



## 3.6 Residue numbering notes (for the parser)

- Store residues as **integers**, compared on `(chain, residue_number)`. Biopython
exposes these as `residue.get_id()[1]` (number) and `residue.parent.id` (chain).
- **Author numbering ≠ sequential position.** PDB residue numbers are chosen by the
depositor and often don't start at 1 (e.g. a structure resolving residues 28–291 numbers
its first residue 28). `target_site_residues: [42]` means "the residue **labelled** 42,"
not "the 42nd residue in the file." Biopython preserves author numbering, so published
residue tables (PDBsum, UniProt) line up — do not re-index.
- **Insertion codes** (e.g. `100`, `100A`) are rare and almost certainly absent in demo
targets. Assert none are present rather than silently mismatching; Biopython's residue
ID is a tuple `(hetflag, resseq, icode)`.

---



# 4. Shared measurement layer (mode-independent)

All modes draw on the **same two measurements**, computed once per run. Modes differ
only in *which* element they compare and *what shape logic* they apply (§5). This is
what keeps the build small — one measurement layer, three (two built) comparison
functions over it.

## 4.1 Protein side — cavity/surface geometry (once per target)

Run fpocket on the prepared single-chain structure (original spec §5.4), then compute
descriptors for each detected pocket:


| Descriptor                               | Source                                                                                                           | Used by                                             |
| ---------------------------------------- | ---------------------------------------------------------------------------------------------------------------- | --------------------------------------------------- |
| `volume_A3`                              | fpocket `_info.txt`                                                                                              | Mode A                                              |
| `d_pocket_A` (characteristic dimension)  | SVD of alpha-sphere coords (spec §5.4)                                                                           | Mode A                                              |
| `d_equiv_A` (equivalent-sphere diameter) | from volume (spec §5.4)                                                                                          | cross-check                                         |
| `elongation` (NEW)                       | ratio of largest to smallest principal-axis extent of the alpha-sphere cloud (from the same SVD — `sv[0]/sv[2]`) | distinguishes channel vs round pocket vs flat patch |
| `planarity_A` (NEW)                      | smallest principal-axis extent (`2*sv[2]/sqrt(k)`); small = flat/planar patch                                    | Mode C                                              |
| `patch_area_A2` (NEW)                    | solvent-accessible surface area of the selected binding-site residues, via freeSASA                              | Mode C                                              |


The two SVD-derived descriptors (`elongation`, `planarity_A`) are essentially free —
they come from the **same** singular values already computed for `d_pocket_A`. Only
`patch_area_A2` needs a new tool call (freeSASA), and only for surface mode.

**New tool:** `freesasa` (`pip install freesasa`) — computes solvent-accessible surface
area per residue. Confirm it installs cleanly at build kickoff; it is a small, standard,
pure-computation library with no external dependencies.

## 4.2 RNA side — element geometry (per candidate)

All already computed in Tier 1 via forgi (original spec §4.6–4.7). No new RNA-side
computation, just exposing existing values to the comparison functions:


| Descriptor                      | Source                                                                                                                      | Used by         |
| ------------------------------- | --------------------------------------------------------------------------------------------------------------------------- | --------------- |
| `loop_nt_median`, `loop_nt_p90` | ensemble loop-size sampling (spec §4.7)                                                                                     | Mode A          |
| `longest_stem_bp`               | forgi stem lengths (spec §4.6)                                                                                              | Mode B (future) |
| `footprint_nt` (NEW, cheap)     | total sequence length as a proxy for overall molecular footprint; refine later with radius-of-gyration proxy if time allows | Mode C          |




## 4.3 Unit conversions (documented constants)

Keep every physical constant in config with a comment on what it assumes. These are
coarse by design.

```yaml
geometry:
  a_per_nt_ss: 6.0        # extended single-stranded rise per nucleotide, Å (Mode A loops)
  a_per_bp_helix: 2.8     # A-form helix rise per base pair, Å (Mode B, future)
  footprint_scale: 1.0    # Mode C footprint scaling, Å per nt (coarse)
```

---



# 5. Binding modes

Each mode is a **comparison function** over the shared measurement layer. Same output
shape for all modes: a geometric agreement score in (0,1], then banded by
shuffled-control percentiles (unchanged machinery). The mode is recorded in the
artifact as `binding_mode`, and the explanation text and dashboard adapt their wording
to the mode.

## 5.1 Mode A — Pocket engagement (BUILD)

- **Premise (shown in UI/paper):** the aptamer presents a flexible loop / unpaired region
that engages a concave pocket on the target.
- **Protein measurement:** pocket characteristic dimension `d_pocket_A` at the selected
residues.
- **RNA measurement:** largest accessible loop extent = `a_per_nt_ss * loop_nt_median`.
- **Comparison:** Gaussian agreement on the size difference (original spec §5.6):
`exp(-(d_apt - d_pocket)^2 / (2 * sigma^2))`.
- **Caveat surfaced:** if the pocket is oddly shaped (`d_pocket_A` and `d_equiv_A`
disagree by > ~2×, i.e. high `elongation`), warn that the geometric comparison assumes
a roughly convex pocket. (This warning already exists — keep it.)
- **Demo target:** NDM-1 (`3SPU`), an enzyme with a real catalytic pocket.



## 5.2 Mode C — Surface-patch recognition (BUILD)

- **Premise (shown in UI/paper):** the aptamer covers a flat-ish surface patch via shape
and charge complementarity, rather than plugging a pocket.
- **Protein measurement:** at the selected residues — `patch_area_A2` (freeSASA),
`planarity_A` (flat vs curved), and the electrostatic character of the patch
(mean potential over the patch, from the APBS grid already planned in spec §5.7).
- **RNA measurement:** `footprint_nt` → footprint area proxy; plus RNA carries a uniform
negative charge (no per-candidate charge computation needed — it's a property of the
backbone).
- **Comparison — a small composite (still coarse, still one score):**
  1. **Size/coverage:** is the RNA footprint large enough to cover the patch area?
    Gaussian agreement on (RNA footprint area − patch area), or a coverage ratio.
  2. **Charge complementarity:** RNA is negative; a positively charged patch is
    favourable, a negative one unfavourable. Map mean patch potential → a
     complementarity factor in (0,1].
  3. Combine (1) and (2) into a single geometric-agreement score (documented weights,
    equal by default), then band as usual.
- **Caveat surfaced:** planarity is coarse; a highly curved "patch" is flagged, and
charge complementarity is a directional signal, not an affinity estimate.
- **Demo target:** IGFBP3 (`7WRQ`, chain per §3.4), a binding/transport protein with an
interface rather than a catalytic pocket.



## 5.3 Mode B — Groove/helical presentation (FUTURE WORK)

Described in the paper as a supported future mode: the aptamer lies a helical stem
(`longest_stem_bp` → length via `a_per_bp_helix`) along an elongated surface channel
(high-`elongation` pocket). Comparison would check helix length vs channel length and
channel width vs helix diameter (~20 Å). **Not implemented for this paper.**

## 5.4 Mode dispatch

```python
def tier2_score(candidate, target, mode, config):
    if mode == "pocket":
        return score_pocket(candidate, target, config)   # §5.1
    elif mode == "surface":
        return score_surface(candidate, target, config)  # §5.2
    else:
        raise ValueError(f"Unsupported binding mode: {mode}")  # 'groove' = future
```

The banding (§5, shuffled-control percentiles), the "annotates never reorders" merge,
and the dashboard all consume `tier2_score` / `tier2_band` **identically regardless of
mode**. Only the wording adapts.

---



# BAND 3



# 6. Polish / defer if tight

- **Electrostatics refinement (Mode C):** if the APBS integration for the patch-potential
signal threatens the Band 2 timeline, ship Mode C with the size/coverage signal only
and add charge complementarity as a fast-follow. Report electrostatics as "computed
where available."
- **RNA footprint refinement:** replace the length-proxy `footprint_nt` with a
radius-of-gyration-style proxy derived from the forgi element graph, if time allows.
- `footprint`**/**`planarity` **validation:** sanity-check the new descriptors against the
two demo targets before relying on them in evaluation figures.

---



# 7. What is explicitly unchanged (contain the blast radius)

So you know what **not** to touch:

- **Tier 1** — all of it (ingest, features, ensemble sampling, shuffled controls,
rank-aggregation composite, ranking).
- **Banding** — strong/moderate/weak from shuffled-control percentiles, recomputed per
target.
- **Ranking rule** — Tier 1 orders; Tier 2 annotates; Tier 2 never reorders;
"not_evaluated" ≠ "weak".
- **Spearman diagnostic** (spec §6.4).
- **Dashboard layout**, run-artifact schema (add the new fields, don't restructure),
explanation-panel mechanism (rule-based; add mode-aware rules), structure diagrams.
- **DEMO/synthetic option** and all provenance/eligibility machinery (except the two
bug fixes in §1.6–1.7).
- **Evaluation design** (spec §8), except that the generation-condition story now runs
on IGFBP3 (FM-generated) + validated aptamers + random RNA (§2.4), and Mode C enables
IGFBP3 as the surface-mode demonstration. NDM-1 remains the Mode A headline; KPC-2 is no
longer part of the generation axis (may survive only as an optional target-swap
illustration if time allows).

---



# 8. Still open for the team


| #   | Item                                                   | Note                                                                                                           |
| --- | ------------------------------------------------------ | -------------------------------------------------------------------------------------------------------------- |
| 1   | freeSASA install + API check                           | Confirm at kickoff; small standard library, low risk                                                           |
| 2   | AlphaFold DB URL version suffix                        | Confirm current `model_v4` (or later) at build time                                                            |
| 3   | Mode C composite weights                               | Default equal; revisit if one signal dominates                                                                 |
| 4   | Whether Mode C charge signal ships in Band 2 or Band 3 | Capacity call for Ali & Mohammad                                                                               |
| 5   | freeSASA patch-area definition                         | Whole selected-residue set vs pocket-overlapping subset — decide with Laura once first numbers are visible     |
| 6   | Mode A (pocket) positive control                       | Unmodified, binary, single-chain, < ~600 Å² interface — Laura checking corpus first, else narrow search (§9.2) |
| 7   | Evaluation positive-control pair list                  | Finalises with the corpus (next day); not blocking — build the scoring path, Laura supplies pairs (§9.2)       |


---



# 9. Demo target inputs and evaluation controls

This section reports the exact user inputs (per §3.2) for the headline demo targets, and specifies how the evaluation positive controls are handled. **Nothing here is a build input** — the pipeline is built to the target-input contract (§3.2) and these
configs are supplied at run time. Build to accept them; do not hardcode any residue list.

> **⚠️ Residue numbers below are placeholders pending Laura's confirmation.** Do not treat
> any `[TBD]` as final, and do not substitute guessed positions — an incorrect residue
> number silently points the tool at the wrong site. Laura supplies the confirmed
> integers (chain-specific author numbering, per §3.6) from each structure before they
> are used in runs.



## 9.1 Headline demo targets



### IGFBP3 — Mode C (surface), FM-generated demo

```yaml
target_name: IGFBP3
source: pdb
id: 7WRQ
chain: [TBD — the IGFBP3 chain in the 7WRQ complex; confirm from the entry]
binding_mode: surface
target_site_residues: [TBD — IGFBP3 residues at the IGF1 interface;
                       from PDBsum interface table or computed as residues
                       within ~4–5 Å of the IGF1 chain (§3.4)]
```

Notes: 7WRQ is a **ternary complex** (IGFBP3 + IGF1 + ALS). Chain selection isolates
IGFBP3; the interface residues are those contacting IGF1, used to locate the patch
**before** the partner chains are stripped (§3.4).

### NDM-1 — Mode A (pocket)

```yaml
target_name: NDM-1
source: pdb
id: 3SPU
chain: A
binding_mode: pocket
target_site_residues: [TBD — zinc-coordinating active-site residues; Laura to supply]
```

Notes: keep `strip_hetatm: false` so the functional zinc ion is retained (it shapes the
cavity, §3.3).

## 9.2 Evaluation positive controls (finalise with the corpus — not blocking)

Positive controls are validated real aptamer–target pairs with solved structures, scored
through the tool to confirm the geometric-agreement score behaves sensibly on known
binders. They are **evaluation fixtures, not headline demos**, so the bar on "unmodified"
and "clean single chain" is relaxed — a modified aptamer is acceptable here because it is
a test point, not a showcase.

- **Source:** drawn from the reference corpus's labelled aptamer–target pairs that have
solved PDB structures. **This list finalises when the corpus does (expected next day)
and does not block implementation** — the engineers build the scoring path; Laura
supplies the specific pairs.
- **Candidate Mode C control — Toggle-25t / human α-thrombin (PDB** `3DD2`**).** A 1.9 Å
binary complex with a quantified interface (buried area ≈ 1,193 Å², shape
complementarity ≈ 0.75) and mutagenesis-confirmed contact residues — an excellent
surface-mode control. Caveat recorded for honesty: the aptamer carries 2′-fluoro
pyrimidines, so it is a modified molecule scored by an unmodified-RNA yardstick;
acceptable for a control, and noted as such.
- **Mode A (pocket) control — placeholder, pending search/corpus.** Still to be selected:
an **unmodified** RNA aptamer in a **binary, single-chain** complex, ≤2.5 Å, with a
**small buried interface (< ~600 Å²)** — the signature of genuine pocket binding
regardless of the authors' wording — and ideally mutagenesis-confirmed interface
residues. Laura is checking whether a corpus pair already meets this filter before
searching externally.



## 9.3 What Laura delivers

For each headline target: confirmed `id`, `chain`, and the integer `target_site_residues`
list in the correct chain's author numbering. For evaluation: the finalised positive-
control pair list (with structures) once the corpus is locked. Engineers treat the
pipeline as ready to accept these; **no residue list or target is hardcoded** — all live
in per-target config (§3.2), supplied at run time.