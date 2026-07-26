"""Generate synthetic development data.

THIS IS NOT SCIENTIFIC DATA. It exists so the pipeline can be built and tested
before the curated validated-aptamer corpus is delivered. Runs scored against
it are stamped `publication_eligible: false` and must never back a paper claim.

    python scripts/make_dev_data.py

writes:
    data/corpus/dev_placeholder_corpus.csv   synthetic "reference" aptamers
    data/demo_candidates.csv                 mixed-quality candidate batch
"""

from __future__ import annotations

import csv
import random
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
COMPLEMENT = {"A": "U", "U": "A", "G": "C", "C": "G"}
SEED = 20260726


def stem_loop(rng: random.Random, stem_bp: int, loop_nt: int, tail_nt: int = 0) -> str:
    """A sequence that folds into one hairpin: stem, loop, reverse-complement."""
    stem = "".join(rng.choices("ACGU", weights=[1, 1.4, 1.4, 1], k=stem_bp))
    loop = "".join(rng.choices("ACGU", k=loop_nt))
    closing = "".join(COMPLEMENT[c] for c in reversed(stem))
    tail = "".join(rng.choices("ACGU", k=tail_nt))
    return stem + loop + closing + tail


def degrade(rng: random.Random, seq: str, n_mutations: int) -> str:
    """Break base pairs so the sequence folds less cleanly."""
    chars = list(seq)
    for pos in rng.sample(range(len(chars)), min(n_mutations, len(chars))):
        chars[pos] = rng.choice([c for c in "ACGU" if c != chars[pos]])
    return "".join(chars)


def random_rna(rng: random.Random, length: int) -> str:
    return "".join(rng.choices("ACGU", k=length))


def make_corpus(rng: random.Random, n: int) -> list[dict[str, str]]:
    rows = []
    targets = [("NDM-1", "3SPU"), ("thrombin", "4DII"), ("VEGF", "1FLT"), ("lysozyme", "1DPX")]
    for i in range(1, n + 1):
        stem_bp = rng.randint(5, 12)
        loop_nt = rng.randint(4, 12)
        tail = rng.randint(0, 14)
        seq = stem_loop(rng, stem_bp, loop_nt, tail)
        if rng.random() < 0.35:  # some corpus entries are two-hairpin
            seq = seq + stem_loop(rng, rng.randint(4, 8), rng.randint(4, 9))
        seq = seq[:100]
        if len(seq) < 20:
            continue
        name, pdb = rng.choice(targets)
        rows.append(
            {
                "id": f"ref{i:05d}",
                "sequence": seq,
                "target_name": name,
                "target_pdb_id": pdb,
                "source_reference": "SYNTHETIC-PLACEHOLDER",
            }
        )
    return rows


def make_candidates(rng: random.Random, n: int) -> list[dict[str, str]]:
    rows = []
    for i in range(1, n + 1):
        draw = rng.random()
        if draw < 0.4:                      # well-formed, aptamer-like
            seq = stem_loop(rng, rng.randint(6, 12), rng.randint(5, 11), rng.randint(0, 10))
        elif draw < 0.7:                    # partially broken
            seq = degrade(rng, stem_loop(rng, rng.randint(6, 11), rng.randint(5, 12), 6),
                          rng.randint(2, 6))
        else:                               # unstructured
            seq = random_rna(rng, rng.randint(28, 70))
        seq = seq[:100]
        if len(seq) < 20:
            continue
        rows.append({"id": f"cand{i:04d}", "sequence": seq})
    return rows


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {len(rows):>5} rows  {path.relative_to(REPO_ROOT)}")


def main() -> None:
    rng = random.Random(SEED)
    write_csv(REPO_ROOT / "data" / "corpus" / "dev_placeholder_corpus.csv", make_corpus(rng, 800))
    write_csv(REPO_ROOT / "data" / "demo_candidates.csv", make_candidates(rng, 200))


if __name__ == "__main__":
    main()
