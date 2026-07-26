"""Dinucleotide-preserving shuffled controls (spec §4.8).

A high structural score could simply reflect letter composition: GC-rich
strings fold tightly regardless of how the letters are arranged. Comparing each
candidate against shuffles of *itself* that preserve dinucleotide composition
shows whether the score reflects arrangement rather than composition.

k=2 preserves every adjacent letter-pair frequency, which is stricter and more
honest than preserving single-letter frequencies alone.

Statistics (deviating from the spec's `percentile >= 0.95` rule, deliberately):
beating 19 of 20 controls is a descriptive 95% win rate, not a 0.05-level test
— with the candidate and M controls pooled, the two best positions already
account for 2/(M+1) of the pooled ranks. We therefore report both the
descriptive percentile and the exact Monte-Carlo p-value

    p = (1 + #{control >= real}) / (M + 1)

and base `shuffle_pass` on p <= alpha. With M=20 and alpha=0.05 that requires
strictly beating all 20 controls (p = 1/21 = 0.0476). Publication runs should
use M=99.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import ushuffle

from ..errors import AptaRankError


class ShuffleCorruptionError(AptaRankError):
    """ushuffle produced something that is not a permutation of its input."""


@dataclass(frozen=True)
class ShuffleOutcome:
    """Result of comparing one candidate against its own shuffled controls."""

    n_shuffles: int
    wins: int
    ties: int
    percentile: float
    p_value: float
    passed: bool
    real_score: float
    median_score: float
    margin: float

    def to_dict(self) -> dict[str, float | int | bool]:
        return {
            "pass": self.passed,
            "percentile": self.percentile,
            "p_value": self.p_value,
            "n_shuffles": self.n_shuffles,
            "wins": self.wins,
            "ties": self.ties,
            "real_score": self.real_score,
            "median_score": self.median_score,
            "margin": self.margin,
        }


def generate_shuffles(sequence: str, n: int, k: int = 2, seed: int | None = None) -> list[str]:
    """`n` dinucleotide-preserving shuffles of `sequence`.

    Seeded per candidate (not once per process) so that results do not depend
    on how a process pool happens to distribute work.

    Two safety measures around the ushuffle C extension:

    * `ushuffle.Shuffler` stores a raw `char*` to the bytes it is constructed
      with and does not hold a reference to them. Passing `sequence.encode()`
      inline lets that temporary be collected, after which every subsequent
      `shuffle()` returns bytes from reused memory — silently corrupted
      controls, not an exception. `payload` below keeps it alive.
    * Every shuffle is checked to be a permutation of the input. A corrupted
      control that happened to stay in the ASCII range would otherwise sail
      through and quietly invalidate the §4.8 statistic.
    """
    if n <= 0:
        return []
    if seed is not None:
        ushuffle.set_seed(int(seed) % (2**31 - 1))

    payload = sequence.encode()  # must outlive the Shuffler — see above
    shuffler = ushuffle.Shuffler(payload, k)
    expected = kmer_counts(sequence, k)

    out = []
    for _ in range(n):
        shuffled = shuffler.shuffle().decode("ascii", errors="replace")
        # Checking the k-mer multiset, not just the letter multiset: preserved
        # dinucleotide composition is exactly the property §4.8 claims, so it
        # is the property worth asserting rather than assuming.
        if len(shuffled) != len(sequence) or kmer_counts(shuffled, k) != expected:
            raise ShuffleCorruptionError(
                f"ushuffle did not preserve {k}-let composition "
                f"({sequence!r} -> {shuffled!r})"
            )
        out.append(shuffled)
    return out


def kmer_counts(sequence: str, k: int) -> dict[str, int]:
    """Multiset of overlapping k-mers — the quantity a k-let shuffle preserves."""
    counts: dict[str, int] = {}
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i : i + k]
        counts[kmer] = counts.get(kmer, 0) + 1
    return counts


def evaluate(real_score: float, shuffle_scores: Sequence[float], alpha: float) -> ShuffleOutcome:
    """Compare a candidate's structural sub-score against its controls."""
    scores = [float(s) for s in shuffle_scores]
    m = len(scores)
    if m == 0:
        return ShuffleOutcome(
            n_shuffles=0, wins=0, ties=0, percentile=float("nan"), p_value=float("nan"),
            passed=False, real_score=real_score, median_score=float("nan"),
            margin=float("nan"),
        )

    wins = sum(1 for s in scores if real_score > s)
    ties = sum(1 for s in scores if real_score == s)
    percentile = (wins + 0.5 * ties) / m
    p_value = (1 + sum(1 for s in scores if s >= real_score)) / (m + 1)
    median = _median(scores)

    return ShuffleOutcome(
        n_shuffles=m,
        wins=wins,
        ties=ties,
        percentile=percentile,
        p_value=p_value,
        passed=p_value <= alpha,
        real_score=real_score,
        median_score=median,
        margin=real_score - median,
    )


def _median(values: Sequence[float]) -> float:
    ordered = sorted(values)
    n = len(ordered)
    mid = n // 2
    return ordered[mid] if n % 2 else 0.5 * (ordered[mid - 1] + ordered[mid])
