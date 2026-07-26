"""Rule-based explanations (spec §7.3).

Rule-based, not an LLM call: deterministic, auditable, no inference dependency,
and every sentence is guaranteed to be supported by a number that is already
stored in the artifact.

Two invariants this module must never break:

1. Every number in the text is read from the candidate record — never
   recomputed, never rounded differently than the table shows.
2. Cautionary rules are first-class. A candidate that fails its shuffled
   control says so, in its own explanation, before any praise is added.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Mapping

from .. import TIER2_CAVEAT

CAUTION = "caution"
POSITIVE = "positive"
NEUTRAL = "neutral"


@dataclass(frozen=True)
class Rule:
    name: str
    kind: str
    priority: int                      # lower fires first
    condition: Callable[[dict[str, Any]], bool]
    template: str
    chip: str | None = None


def _crit(rec: Mapping[str, Any], name: str, field: str = "score") -> float | None:
    entry = (rec.get("criteria") or {}).get(name) or {}
    value = entry.get(field)
    return None if value is None else float(value)


def _flatten(rec: Mapping[str, Any]) -> dict[str, Any]:
    """Flat view of the stored record, for `str.format` and rule conditions."""
    criteria = rec.get("criteria") or {}
    elements = rec.get("elements") or {}
    shuffle = rec.get("shuffle") or {}
    tier2 = rec.get("tier2") or {}
    flat: dict[str, Any] = {
        "candidate_id": rec.get("candidate_id"),
        "rank": rec.get("rank"),
        "length": rec.get("length"),
        "duplicate_count": rec.get("duplicate_count") or 1,
        "tier1_score": rec.get("tier1_score"),
        "batch_rank_fraction": rec.get("batch_rank_fraction"),
        "band": tier2.get("band", "not_evaluated"),
        "tier2_status": tier2.get("status", "not_evaluated"),
    }
    for name, entry in criteria.items():
        flat[name] = entry.get("value")
        score = entry.get("score")
        flat[f"{name}_score"] = score
        # `_pct` reads as "how far up the corpus this sits"; `_top_pct` is its
        # complement, for sentences phrased as "in the top N%".
        flat[f"{name}_pct"] = None if score is None else 100.0 * float(score)
        flat[f"{name}_top_pct"] = None if score is None else 100.0 * (1.0 - float(score))
    flat.update({k: v for k, v in elements.items()})
    n_hairpins = elements.get("n_hairpins")
    if n_hairpins is not None:
        flat["hairpin_phrase"] = f"{n_hairpins} hairpin loop" + ("s" if n_hairpins != 1 else "")
    flat.update({f"shuffle_{k}": v for k, v in shuffle.items()})
    flat.update({f"tier2_{k}": v for k, v in tier2.items()})
    return flat


def _has(rec: dict[str, Any], key: str) -> bool:
    return rec.get(key) is not None


RULES: tuple[Rule, ...] = (
    # -- cautions first: they are the ones a biologist must not miss ------
    Rule(
        "shuffle_fail", CAUTION, 10,
        lambda r: r.get("shuffle_pass") is False,
        "Does not outscore its own dinucleotide-shuffled controls "
        "(beat {shuffle_wins} of {shuffle_n_shuffles}, p = {shuffle_p_value:.3f}), so this "
        "structural score may reflect nucleotide composition rather than arrangement.",
        chip="shuffle control failed",
    ),
    Rule(
        "poorly_defined_fold", CAUTION, 20,
        lambda r: _has(r, "ensemble_defect_score") and r["ensemble_defect_score"] < 0.3,
        "Does not settle on one shape — ensemble defect {ensemble_defect:.3f} is in the "
        "bottom {ensemble_defect_pct:.0f}% of the reference corpus.",
        chip="unstable fold",
    ),
    Rule(
        "barely_pairs", CAUTION, 30,
        lambda r: _has(r, "stem_fraction") and r["stem_fraction"] < 0.2,
        "Only {stem_fraction:.0%} of positions are paired, so the candidate has little "
        "defined structure to present to a target.",
        chip="little structure",
    ),
    Rule(
        "atypical_composition", CAUTION, 40,
        lambda r: _has(r, "gc_fraction_score") and r["gc_fraction_score"] < 0.25,
        "GC content {gc_fraction:.0%} sits in the tail of the corpus distribution "
        "(typicality {gc_fraction_score:.2f}).",
        chip="atypical GC",
    ),
    Rule(
        "uncertain_loop_size", CAUTION, 50,
        lambda r: _has(r, "loop_nt_iqr") and _has(r, "loop_nt_median")
        and r["loop_nt_median"] > 0 and r["loop_nt_iqr"] >= 0.5 * r["loop_nt_median"],
        "Largest-loop size is itself uncertain across the sampled ensemble "
        "(median {loop_nt_median:.0f} nt, IQR {loop_nt_iqr:.0f} nt), so the "
        "target-aware comparison below rests on a soft number.",
        chip="loop size uncertain",
    ),
    Rule(
        "tier2_weak", CAUTION, 60,
        lambda r: r.get("band") == "weak",
        "Largest accessible loop spans approximately {tier2_d_apt_A:.0f} Å, a poor "
        "geometric match to the {tier2_d_pocket_A:.1f} Å cavity detected on the target.",
        chip="geometry: weak",
    ),

    # -- positives --------------------------------------------------------
    Rule(
        "dominant_fold", POSITIVE, 110,
        lambda r: _has(r, "ensemble_defect_score") and r["ensemble_defect_score"] > 0.9,
        "Folds into a single dominant structure — ensemble defect {ensemble_defect:.3f} "
        "sits in the top {ensemble_defect_top_pct:.0f}% of the reference corpus.",
        chip="well-defined fold",
    ),
    Rule(
        "tier2_strong", POSITIVE, 120,
        lambda r: r.get("band") == "strong",
        "Largest accessible loop spans approximately {tier2_d_apt_A:.0f} Å, compatible "
        "with the {tier2_d_pocket_A:.1f} Å cavity detected on the target.",
        chip="geometry: strong",
    ),
    Rule(
        "shuffle_pass", POSITIVE, 130,
        lambda r: r.get("shuffle_pass") is True,
        "Outscores its dinucleotide-shuffled counterparts by {shuffle_margin:.2f} on the "
        "structural criteria (p = {shuffle_p_value:.3f}), so the ranking is not driven "
        "by composition alone.",
        chip="beats shuffled controls",
    ),
    Rule(
        "stable_for_length", POSITIVE, 140,
        lambda r: _has(r, "mfe_norm_score") and r["mfe_norm_score"] > 0.85,
        "Unusually stable for its length at {mfe_norm:.3f} kcal/mol per nucleotide "
        "(top {mfe_norm_top_pct:.0f}% of the corpus).",
        chip="stable",
    ),
    Rule(
        "typical_architecture", POSITIVE, 150,
        lambda r: _has(r, "stem_fraction_score") and r["stem_fraction_score"] > 0.8,
        "Stem/loop balance is typical of validated aptamers "
        "({stem_fraction:.0%} paired, {hairpin_phrase}).",
        chip="typical architecture",
    ),
    Rule(
        "certain_fold", POSITIVE, 160,
        lambda r: _has(r, "positional_entropy_mean_score")
        and r["positional_entropy_mean_score"] > 0.85,
        "Positional entropy averages {positional_entropy_mean:.2f}, so the pairing is "
        "well determined along the whole sequence.",
        chip="low uncertainty",
    ),
    Rule(
        "tier2_moderate", NEUTRAL, 170,
        lambda r: r.get("band") == "moderate",
        "Largest accessible loop spans approximately {tier2_d_apt_A:.0f} Å against a "
        "{tier2_d_pocket_A:.1f} Å cavity — a moderate geometric match relative to "
        "shuffled controls.",
        chip="geometry: moderate",
    ),

    # -- neutral / bookkeeping -------------------------------------------
    Rule(
        "tier2_not_evaluated", NEUTRAL, 200,
        lambda r: r.get("tier2_status") == "not_evaluated",
        "No target-aware evidence was computed for this candidate: it falls outside "
        "the Tier 2 cut.",
        chip="Tier 2 not evaluated",
    ),
    Rule(
        "tier2_not_run", NEUTRAL, 201,
        lambda r: r.get("tier2_status") == "not_run",
        "This run scored intrinsic quality only — no target was supplied, so no "
        "target-aware evidence exists for any candidate.",
        chip="Tier 1 only",
    ),
    Rule(
        "duplicated", NEUTRAL, 210,
        lambda r: (r.get("duplicate_count") or 1) > 1,
        "Submitted {duplicate_count} times in this batch.",
        chip="duplicated input",
    ),
)


def explain(
    record: Mapping[str, Any],
    max_sentences: int = 4,
    max_cautions: int = 2,
    include_caveat: bool = True,
) -> dict[str, Any]:
    """Render one candidate's explanation paragraph and evidence chips.

    Cautions are reserved slots rather than merely high priority: a candidate
    with five glowing positives must still surface its failed shuffle control.
    """
    flat = _flatten(record)

    # A rule only counts as fired if its condition holds *and* every number its
    # sentence quotes is present in the record. Rendering first means an
    # unrenderable rule does not silently consume one of the four slots.
    rendered: list[tuple[Rule, str]] = []
    for rule in RULES:
        if not _safe(rule.condition, flat):
            continue
        try:
            rendered.append((rule, rule.template.format_map(_FormatDefaults(flat))))
        except (KeyError, ValueError, TypeError):
            continue

    cautions = [item for item in rendered if item[0].kind == CAUTION][:max_cautions]
    others = [item for item in rendered if item[0].kind != CAUTION]
    selected = sorted(
        cautions + others[: max(0, max_sentences - len(cautions))],
        key=lambda item: item[0].priority,
    )

    sentences = [sentence for _rule, sentence in selected]
    chips = [
        {"label": rule.chip, "kind": rule.kind} for rule, _s in selected if rule.chip
    ]

    text = " ".join(sentences)
    if include_caveat:
        text = f"{text} {TIER2_CAVEAT}".strip()

    return {
        "text": text,
        "chips": chips,
        "rules_fired": [rule.name for rule, _s in selected],
    }


class _FormatDefaults(dict):
    """Formatting a missing or None value raises, which drops the rule."""

    def __missing__(self, key: str):  # pragma: no cover - defensive
        raise KeyError(key)

    def __getitem__(self, key: str):
        value = super().__getitem__(key)
        if value is None:
            raise KeyError(key)
        return value


def _safe(condition: Callable[[dict[str, Any]], bool], flat: dict[str, Any]) -> bool:
    try:
        return bool(condition(flat))
    except (TypeError, KeyError):
        return False
