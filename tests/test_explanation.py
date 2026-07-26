from __future__ import annotations

from aptarank import TIER2_CAVEAT
from aptarank.artifacts.explanation import explain


def record(**overrides):
    base = {
        "candidate_id": "c00001",
        "rank": 1,
        "length": 40,
        "duplicate_count": 1,
        "tier1_score": 0.87,
        "criteria": {
            "mfe_norm": {"value": -0.448, "score": 0.91},
            "ensemble_defect": {"value": 0.0092, "score": 0.96},
            "positional_entropy_mean": {"value": 0.13, "score": 0.88},
            "stem_fraction": {"value": 0.65, "score": 0.78},
            "gc_fraction": {"value": 0.58, "score": 0.64},
        },
        "elements": {"n_hairpins": 1, "loop_nt_median": 7, "loop_nt_iqr": 1},
        "shuffle": {"pass": True, "p_value": 0.048, "margin": 0.31,
                    "n_shuffles": 20, "wins": 20},
        "tier2": {"status": "not_evaluated", "band": "not_evaluated"},
    }
    base.update(overrides)
    return base


def test_caveat_is_always_appended():
    assert explain(record())["text"].endswith(TIER2_CAVEAT)


def test_top_percent_wording_is_the_complement_of_the_score():
    """score 0.96 means better than 96% of the corpus, i.e. the top 4%."""
    text = explain(record())["text"]
    assert "top 4% of the reference corpus" in text


def test_failed_shuffle_is_reported_even_among_positives():
    rec = record(shuffle={"pass": False, "p_value": 0.42, "margin": -0.05,
                          "n_shuffles": 20, "wins": 11})
    out = explain(rec)
    assert "shuffle_fail" in out["rules_fired"]
    assert "does not outscore" in out["text"].lower()
    assert any(chip["kind"] == "caution" for chip in out["chips"])


def test_numbers_come_from_the_record_verbatim():
    rec = record()
    text = explain(rec)["text"]
    assert "0.009" in text                     # ensemble defect, 3 dp
    assert "0.31" in text                      # shuffle margin, 2 dp
    assert "p = 0.048" in text


def test_rules_needing_absent_numbers_do_not_fire():
    rec = record(tier2={"status": "evaluated", "band": "strong"})  # no d_apt / d_pocket
    out = explain(rec)
    assert "tier2_strong" not in out["rules_fired"]
    assert "Å" not in out["text"]


def test_tier2_geometry_sentence_uses_stored_dimensions():
    rec = record(
        tier2={"status": "evaluated", "band": "strong", "d_apt_A": 21.0,
               "d_pocket_A": 18.4, "difference_A": 2.6}
    )
    text = explain(rec)["text"]
    assert "21 Å" in text and "18.4 Å" in text


def test_explanation_is_deterministic():
    rec = record()
    assert explain(rec) == explain(rec)
