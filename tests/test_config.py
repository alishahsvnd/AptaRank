from __future__ import annotations

import pytest

from aptarank.config import load_config
from aptarank.errors import ConfigError


def test_defaults_load_and_validate():
    cfg = load_config()
    assert cfg.get("run.seed") == 42
    assert cfg.get("tier1.composite.method") == "corpus_weighted_mean"
    assert set(cfg.active_criteria()) == {
        "mfe_norm", "ensemble_defect", "positional_entropy_mean",
        "stem_fraction", "gc_fraction",
    }


def test_cli_set_overrides_are_typed():
    cfg = load_config(cli_sets=["tier1.shuffle.n_shuffles=99", "run.mode=fast"])
    assert cfg.get("tier1.shuffle.n_shuffles") == 99
    assert cfg.is_fast
    assert not cfg.shuffles_enabled  # fast mode skips controls


def test_unreachable_alpha_is_rejected():
    """20 controls can never reach p <= 0.01; say so instead of always failing."""
    with pytest.raises(ConfigError, match="can never reach alpha"):
        load_config(cli_sets=["tier1.shuffle.alpha=0.01"])


def test_zero_weights_everywhere_is_rejected():
    with pytest.raises(ConfigError, match="no criterion with non-zero weight"):
        load_config(
            overrides={"tier1": {"weights": {k: 0.0 for k in (
                "mfe_norm", "ensemble_defect", "positional_entropy_mean",
                "stem_fraction", "gc_fraction")}}}
        )


def test_unknown_criterion_is_rejected():
    with pytest.raises(ConfigError, match="unknown criteria"):
        load_config(overrides={"tier1": {"weights": {"not_a_feature": 1.0}}})


def test_scoring_signature_tracks_weights_but_not_output_paths():
    base = load_config()
    same = load_config(overrides={"output": {"dir": "somewhere/else"}})
    different = load_config(overrides={"tier1": {"weights": {"gc_fraction": 0.5}}})
    assert base.scoring_signature() == same.scoring_signature()
    assert base.scoring_signature() != different.scoring_signature()
