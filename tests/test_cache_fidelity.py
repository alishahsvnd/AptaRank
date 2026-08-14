"""A cached corpus must score exactly like the build that produced it.

The bug this pins: `to_csv`/`read_csv` is not bit-exact, and `empirical_cdf`
counts *exact* ties with a mid-rank. `mfe_norm` is an energy in units of 0.01
divided by a small integer, so the same value recurs hundreds of times across a
real corpus — nudging those off the tie by one ULP moved a candidate's criterion
score by 0.008 and its composite by 0.0016.

The visible symptom was the worst kind: the first run against a new reference
library scored differently from every later run against the same library, and
nothing said so.
"""

from __future__ import annotations

import numpy as np
import pytest

from aptarank.config import load_config
from aptarank.tier1 import corpus as corpus_mod


@pytest.fixture
def dev_corpus_config(tmp_path, mini_corpus_path):
    return load_config(
        overrides={
            "corpus": {
                "path": str(mini_corpus_path),
                "is_placeholder": True,
                "allow_placeholder": True,
                "cache_dir": str(tmp_path / "cache"),
            },
            "tier1": {"parallel": {"workers": 1}},
        }
    )


def test_a_fresh_build_and_a_cache_load_give_the_same_numbers(dev_corpus_config):
    fresh, _info = corpus_mod.build_or_load(dev_corpus_config)      # computes
    cached, _info = corpus_mod.build_or_load(dev_corpus_config)     # reads back

    assert len(fresh) == len(cached)
    for column in ("mfe_norm", "ensemble_defect", "positional_entropy_mean",
                   "stem_fraction", "gc_fraction", "radius_of_gyration_A"):
        assert np.array_equal(
            fresh[column].to_numpy(dtype=float), cached[column].to_numpy(dtype=float)
        ), f"{column} differs between a fresh build and its own cache"


def test_the_reference_distributions_are_identical(dev_corpus_config):
    """What actually matters: the same candidate must get the same percentile."""
    fresh, info = corpus_mod.build_or_load(dev_corpus_config)
    cached, _ = corpus_mod.build_or_load(dev_corpus_config)
    criteria = dev_corpus_config.active_criteria()

    a = corpus_mod.reference_distributions(fresh, info, criteria)
    b = corpus_mod.reference_distributions(cached, info, criteria)
    for name in criteria:
        probe = np.unique(a.values[name])[:20]
        assert np.array_equal(a.score(name, probe), b.score(name, probe)), name


def test_identifiers_stay_identifiers_through_the_cache(tmp_path, mini_corpus_path):
    """A library numbered 1, 2, 3 must not come back as integers.

    Per-control seeds are derived from the candidate id, so an id that changes
    type between a fresh build and a cache load changes the shuffled controls
    the bank is built from.
    """
    import pandas as pd

    numeric = tmp_path / "numeric_ids.csv"
    frame = pd.read_csv(mini_corpus_path, dtype=str, keep_default_na=False)
    frame["id"] = [str(i + 1) for i in range(len(frame))]
    frame.to_csv(numeric, index=False)

    cfg = load_config(
        overrides={
            "corpus": {"path": str(numeric), "is_placeholder": True,
                       "allow_placeholder": True, "cache_dir": str(tmp_path / "c")},
            "tier1": {"parallel": {"workers": 1}},
        }
    )
    fresh, _ = corpus_mod.build_or_load(cfg)
    cached, _ = corpus_mod.build_or_load(cfg)
    assert fresh["candidate_id"].tolist() == cached["candidate_id"].tolist()
    assert all(isinstance(v, str) for v in cached["candidate_id"])


def test_the_float_format_round_trips_exactly(tmp_path):
    """The mechanism, pinned directly."""
    import pandas as pd

    rng = np.random.default_rng(0)
    values = pd.DataFrame(
        {"x": rng.normal(size=500) / 7.0, "candidate_id": [str(i) for i in range(500)]}
    )
    path = tmp_path / "round_trip.csv"
    corpus_mod.write_feature_cache(values, path)
    back = corpus_mod.read_feature_cache(path)
    assert np.array_equal(values["x"].to_numpy(), back["x"].to_numpy())

    # And the default formatting is what the fix exists to avoid.
    naive = tmp_path / "naive.csv"
    values.to_csv(naive, index=False)
    assert not np.array_equal(
        values["x"].to_numpy(), pd.read_csv(naive)["x"].to_numpy()
    )
