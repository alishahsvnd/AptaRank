from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from aptarank.evaluation import groups as groups_mod
from aptarank.evaluation import stats
from aptarank.tier1.shuffles import kmer_counts


def corpus_frame(n: int = 40) -> pd.DataFrame:
    rng = np.random.default_rng(3)
    rows = []
    for i in range(n):
        length = int(rng.integers(25, 60))
        rows.append(
            {
                "candidate_id": f"ref{i:03d}",
                "sequence": "".join(rng.choice(list("ACGU"), size=length)),
                "target_pdb_id": ["3SPU", "4DII", "1FLT", "1DPX"][i % 4],
                "target_name": ["NDM-1", "thrombin", "VEGF", "lysozyme"][i % 4],
            }
        )
    return pd.DataFrame(rows)


def test_random_control_matches_length_and_letter_composition():
    """An unmatched control would let the tool win trivially on composition."""
    corpus = corpus_frame(200)
    random_frame = groups_mod.random_rna_matched(corpus["sequence"], 400, seed=1)

    real_lengths = np.array([len(s) for s in corpus["sequence"]])
    fake_lengths = np.array([len(s) for s in random_frame["sequence"]])
    assert abs(real_lengths.mean() - fake_lengths.mean()) < 3

    def freqs(seqs):
        joined = "".join(seqs)
        return np.array([joined.count(c) / len(joined) for c in "ACGU"])

    assert np.allclose(freqs(corpus["sequence"]), freqs(random_frame["sequence"]), atol=0.03)


def test_hard_negative_preserves_dinucleotides_of_its_source():
    corpus = corpus_frame(20)
    shuffled = groups_mod.shuffled_negatives(corpus, seed=1)
    assert len(shuffled) == len(corpus)
    for source, control in zip(corpus.itertuples(), shuffled.itertuples()):
        assert control.source_id == source.candidate_id
        assert kmer_counts(control.sequence, 2) == kmer_counts(source.sequence, 2)


def test_folds_never_split_a_target_across_folds():
    """The property that makes E1 non-circular."""
    corpus = corpus_frame(40)
    folds = groups_mod.target_folds(corpus, n_folds=4, seed=7)
    # The biological target name is preferred over the PDB ID: two structures
    # can be the same protein, and grouping by structure would leak.
    assert folds and all(f["grouping"] == "target_name" for f in folds)

    seen: dict[str, int] = {}
    for fold in folds:
        for target in corpus.iloc[fold["test_index"]]["target_name"]:
            assert seen.setdefault(target, fold["fold"]) == fold["fold"]

    covered = sorted(i for fold in folds for i in fold["test_index"])
    assert covered == list(range(len(corpus)))


def test_canonical_target_prefers_the_protein_over_the_structure():
    """Aliased PDB entries of one protein must land in the same fold."""
    corpus = pd.DataFrame(
        {
            "candidate_id": ["a", "b", "c"],
            "sequence": ["ACGU" * 6] * 3,
            "target_pdb_id": ["3SPU", "4EYL", "1FLT"],   # first two are both NDM-1
            "target_name": ["NDM-1", "ndm 1", "VEGF"],
        }
    )
    keys, grouping = groups_mod.canonical_target(corpus)
    assert grouping == "target_name"
    assert keys.iloc[0] == keys.iloc[1] != keys.iloc[2]


def test_folds_fall_back_and_say_so_without_labels():
    corpus = corpus_frame(20).drop(columns=["target_pdb_id", "target_name"])
    folds = groups_mod.target_folds(corpus, n_folds=3, seed=7)
    assert all("no target labels" in f["grouping"] for f in folds)


def test_generated_group_drops_sequences_present_in_the_corpus(tmp_path):
    """A memorised training sequence must not be scored as if it were novel."""
    corpus = corpus_frame(10)
    memorised = corpus.iloc[0]["sequence"]
    path = tmp_path / "generated.csv"
    path.write_text(
        "id,sequence\n"
        f"g1,{memorised}\n"
        f"g2,{'ACGU' * 8}\n",
        encoding="utf-8",
    )
    built = groups_mod.build_groups(corpus, seed=1, generated_path=path)
    assert built.provenance["generated_dropped_as_corpus_duplicates"] == 1
    assert memorised not in set(built["generated"]["sequence"])


def test_probability_of_superiority_is_auc():
    assert stats.probability_of_superiority([3, 4, 5], [0, 1, 2]) == pytest.approx(1.0)
    assert stats.probability_of_superiority([0, 1, 2], [3, 4, 5]) == pytest.approx(0.0)
    assert stats.probability_of_superiority([1, 2, 3], [1, 2, 3]) == pytest.approx(0.5)


def test_compare_groups_reports_effect_size_with_a_ci():
    rng = np.random.default_rng(0)
    result = stats.compare_groups(rng.normal(1, 1, 200), rng.normal(0, 1, 200), n_boot=200)
    assert result["auc"] > 0.6
    lo, hi = result["auc_ci95"]
    assert lo < result["auc"] < hi


def test_auc_ci_widens_when_bootstrapping_units_not_pairs():
    """Each aptamer contributes several decoys; they are not independent."""
    rng = np.random.default_rng(1)
    positives = [[float(rng.normal(0.7, 0.1))] for _ in range(12)]
    negatives = [list(rng.normal(0.5, 0.1, 5)) for _ in range(12)]

    by_unit = stats.bootstrap_auc_by_unit(positives, negatives, n_boot=300, seed=0)
    pooled = stats.compare_groups(
        [p[0] for p in positives], [x for n in negatives for x in n], n_boot=300
    )
    by_unit_width = by_unit["ci95"][1] - by_unit["ci95"][0]
    pooled_width = pooled["auc_ci95"][1] - pooled["auc_ci95"][0]
    assert by_unit_width > pooled_width


def test_retrieval_metrics():
    perfect = stats.mean_reciprocal_rank([0.9, 0.8], [[0.1, 0.2], [0.3]])
    assert perfect["top1_accuracy"] == pytest.approx(1.0)
    assert perfect["mrr"] == pytest.approx(1.0)

    worst = stats.mean_reciprocal_rank([0.1], [[0.9, 0.8]])
    assert worst["top1_accuracy"] == pytest.approx(0.0)
    assert worst["mrr"] == pytest.approx(1 / 3)


def test_clustered_bootstrap_widens_the_ci_for_correlated_positives():
    """Aptamers sharing a target are not independent observations."""
    rng = np.random.default_rng(4)
    values, clusters = [], []
    for target in range(6):                       # 6 targets, 20 near-relatives each
        centre = rng.normal(0.75, 0.12)
        for _ in range(20):
            values.append(float(rng.normal(centre, 0.01)))
            clusters.append(f"T{target}")
    controls = rng.normal(0.5, 0.12, 120)

    naive = stats.compare_groups(values, controls, n_boot=400, seed=0)
    clustered = stats.compare_groups_clustered(values, clusters, controls, n_boot=400, seed=0)

    assert clustered["cluster_bootstrap"] is True
    assert clustered["n_clusters_a"] == 6
    naive_width = naive["auc_ci95"][1] - naive["auc_ci95"][0]
    clustered_width = clustered["auc_ci95"][1] - clustered["auc_ci95"][0]
    assert clustered_width > naive_width


def test_paired_bootstrap_reports_direction_and_spread():
    result = stats.paired_bootstrap([0.1, 0.2, 0.05, -0.02, 0.15], n_boot=300)
    assert result["n"] == 5
    assert result["fraction_positive"] == pytest.approx(0.8)
    assert result["ci95"][0] < result["mean"] < result["ci95"][1]
