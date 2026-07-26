"""Bands must depend on the target and the fixed control bank — nothing else."""

from __future__ import annotations

import numpy as np
import pytest

from aptarank.tier2 import calibration
from aptarank.tier2.calibration import CalibrationBank


def bank(loops) -> CalibrationBank:
    return CalibrationBank(
        bank_id="test-bank", loop_nt_median=np.array(loops, dtype=float), meta={}
    )


def test_control_percentile_is_one_sided_and_higher_is_better():
    """A smaller mismatch must score higher: the percentile counts controls worse."""
    # 6*sqrt(L) = 12, 24, 30, 36 against an 18 Å pocket -> mismatches 6, 6, 12, 18
    dist = calibration.target_distribution(
        bank([4, 16, 25, 36]), d_pocket_A=18.0, a_per_nt=6.0, flex_c=1.0, sigma_A=6.0
    )
    assert dist.percentile(0.0)[0] == pytest.approx(1.0)     # beats every control
    assert dist.percentile(1000.0)[0] == pytest.approx(0.0)  # beaten by every control
    assert dist.percentile(6.0)[0] == pytest.approx(0.75)    # 2 worse, 2 tied: (2+1)/4
    # Monotone: a worse mismatch can never score higher.
    values = [dist.percentile(m)[0] for m in (0.0, 3.0, 6.0, 12.0, 18.0, 30.0)]
    assert values == sorted(values, reverse=True)


def test_percentile_uses_midranks_for_ties():
    dist = calibration.ControlDistribution(
        sorted_absolute_mismatch_A=np.array([1.0, 2.0, 2.0, 3.0]),
        bank_id="b", d_pocket_A=10.0, descriptor="flexible", n=4,
    )
    # one control worse (3.0), two tied at 2.0  ->  (1 + 0.5*2)/4
    assert dist.percentile(2.0)[0] == pytest.approx(0.5)


def test_bands_come_from_the_control_percentile_not_the_gaussian():
    assert calibration.assign_band(0.99, 0.75, 0.95) == "strong"
    assert calibration.assign_band(0.95, 0.75, 0.95) == "strong"
    assert calibration.assign_band(0.80, 0.75, 0.95) == "moderate"
    assert calibration.assign_band(0.10, 0.75, 0.95) == "weak"
    assert calibration.assign_band(None, 0.75, 0.95) == "not_evaluated"


def test_sigma_cannot_move_a_candidate_between_bands():
    """The Gaussian is a display transform; the band must not depend on it."""
    controls = bank(np.arange(4, 60, dtype=float))
    records = []
    for sigma in (1.0, 6.0, 30.0):
        dist_f = calibration.target_distribution(controls, 20.0, 6.0, 1.0, sigma, "flexible")
        dist_e = calibration.target_distribution(controls, 20.0, 6.0, 1.0, sigma, "extended")
        records.append(
            calibration.score_candidates(
                [9, 16, 49], 20.0, dist_f, dist_e, 6.0, 1.0, sigma, 0.75, 0.95
            )
        )
    bands = [[r["band"] for r in group] for group in records]
    assert bands[0] == bands[1] == bands[2]


def test_bands_change_when_the_target_changes():
    """The evidence behind the target-swappability claim (E4)."""
    controls = bank(np.arange(4, 60, dtype=float))
    loops = [4.0, 25.0, 100.0]

    def bands_for(d_pocket):
        dist_f = calibration.target_distribution(controls, d_pocket, 6.0, 1.0, 6.0, "flexible")
        dist_e = calibration.target_distribution(controls, d_pocket, 6.0, 1.0, 6.0, "extended")
        return [
            r["band"]
            for r in calibration.score_candidates(
                loops, d_pocket, dist_f, dist_e, 6.0, 1.0, 6.0, 0.75, 0.95
            )
        ]

    assert bands_for(12.0) != bands_for(58.0)


def test_a_candidate_with_no_contact_loop_is_not_scored_as_zero_sized():
    dist = calibration.target_distribution(bank([9, 16, 25]), 18.0, 6.0, 1.0, 6.0)
    records = calibration.score_candidates([0.0], 18.0, dist, dist, 6.0, 1.0, 6.0, 0.75, 0.95)
    assert records[0]["status"] == "not_evaluable_no_contact_loop"
    assert records[0]["band"] == "not_evaluated"


def test_flexible_is_primary_and_extended_is_kept_as_sensitivity():
    dist_f = calibration.target_distribution(bank([9, 16, 25]), 18.0, 6.0, 1.0, 6.0, "flexible")
    dist_e = calibration.target_distribution(bank([9, 16, 25]), 18.0, 6.0, 1.0, 6.0, "extended")
    record = calibration.score_candidates(
        [16.0], 18.0, dist_f, dist_e, 6.0, 1.0, 6.0, 0.75, 0.95
    )[0]
    assert record["d_apt_A"] == record["d_apt_flexible_A"] == pytest.approx(24.0)
    assert record["d_apt_extended_A"] == pytest.approx(96.0)
    assert record["score"] == record["control_percentile_flexible"]
    assert "control_percentile_extended" in record


def test_thresholds_expose_band_boundaries_in_angstrom():
    dist = calibration.target_distribution(
        bank(np.arange(4, 60, dtype=float)), 20.0, 6.0, 1.0, 6.0
    )
    thresholds = calibration.thresholds(dist, 0.75, 0.95)
    assert thresholds["mismatch_at_strong_A"] < thresholds["mismatch_at_moderate_A"]
    assert thresholds["n_controls"] == dist.n
