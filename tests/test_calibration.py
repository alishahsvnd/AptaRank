"""Bands must depend on the target, the mode and the fixed control bank — nothing else."""

from __future__ import annotations

import numpy as np
import pytest

from aptarank.tier2 import calibration, modes
from aptarank.tier2.calibration import CalibrationBank

POCKET_PARAMS = {
    "a_per_nt_ss": 6.0, "flex_c": 1.0, "sigma_A": 6.0,
    "primary_descriptor": "flexible", "footprint_scale": 1.0,
}
SURFACE_PARAMS = {
    "a_per_nt_ss": 6.0, "footprint_scale": 1.0, "sigma_A2": 600.0,
    "footprint_model": "radius_of_gyration",
    "charge_enabled": True, "charge_scale_kT_per_e": 5.0,
    "weights": {"size_coverage": 1.0, "charge_complementarity": 1.0},
}
LENGTH_PARAMS = {**SURFACE_PARAMS, "footprint_model": "length"}


def bank(loops, lengths=None) -> CalibrationBank:
    loops = np.array(loops, dtype=float)
    return CalibrationBank(
        bank_id="test-bank",
        loop_nt_median=loops,
        length=np.array(lengths if lengths is not None else loops, dtype=float),
        meta={},
    )


def pocket_dist(loops, d_pocket_A, params=None):
    return calibration.target_distribution(
        bank(loops), modes.POCKET, {"d_pocket_A": d_pocket_A}, params or POCKET_PARAMS
    )


def pocket_records(loops, d_pocket_A, params=None, moderate=0.75, strong=0.95):
    params = params or POCKET_PARAMS
    target = {"d_pocket_A": d_pocket_A}
    controls = bank(np.arange(4, 60, dtype=float))
    distribution = calibration.target_distribution(controls, modes.POCKET, target, params)
    secondary = calibration.secondary_distributions(controls, modes.POCKET, target, params)
    return calibration.score_candidates(
        loops, modes.POCKET, target, params, distribution, moderate, strong, secondary
    )


def test_control_percentile_is_one_sided_and_higher_is_better():
    """A smaller mismatch must score higher: the percentile counts controls worse."""
    # 6*sqrt(L) = 12, 24, 30, 36 against an 18 Å pocket -> mismatches 6, 6, 12, 18
    dist = pocket_dist([4, 16, 25, 36], 18.0)
    assert dist.percentile(0.0)[0] == pytest.approx(1.0)     # beats every control
    assert dist.percentile(1000.0)[0] == pytest.approx(0.0)  # beaten by every control
    assert dist.percentile(6.0)[0] == pytest.approx(0.75)    # 2 worse, 2 tied: (2+1)/4
    # Monotone: a worse mismatch can never score higher.
    values = [dist.percentile(m)[0] for m in (0.0, 3.0, 6.0, 12.0, 18.0, 30.0)]
    assert values == sorted(values, reverse=True)


def test_percentile_uses_midranks_for_ties():
    dist = calibration.ControlDistribution(
        sorted_disagreement=np.array([1.0, 2.0, 2.0, 3.0]),
        bank_id="b", mode="pocket", units="A", n=4,
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
    bands = []
    for sigma in (1.0, 6.0, 30.0):
        params = {**POCKET_PARAMS, "sigma_A": sigma}
        bands.append([r["band"] for r in pocket_records([9, 16, 49], 20.0, params)])
    assert bands[0] == bands[1] == bands[2]


def test_bands_change_when_the_target_changes():
    """The evidence behind the target-swappability claim (E4)."""
    loops = [4.0, 25.0, 100.0]
    assert (
        [r["band"] for r in pocket_records(loops, 12.0)]
        != [r["band"] for r in pocket_records(loops, 58.0)]
    )


def test_a_candidate_with_no_contact_loop_is_not_scored_as_zero_sized():
    record = pocket_records([0.0], 18.0)[0]
    assert record["status"] == "not_evaluable_no_contact_loop"
    assert record["band"] == "not_evaluated"


def test_flexible_is_primary_and_extended_is_kept_as_sensitivity():
    record = pocket_records([16.0], 18.0)[0]
    assert record["d_apt_A"] == record["d_apt_flexible_A"] == pytest.approx(24.0)
    assert record["d_apt_extended_A"] == pytest.approx(96.0)
    assert record["score"] == record["control_percentile_flexible"]
    assert "control_percentile_extended" in record


def test_thresholds_expose_band_boundaries_in_the_modes_units():
    dist = pocket_dist(np.arange(4, 60, dtype=float), 20.0)
    thresholds = calibration.thresholds(dist, 0.75, 0.95)
    assert thresholds["disagreement_at_strong"] < thresholds["disagreement_at_moderate"]
    assert thresholds["units"] == "A"
    assert thresholds["n_controls"] == dist.n
    # Pocket mode keeps the names earlier artifacts and the dashboard use.
    assert thresholds["mismatch_at_strong_A"] == thresholds["disagreement_at_strong"]


# -- surface mode --------------------------------------------------------


def surface_target(area=2500.0, potential=7.0):
    return {"patch_area_A2": area, "mean_potential_kT_per_e": potential,
            "planarity_A": 10.0, "elongation": 1.8, "shape_warning": False}


def surface_records(values, target=None, params=None, moderate=0.75, strong=0.95):
    params = params or SURFACE_PARAMS
    target = target or surface_target()
    # Controls spanning plausible molecular sizes: Rg 8-40 A, lengths 20-100 nt.
    controls = bank(np.arange(20, 100, dtype=float),
                    lengths=np.arange(20, 100, dtype=float))
    controls.rg_median_A = np.linspace(8.0, 40.0, 80)
    distribution = calibration.target_distribution(controls, modes.SURFACE, target, params)
    return calibration.score_candidates(
        values, modes.SURFACE, target, params, distribution, moderate, strong
    )


def test_surface_mode_compares_the_candidates_footprint_against_the_patch():
    # pi*Rg^2: an Rg of 28.2 A covers 2500 A^2, exactly the patch.
    records = surface_records([28.2, 10.0])
    assert records[0]["footprint_area_A2"] == pytest.approx(2498.0, abs=5.0)
    assert records[0]["absolute_mismatch_A2"] < records[1]["absolute_mismatch_A2"]
    assert records[0]["band"] == "strong"
    assert records[0]["disagreement_units"] == "A^2"
    assert records[0]["footprint_model"] == "radius_of_gyration"


def test_the_length_model_is_still_available_and_reads_nucleotide_counts():
    # 36 A^2 per nt: a 70-mer covers 2520 A^2, almost exactly the 2500 A^2 patch.
    record = surface_records([70.0], params=LENGTH_PARAMS)[0]
    assert record["footprint_area_A2"] == pytest.approx(2520.0)
    assert record["footprint_model"] == "length"
    assert record["footprint_nt"] == 70.0


def test_the_two_footprint_models_read_different_columns():
    """Which column the bank is asked for follows the configured model."""
    assert modes.descriptor_column(modes.SURFACE, SURFACE_PARAMS) == "rg_median_A"
    assert modes.descriptor_column(modes.SURFACE, LENGTH_PARAMS) == "length"
    assert modes.descriptor_column(modes.POCKET, POCKET_PARAMS) == "loop_nt_median"


def test_the_charge_term_cannot_reorder_or_reband_candidates():
    """It is a property of the target, identical for every candidate (§5.2)."""
    radii = [10.0, 18.0, 28.0, 38.0]
    hospitable = surface_records(radii, surface_target(potential=8.0))
    repulsive = surface_records(radii, surface_target(potential=-8.0))

    assert [r["band"] for r in hospitable] == [r["band"] for r in repulsive]
    assert (
        [r["control_percentile"] for r in hospitable]
        == [r["control_percentile"] for r in repulsive]
    )
    # It does move the reported agreement score, in the direction it should.
    assert hospitable[0]["geometric_agreement_score"] > repulsive[0]["geometric_agreement_score"]
    assert hospitable[0]["charge_is_target_level"] is True


def test_disabling_the_charge_term_leaves_size_agreement_alone():
    params = {**SURFACE_PARAMS, "charge_enabled": False}
    record = surface_records([28.2], params=params)[0]
    assert record["charge_complementarity"] is None
    assert record["charge_status"] == "not_computed"
    assert record["geometric_agreement_score"] == pytest.approx(
        record["size_coverage_agreement"]
    )


def test_surface_bands_change_with_the_patch():
    radii = [12.0, 30.0]
    small = [r["band"] for r in surface_records(radii, surface_target(area=900.0))]
    large = [r["band"] for r in surface_records(radii, surface_target(area=3400.0))]
    assert small != large


def test_a_mode_reads_its_own_column_from_the_bank():
    controls = bank(loops=[10, 20, 30], lengths=[40, 60, 80])
    controls.rg_median_A = np.array([11.0, 15.0, 19.0])
    assert list(controls.descriptor(modes.POCKET, POCKET_PARAMS)) == [10, 20, 30]
    assert list(controls.descriptor(modes.SURFACE, LENGTH_PARAMS)) == [40, 60, 80]
    assert list(controls.descriptor(modes.SURFACE, SURFACE_PARAMS)) == [11.0, 15.0, 19.0]


def test_a_bank_without_the_column_a_mode_needs_says_so():
    controls = bank(loops=[10, 20, 30])
    controls.rg_median_A = None
    with pytest.raises(Exception, match="rg_median_A"):
        controls.descriptor(modes.SURFACE, SURFACE_PARAMS)
