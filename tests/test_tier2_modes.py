"""Binding modes: the comparison must adapt to the mode, and only to the mode."""

from __future__ import annotations

import math

import pytest

from aptarank.config import load_config
from aptarank.errors import ConfigError, TargetError
from aptarank.tier2 import modes
from aptarank.tier2.geometry import aptamer_footprint_area, area_compatibility, elongation, planarity

POCKET_TARGET = {"d_pocket_A": 24.0}
SURFACE_TARGET = {"patch_area_A2": 2500.0, "mean_potential_kT_per_e": 6.0}


def params(mode: str, **overrides):
    cfg = load_config(overrides={"tier2": {"binding_mode": mode, **overrides}})
    return modes.parameters(cfg, mode)


# -- dispatch ------------------------------------------------------------


def test_groove_is_rejected_as_future_work_not_as_a_typo():
    """The paper describes it; asking for it must say so rather than 'unknown'."""
    with pytest.raises(TargetError, match="future work"):
        modes.check_mode("groove")
    with pytest.raises(ConfigError, match="future work"):
        load_config(overrides={"tier2": {"binding_mode": "groove"}})


def test_an_unsupported_mode_is_refused():
    with pytest.raises(TargetError, match="unsupported binding mode"):
        modes.check_mode("docking")


def test_each_mode_reads_its_own_rna_descriptor():
    assert modes.descriptor_column(modes.POCKET, params("pocket")) == "loop_nt_median"
    assert modes.descriptor_column(modes.SURFACE, params("surface")) == "rg_median_A"
    assert modes.descriptor_column(
        modes.SURFACE, params("surface", surface={"footprint_model": "length"})
    ) == "length"


def test_modes_report_disagreement_in_their_own_units():
    pocket = modes.compare(modes.POCKET, 16.0, POCKET_TARGET, params("pocket"))
    surface = modes.compare(modes.SURFACE, 28.0, SURFACE_TARGET, params("surface"))
    assert modes.MISMATCH_UNITS[modes.POCKET] == "A"
    assert modes.MISMATCH_UNITS[modes.SURFACE] == "A^2"
    # Both produce the same shape, so nothing downstream needs to know which ran.
    for result in (pocket, surface):
        assert {"disagreement", "agreement", "fields"} <= set(result)
        assert result["disagreement"] >= 0
        assert 0 < result["agreement"] <= 1


# -- pocket mode ---------------------------------------------------------


def test_pocket_mode_compares_loop_reach_against_cavity_width():
    result = modes.compare(modes.POCKET, 16.0, POCKET_TARGET, params("pocket"))
    fields = result["fields"]
    assert fields["d_apt_flexible_A"] == pytest.approx(24.0)   # 6 * sqrt(16)
    assert fields["d_apt_extended_A"] == pytest.approx(96.0)   # 6 * 16
    assert result["disagreement"] == pytest.approx(0.0)


# -- surface mode --------------------------------------------------------


def test_surface_mode_compares_footprint_against_patch_area():
    """Default model: the folded molecule's size, as pi*Rg^2."""
    result = modes.compare(modes.SURFACE, 28.2, SURFACE_TARGET, params("surface"))
    fields = result["fields"]
    assert fields["footprint_area_A2"] == pytest.approx(math.pi * 28.2 ** 2)
    assert fields["patch_area_A2"] == 2500.0
    assert fields["radius_of_gyration_A"] == 28.2
    assert result["disagreement"] < 10.0
    assert fields["coverage_ratio"] == pytest.approx(fields["footprint_area_A2"] / 2500.0)


def test_the_length_model_remains_available_for_comparison():
    result = modes.compare(
        modes.SURFACE, 70.0, SURFACE_TARGET,
        params("surface", surface={"footprint_model": "length"}),
    )
    assert result["fields"]["footprint_area_A2"] == pytest.approx(2520.0)
    assert result["fields"]["footprint_nt"] == 70.0


def test_the_charge_term_is_recorded_as_a_target_level_property():
    """It is the same for every candidate, and the record has to say so."""
    fields = modes.compare(modes.SURFACE, 28.0, SURFACE_TARGET, params("surface"))["fields"]
    assert fields["charge_is_target_level"] is True
    assert "cannot reorder" in fields["charge_note"]

    # Two different candidates, one target: identical charge term.
    other = modes.compare(modes.SURFACE, 12.0, SURFACE_TARGET, params("surface"))["fields"]
    assert fields["charge_complementarity"] == other["charge_complementarity"]


def test_the_band_quantity_ignores_the_charge_term_entirely():
    """`disagreement` drives the band, so it must not move with the charge."""
    hospitable = modes.compare(
        modes.SURFACE, 20.0, {**SURFACE_TARGET, "mean_potential_kT_per_e": 9.0},
        params("surface"),
    )
    repulsive = modes.compare(
        modes.SURFACE, 20.0, {**SURFACE_TARGET, "mean_potential_kT_per_e": -9.0},
        params("surface"),
    )
    assert hospitable["disagreement"] == repulsive["disagreement"]
    assert hospitable["agreement"] > repulsive["agreement"]


def test_a_neutral_patch_gives_exactly_a_half_not_a_flattering_default():
    from aptarank.tier2.apbs import charge_complementarity

    assert charge_complementarity(0.0, 5.0)["factor"] == pytest.approx(0.5)
    assert charge_complementarity(None)["factor"] is None
    # An extreme potential must not overflow the logistic.
    assert charge_complementarity(-1e6, 1.0)["factor"] == 0.0
    assert charge_complementarity(1e6, 1.0)["factor"] == 1.0


def test_surface_mode_refuses_a_target_with_no_measured_patch():
    bundle = {"binding_mode": "surface", "patch": None, "electrostatics": {}}
    with pytest.raises(TargetError, match="without a measured binding-site patch"):
        modes.target_measurement(bundle, modes.SURFACE)


# -- descriptors ---------------------------------------------------------


def test_footprint_area_is_linear_in_length_and_documented():
    assert aptamer_footprint_area(10, 6.0, 1.0) == pytest.approx(360.0)
    assert aptamer_footprint_area(20, 6.0, 1.0) == pytest.approx(720.0)
    assert math.isnan(aptamer_footprint_area(0, 6.0, 1.0))
    assert math.isnan(aptamer_footprint_area(float("nan"), 6.0, 1.0))


def test_area_compatibility_matches_the_pocket_shape_of_the_same_idea():
    result = area_compatibility(1000.0, 1500.0, 600.0)
    assert result["signed_mismatch_A2"] == -500.0
    assert result["absolute_mismatch_A2"] == 500.0
    assert 0 < result["geometric_score"] < 1


def test_elongation_and_planarity_describe_the_shapes_they_claim_to():
    ball = [3.0, 3.0, 3.0]
    channel = [30.0, 3.0, 2.0]
    assert elongation(ball) == pytest.approx(1.0)
    assert elongation(channel) > 10
    # A flat cloud has a small extent along its thinnest axis.
    assert planarity([10.0, 10.0, 0.5], 100) < planarity([10.0, 10.0, 8.0], 100)
    assert math.isnan(planarity([1.0, 1.0], 10))


def test_a_missing_descriptor_is_absence_not_zero():
    assert modes.is_evaluable(12.0) is True
    assert modes.is_evaluable(0.0) is False
    assert modes.is_evaluable(float("nan")) is False
    assert modes.is_evaluable(None) is False
    assert modes.not_evaluable_status(modes.POCKET) == "not_evaluable_no_contact_loop"
    assert modes.not_evaluable_status(modes.SURFACE) == "not_evaluable_no_footprint"
