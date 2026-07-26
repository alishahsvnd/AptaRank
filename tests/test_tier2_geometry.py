"""Geometry is where a silent error becomes a wrong number in every band."""

from __future__ import annotations

import math

import numpy as np
import pytest

from aptarank.errors import TargetError
from aptarank.tier2.geometry import (
    aptamer_dimensions,
    compatibility,
    pocket_geometry,
)


def sphere_points(n: int = 64, radius: float = 5.0) -> np.ndarray:
    """Fibonacci sphere — an isotropic cloud with no preferred axis."""
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    theta = np.pi * (1 + 5**0.5) * i
    return radius * np.column_stack(
        [np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)]
    )


def test_isotropic_cloud_has_near_equal_extents():
    points = sphere_points()
    geom = pocket_geometry(points, np.full(len(points), 1.0), volume_A3=523.6)
    extents = geom.envelope_extents_A
    assert max(extents) / min(extents) < 1.3


def test_elongated_cloud_reports_its_long_axis():
    rng = np.random.default_rng(0)
    points = rng.normal(0, 1, (200, 3)) * np.array([10.0, 1.0, 1.0])
    geom = pocket_geometry(points, np.full(len(points), 1.0), volume_A3=None)
    assert geom.primary_axis_index == 0
    assert geom.envelope_extents_A[0] > 5 * geom.envelope_extents_A[1]


def test_rotation_does_not_change_the_dimension():
    """d_pocket must be a property of the cavity, not of the PDB's frame."""
    rng = np.random.default_rng(1)
    points = rng.normal(0, 1, (200, 3)) * np.array([8.0, 3.0, 1.0])
    radii = np.full(len(points), 1.5)
    plain = pocket_geometry(points, radii, volume_A3=300.0)

    angle = 0.7
    rot = np.array(
        [[math.cos(angle), -math.sin(angle), 0],
         [math.sin(angle), math.cos(angle), 0],
         [0, 0, 1]]
    )
    rotated = pocket_geometry(points @ rot.T + np.array([50, -20, 7]), radii, volume_A3=300.0)
    assert rotated.d_pocket_A == pytest.approx(plain.d_pocket_A, rel=1e-9)


def test_quantile_extent_resists_a_single_outlier():
    """A min-max extent would be inflated by one stray alpha sphere."""
    points = np.vstack([sphere_points(60, 5.0), np.array([[500.0, 0.0, 0.0]])])
    radii = np.full(len(points), 1.0)
    geom = pocket_geometry(points, radii, volume_A3=523.6)
    assert geom.d_pocket_A < 30.0


def test_envelope_is_wider_than_centres_by_roughly_the_radii():
    points = sphere_points(80, 6.0)
    radii = np.full(len(points), 3.0)
    geom = pocket_geometry(points, radii, volume_A3=904.8)
    assert geom.d_pocket_A > geom.d_pocket_centres_A
    assert geom.d_pocket_A - geom.d_pocket_centres_A == pytest.approx(6.0, abs=1.5)


def test_rms_spread_is_stored_but_is_not_d_pocket():
    """The spec's 2*sv/sqrt(n) is twice the RMS spread, not a diameter."""
    points = sphere_points(80, 6.0)
    geom = pocket_geometry(points, np.full(len(points), 2.0), volume_A3=904.8)
    assert geom.rms_spreads_A[0] < geom.d_pocket_A
    assert len(geom.rms_spreads_A) == 3


def test_degenerate_input_is_rejected_not_guessed():
    with pytest.raises(TargetError, match="at least"):
        pocket_geometry([[0, 0, 0], [1, 1, 1]], [1.0, 1.0], volume_A3=10.0)
    with pytest.raises(TargetError, match="positive"):
        pocket_geometry(sphere_points(8), np.zeros(8), volume_A3=10.0)
    with pytest.raises(TargetError, match="non-finite"):
        bad = sphere_points(8)
        bad[0, 0] = np.nan
        pocket_geometry(bad, np.full(8, 1.0), volume_A3=10.0)


def test_shape_warning_fires_for_a_long_thin_cavity():
    rng = np.random.default_rng(2)
    points = rng.normal(0, 1, (300, 3)) * np.array([25.0, 1.0, 1.0])
    geom = pocket_geometry(points, np.full(len(points), 1.0), volume_A3=50.0)
    assert geom.shape_warning


def test_aptamer_dimensions_flexible_is_primary_and_shorter():
    dims = aptamer_dimensions(loop_nt_median=16, a_per_nt=6.0, flex_c=1.0)
    assert dims["extended"] == pytest.approx(96.0)      # contour length, 6 * 16
    assert dims["flexible"] == pytest.approx(24.0)      # 6 * sqrt(16)
    assert dims["flexible"] < dims["extended"]


def test_compatibility_is_symmetric_and_peaks_at_a_perfect_match():
    exact = compatibility(18.0, 18.0, sigma_A=6.0)
    assert exact["geometric_score"] == pytest.approx(1.0)
    assert exact["absolute_mismatch_A"] == 0.0

    over = compatibility(24.0, 18.0, sigma_A=6.0)
    under = compatibility(12.0, 18.0, sigma_A=6.0)
    assert over["geometric_score"] == pytest.approx(under["geometric_score"])
    assert over["signed_mismatch_A"] == pytest.approx(-under["signed_mismatch_A"])
