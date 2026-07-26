"""Pocket and aptamer geometry (spec §5.4–5.6).

Every quantity here is a deliberately coarse, transparent approximation. That
is the design: one nameable simplification whose failure mode we can state,
rather than a stack of poorly-characterised error sources.

Two corrections to the spec, both flagged in the README:

* The spec's `2 * singular_value / sqrt(n)` is twice the RMS spread of the
  alpha-sphere centres, not a diameter. We compute a robust projected extent
  instead, and keep the RMS spreads alongside it for comparison.
* The spec makes contour length (`6 Å × L`) the primary aptamer descriptor.
  Contour length assumes a fully extended loop and therefore over-penalises
  long loops, which bend; a flexible chain's end-to-end distance scales closer
  to √L. The √L proxy is primary here and contour length is the upper-bound
  sensitivity check.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Sequence

import numpy as np

from ..errors import TargetError

GEOMETRY_ALGORITHM = "pca-robust-envelope-v1"
LOWER_Q, UPPER_Q = 0.05, 0.95
MIN_ALPHA_SPHERES = 4


@dataclass
class PocketGeometry:
    """Robust dimensions of one cavity, derived from its alpha spheres."""

    algorithm: str
    n_alpha_spheres: int
    centroid_A: list[float]
    pca_axes: list[list[float]]
    singular_values_A: list[float]
    rms_spreads_A: list[float]
    centre_extents_A: list[float]
    envelope_extents_A: list[float]
    d_pocket_centres_A: float
    d_pocket_A: float                 # primary: radius-aware robust envelope
    d_equiv_A: float
    primary_axis_index: int
    envelope_to_equiv_ratio: float
    shape_warning: bool
    quantiles: dict[str, float] = field(
        default_factory=lambda: {"lower": LOWER_Q, "upper": UPPER_Q}
    )

    def to_dict(self) -> dict:
        return {
            "algorithm": self.algorithm,
            "n_alpha_spheres": self.n_alpha_spheres,
            "centroid_A": self.centroid_A,
            "pca_axes": self.pca_axes,
            "singular_values_A": self.singular_values_A,
            "rms_spreads_A": self.rms_spreads_A,
            "centre_extents_A": self.centre_extents_A,
            "envelope_extents_A": self.envelope_extents_A,
            "d_pocket_centres_A": self.d_pocket_centres_A,
            "d_pocket_A": self.d_pocket_A,
            "d_equiv_A": self.d_equiv_A,
            "primary_axis_index": self.primary_axis_index,
            "envelope_to_equiv_ratio": self.envelope_to_equiv_ratio,
            "shape_warning": self.shape_warning,
            "quantiles": {**self.quantiles, "method": "linear"},
        }


def pocket_geometry(
    centres: Sequence[Sequence[float]],
    radii: Sequence[float],
    volume_A3: float | None,
    shape_warning_ratio: float = 2.0,
) -> PocketGeometry:
    """Robust principal-axis dimensions of a cavity.

    The extent along each principal axis is a 5–95% quantile range rather than
    a min-max, so one stray alpha sphere cannot inflate the cavity. The
    radius-aware variant projects sphere surfaces (`z ± r`) rather than centres,
    because what a loop has to reach across is the cavity, not the set of
    points at its middle.
    """
    xyz = np.asarray(centres, dtype=float)
    r = np.asarray(radii, dtype=float)

    if xyz.ndim != 2 or xyz.shape[1] != 3:
        raise TargetError(f"alpha-sphere centres must be (n, 3), got {xyz.shape}")
    if xyz.shape[0] != r.shape[0]:
        raise TargetError(
            f"{xyz.shape[0]} alpha-sphere centres but {r.shape[0]} radii"
        )
    if xyz.shape[0] < MIN_ALPHA_SPHERES:
        raise TargetError(
            f"pocket has only {xyz.shape[0]} alpha spheres; need at least "
            f"{MIN_ALPHA_SPHERES} for a meaningful extent"
        )
    if not np.isfinite(xyz).all() or not np.isfinite(r).all():
        raise TargetError("alpha-sphere coordinates or radii contain non-finite values")
    if (r <= 0).any():
        raise TargetError("alpha-sphere radii must be positive")

    centroid = xyz.mean(axis=0)
    centred = xyz - centroid
    _u, sv, vt = np.linalg.svd(centred, full_matrices=False)

    # Canonical axis signs, so two runs on the same pocket store identical axes.
    axes = np.array([_canonical_sign(axis) for axis in vt])
    projections = centred @ axes.T           # (n, 3)

    centre_extents, envelope_extents = [], []
    for j in range(3):
        z = projections[:, j]
        centre_extents.append(
            float(np.quantile(z, UPPER_Q) - np.quantile(z, LOWER_Q))
        )
        lower = float(np.quantile(z - r, LOWER_Q))
        upper = float(np.quantile(z + r, UPPER_Q))
        envelope_extents.append(upper - lower)

    primary_axis = int(np.argmax(envelope_extents))
    d_pocket = float(envelope_extents[primary_axis])
    d_equiv = (
        2.0 * (3.0 * float(volume_A3) / (4.0 * math.pi)) ** (1.0 / 3.0)
        if volume_A3 and volume_A3 > 0
        else float("nan")
    )
    ratio = d_pocket / d_equiv if d_equiv and math.isfinite(d_equiv) else float("nan")

    return PocketGeometry(
        algorithm=GEOMETRY_ALGORITHM,
        n_alpha_spheres=int(xyz.shape[0]),
        centroid_A=[float(v) for v in centroid],
        pca_axes=[[float(v) for v in axis] for axis in axes],
        singular_values_A=[float(v) for v in sv],
        rms_spreads_A=[float(v) / math.sqrt(xyz.shape[0]) for v in sv],
        centre_extents_A=centre_extents,
        envelope_extents_A=envelope_extents,
        d_pocket_centres_A=float(max(centre_extents)),
        d_pocket_A=d_pocket,
        d_equiv_A=d_equiv,
        primary_axis_index=primary_axis,
        envelope_to_equiv_ratio=float(ratio),
        # An oddly-shaped cavity is flagged, not silently scored: the whole
        # comparison assumes "a loop reaches across a roughly convex pocket".
        shape_warning=bool(
            math.isfinite(ratio) and (ratio > shape_warning_ratio or ratio < 1 / shape_warning_ratio)
        ),
    )


def _canonical_sign(axis: np.ndarray) -> np.ndarray:
    """SVD axis signs are arbitrary; fix them so the bundle is reproducible."""
    dominant = int(np.argmax(np.abs(axis)))
    return -axis if axis[dominant] < 0 else axis


# -- aptamer descriptors (spec §5.5) ------------------------------------


def aptamer_dimensions(
    loop_nt_median: float, a_per_nt: float, flex_c: float
) -> dict[str, float]:
    """Convert a loop size in nucleotides to a reach in Ångström.

    `flexible` is primary: a single-stranded loop is not a rigid rod, and its
    end-to-end distance scales closer to √L than to L. `extended` is contour
    length — the physical upper bound — kept as a sensitivity check.
    """
    length = float(loop_nt_median)
    return {
        "flexible": a_per_nt * flex_c * math.sqrt(length),
        "extended": a_per_nt * length,
    }


def compatibility(d_apt: float, d_pocket: float, sigma_A: float) -> dict[str, float]:
    """Signed and absolute mismatch, plus the Gaussian display score.

    The Gaussian is a smooth display transform so near-misses degrade
    gracefully. The *evidence* is the absolute mismatch and its position in the
    control distribution — the band is derived from that, not from this number,
    so `sigma` cannot move a candidate between bands.
    """
    signed = float(d_apt) - float(d_pocket)
    absolute = abs(signed)
    return {
        "signed_mismatch_A": signed,
        "absolute_mismatch_A": absolute,
        "geometric_score": math.exp(-(absolute**2) / (2.0 * float(sigma_A) ** 2)),
    }
