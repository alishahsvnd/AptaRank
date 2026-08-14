"""Binding modes: which geometry gets compared (refinements §5).

The old Tier 2 assumed one mechanism — an unpaired loop plugging a cavity — and
presented it as universal. It is not: aptamers also lie a helix along a groove,
or cover a flat surface patch through shape and charge complementarity. Baking
in the pocket model and calling the answer "compatibility" overstated what the
comparison knew.

The tool does **not** infer the mechanism; inferring it would be a research
project. The expert asserts the binding mode, and the comparison adapts to it.
That division of labour is the point: the biologist supplies the judgement, the
tool supplies fast, reproducible, coarse geometry.

Every mode returns the same shape, so nothing downstream has to know which one
ran:

    disagreement    a scalar, lower is better, in the mode's own units
    agreement       a display score in (0, 1]
    fields          mode-specific numbers, all recorded in the artifact

`disagreement` is what the shuffled-control percentile and therefore the band is
computed from. It is deliberately a raw physical mismatch rather than the
Gaussian display score: the Gaussian underflows and depends on `sigma`, and a
display parameter must never be able to move a candidate between bands.
"""

from __future__ import annotations

import math
from typing import Any, Mapping

from ..config import BINDING_MODES, FUTURE_BINDING_MODES, Config
from ..errors import TargetError
from . import apbs
from .geometry import (
    aptamer_dimensions,
    aptamer_footprint_area,
    area_compatibility,
    compatibility,
    footprint_area_from_radius,
)

POCKET = "pocket"
SURFACE = "surface"

#: The candidate/control column each mode compares. Every one already exists in
#: the Tier 1 feature table and in the calibration bank — no new RNA-side
#: computation at scoring time, just a different column (refinements §4.2).
#: Surface mode's column depends on which footprint model is configured.
DESCRIPTOR_COLUMN = {POCKET: "loop_nt_median", SURFACE: "rg_median_A"}
SURFACE_COLUMN_BY_MODEL = {"radius_of_gyration": "rg_median_A", "length": "length"}


def descriptor_column(mode: str, params: Mapping[str, Any]) -> str:
    """Which column this mode reads, given its parameters."""
    if mode == SURFACE:
        model = params.get("footprint_model", "radius_of_gyration")
        return SURFACE_COLUMN_BY_MODEL.get(model, "rg_median_A")
    return DESCRIPTOR_COLUMN[mode]


MISMATCH_UNITS = {POCKET: "A", SURFACE: "A^2"}

PREMISE = {
    POCKET: (
        "the aptamer presents a flexible loop or unpaired region that engages a "
        "concave pocket on the target"
    ),
    SURFACE: (
        "the aptamer covers a flat-ish surface patch through shape and charge "
        "complementarity, rather than plugging a pocket"
    ),
}

#: What the UI calls each mode.
LABEL = {POCKET: "Pocket engagement", SURFACE: "Surface-patch recognition"}


def check_mode(mode: str) -> str:
    if mode in FUTURE_BINDING_MODES:
        raise TargetError(
            f"binding mode {mode!r} is described in the paper as future work and "
            f"is not implemented. Available modes: {BINDING_MODES}."
        )
    if mode not in BINDING_MODES:
        raise TargetError(f"unsupported binding mode: {mode!r}")
    return mode


def parameters(cfg: Config, mode: str) -> dict[str, Any]:
    """The constants this mode uses, exactly as recorded in the artifact."""
    check_mode(mode)
    geometry = {
        "a_per_nt_ss": float(cfg.get("tier2.geometry.a_per_nt_ss")),
        "footprint_scale": float(cfg.get("tier2.geometry.footprint_scale")),
    }
    if mode == POCKET:
        return {
            **geometry,
            "flex_c": float(cfg.get("tier2.flex_c")),
            "sigma_A": float(cfg.get("tier2.sigma_A")),
            "primary_descriptor": cfg.get("tier2.primary_descriptor"),
        }
    weights = cfg.get("tier2.surface.weights")
    return {
        **geometry,
        "footprint_model": cfg.get("tier2.surface.footprint_model"),
        "sigma_A2": float(cfg.get("tier2.surface.sigma_A2")),
        "charge_enabled": bool(cfg.get("tier2.surface.charge.enabled")),
        "charge_scale_kT_per_e": float(cfg.get("tier2.surface.charge.scale_kT_per_e")),
        "weights": {k: float(v) for k, v in weights.items()},
    }


def target_measurement(bundle: Mapping[str, Any], mode: str) -> dict[str, Any]:
    """The target-side numbers this mode compares against.

    Raises rather than substituting a default: a surface-mode run against a
    bundle with no measured patch has nothing to compare, and inventing a number
    there would produce confident output from missing evidence.
    """
    check_mode(mode)
    from . import bundle as bundle_mod

    if mode == POCKET:
        pocket = bundle_mod.selected_pocket(bundle)
        geometry = pocket["geometry"]
        return {
            "d_pocket_A": float(geometry["d_pocket_A"]),
            "d_equiv_A": geometry.get("d_equiv_A"),
            "elongation": geometry.get("elongation"),
            "shape_warning": bool(geometry.get("shape_warning")),
            "pocket_index": pocket["index"],
        }

    patch = bundle.get("patch")
    if not patch or patch.get("patch_area_A2") in (None, 0):
        raise TargetError(
            "this target was prepared without a measured binding-site patch, so "
            "surface mode has nothing to compare. Rebuild the target with "
            "binding-site residues, or score it in pocket mode."
        )
    electrostatics = (bundle.get("electrostatics") or {}).get("sampling") or {}
    return {
        "patch_area_A2": float(patch["patch_area_A2"]),
        "planarity_A": patch.get("planarity_A"),
        "elongation": patch.get("elongation"),
        "n_residues": patch.get("n_residues"),
        "shape_warning": bool(patch.get("shape_warning")),
        "mean_potential_kT_per_e": electrostatics.get("mean_potential_kT_per_e"),
        "electrostatics_status": (bundle.get("electrostatics") or {}).get("status"),
    }


def compare(
    mode: str, value: float, target: Mapping[str, Any], params: Mapping[str, Any]
) -> dict[str, Any]:
    """Compare one candidate descriptor against the target, in `mode`'s terms."""
    check_mode(mode)
    return _POCKET_COMPARE(value, target, params) if mode == POCKET else _SURFACE_COMPARE(
        value, target, params
    )


def _POCKET_COMPARE(
    loop_nt_median: float, target: Mapping[str, Any], params: Mapping[str, Any]
) -> dict[str, Any]:
    """Mode A: does the largest accessible loop span the cavity? (§5.1)"""
    loop = float(loop_nt_median)
    d_pocket = float(target["d_pocket_A"])
    dims = aptamer_dimensions(loop, params["a_per_nt_ss"], params["flex_c"])
    flexible = compatibility(dims["flexible"], d_pocket, params["sigma_A"])
    extended = compatibility(dims["extended"], d_pocket, params["sigma_A"])
    primary = flexible if params["primary_descriptor"] == "flexible" else extended

    return {
        "disagreement": primary["absolute_mismatch_A"],
        "agreement": primary["geometric_score"],
        "fields": {
            "loop_nt_median": loop,
            "d_pocket_A": d_pocket,
            "d_apt_flexible_A": dims["flexible"],
            "signed_mismatch_flexible_A": flexible["signed_mismatch_A"],
            "absolute_mismatch_flexible_A": flexible["absolute_mismatch_A"],
            "geometric_score_flexible": flexible["geometric_score"],
            "d_apt_extended_A": dims["extended"],
            "signed_mismatch_extended_A": extended["signed_mismatch_A"],
            "absolute_mismatch_extended_A": extended["absolute_mismatch_A"],
            "geometric_score_extended": extended["geometric_score"],
            # Display aliases, primary descriptor only.
            "d_apt_A": dims["flexible"],
            "difference_A": flexible["signed_mismatch_A"],
        },
    }


def _SURFACE_COMPARE(
    descriptor: float, target: Mapping[str, Any], params: Mapping[str, Any]
) -> dict[str, Any]:
    """Mode C: is the candidate the right size for the patch, and is the patch
    hospitable to a negatively charged backbone? (§5.2)

    "The right size" is read from the fold by default — the molecule's radius of
    gyration over its element graph — rather than from nucleotide count alone.

    Two signals, combined with documented weights. They are not equivalent in
    kind, and the difference is stated everywhere this number is shown: the
    size term varies per candidate, while the charge term is a property of the
    target and is therefore **identical for every candidate**. The charge term
    can move the reported score; it cannot reorder candidates and cannot change
    a band. That is why `disagreement` — what the band is computed from — is the
    area mismatch alone.
    """
    value = float(descriptor)
    patch_area = float(target["patch_area_A2"])
    model = params.get("footprint_model", "radius_of_gyration")
    if model == "length":
        area = aptamer_footprint_area(
            value, params["a_per_nt_ss"], params["footprint_scale"]
        )
    else:
        area = footprint_area_from_radius(value, params["footprint_scale"])
    size = area_compatibility(area, patch_area, params["sigma_A2"])

    charge = apbs.charge_complementarity(
        target.get("mean_potential_kT_per_e") if params["charge_enabled"] else None,
        params["charge_scale_kT_per_e"],
    )
    weights = dict(params["weights"])
    terms = {"size_coverage": size["geometric_score"]}
    if charge["factor"] is not None:
        terms["charge_complementarity"] = charge["factor"]
    used = {k: weights.get(k, 0.0) for k in terms}
    total = sum(used.values())
    agreement = (
        sum(terms[k] * used[k] for k in terms) / total if total > 0 else float("nan")
    )

    return {
        "disagreement": size["absolute_mismatch_A2"],
        "agreement": agreement,
        "fields": {
            "footprint_model": model,
            "footprint_descriptor": value,
            "radius_of_gyration_A": value if model != "length" else None,
            "footprint_nt": value if model == "length" else None,
            "footprint_area_A2": area,
            "patch_area_A2": patch_area,
            "signed_mismatch_A2": size["signed_mismatch_A2"],
            "absolute_mismatch_A2": size["absolute_mismatch_A2"],
            "coverage_ratio": size["coverage_ratio"],
            "size_coverage_agreement": size["geometric_score"],
            "charge_complementarity": charge["factor"],
            "charge_status": charge["status"],
            "patch_mean_potential_kT_per_e": charge["mean_potential_kT_per_e"],
            "composite_weights": used,
            # Stated in the record itself, not only in the docs: whoever reads
            # this JSON without the paper must still learn it.
            "charge_is_target_level": True,
            "charge_note": (
                "identical for every candidate against this target: it shifts "
                "the score but cannot reorder candidates or change a band"
            ),
        },
    }


def is_evaluable(value: float) -> bool:
    """A missing descriptor is an absence of evidence, not a zero-sized aptamer."""
    try:
        number = float(value)
    except (TypeError, ValueError):
        return False
    return math.isfinite(number) and number > 0


def not_evaluable_status(mode: str) -> str:
    return (
        "not_evaluable_no_contact_loop" if mode == POCKET
        else "not_evaluable_no_footprint"
    )
