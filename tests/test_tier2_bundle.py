from __future__ import annotations

import json

import pytest

from aptarank.errors import TargetError
from aptarank.tier2 import bundle as bundle_mod


def test_bundle_validates_and_round_trips(synthetic_bundle, tmp_path):
    path = bundle_mod.write(synthetic_bundle, tmp_path / "bundles")
    loaded = bundle_mod.load(path)
    assert loaded["bundle_id"] == synthetic_bundle["bundle_id"]
    assert loaded["selection"]["selected_pocket_index"] == 1


def test_bundle_id_ignores_timestamps_but_tracks_science(synthetic_bundle):
    """Two identical builds must agree; a changed measurement must not."""
    same = dict(synthetic_bundle)
    same["created_utc"] = "1999-01-01T00:00:00Z"
    assert bundle_mod.compute_bundle_id(same) == synthetic_bundle["bundle_id"]

    changed = json.loads(json.dumps(synthetic_bundle))
    changed["pockets"][0]["geometry"]["d_pocket_A"] += 1.0
    assert bundle_mod.compute_bundle_id(changed) != synthetic_bundle["bundle_id"]


def test_tampered_bundle_is_rejected_on_load(synthetic_bundle, tmp_path):
    path = bundle_mod.write(synthetic_bundle, tmp_path / "bundles")
    data = json.loads(path.read_text(encoding="utf-8"))
    data["pockets"][0]["geometry"]["d_pocket_A"] = 999.0
    path.write_text(json.dumps(data), encoding="utf-8")

    with pytest.raises(TargetError, match="does not match its contents"):
        bundle_mod.load(path)


def test_bundle_carries_every_pocket_not_just_the_selected_one(synthetic_bundle):
    assert len(synthetic_bundle["pockets"]) == 2
    assert len(synthetic_bundle["selection"]["pocket_evidence"]) == 2


def test_alpha_spheres_are_inline_so_geometry_is_recomputable(synthetic_bundle):
    from aptarank.tier2.geometry import pocket_geometry

    pocket = bundle_mod.selected_pocket(synthetic_bundle)
    spheres = pocket["alpha_spheres"]
    assert spheres and all("center_A" in s and "radius_A" in s for s in spheres)

    recomputed = pocket_geometry(
        [s["center_A"] for s in spheres],
        [s["radius_A"] for s in spheres],
        pocket["fpocket"]["volume_A3"],
    )
    assert recomputed.d_pocket_A == pytest.approx(pocket["geometry"]["d_pocket_A"])


def test_missing_selection_is_rejected(synthetic_bundle):
    broken = json.loads(json.dumps(synthetic_bundle))
    broken["selection"]["selected_pocket_index"] = 99
    with pytest.raises(TargetError, match="not present in the bundle"):
        bundle_mod.validate(broken)


def test_summary_exposes_what_the_artifact_needs(synthetic_bundle):
    summary = bundle_mod.summary(synthetic_bundle)
    assert summary["pdb_id"] == "TEST"
    assert summary["n_pockets"] == 2
    assert summary["pocket_selection"] == "target_site_overlap"
    assert summary["selected_pocket"]["d_pocket_A"] > 0
    assert summary["electrostatics_status"] == "skipped"


def test_find_reports_a_useful_error_when_no_bundle_exists(tmp_path):
    with pytest.raises(TargetError, match="Prepare one with"):
        bundle_mod.find(tmp_path, "3SPU")


def test_the_bundle_id_survives_fpockets_monte_carlo_volume(tmp_path):
    """fpocket re-estimates cavity volume stochastically and offers no seed.

    Two honest builds of identical evidence must still agree on the id, or the
    id certifies nothing. The volume itself is kept as measured.
    """
    from tests.conftest import make_synthetic_bundle

    bundle = make_synthetic_bundle(tmp_path)
    original = bundle["bundle_id"]

    # Same cavity, measured again: volume wobbles a few percent, and everything
    # derived from it moves with it.
    for pocket in bundle["pockets"]:
        pocket["fpocket"]["volume_A3"] *= 1.04
        pocket["geometry"]["d_equiv_A"] *= 1.013
        pocket["geometry"]["envelope_to_equiv_ratio"] *= 0.987
    assert bundle_mod.compute_bundle_id(bundle) == original


def test_a_changed_alpha_sphere_still_changes_the_id(tmp_path):
    """The exclusion must be narrow: real evidence stays covered."""
    from tests.conftest import make_synthetic_bundle

    bundle = make_synthetic_bundle(tmp_path)
    original = bundle["bundle_id"]
    bundle["pockets"][0]["alpha_spheres"][0]["radius_A"] += 0.5
    assert bundle_mod.compute_bundle_id(bundle) != original

    # As does the quantity every pocket-mode band is computed from.
    bundle = make_synthetic_bundle(tmp_path)
    bundle["pockets"][0]["geometry"]["d_pocket_A"] += 0.1
    assert bundle_mod.compute_bundle_id(bundle) != bundle["bundle_id"]


def test_a_bundle_declares_which_fields_its_id_does_not_cover(tmp_path):
    from tests.conftest import make_synthetic_bundle

    bundle = make_synthetic_bundle(tmp_path)
    assert "pockets[].fpocket.volume_A3" in bundle["nondeterministic_fields"]


def test_dropping_a_nondeterministic_field_is_not_the_same_as_carrying_one(tmp_path):
    """Blanked, not deleted: an omitted volume must not hash as a present one."""
    from tests.conftest import make_synthetic_bundle

    bundle = make_synthetic_bundle(tmp_path)
    stripped = json.loads(json.dumps(bundle))
    for pocket in stripped["pockets"]:
        pocket["fpocket"].pop("volume_A3")
    assert bundle_mod.compute_bundle_id(stripped) != bundle["bundle_id"]
