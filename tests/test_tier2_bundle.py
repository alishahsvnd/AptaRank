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
    assert summary["pocket_selection"] == "active_site_overlap"
    assert summary["selected_pocket"]["d_pocket_A"] > 0
    assert summary["electrostatics_status"] == "skipped"


def test_find_reports_a_useful_error_when_no_bundle_exists(tmp_path):
    with pytest.raises(TargetError, match="Build one on Linux"):
        bundle_mod.find(tmp_path, "3SPU")
