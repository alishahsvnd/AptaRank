"""The target-input contract (§3.2) and the key renames it introduced."""

from __future__ import annotations

import pytest

from aptarank.config import load_config, load_target_spec, parse_target_spec
from aptarank.errors import ConfigError

IGFBP3 = """
target_name: IGFBP3
source: pdb
id: 7WRQ
chain: B
binding_mode: surface
partner_chain: C
strip_hetatm: true
target_site_residues: [7, 8, 9, 12]
"""


def test_a_target_description_becomes_config():
    spec = parse_target_spec(IGFBP3, "igfbp3.txt")
    target = spec["tier2"]["target"]
    assert spec["tier2"]["enabled"] is True
    assert spec["tier2"]["binding_mode"] == "surface"
    assert target["id"] == "7WRQ"
    assert target["chain"] == "B"
    assert target["name"] == "IGFBP3"
    assert target["partner_chains"] == ["C"]      # single value becomes a list
    assert target["strip_hetatm"] is True
    assert target["target_site_residues"] == [7, 8, 9, 12]


def test_the_spec_resolves_into_a_valid_config():
    cfg = load_config(overrides=parse_target_spec(IGFBP3, "igfbp3.txt"))
    assert cfg.get("tier2.binding_mode") == "surface"
    assert cfg.get("tier2.target.id") == "7WRQ"
    assert cfg.get("tier2.target.target_site_residues")[:2] == [7, 8]


def test_a_misspelled_key_is_refused_rather_than_ignored():
    """A silently dropped residue list would point at a site nobody chose."""
    with pytest.raises(ConfigError, match="unrecognised key"):
        parse_target_spec("id: 7WRQ\ntarget_site_residue: [7]", "typo.txt")


def test_residues_must_be_plain_integers():
    with pytest.raises(ConfigError, match="write 42, not"):
        parse_target_spec("id: 7WRQ\ntarget_site_residues: [K42]", "letters.txt")


def test_a_description_without_an_id_is_refused():
    with pytest.raises(ConfigError, match="must set `id`"):
        parse_target_spec("chain: B\nbinding_mode: pocket", "no_id.txt")


def test_prose_is_not_a_target_description():
    with pytest.raises(ConfigError, match="key: value"):
        parse_target_spec("please use the IGFBP3 structure", "prose.txt")


def test_a_missing_file_says_so(tmp_path):
    with pytest.raises(ConfigError, match="target file not found"):
        load_target_spec(tmp_path / "nope.txt")


# -- renamed keys --------------------------------------------------------


def test_legacy_keys_are_rewritten_and_the_rewrite_is_recorded():
    cfg = load_config(
        cli_sets=[
            "tier2.target.pdb_id=3SPU",
            "tier2.target.active_site_residues=[120,122]",
            "tier2.a_per_nt=6.5",
        ]
    )
    assert cfg.get("tier2.target.id") == "3SPU"
    assert cfg.get("tier2.target.target_site_residues") == [120, 122]
    assert cfg.get("tier2.geometry.a_per_nt_ss") == 6.5
    renames = [s for s in cfg.sources if s.startswith("<renamed")]
    assert len(renames) == 3


def test_setting_both_the_old_and_new_name_is_an_error():
    with pytest.raises(ConfigError, match="remove 'tier2.target.pdb_id'"):
        load_config(
            overrides={"tier2": {"target": {"pdb_id": "3SPU", "id": "7WRQ"}}}
        )


def test_a_chain_cannot_be_its_own_binding_partner():
    with pytest.raises(ConfigError, match="cannot be its own binding partner"):
        load_config(
            overrides={"tier2": {"target": {"chain": "B", "partner_chains": ["B"]}}}
        )


def test_target_source_is_restricted_to_what_can_be_fetched():
    with pytest.raises(ConfigError, match="must be one of"):
        load_config(overrides={"tier2": {"target": {"source": "swissmodel"}}})


def test_surface_weights_must_leave_at_least_one_signal():
    with pytest.raises(ConfigError, match="no signal with a positive weight"):
        load_config(
            overrides={"tier2": {"surface": {"weights": {
                "size_coverage": 0.0, "charge_complementarity": 0.0}}}}
        )
