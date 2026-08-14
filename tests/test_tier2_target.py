"""Target preparation: chain selection, the partner ordering, residue numbering."""

from __future__ import annotations

import pytest

from aptarank.errors import TargetError
from aptarank.tier2 import target as target_mod

ATOM = (
    "ATOM  {serial:5d}  {name:<3s} ALA {chain}{resseq:4d}    "
    "{x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           {element:>2s}\n"
)
OFFSETS = {"N": (-1.2, 0.0, 0.0), "CA": (0.0, 0.0, 0.0), "C": (1.2, 0.4, 0.0),
           "O": (1.4, 1.6, 0.0), "CB": (0.0, -0.9, 1.3)}


def write_complex(path, first_residue: int = 40) -> None:
    """Two protein chains: B is the target, C sits against its far end.

    B is numbered from `first_residue` on purpose — author numbering rarely
    starts at 1, and the residue selector must use the label, not the position.
    """
    lines, serial = [], 1
    for chain, count, origin, start in (("B", 26, 0.0, first_residue), ("C", 22, 4.0, 1)):
        for index in range(count):
            for name, (dx, dy, dz) in OFFSETS.items():
                lines.append(
                    ATOM.format(
                        serial=serial, name=name, chain=chain, resseq=start + index,
                        x=index * 3.8 + dx, y=dy + origin, z=dz, element=name[0],
                    )
                )
                serial += 1
    path.write_text("".join(lines) + "END\n", encoding="utf-8")


@pytest.fixture
def prepared_complex(tmp_path, monkeypatch):
    """`prepare` against a local two-chain structure, with no network."""
    structure = tmp_path / "TEST.pdb"
    write_complex(structure)

    def fake_fetch(identifier, cache_dir, source="pdb"):
        return {
            "identifier": identifier.upper(), "source": source,
            "structure_kind": "predicted" if source == "alphafold" else "experimental",
            "url": f"test://{identifier}", "format": "pdb", "path": str(structure),
            "sha256": "0" * 64, "size_bytes": structure.stat().st_size,
        }

    monkeypatch.setattr(target_mod, "fetch_structure", fake_fetch)

    def run(**kwargs):
        return target_mod.prepare(
            identifier="TEST",
            cache_dir=tmp_path / "cache",
            work_dir=tmp_path / "work",
            **kwargs,
        )

    return run


def test_only_the_requested_chain_survives(prepared_complex):
    prepared = prepared_complex(chain_id="B", partner_chains=["C"])
    written = prepared.path.read_text(encoding="utf-8")
    assert prepared.chain_id == "B"
    assert " B  " in written.replace("ALA", "ALA")   # chain B atoms are present
    assert not any(line[21:22] == "C" for line in written.splitlines()
                   if line.startswith("ATOM"))
    assert prepared.applied["chains_removed"] == ["C"]
    assert prepared.applied["was_multi_chain"] is True


def test_the_partner_confirms_the_site_before_it_is_stripped(prepared_complex):
    """§3.4: the partner says where to look, then it goes away."""
    prepared = prepared_complex(
        chain_id="B", partner_chains=["C"], target_site_residues=[40, 41, 42]
    )
    evidence = prepared.partner_evidence
    assert evidence["computed"] is True
    assert evidence["partner_chains"] == ["C"]
    assert evidence["n_interface_residues"] > 0
    # The configured residues are compared against the measured interface, and
    # any disagreement is reported rather than silently corrected.
    assert set(evidence["configured_site_residues"]) == {40, 41, 42}
    assert "configured_not_in_interface" in evidence


def test_a_multi_chain_structure_carries_the_induced_fit_caveat(prepared_complex):
    """And states it precisely: the coordinates are still the bound-state ones.

    The partner is deleted, not relaxed away, so calling the result "the unbound
    geometry" overstated what happened to the atoms.
    """
    prepared = prepared_complex(chain_id="B", partner_chains=["C"])
    warning = " ".join(prepared.applied["warnings"])
    assert "bound-state conformation" in warning
    assert "not a relaxed unbound structure" in warning


def test_author_numbering_is_used_verbatim(prepared_complex):
    """Residue 40 means the residue labelled 40, not the 40th in the file."""
    prepared = prepared_complex(chain_id="B", target_site_residues=[40])
    assert [r.residue_number for r in prepared.site_residues] == [40]


def test_a_residue_that_is_not_in_the_chain_is_a_loud_failure(prepared_complex):
    with pytest.raises(TargetError, match="not present in this chain"):
        prepared_complex(chain_id="B", target_site_residues=[7])


def test_a_non_integer_residue_says_to_drop_the_letter(prepared_complex):
    with pytest.raises(TargetError, match="write 42, not"):
        prepared_complex(chain_id="B", target_site_residues=["K42"])


def test_an_absent_partner_chain_is_refused(prepared_complex):
    with pytest.raises(TargetError, match="not in the structure"):
        prepared_complex(chain_id="B", partner_chains=["Z"])


def test_a_missing_chain_lists_what_is_available(prepared_complex):
    with pytest.raises(TargetError, match="available"):
        prepared_complex(chain_id="Q")


def test_alphafold_targets_are_recorded_as_predicted(tmp_path, monkeypatch):
    """AlphaFold models are single-chain A, and must be marked predicted."""
    structure = tmp_path / "AF.pdb"
    lines, serial = [], 1
    for index in range(26):
        for name, (dx, dy, dz) in OFFSETS.items():
            lines.append(ATOM.format(serial=serial, name=name, chain="A",
                                     resseq=index + 1, x=index * 3.8 + dx, y=dy,
                                     z=dz, element=name[0]))
            serial += 1
    structure.write_text("".join(lines) + "END\n", encoding="utf-8")

    monkeypatch.setattr(
        target_mod, "fetch_structure",
        lambda identifier, cache_dir, source="pdb": {
            "identifier": identifier, "source": source, "structure_kind": "predicted",
            "url": "test://af", "format": "pdb", "path": str(structure),
            "sha256": "0" * 64, "size_bytes": 1,
        },
    )
    prepared = target_mod.prepare(
        identifier="P17936", source="alphafold", chain_id="A",
        cache_dir=tmp_path / "cache", work_dir=tmp_path / "work",
    )
    assert prepared.structure_kind == "predicted"
    assert prepared.applied["was_multi_chain"] is False


def test_alphafold_has_only_chain_a(prepared_complex):
    with pytest.raises(TargetError, match="single chain"):
        prepared_complex(source="alphafold", chain_id="B")


def test_a_pdb_id_must_look_like_one(tmp_path):
    with pytest.raises(TargetError, match="4-character PDB ID"):
        target_mod.fetch_structure("7WRQXX", tmp_path, source="pdb")


def test_an_unknown_source_is_refused(tmp_path):
    with pytest.raises(TargetError, match="unknown target source"):
        target_mod.fetch_structure("7WRQ", tmp_path, source="modelarchive")
