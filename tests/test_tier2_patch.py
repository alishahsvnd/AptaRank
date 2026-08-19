"""Surface-patch measurement: the numbers, and the mistakes it must catch."""

from __future__ import annotations

import pytest

from aptarank.errors import TargetError
from aptarank.tier2 import patch as patch_mod

pytest.importorskip("freesasa", reason="surface mode needs freeSASA")

#: A short poly-alanine helix-ish chain. Real atom names matter: freeSASA
#: classifies by residue and atom name, and an invented name would be dropped.
ATOM_TEMPLATE = (
    "ATOM  {serial:5d}  {name:<3s} ALA A{resseq:4d}    "
    "{x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           {element:>2s}\n"
)
OFFSETS = {
    "N": (-1.2, 0.0, 0.0), "CA": (0.0, 0.0, 0.0), "C": (1.2, 0.4, 0.0),
    "O": (1.4, 1.6, 0.0), "CB": (0.0, -0.9, 1.3),
}


def write_chain(path, n_residues: int = 12, spacing: float = 3.8) -> None:
    lines, serial = [], 1
    for index in range(n_residues):
        base_x = index * spacing
        for name, (dx, dy, dz) in OFFSETS.items():
            lines.append(
                ATOM_TEMPLATE.format(
                    serial=serial, name=name, resseq=index + 1,
                    x=base_x + dx, y=dy, z=dz, element=name[0],
                )
            )
            serial += 1
    path.write_text("".join(lines) + "END\n", encoding="utf-8")


@pytest.fixture
def chain(tmp_path):
    path = tmp_path / "chain.pdb"
    write_chain(path)
    return path


def test_patch_area_is_the_sum_of_its_residue_areas(chain):
    geometry = patch_mod.patch_geometry(chain, [3, 4, 5])
    assert geometry.n_residues == 3
    assert geometry.area_A2 == pytest.approx(sum(geometry.per_residue_area_A2.values()))
    assert geometry.area_A2 > 0
    # The patch is part of a larger chain, so it cannot be all of the surface.
    assert geometry.area_A2 < geometry.total_chain_area_A2


def test_a_bigger_residue_set_measures_a_bigger_patch(chain):
    small = patch_mod.patch_geometry(chain, [3, 4])
    large = patch_mod.patch_geometry(chain, [3, 4, 5, 6, 7])
    assert large.area_A2 > small.area_A2
    assert large.n_atoms > small.n_atoms


def test_an_extended_patch_is_flagged_as_elongated(chain):
    """A run of residues along a straight chain is a groove, not a round patch."""
    geometry = patch_mod.patch_geometry(chain, list(range(1, 13)))
    assert geometry.elongation > 3
    assert geometry.planarity_A < 12.0


def test_a_residue_that_is_not_there_is_an_error_not_a_zero(chain):
    with pytest.raises(TargetError, match="absent from the prepared structure"):
        patch_mod.patch_geometry(chain, [3, 999])


def test_surface_mode_without_residues_says_what_is_missing(chain):
    with pytest.raises(TargetError, match="needs binding-site residues"):
        patch_mod.patch_geometry(chain, [])


def test_an_unknown_area_definition_is_refused(chain):
    with pytest.raises(TargetError, match="unknown patch area definition"):
        patch_mod.patch_geometry(chain, [3], definition="whatever_looks_best")


def test_a_buried_residue_is_reported_rather_than_silently_summed(tmp_path):
    """Zero exposed area almost always means the wrong chain or numbering."""
    path = tmp_path / "buried.pdb"
    write_chain(path, n_residues=12)
    geometry = patch_mod.patch_geometry(
        path, [3, 4], buried_threshold_A2=10_000.0   # force everything "buried"
    )
    assert geometry.buried_residue_numbers == [3, 4]


def test_the_dict_form_carries_everything_the_bundle_needs(chain):
    payload = patch_mod.patch_geometry(chain, [3, 4, 5]).to_dict()
    for key in ("patch_area_A2", "planarity_A", "elongation", "n_residues",
                "residue_numbers", "per_residue_area_A2", "definition",
                "buried_residue_numbers", "shape_warning"):
        assert key in payload
    # Keys are JSON-safe: residue numbers become strings, not ints.
    assert all(isinstance(k, str) for k in payload["per_residue_area_A2"])


def test_overlapping_residues_are_reported_for_the_alternative_definition():
    from aptarank.tier2.fpocket import Residue

    lining = [Residue("A", 3, "", "ALA"), Residue("A", 9, "", "ALA")]
    assert patch_mod.overlapping_residue_numbers(lining, [3, 4, 5]) == [3]
