"""The fpocket parser cannot be exercised against the real tool on Windows.

These tests pin the parser's behaviour against hand-built fixtures and, above
all, against malformed input: the danger with a whitespace-sensitive block
format is not a crash, it is a plausible-looking wrong number. The pinned Linux
CI job (.github/workflows/target-bundle.yml) additionally runs the real fpocket
and checks the parse against its output.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from aptarank.errors import ExternalToolError
from aptarank.tier2 import fpocket

FIXTURE_OUT = Path(__file__).parent / "fixtures" / "fpocket" / "demo_clean_out"


def test_parses_every_pocket_and_required_metric():
    pockets = fpocket.parse_info((FIXTURE_OUT / "demo_clean_info.txt").read_text())
    assert [p.index for p in pockets] == [1, 2]
    assert pockets[0].score == pytest.approx(0.412)
    assert pockets[0].druggability == pytest.approx(0.783)
    assert pockets[0].volume_A3 == pytest.approx(412.0)
    assert pockets[0].n_alpha_spheres_reported == 6
    assert pockets[1].metrics["Charge score"] == pytest.approx(-1.0)


def test_labels_without_a_space_before_the_colon_still_parse():
    """fpocket is inconsistent: 'Volume score:' has no space, 'Score :' does."""
    pockets = fpocket.parse_info((FIXTURE_OUT / "demo_clean_info.txt").read_text())
    assert pockets[0].metrics["Volume score"] == pytest.approx(4.1)


def test_crlf_line_endings_parse_identically():
    text = (FIXTURE_OUT / "demo_clean_info.txt").read_text()
    assert [p.metrics for p in fpocket.parse_info(text.replace("\n", "\r\n"))] == [
        p.metrics for p in fpocket.parse_info(text)
    ]


def test_unparsed_line_is_an_error_not_a_shrug():
    text = "Pocket 1 :\n\tScore : \t0.4\n\tthis line has no colon\n"
    with pytest.raises(ExternalToolError, match="unparsed line 3"):
        fpocket.parse_info(text)


def test_missing_required_field_is_an_error():
    text = "Pocket 1 :\n\tScore : \t0.4\n\tVolume : \t100.0\n"
    with pytest.raises(ExternalToolError, match="Number of Alpha Spheres"):
        fpocket.parse_info(text)


def test_field_before_any_pocket_header_is_an_error():
    with pytest.raises(ExternalToolError, match="before any"):
        fpocket.parse_info("\tScore : \t0.4\n")


def test_empty_output_is_an_error():
    with pytest.raises(ExternalToolError, match="no pockets"):
        fpocket.parse_info("\n\n")


def test_pqr_coordinates_and_radii():
    spheres = fpocket.parse_pqr((FIXTURE_OUT / "pockets" / "pocket1_vert.pqr").read_text())
    assert len(spheres) == 6
    assert spheres[0].center_A == (0.0, 0.0, 0.0)
    assert spheres[2].center_A == (12.0, 0.0, 0.0)
    assert spheres[2].radius_A == pytest.approx(3.2)
    assert spheres[2].kind == "polar"
    assert spheres[0].kind == "apolar"


def test_pqr_handles_negative_and_scientific_numbers():
    line = "ATOM      1 APOL STP C   1     -12.500  1.0e1  -3.250  0.0  3.500"
    spheres = fpocket.parse_pqr(line)
    assert spheres[0].center_A == (-12.5, 10.0, -3.25)
    assert spheres[0].radius_A == pytest.approx(3.5)


def test_malformed_pqr_line_is_an_error():
    with pytest.raises(ExternalToolError, match="malformed PQR"):
        fpocket.parse_pqr("ATOM      1 APOL STP C   1     x y z 0.0 3.5")


def test_lining_residues_use_fixed_columns_and_deduplicate():
    residues = fpocket.parse_lining_residues(
        (FIXTURE_OUT / "pockets" / "pocket1_atm.pdb").read_text()
    )
    keys = {r.key() for r in residues}
    assert keys == {("A", 120, ""), ("A", 124, ""), ("A", 208, ""), ("A", 301, "")}
    zinc = next(r for r in residues if r.residue_number == 301)
    assert zinc.residue_name == "ZN" and zinc.record_type == "HETATM"


def test_five_digit_residue_numbers_do_not_run_into_the_chain_id():
    """Whitespace splitting breaks here; fixed-column parsing does not."""
    line = "ATOM      1  CA  ALA A1234       0.000   0.000   0.000  1.00 20.00           C"
    residues = fpocket.parse_lining_residues(line)
    assert residues[0].key() == ("A", 1234, "")


def test_load_pockets_cross_checks_the_three_sources():
    pockets = fpocket.load_pockets(FIXTURE_OUT, "demo_clean")
    assert [p.index for p in pockets] == [1, 2]
    assert len(pockets[0].alpha_spheres) == pockets[0].n_alpha_spheres_reported
    assert len(pockets[1].alpha_spheres) == 5


def test_sphere_count_disagreement_is_an_error(tmp_path):
    """info.txt saying 6 spheres while the PQR has 5 must not pass silently."""
    import shutil

    staged = tmp_path / "demo_clean_out"
    shutil.copytree(FIXTURE_OUT, staged)
    pqr = staged / "pockets" / "pocket1_vert.pqr"
    pqr.write_text("\n".join(pqr.read_text().splitlines()[:-1]) + "\n")

    with pytest.raises(ExternalToolError, match="reports 6 alpha spheres"):
        fpocket.load_pockets(staged, "demo_clean")


def test_missing_vertex_file_is_an_error(tmp_path):
    import shutil

    staged = tmp_path / "demo_clean_out"
    shutil.copytree(FIXTURE_OUT, staged)
    (staged / "pockets" / "pocket2_vert.pqr").unlink()

    with pytest.raises(ExternalToolError, match="inconsistent"):
        fpocket.load_pockets(staged, "demo_clean")
