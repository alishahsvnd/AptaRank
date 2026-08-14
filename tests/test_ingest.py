from __future__ import annotations

import pytest

from aptarank.errors import InputError
from aptarank.ingest import ingest, normalise_sequence, validate_sequence


def write(tmp_path, name, text):
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return path


def test_normalisation_handles_dna_style_and_whitespace():
    assert normalise_sequence(" gttc cau\ngg ") == "GUUCCAUGG"


def test_validation_reasons_are_specific():
    assert validate_sequence("ACGUN" * 5, 20, 100) == "invalid character 'N'"
    assert validate_sequence("ACGU", 20, 100) == "length 4 below minimum 20"
    assert validate_sequence("A" * 200, 20, 100) == "length 200 above maximum 100"
    assert validate_sequence("ACGU" * 6, 20, 100) is None


def test_txt_ingest_rejects_without_dropping_silently(tmp_path):
    path = write(
        tmp_path,
        "c.txt",
        "\n".join(["ACGU" * 6, "ACGUN" * 5, "AC", "ACGU" * 6]),
    )
    result = ingest(path, 20, 100)

    assert result.n_submitted == 4
    assert result.n_valid == 1                      # two rejected, one duplicate merged
    assert result.n_rejected == 2
    assert {r["reason"] for r in result.rejections} == {
        "invalid character 'N'",
        "length 2 below minimum 20",
    }
    assert result.candidates.iloc[0]["duplicate_count"] == 2


def test_csv_ingest_uses_supplied_ids(tmp_path):
    path = write(
        tmp_path,
        "c.csv",
        "id,sequence\nalpha,{0}\nbeta,{1}\n".format("ACGU" * 6, "GCAU" * 7),
    )
    result = ingest(path, 20, 100)
    assert list(result.candidates["candidate_id"]) == ["alpha", "beta"]


def test_csv_without_sequence_column_is_an_error(tmp_path):
    path = write(tmp_path, "c.csv", "id,seq\nalpha,ACGUACGUACGUACGUACGUACGU\n")
    with pytest.raises(InputError, match="requires a 'sequence' column"):
        ingest(path, 20, 100)


def test_fasta_ingest(tmp_path):
    path = write(
        tmp_path,
        "c.fasta",
        f">apt1 description here\n{'ACGU' * 6}\n>apt2\n{'GC' * 12}\n",
    )
    result = ingest(path, 20, 100)
    assert list(result.candidates["candidate_id"]) == ["apt1", "apt2"]


def test_all_invalid_input_fails_loudly(tmp_path):
    path = write(tmp_path, "c.txt", "NNNN\nXXXX\n")
    with pytest.raises(InputError, match="all rejected"):
        ingest(path, 20, 100)


def test_the_original_upload_name_is_kept_beside_the_content_hash(tmp_path):
    """§1.6: uploads are stored under a hash, but the provenance panel has to
    show the user a filename they recognise."""
    from aptarank.provenance import write_origin

    staged = tmp_path / "candidates_9f2c1ab77c10.csv"
    staged.write_text("sequence\n" + "ACGUACGUACGUACGUACGUACGU\n", encoding="utf-8")
    write_origin(staged, "my first batch.csv")

    result = ingest(staged, min_length=20, max_length=100)
    assert result.original_filename == "my first batch.csv"
    assert result.summary()["original_filename"] == "my first batch.csv"
    # The real path is still recorded: the hash is what identifies the content.
    assert result.summary()["filename"].endswith("candidates_9f2c1ab77c10.csv")


def test_a_file_with_no_origin_record_falls_back_to_its_own_name(tmp_path):
    path = tmp_path / "plain.csv"
    path.write_text("sequence\n" + "ACGUACGUACGUACGUACGUACGU\n", encoding="utf-8")
    result = ingest(path, min_length=20, max_length=100)
    assert result.summary()["original_filename"] == "plain.csv"


def test_a_tool_version_is_the_line_with_a_version_in_it(monkeypatch):
    """fpocket answers --version with a banner and an error before naming itself,
    and APBS puts its banner on stdout with the version on stderr. Taking the
    first line of whichever stream spoke recorded '***** POCKET HUNTING BEGINS
    *****' as the fpocket version."""
    import subprocess

    from aptarank import provenance

    class _Result:
        def __init__(self, out, err=""):
            self.stdout, self.stderr, self.returncode = out, err, 0

    cases = {
        "fpocket": (
            "***** POCKET HUNTING BEGINS ***** \n! Invalid pdb name given.\n\n"
            ":||: \x1b[1mfpocket 4.0\x1b[0m :||:\n", "",
        ),
        "apbs": ("\n\n---------------------\n    APBS -- Solver\n", "APBS 3.4.1\n"),
    }
    expected = {"fpocket": "fpocket 4.0", "apbs": "APBS 3.4.1"}

    for tool, (out, err) in cases.items():
        monkeypatch.setattr(
            subprocess, "run", lambda *a, _o=out, _e=err, **k: _Result(_o, _e)
        )
        assert provenance._cli_version([tool, "--version"]) == expected[tool]


def test_a_missing_tool_is_recorded_as_absent_not_as_a_guess(monkeypatch):
    import subprocess

    from aptarank import provenance

    def boom(*args, **kwargs):
        raise OSError("no such file")

    monkeypatch.setattr(subprocess, "run", boom)
    assert provenance._cli_version(["nope", "--version"]) is None
