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
