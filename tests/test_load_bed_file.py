import pytest
from pathlib import Path

import pandas as pd

from methylome_miner.loading import read_bed_file, has_valid_structure


@pytest.mark.parametrize(
    ["file_path", "expected_is_valid", "bed_pickle"],
    [
        (Path("tests/correct_bed_01.bed"), True, "tests/correct_bed_01.pkl"),
        (Path("tests/incorrect_bed_01_sep.bed"), True, "tests/incorrect_bed_01_sep.pkl")
    ]
)
def test_read_bed_file(file_path, expected_is_valid, bed_pickle):
    bed_df, is_valid = read_bed_file(file_path)
    expected_bed_df = pd.read_pickle(bed_pickle)

    assert bed_df.equals(expected_bed_df)
    assert is_valid == expected_is_valid


@pytest.mark.parametrize(
    ["file_path", "expected_has_valid_structure"],
    [
        (Path("tests/correct_bed_01.bed"), True),
        (Path("tests/correct_bed_02_strand+.bed"), True),
        (Path("tests/correct_bed_03_strand-.bed"), True),
        (Path("tests/correct_bed_04_all_methylations.bed"), True),
        (Path("tests/incorrect_bed_01_sep.bed"), False),
        (Path("tests/incorrect_bed_02_col0_contig.bed"), False),
        (Path("tests/incorrect_bed_03_col1_index.bed"), False),
        (Path("tests/incorrect_bed_04_col2_position.bed"), False),
        (Path("tests/incorrect_bed_05_col3_methylation.bed"), False),
        (Path("tests/incorrect_bed_06_col4_coverage.bed"), False),
        (Path("tests/incorrect_bed_07_col9_stats.bed"), False),
    ]
)
def test_has_valid_structure(file_path, expected_has_valid_structure):
    bed_df, is_valid_file = read_bed_file(file_path)

    assert has_valid_structure(bed_df) == expected_has_valid_structure


# bed_df.to_pickle("incorrect_bed_01_sep.pkl")
