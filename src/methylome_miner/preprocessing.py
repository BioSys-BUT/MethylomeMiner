"""
Pre-process bedmethyle file
"""
import statistics
from pathlib import Path

import pandas as pd

from .loading import parse_bed_file


def calculate_med_coverage(bed_dir_path):
    """
    Calculate median coverage for all analyzed samples
    """
    bed_files = list(Path(bed_dir_path).glob("*.bed"))
    medians = []
    for bed_file in bed_files:
        bed_df_part1 = pd.read_csv(bed_file, sep='\t', header=None, engine="pyarrow")
        medians.append(bed_df_part1[4].median())
    return round(statistics.median(medians))


def get_list_of_methylated_positions(bed_file_path, bed_dir_path,
                                    write_to_file=True, output_file_name=None, file_format="json",
                                    min_coverage=None, min_percent_modified=90, strand=None):

    bed_df = parse_bed_file(bed_file_path)

    if min_coverage is None:
        min_coverage = calculate_med_coverage(bed_dir_path)
        print(f"Minimum coverage for methylated position: {min_coverage}")

    if min_percent_modified < 0 or min_percent_modified > 100:
        print("Invalid input, min_percent_modified must be between 0 and 100.")
        return None

    strand_condition = True if strand is None else (bed_df["strand"] == strand)

    df_filter = ((bed_df["score"] >= min_coverage) & (bed_df["percent_modified"] >= min_percent_modified)
                 & strand_condition)

    bed_df_filtered = bed_df[df_filter]

    if bed_df_filtered.empty:
        print("No methylated positions found with the specified conditions.")
        return None
    else:
        if write_to_file:
            if output_file_name is None:
                output_file_name = Path.cwd() / (bed_file_path.stem + "_filtered")
            write_bed_file(bed_df_filtered, output_file_name=output_file_name, file_format=file_format)

        return bed_df_filtered


def write_bed_file(bed_df, output_file_name, file_format="json"):
    output_file_path = output_file_name.with_suffix("." + file_format)

    match file_format:
        case "json":
            pd.DataFrame.to_json(bed_df, output_file_path)
        case "csv":
            bed_df.to_csv(output_file_path, index=False)
        case "tsv":
            bed_df.to_csv(output_file_path, index=False, sep="\t")
        case "bed":
            print("To be implemented.")
        case _:
            print("Unsupported file format. Choose from: json, csv, tsv or bed.")


if __name__ == "__main__":
    # file = Path("input_bed_files/KP825_b53_4mC_5mC_6mA_calls_modifications_whole_run_aligned_sorted_pileup_SUP_qscore.bed")
    # df = parse_bed_file(file)

    bed_files_dir = Path("input_bed_files")

    for bed_file_path in Path(bed_files_dir).glob("*.bed"):
        methylated_positions = get_list_of_methylated_positions(bed_file_path, bed_files_dir, min_coverage=37)
