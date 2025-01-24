"""
Pre-process bedmethyle file
"""
import glob
import json
import os
import statistics
from pathlib import Path

import pandas as pd

from loading import parse_bed_file


def calculate_med_coverage(path_to_dir_with_bed_files):
    """
    Calculate median coverage for all analyzed samples
    """
    bed_files = glob.glob(os.path.join(path_to_dir_with_bed_files, '*.bed'))
    medians = []
    for bed_file in bed_files:
        bed_df_part1 = pd.read_csv(bed_file, sep='\t', header=None, engine="pyarrow")
        medians.append(bed_df_part1[4].median())
    return round(statistics.median(medians))


def get_list_of_methylated_position(file_path, path_to_dir_with_bed_files, strand, file_name, min_percent_modified=90, min_coverage=None):
    """
    Get list of filtred methylated positions
    :param file_path:
    :param min_percent_modified:
    :param strand:
    :param file_name:
    :param min_coverage:
    :return:
    """

    processed_bed_df = parse_bed_file(file_path)

    if min_coverage is None:
        min_coverage = calculate_med_coverage(path_to_dir_with_bed_files)
        print(f"Minimum coverage for methylated position: {min_coverage}")

    if isinstance(min_percent_modified, (int, float)):
        if min_percent_modified < 0 or min_percent_modified > 100:
            print("Invalid input, min_percent_modified must be between 0 and 100.")
            return

    strand_condition = True if strand is None else (processed_bed_df['strand'] == strand)

    condition = (processed_bed_df['score'] >= min_coverage) & (
                processed_bed_df['percent_modified'] >= min_percent_modified) & strand_condition

    selected = processed_bed_df[condition]

    if selected.empty:
        print("No methylated positions found with the specified conditions.")
        return

    filename = file_name.split('.')[0]
    output_csv_path = os.path.join(output_dir, f"{filename}_filtered.csv")
    output_json_path = os.path.join(output_dir2, f'{filename}_methylated_positions.json')

    selected.to_csv(output_csv_path, index=False, sep='\t')

    methylated_positions = selected.groupby(['contig', 'methylation type', 'strand'])['index'].apply(list).to_dict()
    methylated_positions_str_keys = {str(key): value for key, value in methylated_positions.items()}

    with open(output_json_path, 'w') as json_file:
        json.dump(methylated_positions_str_keys, json_file, indent=4)

    return methylated_positions, strand


if __name__ == "__main__":
    file = Path("input_bed_files/KP825_b53_4mC_5mC_6mA_calls_modifications_whole_run_aligned_sorted_pileup_SUP_qscore.bed")
    df = load_bed_file(file)
    dir_with_bed_files = Path("input_bed_files")
    # min_coverage_for_filtering = calculate_med_coverage(dir_with_bed_files)

    # output_dir = r'output/csv_files'
    # output_dir2 = r'output/json_files'
    # os.makedirs(output_dir, exist_ok=True)
    # os.makedirs(output_dir2, exist_ok=True)

    for file_name in os.listdir(dir_with_bed_files):
        if file_name.endswith('.bed'):

            full_file_path_to_bed = os.path.join(dir_with_bed_files, file_name)

            methylated_positions, strand = get_list_of_methylated_position(
                full_file_path_to_bed,
                dir_with_bed_files,
                min_percent_modified=90,
                min_coverage=37,
                strand=None,
                file_name=file_name
            )
