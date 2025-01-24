"""
Module to check and verify input files.
TODO:
Move all checks of files and dirs here
Move all checks of datatypes here
Test try-except in read_bed_file
Is dir with all bed files necessary?
"""
import pandas as pd


METHYLATIONS_KEY = {"21839": "4mC", "m": "5mC", "h": "5hmCG", "a": "6mA"}

# Needs optimization of data types
BED_STRUCTURE = {
    "reference_seq": str,
    "start_index": int,
    "end_index": int,
    "modified_base_code": str,
    "score": int,
    "strand": str,
    "start_position": int,
    "end_position": int,
    "color": str,
    "N_valid_cov": int,
    "percent_modified": float,
    "N_mod": int,
    "N_canonical": int,
    "N_other_mod": int,
    "N_delete": int,
    "N_fail": int,
    "N_diff": int,
    "N_nocall": int,
}


# def check_dir_availability(path_to_dir_with_bed_files):
#     # for get_list_of_methylated_positions resp. calculate_med_coverage
#     if not os.path.isdir(path_to_dir_with_bed_files):
#         print("The path you provided is incorrect or does not exist. Exiting the program.")
#         return


def parse_bed_file(bed_file_path):
    """
    Check structure of bedmethyl file
    :param bed_file_path:
    :return:
    """
    bed_df_part1 = pd.read_csv(bed_file_path, sep="\t", header=None, engine="pyarrow")
    if bed_df_part1.shape[1] != 10:
        print("Invalid number of tab-separated columns in BedMethyl file. Please check the file.")
        return None
    else:
        bed_df_part2 = bed_df_part1[9].str.split(" ", expand=True)
        if bed_df_part2.shape[1] != 9:
            print("Invalid number of space-separated columns in BedMethyl file. Please check the file.")
            return None
        else:
            bed_df = pd.concat([bed_df_part1, bed_df_part2], axis=1).drop(9, axis=1)
            bed_df.columns = list(BED_STRUCTURE.keys())
            bed_df = bed_df.astype(BED_STRUCTURE)

            # Allowed methylations
            for methylation_type in bed_df["modified_base_code"].unique():
                if methylation_type not in METHYLATIONS_KEY:
                    print(f"Methylation {methylation_type} is not supported. Please check the file.")
                    return None
            bed_df["modified_base_code"] = bed_df["modified_base_code"].replace(METHYLATIONS_KEY)

            # Strand
            for strand in bed_df["strand"].unique():
                if strand not in ("+", "-", "."):
                    print(
                        f"BedMethyl file must contain strand information '+', '-' or '.' in column 6."
                        f" Please check the file.")
                    return None

    # bed_df = pd.read_csv(bed_file, sep="\\s+", header=None)
    # if not pd.api.types.is_string_dtype(bed_df[col_index].dtype):
    return bed_df
