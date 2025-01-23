"""
Module to check and verify input files.
TODO:
Move all checks of files and dirs here
Move all checks of datatypes here
Test try-except in read_bed_file
Is dir with all bed files necessary?
"""
import pandas as pd

ALLOWED_METHYLATION_TYPES = ["21839", "m", "h", "a"]


# def check_dir_availability(path_to_dir_with_bed_files):
#     # for get_list_of_methylated_positions resp. calculate_med_coverage
#     if not os.path.isdir(path_to_dir_with_bed_files):
#         print("The path you provided is incorrect or does not exist. Exiting the program.")
#         return


def read_bed_file(bed_file):
    """
    Check structure of bedmethyl file
    :param bed_file:
    :return:
    """
    try:
        bed_df = pd.read_csv(bed_file, sep="\t", header=None, engine="pyarrow")
        return bed_df, True
    except Exception as e:
        print(f"Could not load BED file: {e}")
        return None, False


def has_valid_structure(bed_df):

    # Correct number of columns
    if not bed_df.shape[1] == 10:
        print("BedMethyl file does not have 10 columns. Please check the file.")
        return False

    # Strings - 0: contig, 3: methylation type, 9: statistics
    for col_index in (0, 3, 9):
        if not pd.api.types.is_string_dtype(bed_df[col_index].dtype):
            print(f"BedMethyl file must contain strings in column {col_index + 1}. Please check the file.")
            return False

    # Allowed methylations
    for methylation_type in bed_df[3].unique():
        if methylation_type not in ALLOWED_METHYLATION_TYPES:
            print(f"Methylation {methylation_type} is not supported. Please check the file.")
            return False

    # Integers - 1: index, 2: position, 4: coverage
    for col_index in (1, 2, 4):
        if not pd.api.types.is_integer_dtype(bed_df[col_index].dtype):
            print(f"BedMethyl file must contain integers in column {col_index + 1}. Please check the file.")
            return False

    # Strand
    for strand in bed_df[5].unique():
        if strand not in ("+", "-"):
            print(f"BedMethyl file must contain strand information '+' or '-' in column 6. Please check the file.")
            return False

    # Statistics
    if not bed_df[9].apply(lambda x: len(x.split()) == 9).all():
        print("Each value in column 10 must contain exactly 9 parts separated by spaces.")
        return False

    return True


def load_bed_file(file_path):
    bed_df, is_valid = read_bed_file(file_path)
    if is_valid:
        if has_valid_structure(bed_df):
            return bed_df
            # return bed_df[[0, 1, 2, 3, 4, 5, 9]] # proto to nefunguje!!!!

    return None
