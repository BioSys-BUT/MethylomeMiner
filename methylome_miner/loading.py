"""
Module to check and verify input files.
TODO:
Move all checks of files and dirs here
Move all checks of datatypes here
Test try-except in read_bed_file
Is dir with all bed files necessary?
"""
from io import StringIO
from pathlib import Path

import pandas as pd

from BCBio import GFF
from Bio import SeqIO

# from mm_analysis_of_coding_and_noncoding_regions_server import process_annotation_file

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


def parse_annotation_file(annotation_file_path):
    # Are custom qualifiers set to the right value? locus_tag is in both
    with open(annotation_file_path) as annot_file:
        match annotation_file_path.suffix:
            case ".gbk":
                list_of_records = list(SeqIO.parse(annot_file, "genbank"))
                custom_qualifier = "locus_tag"
            case ".gff":
                new_annot_file = check_and_fix_gff(annot_file)
                list_of_records = list(GFF.parse(new_annot_file))
                custom_qualifier = "ID"
            case _:
                raise ValueError(f"Invalid file format: {annotation_file_path}. Allowed formats are: .gbk, .gff")

    cds_data = []
    for record in list_of_records:
        for feature in record.features:
            if feature.type == "CDS":
                cds_data.append({
                    "contig": record.id,
                    "gene_id": feature.qualifiers[custom_qualifier][0],
                    "start": feature.location.start,
                    "end": feature.location.end,
                    "strand": "+" if feature.location.strand == 1 else "-",
                    "product": feature.qualifiers["product"][0],
                    "note": feature.qualifiers["note"][0]
                })

    return pd.DataFrame(cds_data)


def check_and_fix_gff(gff_file):
    fixed_lines = []
    for line_number, line in enumerate(gff_file, start=1):
        if "<1" in line:
            fixed_line = line.replace("<1", "1") # before line.replace("<1","0")
            fixed_lines.append(fixed_line)
            print(f"Line {line_number}: Replaced '<1' with '1'")
        else:
            fixed_lines.append(line)

    fixed_lines_str = "".join(fixed_lines)
    return StringIO(fixed_lines_str)
        # fixed_gff_file = f"{os.path.basename(gff_file).split(".")[0]}_fixed.gff"
        # with open(fixed_gff_file, 'w') as file:
        #     file.writelines(fixed_lines)
        # return fixed_gff_file


if __name__ == "__main__":
    file1 = Path("D:\OneDrive - VUT\_AZV_Helca\methylome\gbk_gff_data\gff_data_uncorrected\KP1622_genome.gff")
    file2 = Path("D:\OneDrive - VUT\_AZV_Helca\methylome\gbk_gff_data\genbank_data\KP1622_genome.gbk")
    cds1 = parse_annotation_file(file1)
    cds2 = parse_annotation_file(file2)
    # cds3 = process_annotation_file(str(file1))
    # cds4 = process_annotation_file(str(file2))
    file3 = Path("D:\OneDrive - VUT\_AZV_Helca\methylome\input_roary_output_files_corrected\KP1344_genome.gff")
    cds3 = parse_annotation_file(file3)
    print("Done.")
