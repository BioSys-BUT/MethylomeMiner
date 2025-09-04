"""
Module to check and verify input files
"""
from io import StringIO
from pathlib import Path

import pandas as pd

from BCBio import GFF
from Bio import SeqIO

METHYLATIONS_KEY = {"21839": "4mC", "m": "5mC", "h": "5hmCG", "a": "6mA"}

BED_STRUCTURE = {
    "reference_seq": "category",
    "start_index": "int32",
    "end_index": "int32",
    "modified_base_code": "category",
    "score": "int16",  # coverage
    "strand": "category",
    "start_position": "int32",
    "end_position": "int32",
    "color": str,
    "N_valid_cov": "int16",
    "percent_modified": "float16",
    "N_mod": "int16",
    "N_canonical": "int16",
    "N_other_mod": "int16",
    "N_delete": "int16",
    "N_fail": "int16",
    "N_diff": "int16",
    "N_nocall": "int16",
}

ANNOTATION_STRUCTURE = {
    "record_id": "category",
    "gene_id": str,
    "start": "int32",
    "end": "int32",
    "strand": "category",
    "product": str,
    "note": str
}


def parse_bed_file(bed_file_path):
    """
    Check structure of bedMethyl file and load the content

    :param Path bed_file_path: Path to a bedMethyl file with genome-wide single-base methylation data.
    :rtype pd.DataFrame:
    :return: Loaded extended bedMethyl table.
    """
    match bed_file_path.suffix:
        case ".json":
            bed_df = pd.read_json(bed_file_path)
        case ".csv":
            bed_df = pd.read_csv(bed_file_path)
        case ".tsv":
            bed_df = pd.read_csv(bed_file_path, sep="\t")
        case ".bed":
            bed_df_part1 = pd.read_csv(bed_file_path, sep="\t", header=None, engine="pyarrow")
            if bed_df_part1.shape[1] != 10:
                print("Invalid number of tab-separated columns in BedMethyl file. Please check the file.")
                return None
            else:
                bed_df_part2_str = StringIO(bed_df_part1[9].str.cat(sep="\n"))
                bed_df_part2 = pd.read_csv(bed_df_part2_str, sep=" ", header=None, engine="pyarrow")
                if bed_df_part2.shape[1] != 9:
                    print("Invalid number of space-separated columns in BedMethyl file. Please check the file.")
                    return None
                else:
                    bed_df = pd.concat([bed_df_part1, bed_df_part2], axis=1).drop(9, axis=1)
                    bed_df.columns = list(BED_STRUCTURE.keys())
                    bed_df = bed_df.astype(BED_STRUCTURE)
                    bed_df["modified_base_code"] = bed_df["modified_base_code"].cat.rename_categories(METHYLATIONS_KEY)

                    # Strand
                    for strand in bed_df["strand"].cat.categories:
                        if strand not in ("+", "-", "."):
                            print(
                                f"BedMethyl file must contain strand information '+', '-' or '.' in column 6."
                                f" Please check the file.")
                            return None
        case _:
            print("Invalid bedMethyl file format. Valid options are: 'json', 'csv', 'tsv', and 'bed'. "
                  "Please check the file.")
            return None

    return bed_df


def parse_annotation_file(annotation_file_path):
    """
    Parse genome annotation file into custom DataFrame with only CDS features

    Only selected information are stored in the DataFrame:
    - record id
    - feature qualifiers
    - position of the start of the feature
    - position of the end of the feature
    - strand ('+', '-') where the feature is located
    - description of the product of the feature
    - additional notes for the feature

    :param Path annotation_file_path: Path to a file with genome annotation in 'gff' (v3) or 'gbk' file format.
    :rtype pd.DataFrame:
    :return: Table of genome annotation of CDS
    """
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
                raise ValueError(f"Invalid file format: {annotation_file_path}. Allowed formats are: 'gbk', 'gff'")

    cds_data = []
    if not list_of_records:
        print("No annotation records were found.")
        return None
    else:
        for record in list_of_records:
            for feature in record.features:
                if feature.type == "CDS":
                    cds_data.append({
                        "record_id": record.id,
                        "gene_id": feature.qualifiers[custom_qualifier][0],
                        "start": feature.location.start,
                        "end": feature.location.end,
                        "strand": "+" if feature.location.strand == 1 else "-",
                        "product": feature.qualifiers["product"][0],
                        "note": feature.qualifiers["note"][0]
                    })
        cds_df = pd.DataFrame(cds_data)
        cds_df = cds_df.astype(ANNOTATION_STRUCTURE)
        cds_df = cds_df.sort_values(by=["record_id", "start"])
        return cds_df


def check_and_fix_gff(gff_file):
    """
    Fix annotation in gff file format that starts with '<1' instead of '1'.

    :param gff_file: Annotation file opened for reading.
    :rtype str:
    :return: String with corrected annotation ready to be parsed by standard parsing function.
    """
    fixed_lines = []
    for line_number, line in enumerate(gff_file, start=1):
        if "<1" in line:
            fixed_line = line.replace("<1", "1")  # before line.replace("<1","0")
            fixed_lines.append(fixed_line)
            # print(f"Line {line_number}: Replaced '<1' with '1'")
        else:
            fixed_lines.append(line)

    fixed_lines_str = "".join(fixed_lines)
    return StringIO(fixed_lines_str)


def pair_bed_and_annot_files(bed_dir, annot_dir):
    """
    Find matching pairs of bedMethyl and annotation files based on their name

    Match is found based on the prefix of the files. Prefix is every character in front of the first underscore '_'
    character.

    :param bed_dir: Path to a directory with bedMethyl files.
    :param annot_dir: Path to a directory with genome annotations in 'gff' (v3) or 'gbk' file format.
    :rtype dict:
    :return: Prefix oriented dictionary contains a pair of Paths to bedMethyl file and genome annotation file.
    """
    bed_files = list(Path(bed_dir).glob("*.bed"))
    annot_gff_files = list(Path(annot_dir).glob("*.gff"))
    annot_gbk_files = list(Path(annot_dir).glob("*.gbk"))

    bed_group = {file.stem.split("_")[0].split(".")[0]: file for file in bed_files}
    if annot_gff_files:
        annot_group = {file.stem.split("_")[0].split(".")[0]: file for file in annot_gff_files}
    elif annot_gbk_files:
        annot_group = {file.stem.split("_")[0].split(".")[0]: file for file in annot_gbk_files}
    else:
        return None

    paired_files = {}
    for prefix, bed_file in bed_group.items():
        if prefix in annot_group:
            paired_files[prefix] = {"bed_file": bed_file, "annot_file": annot_group[prefix]}
    if paired_files:
        return paired_files
