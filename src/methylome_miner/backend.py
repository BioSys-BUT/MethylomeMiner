"""
Module to filter raw bedMethyl files, store filtered files in various file formats and create panmethylome
"""

from pathlib import Path
import statistics
import sys

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap
import numpy as np

from .loading import METHYLATIONS_KEY, parse_bed_file, parse_annotation_file, pair_bed_and_annot_files


FILTERED_BED_FILE_NAME = "filtered"
CODING_METHYLATIONS_FILE_NAME = "coding"
NON_CODING_METHYLATIONS_FILE_NAME = "non_coding"
ANNOT_WITH_CODING_METHYLATIONS_FILE_NAME = "all_annot_with_coding_methylations"
STATS_FILE_NAME = "methylations_statistics"
STATS_REF_SEQS_FILE_NAME = "methylations_statistics_per_ref"


def calculate_median_coverage(bed_dir_path):
    """
    Calculate median coverage from all available bedMethyl files

    :param Path bed_dir_path: Path to a directory with bedMethyl files.
    :rtype: int
    :return: An integer value of median coverage across all bedMethyl files in input directory.
    """
    bed_files = list(bed_dir_path.glob("*.bed"))

    if not bed_files:
        print(f"Median coverage cannot be calculated. No bedMethyl files found at: {bed_dir_path}")
        return None

    medians = []
    for bed_file in bed_files:
        bed_df = pd.read_csv(bed_file, sep="\t", header=None, engine="pyarrow")
        medians.append(bed_df[4].median())
    return round(statistics.median(medians))


def filter_methylations(bed_file_path, bed_dir_path, min_coverage=None, min_percent_modified=90):
    """
    Filter bedMethyl table based on minimum/median coverage and percent of modified bases.

    Median coverage is calculated from all available bedMethyl files, if minimum coverage is not provided.

    :param Path bed_file_path: Path to a bedMethyl file with genome-wide single-base methylation data.
    :param Path bed_dir_path: Path to a directory with bedMethyl files.
    :param int min_coverage: An integer value of minimum coverage for modified position to be kept.
    :param float min_percent_modified: Minimum required percentage of reads supporting base modification.
    :rtype: pd.DataFrame
    :return: Filtered bedMethyl table
    """

    bed_df = parse_bed_file(bed_file_path)

    if bed_df is None:
        return None

    if min_coverage is None and bed_dir_path is None:
        print("Provide threshold for filtering, use '--min_coverage' or '--input_bed_dir' parameter.")
        return None
    elif bed_dir_path is not None and min_coverage is None:
        min_coverage = calculate_median_coverage(bed_dir_path)
        if min_coverage is None:
            return None
        else:
            print(f"Minimum coverage value not provided. Calculated median coverage {min_coverage} used instead.")

    df_filter = ((bed_df["score"] >= min_coverage) & (bed_df["percent_modified"] >= min_percent_modified))

    bed_df_filtered = bed_df[df_filter]

    if bed_df_filtered.empty:
        print(f"No modified positions found with the specified conditions: "
              f"coverage {min_coverage}, percent modified {min_percent_modified} %.")
        return None
    else:
        return bed_df_filtered


def write_df_to_file(bed_df, file_path):
    """
    Write pandas DataFrame to file of requested file format (options are 'json', 'csv', 'tsv', 'bed').

    :param pd.DataFrame bed_df: DataFrame to be saved.
    :param Path file_path: Path to new file.
    """
    match file_path.suffix:
        case ".json":
            pd.DataFrame.to_json(bed_df, file_path)
        case ".csv":
            bed_df.to_csv(file_path, index=False)
        case ".tsv":
            bed_df.to_csv(file_path, index=False, sep="\t")
        case ".bed":
            to_bed(bed_df, file_path)
        case _:
            print("Unsupported file format. Choose from: 'json', 'csv', 'tsv' or 'bed'.")


def to_bed(bed_df, file_path):
    """
    Save pandas DataFrame as bedMethyl file.

    File structure is the same as the output bedMethyl file from Modkit tool.
    The bedMethyl table contains eighteen columns, but the bedMethyl file has only nine columns that are
    tab-delimited. The remaining columns are stored within the last column and they are space-delimited.

    :param pd.DataFrame bed_df: BedMethyl table to be saved.
    :param Path file_path: Path to new bedMethyl file.
    """
    tab_part = bed_df.iloc[:, :10].astype(str).agg("\t".join, axis=1)
    space_part = bed_df.iloc[:, 10:].astype(str).agg(" ".join, axis=1)

    with open(file_path, mode="w") as file:
        file.write("\n".join(tab_part + " " + space_part) + "\n")


def fill_annot(df, prefix, annot):
    """
    Add columns with annotation to a table with base modifications

    :param pd.DataFrame df: Table of base modifications.
    :param str prefix: Column prefix to distinguish previous/next annotation information.
    :param pd.DataFrame annot: Table with genome feature annotation.
    :rtype pd.DataFrame:
    :return: Table of base modifications with added feature annotation.
    """
    if annot:
        df[f"{prefix}_gene_reference_seq"] = annot.record_id
        df[f"{prefix}_gene_id"] = annot.gene_id
        df[f"{prefix}_gene_start"] = annot.start
        df[f"{prefix}_gene_end"] = annot.end
        df[f"{prefix}_gene_strand"] = annot.strand
        df[f"{prefix}_gene_product"] = annot.product
    else:
        # For the case of modifications in front of the first coding feature and after the last coding feature
        for col in ["reference_seq", "id", "start", "end", "strand", "product"]:
            df[f"{prefix}_gene_{col}"] = ""
    return df


def add_annot_to_non_coding(non_coding_df, prev_annot, next_annot, annot_position="middle"):
    """
    Add annotation to a table of non-coding modifications based on their position within the genome

    Modifications in front of the first coding feature have empty annotation marked as 'prev' (previous).
    Modifications after the last coding feature have empty annotation marked as 'next'.

    :param pd.DataFrame non_coding_df: Table of non-coding base modifications.
    :param pd.DataFrame prev_annot: Table of the previous annotation.
    :param pd.DataFrame next_annot: Table of the next annotation.
    :param str annot_position: Position of the modifications within the genome annotation ('start', 'middle', 'end').
    :rtype pd.DataFrame:
    :return: Table of non-coding base modifications with added annotation information.
    """
    match annot_position:
        case "start":
            fill_annot(non_coding_df, "prev", None)
            fill_annot(non_coding_df, "next", next_annot)
        case "middle":
            fill_annot(non_coding_df, "prev", prev_annot)
            fill_annot(non_coding_df, "next", next_annot)
        case "end":
            fill_annot(non_coding_df, "prev", prev_annot)
            fill_annot(non_coding_df, "next", None)

    return non_coding_df


def add_annot_to_coding(coding_df, annot):
    """
    Add annotation to a table of coding modifications.

    :param pd.DataFrame coding_df: Table of coding base modifications.
    :param pd.DataFrame annot: Table with annotation.
    :rtype pd.DataFrame:
    :return: Table of coding base modifications with added annotation information.
    """
    coding_df["gene_reference_seq"] = annot.record_id
    coding_df["gene_id"] = annot.gene_id
    coding_df["gene_start"] = annot.start
    coding_df["gene_end"] = annot.end
    coding_df["gene_strand"] = annot.strand
    coding_df["gene_product"] = annot.product
    return coding_df


def add_coding_to_annot(annot, methylations_df=None):
    """
    Add coding base modifications to annotation as a new columns (one column per modification type)

    :param pd.DataFrame annot: Table with annotation.
    :param pd.DataFrame methylations_df: Table of coding base modifications.
    :rtype pd.DataFrame:
    :return: Table with annotation and coding base modifications.
    """
    new_annot = pd.DataFrame([annot])
    # Remove 'm' from the names of columns with modification name
    new_col_names = new_annot.columns[:7].union([cat[1:] for cat in new_annot.columns[7:]], sort=False)
    new_annot.columns = new_col_names

    # Add empty lists to each modification type
    for name in new_col_names[7:]:
        new_annot[name] = [[]]

    # Add coding modifications to annotation
    if methylations_df is not None:
        for group, methylations_group in methylations_df.groupby("modified_base_code", observed=False):
            new_annot.at[0, group] = methylations_group["start_index"].tolist()

    return new_annot


def sort_coding_non_coding_methylations(bed_df, annot_df):
    """
    Sort base modifications based on annotation to coding and non-coding types and add coding to annotation.

    'Coding' modifications are bases within a coding region according to the annotation. Each modification receive
    information about corresponding coding region.

    'Non-coding' modifications are in the intergenic region. Each modification receive information about the previous
    and the next coding region.

    :param pd.DataFrame bed_df: Table with base modification information.
    :param pd.DataFrame annot_df: Table with genome annotation.
    :rtype (pd.DataFrame, pd.DataFrame, pd.DataFrame):
    :return: coding_df, non_coding_df, new_annot_df

        - coding_df - Table of all coding modifications with coding feature annotation.
        - non_coding_df - Table of all non-coding modifications with information about neighboring coding features.
        - new_annot_df - Table of all annotation features with information about found coding modifications.
    """
    bed_df = bed_df.reset_index(drop=True)

    coding_parts = []
    non_coding_parts = []
    new_annot_parts = []
    next_methylation_search_idx = 0

    # Iterate over annotation features
    for idx, annot in enumerate(annot_df.itertuples(index=False)):

        # Get modifications that were not sorted yet
        bed_slice = bed_df.iloc[next_methylation_search_idx:]

        # Find modifications that belong to current coding feature of annotation
        coding_part = (
            bed_slice
            .loc[lambda df: df["start_index"] <= annot.end]
            .loc[lambda df: annot.start <= df["start_index"]]
            .loc[lambda df: df["strand"] == annot.strand]
            .loc[lambda df: df["reference_seq"] == annot.record_id]
        ).copy()

        # Get index were coding modifications start
        coding_start_index = coding_part.index.min() if not coding_part.empty else bed_slice.index.min()

        # Non-coding modifications:
        # Modifications are in front of the first coding feature
        if idx == 0:
            non_coding_part = bed_df.loc[lambda df: df["start_index"] < annot.start].copy()
            if not non_coding_part.empty:
                non_coding_part = add_annot_to_non_coding(non_coding_part, "", annot, "start")
                next_methylation_search_idx = non_coding_part.index.max() + 1
                non_coding_parts.append(non_coding_part)

        # Modifications are in front of some coding feature
        elif coding_start_index > next_methylation_search_idx:
            non_coding_part = bed_df.iloc[next_methylation_search_idx:coding_start_index].copy()
            if not non_coding_part.empty:
                non_coding_part = add_annot_to_non_coding(non_coding_part, prev_annot, annot)
                next_methylation_search_idx = non_coding_part.index.max() + 1
                non_coding_parts.append(non_coding_part)

        # Modifications are after the last coding feature
        if idx == annot_df.shape[0] - 1:
            tail_start_idx = coding_part.index.max() + 1 if not coding_part.empty else next_methylation_search_idx
            non_coding_part = bed_df.loc[tail_start_idx:].copy()
            if not non_coding_part.empty:
                non_coding_part = add_annot_to_non_coding(non_coding_part, prev_annot, annot, "end")
                non_coding_parts.append(non_coding_part)

        # Add information from annotation to modifications
        if not coding_part.empty:
            new_annot_parts.append(add_coding_to_annot(annot, coding_part))
            coding_part = add_annot_to_coding(coding_part, annot)
            coding_parts.append(coding_part)

            next_methylation_search_idx = coding_part.index.max() + 1
        else:
            new_annot_parts.append(add_coding_to_annot(annot))

        prev_annot = annot

    # Concat all partial DataFrames
    coding_df = pd.concat(coding_parts)
    non_coding_df = pd.concat(non_coding_parts)
    new_annot_df = pd.concat(new_annot_parts)

    return coding_df, non_coding_df, new_annot_df

def sort_methylations(bed_df, annot_df):
    """
    Sort base modifications according to annotation consecutively by reference sequence name and strand.

    :param pd.DataFrame bed_df: Table of base modifications.
    :param pd.DataFrame annot_df: Table with annotation.
    :rtype (pd.DataFrame, pd.DataFrame, pd.DataFrame):
    :return: coding_df, non_coding_df, new_annot_df

        - coding_df - Table of all coding modifications with coding feature annotation.
        - non_coding_df - Table of all non-coding modifications with information about neighboring coding features.
        - new_annot_df - Table of all annotation features with information about found coding modifications.
    """
    # Group bedMethyl table by reference sequence and strand
    bed_df_grouped = bed_df.groupby(["reference_seq", "strand"], observed=False)
    # Group annotation table by sequence record id and strand
    annot_df_grouped = annot_df.groupby(["record_id", "strand"], observed=False)

    all_coding = []
    all_non_coding = []
    new_annots = []
    # Get all base modifications present the bedMethyl table
    methylation_types = sorted(bed_df["modified_base_code"].cat.categories)
    # Add 'm' to the start of new column names (methylations types)
    # Have to by done, because of itertuples method, that create a Pandas object, where column names cannot start
    # with a number.
    methylation_types_edited = ["m" + cat for cat in methylation_types]

    # Go through each group and sort modifications
    for group, bed_df_group in bed_df_grouped:
        # Add new columns named after base modifications to annotation
        annot_df_group = annot_df_grouped.get_group(group)
        annot_df_group = annot_df_group.reindex(columns=annot_df_group.columns.to_list() + methylation_types_edited)
        # annot_df_group = pd.concat([annot_df_group, pd.DataFrame(columns=methylation_types, dtype="object")])

        coding_df, non_coding_df, new_annot_df = sort_coding_non_coding_methylations(bed_df_group, annot_df_group)

        all_coding.append(coding_df)
        all_non_coding.append(non_coding_df)
        new_annots.append(new_annot_df)

    all_coding_df = pd.concat(all_coding).sort_values(by=["reference_seq", "start_index"])
    all_non_coding_df = pd.concat(all_non_coding).sort_values(by=["reference_seq", "start_index"])
    new_annot_df = pd.concat(new_annots).sort_values(by=["record_id", "start"])

    return all_coding_df, all_non_coding_df, new_annot_df


def get_methylations_stats(coding_df, non_coding_df, file_path):
    """
    Get counts statistics per methylation types and per reference sequence in annotation

    :param pd.DataFrame coding_df: Table of all coding modifications with coding feature annotation.
    :param pd.DataFrame non_coding_df: Table of all non-coding modifications with coding feature annotation.
    :param Path file_path: Path of output files of MethylomeMiner
    """

    # Counts statistics per methylation type
    coding_stats = coding_df["modified_base_code"].value_counts()
    non_coding_stats = non_coding_df["modified_base_code"].value_counts()

    methylations = sorted(
        set(coding_df["modified_base_code"].unique()).union(non_coding_df["modified_base_code"].unique())
    )
    ref_seqs = sorted(set(coding_df["reference_seq"].unique()).union(non_coding_df["reference_seq"].unique()))

    stats = {"genome": file_path.stem}
    total_count = 0
    for methylation in methylations:
        stats[f"{methylation}_coding_count"] = coding_stats[methylation]
        stats[f"{methylation}_non_coding_count"] = non_coding_stats[methylation]
        stats[f"{methylation}_total_count"] = coding_stats[methylation] + non_coding_stats[methylation]
        total_count = total_count + stats[f"{methylation}_total_count"]
    stats["total_count"] = total_count
    stats_df = pd.DataFrame([stats])

    # Count statistics per methylation type and per reference sequence in annotation
    ref_seq_coding_stats = coding_df.groupby("reference_seq", observed=False)["modified_base_code"].value_counts()
    ref_seq_non_coding_stats = non_coding_df.groupby("reference_seq", observed=False)["modified_base_code"].value_counts()

    ref_seq_stats = []
    for ref_seq in ref_seqs:
        row = {"genome": file_path.stem, "reference_seq": ref_seq}
        total_count = 0
        for methylation in methylations:
            row[f"{methylation}_coding_count"] = ref_seq_coding_stats[ref_seq][methylation]
            row[f"{methylation}_non_coding_count"] = ref_seq_non_coding_stats[ref_seq][methylation]
            row[f"{methylation}_total_count"] = ref_seq_coding_stats[ref_seq][methylation] + \
                                                ref_seq_non_coding_stats[ref_seq][methylation]
            total_count = total_count + row[f"{methylation}_total_count"]
        row["total_count"] = total_count
        ref_seq_stats.append(row)
    ref_seq_stats_df = pd.DataFrame(ref_seq_stats)

    # Write statistics to CSV files
    write_df_to_file(stats_df, file_path.with_stem(f"{file_path.stem}_{STATS_FILE_NAME}.csv"))
    write_df_to_file(ref_seq_stats_df, file_path.with_stem(f"{file_path.stem}_{STATS_REF_SEQS_FILE_NAME}.csv"))


def _mine_methylations(input_bed_file, input_annot_file, input_bed_dir, min_coverage, min_percent_modified,
                       work_dir, file_name, write_filtered_bed, filtered_bed_format, split_by_reference,
                       write_all=True):
    """
    Filter modified bases stored in bedMethyl file and sort them according to annotation into coding and non-coding.

    Filtration is performed based on requested coverage and the percent of modified bases.
    Sorting is conducted based on provided annotation of the genome. Modified bases are sorted into coding
    (modification is within coding region) and non-coding (modification is in intergenic region) groups.


    :param Path input_bed_file: Path to a bedMethyl file with genome-wide single-base methylation data.
    :param Path input_annot_file: Path to a file with genome annotation in 'gff' (v3) or 'gbk' file format.
    :param Path input_bed_dir: Path to a directory with bedMethyl files.
    :param int min_coverage: An integer value of minimum coverage for modified position to be kept.
    :param float min_percent_modified: Minimum required percentage of reads supporting base modification. Default: 90
    :param Path work_dir: Path to directory for MethylomeMiner outputs. Default: MethylomeMiner_output
    :param str file_name: Custom name for MethylomeMiner outputs.
    :param bool write_filtered_bed: Write filtered bedMethyl file to a new file. Default: False
    :param str filtered_bed_format: File format for filtered bedMethyl file. Default: csv
    :param bool split_by_reference: Write all outputs to separate files based on reference sequence.
    :param bool write_all: True to write all DataFrames (coding, non-coding, extended annotation) to CSV files.
    :rtype pd.DataFrame:
    :return: DataFrame of extended annotation with added base modification positions.
    """
    print(f"Mining of {input_bed_file} started.")

    # Filter out invalid base modifications
    filtered_bed_df = filter_methylations(
        input_bed_file,
        input_bed_dir,
        min_coverage=min_coverage,
        min_percent_modified=min_percent_modified,
    )

    if filtered_bed_df is None:
        return None
    else:
        # Check (and create) working directory for MethylomeMiner
        if not work_dir.exists():
            work_dir.mkdir(parents=True, exist_ok=True)

        # Check (and create) name of output files of MethylomeMiner
        if file_name is None:
            file_path = Path(work_dir, input_bed_file.stem.split("_")[0].split(".")[0])
        else:
            file_path = Path(work_dir, file_name)

        # Save filtered bedMethyl table if requested
        if write_filtered_bed:
            bed_file_path = file_path.with_stem(f"{file_path.stem}_{FILTERED_BED_FILE_NAME}.{filtered_bed_format}")
            write_df_to_file(filtered_bed_df, bed_file_path)

        # Parse genome annotation file
        annot_df = parse_annotation_file(input_annot_file)
        if annot_df is not None:

            # According to annotation sort modifications into coding and non-coding groups
            coding_df, non_coding_df, new_annot_df = sort_methylations(filtered_bed_df, annot_df)

            # Get basic methylations statistics
            get_methylations_stats(coding_df, non_coding_df, file_path)

            # Split outputs by reference sequences from annotation
            if split_by_reference:

                reference_seqs = new_annot_df["record_id"].astype("category").cat.categories
                for ref_seq in reference_seqs:
                    coding_df_part = coding_df[coding_df["reference_seq"] == ref_seq]
                    non_coding_df_part = non_coding_df[non_coding_df["reference_seq"] == ref_seq]
                    new_annot_df_part = new_annot_df[new_annot_df["record_id"] == ref_seq]

                    write_df_to_file(coding_df_part, file_path.with_stem(
                        f"{file_path.stem}_{ref_seq}_{CODING_METHYLATIONS_FILE_NAME}.csv"))
                    write_df_to_file(non_coding_df_part, file_path.with_stem(
                        f"{file_path.stem}_{ref_seq}_{NON_CODING_METHYLATIONS_FILE_NAME}.csv"))
                    write_df_to_file(new_annot_df_part, file_path.with_stem(
                        f"{file_path.stem}_{ref_seq}_{ANNOT_WITH_CODING_METHYLATIONS_FILE_NAME}.csv"))

            else:
                # Store all results
                if write_all:
                    write_df_to_file(coding_df, file_path.with_stem(
                        f"{file_path.stem}_{CODING_METHYLATIONS_FILE_NAME}.csv"))
                    write_df_to_file(non_coding_df, file_path.with_stem(
                        f"{file_path.stem}_{NON_CODING_METHYLATIONS_FILE_NAME}.csv"))

                # Write always because of panmethylome
                write_df_to_file(new_annot_df, file_path.with_stem(
                    f"{file_path.stem}_{ANNOT_WITH_CODING_METHYLATIONS_FILE_NAME}.csv"))

        print("Methylome mining done.")
        return new_annot_df


def get_pangenome_genes(roary_df, annot_df, genome_name_abbr, work_dir, annot_file):
    """
    Get list of pangenome genes for one genome, verify genes' names in extended annotation and if needed,
    use new annotation (e.g. fixed annotation from Roary)

    :param pd.DataFrame roary_df: Results of pangenome analysis by Roary.
    :param pd.DataFrame annot_df: Extended annotation with added base modification positions.
    :param str genome_name_abbr: Prefix or genome name used to distinguish samples in pangenome analysis.
    :param Path work_dir: Path to directory for (Pan)MethylomeMiner outputs.
    :param Path annot_file: Path to a file with genome annotation in 'gff' (v3) or 'gbk' file format.
    :return: (new)valid_genes, (new)annot_df

        - (new)valid_genes - a pd.Series with names of genes included in pangenome analysis for anqlyzed genome.
        - (new)annot_df - a pd.DataFrame with extended annotation with added base modification positions.
    """
    valid_genes = roary_df[pd.notna(roary_df) & roary_df.isin(annot_df["gene_id"])]

    if len(valid_genes) / len(annot_df["gene_id"]) < 0.9:
        print(
            f"Genes' annotations are not matching between input annotation file and Roary file. "
            f"Running MethylomeMiner for new annotation file {annot_file}."
        )

        # Load new annotation file
        new_annot_df = parse_annotation_file(annot_file)

        # Check intersection of Roary file and new annotation
        new_valid_genes = roary_df[pd.notna(roary_df) & roary_df.isin(new_annot_df["gene_id"])]

        if len(new_valid_genes) / len(new_annot_df["gene_id"]) > 0.9:
            # Load filtered bedMethyl file
            filtered_bed_files_options = [
                Path(work_dir, f"{genome_name_abbr}_{FILTERED_BED_FILE_NAME}.{extension}")
                for extension in ("json", "tsv", "csv", "bed")
            ]
            available_files = [file for file in filtered_bed_files_options if file.exists()]
            filtered_bed_df = parse_bed_file(available_files[0])

            # Sort methylation according to new annotation
            coding_df, non_coding_df, new_annot_df = sort_methylations(filtered_bed_df, new_annot_df)
            write_df_to_file(new_annot_df,
                             Path(work_dir, f"{genome_name_abbr}_{ANNOT_WITH_CODING_METHYLATIONS_FILE_NAME}.csv"))

            return new_valid_genes, new_annot_df

        else:
            print(f"Genes' annotations are still not matching between Roary and annotation file {annot_file}.")
            return None, None
    else:
        return valid_genes, annot_df


def get_panmethylome(roary_output_file, files_pairs, work_dir, matrix_values):
    """
    Create panmethylome for all base modification types

    :param Path roary_output_file: Path to output file from Roary tool named 'gene_presence_absence.csv'.
    :param Path work_dir: Path to directory for (Pan)MethylomeMiner outputs.
    :param dict files_pairs: Prefix oriented dictionary for pairs of Paths to bedMethyl file and genome annotation file.
    :param str matrix_values: Type of values in the output panmethylome matrix. Options: 'presence': '0' value for
         no detected base modifications, '1' value for detected base modification, 'positions': a list of exact 
         locations of base modifications within a pangenome genes. Default is 'presence' option.
    :rtype dict:
    :return: Dictionary with keys corresponding to found modification types and values are DataFrames, where rows
         are pangenome genes according to Roary and columns are analyzed genomes. Values in the DataFrame are
         as requested through `matrix_values` parameter.
    """
    # Load results from Roary
    roary_df = pd.read_csv(roary_output_file)

    # Prepare output dictionary with DataFrames 
    panmethylomes = {methylation: roary_df["Gene"].to_frame() for methylation in METHYLATIONS_KEY.values()}
    
    # Fill DataFrames based on results from MethylomeMiner, more precisely files with extended annotation
    for genome_name_abbr in files_pairs:

        annot_df = pd.read_csv(Path(work_dir, f"{genome_name_abbr}_{ANNOT_WITH_CODING_METHYLATIONS_FILE_NAME}.csv"))
        methylation_types = annot_df.iloc[:, 7:].columns

        # Get column from Roary file that corresponds to current extended annotation DataFrame
        col_names = list(roary_df.filter(regex=genome_name_abbr).columns)
        if len(col_names) != 1:
            print("Invalid column names in Roary file.")
            return None

        # Find pangenome genes in the extended annotation
        valid_genes, annot_df = get_pangenome_genes(roary_df[col_names[0]], annot_df, genome_name_abbr, work_dir,
                                                    files_pairs[genome_name_abbr]["annot_file"])
        if valid_genes is None:
            return None

        # For each extended annotation file, add found modifications to correct pangenome gene
        for methylation, panmethylome in panmethylomes.items():

            panmethylome[genome_name_abbr] = pd.Series()

            if methylation in methylation_types:

                # Add found base modifications' positions to pangenome genes
                panmethylome.loc[valid_genes.index, genome_name_abbr] = valid_genes.map(
                    annot_df.set_index("gene_id")[methylation])

                # Change lists of positions to '0' and '1' values, if requested
                if matrix_values == "presence":
                    panmethylome.loc[panmethylome[genome_name_abbr].astype(str) == "[]", genome_name_abbr] = 0
                    panmethylome.loc[(panmethylome[genome_name_abbr] != 0) & ~panmethylome[
                        genome_name_abbr].isna(), genome_name_abbr] = 1
                elif matrix_values == "positions":
                    continue
                else:
                    print("Invalid matrix values requested. Only 'presence' and 'positions' options are possible.")
                    return None

    return panmethylomes


def get_panmethylations_stats(work_dir, prefixes):
    """
    Get counts statistics per methylation types and per reference sequence for all genomes in pangenome analysis

    :param Path work_dir: Path to directory for (Pan)MethylomeMiner outputs.
    :param list prefixes: Genomes names derived from file names' prefixes
    """

    stats = []
    ref_seq_stats = []
    for prefix in prefixes:
        stats.append(pd.read_csv(Path(work_dir, f"{prefix}_{STATS_FILE_NAME}.csv")))
        ref_seq_stats.append(pd.read_csv(Path(work_dir, f"{prefix}_{STATS_REF_SEQS_FILE_NAME}.csv")))
    stats_df = pd.concat(stats)
    ref_seq_stats_df = pd.concat(ref_seq_stats)

    write_df_to_file(stats_df, Path(work_dir, f"panmethylome_statistics.csv"))
    write_df_to_file(ref_seq_stats_df, Path(work_dir, f"panmethylome_statistics_per_ref.csv"))


def _mine_panmethylations(input_bed_dir, input_annot_dir, roary_file, min_coverage, min_percent_modified,
                          matrix_values, work_dir, write_all_results, heatmap_type, heatmap_file_format,
                          heatmap_min_percent_presence):
    """
    Create panmethylome from bedMethyl files, genome annotation and Roary output.

    :param Path input_bed_dir: Path to a directory with bedMethyl files.
    :param Path input_annot_dir: Path to a directory with genome annotations in 'gff' (v3) or 'gbk' file format.
    :param Path roary_file: Path to output file from Roary tool named 'gene_presence_absence.csv'.
    :param int min_coverage: An integer value of minimum coverage for modified position to be kept.
    :param float min_percent_modified: Minimum required percentage of reads supporting base modification. Default: 90
    :param str matrix_values: Type of values in the output panmethylome matrix. Options: 'presence': '0' value for
         no detected base modifications, '1' value for detected base modification, 'positions': a list of exact 
         locations of base modifications within a pangenome gene. Default is 'presence' option.
    :param Path work_dir: Path to directory for (Pan)MethylomeMiner outputs. Default: MethylomeMiner_output
    :param bool write_all_results: Write all results from (Pan)MethylomeMiner to files. Default: False
    :param str heatmap_type: Choose which heatmap(s) should be created (compact, full, both or none).
        Default: None - no heatmap is created.
    :param str heatmap_file_format: File format of heatmap(s) to be created. Default: 'png'.
    :param float heatmap_min_percent_presence: Minimum required percentage of genomes where are genes present.
    """
    print("Panmethylome mining started.")

    # For each input bed file find matching annotation file based on the prefix in the name of the files
    files_pairs = pair_bed_and_annot_files(input_bed_dir, input_annot_dir)
    if files_pairs is None:
        print("No matching bedMethyl and annotation files found.")

    for prefix, files in files_pairs.items():
        miner_file_exists = {
            "coding": Path(work_dir, f"{prefix}_{CODING_METHYLATIONS_FILE_NAME}.csv").exists(),
            "non_coding": Path(work_dir, f"{prefix}_{NON_CODING_METHYLATIONS_FILE_NAME}.csv").exists(),
            "annot_coding": Path(work_dir, f"{prefix}_{ANNOT_WITH_CODING_METHYLATIONS_FILE_NAME}.csv").exists(),
            "filtered_bed": any(
                [Path(work_dir, f"{prefix}_{FILTERED_BED_FILE_NAME}.{extension}").exists()
                 for extension in ("json", "tsv", "csv", "bed")]
            ),
            "stats": Path(work_dir, f"{prefix}_{STATS_FILE_NAME}.csv").exists(),
            "stats_ref":  Path(work_dir, f"{prefix}_{STATS_REF_SEQS_FILE_NAME}.csv").exists(),
        }
        # Check if it is necessary to run MethylomeMiner
        if ((write_all_results and not all(miner_file_exists.values())) or
                (not write_all_results and (not miner_file_exists["annot_coding"] or not miner_file_exists["stats"] or
                                            not miner_file_exists["stats_ref"]))):
            # Before running MethylomeMiner, check is minimum coverage was set, if not calculate median coverage
            # and use it for filtration
            if min_coverage is None:
                min_coverage = calculate_median_coverage(input_bed_dir)
                print(
                    f"Minimum coverage value not provided. Calculated median coverage {min_coverage} used instead.")

            # Run MethylomeMiner
            _mine_methylations(
                input_bed_file=files["bed_file"],
                input_annot_file=files["annot_file"],
                input_bed_dir=input_bed_dir,
                min_coverage=min_coverage,
                min_percent_modified=min_percent_modified,
                work_dir=work_dir,
                file_name=prefix,
                write_filtered_bed=write_all_results,
                filtered_bed_format="bed",
                write_all=write_all_results,
                split_by_reference=False,
            )

    # Create panmethylomes for all present base modifications
    panmethylomes = get_panmethylome(roary_file, files_pairs, work_dir, matrix_values)

    if panmethylomes is not None:
        # For each base modification save individual output files
        for methylation, panmethylome in panmethylomes.items():
            if not panmethylome.iloc[:, 1:].isnull().values.all():
                # Panmethylome matrix
                write_df_to_file(panmethylome, Path(work_dir, f"{methylation}_panmethylome_{matrix_values}.csv"))
                # Heatmap(s)
                if heatmap_type and matrix_values == "presence":
                    generate_heatmaps(panmethylome, methylation, work_dir, type=heatmap_type,
                                      file_format=heatmap_file_format,
                                      min_percent_presence=heatmap_min_percent_presence)

        # Get panmethylome statistics
        get_panmethylations_stats(work_dir, files_pairs.keys())

    print("Panmethylome mining done.")


def generate_heatmaps(panmethylome_df, methylation, work_dir, type, file_format, min_percent_presence):
    """
    Filter panmethylome and create heatmap(s).

    :param pd.DataFrame panmethylome_df: Panmethylome matrix with 'presence' values.
    :param str methylation: Name of the methylation for which panmethylome was created.
    :param Path work_dir: Path to directory for (Pan)MethylomeMiner outputs.
    :param str type: Choose which heatmap(s) of panmethylome should be created.
        Options:
         'compact': heatmap of panmethylome in lower resolution without genes' names,
         'full': heatmap of panmethylome in full resolution with genes' names,
         'both': both 'compact' and 'full' heatmaps are created.
    :param str file_format: File format of heatmap(s) to be created.
    :param float min_percent_presence: Minimum required percentage of genomes where are genes present.
    """
    panmethylome_df.index = panmethylome_df["Gene"]
    panmethylome_df = panmethylome_df.iloc[:, 1:].astype(float)

    min_presence = int(np.floor(panmethylome_df.shape[1] * (min_percent_presence / 100)))

    filtered_panmethylome_df = panmethylome_df[panmethylome_df.notna().sum(axis=1) >= min_presence]

    if filtered_panmethylome_df.empty:
        print(f"No genes passed the presence threshold ({min_percent_presence}%). Heatmaps will not be generated.")
        return

    if type == "compact" or type == "full":
        create_heatmap(filtered_panmethylome_df, methylation, work_dir, type, file_format, min_percent_presence)
    else:
        create_heatmap(filtered_panmethylome_df, methylation, work_dir, "compact", file_format, min_percent_presence)
        create_heatmap(filtered_panmethylome_df, methylation, work_dir, "full", file_format, min_percent_presence)


def create_heatmap(filtered_panmethylome_df, methylation, work_dir, type, file_format, min_percent_presence):
    """
    Create panmethylome heatmap.

    :param pd.DataFrame filtered_panmethylome_df: Filtered panmethylome matrix based on minimum percent of presence.
    :param str methylation: Name of the methylation for which panmethylome was created.
    :param Path work_dir: Path to directory for (Pan)MethylomeMiner outputs.
    :param type: Choose which heatmap(s) of panmethylome should be created.
        Options:
         'compact': heatmap of panmethylome in lower resolution without genes' names,
         'full': heatmap of panmethylome in full resolution with genes' names,
         'both': both 'compact' and 'full' heatmaps are created.
    :param str file_format: File format of heatmap(s) to be created.
    :param float min_percent_presence: Minimum required percentage of genomes where are genes present.
    """
    sys.setrecursionlimit(5000)

    blue = (0/255, 114/255, 178/255)
    yellow = (253/255, 179/255, 56/255)
    gray = (204/255, 204/255, 204/255)
    cmap = ListedColormap([blue, yellow])

    if type == "compact":
        bbox_anchor = (1.05, 0.65)
        figsize = (20, 20)
        dpi = 150
        show_y_axis = False
    else:
        bbox_anchor = (1.2, 0.65)
        figsize = (max(10, filtered_panmethylome_df.shape[1] * 0.4), max(10, filtered_panmethylome_df.shape[0] * 0.4))
        dpi = 300
        show_y_axis = True

    plt.ioff()
    sns.set(font_scale=1.2)

    g = sns.clustermap(filtered_panmethylome_df,
                       metric="hamming",
                       method="average",
                       cmap=cmap,
                       figsize=figsize,
                       yticklabels=show_y_axis,
                       xticklabels=True,
                       cbar_kws={"label": "Methylation"},
                       mask=filtered_panmethylome_df.isna(),
                       dendrogram_ratio=(0.1, 0.05))

    g.cax.set_visible(False)
    g.fig.subplots_adjust(left=0.05, right=0.85)

    legend_elements = [
        mpatches.Patch(color=blue, label="non-methylated"),
        mpatches.Patch(color=yellow, label="methylated")
    ]
    if min_percent_presence < 100:
        legend_elements.append(mpatches.Patch(color=gray, label="missing gene"))

    g.ax_heatmap.legend(handles=legend_elements,
                        loc="center left",
                        bbox_to_anchor=bbox_anchor,
                        frameon=False,
                        fontsize=14)

    g.ax_heatmap.set_xlabel("Genome", fontsize=14, labelpad=10)

    gene_count = filtered_panmethylome_df.shape[0]
    ylabel_text = f"Genes that are presented at {min_percent_presence}% of the analyzed genomes (n={gene_count})"
    g.ax_heatmap.set_ylabel(ylabel_text, fontsize=14, labelpad=20)

    if not show_y_axis:
        g.ax_heatmap.set_yticks([])
        g.ax_heatmap.set_yticklabels([])

    output_path = Path(work_dir, f"{methylation}_panmethylome_presence_{type}_heatmap.{file_format}")
    g.savefig(output_path, dpi=dpi, bbox_inches="tight")
