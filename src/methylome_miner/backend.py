"""
Pre-processing module filter raw bedMethyle files and store filtered files in various file formats
"""

from pathlib import Path
import statistics

import pandas as pd

from .loading import METHYLATIONS_KEY, parse_bed_file


def calculate_median_coverage(bed_dir_path):
    """
    Calculate median coverage from all available bedMethyl files

    :param Path bed_dir_path: Path to a directory with bedMethyl files.
    :rtype: int
    :return: An integer value of median coverage across all bedMethyl files.
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

    :param bed_file_path: Path to a bedMethyl file with information about modified and unmodified bases.
    :param bed_dir_path: Path to a directory with bedMethyl files.
    :param min_coverage: An integer value of minimum coverage for modified position to be kept.
    :param min_percent_modified: A minimum percent of modified base occurrence.
    :rtype: pd.DataFrame
    :return: Filtered bedMethyl table
    """

    bed_df = parse_bed_file(bed_file_path)

    if bed_df is None:
        return None

    if min_coverage is None:
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
    Write pandas DataFrame to file of requested file format (options are '.json', '.csv', '.tsv', '.bed').

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
            print("Unsupported file format. Choose from: '.json', '.csv', '.tsv' or '.bed'.")


def to_bed(bed_df, file_path):
    """
    Save pandas DataFrame as bedMethyl file.

    File structure is the same as the output bedMethyl file from Modkit tool.
    The bedMethyl table contains eighteen columns, but the bedMethyl file has only nine columns that are
    tab-delimited. The remaining columns are stored within the last column and they are space-delimited.

    :param pd.DataFrame bed_df: BedMethyl table to be saved.
    :param Path file_path: Path to new bedMethyl file.
    """
    tab_part = bed_df.iloc[:, :9].astype(str).agg("\t".join, axis=1)
    space_part = bed_df.iloc[:, 9:].astype(str).agg(" ".join, axis=1)

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
    :param methylations_df: Table of coding base modifications.
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

    TODO: Overlapping genes (remove next_methylation_search_idx?) and ?

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

    TODO: Can be sort_methylations and sort_coding_non_coding_methylations merged together?
    :param pd.DataFrame bed_df: Table of base modifications.
    :param annot_df: Table with annotation.
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
    methylation_types = ["m" + cat for cat in bed_df["modified_base_code"].astype("category").cat.categories]

    # Go through each group and sort modifications
    for group, bed_df_group in bed_df_grouped:
        # Add new columns named after base modifications to annotation
        annot_df_group = annot_df_grouped.get_group(group)
        annot_df_group = pd.concat([annot_df_group, pd.DataFrame(columns=methylation_types, dtype="object")])

        coding_df, non_coding_df, new_annot_df = sort_coding_non_coding_methylations(bed_df_group, annot_df_group)

        all_coding.append(coding_df)
        all_non_coding.append(non_coding_df)
        new_annots.append(new_annot_df)

    all_coding_df = pd.concat(all_coding).sort_values(by=["reference_seq", "start_index"])
    all_non_coding_df = pd.concat(all_non_coding).sort_values(by=["reference_seq", "start_index"])
    new_annot_df = pd.concat(new_annots).sort_values(by=["record_id", "start"])

    return all_coding_df, all_non_coding_df, new_annot_df


def get_core_methylome(roary_output_file, miner_output_dir="MethylomeMiner_output", matrix_values="presence"):
    roary_df = pd.read_csv(roary_output_file)

    core_methylomes = {methylation: roary_df["Gene"].to_frame() for methylation in METHYLATIONS_KEY.values()}

    for coding_df_file in list(Path(miner_output_dir).glob("*_all_annot_with_methylations.csv")):
        coding_df = pd.read_csv(coding_df_file, delimiter=";")
        methylation_types = coding_df.iloc[:, 7:].columns

        genome_name_abbr = coding_df_file.stem.split("_")[0]

        for methylation, core_methylome in core_methylomes.items():
            if methylation in methylation_types:
                core_methylome[genome_name_abbr] = pd.Series()

                genes = roary_df[f"{genome_name_abbr}_genome"]
                valid_genes = genes[pd.notna(genes) & genes.isin(coding_df["gene_id"])]
                core_methylome.loc[valid_genes.index, genome_name_abbr] = valid_genes.map(coding_df.set_index("gene_id")[methylation])

                if matrix_values == "presence":
                    core_methylome.loc[core_methylome[genome_name_abbr] == "[]", genome_name_abbr] = 0
                    core_methylome.loc[(core_methylome[genome_name_abbr] != 0) & ~core_methylome[genome_name_abbr].isna(), genome_name_abbr] = 1
                elif matrix_values == "positions":
                    continue
                else:
                    return None
            else:
                core_methylome[genome_name_abbr] = pd.Series()
    return core_methylomes
