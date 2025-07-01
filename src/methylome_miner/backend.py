"""
Pre-processing module filter raw bedMethyle files and store filtered files in various file formats
"""

from pathlib import Path
import statistics
import time

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


def write_bed_file(bed_df, file_path, file_format):
    """
    Write pandas DataFrame to file of requested file format.

    :param pd.DataFrame bed_df: BedMethyl table to be saved.
    :param Path file_path: Path to new file.
    :param str file_format: File format of the new file (options are 'json', 'csv', 'tsv', 'bed').
    """
    output_file_path = file_path.with_stem(file_path.stem + f"_filtered.{file_format}")

    match file_format:
        case "json":
            pd.DataFrame.to_json(bed_df, output_file_path)
        case "csv":
            bed_df.to_csv(output_file_path, index=False)
        case "tsv":
            bed_df.to_csv(output_file_path, index=False, sep="\t")
        case "bed":
            to_bed(bed_df, output_file_path)
        case _:
            print("Unsupported file format. Choose from: json, csv, tsv or bed.")


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


def sort_coding_non_coding_methylations(bed_df, annot_df, all_annot=False):
    bed_df = bed_df.reset_index(drop=True)

    coding_parts = []
    non_coding_parts = []
    new_annot_parts = []
    next_methylation_search_idx = 0

    start = time.time()

    for idx, annot in enumerate(annot_df.itertuples(index=False)):

        bed_slice = bed_df.iloc[next_methylation_search_idx:]

        coding_part = (
            bed_slice
            .loc[lambda df: df["start_index"] <= annot.end]
            .loc[lambda df: annot.start <= df["start_index"]]
            .loc[lambda df: df["strand"] == annot.strand]
            .loc[lambda df: df["reference_seq"] == annot.record_id]
        ).copy()

        coding_start_index = coding_part.index.min() if not coding_part.empty else bed_slice.index.min()

        if idx == 0:
            non_coding_part = bed_df.loc[lambda df: df["start_index"] < annot.start].copy()
            if not non_coding_part.empty:
                non_coding_part = add_annot_to_non_coding(non_coding_part, "", annot, "start")
                next_methylation_search_idx = non_coding_part.index.max() + 1
                non_coding_parts.append(non_coding_part)

        elif coding_start_index > next_methylation_search_idx:
            non_coding_part = bed_df.iloc[next_methylation_search_idx:coding_start_index].copy()
            if not non_coding_part.empty:
                non_coding_part = add_annot_to_non_coding(non_coding_part, prev_annot, annot)
                next_methylation_search_idx = non_coding_part.index.max() + 1
                non_coding_parts.append(non_coding_part)

        if idx == annot_df.shape[0] - 1:
            tail_start_idx = coding_part.index.max() + 1 if not coding_part.empty else next_methylation_search_idx
            non_coding_part = bed_df.loc[tail_start_idx:].copy()
            if not non_coding_part.empty:
                non_coding_part = add_annot_to_non_coding(non_coding_part, prev_annot, annot, "end")
                non_coding_parts.append(non_coding_part)

        if not coding_part.empty:
            new_annot_parts.append(add_coding_to_annot(annot, coding_part))
            coding_part = add_annot_to_coding(coding_part, annot)
            coding_parts.append(coding_part)

            next_methylation_search_idx = coding_part.index.max() + 1
        elif all_annot:
            new_annot_parts.append(add_coding_to_annot(annot))

        prev_annot = annot

    coding_df = pd.concat(coding_parts)
    non_coding_df = pd.concat(non_coding_parts)
    new_annot_df = pd.concat(new_annot_parts)

    end = time.time()
    # print(end - start)

    # all_df = pd.concat([coding_df, non_coding_df])
    # all_df = all_df.sort_values(["reference_seq", "start_index"])
    # my_diff = pd.merge(bed_df, all_df, how="outer", indicator=True)

    return coding_df, non_coding_df, new_annot_df


def fill_annot(df, prefix, annot):
    if annot:
        df[f"{prefix}_gene_reference_seq"] = annot.record_id
        df[f"{prefix}_gene_id"] = annot.gene_id
        df[f"{prefix}_gene_start"] = annot.start
        df[f"{prefix}_gene_end"] = annot.end
        df[f"{prefix}_gene_strand"] = annot.strand
        df[f"{prefix}_gene_product"] = annot.product
    else:
        for col in ["reference_seq", "id", "start", "end", "strand", "product"]:
            df[f"{prefix}_gene_{col}"] = ""
    return df


def add_annot_to_non_coding(non_coding_df, prev_annot, next_annot, annot_position="middle"):

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
    coding_df["gene_reference_seq"] = annot.record_id
    coding_df["gene_id"] = annot.gene_id
    coding_df["gene_start"] = annot.start
    coding_df["gene_end"] = annot.end
    coding_df["gene_strand"] = annot.strand
    coding_df["gene_product"] = annot.product
    return coding_df


def add_coding_to_annot(annot, methylations_df=None):
    new_annot = pd.DataFrame([annot])
    new_col_names = new_annot.columns[:7].union([cat[1:] for cat in new_annot.columns[7:]], sort=False)
    new_annot.columns = new_col_names

    for name in new_col_names[7:]:
        new_annot[name] = [[]]

    if methylations_df is not None:
        for group, methylations_group in methylations_df.groupby("modified_base_code", observed=False):
            new_annot.at[0, group] = methylations_group["start_index"].tolist()

    return new_annot


def sort_methylations(bed_df, annot_df, all_annot=False):

    bed_df_grouped = bed_df.groupby(["reference_seq", "strand"], observed=False)
    annot_df_grouped = annot_df.groupby(["record_id", "strand"], observed=False)
    all_coding = []
    all_non_coding = []
    new_annots = []
    methylation_types = ["m" + cat for cat in bed_df["modified_base_code"].astype("category").cat.categories]
    for group, bed_df_group in bed_df_grouped:
        annot_df_group = annot_df_grouped.get_group(group)
        annot_df_group = pd.concat([annot_df_group, pd.DataFrame(columns=methylation_types, dtype="object")])

        coding_df, non_coding_df, new_annot_df = sort_coding_non_coding_methylations(bed_df_group, annot_df_group, all_annot=all_annot)

        all_coding.append(coding_df)
        all_non_coding.append(non_coding_df)
        new_annots.append(new_annot_df)

    all_coding_df = pd.concat(all_coding).sort_values(by=["reference_seq", "start_index"])
    all_non_coding_df = pd.concat(all_non_coding).sort_values(by=["reference_seq", "start_index"])
    new_annot_df = pd.concat(new_annots).sort_values(by=["record_id", "start"])
    # all_df = pd.concat([all_coding_df, all_non_coding_df])

    # all_df = all_df.sort_values(["reference_seq", "start_index"])
    # diff = pd.merge(bed_df, all_df.iloc[:, :18], how="outer", indicator=True)
    return all_coding_df, all_non_coding_df, new_annot_df


def write_df_to_file(methylations_df, file_name, file_type="csv"):
    Path(file_name).parent.mkdir(parents=True, exist_ok=True)

    match file_type:
        case "csv":
            methylations_df.to_csv(file_name, sep=";", index=False)
        case "json":
            methylations_df.to_json(file_name, indent=4)


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


if __name__ == "__main__":
    # bed_df1 = get_list_of_methylated_positions(Path(r"input_bed_files/KP825_b53_4mC_5mC_6mA_calls_modifications_whole_run_aligned_sorted_pileup_SUP_qscore.bed"), Path("input_bed_files"))
    # bed_df1 = pd.read_json(Path(r"KP825_b53_4mC_5mC_6mA_calls_modifications_whole_run_aligned_sorted_pileup_SUP_qscore_filtered2.json"))
    # annot_df1 = parse_annotation_file(Path(r"input_roary_output_files_corrected\KP825_genome.gff"))
    # coding_df, non_coding_df, new_annot_coding_df = sort_methylations(bed_df1, annot_df1)

    print(get_core_methylome(Path("roary_output", "gene_presence_absence.csv"), matrix_values="positions"))

