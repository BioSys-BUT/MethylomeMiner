from pathlib import Path

import pandas as pd

from loading import parse_annotation_file

import time


def sort_coding_non_coding_methylations(bed_df, annot_df, all_coding=False):
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
        elif all_coding:
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
    new_annot.columns = new_annot.columns[:7].union([cat[1:] for cat in new_annot.columns[7:]], sort=False)
    if methylations_df is not None:
        for group, methylations_group in methylations_df.groupby("modified_base_code", observed=False):
            new_annot[group] = new_annot[group].astype("object")
            new_annot.at[0, group] = methylations_group["start_index"].values.tolist()
    return new_annot


def sort_methylations(bed_df, annot_df):

    bed_df_grouped = bed_df.groupby(["reference_seq", "strand"], observed=False)
    annot_df_grouped = annot_df.groupby(["record_id", "strand"], observed=False)
    all_coding = []
    all_non_coding = []
    new_annots = []
    methylation_types = ["m" + cat for cat in bed_df["modified_base_code"].astype("category").cat.categories]
    for group, bed_df_group in bed_df_grouped:
        annot_df_group = annot_df_grouped.get_group(group)
        annot_df_group = pd.concat([annot_df_group, pd.DataFrame(columns=methylation_types, dtype="object")])

        coding_df, non_coding_df, new_annot_df = sort_coding_non_coding_methylations(bed_df_group, annot_df_group)

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


def pair_json_with_annotation(json_dir, annotation_dir):

    json_files = list(Path(json_dir).glob("*.json"))
    gff_files = list(Path(annotation_dir).glob("*.gff"))

    json_groups = {file.stem.split("_")[0].split(".")[0]: file for file in json_files}
    gff_groups = {file.stem.split("_")[0].split(".")[0]: file for file in gff_files}

    paired_files = []
    for prefix, json_file in json_groups.items():
        if prefix in gff_groups:
            paired_files.append((prefix, json_file, gff_groups[prefix]))
    return paired_files


if __name__ == "__main__":
    # bed_df1 = get_list_of_methylated_positions(Path(r"input_bed_files/KP825_b53_4mC_5mC_6mA_calls_modifications_whole_run_aligned_sorted_pileup_SUP_qscore.bed"), Path("input_bed_files"))
    bed_df1 = pd.read_json(Path(r"KP825_b53_4mC_5mC_6mA_calls_modifications_whole_run_aligned_sorted_pileup_SUP_qscore_filtered2.json"))
    annot_df1 = parse_annotation_file(Path(r"input_roary_output_files_corrected\KP825_genome.gff"))
    coding_df, non_coding_df, new_annot_coding_df = sort_methylations(bed_df1, annot_df1)
