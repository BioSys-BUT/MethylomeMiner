import json

from pathlib import Path

import pandas as pd

from loading import parse_annotation_file


def get_list_of_methylated_genes(data_frame, strand, list_of_methylated_positions):

    rows = pd.DataFrame()
    for contig_info, values in list_of_methylated_positions.items():
        contig_name, methylation_type, meth_strand = contig_info
        for value in values:
            row = data_frame.loc[(data_frame["contig"] == contig_name) & (data_frame["start"] <= value) & (value <= data_frame["end"])]
            if not row.empty:
                rows = pd.concat([rows, pd.DataFrame(
                      {"contig": contig_name, "gene_id": row["gene_id"], "start": row["start"], "end": row["end"], "position": value,
                      "methylation_type": methylation_type, "strand": row["strand"], "product": row["product"]})])

    if not rows.empty:
        df = rows
        df_grouped = df.groupby(['contig', "gene_id", 'start', 'end', 'methylation_type', 'strand', 'product']).agg(
            count=('position', 'size'),
            positions=('position', lambda x: ','.join(map(str, sorted(x))))
        ).reset_index()

        existing_positions = set(zip(df['position'], df['methylation_type'], df['contig'], df['strand']))
    else:
        df_grouped = pd.DataFrame(columns=['contig', 'gene_id', 'start', 'end', 'methylation_type', 'strand', 'product', 'count', 'positions'])
        existing_positions = set()

    # Převést list_of_methylated_positions na množinu
    methylated_positions_set = set(
        (value, key[1], key[0], key[2])
        for key, values in list_of_methylated_positions.items()
        for value in values
    )

    # Najít rozdíly mezi existing_positions a methylated_positions_set
    missing_in_existing = methylated_positions_set - existing_positions

    grouped = data_frame.groupby('contig')

    # vytvorim skupiny dle jednotlivych contigu
    grouped_missing = {}
    for value, key, contig, orientation in missing_in_existing:
        if contig not in grouped_missing:
            grouped_missing[contig] = []
        grouped_missing[contig].append((value, key, contig, orientation))

    # v ramci contigu seradim od nejmensiho po nejvetsi
    for contig in grouped_missing:
        grouped_missing[contig] = sorted(grouped_missing[contig], key=lambda x: x[0])

    previous_row = None
    gene_near = []
    for contig_name, group in grouped:
        if contig_name not in grouped_missing:
            continue
        for value, key, contig, orientation in grouped_missing[contig_name]:
            for row in group.itertuples():
                if value < row.start and row.Index == group.index[0]:
                    gene_near.append((contig, key, value, 'before first gene', 0, 0, "-", row.gene_id,
                                      row.start, row.end, row.product, row.strand))
                    break
                elif value < row.start:
                    gene_near.append((contig, key, value, previous_row.gene_id, previous_row.start, previous_row.end, previous_row.product,
                                      row.gene_id, row.start, row.end, row.product, row.strand))
                    break
                elif row.Index == group.index[-1]:
                    gene_near.append((contig, key, value, row.gene_id, row.start, row.end,  row.product,
                                      "behind last gene", 0, 0, "-", row.strand))
                    break
                previous_row = row


    df_missing_methylation = pd.DataFrame(gene_near, columns=[
        "contig", "methylation", "position", "previous_gene_id", "previous_gene_start", "previous_gene_end", "previous_gene_product",
        "next_gene_id", "next_gene_start", "next_gene_end", "next_gene_product", "strand"])

    return df_grouped, df_missing_methylation

def generate_coding_and_noncoding_outputs(df_grouped_leading, df_grouped_complement, df_missing_methylation_leading, df_missing_methylation_complement, strand, output_file_coding, output_file_non_coding):
    if strand == None:
        combined_df_coding = pd.concat([df_grouped_leading, df_grouped_complement])
        combined_df_coding_sorted = combined_df_coding.sort_values(by=['contig', 'start'])
        combined_df_coding_sorted.to_csv(output_file_coding, sep='\t', index=False)

        combined_df_non_coding = pd.concat([df_missing_methylation_leading, df_missing_methylation_complement])
        combined_df_non_coding_sorted = combined_df_non_coding.sort_values(by=['contig', 'position'])
        combined_df_non_coding_sorted.to_csv(output_file_non_coding, sep='\t', index=False)

    elif strand == '+':
        df_grouped_leading.to_csv(output_file_coding, sep='\t', index=False)
        df_missing_methylation_leading.to_csv(output_file_non_coding, sep='\t', index=False)

    elif strand == '-':
        df_grouped_complement.to_csv(output_file_coding, sep='\t', index=False)
        df_missing_methylation_complement.to_csv(output_file_non_coding, sep='\t', index=False)

    else:
        raise ValueError(f"Invalid strand: {strand}.")

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
    # file = Path("D:\OneDrive - VUT\_AZV_Helca\methylome\input_roary_output_files_corrected\KP825_genome.gff")
    # annot = parse_annotation_file(file)
    methylation_dir = Path.cwd()
    # annotation_dir = Path.cwd() / "input_roary_output_files_corrected"
    annotation_dir = Path.cwd() / "gbk_gff_data" / "gff_data_uncorrected"
    paired_files = pair_json_with_annotation(methylation_dir, annotation_dir)
    df_cds_info = []
    for files in paired_files:
        with open(files[1]) as json_file:
            methylated_positions_str = json.load(json_file)
        methylated_positions = {eval(key): value for key, value in methylated_positions_str.items()}
        # print(files[2])
        df_cds_info = parse_annotation_file(files[2])

        df_cds_info_leading = df_cds_info[df_cds_info["strand"] == "+"]
        methylated_positions_leading = {key: values for key, values in methylated_positions.items() if key[2] == '+'}
        df_grouped_leading, df_missing_methylation_leading = get_list_of_methylated_genes(df_cds_info_leading, '+', methylated_positions_leading)

        # df_cds_info_complement = df_cds_info[df_cds_info['strand'] == '-']
        # methylated_positions_complement = {key: values for key, values in methylated_positions.items() if key[2] == '-'}
        # df_grouped_complement, df_missing_methylation_complement = get_list_of_methylated_genes(df_cds_info_complement, '-',
        #                                                                                         methylated_positions_complement)
        #
        # ## choose strand
        # strand = None
        # output_dir = r'/home/user_pool_1/nykrynova/methylation_KP_Dorado/ESCMID_version/MM_output/coding_and_noncoding'
        # os.makedirs(output_dir, exist_ok=True)
        #
        # output_coding = os.path.join(output_dir, f"{files[0]}_coding.csv")
        # output_noncoding = os.path.join(output_dir, f"{files[0]}_noncoding.csv")
        # generate_coding_and_noncoding_outputs(df_grouped_leading, df_grouped_complement, df_missing_methylation_leading,
        #                                       df_missing_methylation_complement, strand, output_coding, output_noncoding)
    print("")