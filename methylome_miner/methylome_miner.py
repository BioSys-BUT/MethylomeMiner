"""
MethylomeMiner module contains functions to run the tool to detect methylations
"""
from pathlib import Path

import click

from methylome.methylome_miner.loading import parse_annotation_file
from methylome.methylome_miner.backend import sort_methylations, write_df_to_file
from methylome.methylome_miner.preprocessing import get_list_of_methylated_positions


@click.command()
@click.option(
    '--input_bed_file',
    required=True,
    type=click.Path(
        exists=True, file_okay=True, dir_okay=False,
        readable=True, resolve_path=True),
    help='Path to a BED file with methylation data.',
)
@click.option(
    '--input_bed_dir',
    required=True,
    type=click.Path(
        exists=True, dir_okay=True,
        readable=True, resolve_path=True),
    help='Path to a dir with BED files with methylation data.',
)
@click.option(
    '--input_annot_file',
    required=True,
    type=click.Path(
        exists=True, file_okay=True, dir_okay=False,
        readable=True, resolve_path=True),
    help='Path to a file with genome annotation.',
)
@click.option(
    '--output_files_prefix',
    required=True,
    type=click.Path(
        file_okay=True, dir_okay=False, writable=True,
        resolve_path=True, path_type=Path),
    help='Path to files with detected methylations.',
)
@click.option(
    '--min_coverage',
    required=False,
    default=None,
    type=int,
    help='Minimum coverage for methylated position.',
)
@click.option(
    '--min_percent_modified',
    required=False,
    default=90,
    type=click.IntRange(0, 100),
    help='Minimum percent modified.',
)
def mine_methylations(input_bed_file, input_bed_dir, input_annot_file, output_files_prefix, min_coverage, min_percent_modified):

    bed_df = get_list_of_methylated_positions(
        input_bed_file,
        input_bed_dir,
        write_to_file=False,
        min_coverage=min_coverage,
        min_percent_modified=min_percent_modified
    )
    if bed_df is not None:
        annot_df = parse_annotation_file(Path(input_annot_file))
        coding_df, non_coding_df, new_annot_df = sort_methylations(bed_df, annot_df)
        write_df_to_file(coding_df, Path(output_files_prefix.with_stem(output_files_prefix.stem + "_coding.csv")))
        write_df_to_file(non_coding_df, Path(output_files_prefix.with_stem(output_files_prefix.stem + "_non_coding.csv")))
        write_df_to_file(new_annot_df, Path(output_files_prefix.with_stem(output_files_prefix.stem + "_annot_with_methylations.csv")))


def main():
    mine_methylations()


if __name__ == "__main__":
    main()
