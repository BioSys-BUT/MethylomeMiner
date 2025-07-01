"""
MethylomeMiner module contains functions to run the tool to filter and sort methylations and core methylations
"""
from pathlib import Path

import click

from .loading import parse_annotation_file, pair_bed_and_annot_files
from .backend import sort_methylations, write_df_to_file, get_core_methylome, filter_methylations, calculate_median_coverage, write_bed_file


@click.command()
@click.option(
    "--input_bed_dir",
    required=True,
    type=click.Path(
        exists=True, dir_okay=True,
        readable=True, resolve_path=True, path_type=Path),
    help="Path to a directory with bedMethyl files.",
)
@click.option(
    "--input_annot_dir",
    required=True,
    type=click.Path(
        exists=True, dir_okay=True,
        readable=True, resolve_path=True, path_type=Path),
    help="Path to a directory with genome annotations in *.gff (v3) or *.gbk file format.",
)
@click.option(
    "--roary_file",
    required=True,
    type=click.Path(
        exists=True, file_okay=True, dir_okay=False,
        readable=True, resolve_path=True, path_type=Path),
    help="Path to output file from Roary tool named 'gene_presence_absence.csv'.",
)
@click.option(
    "--work_dir",
    required=False,
    default="MethylomeMiner_output",
    type=click.Path(
        dir_okay=True, writable=True,
        resolve_path=True, path_type=Path),
    help="Path to directory for MethylomeMiner and CoreMethylomeMiner outputs.\n"
         "If not provided, 'MethylomeMiner_output' folder will be created in the current working directory.",
)
@click.option(
    "--min_coverage",
    required=False,
    default=None,
    type=int,
    help="An integer value of minimum coverage for modified position to be kept.\n"
         "This option is used only if MethylomeMiner results was not used for input bedMethyl files"
         "and the results are not present in the 'work_dir' folder.\n"
         "If value is not provided, median coverage will be calculated from files in the 'input_bed_dir' folder.\n",
)
@click.option(
    "--min_percent_modified",
    required=False,
    default=90,
    type=click.FloatRange(0, 100),
    help="A minimum percent of modified base occurrence.\n"
         "This option is used only if MethylomeMiner results was not used for input bedMethyl files"
         "and the results are not present in the 'work_dir' folder.\n"
         "A float value between 0 and 100 inclusive. Default value is 90 %.",
)
@click.option(
    "--matrix_values",
    required=False,
    default="presence",
    help="Type of values in the output core methylome matrix.\n"
         "Options:\n"
         "'presence': '0' value for no detected base modifications, '1' value for detected base modification,\n"
         "'positions': a list of exact location of base modification within a core gene.\n"
         "Default is 'presence' option.",
)
@click.option(
    "--write_all",
    required=False,
    default=False,
    help="Write all results from MethylomeMiner to files: methylations sorted to coding and non-coding groups, "
         "genome annotation with found methylations in coding regions.\n"
         "Default is 'False', set to 'True' to save all results to files.",
)
def mine_core_methylations(input_bed_dir, input_annot_dir, roary_file, work_dir, min_coverage, min_percent_modified, matrix_values,
                        write_all):
    """
    Create core methylome from bedMethyl files, genome annotation and Roary output.

    :param str bed_dir: Path to a directory with bedMethyl files.
    :param str annot_dir: Path to a directory with genome annotations in *.gff (v3) or *.gbk file format.
    :param str roary_file: Path to output file from Roary tool named 'gene_presence_absence.csv'.
    :param str work_dir: Path to directory for (Core)MethylomeMiner outputs. Default: MethylomeMiner_output
    :param int min_coverage: An integer value of minimum coverage for modified position to be kept.
    :param float min_percent_modified: A minimum percent of modified base occurrence. Default: 90
    :param str matrix_values: Type of values in the output core methylome matrix. Default: presence
    :param bool write_all: Write all results from MethylomeMiner to files. Default: False
    """

    files_pairs = pair_bed_and_annot_files(input_bed_dir, input_annot_dir)

    for prefix, files in files_pairs.items():
        if (write_all and (not Path(work_dir, prefix + "_coding.csv").exists() or
                           not Path(work_dir, prefix + "_non_coding.csv").exists() or
                           not Path(work_dir, prefix + "_all_annot_with_methylations.csv").exists()
        ) or (not write_all and not Path(work_dir, prefix + "_all_annot_with_methylations.csv").exists())):
                bed_file = files["bed_file"]
                annot_file = files["annot_file"]
                if min_coverage is None:
                    min_coverage = calculate_median_coverage(input_bed_dir)
                    print(
                        f"Minimum coverage value not provided. Calculated median coverage {min_coverage} used instead.")
                bed_df = filter_methylations(
                    bed_file,
                    input_bed_dir,
                    min_coverage=min_coverage,
                    min_percent_modified=min_percent_modified
                )
                if write_all:
                    # save filtered bed files + make statistics
                    write_bed_file(bed_df, Path(work_dir, prefix), file_format="csv")

                if bed_df is not None:
                    annot_df = parse_annotation_file(Path(annot_file))
                    coding_df, non_coding_df, new_annot_df = sort_methylations(bed_df, annot_df, all_annot=True)
                    write_df_to_file(new_annot_df, Path(work_dir, prefix + "_all_annot_with_methylations.csv"))
                    if write_all:
                        write_df_to_file(coding_df, Path(work_dir, prefix + "_coding.csv"))
                        write_df_to_file(non_coding_df, Path(work_dir, prefix + "_non_coding.csv"))

    core_methylomes = get_core_methylome(roary_file, miner_output_dir=work_dir, matrix_values=matrix_values)

    for methylation, methylome in core_methylomes.items():
        if not methylome.iloc[:, 1:].isnull().values.all():
            write_df_to_file(methylome, Path(work_dir, methylation + "_core_methylome_" + matrix_values +".csv"))

    print("Core methylome mining done.")
