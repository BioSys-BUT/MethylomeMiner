"""
MethylomeMiner module contains functions to run the tool to filter and sort methylations and core methylations
"""
from pathlib import Path

import click

from .loading import parse_annotation_file, pair_bed_and_annot_files
from .backend import sort_methylations, write_df_to_file, get_core_methylome
from .preprocessing import get_list_of_methylated_positions, calculate_med_coverage


class Mutex(click.Option):
    """
    Edited click Option class to make mutually exclusive options work
    """

    def __init__(self, *args, **kwargs):
        self.not_required_if = kwargs.pop("not_required_if")

        assert self.not_required_if, "'not_required_if' parameter required"
        kwargs["help"] = (kwargs.get("help", "") + "Option is mutually exclusive with " + ", ".join(
            self.not_required_if) + ".").strip()
        super(Mutex, self).__init__(*args, **kwargs)

    def handle_parse_result(self, ctx, opts, args):
        current_opt = self.name in opts
        for mutex_opt in self.not_required_if:
            if mutex_opt in opts:
                if current_opt:
                    raise click.UsageError(
                        "Illegal usage: '" + str(self.name) + "' is mutually exclusive with '" + str(mutex_opt) + "'.")
                else:
                    self.prompt = None
        return super(Mutex, self).handle_parse_result(ctx, opts, args)


@click.command()
@click.option(
    "--input_bed_file",
    required=True,
    type=click.Path(
        exists=True, file_okay=True, dir_okay=False,
        readable=True, resolve_path=True),
    help="Path to a bedMethyl file with information about modified and unmodified bases.",
)
@click.option(
    "--input_annot_file",
    required=True,
    type=click.Path(
        exists=True, file_okay=True, dir_okay=False,
        readable=True, resolve_path=True),
    help="Path to a file with genome annotation in *.gff (v3) or *.gbk file format.",
)
@click.option(
    "--input_bed_dir",
    required=True, cls=Mutex, not_required_if=["min_coverage"],
    type=click.Path(
        exists=True, dir_okay=True,
        readable=True, resolve_path=True),
    help="Path to a directory with bedMethyl files.\n"
         "This option is mutually exclusive with 'minimum_coverage'.\n"
         "This folder will be used to calculated median coverage from all present bedMethyl files.",
)
@click.option(
    "--min_coverage",
    required=True, cls=Mutex, not_required_if=["input_bed_dir"],
    default=None,
    type=int,
    help="An integer value of minimum coverage for modified position to be kept.\n"
         "This option is mutually exclusive with 'input_bed_dir'.",
)
@click.option(
    "--min_percent_modified",
    required=False,
    default=90,
    type=click.FloatRange(0, 100),
    help="A minimum percent of modified base occurrence.\n"
         "A float value between 0 and 100 inclusive. Default value is 90 %.",
)
@click.option(
    "--work_dir",
    required=False,
    default="MethylomeMiner_output",
    type=click.Path(
        dir_okay=True, writable=True,
        resolve_path=True, path_type=Path),
    help="Path to directory for MethylomeMiner outputs.\n"
         "If not provided, 'MethylomeMiner_output' folder will be created in the current working directory.",
)
@click.option(
    "--file_name",
    required=False,
    default=None,
    help="Custom name for MethylomeMiner outputs.",
)
def mine_methylations(bed_file, annot_file, bed_dir, min_coverage, min_percent_modified, work_dir, file_name):
    """
    Filter modified bases stored in bedMethyl file and sort them according to annotation into coding and non-coding.

    Filtration is performed based on requested coverage and the percent of modified bases.
    Sorting is conducted based on provided annotation of the genome. Modified bases are sorted into coding
    (modification is within coding region) and non-coding (modification is in intergenic region) groups.


    :param str bed_file: Path to a bedMethyl file with information about modified and unmodified bases.
    :param str annot_file: Path to a file with genome annotation in *.gff (v3) or *.gbk file format.
    :param str bed_dir: Path to a directory with bedMethyl files.
    :param int min_coverage: An integer value of minimum coverage for modified position to be kept.
    :param float min_percent_modified: A minimum percent of modified base occurrence. Default: 90
    :param str work_dir: Path to directory for MethylomeMiner outputs. Default: MethylomeMiner_output
    :param str file_name: Custom name for MethylomeMiner outputs.
    """

    bed_df = get_list_of_methylated_positions(
        bed_file,
        bed_dir,
        write_to_file=False,
        min_coverage=min_coverage,
        min_percent_modified=min_percent_modified
    )

    if bed_df is not None:
        annot_df = parse_annotation_file(Path(annot_file))
        coding_df, non_coding_df, new_annot_df = sort_methylations(bed_df, annot_df)

        if file_name is None:
            file_name = bed_file.stem.split("_")[0].split(".")[0]
        if not work_dir.exists():
            work_dir.mkdir(parents=True, exist_ok=True)
        write_df_to_file(coding_df, Path(work_dir, file_name + "_coding.csv"))
        write_df_to_file(non_coding_df, Path(work_dir, file_name + "_non_coding.csv"))
        write_df_to_file(new_annot_df, Path(work_dir, file_name + "_annot_with_methylations.csv"))

    print("Methylome mining done.")


@click.command()
@click.option(
    "--input_bed_dir",
    required=True,
    type=click.Path(
        exists=True, dir_okay=True,
        readable=True, resolve_path=True),
    help="Path to a directory with bedMethyl files.",
)
@click.option(
    "--input_annot_dir",
    required=True,
    type=click.Path(
        exists=True, dir_okay=True,
        readable=True, resolve_path=True),
    help="Path to a directory with genome annotations in *.gff (v3) or *.gbk file format.",
)
@click.option(
    "--roary_file",
    required=True,
    type=click.Path(
        exists=True, file_okay=True, dir_okay=False,
        readable=True, resolve_path=True),
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
def mine_core_methylome(bed_dir, annot_dir, roary_file, work_dir, min_coverage, min_percent_modified, matrix_values,
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

    files_pairs = pair_bed_and_annot_files(bed_dir, annot_dir)

    for prefix, files in files_pairs.items():
        if (write_all and (not Path(work_dir, prefix + "_coding.csv").exists() or
                           not Path(work_dir, prefix + "_non_coding.csv").exists() or
                           not Path(work_dir, prefix + "_all_annot_with_methylations.csv").exists()
        ) or (not write_all and not Path(work_dir, prefix + "_all_annot_with_methylations.csv").exists())):
                bed_file = files["bed_file"]
                annot_file = files["annot_file"]
                if min_coverage is None:
                    min_coverage = calculate_med_coverage(bed_dir)
                bed_df = get_list_of_methylated_positions(
                    bed_file,
                    bed_dir,
                    write_to_file=False,
                    min_coverage=min_coverage,
                    min_percent_modified=min_percent_modified
                )

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
