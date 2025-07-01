"""
MethylomeMiner module contains functions to run the tool to filter and sort methylations and core methylations
"""
from pathlib import Path

import click

from .loading import parse_annotation_file
from .backend import sort_methylations, write_df_to_file, filter_methylations


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
    type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True, path_type=Path),
    help="Path to a bedMethyl file with information about modified and unmodified bases.",
)
@click.option(
    "--input_annot_file",
    required=True,
    type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True, path_type=Path),
    help="Path to a file with genome annotation in *.gff (v3) or *.gbk file format.",
)
@click.option(
    "--input_bed_dir",
    cls=Mutex, not_required_if=["min_coverage"],
    type=click.Path(exists=True, file_okay=False, dir_okay=True, readable=True, resolve_path=True, path_type=Path),
    help="Path to a directory with bedMethyl files.\n"
         "This option is mutually exclusive with 'min_coverage'.\n"
         "This folder will be used to calculated median coverage from all present bedMethyl files.",
)
@click.option(
    "--min_coverage",
    cls=Mutex, not_required_if=["input_bed_dir"], default=None, type=int,
    help="An integer value of minimum coverage for modified position to be kept.\n"
         "This option is mutually exclusive with 'input_bed_dir'.",
)
@click.option(
    "--min_percent_modified",
    required=False, default=90, type=click.FloatRange(0, 100),
    help="A minimum percent of modified base occurrence.\n"
         "A float value between 0 and 100 inclusive. Default value is 90 %.",
)
@click.option(
    "--work_dir",
    required=False,
    default=Path(Path.cwd(), "MethylomeMiner_output"),
    type=click.Path(file_okay=False, dir_okay=True, writable=True, resolve_path=True, path_type=Path),
    help="Path to directory for MethylomeMiner outputs.\n"
         "If not provided, 'MethylomeMiner_output' folder will be created in the current working directory.",
)
@click.option(
    "--file_name",
    required=False,
    default=None,
    help="Custom name for MethylomeMiner outputs.",
)
@click.option(
    "--write_filtered_bed",
    required=False,
    is_flag=True,
    help="Write filtered bedMethyl file to a new file.",
)
@click.option(
    "--filtered_bed_format",
    required=False,
    default="csv",
    help="Choose filtered bedMethyl file format.\n"
         "Options: 'json', 'csv', 'tsv', 'bed'.\nDefault format is 'csv'.",
)
def mine_methylations(input_bed_file, input_annot_file, input_bed_dir, min_coverage, min_percent_modified,
                      work_dir, file_name, write_filtered_bed, filtered_bed_format):
    """
    Filter modified bases stored in bedMethyl file and sort them according to annotation into coding and non-coding.

    Filtration is performed based on requested coverage and the percent of modified bases.
    Sorting is conducted based on provided annotation of the genome. Modified bases are sorted into coding
    (modification is within coding region) and non-coding (modification is in intergenic region) groups.


    :param Path input_bed_file: Path to a bedMethyl file with information about modified and unmodified bases.
    :param Path input_annot_file: Path to a file with genome annotation in '.gff' (v3) or '.gbk' file format.
    :param Path input_bed_dir: Path to a directory with bedMethyl files.
    :param int min_coverage: An integer value of minimum coverage for modified position to be kept.
    :param float min_percent_modified: A minimum percent of modified base occurrence. Default: 90
    :param Path work_dir: Path to directory for MethylomeMiner outputs. Default: MethylomeMiner_output
    :param str file_name: Custom name for MethylomeMiner outputs.
    :param bool write_filtered_bed: Write filtered bedMethyl file to a new file. Default: False
    :param str filtered_bed_format: File format for filtered bedMethyl file. Default: csv
    """
    print(f"Mining of {input_bed_file} started.")

    # Filter out invalid base modifications
    filtered_bed_df = filter_methylations(
        input_bed_file,
        input_bed_dir,
        min_coverage=min_coverage,
        min_percent_modified=min_percent_modified,
    )

    if filtered_bed_df is not None:
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
            file_path = file_path.with_stem(file_path.stem + f"_filtered.{filtered_bed_format}")
            write_df_to_file(filtered_bed_df, file_path)

        # Parse genome annotation file
        annot_df = parse_annotation_file(input_annot_file)
        if annot_df is not None:

            # According to annotation sort modifications into coding and non-coding groups
            coding_df, non_coding_df, new_annot_df = sort_methylations(filtered_bed_df, annot_df)

            # Store all results
            write_df_to_file(coding_df, file_path.with_stem(file_path.stem + "_coding.csv"))
            write_df_to_file(non_coding_df, file_path.with_stem(file_path.stem + "_non_coding.csv"))
            write_df_to_file(new_annot_df, file_path.with_stem(file_path.stem + "_all_annot_with_methylations.csv"))

    print("Methylome mining done.")
