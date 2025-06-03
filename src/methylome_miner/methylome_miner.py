"""
MethylomeMiner module contains functions to run the tool to detect methylations
"""
from pathlib import Path

import click

from .loading import parse_annotation_file, pair_bed_and_annot_files
from .backend import sort_methylations, write_df_to_file, get_core_methylome
from .preprocessing import get_list_of_methylated_positions


class Mutex(click.Option):
    def __init__(self, *args, **kwargs):
        self.not_required_if:list = kwargs.pop("not_required_if")

        assert self.not_required_if, "'not_required_if' parameter required"
        kwargs["help"] = (kwargs.get("help", "") + "Option is mutually exclusive with " + ", ".join(self.not_required_if) + ".").strip()
        super(Mutex, self).__init__(*args, **kwargs)

    def handle_parse_result(self, ctx, opts, args):
        current_opt:bool = self.name in opts
        for mutex_opt in self.not_required_if:
            if mutex_opt in opts:
                if current_opt:
                    raise click.UsageError("Illegal usage: '" + str(self.name) + "' is mutually exclusive with " + str(mutex_opt) + ".")
                else:
                    self.prompt = None
        return super(Mutex, self).handle_parse_result(ctx, opts, args)


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
    '--input_annot_file',
    required=True,
    type=click.Path(
        exists=True, file_okay=True, dir_okay=False,
        readable=True, resolve_path=True),
    help='Path to a file with genome annotation.',
)
@click.option(
    '--input_bed_dir',
    required=True, cls=Mutex, not_required_if="min_coverage",  # Funguje Mutex? A co když jsou obě varianty?
    type=click.Path(
        exists=True, dir_okay=True,
        readable=True, resolve_path=True),
    help='Path to a dir with BED files with methylation data.',
)
@click.option(
    '--min_coverage',
    required=True, cls=Mutex, not_required_if="input_bed_dir",
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
@click.option(
    '--work_dir',
    required=False,
    default="MethylomeMiner_output",
    type=click.Path(
        dir_okay=True, writable=True,
        resolve_path=True, path_type=Path),
    help='Path to files with detected methylations.',
)
@click.option(
    '--file_name',
    required=False,
    default="None",
    help='New file name for files with detected methylations.',
)
def mine_methylations(input_bed_file, input_annot_file, input_bed_dir, min_coverage, min_percent_modified, work_dir, file_name):

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

        if file_name is None:
            file_name = input_bed_file.stem.split("_")[0].split(".")[0]
        if not work_dir.exists():
            work_dir.mkdir(parents=True, exist_ok=True)
        write_df_to_file(coding_df, Path(work_dir, file_name + "_coding.csv"))
        write_df_to_file(non_coding_df, Path(work_dir, file_name + "_non_coding.csv"))
        write_df_to_file(new_annot_df, Path(work_dir, file_name + "_annot_with_methylations.csv"))
    print("Methylome mining done.")


@click.command()
@click.option(
    '--input_bed_dir',
    required=True,
    type=click.Path(
        exists=True, dir_okay=True,
        readable=True, resolve_path=True),
    help='Path to a dir with BED files with methylation data.',
)
@click.option(
    '--input_annot_dir',
    required=False,
    type=click.Path(
        exists=True, dir_okay=True,
        readable=True, resolve_path=True),
    help='Path to a dir with annotation files for BED files.',
)
@click.option(
    "--roary_output_file",
    required=True,
    type=click.Path(
        exists=True, file_okay=True, dir_okay=False,
        readable=True, resolve_path=True),
    help='Path to a output file ("gene_presence_absence.csv") from Roary tool. Required only for core methylome analysis.',
)
@click.option(
    '--work_dir',
    required=False,
    default="MethylomeMiner_output",
    type=click.Path(
        dir_okay=True, writable=True,
        resolve_path=True, path_type=Path),
    help='Path to files with detected core methylations.',
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
@click.option(
    '--matrix_values',
    required=False,
    default="presence",
    help='Select type of values in the core methylome matrix. Default "presence".'
         '"presence" - 0 for no methylation, 1 for any methylation; '
         '"positions" - values of exact location of present methylations in a gene.',
)
def mine_core_methylome(input_bed_dir, input_annot_dir, roary_output_file, work_dir, min_coverage, min_percent_modified, matrix_values):
    files_pairs = pair_bed_and_annot_files(input_bed_dir, input_annot_dir)
    for prefix, files in files_pairs.items():
        if not Path(work_dir, prefix + "_annot_with_methylations.csv").exists():
            bed_file = files["bed_file"]
            annot_file = files["annot_file"]
            bed_df = get_list_of_methylated_positions(
                bed_file,
                input_bed_dir,
                write_to_file=False,
                min_coverage=min_coverage,
                min_percent_modified=min_percent_modified
            )

            if bed_df is not None:
                annot_df = parse_annotation_file(Path(annot_file))
                coding_df, non_coding_df, new_annot_df = sort_methylations(bed_df, annot_df)
                write_df_to_file(new_annot_df, Path(work_dir, prefix + "_annot_with_methylations.csv"))
    core_methylomes = get_core_methylome(roary_output_file, miner_output_dir=work_dir, matrix_values=matrix_values)

    for methylation, methylome in core_methylomes.items():
        if not methylome.iloc[:, 1:].isnull().values.all():
            write_df_to_file(methylome, Path(work_dir, methylation + "_core_methylome.csv"))
    print("Core methylome mining done.")


if __name__ == "__main__":
    mine_methylations()
    # mine_core_methylome()
