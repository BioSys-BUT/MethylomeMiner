"""
MethylomeMiner module contains functions to run the tool to detect methylations
"""
from pathlib import Path

import click

from loading import load_bed_file


@click.command()
@click.option(
    '-i', '--input_file',
    required=True,
    type=click.Path(
        exists=True, file_okay=True, dir_okay=False,
        readable=True, resolve_path=True),
    help='Path to a BED file with methylation data.',
)
@click.option(
    '-o', '--output_file',
    required=True,
    type=click.Path(
        file_okay=True, dir_okay=False, writable=True,
        resolve_path=True, path_type=Path),
    help='Path to a new JSON file with detected methylations.',
)
def mine_methylations(input_file, output_file):
    bed_df, is_valid = load_bed_file(input_file)
    print(is_valid)


def main():
    mine_methylations()


if __name__ == "__main__":
    main()
