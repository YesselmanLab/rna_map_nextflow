import click
import os
import glob
import pandas as pd
import subprocess

from rna_map_nextflow.logger import get_logger, setup_applevel_logger
from rna_map_nextflow.util import random_string
from rna_map_nextflow.fastq import PairedFastqFiles, FastqFile
from rna_map_nextflow.demultiplex import run_demultiplex
from rna_map_nextflow.int_demultiplex import run_int_demultiplex, join_int_demult_files
from rna_map_nextflow.run_rna_map import run_rna_map as rna_map
from rna_map_nextflow.run_rna_map import combine_outputs

log = get_logger("CLI")


@click.group()
def cli():
    pass


@cli.command()
@click.argument("csv")
@click.argument("r1_path", type=click.Path(exists=True))
@click.argument("r2_path", type=click.Path(exists=True))
@click.option("--dtype", type=str, default="sabre")
@click.option("--output-dir", type=click.Path(exists=True), default=None)
@click.option("--debug", is_flag=True)
def demultiplex(csv, r1_path, r2_path, dtype, output_dir, debug):
    """
    demultiplexes paired fastq files given 3' end barcodes
    """
    setup_applevel_logger(file_name="demultiplex.log", is_debug=debug)
    if output_dir is None:
        output_dir = os.getcwd()
    paired_fastqs = PairedFastqFiles(FastqFile(r1_path), FastqFile(r2_path))
    df = pd.read_csv(csv)
    run_demultiplex(df, paired_fastqs, dtype, output_dir)


@cli.command()
@click.argument("barcode", type=str)
@click.argument("input_dir", type=click.Path(exists=True))
@click.argument("result_dirs", nargs=-1, type=click.Path(exists=True))
def combine_rna_map_outputs(barcode, input_dir, result_dirs):
    setup_applevel_logger()
    combine_outputs(barcode, input_dir, result_dirs)


@cli.command()
@click.argument("input_dir", type=click.Path(exists=True))
@click.argument("fastq_dir", type=click.Path(exists=True))
@click.option("-o", "--output-dir", type=click.Path(exists=True), default=None)
def int_demultiplex(input_dir, fastq_dir, output_dir):
    """
    uses barcode_demultiplex to demultiplex based on internal helices
    """
    setup_applevel_logger()
    run_int_demultiplex(input_dir, fastq_dir, output_dir)


@cli.command()
@click.argument("input_dir", type=click.Path(exists=True))
@click.argument("fastq_dir", type=click.Path(exists=True))
@click.option("-o", "--output-dir", type=click.Path(exists=True), default=None)
def run_rna_map(input_dir, fastq_dir, output_dir):
    setup_applevel_logger()
    if output_dir is None:
        output_dir = os.getcwd()
    rna_map(input_dir, fastq_dir, output_dir)


@cli.command()
@click.argument("barcode", type=str)
@click.argument("zip_files", nargs=-1, type=click.Path(exists=True))
@click.option("--threads", type=int, default=1)
@click.option("-o", "--output-dir", type=click.Path(exists=True), default=None)
def join_int_demultiplex_zips(barcode, zip_files, threads, output_dir):
    """
    joins int_demultiplex result zip files together
    """
    setup_applevel_logger()
    if output_dir is None:
        output_dir = os.getcwd()
    join_int_demult_files(barcode, zip_files, threads, output_dir)


@cli.command()
@click.argument("r1_path", type=click.Path(exists=True))
@click.argument("r2_path", type=click.Path(exists=True))
@click.argument("num_chunks", type=int)
@click.option("--output-dir", type=click.Path(exists=True), default=None)
@click.option("--debug", is_flag=True)
def split_fastqs(r1_path, r2_path, num_chunks, output_dir, debug):
    setup_applevel_logger()
    if output_dir is None:
        output_dir = os.getcwd()
    r1_output_files = ""
    r2_output_files = ""
    for i in range(1, num_chunks + 1):
        uid = random_string(10)
        r1_output_files += f"-o test_R1_{uid}_.fastq.gz "
        r2_output_files += f"-o test_R2_{uid}_.fastq.gz "

    if debug:
        log.info(f"fastqsplitter -i {r1_path} {r1_output_files}")
        log.info(f"fastqsplitter -i {r2_path} {r2_output_files}")
    subprocess.call(f"fastqsplitter -i {r1_path} {r1_output_files}", shell=True)
    subprocess.call(f"fastqsplitter -i {r2_path} {r2_output_files}", shell=True)


if __name__ == "__main__":
    cli()
