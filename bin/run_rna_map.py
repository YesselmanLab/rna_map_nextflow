#!/usr/bin/env python

import subprocess
import click
import os
import random
import string
import glob
import shutil
import pandas as pd
from pathlib import Path


def random_string(length):
    return "".join(random.choices(string.ascii_letters, k=length))


@click.command()
@click.argument("base_dir", type=click.Path(exists=True))
@click.argument("fastq_dir", type=click.Path(exists=True))
def main(base_dir, fastq_dir):
    """
    main function for script
    """
    fastq_dir = Path(fastq_dir)
    full_barcode = fastq_dir.stem
    construct_barcode, seq_barcode = full_barcode.split("-")
    df = pd.read_csv(f"/inputs/data.csv")
    construct_row = df[df["barcode_seq"] == construct_barcode].iloc[0]
    construct_code = construct_row["code"]
    df_barcodes = pd.read_json(f"/inputs/barcode_jsons/{construct_code}.json")
    seq_row = df_barcodes[df_barcodes["full_barcode"] == seq_barcode].iloc[0]
    f = open("test.fasta", "w")
    f.write(f">{seq_row['name']}\n")
    f.write(f"{seq_row['sequence'].replace('U', 'T')}\n")
    f.close()
    f = open("test.csv", "w")
    f.write(f"name,sequence,structure\n")
    f.write(
        f"{seq_row['name']},{seq_row['sequence'].replace('U', 'T')},{seq_row['structure']}\n"
    )
    f.close()
    mate_1_path = glob.glob(f"{fastq_dir}/*_mate1.fastq.gz")[0]
    mate_2_path = glob.glob(f"{fastq_dir}/*_mate2.fastq.gz")[0]
    cmd = (
        f"rna-map -fa test.fasta -fq1 {mate_1_path} -fq2 {mate_2_path} "
        f"--dot-bracket test.csv --param-preset 'barcoded-library' "
    )
    subprocess.run(cmd, shell=True)
    shutil.rmtree("log")
    shutil.rmtree("input")
    shutil.rmtree("output/Mapping_Files")
    rand_str = random_string(10)
    shutil.move("output/BitVector_Files", f"output-{construct_barcode}-{rand_str}")
    shutil.rmtree("output")


# pylint: disable=no-value-for-parameter
if __name__ == "__main__":
    main()
