#!/usr/bin/env python

import shutil
import random
import string
import re
import os
import subprocess
import zipfile
import click
import pandas as pd
from pathlib import Path


def flatten_directory(input_directory, output_zip):
    with zipfile.ZipFile(output_zip, "w") as zip_ref:
        for root, _, files in os.walk(input_directory):
            for file in files:
                file_path = os.path.join(root, file)
                # Save the file in the zip with only its base name
                zip_ref.write(file_path, os.path.basename(file))


def random_string(length):
    return "".join(random.choices(string.ascii_letters, k=length))


def find_paired_fastq_files(directory):
    all_fastq_files = list(Path(directory).glob("**/*.fastq.gz"))
    paired_files = {}

    for fastq_file in all_fastq_files:
        basename = fastq_file.name
        match = re.match(r"(.+)_S\d+_L\d+_R([12])_\d+\.fastq\.gz", basename)

        if match:
            sample_name = match.group(1)
            read_pair = int(match.group(2))

            if sample_name in paired_files:
                paired_files[sample_name][read_pair - 1] = fastq_file
            else:
                paired_files[sample_name] = [None, None]
                paired_files[sample_name][read_pair - 1] = fastq_file

    # Remove incomplete pairs
    paired_files = {key: value for key, value in paired_files.items() if all(value)}

    return paired_files


@click.command()
@click.argument("csv", type=click.Path(exists=True))
@click.argument("dir_path", type=click.Path(exists=True))
def main(csv, dir_path):
    """
    main function for script
    """
    # get fastq files
    paired_fastqs = find_paired_fastq_files(dir_path)
    if len(paired_fastqs) != 1:
        raise ValueError("No fastq files found in directory")
    paired_fastqs = list(paired_fastqs.values())[0]
    R1_path = paired_fastqs[0]
    R2_path = paired_fastqs[1]
    # get the construct information based on barcode_seq in folder name
    df = pd.read_csv(csv)
    barcode_seq = Path(dir_path).stem
    df_sub = df.loc[df["barcode_seq"] == barcode_seq]
    if df_sub.empty:
        raise ValueError("No barcode found in csv file")
    row = df_sub.iloc[0]
    unique_code = random_string(10)
    # data_path = f"/tmp/bc-path/{unique_code}"
    data_path = "data"
    # cmd = (
    #    f"barcode_demultiplex -csv /inputs/rnas/{row['code']}.csv "
    #    f"-fq1 {R2_path} -fq2 {R1_path} --data-path {data_path} {row['demult_args']}"
    # )
    cmd = (
        f"barcode_demultiplex -csv ../../../inputs/rnas/{row['code']}.csv "
        f"-fq1 {R2_path} -fq2 {R1_path} --data-path {data_path} {row['demult_args']}"
    )
    subprocess.call(cmd, shell=True)
    # subprocess.call(
    #    f"zip -r {barcode_seq}_.{unique_code}.demultiplexed.zip {data_path}", shell=True
    # )
    flatten_directory(data_path, f"{barcode_seq}_.{unique_code}.demultiplexed.zip")
    # shutil.rmtree(data_path)


# pylint: disable=no-value-for-parameter
if __name__ == "__main__":
    main()
