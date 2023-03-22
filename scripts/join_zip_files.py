import shutil
import zipfile
import re
import glob
import gzip
import subprocess
import os
import click
import pandas as pd
from pathlib import Path
from collections import defaultdict


def group_fastq_files(fastq_files):
    grouped_files = defaultdict(lambda: [None, None])

    for file in fastq_files:
        basename = os.path.basename(file)
        base = basename.split(".")[0]
        prefix = "_".join(base.split("_")[0:2])
        if base.endswith("mate1"):
            grouped_files[prefix][0] = file
        elif base.endswith("mate2"):
            grouped_files[prefix][1] = file

    return dict(grouped_files)


def combine_gzipped_fastq(input_files, output_file):
    with gzip.open(output_file, "wb") as output_gz:
        for input_file in input_files:
            with gzip.open(input_file, "rb") as input_gz:
                for line in input_gz:
                    output_gz.write(line)


@click.command()
@click.argument("barcode", type=str)
@click.argument("zip_files", nargs=-1, type=click.Path(exists=True))
def main(barcode, zip_files):
    """
    main function for script
    """
    files = []
    for zip_file in zip_files:
        with zipfile.ZipFile(zip_file, "r") as zip_ref:
            for file_info in zip_ref.infolist():
                if file_info.filename not in files:
                    files.append(file_info.filename)
    pairs = group_fastq_files(files)
    outdir = f"output"
    os.makedirs(outdir, exist_ok=True)
    for name, pair in pairs.items():
        count = 0
        for zip_file in zip_files:
            with zipfile.ZipFile(zip_file, "r") as zip_ref:
                zip_ref.extract(pair[0], f"{count}")
                zip_ref.extract(pair[1], f"{count}")
            count += 1
        os.makedirs(f"{outdir}/{barcode}-{name}", exist_ok=True)
        mate1_files = glob.glob(f"*/{pair[0]}")
        mate2_files = glob.glob(f"*/{pair[1]}")
        combine_gzipped_fastq(mate1_files, f"{outdir}/{barcode}-{name}/{pair[0]}")
        combine_gzipped_fastq(mate2_files, f"{outdir}/{barcode}-{name}/{pair[1]}")
        subprocess.call(f"rm -r */{pair[0]}", shell=True)
        subprocess.call(f"rm -r */{pair[1]}", shell=True)
    count = 0
    for zip_file in zip_files:
        shutil.rmtree(f"{count}")
        count += 1


# pylint: disable=no-value-for-parameter
if __name__ == "__main__":
    main()
