#!/usr/bin/env python

import os
import sys
import click
import glob
import subprocess


@click.command()
@click.argument("fastq_dir", type=click.Path(exists=True))
@click.argument("chunks", type=int)
def main(fastq_dir, chunks):
    """
    main function for script
    """

    r1_output_files = " ".join(
        [f"-o test_R1_{i:03d}.fastq.gz" for i in range(1, chunks + 1)]
    )
    r2_output_files = " ".join(
        [f"-o test_R2_{i:03d}.fastq.gz" for i in range(1, chunks + 1)]
    )

    R1_path = glob.glob(f"{fastq_dir}/*R1*.fastq.gz")[0]
    R2_path = glob.glob(f"{fastq_dir}/*R2*.fastq.gz")[0]

    subprocess.call(f"fastqsplitter -i {R1_path} {r1_output_files}", shell=True)
    subprocess.call(f"fastqsplitter -i {R2_path} {r2_output_files}", shell=True)


# pylint: disable=no-value-for-parameter
if __name__ == "__main__":
    main()
