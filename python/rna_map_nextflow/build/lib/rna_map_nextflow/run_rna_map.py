import pandas as pd
from pathlib import Path
import glob
import shutil

import rna_map
from rna_map.parameters import get_preset_params
from rna_map.logger import get_logger, setup_applevel_logger

from rna_map_nextflow.util import random_string


def write_fasta_file(row):
    f = open("test.fasta", "w")
    f.write(f">{row['name']}\n")
    f.write(f"{row['sequence'].replace('U', 'T')}\n")
    f.close()


def write_dotbracket_file(row):
    f = open("test.csv", "w")
    f.write(f"name,sequence,structure\n")
    f.write(f"{row['name']},{row['sequence'].replace('U', 'T')},{row['structure']}\n")
    f.close()


def run_rna_map(input_dir, fastq_dir):
    """
    main function for script
    """
    df = pd.read_csv(f"{input_dir}/data.csv")
    # parse barcode information in directory name
    # in the format of construct_barcode-seq_barcode
    construct_barcode, seq_barcode = fastq_dir.split("-")
    construct_row = df[df["barcode_seq"] == construct_barcode].iloc[0]
    construct_code = construct_row["code"]
    df_barcodes = pd.read_json(f"{input_dir}/barcode_jsons/{construct_code}.json")
    seq_row = df_barcodes[df_barcodes["full_barcode"] == seq_barcode].iloc[0]
    write_fasta_file(seq_row)
    write_dotbracket_file(seq_row)
    mate_1_path = glob.glob(f"{fastq_dir}/*_mate1.fastq.gz")[0]
    mate_2_path = glob.glob(f"{fastq_dir}/*_mate2.fastq.gz")[0]
    # run rna_map pipeline
    setup_applevel_logger()
    rna_map.run.run(
        "test.fasta",
        mate_1_path,
        mate_2_path,
        "test.csv",
        get_preset_params("barcoded-library"),
    )
    # clean up unncessary files
    shutil.rmtree("log")
    shutil.rmtree("input")
    shutil.rmtree("output/Mapping_Files")
    # move to unique directory name to avoid collisions in nextflow
    rand_str = random_string(10)
    shutil.move("output/BitVector_Files", f"output-{construct_barcode}-{rand_str}")
    shutil.rmtree("output")
