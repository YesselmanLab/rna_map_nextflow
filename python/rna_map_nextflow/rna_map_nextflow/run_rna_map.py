import glob
import os
import shutil
import glob
import pandas as pd
from pathlib import Path

import rna_map
from rna_map.parameters import get_preset_params
from rna_map.logger import setup_applevel_logger
from rna_map.mutation_histogram import merge_mut_histo_files

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


def get_file_size(file_path):
    file_path = os.path.realpath(file_path)
    return os.path.getsize(file_path)


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
    fsize_1 = get_file_size(mate_1_path)
    fsize_2 = get_file_size(mate_2_path)
    rand_str = random_string(10)
    if fsize_1 < 100 or fsize_2 < 100:
        print(f"skipping {fastq_dir} because file size is too small")
        os.mkdir(f"output-{construct_barcode}-{rand_str}")
        return
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
    shutil.move("output/BitVector_Files", f"output-{construct_barcode}-{rand_str}")
    shutil.rmtree("output")


def keep_file(file_path):
    _, extension = os.path.splitext(file_path)
    return extension in {".html", ".txt"}


def combine_outputs(barcode, input_dir, result_dirs):
    """
    combines all of the mutation histogram files into one file
    """
    df = pd.read_csv(f"{input_dir}/data.csv")
    row = df.loc[df["barcode_seq"] == barcode].iloc[0]
    dir_name = row["construct"] + "_" + row["code"] + "_" + row["data_type"]
    dir_path = f"{dir_name}/output/BitVector_Files/"
    os.makedirs(dir_path, exist_ok=True)
    mhs_files = []
    for result_dir in result_dirs:
        files = glob.glob(result_dir + "/*")
        for file in files:
            if keep_file(file):
                # since files are symlinked, we need to copy the content
                f = open(file, "r")
                content = f.read()
                f.close()
                f = open(os.path.join(dir_path, os.path.basename(file)), "w")
                f.write(content)
                f.close()
            elif file.endswith(".p"):
                mhs_files.append(file)
    merge_mut_histo_files(mhs_files, dir_path)
