import glob
import os
import gzip
import zipfile
import pandas as pd
import shutil
import subprocess
from pathlib import Path
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor

from barcode_demultiplex.demultiplex import demultiplex

from rna_map_nextflow.util import random_string, flatten_directory


def run_int_demultiplex(input_dir, fastq_dir, output_dir=None):
    if output_dir is None:
        output_dir = os.getcwd()
    df = pd.read_csv(f"{input_dir}/data.csv")
    fastq1_path = glob.glob(os.path.join(fastq_dir, "*R1*.fastq.gz"))[0]
    fastq2_path = glob.glob(os.path.join(fastq_dir, "*R2*.fastq.gz"))[0]
    barcode_seq = Path(fastq_dir).stem
    df_sub = df.loc[df["barcode_seq"] == barcode_seq]
    if df_sub.empty:
        raise ValueError("No barcode found in csv file")
    row = df_sub.iloc[0]
    # get helices from commandline
    helices = []
    args = row["demult_args"].split()
    for i in range(0, len(args)):
        if args[i] == "--helix" or args[i] == "-helix":
            helices.append([int(args[i + 1]), int(args[i + 2]), int(args[i + 3])])
    unique_code = random_string(10)
    data_path = f"{output_dir}/{unique_code}"
    df_seq = pd.read_csv(f"{input_dir}/rnas/{row['code']}.csv")
    demultiplex(df_seq, Path(fastq2_path), Path(fastq1_path), helices, data_path)
    flatten_directory(data_path, f"{barcode_seq}_.{unique_code}.demultiplexed.zip")
    shutil.rmtree(data_path)


# join zip results ##############################################################


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


def process_pair(name_pair, barcode, zip_files, outdir, tmp_dir):
    name, pair = name_pair
    count = 0
    for zip_file in zip_files:
        with zipfile.ZipFile(zip_file, "r") as zip_ref:
            zip_ref.extract(pair[0], f"{tmp_dir}/{count}")
            zip_ref.extract(pair[1], f"{tmp_dir}/{count}")
        count += 1
    os.makedirs(f"{outdir}/{barcode}-{name}", exist_ok=True)
    mate1_files = glob.glob(f"{tmp_dir}/*/{pair[0]}")
    mate2_files = glob.glob(f"{tmp_dir}/*/{pair[1]}")
    combine_gzipped_fastq(mate1_files, f"{outdir}/{barcode}-{name}/{pair[0]}")
    combine_gzipped_fastq(mate2_files, f"{outdir}/{barcode}-{name}/{pair[1]}")
    subprocess.call(f"rm -r {tmp_dir}/*/{pair[0]}", shell=True)
    subprocess.call(f"rm -r {tmp_dir}/*/{pair[1]}", shell=True)


def join_int_demult_files(barcode, zip_files, threads=1, tmp_dir=None):
    if tmp_dir is None:
        tmp_dir = os.getcwd()
    files = []
    for zip_file in zip_files:
        with zipfile.ZipFile(zip_file, "r") as zip_ref:
            for file_info in zip_ref.infolist():
                if file_info.filename not in files:
                    files.append(file_info.filename)
    pairs = group_fastq_files(files)
    outdir = f"output"
    os.makedirs(outdir, exist_ok=True)
    tmp_dir = tmp_dir + "/" + random_string(10)
    os.makedirs(tmp_dir, exist_ok=True)

    with ProcessPoolExecutor(max_workers=threads) as executor:
        executor.map(
            process_pair,
            pairs.items(),
            [barcode] * len(pairs),
            [zip_files] * len(pairs),
            [outdir] * len(pairs),
            [tmp_dir] * len(pairs),
        )

    count = 0
    for zip_file in zip_files:
        shutil.rmtree(f"{tmp_dir}/{count}")
        count += 1
