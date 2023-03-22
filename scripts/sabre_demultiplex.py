import click
import pandas as pd
import subprocess
import os
import shutil
import gzip


def generate_barcode_file(df, fname="barcode.txt"):
    seen = []
    warning = False
    s = ""
    for i, row in df.iterrows():
        if row["barcode"] in seen:
            continue
        s += f"{row['barcode_seq']}\t"
        s += f"{row['barcode_seq']}/test_S1_L001_R1_001.fastq\t"
        s += f"{row['barcode_seq']}/test_S1_L001_R2_001.fastq\n"
        # os.makedirs(f"{row['barcode_seq']}", exist_ok=True)
        seen.append(row["barcode"])
    with open(fname, "w", encoding="utf8") as f:
        f.write(s)


def gzip_files(directory):
    for root, dirs, files in os.walk(directory):
        for file in files:
            if not file.endswith(".gz"):  # Ignore already compressed files
                file_path = os.path.join(root, file)
                compressed_file_path = f"{file_path}.gz"
                with open(file_path, "rb") as f_in:
                    with gzip.open(compressed_file_path, "wb") as f_out:
                        shutil.copyfileobj(f_in, f_out)

                os.remove(file_path)  # Remove the original file


@click.command()
@click.argument("csv", type=click.Path(exists=True))
@click.argument("r1_path", type=click.Path(exists=True))
@click.argument("r2_path", type=click.Path(exists=True))
def main(csv, r1_path, r2_path):
    """
    main function for script
    """
    df = pd.read_csv(csv)
    generate_barcode_file(df)
    os.makedirs("NC", exist_ok=True)
    for _, row in df.iterrows():
        os.makedirs(row["barcode_seq"], exist_ok=True)
    cmd = (
        f"sabre pe -f {r1_path} -r {r2_path} -b barcode.txt -m 4 "
        f"-u NC/test_S1_L001_R1_001.fastq -w NC/test_S1_L001_R2_001.fastq"
    )
    subprocess.call(cmd, shell=True)
    shutil.rmtree("NC")
    for _, row in df.iterrows():
        gzip_files(row["barcode_seq"])


# pylint: disable=no-value-for-parameter
if __name__ == "__main__":
    main()
