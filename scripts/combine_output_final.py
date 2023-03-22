import subprocess
import click
import os
import random
import string
import glob
import shutil
import pandas as pd
from pathlib import Path


from rna_map.mutation_histogram import (
    get_mut_histos_from_pickle_file,
    merge_all_merge_mut_histo_dicts,
    write_mut_histos_to_pickle_file,
    write_mut_histos_to_json_file,
)


@click.command()
@click.argument("csv", type=click.Path(exists=True))
@click.argument("barcode", type=str)
@click.argument("result_dirs", nargs=-1, type=click.Path(exists=True))
def main(csv, barcode, result_dirs):
    """
    main function for script
    """
    df = pd.read_csv(csv)
    row = df.loc[df["barcode_seq"] == barcode].iloc[0]
    dir_name = row["construct"] + "_" + row["code"] + "_" + row["data_type"]
    dir_path = f"{dir_name}/output/BitVector_Files/"
    os.makedirs(dir_path, exist_ok=True)
    mhs = []
    for result_dir in result_dirs:
        files = glob.glob(result_dir + "/*")
        for file in files:
            if file.endswith(".p"):
                mhs.append(get_mut_histos_from_pickle_file(file))
            # Not sure if this is necessary but read in the file contents instead
            # of just copying it to ensure its not the symlink
            elif file.endswith(".html") or file.endswith(".txt"):
                f = open(file, "r")
                content = f.read()
                f.close()
                f = open(os.path.join(dir_path, os.path.basename(file)), "w")
                f.write(content)
                f.close()
    merged_mhs = merge_all_merge_mut_histo_dicts(mhs)
    write_mut_histos_to_pickle_file(merged_mhs, dir_path + "mutation_histos.p")
    write_mut_histos_to_json_file(merged_mhs, dir_path + "mutation_histos.json")


# pylint: disable=no-value-for-parameter
if __name__ == "__main__":
    main()
