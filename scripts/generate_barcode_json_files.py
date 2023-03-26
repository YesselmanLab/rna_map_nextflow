import click
import os
import pandas as pd
from pathlib import Path

from barcode_demultiplex.demultiplex import find_helix_barcodes


@click.command()
@click.argument("csv", type=click.Path(exists=True))
def main(csv):
    """
    main function for script
    """
    df = pd.read_csv(csv)
    seq_path = os.environ["SEQPATH"]
    os.makedirs("barcode_jsons", exist_ok=True)
    for i, row in df.iterrows():
        df_seq = pd.read_csv(f"{seq_path}/rna/{row['code']}.csv")
        helices = []
        args = row["demult_args"].split()
        for i in range(0, len(args)):
            if args[i] == "--helix" or args[i] == "-helix":
                helices.append([int(args[i + 1]), int(args[i + 2]), int(args[i + 3])])
        df_barcodes = find_helix_barcodes(df_seq, helices)
        df_barcodes.to_json(f"inputs/barcode_jsons/{row['code']}.json", orient="records")


# pylint: disable=no-value-for-parameter
if __name__ == "__main__":
    main()
