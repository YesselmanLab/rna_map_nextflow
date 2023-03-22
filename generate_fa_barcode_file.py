import pandas as pd
import click

def generate_barcode_file(df):
    """
    generates a .fa file for novobarcode
    :param df: a dataframe with that contains informmation of barcode sequences
    :return: None
    """
    expects = ["barcode", "barcode_seq", "construct"]
    for e in expects:
        if e not in df:
            raise ValueError(f"{e} column is required in the csv")
    s = "Distance\t4\nFormat\t5\n"
    seen = []
    for i, row in df.iterrows():
        if row["barcode"] in seen:
            continue
        s += f"{row['barcode']}\t{row['barcode_seq']}\n"
        seen.append(row["barcode"])
    f = open("rtb_barcodes.fa", "w")
    f.write(s)
    f.close()

@click.command()
@click.argument("data", type=click.Path(exists=True))
def main(data):
    """
    main function for script
    """
    df = pd.read_csv(data)
    generate_barcode_file(df)

# pylint: disable=no-value-for-parameter
if __name__ == '__main__':
    main()
