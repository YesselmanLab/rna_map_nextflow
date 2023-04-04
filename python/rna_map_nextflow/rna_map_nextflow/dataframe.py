"""
functions for handling dataframes of data
"""

import pandas as pd
from seq_tools.dataframe import trim as seq_ss_trim


def validate_data(df):
    """
    validates the data in the dataframe
    :param df: a dataframe with data
    :return: None
    """
    expects = ["barcode", "barcode_seq", "construct"]
    for e in expects:
        if e not in df:
            raise ValueError(f"{e} column is required in the dataframe")


def trim(df: pd.DataFrame, start: int, end: int) -> pd.DataFrame:
    """
    trims the dataframe to the given start and end
    :param df: a dataframe with data
    :param start: the start index
    :param end: the end index
    :return: a trimmed dataframe
    """
    df = seq_ss_trim(df, start, end)
    df["data"] = df["data"].apply(lambda x: x[start:end])
    return df


def check_if_columns_exist(df: pd.DataFrame, cols) -> None:
    """
    validates the data in the dataframe
    :param df: a dataframe with data
    :return: None
    """
    for e in cols:
        if e not in df:
            raise ValueError(
                f"{e} column is required in the dataframe! colums are "
                f"{list(df.columns)} "
            )
