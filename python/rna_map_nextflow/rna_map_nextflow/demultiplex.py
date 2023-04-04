"""
handles demultiplexing of fastq files currently using novobaracode / sabre
"""
import os
import shutil
import subprocess
import pandas as pd
from tabulate import tabulate
from pathlib import Path
import gzip

from rna_map_nextflow.dataframe import check_if_columns_exist
from rna_map_nextflow.logger import get_logger
from rna_map_nextflow.fastq import PairedFastqFiles

log = get_logger("DEMULTIPLEX")


def run_demultiplex(df, paired_fqs, demultiplex_type, demultiplex_path, params=None):
    """
    demultiplexes fastq files
    :param df: dataframe containing sample information
    :param paired_fqs: PairedFastqFiles object that contains the paths to
    paired read fastq files
    :param demultiplex_path: path to the directory where demultiplexed fastq
    files will be stored
    :param params: demultiplexing parameters
    :return: None
    """
    if params is None:
        params = {}
    if demultiplex_type == "novobarcode":
        demultiplexer = NovobarcodeDemultiplexer()
    elif demultiplex_type == "sabre":
        demultiplexer = SabreDemultiplexer()
    else:
        raise ValueError(f"Demultiplex type {demultiplex_type} not supported")
    demultiplexer.setup(params)
    demultiplexer.run(df, paired_fqs, demultiplex_path)


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


class Demultiplexer:
    """
    An abstract class for demultiplexing fastq files
    """

    def setup(self, params) -> None:
        """
        setup demultiplexer
        :param params: setup parameters which are in the format of a dictionary
        :return: None
        """
        self._params = params

    def run(
        self,
        df: pd.DataFrame,
        paired_fqs: PairedFastqFiles,
        demultiplex_path,
    ):
        pass

    def _uncompress_fastqs(self, paired_fqs: PairedFastqFiles) -> None:
        """
        handles unzipping of fastq files if the imput fastq files are compressed
        :param paired_fqs: PairedFastqFiles object that contains the paths to
        paired read fastq files
        :return: None
        """
        log.info("unzipping fastq files")
        shutil.copy(paired_fqs.read_1.path, "test_S1_L001_R1_001.fastq.gz")
        shutil.copy(paired_fqs.read_2.path, "test_S1_L001_R2_001.fastq.gz")
        subprocess.call(f"gunzip --force *.gz", shell=True)
        log.info("read1 fastq -> test_S1_L001_R1_001.fastq")
        log.info("read2 fastq -> test_S1_L001_R2_001.fastq")

    def _prepare_fastq_files(self, paired_fqs: PairedFastqFiles):
        """
        copies and decompresses fastq files to the current working directory
        :param paired_fqs: PairedFastqFiles object that contains the paths to
        paired read fastq files
        :return:
        """
        if paired_fqs.is_compressed():
            self._uncompress_fastqs(paired_fqs)
        else:
            shutil.copy2(paired_fqs.read_1.path, "test_S1_L001_R1_001.fastq")
            shutil.copy2(paired_fqs.read_2.path, "test_S1_L001_R2_001.fastq")
            log.info(f"copying {paired_fqs.read_1.path} -> test_S1_L001_R1_001.fastq")
            log.info(f"copying {paired_fqs.read_2.path} -> test_S1_L001_R2_001.fastq")


class NovobarcodeDemultiplexer(Demultiplexer):
    def run(
        self,
        df: pd.DataFrame,
        paired_fqs: PairedFastqFiles,
        demultiplex_path,
    ) -> None:
        """ """
        if not os.path.isdir(demultiplex_path):
            log.error(f"{demultiplex_path} does not exist")
            exit()
        os.chdir(demultiplex_path)
        self._prepare_fastq_files(paired_fqs)
        log.info("preparing rtb_barcodes.fa file for demultiplexing")
        self.__generate_barcode_file(df)
        output = subprocess.check_output(
            "novobarcode -b rtb_barcodes.fa -f test_S1_L001_R1_001.fastq "
            "test_S1_L001_R2_001.fastq",
            shell=True,
        )
        output = output.decode("UTF-8")
        log.info(f"output from novobarcode:\n{output}")
        df_demult = self.__parse_novobarcode_stdout(output)
        log.info(f"total number of reads: {df_demult['count'].sum()}")
        log.info(
            f"total number of data reads: "
            f"{df_demult['count'].sum() - df_demult.iloc[-1]['count']}"
        )
        if not os.path.isdir("NC"):
            log.error("nothing was generated STOPING NOW!")
            exit()
        if self._params["delete_fastqs"]:
            log.info("deleting copied fastq files")
            os.remove("test_S1_L001_R1_001.fastq")
            os.remove("test_S1_L001_R2_001.fastq")
        if self._params["delete_non_barcoded"]:
            log.info("deleting reads that do not have a barcode")
            shutil.rmtree("NC")

    def __parse_novobarcode_stdout(self, output) -> pd.DataFrame:
        """
        parses the demultiplex log file
        :param output: the output of the demultiplex command
        :return:
        """
        lines = output.split("\n")
        data = []
        for l in lines:
            if l.startswith("#"):
                continue
            spl = l.split()
            # check if last row has number of reads
            try:
                float(spl[-1])
            except:
                continue
            data.append(spl)
        df = pd.DataFrame(data, columns="id,tag,count".split(","))
        df["count"] = pd.to_numeric(df["count"])
        df.to_csv("demultiplex.csv", index=False)
        return df

    def __generate_barcode_file(self, df, fname="rtb_barcodes.fa"):
        """
        generates a .fa file for novobarcode
        :param df: a dataframe with that contains informmation of barcode
        sequences
        :param fname: the name of the file to write to
        :return: None
        """
        expects = ["barcode", "barcode_seq", "construct"]
        check_if_columns_exist(df, expects)
        s = "Distance\t4\nFormat\t5\n"
        seen = []
        warning = False
        log.info(
            "constructs:\n\n"
            + tabulate(
                df[expects],
                expects,
                tablefmt="github",
                showindex=False,
            )
            + "\n"
        )
        for i, row in df.iterrows():
            if row["barcode"] in seen:
                log.warning(
                    f"{row['barcode']} has been used more than once this may be an "
                    f"issue"
                )
                warning = True
                continue
            s += f"{row['barcode']}\t{row['barcode_seq']}\n"
            seen.append(row["barcode"])
        log.info(f"{len(seen)} unique barcodes found from csv file")
        if not warning:
            log.info("no barcode conflicts detected")
        with open(fname, "w", encoding="utf8") as f:
            f.write(s)


class SabreDemultiplexer(Demultiplexer):
    def run(
        self,
        df: pd.DataFrame,
        paired_fqs: PairedFastqFiles,
        demultiplex_path,
    ):
        if not os.path.isdir(demultiplex_path):
            log.error(f"{demultiplex_path} does not exist")
            exit()
        os.chdir(demultiplex_path)
        log.info("preparing barcodes.txt file for demultiplexing")
        self.__generate_barcode_file(df)
        r1_path = paired_fqs.read_1.path
        r2_path = paired_fqs.read_2.path
        output = subprocess.check_output(
            f" sabre pe -f {r1_path} -r {r2_path} -b barcode.txt "
            f"-u NC/test_S1_L001_R1_001.fastq -w NC/test_S1_L001_R2_001.fastq -m 4",
            shell=True,
        )
        output = output.decode("UTF-8")
        log.info(f"output from sabre:\n{output}")
        for _, row in df.iterrows():
            gzip_files(row["barcode_seq"])

    def __generate_barcode_file(self, df, fname="barcode.txt"):
        expects = ["barcode", "barcode_seq", "construct"]
        check_if_columns_exist(df, expects)
        seen = []
        warning = False
        log.info(
            "constructs:\n\n"
            + tabulate(
                df[expects],
                expects,
                tablefmt="github",
                showindex=False,
            )
            + "\n"
        )
        s = ""
        for i, row in df.iterrows():
            if row["barcode"] in seen:
                log.warning(
                    f"{row['barcode']} has been used more than once this may "
                    f"be an issue"
                )
                warning = True
                continue
            s += f"{row['barcode_seq']}\t"
            s += f"{row['barcode_seq']}/test_S1_L001_R1_001.fastq\t"
            s += f"{row['barcode_seq']}/test_S1_L001_R2_001.fastq\n"
            os.makedirs(f"{row['barcode_seq']}", exist_ok=True)
            seen.append(row["barcode"])
        os.makedirs("NC", exist_ok=True)
        log.info(f"{len(seen)} unique barcodes found from csv file")
        if not warning:
            log.info("no barcode conflicts detected")
        with open(fname, "w", encoding="utf8") as f:
            f.write(s)
