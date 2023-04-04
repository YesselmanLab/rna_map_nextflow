import os
import glob
from dataclasses import dataclass


@dataclass(frozen=True, order=True)
class FastqFile:
    """
    Holds the path to a fastq file
    """

    path: str

    def __post_init__(self):
        """
        Check that the file exists
        """
        if not os.path.exists(self.path):
            raise FileNotFoundError(f"File {self.path} does not exist")
        if not self.is_r1() and not self.is_r2():
            raise ValueError(
                f"File {self.path} is not R1 or R2 please check the file name "
                f"must have either _R1_ or _R2_ in the name"
            )

    def is_compressed(self):
        """
        Check if the file is compressed
        """
        if self.path.endswith(".gz"):
            return True
        return False

    def is_r1(self):
        """
        Check if the file is R1
        """
        if "_R1_" in self.path:
            return True
        return False

    def is_r2(self):
        """
        Check if the file is R2
        """
        if "_R2_" in self.path:
            return True
        return False


@dataclass(frozen=True, order=True)
class PairedFastqFiles:
    """
    Holds the paths to paired fastq files
    """

    read_1: FastqFile
    read_2: FastqFile

    def is_compressed(self):
        """
        Check if the files are compressed
        """
        if self.read_1.is_compressed() and self.read_2.is_compressed():
            return True
        return False


def get_paired_fastqs(dir_path: str) -> PairedFastqFiles:
    """
    Get the paired fastq files from a directory
    :dir_path: path to directory
    :return: list of paired fastq files
    """
    if os.path.isdir(dir_path):
        f1_paths = glob.glob(os.path.join(dir_path, "*_R1_*"))
        f2_paths = glob.glob(os.path.join(dir_path, "*_R2_*"))
    else:
        f1_paths = glob.glob(dir_path + "*_R1_*")
        f2_paths = glob.glob(dir_path + "*_R2_*")
    if len(f1_paths) == 0 or len(f2_paths) == 0:
        raise FileNotFoundError(
            f"Could not find paired fastq files in {dir_path}"
        )
    if len(f1_paths) > 1 or len(f2_paths) > 1:
        raise ValueError(f"Found more than one fastq file in {dir_path}")
    return PairedFastqFiles(FastqFile(f1_paths[0]), FastqFile(f2_paths[0]))
