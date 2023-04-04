import os
import random
import string
import zipfile


def flatten_directory(input_directory, output_zip):
    with zipfile.ZipFile(output_zip, "w") as zip_ref:
        for root, _, files in os.walk(input_directory):
            for file in files:
                file_path = os.path.join(root, file)
                # Save the file in the zip with only its base name
                zip_ref.write(file_path, os.path.basename(file))


def random_string(length):
    return "".join(random.choices(string.ascii_letters, k=length))
