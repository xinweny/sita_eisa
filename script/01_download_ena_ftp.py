#### Packages ####
import argparse
import os
import re
import pandas as pd
import glob

#### Functions ####
def get_acc_names(metadata):
    acc_names = {}

    for acc in list(metadata['run_accession']):
        acc_names[acc] = metadata.loc[metadata['run_accession'] == acc, 'sample_title'].iloc[0]
    
    return acc_names

def aspera_download(acc, acc_names, metadata, paired, ext):
    name = acc_names[acc]

    if paired:
        ftps = metadata.loc[metadata['run_accession'] == acc, 'fastq_ftp'].iloc[0].split(';')

        sys_value = 0

        for i, ftp in enumerate(ftps):
            sys_val = os.system(f"wget -P ./fastq {ftp} && mv ./fastq/{acc}_{i + 1}.fastq.gz ./fastq/{name}_{i + 1}.{ext}")

            sys_value += sys_val
    else:
        ftp = metadata.loc[metadata['run_accession'] == acc, 'fastq_ftp'].iloc[0]

        sys_val = os.system(f"wget -P ./fastq {ftp} && mv ./fastq/{acc}.fastq.gz ./fastq/{name}.{ext}")

        sys_value = sys_val

    return sys_value

def is_paired(metadata, element):
    return metadata.loc[metadata['run_accession'] == element, 'library_layout'].iloc[0] == "PAIRED"

def check_exitcodes(sys_values):
    if sys_values == {}:
        return True

    if sum(sys_values.values()) > 0:
        print("WARNING: Some files not fully downloaded:")
        for acc, value in sys_values.items():
            if value != 0:
                print(f"{acc}, STATUS: {value}")

        return False
    else:
        print("All files downloaded completely.")

        return True

def main():
    # Parser
    parser = argparse.ArgumentParser(description='Download fastq.gz files from ENA')
    parser.add_argument('-m', required=True, help='Path to ENA run table metadata.')
    parser.add_argument('-e', nargs='?', const='fq.gz', default='fq.gz', help='Extension name (default: fq.gz)')

    args = parser.parse_args()
    metadata_path = args.m
    ext = args.e

    metadata_df = pd.read_csv(metadata_path, header=0, sep='\t')

    acc_names = get_acc_names(metadata_df)

    # Set up file architecture
    os.system("mkdir fastq")

    not_all_downloaded = True
    sys_values = {}

    while not_all_downloaded:
        # Download .fq.gz with aspera
        for acc, name in acc_names.items():
            print(f"Processing {acc}...")
            num_files = 2 if is_paired(metadata_df, acc) else 1

            if len([file for file in glob.glob(f"./fastq/{name}*.{ext}")]) == num_files:
                print(f"{name}*.{ext} already downloaded. Skipping...")
                continue

            print(f"Downloading {acc}...")
            sys_value = aspera_download(acc=acc, acc_names=acc_names, metadata=metadata_df, paired=is_paired(metadata_df, acc), ext=ext)

            sys_values[acc] = sys_value

            not_all_downloaded = not check_exitcodes(sys_values)

#### Execute code ####
if __name__ == "__main__":
    main()
