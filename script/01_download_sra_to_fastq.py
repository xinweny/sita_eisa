import argparse
import os
import pandas as pd

def mkdir(directory):
    if not os.path.exists(directory):
        os.system(f"mkdir -p {directory}")

    return directory

# ftp://ftp.sra.ebi.ac.uk/vol1/srr/SRR602/002/SRR6026772 - URL architecture for download with wget

def main():
    parser = argparse.ArgumentParser(description='Script to download SRA files')
    parser.add_argument('-i', required=True, help='path to SRR accession list')
    parser.add_argument('-m', required=True, help='path to SRA run table metadata')

    args = parser.parse_args()
    input_path = args.i
    metadata_path = args.m

    os.system("mkdir bam processed cache fastq script")

    with open(str(input_path), 'r') as acc_list:
        srrs = acc_list.readlines()

    metadata_df = pd.read_csv(metadata_path, header=True, sep=',')

    os.chdir("fastq")

    for i, srr in enumerate(srrs):
        print(f"({i + 1}/{len(srrs)}) Downloading {srr}...")
        os.system(f"prefetch -v {srr}")

        print(f"Converting {srr} to .fastq...")
        os.system(f"fastq-dump {srr}.sra --skip-technical --readids")

        sample_name = metadata_df.loc[metadata_df['Run'] == srr, 'Sample Name'].iloc[0]
        os.system(f"mv {srr}.fastq {sample_name}.fastq")

# Execute code
if __name__ == "__main__":
    main()
