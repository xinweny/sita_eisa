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

    os.system("mkdir bam processed cache fastq")

    with open(str(input_path), 'r') as acc_list:
        srrs = acc_list.read().splitlines()

    metadata_df = pd.read_csv(metadata_path, header=0, sep=',')

    os.chdir("fastq")

    for i, srr in enumerate(srrs):
        print(f"({i + 1}/{len(srrs)}) Downloading {srr}...")
        os.system(f"prefetch -v -o ./{srr}.sra {srr}")

        print(f"Converting {srr} to .fastq...")
        library_layout = metadata_df.loc[metadata_df['Run'] == srr, 'LibraryLayout'].iloc[0]
        sample_name = metadata_df.loc[metadata_df['Run'] == srr, 'Sample Name'].iloc[0]

        if library_layout == "PAIRED":
            os.system(f"fastq-dump --split-files {srr}.sra --skip-technical --readids")
            os.system(f"mv -v {srr}_1.fastq {sample_name}_1.fastq")
            os.system(f"mv -v {srr}_2.fastq {sample_name}_2.fastq")
        else:
            os.system(f"fastq-dump {srr}.sra --skip-technical --readids")
            os.system(f"mv -v {srr}.fastq {sample_name}.fastq")

        os.system(f"rm {srr}.sra")

# Execute code
if __name__ == "__main__":
    main()
