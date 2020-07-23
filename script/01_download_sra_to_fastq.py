import argparse
import os

def mkdir(directory):
    if not os.path.exists(directory):
        os.system(f"mkdir -p {directory}")

    return new_dir

# ftp://ftp.sra.ebi.ac.uk/vol1/srr/SRR602/002/SRR6026772 - URL architecture for download with wget

def main():
    parser = argparse.ArgumentParser(description='Script to download SRA files')
    parser.add_argument('-i', required=True, help='path to SRR accession list')
    parser.add_argument('-o', default='./', help='output directory path')

    args = parser.parse_args()
    outdir = mkdir(args.o)

    with open(str(args.i), 'r') as acc_list:
        srrs = acc_list.readlines()

    os.chdir(outdir)

    for i, srr in enumerate(srrs):
        print(f"({i + 1}/{len(srrs)}) Downloading {srr}...")
        os.system(f"prefetch -v {srr}")

        print(f"Converting {srr} to .fastq...")
        os.system(f"fastq-dump {srr}.sra --skip-technical --readids")

# Execute code
if __name__ == "__main__":
    main()
