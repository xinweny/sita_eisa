import argparse
import os

def main():
    parser = argparse.ArgumentParser(description='Script to download SRA files')
    parser.add_argument('-i', help='path to SRR accession list')
    args = parser.parse_args()

    with open(str(args.i), 'r') as acc_list:
        srrs = acc_list.readlines()

    for i, srr in enumerate(srrs):
        print(f"({i + 1}/{len(srrs)}) Downloading {srr}...")
        os.system(f"prefetch -v {srr}")

        print(f"Converting {srr} to .fastq...")
        os.system(f"fastq-dump {srr}.sra --skip-technical --readids")

# Execute code
if __name__ == "__main__":
    main()
