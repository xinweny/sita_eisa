import argparse
import pandas as pd
import os, glob, sys

def main():
    parser = argparse.ArgumentParser(description='Generate QuasR-format sample file.')
    parser.add_argument('-e', required=True,
                        help="Extension (e.g. fq, fastq)")

    args = parser.parse_args()
    ext = args.e
    metadata_path = args.m

    metadata_df = pd.read_csv(metadata_path, header=0, sep=',')
    is_paired = True if metadata_df['LibraryLayout'][0] == "PAIRED" else False

    if is_paired:
        sample_file = pd.DataFrame(columns=['FileName1', 'FileName2', 'SampleName'])

        pair1_files = [file for file in glob.glob(f"fastq/*_1.{ext}")]
        pair2_files = [file for file in glob.glob(f"fastq/*_2.{ext}")]

        if len(pair1_files) != len(pair2_files):
            sys.exit("Not all samples found in pairs.")

        pair1_files.sort()
        pair2_files.sort()

        sample_names = [file.split("/")[-1][:-(3 + len(ext))] for file in pair1_files]

        sample_file['FileName1'] = pair1_files
        sample_file['FileName2'] = pair2_files
        sample_file['SampleName'] = sample_names

    else:
        sample_file = pd.DataFrame(columns=['FileName', 'SampleName'])

        files = [file for file in glob.glob(f"fastq/*.{ext}")]
        files.sort()

        sample_names = [file.split("/")[-1][:-(1 + len(ext))] for file in files]

        sample_file['FileName'] = files
        sample_file['SampleName'] = sample_names

    sample_file.to_csv('SampleFile.txt', header=True, index=False, sep='\t')

if __name__ == "__main__":
    main()
