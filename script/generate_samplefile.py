import argparse
import pandas as pd
import os, glob, sys

def make_samplefile(ext, layout):
    if layout == "PAIRED":
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

        files = [file for file in glob.glob(f"fastq/*.{ext}") if not any(end in file for end in ['_1.fastq', '_2.fastq'])]
        files.sort()

        sample_names = [file.split("/")[-1][:-(1 + len(ext))] for file in files]

        sample_file['FileName'] = files
        sample_file['SampleName'] = sample_names

    return sample_file

def main():
    parser = argparse.ArgumentParser(description='Generate QuasR-format sample file.')
    parser.add_argument('-m', required=True,
                        help="Path to SraRunTable")
    parser.add_argument('-e', required=True,
                        help="Extension (e.g. fq, fastq)")

    args = parser.parse_args()
    ext = args.e
    metadata_path = args.m

    metadata_df = pd.read_csv(metadata_path, header=0, sep=',')
    layout = list(metadata_df['LibraryLayout'].unique())

    if len(layout) == 1:
        sample_file = make_samplefile(ext, layout[0])

        sample_file.to_csv(f"SampleFile_{layout[0]}.txt", header=True, index=False, sep='\t')
    else:
        for l in layout:
            sample_file = make_samplefile(ext, l)
            sample_file.to_csv(f"SampleFile_{l}.txt", header=True, index=False, sep='\t')


if __name__ == "__main__":
    main()
