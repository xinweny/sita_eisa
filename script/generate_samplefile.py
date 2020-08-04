import argparse
import pandas as pd
import os, glob, sys

def main():
    parser = argparse.ArgumentParser(description='Generate QuasR-format sample file.')
    parser.add_argument('--paired', default=False, action='store_true',
                        help="""Flag to indicate paired-end RNAseq data.
File names for each read should end in _1 and _2 for each paired sample.""")
    parser.add_argument('--single', dest='paired', action='store_false',
                        help="""Flag to indicate single-read RNAseq data.""")
    parser.add_argument('-e',
                        help="Extension (e.g. fq, fastq)")

    args = parser.parse_args()
    is_paired = args.paired
    ext = args.e

    cwd = os.getcwd()

    if is_paired:
        sample_file = pd.DataFrame(columns=['FileName1', 'FileName2', 'SampleName'])

        pair1_files = [file for file in glob.glob(f"*_1.{ext}")]
        pair2_files = [file for file in glob.glob(f"*_2.{ext}")]

        if len(pair1_files) != len(pair2_files):
            sys.exit("Not all samples found in pairs.")

        pair1_files.sort()
        pair2_files.sort()

        sample_names = [file[:-(3 + len(ext))] for file in pair1_files]

        sample_file['FileName1'] = pair1_files
        sample_file['FileName2'] = pair2_files
        sample_file['SampleName'] = sample_names

    else:
        sample_file = pd.DataFrame(columns=['FileName', 'SampleName'])

        files = [file for file in glob.glob(f"*.{ext}")]
        files.sort()

        sample_names = [file[:-(1 + len(ext))] for file in files]

        sample_file['FileName'] = files
        sample_file['SampleName'] = sample_names

    sample_file.to_csv('SampleFile.txt', header=True, index=False, sep='\t')

if __name__ == "__main__":
    main()
