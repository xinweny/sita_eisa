import argparse
import pandas as pd
import os, glob, sys

def make_samplefile(ext, layout):
    if layout == "PAIRED":
        sample_file = pd.DataFrame(columns=['FileName1', 'FileName2', 'SampleName'])

        pair1_files = [file for file in glob.glob(f"fastq_trimmed/*_1_val_1.{ext}")]
        pair2_files = [file for file in glob.glob(f"fastq_trimmed/*_2_val_2.{ext}")]

        if len(pair1_files) != len(pair2_files):
            sys.exit("Not all samples found in pairs.")

        pair1_files.sort()
        pair2_files.sort()

        sample_names = [file.split("/")[-1][:-(9 + len(ext))] for file in pair1_files]

        sample_file['FileName1'] = pair1_files
        sample_file['FileName2'] = pair2_files
        sample_file['SampleName'] = sample_names
    else:
        sample_file = pd.DataFrame(columns=['FileName', 'SampleName'])

        files = [file for file in glob.glob(f"fastq_trimmed/*_trimmed.{ext}") if not any(end in file for end in ['_1.fastq', '_2.fastq'])]
        files.sort()

        sample_names = [file.split("/")[-1][:-(9 + len(ext))] for file in files]

        sample_file['FileName'] = files
        sample_file['SampleName'] = sample_names

    return sample_file

def main():
    parser = argparse.ArgumentParser(description='Generate QuasR-format sample file.')
    parser.add_argument('-e', required=True,
                        help="Extension (e.g. fq, fastq, fastq.gz)")
    parser.add_argument('-l', required=True,
                        help="Library layout (SINGLE or PAIRED)")
    parser.add_argument('-w', required=True,
                        help="Working directory")

    args = parser.parse_args()
    ext = args.e

    os.chdir(args.w)

    layout = args.l.upper()

    sample_file = make_samplefile(ext, layout)

    sample_file.to_csv(f"SampleFile_{layout}.txt", header=True, index=False, sep='\t')

main()
