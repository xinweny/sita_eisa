#### Packages ####
import argparse
import os
import re
import pandas as pd
import glob

#### Functions ####
def aspera_download(err, metadata, paired=False):
    name = metadata.loc[metadata['run_accession'] == err, 'sample_title'].iloc[0]

    if paired:
        asperas = metadata.loc[metadata['run_accession'] == err, 'fastq_aspera'].iloc[0].split(';')

        sys_value = 0

        for i, aspera in enumerate(asperas):
            sys_val = os.system(f"ascp -QT -l 1000m -P33001 -i ~/.aspera/cli/etc/asperaweb_id_dsa.openssh \
                        era-fasp@{aspera} \
                        ./{err}_{name}_{i + 1}.fq.gz")

            sys_value += sys_val
    else:
        aspera = metadata.loc[metadata['run_accession'] == err, 'fastq_aspera'].iloc[0]

        sys_val = os.system(f"ascp -QT -l 1000m -P33001 -i ~/.aspera/cli/etc/asperaweb_id_dsa.openssh \
                    era-fasp@{aspera} \
                    ./{err}_{name}.fq.gz")

        sys_value = sys_val

    return sys_value

def is_paired(metadata, element):
    return metadata.loc[metadata['run_accession'] == element, 'library_layout'].iloc[0] == "PAIRED"

def check_exitcodes(sys_values):
    if sys_values == {}:
        return True

    if sum(sys_values.values()) > 0:
        print("WARNING: Some files not fully downloaded:")
        for err, value in sys_values.items():
            if value != 0:
                print(f"{err}, STATUS: {value}")

        return False
    else:
        print("All files downloaded completely.")

        return True

def main():
    # Parser
    parser = argparse.ArgumentParser(description='Download fq.gz files from ENA')
    parser.add_argument('-m', required=True, help='Path to ENA run table metadata.')

    args = parser.parse_args()
    metadata_path = args.m

    metadata_df = pd.read_csv(metadata_path, header=0, sep='\t')
    run_accs = metadata_df['run_accession']

    # Set up file architecture
    os.system("mkdir fastq")
    os.chdir("fastq")

    not_all_downloaded = True
    sys_values = {}

    while not_all_downloaded:
        # Download .fq.gz with aspera
        for i, acc in enumerate(run_accs):
            print(f"({i + 1}/{len(run_accs)}) Processing {acc}...")

            if len([file for file in glob.glob(f"{acc}*.fq.gz")]) > 0:
                print(f"{acc}*.fq.gz already downloaded. Skipping...")
                continue

            print(f"Downloading {acc}...")
            sys_value = aspera_download(acc, metadata_df, paired=is_paired(metadata_df, acc))

            sys_values[acc] = sys_value

            not_all_downloaded = not check_exitcodes(sys_values)

#### Execute code ####
if __name__ == "__main__":
    main()
