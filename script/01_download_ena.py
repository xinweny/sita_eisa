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

def ftp_download(acc, metadata, paired, ext):
    if paired:
        ftps = metadata.loc[metadata['run_accession'] == acc, 'fastq_ftp'].iloc[0].split(';')

        sys_value = 0

        for ftp in ftps:
            sys_val = os.system(f"wget -P ./fastq {ftp}")
            
            if sys_val != 0:
                sys_value = os.WEXITSTATUS(sys_val)

    else:
        ftp = metadata.loc[metadata['run_accession'] == acc, 'fastq_ftp'].iloc[0]

        sys_val = os.system(f"wget -P ./fastq {ftp}")

        sys_value = os.WEXITSTATUS(sys_val)

    return sys_value

def aspera_download(acc, metadata, paired, ext):
    if paired:
        ascps = metadata.loc[metadata['run_accession'] == acc, 'fastq_aspera'].iloc[0].split(';')

        sys_value = 0

        for i, ascp in enumerate(ascps):
            sys_val = os.system(f"ascp -QT -l 1000m -P33001 -i ~/.aspera/cli/etc/asperaweb_id_dsa.openssh \
                        era-fasp@{ascp} \
                        ./fastq/{acc}_{i + 1}.{ext}")

            if sys_val != 0:
                sys_value = os.WEXITSTATUS(sys_val)

    else:
        ascp = metadata.loc[metadata['run_accession'] == acc, 'fastq_aspera'].iloc[0]

        sys_val = os.system(f"ascp -QT -l 1000m -P33001 -i ~/.aspera/cli/etc/asperaweb_id_dsa.openssh \
                    era-fasp@{ascp} \
                    ./fastq/{acc}.{ext}")

        sys_value = os.WEXITSTATUS(sys_val)

    return sys_value

def merge_accs(accs, name, paired, ext):
    if paired:
        for i in [1, 2]:
            accs_filenames = ' '.join([f"./fastq/{acc}_{i}.fastq.gz" for acc in accs])
            os.system(f"cat {accs_filenames} > ./fastq/{name}_{i}.{ext} && rm {accs_filenames}")
    else:
        accs_filenames = ' '.join([f"./fastq/{acc}.fastq.gz" for acc in accs])
        os.system(f"cat {accs_filenames} > ./fastq/{name}.{ext} && rm {accs_filenames}")

def rename_file(acc, name, paired, ext):
    if paired:
        for i in [1, 2]:
            os.system(f"mv ./fastq/{acc}_{i}.fastq.gz ./fastq/{name}_{i}.{ext}")
    else:
        os.system(f"mv ./fastq/{acc}.fastq.gz ./fastq/{name}.{ext}")

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
    parser.add_argument('-e', nargs='?', const='fastq.gz', default='fastq.gz', help='Extension name (default: fastq.gz)')
    parser.add_argument('-d', nargs='?', const='ascp', default='ascp', help='Download type (default: ascp)')

    args = parser.parse_args()
    metadata_path = args.m
    ext = args.e
    download_type = args.d

    metadata_df = pd.read_csv(metadata_path, header=0, sep='\t')

    acc_names = get_acc_names(metadata_df)

    # Set up file architecture
    os.system("mkdir fastq")

    not_all_downloaded = True
    sys_values = {}

    while not_all_downloaded:
        # Download .fq.gz with aspera
        for i, name in enumerate(set(acc_names.values())):
            print(f"({i + 1}/{len(set(acc_names.values()))}) Processing {name}.{ext}...")

            accs = [acc for acc, n in acc_names.items() if n == name]

            paired = is_paired(metadata_df, accs[0])
            num_files = 2 if paired else 1

            format_name = re.sub(r'[\/:,]', "", name).replace(' ', '_')

            if len([file for file in glob.glob(f"./fastq/{format_name}*.{ext}")]) == num_files:
                    print(f"{format_name}*.{ext} already downloaded. Skipping...")
                    continue

            name_sys_value = 0

            for j, acc in enumerate(accs):
                if len([file for file in glob.glob(f"./fastq/{acc}*.{ext}")]) == num_files:
                    print(f"{acc}*.{ext} already downloaded. Skipping...")
                    continue

                print(f"({j + 1}/{len(accs)}) Downloading {acc}...")
                if download_type == 'ascp':
                    sys_value = aspera_download(acc=acc, metadata=metadata_df, paired=paired, ext=ext)
                else:
                    sys_value = ftp_download(acc=acc, metadata=metadata_df, paired=paired, ext=ext)

                sys_values[acc] = sys_value
                name_sys_value += sys_value

            if name_sys_value == 0 and len(accs) > 1:
                print(f"Merging accessions for multiple technical runs of {format_name}.{ext}...")
                merge_accs(accs, name, paired, ext)
            elif name_sys_value == 0:
                print(f"Renaming to {name}.{ext}...")
                rename_file(acc, name, paired, ext)

        not_all_downloaded = not check_exitcodes(sys_values)

#### Execute code ####
if __name__ == "__main__":
    main()
