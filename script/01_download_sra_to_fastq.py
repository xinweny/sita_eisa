import argparse
import os
import pandas as pd
import glob
from requests import get
from bs4 import BeautifulSoup

def scrape_sample_names(gse):
    url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse}"
    response = get(url)

    html_soup = BeautifulSoup(response.text, 'html.parser')
    tds = [td.text for td in html_soup.findAll("td")]
    tables = html_soup.findAll('table')

    if 'NIH grant(s)' in tds:
        start = 19
        end = 21
    else:
        start = 18
        end = 20

    sample_table = tables[start:end]

    gsm_sample = {}

    for table in sample_table:
        for i, row in enumerate(table.findAll("tr")):
            cells = row.findAll("td")
            text = [cell.text for cell in cells]
            gsm_sample[text[0]] = text[1]

    return gsm_sample

def is_paired(metadata, col_name, element):
    return metadata.loc[metadata[col_name] == element, 'LibraryLayout'].iloc[0] == "PAIRED"

# ftp://ftp.sra.ebi.ac.uk/vol1/srr/SRR602/002/SRR6026772 - URL architecture for download with wget
def main():
    parser = argparse.ArgumentParser(description='Script to download SRA files')
    parser.add_argument('-m', required=True, help='path to SRA run table metadata')

    args = parser.parse_args()
    metadata_path = args.m

    metadata_df = pd.read_csv(metadata_path, header=0, sep=',')
    gsms = metadata_df['Sample Name'].unique()

    gsm_samples = scrape_sample_names(os.getcwd().split("/")[-1])

    # Set up file architecture
    os.system("mkdir bam processed cache fastq")
    os.chdir("fastq")

    for i, gsm in enumerate(gsms):
        print(f"({i + 1}/{len(gsms)}) Processing {gsm}...")

        if len([file for file in glob.glob(f"{gsm}_{gsm_samples[gsm]}*.fastq")]) == 0:
            gsm_srrs = metadata_df.loc[metadata_df['Sample Name'] == gsm, 'Run']

            for srr in list(gsm_srrs):
                if len([file for file in glob.glob(f"{srr}*.fastq")]) == 0:
                    print(f"Downloading {srr}...")
                    os.system(f"prefetch -v -o ./{srr}.sra {srr}")

                    print(f"Converting {srr} to .fastq...")
                    if is_paired(metadata_df, 'Run', srr):
                        os.system(f"fastq-dump --split-files {srr}.sra --skip-technical --readids")
                    else:
                        os.system(f"fastq-dump {srr}.sra --skip-technical --readids")

                    os.system(f"rm {srr}.sra")
                else:
                    print(f"{srr}*.fastq already downloaded. Skipping...")

            if len(gsm_srrs) == 1:
                print(f"Renaming fastq files for {gsm}...")
                srr = gsm_srrs.iloc[0]

                if is_paired(metadata_df, 'Run', srr):
                    os.system(f"mv -v '{srr}_1.fastq' '{gsm}_{gsm_samples[gsm]}_1.fastq'")
                    os.system(f"mv -v '{srr}_2.fastq' '{gsm}_{gsm_samples[gsm]}_2.fastq'")
                else:
                    os.system(f"mv -v '{srr}.fastq' '{gsm_samples[gsm]}.fastq'")
            else:
                print(f"Merging technical runs and enaming fastq files for {gsm}...")
                if is_paired(metadata_df, 'Run', srr):
                    fwd_fqs = ' '.join(gsm_srrs + '_1.fastq')
                    rv_fqs = ' '.join(gsm_srrs + '_2.fastq')

                    os.system(f"cat {fwd_fqs} > '{gsm}_{gsm_samples[gsm]}_1.fastq' && rm {fwd_fqs}")
                    os.system(f"cat {rv_fqs} > '{gsm}_{gsm_samples[gsm]}_2.fastq' && rm {rv_fqs}")
                else:
                    srr_fqs = ' '.join(gsm_srrs + '.fastq')

                    os.system(f"cat {srr_fqs} > '{gsm}_{gsm_samples[gsm]}.fastq' && rm {srr_fqs}")

        else:
            print(f"{gsm}*.fastq already downloaded. Skipping...")

# Execute code
if __name__ == "__main__":
    main()
