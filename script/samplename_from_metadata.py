# Packages
import pandas as pd
import argparse
import os, sys
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

def main():
    parser = argparse.ArgumentParser(description='Script to convert sample names based on SRA Run Table metadata')
    parser.add_argument('-g', required=True, help='GSE accession number')
    parser.add_argument('-c', required=True, help='columns in metadata table, separated by \",\"')

    args = parser.parse_args()
    gse = args.g
    columns = (args.c).split(',')

    metadata_df = pd.read_csv(f"{gse}/SraRunTable_{gse}.txt", header=0, sep=',')
    gsms = metadata_df['Sample Name'].unique()
    gsm_sample = scrape_sample_names(gse)

    os.chdir(f"{gse}/fastq")

    for gsm in gsms:
        info = [metadata_df.loc[metadata_df['Sample Name'] == gsm, col].iloc[0] for col in columns]
        new_name = '_'.join(info)

        if metadata_df.loc[metadata_df['Sample Name'] == gsm, 'LibraryLayout'].iloc[0] == "PAIRED":
            os.system(f"mv -v '{gsm}_{gsm_sample[gsm]}_1.fastq.gz' '{gsm}_{new_name}_1.fq.gz'")
            os.system(f"mv -v '{gsm}_{gsm_sample[gsm]}_2.fastq.gz' '{gsm}_{new_name}_2.fq.gz'")
        else:
            os.system(f"mv -v '{gsm}_{gsm_sample[gsm]}.fastq.gz' '{gsm}_{new_name}.fq.gz'")

# Execute code
if __name__ == "__main__":
    main()
