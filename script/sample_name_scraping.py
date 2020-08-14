from requests import get
from bs4 import BeautifulSoup
import argparse

parser = argparse.ArgumentParser(description='Script to check GSM-sample name table.')
parser.add_argument('-g', required=True, help='GSE number')

args = parser.parse_args()

url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={args.g}"
response = get(url)

html_soup = BeautifulSoup(response.text, 'html.parser')
tables = html_soup.findAll('table')[18:20]

gsm_sample = {}

for table in tables:
    for i, row in enumerate(table.findAll("tr")):
        cells = row.findAll("td")
        text = [cell.text for cell in cells]
        gsm_sample[text[0]] = text[1]

print(f"No. of samples: {len(gsm_sample)}")

for gsm, name in gsm_sample.items():
    print(f"{gsm} - {name}")
