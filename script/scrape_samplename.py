from requests import get
from bs4 import BeautifulSoup
import argparse

parser = argparse.ArgumentParser(description='Script to check GSM-sample name table.')
parser.add_argument('-g', required=True, help='GSE number')

args = parser.parse_args()

url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={args.g}"
response = get(url)

html_soup = BeautifulSoup(response.text, 'html.parser')
tables = html_soup.findAll('table')

if len(tables) == 24:
    start = 19
    end = 21
elif len(tables) == 23:
    start = 18
    end = 20

sample_table = tables[start:end]

gsm_sample = {}

for table in sample_table:
    for i, row in enumerate(table.findAll("tr")):
        cells = row.findAll("td")
        text = [cell.text for cell in cells]
        gsm_sample[text[0]] = text[1]

counter = 1

for gsm, name in gsm_sample.items():
    print(f"{counter}. {gsm} - {name}")
    counter += 1
