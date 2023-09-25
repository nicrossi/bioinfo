import argparse
import sys

import requests
import datetime
import xml.etree.ElementTree as ElementTree
import re

# Constants
BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
DATABASE = "gene"
ORGANISM = "Homo sapiens"


def search_gene_id(gene_name):
    search_url = f"{BASE_URL}esearch.fcgi?db={DATABASE}&term={gene_name}&retmax=10&sort=relevance&organism={ORGANISM}"

    try:
        response = requests.get(search_url)
        if response.status_code == 200:
            # Parse XML response to get the Gene ID
            root = ElementTree.fromstring(response.text)
            gene_id = root.find(".//Id").text

            return gene_id
        else:
            print(f"[search_gene_id] HTTP request. Status code: {response.status_code}")
    except requests.exceptions.RequestException as e:
        print(f"Request error: {e}")
    except ElementTree.ParseError as e:
        print(f"XML parsing error: {e}")


def fetch_genbank_file(gene, acc_id):
    genbank_url = f"{BASE_URL}efetch.fcgi?db=nucleotide&id={acc_id}&rettype=gb&retmode=text"

    response = requests.get(genbank_url)
    if response.status_code == 200:
        # Get the current timestamp
        timestamp = datetime.datetime.now().strftime("%Y%m%d%H%M%S")
        # Save content to local GenBank file
        with open(f"{gene}-{acc_id}-{timestamp}.gbk", "wb") as genbank_file:
            genbank_file.write(response.content)
        print(f"GenBank file downloaded as {gene}-{acc_id}-{timestamp}.gbk")
    else:
        print(f"[fetch_genbank_file] HTTP request. Status code: {response.status_code}")


def fetch_accession_id(gene_id):
    fetch_url = f"{BASE_URL}efetch.fcgi?db=gene&id={gene_id}"

    response = requests.get(fetch_url)
    if response.status_code == 200:
        return parse_accession_id(response.text)
    else:
        print(f"[fetch_accession_id] HTTP request. Status code: {response.status_code}")


# Parse mRNA accession_id (first match) from gene search query response
def parse_accession_id(data):
    pattern = r'accession "(NM_\d+)"'

    match = re.search(pattern, data)
    if match:
        accession = match.group(1)
        print(f"mRNA --> Accession_id: {accession}")
        return accession
    else:
        print("[parse_accession_id] No match found for accession_id with format (NM_\\d+)")


def main():
    if len(sys.argv) != 2:
        print("----------------------------")
        print("Genbank (.gbk) download tool")
        print("----------------------------\n")
        print("[Use]     python fetch_genbank.py GENE-ID")
        print("[Example] python fetch_genbank.py FBN1\n")
        print("[Output]  Genbank (.gbk) file")
        sys.exit(1)

    print(f"Querying NCBI-Gene Database...")
    parser = argparse.ArgumentParser(description="Search for a gene ID on NCBI-Gene.")
    parser.add_argument("search_term", type=str, help="Gene name to search for")

    args = parser.parse_args()
    print(f"searching for {args.search_term}...")

    try:
        gene_id = search_gene_id(args.search_term)
        if gene_id:
            print(f"Gene: {args.search_term} --> Id: {gene_id}")
            accession_id = fetch_accession_id(gene_id)
            fetch_genbank_file(args.search_term, accession_id)
    except Exception as e:
        print(f"Error: {e}")


if __name__ == '__main__':
    main()
