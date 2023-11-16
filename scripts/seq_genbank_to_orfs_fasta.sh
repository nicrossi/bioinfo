#!/bin/bash

# Check arguments
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <accession_id>"
    exit 1
fi

gene_name="$1"
# Fetch genbank, store the result to use as input later
python fetch_genbank.py "$gene_name"
# Run genbank to fasta on success
if [ $? -eq 0 ]; then
    echo "fetch_genbank.py successfully executed"
    python genbank_to_fasta.py "$gene_name.gbk"
else
    echo "Error: fetch_genbank.py failed"
    exit 1
fi
