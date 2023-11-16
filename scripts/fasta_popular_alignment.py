import os

import requests
from Bio import SeqIO, Entrez
import subprocess
import sys


# Gets the n-Best BLAST search results
def get_best_ids(blast_results_file, n=10):
    top_accession = {}
    with open(blast_results_file, "r") as blast_file:
        for line_n, line in enumerate(blast_file):
            if len(top_accession) == n:
                break
            fields = line.strip().split("\t")
            if len(fields) >= 2 and not fields[1] in top_accession:
                top_accession[fields[1]] = fields
    return top_accession


def fetch_and_concat_sequences(accession_ids, query_seq, local_db=None):
    data_sequences = []
    if query_seq:
        data_sequences.append(efetch_nucleotide_fasta(query_seq))

    # Fetch and store protein sequences for top results
    for accession_id in accession_ids.keys():
        if local_db:
            try:
                # Try fetching from provided database
                data_sequences.append(dbfetch_prot_sequence(local_db, accession_id))
            except subprocess.CalledProcessError as e:
                print(f"Error fetching protein with accession_id: {accession_id} from local db {local_db}:\n {e}")
                data_sequences.append(efetch_prot_sequence(accession_id))
        else:
            # Fetch using Entrez API in case DB fails
            data_sequences.append(efetch_prot_sequence(accession_id))

    output_file = "msa_input_sequences.fasta"
    SeqIO.write(data_sequences, output_file, "fasta")
    return output_file


# Fetch protein sequence with id 'accession_id' using Entrez API
def efetch_prot_sequence(accession_id):
    try:
        print(f"Entrez.efetch(db=\"protein\", id={accession_id}, rettype=\"fasta\", retmode=\"text\")")
        handle = Entrez.efetch(db="protein", id=accession_id, rettype="fasta", retmode="text")
        seq = SeqIO.read(handle, "fasta")
        handle.close()
        return seq
    except Exception as e:
        print(f"Error fetching protein with accession_id {accession_id} using Entrez.efetch: {e}")


# Fetch sequence with id 'accession_id' from database in 'db_path'
def dbfetch_prot_sequence(db_path, accession_id):
    cmd = f"blastdbcmd -db {db_path} -entry {accession_id} -outfmt %f"
    seq = subprocess.check_output(cmd, shell=True, text=True)
    seq = SeqIO.read(seq.splitlines(), "fasta")
    return seq


# Function to perform a multiple sequence alignment MUSCLE
def perform_msa_muscle(input_fasta, output_msa):
    cmd = f"muscle -in {input_fasta} -out {output_msa}"
    subprocess.run(cmd, shell=True)


# Clean up temporary file
def clean_tmp_file(tmp_filename):
    try:
        os.remove(tmp_filename)
    except Exception as e:
        print(f"Error while removing temporary file: {e}")


def perform_msa_clustalo(msa_input, msa_output, email):
    clustalo_script = "clustalo.py"
    command = [
        "python",
        clustalo_script,
        "--email", email,
        "--sequence", msa_input,
        "--outfile", msa_output
    ]
    print("=== Multiple Sequence Alignment ===")
    print(f"input_file: {msa_input}")
    print("Clustal Omega is running, please wait for SUCCESS message...")
    # Execute the clustalo script
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # Wait for the command to complete and capture the output
    stdout, stderr = process.communicate()
    # Check the return code to see if the command was successful
    if process.returncode == 0:
        print("[SUCCESS] Clustal Omega command was successful.")
    else:
        print("[FAILURE] Clustal Omega command failed. Error message:")
        print(stderr.decode('utf-8'))


def efetch_nucleotide_fasta(accession_id):
    try:
        print(f'Entrez.efetch(db="nucleotide", id={accession_id}, rettype="fasta", retmode="text")')
        handle = Entrez.efetch(db="nucleotide", id=accession_id, rettype="fasta", retmode="text")
        record = SeqIO.read(handle, "fasta")
        handle.close()
        return record
    except Exception as e:
        print(f"Error fetching gene with accession_id {accession_id} using Entrez.efetch: {e}")


def main():
    if len(sys.argv) < 3 or not sys.argv[2].lower().endswith(".out"):
        print("Invalid arguments")
        print("---------------------------------------------")
        print("Alignment of top 10 results of a BLAST Report")
        print("-------------------------------------------\n")
        print("[Use]     python msa.py query_sequence_path blast_report.out email path_to_db")
        print("[Example] python msa.py /FBN1-NM_000138.fas blast_results.out example@email.com /ncbi-blast/db/swissprot")
        print("[Output]  Fasta file with the MSA result")

        if not sys.argv[2].lower().endswith(".out"):
            print("File must be of (.out) type")

        sys.exit(1)

    query_sequence = sys.argv[1]
    blast_results_file = sys.argv[2]
    email = sys.argv[3]
    Entrez.email = email
    db_path = None if len(sys.argv) < 5 else sys.argv[4]
    msa_output = blast_results_file.rsplit('.', 1)[0] + "_top_MSA"

    # Get 10 accession ids for best results
    accession_ids_map = get_best_ids(blast_results_file)
    # Fetch protein sequences and merge them to query sequence inta a fasta file
    msa_input = fetch_and_concat_sequences(accession_ids_map, query_sequence, None if db_path is None else db_path)
    # perform MSA
    perform_msa_clustalo(msa_input, msa_output, email)
    print(f"Check results in {msa_output}.out")
    clean_tmp_file(msa_input)


if __name__ == '__main__':
    main()
