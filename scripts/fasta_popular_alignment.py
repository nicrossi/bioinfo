from Bio import AlignIO, SeqIO
import subprocess
import sys


#Gets 10 best BLAST search results
def get_best_ids(blast_results_file, top_ids):
    with open(blast_results_file, "r") as blast_file:
        for line_n, line in enumerate(blast_file):
            if line_n >= 10:
                break
            fields = line.strip().split("\t")
            fasta_id = fields[0]
            top_ids.append(fasta_id)  # Use append to add IDs to the list


# Receives a list of fastas ids and performs an alignment
def align_fastas(top_ids, db_name, out_name):
    # Assuming you want to retrieve the sequences from the BLAST database and align them
    blastdbcmd_cmd = f"blastdbcmd -db {db_name} -entry {' '.join(top_ids)} > {out_name}"
    subprocess.run(blastdbcmd_cmd, shell=True)

    # Perform the multiple sequence alignment using MUSCLE
    muscle_cmd = f"muscle -in {out_name} -out {out_name}_aligned.fasta"
    subprocess.run(muscle_cmd, shell=True)


def main():
    if len(sys.argv) != 3 or not sys.argv[1].lower().endswith(".out"):
        print("FAILED")
        print("----------------------------")
        print("Alignment of top 10 results of a BLAST Report")
        print("----------------------------\n")
        print("[Use]     python fasta_to_blast_report.py blast_report.out path_to_db")
        print("[Example] python fasta_to_blast_report.py blast_results.out ../../ncbi-blast-2.14.1+/data/swissprot")
        print("[Output]  List of top sequence identifiers")
        print("FAILED")

        if not sys.argv[1].lower().endswith(".out"):
            print("File must be of (.out) type")

        sys.exit(1)

    blast_results_file = sys.argv[1]
    db_name = sys.argv[2]
    out_name = blast_results_file.rsplit('.', 1)[0] + "_top_10_fasta_alignment.out"  

    top_ids = []

    get_best_ids(blast_results_file, top_ids)
    align_fastas(top_ids, db_name, out_name)
    
    

if __name__ == '__main__':
    main()





"""
# Function to retrieve sequences from a local BLAST database
def retrieve_sequences_from_db(id_list, db_name, output_fasta):
    cmd = f"blastdbcmd -db {db_name} -entry {' '.join(id_list)} -out {output_fasta}"
    subprocess.run(cmd, shell=True)

# Function to perform a multiple sequence alignment
def perform_msa(input_fasta, output_msa):
    cmd = f"muscle -in {input_fasta} -out {output_msa}"
    subprocess.run(cmd, shell=True)

# List of FASTA sequence IDs you want to align
fasta_ids = ["NM_000138.5_FBN1", "P35555.4", "P98133.2", "Q9TV36.1"]  # Replace with your IDs

# Name of your local BLAST database
db_name = "path/to/your/local/database"

# Output filenames
output_fasta = "sequences_to_align.fasta"
output_msa = "alignment_result.fasta"

# Retrieve sequences from the local BLAST database
retrieve_sequences_from_db(fasta_ids, db_name, output_fasta)

# Perform multiple sequence alignment
perform_msa(output_fasta, output_msa)

print(f"Alignment result saved in {output_msa}")
"""