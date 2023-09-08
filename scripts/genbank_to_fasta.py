import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Validate command line arguments
if len(sys.argv) != 2:
    print("------------------------------------")
    print("Genbank to FASTA translation script")
    print("------------------------------------\n")
    print("[Use]    python genbank_to_fasta.py file_to_translate.gbk\n")
    print("[Input]  Genbank (.gbk) file")
    print("[Output] FASTA (.fas) file")
    sys.exit(1)

# Input and Output files
input_filename = sys.argv[1]
output_filename = input_filename.replace(".gbk", ".fas")

# Read GenBank and translate to FASTA
with open(input_filename, "r") as input_file, open(output_filename, "w") as output_file:
    for record in SeqIO.parse(input_file, "genbank"):
        for feature in record.features:
            if feature.type == "CDS":
                seq = feature.qualifiers["translation"][0]
                seq_obj = Seq(seq)
                # Using gene name (if available) as identifier
                identifier = feature.qualifiers.get("gene", ["unknown_gene"])[0]
                seq_record = SeqRecord(seq_obj, id=f"{record.id}_{identifier}")
                SeqIO.write(seq_record, output_file, "fasta")

print(f"GenBank to FASTA translation complete, output file: {output_filename}")
