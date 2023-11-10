import sys
from Bio import SeqIO
from Bio.Seq import translate


def main():
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
            print(f"Dealing with GenBank record {record.id}")
            seq = record.seq
            orfs = generate_orfs(seq)
            for key, value in orfs.items():
                # Insert a line break every 60 chars (just for a pretty output)
                seq_trans = '\n'.join(str(value[i:i + 60]) for i in range(0, len(value), 60))
                output_file.write(f">{record.id} | frame {key}\n{seq_trans}\n")

    print(f"GenBank to FASTA translation complete, output file: {output_filename}")


# Generates the six possible frames translation per one sequence
def generate_orfs(seq):
    frames = {'+1': [], '+2': [], '+3': [], '-1': [], '-2': [], '-3': []}
    seq_rev = seq[::-1]
    for j in range(0, 3):
        seq_trans = translate(seq[j::])
        seq_rev_trans = translate(seq_rev[j::])
        if j == 0:
            frames['+1'] = seq_trans
            frames['-1'] = seq_rev_trans
        if j == 1:
            frames['+2'] = seq_trans
            frames['-2'] = seq_rev_trans
        if j == 2:
            frames['+3'] = seq_trans
            frames['-3'] = seq_rev_trans

    return frames


if __name__ == '__main__':
    main()
