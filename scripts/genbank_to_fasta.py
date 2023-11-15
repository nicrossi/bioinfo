import sys
from Bio import SeqIO
from Bio.Seq import translate, Seq
from Bio.SeqUtils import six_frame_translations, nt_search


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
    with open(input_filename, "r") as input_file, open(
        output_filename, "w"
    ) as output_file:
        for record in SeqIO.parse(input_file, "genbank"):
            print(f"Dealing with GenBank record {record.id}")
            seq = record.seq
            for orf in find_orfs(seq, record.id):
                output_file.write(f"{orf}\n")

    print(f"GenBank to FASTA translation complete, output file: {output_filename}")


def find_orfs(sequence, record_id):
    """
    Find Open Reading Frames (ORFs) in a DNA sequence.

    Parameters:
    - dna_sequence (sequence, record_id): DNA sequence, Sequence Record id

    Returns:
    - List of translation of all 6 possible ORFs in fasta format
    """

    def process_orf(start_pos, frame_sign):
        orf = sequence[start_pos:].translate(to_stop=True)
        orf_len = len(orf)
        if orf_len >= 30:
            frame = frame_sign * ((abs(start_pos) % 3) + 1)
            seq_trans = "\n".join(str(orf[i : i + 60]) for i in range(0, len(orf), 60))
            frames.append(
                f">{record_id} | Translation ORF {frame:+} | {start_pos}/{start_pos + orf_len} | len {orf_len}\n"
                + str(seq_trans)
            )

    start_positions = list(nt_search(str(sequence), "ATG"))
    start_positions.pop(0)
    reverse_starts = list(nt_search(str(sequence.reverse_complement()), "ATG"))
    reverse_starts.pop(0)

    frames = []

    for start_pos in start_positions:
        process_orf(start_pos, 1)

    for start_pos in reverse_starts:
        process_orf(start_pos, -1)

    return frames


if __name__ == "__main__":
    main()
