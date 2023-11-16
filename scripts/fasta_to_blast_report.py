import os
import sys
import subprocess
import threading
import time
import re


# Run child process to execute the BLAST search
def run_command(cmd):
    result = subprocess.run(cmd, shell=True, text=True, capture_output=True)
    return result


# Print dots to give the user some feedback while searching
def print_dots(exit_signal):
    while not exit_signal.is_set():
        print(".", end="", flush=True)
        time.sleep(2)


# Executes the BLAST search in the provided database and print result on a .out file
def to_blast(input_file, db_name, out_name):
    longest = get_longest_orf(input_file)
    cmd = f"blastp -query {longest} -db {db_name} -out {out_name}.out -outfmt '6 std qlen slen'"

    print("Fasta to Blast report:")
    print(f"Processing input {input_file}")
    print(f"Searching in {db_name} (this could take a while)")

    exit_signal = threading.Event()
    dots_thread = threading.Thread(target=print_dots, args=(exit_signal,))
    dots_thread.start()

    result = run_command(cmd)
    if result.returncode == 0:
        exit_signal.set()
        dots_thread.join()
        clean_tmp_file(longest)
        return 0
    else:
        return -1


def extract_len_from_header(header):
    match = re.search(r'len (\d+)', header)
    if match:
        return int(match.group(1))
    else:
        return None


# Clean up temporary file
def clean_tmp_file(tmp_filename):
    try:
        os.remove(tmp_filename)
    except Exception as e:
        print(f"Error while removing temporary file: {e}")


def get_longest_orf(fastafile):
    """
    Reads a fasta file and save the longest ORF
    translation in a temp fasta file to use on the BLAST report
    """
    sequences = []
    with open(fastafile, "r") as f:
        lines = f.readlines()
        for line in lines:
            sequences.append(line.rstrip("\n"))

    max_len = 0
    max_seq_id = None
    seq_id = []
    seq_dic = {}

    for i in range(len(sequences)):
        if sequences[i][0] == ">":
            seq_id.append(sequences[i])
            curr_len = extract_len_from_header(sequences[i])
            if curr_len > max_len:
                max_len = curr_len
                max_seq_id = sequences[i]

    if max_seq_id:
        seq_id_index = [sequences.index(seq_id[i]) for i in range(len(seq_id))]
        for i in range(len(seq_id_index)):
            if i == (len(seq_id_index) - 1):
                seq_dic[seq_id[i]] = "".join(sequences[seq_id_index[i] + 1:])
            else:
                seq_dic[seq_id[i]] = "".join(sequences[seq_id_index[i] + 1:seq_id_index[i + 1]])

        # Create a temporary FASTA file with the longest sequence
        temp_fasta_file = "temp_longest_orf.fasta"
        with open(temp_fasta_file, "w") as temp_fasta:
            temp_fasta.write(f">{max_seq_id}\n{seq_dic[max_seq_id]}")

        return temp_fasta_file
    else:
        raise Exception(f"Error while reading input file {fastafile}")


def main():
    if len(sys.argv) != 3 or not sys.argv[1].lower().endswith((".fasta", ".fas")):
        print("FAILED")
        print("----------------------------")
        print("Fasta to Blast report tool")
        print("----------------------------\n")
        print("[Use]     python fasta_to_blast_report.py fasta_file.fas path_to_db")
        print("[Example] python fasta_to_blast_report.py fasta_file.fas ../../ncbi-blast-2.14.1+/data/swissprot")
        print("[Output]  .out file containing the BLAST search results")
        print("FAILED")

        if not sys.argv[1].lower().endswith((".fas", ".fasta")):
            print("File must be of Fasta(.fas) type")

        sys.exit(1)

    input_fasta = sys.argv[1]
    out_name = input_fasta.rsplit(".", 1)[0]
    if to_blast(input_fasta, sys.argv[2], out_name) == 0:
        print("\nFasta to Blast report complete!")
        print(f"Check results in {out_name}.out")
    else:
        print(f"There was an error processing {input_fasta}")


if __name__ == "__main__":
    main()
