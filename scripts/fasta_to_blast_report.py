import sys
import subprocess
import threading
import time

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
    cmd = f"blastp -query {input_file} -db {db_name} -out {out_name}.out -outfmt '6 std qlen slen'"

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
        return 0
    else:
        return -1



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
