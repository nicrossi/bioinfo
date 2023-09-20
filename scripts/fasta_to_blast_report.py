import sys
import subprocess

#PROGRAMA QUE EJECUTA UNA BUSQUEDA BLAST CON IMPUT UN ARCHIVO FASTA

def to_blast(input_file, db_name, out_name):
    cmd = f"blastp -query {input_file} -db {db_name} -out {out_name}.out -outfmt '6 std qlen slen'"
    
    print("Ejecutando la busqueda en BLAST (puede tardar unos minutos)")

    subprocess.run(cmd, shell=True)
    print(f"BLAST para {input_file} completo. Resultados en {out_name}.out")

def main():

    if (len(sys.argv) != 3 or not sys.argv[1].lower().endswith((".fasta", ".fas"))):
        print("FAILED")
        print("----------------------------")
        print("Fasta to Blast report tool")
        print("----------------------------\n")
        print("[Use]     python fasta_to_blast_report.py fasta_file.fas path_to_db")
        print("[Example] python fasta_to_blast_report.py fasta_file.fas ../../ncbi-blast-2.14.1+/data/swissprot")
        print("[Output]  .out file containing the BLAST search results")
        print("FAILED")
        
        if(not sys.argv[1].lower().endswith((".fas", ".fasta"))):
            print("File must be of Fasta(.fas) type")

        sys.exit(1)
    

    input_fasta = sys.argv[1]
    out_name = input_fasta.rsplit('.', 1)[0]   

    

    to_blast(input_fasta, sys.argv[2], out_name)

if __name__ == '__main__':
    main()
