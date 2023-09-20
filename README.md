# bioinfo
Bioinformatica 2023 2Q 

# Requirements

-Python
-Blast+
-Muscle
-Blast+ databse already set up


# fetch_genbank 

Searchs for a gene ID on NCBI-Gene.

[Use]       python fetch_genbank.py GENE-ID
[Example]   python fetch_genbank.py FBN1

# genbank_to_fasta 

Translate a .gbk to a .fas file 

[Use]       python genbank_to_fasta.py file_to_translate.gbk
[Example]   python genbank_to_fasta.py gene=FBN1-NM_000138-20230919203303.gbk

# fasta_to_blast_report 

Perfoms a Blast search with a .fas imput file

[Use]       python fasta_to_blast_report.py fasta_file.fas path_to_db
[Example]   python fasta_to_blast_report.py gene=FBN1-NM_000138-20230919203303.fas ../../ncbi-blast-2.14.1+/data/swissprot