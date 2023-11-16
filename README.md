# ITBA BioInformática
Trabajo Práctico 2023 2Q 

# Requirements
-Python \
-Blast+ \
-Muscle \
-Blast+ database already set up

-Emboss 

# SCRIPTS

## 1. Bash script for full translation into a fasta file
Fetch genbank file, translate all ORFs and writes result into a FASTA file.
<pre>
[Use]       python genbank_to_fasta.py file_to_translate.gbk
[Example]   python genbank_to_fasta.py FBN1-NM_000138-20230919203303.gbk
</pre>

### fetch_genbank
Fetches the Genbank file from NCBI Gene-Database using the NCBI API. Stores the result
in a .gbk file
<pre>
[Use]       python fetch_genbank.py GENE-ID
[Example]   python fetch_genbank.py FBN1 (output FBN1.gbk)
</pre>

### genbank_to_fasta
Translate a Genbank (.gbk) file to a Fasta (.fas) file 
<pre>
[Use]       /bin/bash ./seq_genbank_to_orfs_fasta.sh GENE-ID
[Example]   /bin/bash ./seq_genbank_to_orfs_fasta.sh FBN1
</pre>

## fasta_to_blast_report
Performs a Blast search with a .fas input file
<pre>
[Use]       python fasta_to_blast_report.py fasta_file.fas path_to_db
[Example]   python fasta_to_blast_report.py gene=FBN1-NM_000138-20230919203303.fas ../../ncbi-blast-2.14.1+/data/swissprot
</pre>

## fasta_popular_alignment
Performs a Multiple Sequence Alignment (MSA)
<pre>
[Use]       python fasta_popular_alignment.py query_sequence_path blast_report.out email path_to_db
[Example]   python fasta_popular_alignment.py /FBN1-NM_000138.fas blast_results.out example@email.com /ncbi-blast/db/swissprot"
</pre>

## get_orfs_and_domains_from_fasta.sh
Calculate the orfs and domains present in a fasta file
<pre>
[Use]       ./get_orfs_and_domains_from_fasta.sh gene-to-orf-and-fasta.fasta
[Example]   ./get_orfs_and_domains_from_fasta.sh huntington-disease.fasta
</pre>

## get_orgs_from_fasta.sh
Calculate the orfs present in a fasta file, rsults will be in orfs.fasta
<pre>
[Use]       ./get_orfs_from_fasta.sh gene-to-orf.fasta
[Example]   ./get_orfs_from_fasta.sh huntington-disease.fasta
</pre>

## get_domains_from_fasta.sh
Calculates domain from a fasta file
<pre>
[Use]       ./get_domains_from_fasta.sh gene-to-domain.fasta
[Example]   ./get_domains_from_fasta.sh huntington-disease.fasta
</pre>

## fasta_to_primer.py
Calculate 5 primers based on a fasta file and parameters in a .json file (with fields: "min_length", "max_length", "min_gc", "max_gc", "max_tm").
Fasta sequence must only have one record.
<pre>
[Use]       ./fasta_to_primer.py gene.fasta config.json
[Example]   ./fasta_to_primer.py FBQ1-NM_0038255.fasta config.json
</pre>
