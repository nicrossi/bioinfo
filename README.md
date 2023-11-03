# ITBA BioInformática
Trabajo Práctico 2023 2Q 

# Requirements
-Python \
-Blast+ \
-Muscle \
-Blast+ database already set up

-Emboss 

# SCRIPTS
## fetch_genbank
Fetches the Genbank file from NCBI Gene-Database using the NCBI API 
<pre>
[Use]       python fetch_genbank.py GENE-ID
[Example]   python fetch_genbank.py FBN1
</pre>

## genbank_to_fasta
Translate a Genbank (.gbk) file to a Fasta (.fas) file 
<pre>
[Use]       python genbank_to_fasta.py file_to_translate.gbk
[Example]   python genbank_to_fasta.py FBN1-NM_000138-20230919203303.gbk
</pre>

## fasta_to_blast_report
Performs a Blast search with a .fas input file
<pre>
[Use]       python fasta_to_blast_report.py fasta_file.fas path_to_db
[Example]   python fasta_to_blast_report.py gene=FBN1-NM_000138-20230919203303.fas ../../ncbi-blast-2.14.1+/data/swissprot
</pre>

## fasta_popular_alignment
Performs a Multiple Sequence Alignment (MSA)

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
