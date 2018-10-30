# Annovar2VCF
Python code for converting annovar format into VCF format.
Annovar result file is different from generic VCF format, especially on indel notation.
This code easily convert it by deploying perl code provided by annovar developers. [retrieve_seq_from_fasta.pl]

Requirements
- "retrieve_seq_from_fasta.pl" from annovar at PATH environment variable.
- reference fasta file for retrieving reference sequences of indels

Input format
- TAB-delimited table_annovar result file

output format
- VCF format
