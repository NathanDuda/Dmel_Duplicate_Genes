
#!/bin/bash


# run blast
./ncbi-blast-2.14.1+/bin/makeblastdb -in ./Dsim_annotations.fasta -out ./Prot_Blast_Output/Blast_Databases/Dsim_db -dbtype prot 

./ncbi-blast-2.14.1+/bin/blastp -query ./Dsim_annotations.fasta -db ./Prot_Blast_Output/Blast_Databases/Dsim_db -outfmt "6 qseqid sseqid length qstart qend qlen sstart send slen pident evalue" -out ./Prot_Blast_Output/Dsim_blastp.tsv

