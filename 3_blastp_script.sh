#!/bin/bash


names=(
  "A1"
  "A2"
  "A3"
  "A4"
  "A5"
  "A6"
  "A7"
  "AB8"
  "AKA-017"
  "AKA-018"
  "B1"
  "B2"
  "B3"
  "B4"
  "B59"
  "B6"
  "COR-014"
  "COR-018"
  "COR-023"
  "COR-025"
  "GIM-012"
  "GIM-024"
  "I23"
  "ISO-1"
  "JUT-008"
  "JUT-011"
  "KIE-094"
  "LUN-004"
  "LUN-007"
  "MUN-008"
  "MUN-009"
  "MUN-013"
  "MUN-015"
  "MUN-016"
  "MUN-020"
  "N25"
  "ORE"
  "RAL-059"
  "RAL-091"
  "RAL-176"
  "RAL-177"
  "RAL-375"
  "RAL-426"
  "RAL-737"
  "RAL-855"
  "SLA-001"
  "STO-022"
  "T29A"
  "TEN-015"
  "TOM-007"
  "TOM-008"
  "ZH26"
)


for name in "${names[@]}"; do

  makeblastdb -in "./Prot_Fastas/${name}_prot.fasta" -out "./Prot_Blast_Output/Blast_Databases/${name}_db" -dbtype prot 

  blastp -query "./Prot_Fastas/${name}_prot.fasta" -db "./Prot_Blast_Output/Blast_Databases/${name}_db" -outfmt "6 qseqid sseqid length qstart qend qlen sstart send slen pident evalue" -out "./Prot_Blast_Output/${name}_blastp.tsv" &

done

wait





