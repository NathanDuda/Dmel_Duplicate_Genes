#!/bin/bash


# list of names
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


# running augustus 
for name in "${names[@]}"; do
  ./Augustus/src/augustus --species=fly "./Input_Genomes/${name}.fasta" > "./Output_Annotations/${name}_annotations.gff" &
done
wait


# formatting the output files 
for name in "${names[@]}"; do
cp "./Augustus_Output/${name}_annotations.aa" "./Prot_Fastas/${name}_annotations.fasta"
done

cd ./Prot_Fastas
for name in "${names[@]}"; do
sed -i 's/\.t1//g' "./${name}_annotations.fasta"
done


# getting the number of genes per chromosome: 
for name in "${names[@]}"; do
  awk 'BEGIN { count = 0; inside_sequence = 0; }
       /sequence number/ { if (inside_sequence) print count; inside_sequence = 1; count = 0; next; }
       /protein sequence/ { if (inside_sequence) count++; }
       END { if (inside_sequence) print count; }' "${name}_annotations.gff" > "./aa_per_chrom/${name}_aa_per_chrom.txt"
done


