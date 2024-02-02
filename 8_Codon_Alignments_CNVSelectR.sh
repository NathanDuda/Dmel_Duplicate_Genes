




# align proteins of groups of pairs 
for file in ./CNVSelectR/Protein_Sequences/*; do
 filename=$(basename "$file")
 muscle -in "$file" -out "./CNVSelectR/Protein_Alignments/${filename}"
done


# get codon alignment of orthogroups with pal2nal
for file in ./CNVSelectR/Protein_Sequences/*; do
 filename=$(basename "$file")
 pal2nal.pl "./CNVSelectR/Protein_Alignments/${filename}" "./CNVSelectR/Nucleotide_Sequences/${filename}" -output fasta > "./CNVSelectR/Codon_Alignments/${filename}"
done

# combine 


for file in ./Protein_Alignments/*.*; do
    base_name=$(basename "$file" | cut -d '.' -f 1)
    cat "$file" >> "./Combined_Protein_Alignments/${base_name}.fa"
done







