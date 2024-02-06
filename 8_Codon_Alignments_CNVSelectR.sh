

# global alignments:
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

# combine pairs of codon alignments with codon alignments of groupmate pairs 
for file in ./CNVSelectR/Codon_Alignments/*; do
    base_name=$(basename "$file" | sed 's/\(group_[0-9]\+\).*\(\..*\)/\1\2/')
    cat "$file" >> "./CNVSelectR/Combined_Codon_Alignments/${base_name}"
done


######################################

# blastp nucleotides:
# align proteins of groups of pairs 
for file in ./CNVSelectR/Blastp_Dups_Protein_Sequences/*; do
 filename=$(basename "$file")
 muscle -in "$file" -out "./CNVSelectR/Blastp_Dups_Protein_Alignments/${filename}"
done

# get codon alignment of orthogroups with pal2nal
for file in ./CNVSelectR/Blastp_Dups_Protein_Sequences/*; do
 filename=$(basename "$file")
 pal2nal.pl "./CNVSelectR/Blastp_Dups_Protein_Alignments/${filename}" "./CNVSelectR/Blastp_Dups_Nucleotide_Sequences/${filename}" -output fasta > "./CNVSelectR/Blastp_Dups_Codon_Alignments/${filename}"
done

# combine pairs of codon alignments with codon alignments of groupmate pairs 
for file in ./CNVSelectR/Blastp_Dups_Codon_Alignments/*; do
    base_name=$(basename "$file" | sed 's/\(group_[0-9]\+\).*\(\..*\)/\1\2/')
    cat "$file" >> "./CNVSelectR/Blastp_Dups_Combined_Codon_Alignments/${base_name}"
done


