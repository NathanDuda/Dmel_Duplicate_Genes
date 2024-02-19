

# Author: Nathan Duda
# Purpose: Get codon alignments of duplicate pairs from prot and nuc fastas for CNVSelectR 

# connected eq nucleotides:
# align proteins of groups of pairs 
for file in ./CNVSelectR/Connected_Eq_Protein_Sequences/*; do
 filename=$(basename "$file")
 muscle -in "$file" -out "./CNVSelectR/Connected_Eq_Protein_Alignments/${filename}"
done

# get codon alignment of orthogroups with pal2nal
for file in ./CNVSelectR/Connected_Eq_Protein_Sequences/*; do
 filename=$(basename "$file")
 pal2nal.pl "./CNVSelectR/Connected_Eq_Protein_Alignments/${filename}" "./CNVSelectR/Connected_Eq_Nucleotide_Sequences/${filename}" -output fasta > "./CNVSelectR/Connected_Eq_Codon_Alignments/${filename}"
done

# combine pairs of codon alignments with codon alignments of groupmate pairs 
for file in ./CNVSelectR/Connected_Eq_Codon_Alignments/*; do
    base_name=$(basename "$file" | sed 's/\(group_[0-9]\+\).*\(\..*\)/\1\2/')
    cat "$file" >> "./CNVSelectR/Connected_Eq_Combined_Codon_Alignments/${base_name}"
done

