# Title: 6_sequences_of_duplicates.R
# Author: Nathan Duda
# Date: 10/22/2023
# Purpose: Get sequences and information for each blast hit that hit to a gene. 


source("startup.R")

# read in genes with their blast hits 
dups <- read.csv("./Blast_Outputs/mRNA_raw_Duplicates_from_Blast.tsv", sep="")

# format the file 
dups <- dups[,c(1,2,3,12,8,9)] 
colnames(dups) <- c('gn_id','fly','hit_chrom','hit_strand','hit_start_on_chrom','hit_end_on_chrom')

# read in mRNAs 
##mRNA <- read.csv("./Non_Overlapping_mRNA.tsv", sep="")
##mRNA <- mRNA[,c(7,2,1,9,4,5,6,8)]

# get reverse complement of mRNA sequences on the negative strand 
##for (i in 1:nrow(mRNA)) {
##  if (mRNA$strand[i] == '-') {mRNA$nuc_sequence[i] <- seqRFLP::revComp(mRNA$nuc_sequence[i])}
##  else if (mRNA$strand[i] == '+') {mRNA$nuc_sequence[i] <- mRNA$nuc_sequence[i]}}
##colnames(mRNA)[8] <- 'strand_nuc_sequence'

##write.table(mRNA,'mRNA_Non_Overlapping_Strand_seq.tsv')
mRNA <- read.csv("./mRNA_Non_Overlapping_Strand_seq.tsv", sep="")


# merge the two dataframes to get original mRNA and duplicate information in one dataframe
colnames(mRNA) <- c('gn_id','fly','gn_chrom','gn_FBgn','gn_start','gn_end','gn_strand','gn_strand_nuc_sequence')
dups <- merge(dups, mRNA, by = c('gn_id','fly'))
# the removed rows were the original overlapping mRNAs that the blast was run with but were
# cleaned out of the mRNA dataframe used here

# format the merged dataframe 
dups <- dups[,c(2,1,8,7,11,9,10,12,3:6)]

# remove duplicates that are inside their mrna or hits to their mrna 
dups <- dups %>%
  group_by(gn_chrom,hit_chrom, fly) %>%
  mutate(is_between = case_when((hit_start_on_chrom >= gn_start) & (hit_end_on_chrom <= gn_end) ~ T,
                                (hit_start_on_chrom >= gn_start) & (hit_start_on_chrom <= gn_end) ~ T,
                                (hit_end_on_chrom >= gn_start) & (hit_end_on_chrom <= gn_end) ~ T,
                                T ~ F)) %>%
  filter(is_between == F)


# format mRNA dataframe to get information on gene hit 
colnames(mRNA) <- c('gene_hit_id','fly','hit_chrom','gene_hit_FBgn','gene_hit_start','gene_hit_end','gene_hit_strand','gene_hit_strand_nuc_sequence')
mRNA$hit_length <- mRNA$gene_hit_end - mRNA$gene_hit_start

# get information on gene hit when blast hit to a gene or at most 1/4 of the gene lengths away from the starts/ends of the gene
dups <- dups %>%
  inner_join(mRNA, by = c('fly','hit_chrom')) %>%
  filter((abs(hit_start_on_chrom - gene_hit_start) <= (hit_length * (1/4))) & (abs(hit_end_on_chrom - gene_hit_end) <=  (hit_length * (1/4))))

# write to file 
write.table(dups,'Duplicate_Sequences.tsv')
