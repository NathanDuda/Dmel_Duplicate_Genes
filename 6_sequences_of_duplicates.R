# Title: 6_sequences_of_duplicates.R
# Author: Nathan Duda
# Date: 10/8/2023
# Purpose: Extract sequences for each blast hit from the according chromosomes. 


source("startup.R")

# read in genes with their blast hits 
dups <- read.csv("./Blast_Outputs/mRNA_raw_Duplicates_from_Blast.tsv", sep="")

# format the file 
dups <- dups[,c(1,2,3,12,8,9)] 
colnames(dups) <- c('gn_id','fly','hit_chrom','hit_strand','hit_start_on_chrom','hit_end_on_chrom')

# read in mRNAs 
mRNA <- read.csv("./Non_Overlapping_mRNA.tsv", sep="")
mRNA <- mRNA[,c(7,2,1,9,4,5,6,8)]

# get reverse complement of mRNA sequences on the negative strand 
mRNA <- mRNA %>% 
  mutate(nuc_sequence = case_when(strand == '-' ~ seqRFLP::revComp(nuc_sequence),
                                  strand == '+' ~ nuc_sequence))
colnames(mRNA)[8] <- 'strand_nuc_sequence'

#write.table(mRNA,'mRNA_Non_Overlapping_Strand_seq.tsv')
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


# import chromosome sequences
genomes <- readDNAStringSet("./Genome_Assemblies/All_Genome_Assemblies.fasta")

# get sequences of blast hits
dups <- dups %>%
  mutate(nuc_sequence = as.character(subseq(genomes[paste0(' ',fly,'.fasta_',hit_chrom, sep = '')],
                                            start = hit_start_on_chrom, end = hit_end_on_chrom))) %>%
  select(-is_between)

# get reverse complement of duplicate sequences hit to negative strand 
trash <- dups %>% 
  mutate(nuc_sequence = case_when(hit_strand == '-' ~ seqRFLP::revComp(nuc_sequence),
                                  hit_strand == '+' ~ nuc_sequence))
colnames(dups)[13] <- 'hit_strand_nuc_sequence'

# write to file
write.table(dups,'Duplicate_Sequences.tsv')
# both nucleotide columns are the reverse complement where necessary, no further changes are needed


