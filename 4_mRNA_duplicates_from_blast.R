# Title: 3_mRNA_duplicates_from_blast
# Author: Nathan Duda
# Date: 10/6/2023
# Purpose: Format the file with all mRNA sequences to not contain overlapping mRNAs.
#          Use this to clean the overlapping mRNA query_ids from blast outputs.
#          Extract quality blast hits from the mRNA blast output files.
#          Remove overlapping blast hits and write all to a file.

source("startup.R")


# read in all mRNA sequences
all_mRNA <- read.csv("./all_mRNA.tsv", sep="")

# find overlapping mRNAs and mark them with TRUE in is_between and is_between_end columns
nonoverlapping_mrna <- all_mRNA %>%
  group_by(chrom, fly) %>%
  arrange(start) %>%
  mutate(is_between = case_when(
    start >= lag(start) & start <= lag(end) ~ T,
    start >= lead(start) & start <= lead(end) ~ T,
    T ~ F)) %>%
  arrange(end) %>%
  mutate(is_between_end = case_when(
    end >= lag(start) & end <= lag(end) ~ T,
    end >= lead(start) & end <= lead(end) ~ T,
    T ~ F)) %>%
  ungroup()

# extract FBgn into new column 
nonoverlapping_mrna <- nonoverlapping_mrna %>% mutate(FBgn = sub(".*\\b(FBgn\\d{7})\\b.*", "\\1", id))

# get lengths of sequences 
nonoverlapping_mrna$nchar <- nchar(nonoverlapping_mrna$nuc_sequence)

# remove overlapping mrnas and keeping the longest ones 
nonoverlapping_mrna <- nonoverlapping_mrna %>%
  group_by(chrom, fly, FBgn) %>%
  mutate(keep = case_when(is_between == F & is_between_end == F ~ 'yes',
                          is_between == T & is_between_end == T & nchar == max(nchar) ~ 'yes',
                          is_between == F & is_between_end == T & nchar == max(nchar) ~ 'yes',
                          is_between == T & is_between_end == F & nchar == max(nchar) ~ 'yes')) %>%
  ungroup() %>%
  na.omit() 

# keep only one random mrna when multiple are in the same exact spot
nonoverlapping_mrna <- nonoverlapping_mrna %>%
  distinct(chrom,fly,FBgn, .keep_all = T) %>%
  select(-c(is_between,is_between_end, keep)) # remove unnecessary columns

# remove the two mrnas annotated twice on chromosomes different from all other flies: 
counts <- as.data.frame(table(nonoverlapping_mrna$fly, nonoverlapping_mrna$FBgn))
# RAL-091 FBgn0085334 on X chromosome 
# MUN-009 FBgn0000500 on 2R chromosome
nonoverlapping_mrna <- nonoverlapping_mrna %>%
  filter(!(fly == 'RAL-091' & FBgn == 'FBgn0085334' & chrom == 'X')) %>%
  filter(!(fly == 'MUN-009' & FBgn == 'FBgn0000500' & chrom == '2R'))

# write to table 
write.table(nonoverlapping_mrna,'Non_Overlapping_mRNA.tsv')

# these nonoverlapping mRNAs will be used to filter out the overlapping mRNAs used in blast



################ get mRNA duplicates from blast outputs:

all_mRNA_orig <- nonoverlapping_mrna

# read in list of all fly names 
list <- read.table("./Individuals_List.txt", quote="\"", comment.char="")

# create empty dataframe for output
duplicates <- as.data.frame(matrix(ncol = 11, nrow = 0))

for (row_num in 1:nrow(list)){
  
  # get species name 
  name <- list[row_num,1]
  
  # extract exon for given species 
  all_mRNA <- all_mRNA_orig[all_mRNA_orig$fly == paste(name),]
  
  # get blast output file for that species 
  blast_output <- read.delim(paste("./Blast_Outputs/mRNA_blast_output/",name,"_results.txt",sep=''), header=FALSE)
  
  # format the blast output file
  colnames(blast_output) <- c('query_id', 'subject_id', 'perc_identity', 'alig_length', 'q_start', 'q_end', 's_start', 's_end', 'evalue', 'bitscore')
  blast_output[c('individual', 'subject_id')] <- str_split_fixed(blast_output$subject_id, '_', 2)
  blast_output$individual <- gsub('.fasta','',blast_output$individual)
  blast_output$query_id <- gsub(".*_ID=","ID=",blast_output$query_id)
  blast_output <- blast_output[,c('query_id', 'individual', 'subject_id', 'perc_identity', 'alig_length','q_start', 'q_end', 's_start', 's_end', 'evalue', 'bitscore')]
  
  # keep blast hit only if the query_id (mrna) is in the non-overlapping mrna dataframe
  blast_output <- blast_output %>%
    filter(query_id %in% nonoverlapping_mrna$id)
  
  # keep only hits with percentage identity above 90 and alignment length above 300 nucleotides
  blast_output <- blast_output[blast_output$perc_identity > 90,]
  blast_output <- blast_output[blast_output$alig_length > 300,]
  
  # remove overlapping blast hits 
  blast_output <- blast_output %>%
    
    # if start coordinate is larger than end coordinate, flip:
    mutate(s_end_r = case_when(s_end > s_start ~ s_end,
                               s_start > s_end ~ s_start)) %>%
    mutate(s_start = case_when(s_start > s_end ~ s_end,
                               s_end > s_start ~ s_start,)) %>% 
    mutate(s_end = s_end_r) %>% 
    select(-s_end_r) %>% # remove intermediate column
    
    # remove overlapping blast hits 
    group_by(subject_id, individual) %>% # subject_id is chromosome
    arrange(s_start) %>%
    mutate(is_between = case_when(
      s_start >= lag(s_start) & s_start <= lag(s_end) ~ T,
      s_start >= lead(s_start) & s_start <= lead(s_end) ~ T,
      T ~ F)) %>%
    arrange(s_end) %>%
    mutate(is_between_end = case_when(
      s_end >= lag(s_start) & s_end <= lag(s_end) ~ T,
      s_end >= lead(s_start) & s_end <= lead(s_end) ~ T,
      T ~ F)) %>%
    filter(is_between == F & is_between_end == F) %>%
    select(-c(is_between,is_between_end)) %>%
    ungroup()
  
  # append blast output to all duplicates
  duplicates <- rbind(duplicates,blast_output)  
  
  # print progress
  print(row_num)
  
}

# write duplicates to table 
colnames(duplicates) <- colnames(blast_output)
write.table(duplicates, './Blast_Outputs/Duplicates_from_Blast.tsv')

  