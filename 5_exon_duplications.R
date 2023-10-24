

source("startup.R")

# read in all exon sequences
all_exon <- read.csv("./all_exon.tsv", sep="")

# find overlapping exons and mark them with TRUE in is_between and is_between_end columns
nonoverlapping_exon <- all_exon %>%
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

# get lengths of sequences 
##nonoverlapping_exon$nchar <- nchar(nonoverlapping_exon$nuc_sequence)

# remove overlapping exons
##nonoverlapping_exon <- nonoverlapping_exon %>%
##  group_by(chrom, fly, id) %>%
##  mutate(keep = case_when(is_between == F & is_between_end == F ~ 'yes',
##                          is_between == T & is_between_end == T ~ NA,
##                          is_between == F & is_between_end == T ~ NA,
##                          is_between == T & is_between_end == F ~ NA)) %>%
##  ungroup() %>%
##  na.omit() 

# make sure no exons are in the same exact spot
##nonoverlapping_exon <- nonoverlapping_exon %>%
##  distinct(chrom,fly,id,start,end, .keep_all = T) %>%
##  select(-c(is_between,is_between_end, keep)) # remove unnecessary columns

# write to table 
##write.table(nonoverlapping_exon,'Non_Overlapping_exon.tsv')
nonoverlapping_exon <- read.csv("./Non_Overlapping_exon.tsv", sep="")
# these nonoverlapping exons will be used to filter out the overlapping exons used in blast


################ get exon duplicates from blast outputs:

all_exon_orig <- nonoverlapping_exon

# read in list of all fly names 
list <- read.table("./Individuals_List.txt", quote="\"", comment.char="")

# create empty dataframe for output
duplicates <- as.data.frame(matrix(ncol = 12, nrow = 0))

for (row_num in 1:nrow(list)){
  
  # get species name 
  name <- list[row_num,1]
  
  # extract exon for given species 
  all_exon <- all_exon_orig[all_exon_orig$fly == paste(name),]
  
  # get blast output file for that species 
  blast_output <- read.delim(paste("./Blast_Outputs/exon_blast_output/",name,"_results.txt",sep=''), header=FALSE)
  
  # format the blast output file
  colnames(blast_output) <- c('query_id', 'subject_id', 'perc_identity', 'alig_length', 'q_start', 'q_end', 's_start', 's_end', 'evalue', 'bitscore')
  blast_output[c('individual', 'subject_id')] <- str_split_fixed(blast_output$subject_id, '_', 2)
  blast_output$individual <- gsub('.fasta','',blast_output$individual)
  blast_output$query_id <- gsub(".*_Parent=","Parent=",blast_output$query_id)
  blast_output <- blast_output[,c('query_id', 'individual', 'subject_id','perc_identity', 'alig_length','q_start', 'q_end', 's_start', 's_end', 'evalue', 'bitscore')]
  
  # keep blast hit only if the query_id (exon) is in the non-overlapping exon dataframe
  blast_output <- blast_output %>%
    filter(query_id %in% nonoverlapping_exon$id)
  
  # keep only hits with percentage identity above 90 and alignment length above 300 nucleotides
  blast_output <- blast_output[blast_output$perc_identity > 90,]
  blast_output <- blast_output[blast_output$alig_length > 30,]
  
  # remove overlapping blast hits 
  blast_output <- blast_output %>%
    
    # get strand hit is on  
    mutate(strand = case_when(s_end < s_start ~ '-',
                              s_end > s_start ~ '+')) %>%
    
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
write.table(duplicates, './Blast_Outputs/exon_raw_Duplicates_from_Blast.tsv')


















#########################

# read in genes with their blast hits 
dups <- read.csv("./Blast_Outputs/exon_raw_Duplicates_from_Blast.tsv", sep="")

# format the file 
dups <- dups[,c(1,2,3,12,8,9)] 
colnames(dups) <- c('gn_id','fly','hit_chrom','hit_strand','hit_start_on_chrom','hit_end_on_chrom')

# read in exons 
##exon <- read.csv("./Non_Overlapping_exon.tsv", sep="")
##exon <- exon[,c(7,2,1,9,4,5,6,8)]

# get reverse complement of exon sequences on the negative strand 
##for (i in 1:nrow(exon)) {
##  if (exon$strand[i] == '-') {exon$nuc_sequence[i] <- seqRFLP::revComp(exon$nuc_sequence[i])}
##  else if (exon$strand[i] == '+') {exon$nuc_sequence[i] <- exon$nuc_sequence[i]}}
##colnames(exon)[8] <- 'strand_nuc_sequence'

##write.table(exon,'exon_Non_Overlapping_Strand_seq.tsv')
exon <- read.csv("./exon_Non_Overlapping_Strand_seq.tsv", sep="")


# merge the two dataframes to get original exon and duplicate information in one dataframe
colnames(exon) <- c('gn_id','fly','gn_chrom','gn_FBgn','gn_start','gn_end','gn_strand','gn_strand_nuc_sequence')
dups <- merge(dups, exon, by = c('gn_id','fly'))
# the removed rows were the original overlapping exons that the blast was run with but were
# cleaned out of the exon dataframe used here

# format the merged dataframe 
dups <- dups[,c(2,1,8,7,11,9,10,12,3:6)]

# remove duplicates that are inside their exon or hits to their exon 
dups <- dups %>%
  group_by(gn_chrom,hit_chrom, fly) %>%
  mutate(is_between = case_when((hit_start_on_chrom >= gn_start) & (hit_end_on_chrom <= gn_end) ~ T,
                                (hit_start_on_chrom >= gn_start) & (hit_start_on_chrom <= gn_end) ~ T,
                                (hit_end_on_chrom >= gn_start) & (hit_end_on_chrom <= gn_end) ~ T,
                                T ~ F)) %>%
  filter(is_between == F)


# format exon dataframe to get information on gene hit 
colnames(exon) <- c('gene_hit_id','fly','hit_chrom','gene_hit_FBgn','gene_hit_start','gene_hit_end','gene_hit_strand','gene_hit_strand_nuc_sequence')
exon$hit_length <- exon$gene_hit_end - exon$gene_hit_start

# get information on gene hit when blast hit to a gene or at most 1/4 of the gene lengths away from the starts/ends of the gene



#

# separate the dataframes by chromosome and fly because otherwise the dataframes are too big 

##trash <- dups %>%
##  inner_join(exon, by = c('fly','hit_chrom')) %>%
##  filter((abs(hit_start_on_chrom - gene_hit_start) <= (hit_length * (1/4))) & (abs(hit_end_on_chrom - gene_hit_end) <=  (hit_length * (1/4))))


#

dups_orig <- dups
exon_orig <- exon

dups <- dups_orig
exon <- exon_orig

output <- as.data.frame(matrix(ncol = 12, nrow = 0))

for (flyy in (unique(dups_orig$fly))){
  for (chrom in (unique(dups_orig$hit_chrom))){
    dups <- dups_orig[(dups_orig$hit_chrom == chrom) & (dups_orig$fly == flyy),]
    exon <- exon_orig[exon_orig$hit_chrom == chrom & exon_orig$fly == flyy,]
    
    trash <- dups %>%
      inner_join(exon, by = c('hit_chrom','fly')) %>%
      filter((abs(hit_start_on_chrom - gene_hit_start) <= (hit_length * (1/4))) & (abs(hit_end_on_chrom - gene_hit_end) <=  (hit_length * (1/4))))
    output <- rbind(output,trash)  
    
  }
  print(flyy)
}

colnames(output) <- colnames(trash)
#

# remove overlapping exon hits 

output <- output %>%
  group_by(hit_chrom, fly) %>%
  arrange(gene_hit_start) %>%
  mutate(is_between = case_when(
    gene_hit_start >= lag(gene_hit_start) & gene_hit_start <= lag(gene_hit_end) ~ T,
    gene_hit_start >= lead(gene_hit_start) & gene_hit_start <= lead(gene_hit_end) ~ T,
    T ~ F)) %>%
  arrange(gene_hit_end) %>%
  mutate(is_between_end = case_when(
    gene_hit_end >= lag(gene_hit_start) & gene_hit_end <= lag(gene_hit_end) ~ T,
    gene_hit_end >= lead(gene_hit_start) & gene_hit_end <= lead(gene_hit_end) ~ T,
    T ~ F)) %>%
  ungroup()

output <- output %>%
  group_by(hit_chrom, fly) %>%
  mutate(keep = case_when(is_between == F & is_between_end == F ~ 'yes')) %>%
  ungroup() %>%
  na.omit() 

# remove hits to itself
output <- output[output$gn_id != output$gene_hit_id,]

# write to file 
write.table(output,'Exon_Duplicate_Sequences.tsv')



