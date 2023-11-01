# Title: 7_FBgn_attachment.R
# Author: Nathan Duda
# Date: 10/31/2023
# Purpose: Get FBgns for all my annotated genes.
#          Make a list of all FBgns in my annotations, all FBgns in my duplicates,
#          and all FBgns that are not duplicates.

source("startup.R")


# import flybase proteins 
fbgn_prots <- read.csv("./dmel-all-translation-r6.54.fasta", sep="", header = F)

# format the fasta into a dataframe 
colnames(fbgn_prots)[1] <- 'prot_seq'
fbgn_prots <- fbgn_prots %>%
  mutate(Group = cumsum(grepl("^>", prot_seq))) %>%
  group_by(Group) %>%
  mutate(prot = paste(prot_seq, collapse = "")) %>%
  ungroup() %>%
  mutate(prot_seq = prot) %>%
  select(-Group,-prot)

fbgn_prots$FBpp <- substr(fbgn_prots$prot_seq, start = 2, stop = 12)
fbgn_prots$prot_seq <- substring(fbgn_prots$prot_seq, first = 13)


# remove the unnecessary rows 
fbgn_prots <- fbgn_prots %>% mutate_all(~na_if(., "")) %>% na.omit()

# extract the FBgn id 
fbgn_prots$FBgn <- substr(fbgn_prots$V6, start = 8, stop = 18)

# format the dataframe
fbgn_prots <- fbgn_prots[,c(12,13,1,3,5:9)]

# get only the FBgn and FBpp 
fbgn_fbpp <- fbgn_prots[c('FBgn','FBpp')]
rm(fbgn_prots)


### get FBgn for each gn from FBgn blast 

# read in list with all fly names 
fly_name_list <- read.table("./Individuals_List.txt", quote="\"", comment.char="")

# create empty dataframes for output
my_FBgns <- as.data.frame(matrix(ncol=4,nrow=0))

# loop over all flies 
for (row in 1:nrow(fly_name_list)){
  
  # get name of current fly 
  name <- fly_name_list[row,1]
  
  # import blast result
  blastp_FBgn <- read.delim(paste0("C:/Users/17735/Downloads/Dmel_Duplicate_Genes/FBgn_Blast_Output/",name,"_blastp.tsv"), header=FALSE)
  colnames(blastp_FBgn) <- c('qseqid', 'sseqid', 'length', 'bitscore', 'pident', 'evalue')
  
  # keep only the gene with its FBgn and the evalue
  blastp_FBgn <- blastp_FBgn[,c('qseqid','sseqid','evalue')]
  colnames(blastp_FBgn) <- c('gn','FBpp','evalue')
  
  # keep only highest evalue blast
  #blastp_FBgn <- blastp_FBgn %>%
  #  group_by(gn) %>%
  #  filter(evalue == min(evalue)) %>%
  #  ungroup() %>%
  #  distinct(gn,evalue, .keep_all = T) # keep only one if lowest evalue is same for multiple 

  # combine all genes with their FBgns into dataframe 
  blastp_FBgn$fly <- name
  my_FBgns <- rbind(my_FBgns,blastp_FBgn)
  
  print(row)
  
}

# get FBgn for each pp
my_FBgns <- merge(my_FBgns,fbgn_fbpp,by='FBpp')

length(unique(my_FBgns$FBgn))
length(unique(my_FBgns$gn))

length(unique(fbgn_fbpp$FBgn))


# write to table
write.table(my_FBgns,'./gn_FBgn_FBpp.tsv')

# extract the FBgns 
FBgn <- my_FBgns[c('FBgn')]
FBgn <- unique(FBgn)

# write a list of FBgns 
write.table(FBgn, file = './FBgn_list.txt', sep = "\t", row.names = F, col.names = F, quote = F)



#############
# get FBgns for my duplicates 

# import my duplicate genes 
dups <- read.csv("./Duplicate_Proteins.tsv", sep="")

#dup_FBgns_orig <- my_FBgns
dup_FBgns <- dup_FBgns_orig


# merge FBgns with dups dataframe

# remove duplicates caused by multiple FBpp per FBgn
dup_FBgns <- dup_FBgns[!duplicated(dup_FBgns[c('gn','fly','FBgn')]),]
dup_FBgns <- dup_FBgns[c('gn','fly','FBgn')]

colnames(dup_FBgns) <- c('dup1','fly','dup1_FBgn')
dups <- left_join(dups,dup_FBgns,by=c('dup1','fly'))

colnames(dup_FBgns) <- c('dup2','fly','dup2_FBgn')
dups <- left_join(dups,dup_FBgns,by=c('dup2','fly'))


# get a list of FBgns of just the duplicates

dup1_fbgn <- dups[c('dup1_FBgn')]
dup2_fbgn <- dups[c('dup2_FBgn')]

colnames(dup1_fbgn) <- 'FBgn'
colnames(dup2_fbgn) <- 'FBgn'

dup_fbgn <- rbind(dup1_fbgn,dup2_fbgn)

dup_fbgn <- dup_fbgn[!duplicated(dup_fbgn),]
dup_fbgn <- as.data.frame(dup_fbgn)

# write a list of FBgns of duplicates 
write.table(dup_fbgn, file = './Dup_FBgn_list.txt', sep = "\t", row.names = F, col.names = F, quote = F)


# get a list of FBgns of not duplicates 
nondup_fbgn <- dup_fbgn[!(FBgn$FBgn %in% dup_fbgn$dup_fbgn),]
nondup_fbgn <- na.omit(as.data.frame(nondup_fbgn))

# write a list of FBgns of that are not duplicates 
write.table(nondup_fbgn, file = './ NonDup_FBgn_list.txt', sep = "\t", row.names = F, col.names = F, quote = F)

