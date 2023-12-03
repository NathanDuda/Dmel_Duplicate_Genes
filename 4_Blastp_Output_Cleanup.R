# Title: 4_Blastp_Output_Cleanup.R
# Author: Nathan Duda
# Date: 10/29/2023
# Purpose: Format all the blastp outputs to get quality protein duplicates.
#          Combine blastp output with annotation information for each duplicate copy. 


source("startup.R")

# import annotations for all flies 
all_annotations <- read.csv("./Annotations.tsv", sep="")

# keep only the full gene information 
all_annotations <- all_annotations[all_annotations$type=='gene',]
write.table(all_annotations,file='./Annotations_Gene.tsv')
#all_annotations <- read.csv("./Annotations_Gene.tsv", sep="")

# make sure no genes overlap
nonoverlapping_all_annotations <- all_annotations %>%
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
# they dont


# read in list with all fly names 
fly_name_list <- read.table("./Individuals_List.txt", quote="\"", comment.char="")

# create empty dataframe for output
all_dups <- as.data.frame(matrix(ncol=34,nrow=0))

# loop over all flies 
for (row in 1:nrow(fly_name_list)){
  
  # get name of current fly 
  name <- fly_name_list[row,1]
  
  blastp <- read.delim(paste0("./Prot_Blast_Output/",name,"_blastp.tsv"), header=FALSE)
  
  colnames(blastp) <- c('qseqid','sseqid','length','qstart','qend','qlen','sstart','send','slen','pident','evalue')
  
  # remove hits to itself
  blastp <- blastp[blastp$qseqid != blastp$sseqid,]
  
  # filter out hits with percentage identity lower than or equal to 99 and 
  # length lower than or equal to 100 aa
  blastp <- blastp %>% 
    filter((length >= 100) & (pident >= 99))
  
  # keep only the longest hit for same hits 
  blastp <- blastp %>% 
    arrange(length) %>%
    group_by(qseqid, sseqid) %>%
    filter(length == max(length)) %>%
    distinct(qseqid,sseqid, .keep_all = T) %>%
    ungroup() 

  # keep only reciprocal hits 
  blastp <- blastp %>% 
    mutate(q_s = paste0(qseqid,'_',sseqid)) %>%
    mutate(s_q = paste0(sseqid,'_',qseqid)) %>%
    arrange(desc(pident)) %>%
    filter(s_q %in% q_s) 
  
  # keep only one of the reciprocal hits 
  blastp <- blastp %>%
    filter(!(q_s > s_q)) %>%
    select(-q_s,-s_q)
  
  # get annotations for given fly
  annotations <- all_annotations[all_annotations$fly==name,]
  
  # merge blast hits with gene annotations for dup1
  colnames(blastp) <- c('dup1','dup2','match_length','match_dup1_start','match_dup1_end','match_dup1_length','match_dup2_start','match_dup2_end','match_dup2_length','match_pident','match_evalue')
  blastp$dup1 <- paste0(name, '_', blastp$dup1)
  colnames(annotations) <- c('dup1','dup1_type','dup1_start_on_chrom','dup1_end_on_chrom','dup1_n','dup1_strand','dup1_s','dup1_length','dup1_approx_expected_prot_len','dup1_chrom','dup1_prot','dup1_nchar','fly')
  blastp <- merge(blastp,annotations,by='dup1')
  
  # merge blast hits with gene annotations for dup2
  colnames(annotations) <- c('dup2','dup2_type','dup2_start_on_chrom','dup2_end_on_chrom','dup2_n','dup2_strand','dup2_s','dup2_length','dup2_approx_expected_prot_len','dup2_chrom','dup2_prot','dup2_nchar','fly')
  blastp$dup2 <- paste0(name, '_', blastp$dup2)
  blastp <- merge(blastp,annotations,by=c('dup2','fly'))
  
  # combine with all duplicates 
  all_dups <- rbind(all_dups,blastp)
  
  print(row)
}


# write to file 
write.table(all_dups,'./Duplicate_Hits.tsv')


# get duplicate families 

#all_dups <- read.csv("./Duplicate_Hits.tsv", sep="")
#all_dups_orig <- all_dups
#all_dups <- all_dups_orig




fly_name_list <- read.table("./Individuals_List.txt", quote="\"", comment.char="")


dup_families <- as.data.frame(matrix(nrow=0,ncol=3))

# loop over all flies 
for (row in 1:nrow(fly_name_list)){
  
  # get name of current fly
  name <- fly_name_list[row,1]
  
  dup1_dup2 <- all_dups[all_dups$fly==name,c('dup1','dup2')]

  edge_list <- data.frame(source = dup1_dup2$dup1, target = dup1_dup2$dup2)
  
  # Create a graph from the edge list
  graph <- graph.data.frame(edge_list, directed = FALSE)
  
  # Perform community detection using the Louvain algorithm
  community <- cluster_louvain(graph)
  
  # Add the community labels to the nodes
  #V(graph)$community <- membership(community)
  
  # Get the community labels for each node
  node_communities <- data.frame(
    node = community$names,
    community = community$membership
  )
  
  colnames(node_communities) <- c('dup','dup_family')
  node_communities$fly <- name
  
  dup_families <- rbind(dup_families,node_communities)

  print(row)

}




dup1_info <- all_dups[, c('fly', 'dup1', 'dup1_start_on_chrom', 'dup1_end_on_chrom', 'dup1_n', 'dup1_strand', 'dup1_s', 'dup1_length', 'dup1_approx_expected_prot_len', 'dup1_chrom', 'dup1_prot', 'dup1_nchar')]
dup2_info <- all_dups[, c('fly', 'dup2', 'dup2_start_on_chrom', 'dup2_end_on_chrom', 'dup2_n', 'dup2_strand', 'dup2_s', 'dup2_length', 'dup2_approx_expected_prot_len', 'dup2_chrom', 'dup2_prot', 'dup2_nchar')]

colnames(dup1_info) <- c('fly', 'dup', 'dup_start_on_chrom', 'dup_end_on_chrom', 'dup_n', 'dup_strand', 'dup_s', 'dup_length', 'dup_approx_expected_prot_len', 'dup_chrom', 'dup_prot', 'dup_nchar')
colnames(dup2_info) <- c('fly', 'dup', 'dup_start_on_chrom', 'dup_end_on_chrom', 'dup_n', 'dup_strand', 'dup_s', 'dup_length', 'dup_approx_expected_prot_len', 'dup_chrom', 'dup_prot', 'dup_nchar')

dups <- rbind(dup1_info,dup2_info)

# merge with duplicate families
dups <- left_join(dups,dup_families,by=c('fly','dup'))


dups <- dups[!duplicated(dups),]
table(dups$fly)


# remove the crappy assemblies
dups <- dups[!dups$fly %in% c('I23','ZH26','N25','T29A','B59'),]


#
dups$dup_family <- paste0(dups$fly,'_',dups$dup_family)


# write to file 
write.table(dups,'./Duplicate_Proteins.tsv')



