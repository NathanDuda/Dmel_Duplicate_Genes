

source("startup.R")

# import annotations for all flies 
all_annotations <- read.csv("./Dsim/Dsim_Annotations.tsv", sep="")

# keep only the full gene information 
all_annotations <- all_annotations[all_annotations$type=='gene',]
write.table(all_annotations,file='./Dsim/Dsim_Annotations_Gene.tsv')

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



blastp <- read.delim("./Dsim/Prot_Blast_Output/Dsim_blastp.tsv", header=FALSE)

colnames(blastp) <- c('qseqid','sseqid','length','qstart','qend','qlen','sstart','send','slen','pident','evalue')

# remove hits to itself
blastp <- blastp[blastp$qseqid != blastp$sseqid,]

# filter out hits with percentage identity lower than or equal to 99 and 
# length lower than or equal to 100 aa
blastp <- blastp %>% 
  filter((length >= 100) & (pident >= 99))

# keep only the longest hit for same hits 
blastp <- blastp %>% 
  group_by(qseqid, sseqid) %>%
  filter(length == max(length)) %>%
  ungroup() 

# keep only one of the reciprocal hits 
blastp <- blastp %>% 
  mutate(q_s = paste0(qseqid,'_',sseqid)) %>%
  mutate(s_q = paste0(sseqid,'_',qseqid)) %>%
  arrange(desc(pident)) %>%
  filter(!((s_q %in% q_s) & (q_s > s_q))) %>%
  select(-q_s,-s_q)


# merge blast hits with gene annotations for dup1
colnames(blastp) <- c('dup1','dup2','match_length','match_dup1_start','match_dup1_end','match_dup1_length','match_dup2_start','match_dup2_end','match_dup2_length','match_pident','match_evalue')
all_annotations <- all_annotations[,-c(11)] # remove extra gene id. column name "aa"
colnames(all_annotations) <- c('dup1','dup1_type','dup1_start_on_chrom','dup1_end_on_chrom','dup1_n','dup1_strand','dup1_s','dup1_length','dup1_approx_expected_prot_len','dup1_chrom','dup1_prot','dup1_nchar','fly')
blastp <- merge(blastp,all_annotations,by='dup1')

# merge blast hits with gene annotations for dup2
colnames(all_annotations) <- c('dup2','dup2_type','dup2_start_on_chrom','dup2_end_on_chrom','dup2_n','dup2_strand','dup2_s','dup2_length','dup2_approx_expected_prot_len','dup2_chrom','dup2_prot','dup2_nchar','fly')
blastp <- merge(blastp,all_annotations,by=c('dup2','fly'))


# write to file 
write.table(blastp,'./Dsim/Duplicate_Hits.tsv')



# get duplicate families 

all_dups_orig <- blastp
all_dups <- all_dups_orig




edge_list <- data.frame(source = all_dups$dup1, target = all_dups$dup2)

# Create a graph from the edge list
graph <- graph.data.frame(edge_list, directed = FALSE)

# Perform community detection using the Louvain algorithm
community <- cluster_louvain(graph)

# Add the community labels to the nodes
V(graph)$community <- membership(community)

# Get the community labels for each node
node_communities <- data.frame(
  node = community$names,
  community = community$membership)


colnames(node_communities) <- c('dup','dup_family')
node_communities$fly <- 'Dsim'

dup_families <- node_communities




dup1_info <- all_dups[, c('fly', 'dup1', 'dup1_start_on_chrom', 'dup1_end_on_chrom', 'dup1_n', 'dup1_strand', 'dup1_s', 'dup1_length', 'dup1_approx_expected_prot_len', 'dup1_chrom', 'dup1_prot', 'dup1_nchar')]
dup2_info <- all_dups[, c('fly', 'dup2', 'dup2_start_on_chrom', 'dup2_end_on_chrom', 'dup2_n', 'dup2_strand', 'dup2_s', 'dup2_length', 'dup2_approx_expected_prot_len', 'dup2_chrom', 'dup2_prot', 'dup2_nchar')]

colnames(dup1_info) <- c('fly', 'dup', 'dup_start_on_chrom', 'dup_end_on_chrom', 'dup_n', 'dup_strand', 'dup_s', 'dup_length', 'dup_approx_expected_prot_len', 'dup_chrom', 'dup_prot', 'dup_nchar')
colnames(dup2_info) <- c('fly', 'dup', 'dup_start_on_chrom', 'dup_end_on_chrom', 'dup_n', 'dup_strand', 'dup_s', 'dup_length', 'dup_approx_expected_prot_len', 'dup_chrom', 'dup_prot', 'dup_nchar')

dups <- rbind(dup1_info,dup2_info)

# merge with duplicate families
dups <- left_join(dups,dup_families,by=c('fly','dup'))


# write to file 
write.table(dups,'./Dsim/Duplicate_Proteins.tsv')




