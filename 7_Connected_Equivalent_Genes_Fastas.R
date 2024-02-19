

# Author: Nathan Duda
# Purpose: Connect equivalent genes that appear multiple times in the dataframe to each other
#          to get complete equivalent genes.
#          Get prot and nuc fastas of duplicate pairs for CNVSelectR


#########################################################################################################
# on cluster:

library(tidyr)
library(dplyr)
library(igraph)


all_genes <- read.csv("./Equivalent_Genes.tsv", sep="")

dups <- read.csv("./Duplicate_Proteins.tsv", sep="")


# format equivalent genes into pairs 
eq_pairs <- all_genes %>%
  pivot_longer(., cols=c(2:ncol(all_genes))) %>%
  filter(value != 0)


# add dup pairs to equivalent gene pairs 
dups <- dups %>% select(fly, dup, dup_family)

create_pairs <- function(data) {
  pairs <- expand.grid(data$dup, data$dup)
  pairs <- pairs[pairs[,1] != pairs[,2] & as.character(pairs[,1]) < as.character(pairs[,2]),]
  pairs
}

result <- dups %>%
  group_by(dup_family) %>%
  do(create_pairs(.))

dups <- result %>%
  select(gene_group = Var1,
         name = dup_family,
         value = Var2) %>%
  mutate(name = sub("_.*", "", name))

pairs <- rbind(eq_pairs, dups)


# separate multi-copy duplicates into pairs 
pairs <- pairs %>%
  select(gene_group, value) %>%
  separate_rows(value, sep = ",")


# build communities from pair connections 
edge_list <- data.frame(source = pairs$gene_group, target = pairs$value)
graph <- graph.data.frame(edge_list, directed = FALSE)

community <- cluster_louvain(graph) # TRY DIFF CLUSTERING

node_communities <- data.frame(
  node = community$names,
  community = community$membership
)

write.table(node_communities, file = 'connected_dups.tsv')


#########################################################################################################


source("startup.R")


connected_dups <- read.csv("./connected_dups.tsv", sep="")

connected_dups <- connected_dups %>%
  separate(node, sep = '_', into = c('fly','gn'))


community_fly_counts <- as.data.frame(table(connected_dups$fly, connected_dups$community))

anys_to_anys <- community_fly_counts %>%
  group_by(Var2) %>%
  filter(all(Freq %in% c(0, 1, 2)))


one_to_ones <- community_fly_counts %>%
  group_by(Var2) %>%
  filter(all(Freq == 1))
length(unique(one_to_ones$Var2)) # 816 one to ones out of the 25,005 total communities 


one_to_twos <- anys_to_anys %>%
  group_by(Var2) %>%
  filter(any(Freq == 2))
length(unique(one_to_twos$Var2))


one_to_two_communities <- one_to_twos %>% select(Var2) %>% distinct()

one_to_two_genes <- connected_dups %>% 
  filter(community %in% one_to_two_communities$Var2)



########################################################################################

# write fastas of the one to two genes for CNVSelectR

seqs <- read.csv("./Annotations_Gene.tsv", sep="")


pair_groups <- one_to_two_genes %>%
  group_by(fly, community) %>%
  mutate(n = n()) %>%
  filter(n == 2) %>%
  mutate(group = community,
         gene_group = paste0(fly,'_',gn))


pair_groups <- pair_groups[c('group','gene_group')]
prot_nucs <- seqs[c('gene_group','prot','nuc')]

pair_groups <- merge(pair_groups,prot_nucs,by='gene_group')

pair_groups <- pair_groups %>%
  mutate(fly = str_extract(gene_group, "^[^_]+"))

# write nucleotide and protein fasta files for the groups of pairs 
for (group_number in unique(pair_groups$group)) {
  group_rows <- pair_groups %>% filter(group == group_number)
  
  group_num <- as.numeric(group_number)
  group_num_orig <- group_num
  
  for (row_num in 1:nrow(group_rows)) {
    row <- group_rows[row_num,]
    fly_name <- gsub('-','_',row$fly)
    group_num <- paste0(group_num_orig, '_', fly_name)
    
    output_file <- paste0('./CNVSelectR/Connected_Eq_Nucleotide_Sequences/group_',group_num,'.fa')  
    cat(">", row$gene_group, '\n', row$nuc, "\n", file = output_file, append = T, sep = '')
    
    output_file <- paste0('./CNVSelectR/Connected_Eq_Protein_Sequences/group_',group_num,'.fa')  
    cat(">", row$gene_group, '\n', row$prot, "\n", file = output_file, append = T, sep = '')
    
  }
}



