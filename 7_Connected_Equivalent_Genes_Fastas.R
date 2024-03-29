

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


min(table(connected_dups$fly))
max(table(connected_dups$fly))
# each fly has 9,224 to 12,055 genes in a community 


community_fly_counts <- as.data.frame(table(connected_dups$fly, connected_dups$community))

anys_to_anys <- community_fly_counts %>%
  group_by(Var2) %>%
  filter(all(Freq %in% c(0, 1, 2)))
length(unique(anys_to_anys$Var2))
# 24,720 out of 25,005 have all 0, 1, or 2 copies 

one_to_zeros <- community_fly_counts %>%
  group_by(Var2) %>%
  filter(all(Freq %in% c(0, 1)))
length(unique(one_to_zeros$Var2))
# 23,713 out of 25,005 have all 0 or 1 copies 

only_one_fly <- community_fly_counts %>% filter(Freq != 0) 
only_one_fly <- as.data.frame(table(only_one_fly$Var2))
only_one_fly %>% summarize(sum(Freq == 1))
# 97 out of 25,005 have only one fly in the community (one dup)

single_copies <- one_to_zeros %>%
  group_by(Var2) %>%
  summarise(count_ones = sum(Freq == 1))
table(single_copies$count_ones)
# 5,678 communities are single-copy pairs 
# 816 one to ones out of the 25,005 total communities 

same_n_copies <- community_fly_counts %>%
  group_by(Var2) %>%
  summarise(distinct_values = n_distinct(Freq))
sum(same_n_copies$distinct_values == 1)
# 817 communities have the same number of copies in every fly
community_fly_counts %>%
  group_by(Var2) %>%
  filter(all(Freq == 2))
# one community has all flies with 2 copies 

copy_variation <- community_fly_counts %>%
  group_by(Var2) %>%
  summarise(sd = sd(Freq))
ggplot(copy_variation, aes(y = sd, x = 0)) + geom_violin()
median(copy_variation$sd)
# median value of copy number standard deviation is 0.2820567
mean(copy_variation$sd)
# mean value of copy number standard deviation is 0.3020396


one_to_twos <- anys_to_anys %>%
  group_by(Var2) %>%
  filter(any(Freq == 2))
length(unique(one_to_twos$Var2))
# 1007 communities have 0, 1, or 2 copies with at least one 2 


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



