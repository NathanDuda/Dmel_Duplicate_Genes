


source('startup.R')


# read in genes with their seqyences and information 
all_annotations <- read.csv("./Annotations_Gene.tsv", sep="")
all_annotations <- all_annotations[c('gene_group','chrom','prot','nchar','fly')]

# read in global alignment scores
scores <- read.csv("./Global_Alignment_Output/all_global_align.tsv", row.names=NULL, sep="")
scores <- scores %>% select(-row.names) %>% distinct() %>% na.omit()
colnames(scores) <- c('gene_group.x','gene_group.y','score')


# merge each gene's information with their score
colnames(all_annotations) <- c('gene_group.x','chrom.x','prot.x','nchar.x','fly.x')
gene_pairs <- merge(scores,all_annotations,by='gene_group.x')

colnames(all_annotations) <- c('gene_group.y','chrom.y','prot.y','nchar.y','fly.y')
gene_pairs <- merge(gene_pairs,all_annotations,by='gene_group.y')


# get duplicate genes by requiring a protein length of at least 100 aa 
gene_pairs <- gene_pairs %>% 
  mutate(shorter_length = pmin(nchar.x,nchar.y)) %>%
  mutate(max_length_diff = 0) %>% # floor(shorter_length * 0.01) # for 99% match, now its 100%
  mutate(min_score = shorter_length - max_length_diff) %>%
  filter(score >= min_score) %>%
  filter(nchar.x >= 100 & nchar.y >= 100)

# get duplicate gene families 

fly_name_list <- read.table("./Individuals_List_47.txt", quote="\"", comment.char="")
dup_families <- data.frame()

# loop over all flies 
for (row in 1:nrow(fly_name_list)){
  
  # get name of current fly
  name <- fly_name_list[row,1]
  
  dup1_dup2 <- gene_pairs[gene_pairs$fly.x==name,c('gene_group.x','gene_group.y')]
  
  edge_list <- data.frame(source = dup1_dup2$gene_group.x, target = dup1_dup2$gene_group.y)
  
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


# merge gene information with duplicate families
colnames(all_annotations) <- c('dup','chrom','prot','nchar','fly')
dups <- merge(all_annotations,dup_families,by=c('fly','dup'))

# write to file
write.table(dups,file='Duplicate_Gene_Families.tsv')


