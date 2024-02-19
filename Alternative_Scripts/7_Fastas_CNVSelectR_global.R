
# global alignment


# get only duplicate pairs
# align duplicate pairs
# have each file with all equivalent duplicate pair(s)
source("startup.R")


seqs <- read.csv("./Annotations_Gene.tsv", sep="")


all_genes <- read.csv("./Equivalent_Genes.tsv", sep="")

# keep only the first two copies for genes that have more than 2 copies
remove_after_second_comma <- function(x) {
  parts <- unlist(strsplit(x, ",")) 
  if (length(parts) > 2) {result <- paste(parts[1:2], collapse = ",") } 
  else {result <- x}
  return(result) }

all_genes[] <- lapply(all_genes, function(x) sapply(x, remove_after_second_comma))

# count the number of genes 
equiv_counts <- all_genes %>%
  mutate_at(vars(2:48), ~ case_when(.==0~0,
                                    .!=0~str_count(as.character(.), ",") + 1)) 

# keep only duplicate pairs
pairs <- rowSums(equiv_counts[, 2:48] > 2) == 0
equiv_counts <- equiv_counts[pairs, ]


# get the duplicate status of the reference gene using the global alignments
dups <- read.csv("./Duplicate_Proteins.tsv", sep="")
dups <- dups[c('fly','dup','dup_chrom','dup_prot','dup_nchar','dup_family')]
colnames(dups) <- c('fly','dup','chrom','prot','nchar','dup_family')


dup_pairs <- dups %>%
  group_by(fly, dup_family) %>%
  mutate(n_copies = n()) %>%
  filter(n_copies==2)

dup_pairs_combined <- dup_pairs %>%
  group_by(fly, dup_family) %>%
  summarise(dup_pair = paste(dup, collapse = ','), .groups='keep') #%>%
  
  # get the reversed duplicate pair 
  #mutate(dup_pair_reversed = dup_pair) %>%
  #separate(dup_pair_reversed, into = c('second','first'), sep = ",") %>%
  #mutate(dup_pair_reversed = paste0(first,',',second)) %>%
  #select(-first,-second)


# place the gene_group gene into the correct column of the equivalent genes dataframe 
all_genes <- all_genes %>%
  mutate(fly = str_extract(gene_group, "^[^_]+"))

fly_name_list <- read.table("./Individuals_List_47.txt", quote="\"", comment.char="")
for (row in 1:nrow(fly_name_list)){
  fly_name <- fly_name_list[row,1]
  fly_name <- gsub('-','.',fly_name)
  all_genes <- all_genes %>%
    mutate(
      {{fly_name}} := case_when(fly == {{fly_name}} ~ gene_group,
                                T ~ get({{fly_name}})))
}

# make the new gene_group genes and any other missed genes into duplicate pairs
# found by global alignment if case not found by equivalent gene blast
dup_pairs <- merge(dup_pairs,dup_pairs,by=c('fly','dup_family'))
dup_pairs <- dup_pairs[dup_pairs$dup.x != dup_pairs$dup.y,]

all_genes <- all_genes %>%
  select(-gene_group, -fly) %>%
  mutate_at(vars(1:47), ~ case_when(
    . %in% dup_pairs$dup.x ~ paste0(., ',', dup_pairs$dup.y[match(., dup_pairs$dup.x)]),
    TRUE ~ .))

# keep equivalent genes only if they are duplicate pairs found by global alignment
pair_groups <- all_genes %>%
  distinct() %>%
  mutate_at(vars(1:47), ~ifelse(. %in% dup_pairs_combined$dup_pair | . %in% dup_pairs_combined$dup_pair_reversed, ., NA)) %>%
  distinct()


# format pairs
pair_groups <- pair_groups %>%
  mutate(group = 1:nrow(pair_groups)) %>%
  pivot_longer(cols=c(1:47),
               values_to = 'gene_group') %>%
  na.omit() %>%
  separate_rows(gene_group, sep = ',')


pair_groups <- pair_groups[c('group','gene_group')]
prot_nucs <- seqs[c('gene_group','prot','nuc')]

pair_groups <- merge(pair_groups,prot_nucs,by='gene_group')


# remove pairs that contain genes that are a part of multiple duplicate pairs
##pair_groups <- pair_groups %>%
##  group_by(gene_group) %>%
##  mutate(count = n())


#

pair_groups <- pair_groups %>%
  mutate(fly = str_extract(gene_group, "^[^_]+"))



# write nucleotide and protein fasta files for the groups of pairs 
for (group_number in 2:length(unique(pair_groups$group))) {
  group_rows <- pair_groups %>% filter(group == group_number)
  
  group_num <- as.numeric(group_number)
  group_num_orig <- group_num
  
  for (row_num in 1:nrow(group_rows)) {
    row <- group_rows[row_num,]
    fly_name <- gsub('-','_',row$fly)
    group_num <- paste0(group_num_orig, '_', fly_name)
    
    output_file <- paste0('./CNVSelectR/Nucleotide_Sequences/group_',group_num,'.fa')  
    cat(">", row$gene_group, '\n', row$nuc, "\n", file = output_file, append = T, sep = '')
    
    output_file <- paste0('./CNVSelectR/Protein_Sequences/group_',group_num,'.fa')  
    cat(">", row$gene_group, '\n', row$prot, "\n", file = output_file, append = T, sep = '')
    
  }
}




