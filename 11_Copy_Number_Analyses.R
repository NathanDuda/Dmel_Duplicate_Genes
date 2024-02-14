

source("startup.R")



all_genes <- read.csv("./Corrected_Equivalent_Genes.tsv", sep="")

all_genes <- all_genes %>% distinct(.keep_all=T)

# count the number of genes 
equiv_counts <- all_genes %>%
  mutate_at(vars(1:47), ~ case_when(.==0~0,
                                    .!=0~str_count(as.character(.), ",") + 1)) 

# get the one to ones 
one_to_ones <- equiv_counts %>%
  filter_all(all_vars(. %in% c(1)))


row_counts <- apply(equiv_counts, 1, function(row) sum(row == 1) == ncol(equiv_counts))




