

source("startup.R")



# all_genes <- read.csv("./Corrected_Equivalent_Genes.tsv", sep="")
# all_genes <- all_genes %>% distinct(.keep_all=T)

connected_dups <- read.csv("./connected_dups.tsv", sep="")

all_genes <- connected_dups %>%
  mutate(fly = sub("^(.*?)_.*", "\\1", node)) %>%
  group_by(community, fly) %>%
  summarise(gn = paste(node, collapse = ","), .groups = 'keep') %>%
  ungroup() %>%
  pivot_wider(., values_from = 'gn', names_from = 'fly')

# count the number of genes 
equiv_counts_matrix <- all_genes %>%
  mutate_at(vars(2:48), ~ case_when(.==0~0,
                                    .!=0~str_count(as.character(.), ",") + 1)) 

# get the one to ones 
one_to_ones <- equiv_counts_matrix %>%
  select(-community) %>%
  filter_all(all_vars(. %in% c(1)))


# make frequency distributions 
equiv_counts <- equiv_counts_matrix %>%
  mutate_all(., ~replace(., is.na(.), 0)) %>%
  pivot_longer(cols=c(2:ncol(equiv_counts_matrix)))

n_copies_freq_dist <- data.frame()
for (n_copies in c(1:3)){
  n_df <- equiv_counts %>%
    group_by(community) %>%
    summarise(flies_with_n_copies = sum(value == n_copies)) %>%
    mutate(n_copies = n_copies) %>%
    filter(flies_with_n_copies != 0)
  
  n_copies_freq_dist <- rbind(animate_df, n_df)
  print(n_copies)
}

ggplot(n_copies_freq_dist, aes(x = flies_with_n_copies)) +
  geom_histogram(binwidth = 1) +
  facet_wrap(.~n_copies)


# more than 1 copy frequency distribution
equiv_counts %>%
  group_by(community) %>%
  summarise(flies_with_n_copies = sum(value > 1)) %>%
  filter(flies_with_n_copies != 0) %>%
  
  ggplot(aes(x = flies_with_n_copies)) +
  geom_histogram(binwidth = 1)



# make pca plot
equiv_counts <- equiv_counts_matrix %>%
  mutate_all(., ~replace(., is.na(.), 0)) 

pca_result <- prcomp(as.data.frame(t(equiv_counts)))

pca_df <- data.frame(
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  Individuals = names(pca_result$x[, 1]))

ggplot(pca_df, aes(x = PC1, y = PC2, label = Individuals)) +
  geom_point() + 
  geom_text(size = 2, vjust = 2)

pca_plot <- pca_df %>%
  mutate(area = case_when(
    str_detect(Individuals, "B1|RAL|ORE|ISO|A1|A6|I23|B4") ~ 'North_America',
    str_detect(Individuals, "A2|B6") ~ 'South_America',
    str_detect(Individuals, "TEN|GIM|MUN|JUT|LUN|TOM|COR|A3|STO|KIE|A5|B3|SLA|AKA") ~ 'Europe',
    Individuals == "AB8" ~ 'Middle_East',
    Individuals %in% c("A4", 'B2') ~ 'Southish_Africa',  
    Individuals == "A7" ~ 'East_Asia' 
  ))

ggplot(pca_plot, aes(x = PC1, y = PC2, color = continent)) +
  geom_point()


library(ggrepel)
ggplot(t, aes(x = PC1, y = PC2, color = area, label = Individuals)) +
  geom_point(size = 2.5) +
  geom_text_repel(size = 2.5, max.overlaps = 20) +
  theme_bw()

#############################


fly_variances <- apply(equiv_counts, 2, var)
t <- as.data.frame(fly_variances)


gn_variances <- apply(equiv_counts, 1, var)
g <- as.data.frame(gn_variances)




###



t <- equiv_counts %>%
  mutate(n_zero = rowSums(select(., 1:47) == 0)) %>%
  arrange(n_zero) %>%
  summarise(across(everything(), ~unique(.)))



