

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
equiv_counts <- all_genes %>%
  mutate_at(vars(2:48), ~ case_when(.==0~0,
                                    .!=0~str_count(as.character(.), ",") + 1)) 


# get the one to ones 
one_to_ones <- equiv_counts %>%
  select(-community) %>%
  filter_all(all_vars(. %in% c(1)))

# make frequency distributions 


equiv_counts <- equiv_counts %>%
  mutate_all(., ~replace(., is.na(.), 0)) %>%
  pivot_longer(cols=c(2:ncol(equiv_counts)))



animate_df <- data.frame()
for (n_copies in c(1:3)){
  n_df <- equiv_counts %>%
    group_by(community) %>%
    summarise(flies_with_n_copies = sum(value == n_copies)) %>%
    mutate(n_copies = n_copies) %>%
    filter(flies_with_n_copies != 0)
  
  animate_df <- rbind(animate_df, n_df)
  print(n_copies)
}


ggplot(animate_df, aes(x = flies_with_n_copies)) +
  geom_histogram(binwidth = 1) +
  facet_wrap(.~n_copies)


####
  

t <- equiv_counts %>%
  group_by(community) %>%
  summarise(flies_with_n_copies = sum(value == 2)) %>%
  
  ggplot(aes(x = flies_with_n_copies)) +
  geom_histogram(binwidth = 1)





n_more_than_1 <- equiv_counts %>%
  mutate_all(., ~replace(., is.na(.), 0)) %>%
  rowwise() %>%
  mutate(n_more_than_1 = sum(c_across(everything()) > 1))

n_more_than_1 %>% filter(n_more_than_1 != 0) %>%
  ggplot(., aes(x = n_more_than_1)) +
    geom_histogram(binwidth = 1)


equiv_counts <- equiv_counts %>%
  mutate_all(., ~replace(., is.na(.), 0))
  


animate_df <- data.frame()
for (n_copies in c(1:5)){
  n_df <- equiv_counts %>%
    rowwise() %>%
    mutate(flies_with_n_copies = sum(c_across(everything()) == n_copies)) %>% # can change to > 
    filter(flies_with_n_copies != 0) %>%
    mutate(n_copies = n_copies)
  
  animate_df <- rbind(animate_df, n_df)
  print(n_copies)
}


library(gganimate)

animate_df %>% 
  filter(n_copies == 1) %>%
  ggplot(., aes(x = flies_with_n_copies)) +
    geom_histogram(binwidth = 1)



ggplot(animate_df, aes(x = flies_with_n_copies)) +
  geom_histogram(binwidth = 1) +
  transition_time(n_copies)





##### animate



#####



# make a pca 


transposed_df <- as.data.frame(t(equiv_counts))

# Step 2: Perform PCA
pca_result <- prcomp(transposed_df)


# Step 3: Create a dataframe containing the PCA results
pca_df <- data.frame(
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  Individuals = names(pca_result$x[, 1])  # Assuming column names are individuals
)

# Step 4: Create the PCA plot using ggplot2
ggplot(pca_df, aes(x = PC1, y = PC2, label = Individuals)) +
  geom_point() + 
  geom_text(size = 2, vjust = 2)

t <- pca_df %>%
  mutate(area = case_when(
    str_detect(Individuals, "B1|RAL|ORE|ISO|A1|A6|I23|B4") ~ 'North_America',
    str_detect(Individuals, "A2|B6") ~ 'South_America',
    str_detect(Individuals, "TEN|GIM|MUN|JUT|LUN|TOM|COR|A3|STO|KIE|A5|B3|SLA|AKA") ~ 'Europe',
    Individuals == "AB8" ~ 'Middle_East',
    Individuals %in% c("A4", 'B2') ~ 'Southish_Africa',  
    Individuals == "A7" ~ 'East_Asia' 
  ))



ggplot(t, aes(x = PC1, y = PC2, color = continent)) +
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



