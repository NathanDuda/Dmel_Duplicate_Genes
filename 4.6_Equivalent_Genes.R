



source("startup.R")
names <- c(
  'A1','A2','A3','A4','A5','A6','A7','AB8','AKA-017','AKA-018','B1','B2','B3','B4','B6','COR-014',
  'COR-018','COR-023','COR-025','GIM-012','GIM-024','ISO-1','JUT-008','JUT-011','KIE-094','LUN-004',
  'LUN-007','MUN-008','MUN-009','MUN-013','MUN-015','MUN-016','MUN-020','ORE','RAL-059','RAL-091',
  'RAL-176','RAL-177','RAL-375','RAL-426','RAL-737','RAL-855','SLA-001','STO-022','TEN-015','TOM-007',
  'TOM-008')

dups <- read.csv("C:/Users/17735/Downloads/Dmel_Duplicate_Genes/Duplicate_Proteins.tsv", sep="")
dups <- dups[c('dup')]
new_columns_df <- as.data.frame(setNames(rep(list(0), length(names)), names))
dups <- bind_cols(dups, new_columns_df)


all_genes <- read.csv("C:/Users/17735/Downloads/Dmel_Duplicate_Genes/Annotations_Gene.tsv", sep="")
all_genes <- all_genes[c('gene_group')]
new_columns_df <- as.data.frame(setNames(rep(list(0), length(names)), names))
all_genes <- bind_cols(all_genes, new_columns_df)


# make empty dataframe for output
all_eq_genes <- as.data.frame(matrix(ncol=2,nrow=0))
result_df <- data.frame()


p = 0 


# loop over all names 
for (name in names) {
  for (second_name in names) {
    if (name < second_name) {
      blastp <- read.delim(paste0('./Equivalent_Genes_Blastp_Output/Rev_Blastp_Output/q_',second_name,'_db_',name,'_blastp.tsv'))
      blastp <- blastp[,c(1,2)]
      blastp <- blastp[!duplicated(blastp),]
      colnames(blastp) <- c('q','db')
      
      blastp$q  <- paste0(second_name,'_',blastp$q)
      blastp$db <- paste0(name,'_',blastp$db)
      
      
      blastp <- blastp %>%
        filter((q %in% all_genes$gene_group) | (db %in% all_genes$gene_group)) %>%
        mutate(n_duplicated = case_when(
          (q %in% all_genes$gene_group) & !(db %in% all_genes$gene_group) ~ 'q',
          !(q %in% all_genes$gene_group) & (db %in% all_genes$gene_group) ~ 'db',
          (q %in% all_genes$gene_group) & (db %in% all_genes$gene_group) ~ 'two')) %>%
        distinct(
          q = if_else(n_duplicated %in% c('q', 'two'), q, q),
          db = if_else(n_duplicated %in% c('db', 'two'), db, db),
          .keep_all = TRUE
        ) 
      
      
      blastp <- blastp %>% filter(n_duplicated == 'two') %>%
        mutate(q_name =  gsub('-','.',(str_extract(q, "^[^_]+")))) %>%
        mutate(db_name = gsub('-','.',(str_extract(db, "^[^_]+")))) 

      blastp <- blastp %>%
        group_by(db) %>%
        mutate(q_concatenated = paste(q, collapse = ',')) %>%
        mutate(q = if_else(row_number() == 1, q_concatenated, NA_character_)) %>%
        select(-q_concatenated) %>%
        ungroup()
      
      
      all_genes <- all_genes %>%
        mutate(across(all_of(blastp$q_name), 
                      ~ if_else(gene_group %in% blastp$db, 
                                blastp$q[match(gene_group, blastp$db)], 
                                blastp$q[match(gene_group, blastp$db)])))
      
      rm(blastp)
    }
  }
 #rows that start with name 
  all <- all_genes %>%
    filter(str_starts(gene_group, name))
  result_df <- bind_rows(result_df, all)
  
  p = p + 1
  print(p)
}

#write.table(result_df,file='./Equiv_All_Genes.tsv')











source("startup.R")
names <- c(
  'A1','A2','A3','A4','A5','A6','A7','AB8','AKA-017','AKA-018','B1','B2','B3','B4','B6','COR-014',
  'COR-018','COR-023','COR-025','GIM-012','GIM-024','ISO-1','JUT-008','JUT-011','KIE-094','LUN-004',
  'LUN-007','MUN-008','MUN-009','MUN-013','MUN-015','MUN-016','MUN-020','ORE','RAL-059','RAL-091',
  'RAL-176','RAL-177','RAL-375','RAL-426','RAL-737','RAL-855','SLA-001','STO-022','TEN-015','TOM-007',
  'TOM-008')

dups <- read.csv("C:/Users/17735/Downloads/Dmel_Duplicate_Genes/Duplicate_Proteins.tsv", sep="")
dups <- dups[c('dup')]
new_columns_df <- as.data.frame(setNames(rep(list(0), length(names)), names))
dups <- bind_cols(dups, new_columns_df)


all_genes <- read.csv("C:/Users/17735/Downloads/Dmel_Duplicate_Genes/Annotations_Gene.tsv", sep="")
all_genes <- all_genes[c('gene_group')]
new_columns_df <- as.data.frame(setNames(rep(list(0), length(names)), names))
all_genes <- bind_cols(all_genes, new_columns_df)


# make empty dataframe for output
all_eq_genes <- as.data.frame(matrix(ncol=2,nrow=0))








# repeat for reverse 
rev_result_df <- data.frame()
p = 0

for (name in names) {
  for (second_name in names) {
    if (name > second_name) {
      blastp <- read.delim(paste0('./Equivalent_Genes_Blastp_Output/Blastp_Output/q_',second_name,'_db_',name,'_blastp.tsv'))
      blastp <- blastp[,c(1,2)]
      blastp <- blastp[!duplicated(blastp),]
      colnames(blastp) <- c('q','db')
      
      blastp$q  <- paste0(second_name,'_',blastp$q)
      blastp$db <- paste0(name,'_',blastp$db)
      
      
      blastp <- blastp %>%
        filter((q %in% all_genes$gene_group) | (db %in% all_genes$gene_group)) %>%
        mutate(n_duplicated = case_when(
          (q %in% all_genes$gene_group) & !(db %in% all_genes$gene_group) ~ 'q',
          !(q %in% all_genes$gene_group) & (db %in% all_genes$gene_group) ~ 'db',
          (q %in% all_genes$gene_group) & (db %in% all_genes$gene_group) ~ 'two')) %>%
        distinct(
          q = if_else(n_duplicated %in% c('q', 'two'), q, q),
          db = if_else(n_duplicated %in% c('db', 'two'), db, db),
          .keep_all = TRUE
        ) 
      
      
      blastp <- blastp %>% filter(n_duplicated == 'two') %>%
        mutate(q_name =  gsub('-','.',(str_extract(q, "^[^_]+")))) %>%
        mutate(db_name = gsub('-','.',(str_extract(db, "^[^_]+")))) 
      
      blastp <- blastp %>%
        group_by(db) %>%
        mutate(q_concatenated = paste(q, collapse = ',')) %>%
        mutate(q = if_else(row_number() == 1, q_concatenated, NA_character_)) %>%
        select(-q_concatenated) %>%
        ungroup()
      
      
      all_genes <- all_genes %>%
        mutate(across(all_of(blastp$q_name), 
                      ~ if_else(gene_group %in% blastp$db, 
                                blastp$q[match(gene_group, blastp$db)], 
                                blastp$q[match(gene_group, blastp$db)])))
      
      rm(blastp)
    }
  }
  #rows that start with name 
  all <- all_genes %>%
    filter(str_starts(gene_group, name))
  rev_result_df <- bind_rows(rev_result_df, all)
  
  p = p + 1
  print(p)  
}
  

#write.table(rev_result_df,file='./Rev_Equiv_All_Genes.tsv')

  


######################################################################################

# read in output files
all_genes <- read.csv("./Equiv_All_Genes.tsv", sep="")
rev_all_genes <- read.csv("./Rev_Equiv_All_Genes.tsv", sep="")


# change 0's to NA
all_genes <- all_genes %>% mutate_all(~ ifelse(. == 0, NA, .))



t <- all_genes %>%
  merge(.,rev_all_genes,by='gene_group') %>%
  transmute(names[1] = coacross(starts_with(names[1])))

  


t <- merge(all_genes, rev_all_genes, by = 'gene_group')

for (col_name in names) {
  t <- t %>% mutate(!!col_name := coalesce(!!sym(gsub('-','.',paste0(col_name,'.x'))), !!sym(gsub('-','.',paste0(col_name,'.y')))))
}

t <- t %>% select(-ends_with('.x'), -ends_with('.y'))



write.table(t,file='Equivalent_Genes.tsv')



#



equiv_genes <- read.csv("./Equivalent_Genes.tsv", sep="")
dups <- read.csv("./Duplicate_Proteins.tsv", sep="")


# get the number of duplicates each gene has in its family
dups <- dups %>%
  select(dup, dup_family) %>%
  group_by(dup_family) %>%
  mutate(n_copies_in_fam = n())

t <- merge(equiv_genes, dups, by.x = 'gene_group', by.y = 'dup')


t2 <- t %>%
  mutate_at(vars(2:48), ~ str_count(as.character(.), ",")) %>%
  mutate_at(vars(2:48), ~ . + 1)



plotz <- t2 %>%
  mutate_at(vars(2:48), ~ coalesce(., n_copies_in_fam)) %>%
  rowwise() %>%
  mutate(identitcal_n_copies_across = case_when(all(c_across(2:48) == n_copies_in_fam) ~ T,
                                                T ~ F)) %>%
  mutate(duplicated_across = case_when(all(c_across(2:48) > 1) ~ T,
                                       T ~ F))

# n duplicated same amount of copies 

result <- t2 %>%
  mutate_at(vars(2:48), ~ coalesce(., n_copies_in_fam)) %>%
  rowwise() %>%
  mutate(n_same_n_duplicated = sum(across(2:48) == n_copies_in_fam, na.rm = TRUE))


plot <- as.data.frame(table(result$n_same_n_duplicated))


gg <- ggplot(plot,aes(x=Var1,y=Freq)) +
  geom_bar(stat='identity') +
  #scale_x_discrete(limits=factor(1:52)) + # not 52 without crappy assemblies 
  xlab('Number of flies duplicate is present in') +
  ylab('Freq') +
  theme_bw()

###


# n duplicated 

result2 <- t2 %>%
  mutate_at(vars(2:48), ~ coalesce(., n_copies_in_fam)) %>%
  rowwise() %>%
  mutate(n_same_n_duplicated = sum(across(2:48) > 1, na.rm = TRUE))


plot2 <- as.data.frame(table(result2$n_same_n_duplicated))


gg2 <- ggplot(plot2,aes(x=Var1,y=Freq)) +
  geom_bar(stat='identity') +
  #scale_x_discrete(limits=factor(1:52)) + # not 52 without crappy assemblies 
  xlab('Number of flies duplicate is present in') +
  ylab('Freq') +
  theme_bw()



gg
gg2


###################


# repeating same gene for the 50 diff flies 
### but not repeating rows because same gene may hit to multiple others genes 
has_duplicates <- any(duplicated(as.vector(unlist(df))))

############


#


t <- result_df %>%
  mutate_all(~ ifelse(. == 0, 'a', .)) %>%
  rowwise() %>%
  mutate(comma_count = sum(str_detect(c_across(2:48), ",")))

table(t$comma_count)






#



t2 <- result_df %>%
  mutate_all(~ str_count(as.character(.), ",")) 


t2 <- t2 %>%
  rowwise() %>%
  mutate(same_comma_n = case_when(length(unique(c_across(2:48))) == 1 ~ 't',
                                  T ~ 'f')) 



#


#

##################################


dups <- read.csv("C:/Users/17735/Downloads/Dmel_Duplicate_Genes/Equiv_Dups.tsv", sep="")
all_genes <- read.csv("C:/Users/17735/Downloads/Dmel_Duplicate_Genes/Equiv_All_Genes.tsv", sep="")


t <- all_genes %>%
  mutate(n_flies = rowSums(select(., A1:TOM.008)))
t$n_flies <- t$n_flies + 1



table(t$n_flies)


plot_freq_dist <- as.data.frame(table(t$n_flies))

ggplot(plot_freq_dist,aes(x=Var1,y=Freq)) +
  geom_bar(stat='identity') +
  #scale_x_discrete(limits=factor(1:52)) + # not 52 without crappy assemblies 
  xlab('Number of flies duplicate is present in') +
  ylab('Freq') +
  theme_bw()





