



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
      
      
      blastp <- blastp %>% 
        filter(n_duplicated == 'two') %>%
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
    filter(str_starts(gene_group, name)) %>%
    mutate_all(~ ifelse(. == 0, NA, .))
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




t <- full_join(all_genes, rev_all_genes, by = 'gene_group')

for (col_name in names) {
  t <- t %>% mutate(!!col_name := coalesce(!!sym(gsub('-','.',paste0(col_name,'.x'))), !!sym(gsub('-','.',paste0(col_name,'.y')))))
}

t <- t %>% select(-ends_with('.x'), -ends_with('.y'))



write.table(t,file='Equivalent_Genes.tsv')



#




