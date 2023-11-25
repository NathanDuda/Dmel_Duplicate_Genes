



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


# make empty dataframe for output
all_eq_genes <- as.data.frame(matrix(ncol=2,nrow=0))


p = 0 


# loop over all names 
for (name in names) {
  for (second_name in names) {
    if (name == second_name) {break}
    if (name < second_name) {
      blastp <- read.delim(paste0('./Equivalent_Genes_Blastp_Output/Rev_Blastp_Output/q_',second_name,'_db_',name,'_blastp.tsv'))
      blastp <- blastp[,c(1,2)]
      blastp <- blastp[!duplicated(blastp),]
      colnames(blastp) <- c('q','db')
      
      blastp$q  <- paste0(second_name,'_',blastp$q)
      blastp$db <- paste0(name,'_',blastp$db)
    }
    if (name > second_name) {
      blastp <- read.delim(paste0('./Equivalent_Genes_Blastp_Output/Blastp_Output/q_',second_name,'_db_',name,'_blastp.tsv'))
      blastp <- blastp[,c(1,2)]
      blastp <- blastp[!duplicated(blastp),]
      colnames(blastp) <- c('q','db')
      
      blastp$db <- paste0(name,'_',blastp$db)
      blastp$q  <- paste0(second_name,'_',blastp$q)
    }
    
    
    blastp <- blastp %>%
      filter((q %in% dups$dup) | (db %in% dups$dup)) %>%
      mutate(n_duplicated = case_when(
        (q %in% dups$dup) & !(db %in% dups$dup) ~ 'q',
        !(q %in% dups$dup) & (db %in% dups$dup) ~ 'db',
        (q %in% dups$dup) & (db %in% dups$dup) ~ 'two')) %>%
      distinct(
        q = if_else(n_duplicated %in% c('q', 'two'), q, NA),
        db = if_else(n_duplicated %in% c('db', 'two'), db, NA),
        .keep_all = TRUE
      ) %>%
      
      
      filter(n_duplicated == 'two') %>%
        mutate(q_name =  gsub('-','.',(str_extract(q, "^[^_]+")))) %>%
        mutate(db_name = gsub('-','.',(str_extract(db, "^[^_]+")))) %>%
        mutate(across(1, ~ ifelse(duplicated(.), NA, .))) %>%
        mutate(across(db, ~ ifelse(duplicated(.), NA, .)))

    
    
    dups <- dups %>%
      mutate(across(all_of(blastp$db_name), ~if_else(dup %in% blastp$q, . + 1, .))) %>%
      mutate(across(all_of(blastp$q_name), ~if_else(dup %in% blastp$db, . + 1, .)))
    
    p = p + 1
    print(p)
    rm(blastp)
  }}

write.table(dups,file='./Equiv_Dups.tsv')


