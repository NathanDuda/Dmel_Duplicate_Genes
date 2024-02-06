

# require these to be in duplicate genes blastp output 


source("startup.R")
names <- c(
  'A1','A2','A3','A4','A5','A6','A7','AB8','AKA-017','AKA-018','B1','B2','B3','B4','B6','COR-014',
  'COR-018','COR-023','COR-025','GIM-012','GIM-024','ISO-1','JUT-008','JUT-011','KIE-094','LUN-004',
  'LUN-007','MUN-008','MUN-009','MUN-013','MUN-015','MUN-016','MUN-020','ORE','RAL-059','RAL-091',
  'RAL-176','RAL-177','RAL-375','RAL-426','RAL-737','RAL-855','SLA-001','STO-022','TEN-015','TOM-007',
  'TOM-008')

#dups <- read.csv("./Duplicate_Proteins.tsv", sep="")
#dups <- dups[c('dup')]
#new_columns_df <- as.data.frame(setNames(rep(list(0), length(names)), names))
#dups <- bind_cols(dups, new_columns_df)


all_genes <- read.csv("C:/Users/17735/Downloads/Dmel_Duplicate_Genes/Annotations_Gene.tsv", sep="")

all_prot_lengths <- all_genes[c('gene_group','prot')]
all_prot_lengths$nchar <- nchar(all_prot_lengths$prot)
all_prot_lengths <- all_prot_lengths[c('gene_group','nchar')]

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
      blastp <- read.delim(paste0('./Equivalent_Genes_Blastp_Output/Blastp_Output/q_',second_name,'_db_',name,'_blastp.tsv'), header=FALSE)
      
      colnames(blastp) <- c('qseqid', 'sseqid', 'evalue', 'pident', 'bitscore', 'length')
      blastp$qseqid <- paste0(second_name,'_',blastp$qseqid)
      blastp$sseqid <- paste0(name,'_',blastp$sseqid)
      
      rev_blastp <- read.delim(paste0('./Equivalent_Genes_Blastp_Output/Rev_Blastp_Output/q_',name,'_db_',second_name,'_blastp.tsv'), header=FALSE)
      colnames(rev_blastp) <- c('qseqid', 'sseqid', 'evalue', 'pident', 'bitscore', 'length') # PUT 
      rev_blastp$qseqid <- paste0(name,'_',rev_blastp$qseqid)
      rev_blastp$sseqid <- paste0(second_name,'_',rev_blastp$sseqid)
      
      
      # get protein lengths
      colnames(all_prot_lengths) <- c('qseqid','qlen')
      blastp <- merge(blastp,all_prot_lengths,by='qseqid')
      rev_blastp <- merge(rev_blastp,all_prot_lengths, by='qseqid')
      
      colnames(all_prot_lengths) <- c('sseqid','slen')
      blastp <- merge(blastp,all_prot_lengths,by='sseqid')
      rev_blastp <- merge(rev_blastp,all_prot_lengths,by='sseqid')
      
      # keep only quality hits 
      #blastp <- blastp %>% filter((bitscore >= 2*length) & (length>=100) & (pident >= 99))
      #rev_blastp <- rev_blastp %>% filter((bitscore >= 2*length) & (length>=100) & (pident >= 99)) # PUT 
      blastp <- blastp %>% 
        filter(length >= 100 & pident >= 99) %>%
        filter(case_when((qlen>slen  & ((.90*qlen)<length))~T, # require the hit be at least 90% of the longer gene 
                         (slen>qlen  & ((.90*slen)<length))~T,
                         (slen==qlen & ((.90*slen)<length))~T)) # IGNORES SEPARATE BLAST HITS TO SAME GN
      
      
      
      # keep only reciprocal q/db and db/q hits
      blastp$name_second_name <- paste0(blastp$sseqid,'_',blastp$qseqid)
      rev_blastp$name_second_name <- paste0(rev_blastp$qseqid,'_',rev_blastp$sseqid)
      
      blastp <- blastp[blastp$name_second_name %in% rev_blastp$name_second_name,]
      
      # now blastp contains both blastp and rev_blastp info 
      # if n different genes all hit to same 
      # 
      
      t <- blastp[c('qseqid','sseqid')]
      t <- t[!duplicated(t),]
      
      
      qseqid_dups <- t[t$sseqid %in% t$sseqid[duplicated(t$sseqid) | duplicated(t$sseqid, fromLast = TRUE)], ]
      
      q <- qseqid_dups %>%
        mutate(q_dup_group = paste0(second_name,'_',dense_rank(sseqid))) %>%
        mutate(s_dup_group = paste0(name,'_',dense_rank(qseqid))) %>%
        
        group_by(q_dup_group) %>%
        mutate(q_n_in_dup_group = n()) %>%
        ungroup() %>%
        
        group_by(s_dup_group) %>%
        mutate(s_n_in_dup_group = n()) %>%
        ungroup()
      
      
      
      
      blastp <- t[!duplicated(t),]
      colnames(blastp) <- c('q','db')
      
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
      
      
      # for q 
      
      blastp_q <- blastp %>%
        group_by(db) %>%
        arrange(q) %>%
        mutate(q_concatenated = paste(q, collapse = ',')) %>%
        mutate(q = if_else(row_number() == 1, q_concatenated, NA_character_)) %>%
        mutate(q = q_concatenated) %>%
        select(-q_concatenated) %>%
        ungroup()
      
      # for db 
      
      blastp_db <- blastp %>%
        group_by(q) %>%
        arrange(db) %>%
        mutate(db_concatenated = paste(db, collapse = ',')) %>%
        mutate(db = if_else(row_number() == 1, db_concatenated, NA_character_)) %>%
        mutate(db = db_concatenated) %>%
        select(-db_concatenated) %>%
        ungroup()
      
      # put into all_genes 
      
      all_genes <- all_genes %>%
        mutate(across(all_of(blastp_q$q_name), 
                      ~ if_else(gene_group %in% blastp_q$db, 
                                as.character(blastp_q$q[match(gene_group, blastp_q$db)]), 
                                as.character(.)))) %>%
        mutate(across(all_of(blastp_db$db_name), 
                      ~ if_else(gene_group %in% blastp_db$q, 
                                as.character(blastp_db$db[match(gene_group, blastp_db$q)]), 
                                as.character(.))))
      
      
      #rm(blastp)
    }
  }
  #rows that start with name 
  #all <- all_genes %>%
  #  filter(str_starts(gene_group, name)) %>%
  #  mutate_all(~ ifelse(. == 0, NA, .))
  #rev_result_df <- bind_rows(rev_result_df, all)
  
  p = p + 1
  print(p)  
}


all_genes <- all_genes[!grepl("^(I23|ZH26|N25|T29A|B59)", all_genes$gene_group), ]
write.table(all_genes,file='./Equivalent_Genes.tsv')



