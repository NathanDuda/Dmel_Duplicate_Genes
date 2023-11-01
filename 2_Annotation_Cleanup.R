# Title: 2_Annotation_Cleanup_and_Blast_Run.R
# Author: Nathan Duda
# Date: 10/28/2023
# Purpose: Format the annotations into fasta files for each fly.
#          Format all annotations into a table with gene id, prot sequence, chromosome, location, and other information.


source("startup.R")

# read in list with all fly names 
fly_name_list <- read.table("./Individuals_List.txt", quote="\"", comment.char="")

# create empty dataframe for output
all_annotations <- as.data.frame(matrix(ncol=13,nrow=0))

# loop over all flies 
for (row in 1:nrow(fly_name_list)){
  
  # get name of current fly 
  name <- fly_name_list[row,1]
  
  # import the dataframe 
  aa <- read.csv(paste0("./Augustus_Output/aa_fastas/",name,"_aa.fasta"), sep="", header = F)
  
  # format the fasta into a dataframe 
  colnames(aa)[1] <- 'aa'
  aa <- aa %>%
    mutate(Group = cumsum(grepl("^>", aa))) %>%
    group_by(Group) %>%
    mutate(prot = paste(aa, collapse = "")) %>%
    ungroup() %>%
    select(-Group)
  aa[1,1] <- '>'
  aa <- aa[grepl(">", aa$aa), ]
  aa$aa <- '>'
  aa$prot <- gsub('>','',aa$prot)
  
  # get protein length
  aa$nchar <- nchar(aa$prot)
  
  
  # import the number of genes on each chromosome to get the chromosomes of each gene (since theyre in order)
  aa_per_chrom <- read.table(paste0("./Augustus_Output/aa_per_chrom/",name,"_aa_per_chrom.txt"), quote="\"", comment.char="")
  aa_per_chrom[aa_per_chrom$V1 == 0,] <- NA
  aa_per_chrom <- na.omit(aa_per_chrom)
  
  chrom <- c('2L','2R','3L','3R','X')
  
  aa_per_chrom$chrom <- chrom
  aa_per_chrom$V1 <- cumsum(aa_per_chrom$V1)
  
  aa_per_chrom$gene_group <- paste0('g',aa_per_chrom$V1)
  
  
  # import annotation output from augustus
  annotations <- read.delim(paste0("./Augustus_Output/",name,"_annotations.gff"), header=FALSE, comment.char="#")
  
  # format the annotation dataframe
  annotations <- annotations[,-c(1,2)]
  annotations$V9 <- str_extract(annotations$V9, "g\\d+")
  colnames(annotations) <- c('type','start','end','n','strand','s','gene_group')
  
  # get the length of each annotation
  annotations$len <- annotations$end - annotations$start
  
  # get the approximate length of each protein by taking the difference in length between the
  # stop and start codons and subtract out the length of all introns, then dividing by 3
  annotations <- annotations %>%
    mutate(approx_expected_prot_len = case_when((strand=='-' & type=='stop_codon') ~ -start,
                                                (strand=='-' & type=='intron') ~ -len,
                                                (strand=='-' & type=='start_codon') ~ start,
                                                
                                                (strand=='+' & type=='start_codon') ~ -start,
                                                (strand=='+' & type=='intron') ~ -len,
                                                (strand=='+' & type=='stop_codon') ~ start,
                                                T ~ 0)) %>%
    group_by(gene_group) %>%
    mutate(approx_expected_prot_len = sum(approx_expected_prot_len) / 3) %>%
    ungroup()
  
  
  # get the chromosomes for each gene
  annotations <- left_join(annotations,aa_per_chrom,by='gene_group')
  annotations <- annotations %>% fill(chrom, .direction = "up") %>% select(-V1)
  
  
  # get the proteins for each gene 
  aa$gene_group <- paste0('g',1:nrow(aa))
  annotations <- merge(annotations,aa,by='gene_group')
  annotations <- annotations %>% select(-aa)
  
  ####################################
  # write proteins into fasta file 
  
  # Convert the dataframe to FASTA-like format
  annotations_subset <- annotations[c('gene_group','prot')]
  annotations_subset <- annotations_subset[!duplicated(annotations_subset),]
  
  prot_fasta <- character(nrow(annotations_subset) * 2)
  prot_fasta[c(TRUE, FALSE)] <- paste0(">", annotations_subset$gene_group)
  prot_fasta[c(FALSE, TRUE)] <- annotations_subset$prot
  
  # write into fasta file 
  writeLines(prot_fasta,paste0('./Prot_Fastas/',name,'_prot.fasta'))
  
  ####################################
  
  # add fly name to the annotation dataframe
  annotations$fly <- name 
  
  # add annotation to dataframe with all annotations
  all_annotations <- rbind(all_annotations,annotations)

  print(row)
}

# write to table
write.table(all_annotations,'Annotations.tsv')
