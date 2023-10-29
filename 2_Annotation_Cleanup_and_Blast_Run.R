# Title: 2_Annotation_Cleanup_and_Blast_Run.R
# Author: Nathan Duda
# Date: 10/28/2023
# Purpose: Format the annotations into fasta files for each fly.
#          Format all annotations into a table with gene id, prot sequence, chromosome, location, and other information.
#          Run blast on all proteins of flies 


source("startup.R")

# function to write table with gene ids and sequences into fasta format 
format_protein_fasta <- function(gene_data) {
  fasta_lines <- paste(
    ">", gene_data$gene_group,
    "\n", gene_data$prot,
    sep = ""
  )
  return(paste(fasta_lines, collapse = "\n"))
}

# read in list with all fly names 
fly_name_list <- read.table("./Individuals_List.txt", quote="\"", comment.char="")

# create empty dataframe for output
all_annotations <- as.data.frame(matrix(ncol=13,nrow=0))

# loop over all flies 
for (row in 1:nrow(fly_name_list)){
  
  # get name of current fly 
  name <- fly_name_list[row,1]
  
  # import the dataframe 
  aa <- read.csv(paste0("C:/Users/17735/Downloads/Dmel_Duplicate_Genes/Augustus_Output/aa_fastas/",name,"_aa.fasta"), sep="")
  
  # format the fasta into a dataframe 
  colnames(aa) <- 'aa'
  aa <- aa %>%
    mutate(Group = cumsum(grepl("^>", aa))) %>%
    group_by(Group) %>%
    mutate(prot = paste(aa, collapse = "")) %>%
    ungroup() %>%
    select(-Group)
  aa[1,1] <- '>'
  aa <- aa[aa$aa == '>',]
  aa$prot <- gsub('>','',aa$prot)
  
  # get protein length
  aa$nchar <- nchar(aa$prot)
  
  
  # import the number of genes on each chromosome to get the chromosomes of each gene (since theyre in order)
  aa_per_chrom <- read.table(paste0("C:/Users/17735/Downloads/Dmel_Duplicate_Genes/Augustus_Output/aa_per_chrom/",name,"_aa_per_chrom.txt"), quote="\"", comment.char="")
  aa_per_chrom[aa_per_chrom$V1 == 0,] <- NA
  aa_per_chrom <- na.omit(aa_per_chrom)
  
  chrom <- c('2L','2R','3L','3R','X')
  
  aa_per_chrom$chrom <- chrom
  aa_per_chrom$V1 <- cumsum(aa_per_chrom$V1)
  
  aa_per_chrom$gene_group <- paste0('g',aa_per_chrom$V1)
  
  
  # import annotation output from augustus
  annotations <- read.delim(paste0("C:/Users/17735/Downloads/Dmel_Duplicate_Genes/Augustus_Output/",name,"_annotations.gff"), header=FALSE, comment.char="#")
  
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
  prot_fasta <- format_protein_fasta(annotations)
  
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


####################################################################################################################

# run blastp with each flys proteins against itself
for (row in 1:nrow(fly_name_list)){
  name <- fly_name_list[row,1]
  
  file <- paste('./Prot_Fastas/',name,'_prot.fasta',sep='')

  metablastr::blast_protein_to_protein(file,file, search_type = 'protein_to_protein', 
                           output.path = "./Prot_Blast_Output/")
  
}








