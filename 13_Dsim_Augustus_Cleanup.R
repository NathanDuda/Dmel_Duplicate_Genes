# Title: 10_Dsim_Annotation_Cleanup.R
# Author: Nathan Duda
# Date: 11/10/2023
# Purpose: Format all annotations into a table with gene id, prot sequence, chromosome, location, and other information
#          for the Drosophila simulans dataset. 


source("startup.R")

# import the annotations dataframe 
aa <- read.csv("./Dsim/Prot_Fastas/Dsim_annotations.fasta", sep="", header = F)

# format the fasta into a dataframe 
colnames(aa)[1] <- 'aa'
aa <- aa %>%
  mutate(Group = cumsum(grepl("^>", aa))) %>%
  group_by(Group) %>%
  mutate(prot = paste(aa, collapse = "")) %>%
  ungroup() %>%
  select(-Group)
aa <- aa[grepl(">", aa$aa), ]

# remove the gene name from the protein sequence column 
aa$prot <- sub(">g\\d+", "", aa$prot)

# remove the proteins that dont start with M
aa <- aa[grepl("^M", aa$prot), ]

# get protein length
aa$nchar <- nchar(aa$prot)

# import the number of genes on each chromosome to get the chromosomes of each gene (since theyre in order)
aa_per_chrom <- read.table("./Dsim/Output_Annotations/Dsim_aa_per_chrom.txt", quote="\"", comment.char="")
aa_per_chrom[aa_per_chrom$V1 == 0,] <- NA
aa_per_chrom <- na.omit(aa_per_chrom)

chrom <- c('2L','2R','3L','3R','4','X','unk')

# sum up rest of rows to be unplaced genes 
row_7 <- aa_per_chrom %>% dplyr::slice(7:155) %>% summarise(across(everything(), sum))
aa_per_chrom$V1[7] <- row_7

aa_per_chrom <- as.data.frame(aa_per_chrom[c(1:7),])
aa_per_chrom <- as.data.frame(t(aa_per_chrom))
colnames(aa_per_chrom) <- 'V1'

aa_per_chrom$chrom <- chrom
aa_per_chrom$V1 <- cumsum(aa_per_chrom$V1)

aa_per_chrom$gene_group <- paste0('g',aa_per_chrom$V1)

# import annotation output from augustus
annotations <- read.delim("./Dsim/Output_Annotations/Dsim_annotations.gff", header=FALSE, comment.char="#")

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
annotations <- annotations 

# add fly name to the annotation dataframe
annotations$fly <- 'Dsim' 

# write to table
write.table(annotations,'./Dsim/Dsim_Annotations.tsv')