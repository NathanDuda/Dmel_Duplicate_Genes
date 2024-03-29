
# Author: Nathan Duda
# Purpose: Format the AUGUSTUS annotation output into a tidy dataframe 
source("startup.R")


# function to extract sequences using coordinates 
extract_sequence <- function(chrom, start, end, strand, genome) {
  matches <- genome[names(genome) == chrom]
  if (length(matches) > 0) {
    seq <- subseq(matches, start, end)
      if (strand == "-") {seq <- reverseComplement(seq)}
    as.character(seq)
  } 
  else {NA}
}


# read in list with all fly names 
fly_name_list <- read.table("./Individuals_List_47.txt", quote="\"", comment.char="")

# create empty dataframe for output
all_annotations <- as.data.frame(matrix(ncol=13,nrow=0))

# loop over all flies 
for (row in 1:nrow(fly_name_list)){
  
  # get name of current fly 
  name <- fly_name_list[row,1]
  
  # get genome of fly
  genome <- readDNAStringSet(paste0("./Genome_Assemblies/", name, ".fasta"))
  
  # import the dataframe 
  aa <- read.csv(paste0("./Prot_Fastas/",name,"_annotations.fasta"), sep="", header = F)
  
  # format the fasta into a dataframe 
  colnames(aa)[1] <- 'aa'
  aa <- aa %>%
    mutate(Group = cumsum(grepl("^>", aa))) %>%
    group_by(Group) %>%
    mutate(prot = paste(aa, collapse = "")) %>%
    ungroup() %>%
    select(-Group)
  aa <- aa[grepl(">", aa$aa), ]
  aa$prot <- sub(".*\\.t1", "", aa$prot)
  
  
  # get protein length
  aa$nchar <- nchar(aa$prot)
  
  # import the number of genes on each chromosome to get the chromosomes of each gene (since theyre in order)
  aa_per_chrom <- read.table(paste0("./Augustus_Output/aa_per_chrom/",name,"_aa_per_chrom.txt"), quote="\"", comment.char="")
  aa_per_chrom[aa_per_chrom$V1 == 0,] <- NA
  aa_per_chrom <- na.omit(aa_per_chrom)
  
  aa_per_chrom$chrom <- c('2L','2R','3L','3R','X')
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
  
  # get the nucleotide sequences for each gene
  names(genome) <- sub(paste0(' ', name, '.fasta_'),'',names(genome))
  
  nuc_seqs <- annotations %>%
    filter(type=='CDS')
    
  nuc_seqs$nuc <- mapply(
    extract_sequence,
    nuc_seqs$chrom,
    nuc_seqs$start,
    nuc_seqs$end,
    nuc_seqs$strand,
    MoreArgs = list(genome)
  )
  
  
  nuc_seqs_neg <- nuc_seqs %>%
    filter(strand == '-') %>%
    select(gene_group, chrom, nuc, strand, start) %>%
    group_by(gene_group) %>%
    arrange(desc(start)) %>%
    mutate(nuc = paste(nuc, collapse = "")) %>%
    distinct(gene_group, chrom, nuc, strand)
  
  nuc_seqs_pos <- nuc_seqs %>%
    filter(strand == '+') %>%
    select(gene_group, chrom, nuc, strand, start) %>%
    group_by(gene_group) %>%
    arrange(start) %>%
    mutate(nuc = paste(nuc, collapse = "")) %>%
    distinct(gene_group, chrom, nuc, strand)
  
  nuc_seqs <- rbind(nuc_seqs_pos,nuc_seqs_neg)
  
  
  # keep only the gene in the annotations
  annotations <- annotations[annotations$type == 'gene',]
  
  # merge nuc sequences with the gene annotations
  nuc_seqs <- nuc_seqs[c('gene_group','nuc')]
  annotations <- merge(annotations,nuc_seqs,by='gene_group')
  
  # remove stop codons
  annotations$nuc <- substr(annotations$nuc, 1, nchar(annotations$nuc) - 3)
  
  
  # get the proteins for each gene 
  aa$gene_group <- paste0('g',1:nrow(aa))
  annotations <- merge(annotations,aa,by='gene_group')
  annotations <- annotations %>% select(-aa)
  
  # add fly name to the annotation dataframe
  annotations$fly <- name 
  
  # add annotation to dataframe with all annotations
  all_annotations <- rbind(all_annotations,annotations)
  
  print(row)
}

# remove the gene name from the protein sequence column 
all_annotations$prot <- sub(">g\\d+", "", all_annotations$prot)

# remove the proteins that dont start with M
all_annotations <- all_annotations[grepl("^M", all_annotations$prot), ]

# make gn id contain fly name since same id is different gene in different fly
all_annotations$gene_group <- paste0(all_annotations$fly,'_',all_annotations$gene_group)

# write to file
write.table(all_annotations, file = 'Annotations_Gene.tsv')
