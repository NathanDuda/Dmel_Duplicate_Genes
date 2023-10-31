

source("startup.R")


# import flybase proteins 
fbgn_prots <- read.csv("./dmel-all-translation-r6.54.fasta", sep="", header = F)

# format the fasta into a dataframe 
colnames(fbgn_prots)[1] <- 'prot_seq'
fbgn_prots <- fbgn_prots %>%
  mutate(Group = cumsum(grepl("^>", prot_seq))) %>%
  group_by(Group) %>%
  mutate(prot = paste(prot_seq, collapse = "")) %>%
  ungroup() %>%
  mutate(prot_seq = prot) %>%
  select(-Group,-prot)

fbgn_prots$FBpp <- substr(fbgn_prots$prot_seq, start = 2, stop = 12)
fbgn_prots$prot_seq <- substring(fbgn_prots$prot_seq, first = 13)


# remove the unnecessary rows 
fbgn_prots <- fbgn_prots %>% mutate_all(~na_if(., "")) %>% na.omit()

# extract the FBgn id 
fbgn_prots$FBgn <- substr(fbgn_prots$V6, start = 8, stop = 18)

# format the dataframe
fbgn_prots <- fbgn_prots[,c(12,13,1,3,5:9)]

# get only the FBgn and FBpp 
fbgn_fbpp <- fbgn_prots[c('FBgn','FBpp')]
rm(fbgn_prots)


### get FBgn for each gn from FBgn blast 

# read in list with all fly names 
fly_name_list <- read.table("./Individuals_List.txt", quote="\"", comment.char="")

# create empty dataframes for output
my_FBgns <- as.data.frame(matrix(ncol=4,nrow=0))

# loop over all flies 
for (row in 1:nrow(fly_name_list)){
  
  # get name of current fly 
  name <- fly_name_list[row,1]
  
  # import blast result
  blastp_FBgn <- read.delim(paste0("C:/Users/17735/Downloads/Dmel_Duplicate_Genes/FBgn_Blast_Output/",name,"_blastp.tsv"), header=FALSE)
  colnames(blastp_FBgn) <- c('qseqid', 'sseqid', 'length', 'bitscore', 'pident', 'evalue')
  
  # keep only the gene with its FBgn and the evalue
  blastp_FBgn <- blastp_FBgn[,c('qseqid','sseqid','evalue')]
  colnames(blastp_FBgn) <- c('gn','FBpp','evalue')
  
  # keep only highest evalue blast
  blastp_FBgn <- blastp_FBgn %>%
    group_by(gn) %>%
    filter(evalue == min(evalue)) %>%
    ungroup() %>%
    distinct(gn,evalue, .keep_all = T) # keep only one if lowest evalue is same for multiple 

  # combine all genes with their FBgns into dataframe 
  blastp_FBgn$fly <- name
  my_FBgns <- rbind(my_FBgns,blastp_FBgn)
  
  print(row)
  
}

# get FBgn for each pp
my_FBgns <- merge(my_FBgns,fbgn_fbpp,by='FBpp')

# extract the FBgns 
FBgn <- my_FBgns[c('FBgn')]
FBgn <- unique(FBgn)

# write a list of FBgns 
write.table(FBgn, file = './FBgn_list.txt', sep = "\t", row.names = F, col.names = F, quote = F)







###########################################################################################
#

any(duplicated(my_FBgns[c('gn','fly','FBgn')]))

t <- as.data.frame(table(my_FBgns$fly))
any(duplicated(my_FBgns[c('fly','FBgn')])) # T
any(duplicated(my_FBgns[c('fly','gn')]))   # F
# so multiple of my gn's hits to the same FBgn occasionally 

t <- my_FBgns[(duplicated(my_FBgns[c('fly','FBgn')])),]


n_of_my_gns_hit_to_FBgn <- as.data.frame(table(t$fly,t$FBgn))
n_of_my_gns_hit_to_FBgn[n_of_my_gns_hit_to_FBgn$Freq==0,] <- NA
n_of_my_gns_hit_to_FBgn <- na.omit(n_of_my_gns_hit_to_FBgn)

length(unique(my_FBgns$FBgn))


# import my duplicate genes 
dups <- read.csv("./Duplicate_Proteins.tsv", sep="")





# merge FBgns with dups dataframe

dup_FBgns <- dup_FBgns_orig

# remove duplicates caused by multiple FBpp per FBgn
dup_FBgns <- dup_FBgns[!duplicated(dup_FBgns[c('gn','fly','FBgn')]),]
dup_FBgns <- dup_FBgns[c('gn','fly','FBgn')]

colnames(dup_FBgns) <- c('dup1','fly','dup1_FBgn')
t <- merge(dups,dup_FBgns,by=c('dup1','fly'))

colnames(dup_FBgns) <- c('dup2','fly','dup1_FBgn')
t <- merge(dups,dup_FBgns,by=c('dup2','fly'))




flybase <- read_delim("FlyBase_Fields_download.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)



#############
# get FBgn for the duplicates FBpps
t <- merge(my_FBgns,fbgn_fbpp,by='FBpp')


for (row in 1:nrow(fly_name_list)){
  
  # get name of current fly 
  name <- fly_name_list[row,1]
  
  lenn <- as.character(length(unique(t[t$fly == name,]$FBgn)))
  print(paste0(lenn,'   ',name))
  
}

# ~13,400 for each fly 


##

###################



all_genes <- read.csv("./Annotations_Gene.tsv", sep="")

t <- as.data.frame(table(all_genes$fly))

#
drosomics_genes <- read.delim("./Drosomics_Annotated_Genes.gff3", header=FALSE, comment.char="#")

# remove unnecessary characters 
drosomics_genes$V2 <- gsub('.gff3','',drosomics_genes$V2)
drosomics_genes$V2 <- gsub('Liftoff ','',drosomics_genes$V2)
drosomics_genes <- drosomics_genes[,c(1:5,7,9)]
colnames(drosomics_genes) <- c('chrom','fly','type','start','end','strand','id')


fb_drosomics <- drosomics_genes[drosomics_genes$fly == 'FlyBase',]
# keep only types: gene, mRNA, exon, CDS, and keep only chromosomes: 2L, 2R, 3L, 3R, X, and remove flybase fly
drosomics_genes <- drosomics_genes %>% 
  filter(type == 'gene') %>%
  filter(chrom %in% c('2L','2R','3L','3R','X')) %>%
  filter(!fly %in% c('FlyBase'))
#

perfect_my_drosomics_genes <- merge(all_genes,drosomics_genes,
                                    by=c('chrom','fly','start','end','strand'))

#############





all_annotations <- read.csv("./Annotations.tsv", sep="")
cds <- all_annotations[all_annotations$type == 'start_codon',]
t <- as.data.frame(table(cds$fly))


all_genes <- read.csv("./Annotations_Gene.tsv", sep="")
length(unique(all_genes$gene_group))
t <- as.data.frame(table(all_genes$fly))

# import chromosomes for every fly
genomes <- readDNAStringSet("./Genome_Assemblies/All_Genome_Assemblies.fasta")

# get chromosome lengths
chrom_lengths <- data.frame(chrom = names(genomes), length = width(genomes))
chrom_lengths <- chrom_lengths %>% 
  separate(chrom, into = c("fly", "chrom"), sep = ".fasta_")

genome_lengths <- chrom_lengths %>%
  group_by(fly) %>%
  summarise(total_length = sum(length))


############################################################


