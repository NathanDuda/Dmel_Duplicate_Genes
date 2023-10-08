# Title: 2_extract_genes_from_genomes.R
# Author: Nathan Duda
# Date: 10/6/2023
# Purpose: Extract the nucleotide sequences of each gene annotation from full chromosome sequences. 


source("startup.R")

# input chromosome sequences
genomes <- readDNAStringSet("./Genome_Assemblies/All_Genome_Assemblies.fasta")

# read in gff3 file with all annotations 
all_genes <- read.delim("./Genes/All_Genes.gff3", header=FALSE, comment.char="#")

# remove unnecessary characters 
all_genes$V2 <- gsub('.gff3','',all_genes$V2)
all_genes$V2 <- gsub('Liftoff ','',all_genes$V2)
all_genes <- all_genes[,c(1:5,7,9)]
colnames(all_genes) <- c('chrom','fly','type','start','end','strand','id')

# keep only types: gene, mRNA, exon, CDS, and keep only chromosomes: 2L, 2R, 3L, 3R, X, and remove flybase fly
all_genes <- all_genes %>% 
  filter(!type %in% c('ncRNA','snoRNA','pre_miRNA','miRNA','snRNA','tRNA','pseudogene','rRNA','transcript')) %>%
  filter(chrom %in% c('2L','2R','3L','3R','X')) %>%
  filter(!fly %in% c('FlyBase'))

# extract gene sequences 
all_genes <- all_genes %>%
  mutate(nuc_sequence = as.character(subseq(genomes[paste0(' ',fly,'.fasta_',chrom, sep = '')],
                                            start = start, end = end)))

# separate by types and write to file 
mrna <- all_genes[all_genes$type == 'mRNA',]
write.table(mrna,'all_mRNA.tsv')

exon <- all_genes[all_genes$type == 'exon',]
write.table(exon,'all_exon.tsv')

CDS <- all_genes[all_genes$type == 'CDS',]
write.table(CDS,'all_CDS.tsv')

gene <- all_genes[all_genes$type == 'gene',]
write.table(gene,'all_gene.tsv')

