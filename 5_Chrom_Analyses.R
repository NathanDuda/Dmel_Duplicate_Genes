


source("startup.R")

# import my duplicate genes 
dups <- read.csv("C:/Users/17735/Downloads/Dmel_Duplicate_Genes/Duplicate_Proteins.tsv", sep="")


# import chromosomes for every fly
genomes <- readDNAStringSet("./Genome_Assemblies/All_Genome_Assemblies.fasta")

# get chromosome lengths
chrom_lengths <- data.frame(chrom = names(genomes), length = width(genomes))
chrom_lengths <- chrom_lengths %>% 
  separate(chrom, into = c("fly", "chrom"), sep = ".fasta_")

# remove the whitespace in the fly column 
chrom_lengths$fly <- trimws(chrom_lengths$fly)


# import annotations for all flies 
###all_annotations <- read.csv("C:/Users/17735/Downloads/Dmel_Duplicate_Genes/Annotations.tsv", sep="")
# keep only the full gene information 
all_annotations <- all_annotations[all_annotations$type=='gene',]

# get the number of genes on each chromosome for each fly 
chrom_gene_counts <- as.data.frame(table(all_annotations$chrom,all_annotations$fly))
colnames(chrom_gene_counts) <- c('chrom','fly','n_genes')

# merge number of genes with chromosome lengths 
chrom_info <- merge(chrom_lengths,chrom_gene_counts,by=c('chrom','fly'))

# plot the number of genes on a chromosome versus the chromosomes length
ggplot(chrom_info,aes(x=length,y=n_genes,color=chrom)) +
  geom_point() +
  theme_bw()


# chromosomal location heatmap 
chrom_matrix <- as.data.frame(table(dups$dup1_chrom,dups$dup2_chrom))

chrom_matrix <- as.data.frame(pivot_wider(chrom_matrix,names_from = 'Var2', values_from = 'Freq'))
rownames(chrom_matrix) <- chrom_matrix$Var1
chrom_matrix <- chrom_matrix[,-c(1)]

one_sided_chrom_matrix <- chrom_matrix + t(chrom_matrix)

one_sided_chrom_matrix <- one_sided_chrom_matrix
one_sided_chrom_matrix$chrom1 <- rownames(one_sided_chrom_matrix)
one_sided_chrom_matrix <- pivot_longer(one_sided_chrom_matrix,names_to = 'chrom2',values_to = 'n',cols = c(1:5))

one_sided_chrom_matrix <- one_sided_chrom_matrix[!duplicated(one_sided_chrom_matrix$n),]

ggplot(one_sided_chrom_matrix,aes(x=chrom1,y=chrom2,fill=n)) +
  geom_tile() +
  theme_bw()


# frequency plot 
plot_freq_dist <- dups[c('dup1','dup2')]
#plot_freq_dist <- unique(plot_freq_dist)

plot_freq_dist <- as.data.frame(table(plot_freq_dist$dup1))
plot_freq_dist <- as.data.frame(table(plot_freq_dist$Freq))

library(ggplot2)
ggplot(plot_freq_dist,aes(x=Var1,y=Freq)) +
  geom_bar(stat='identity') +
  scale_x_discrete(limits=c(1:52)) +
  xlab('number of flies duplicate is present in') +
  ylab('freq') +
  theme_bw()


# distance when on same chrom
same_chrom <- dups[dups$dup1_chrom == dups$dup2_chrom,]
same_chrom$dist <- same_chrom$dup1_end_on_chrom - same_chrom$dup2_start_on_chrom
same_chrom$dist2 <- same_chrom$dup2_end_on_chrom - same_chrom$dup1_start_on_chrom

same_chrom$dist[same_chrom$dist < 0] <- NA
same_chrom$dist2[same_chrom$dist2 < 0] <- NA

same_chrom$dist <- coalesce(same_chrom$dist,same_chrom$dist2)
same_chrom <- same_chrom %>% select(-dist2)

ggplot(same_chrom,aes(y=dist,x=dup1_chrom,color=dup1_chrom)) +
  geom_boxplot() +
  theme_bw()


# number of duplicates per chromosome in each fly 
dup1 <- dups[,c(3,2,21)]
dup2 <- dups[,c(1,2,32)]

colnames(dup1) <- c('dup','fly','chrom')
colnames(dup2) <- c('dup','fly','chrom')

dup <- rbind(dup1,dup2)

dup <- as.data.frame(table(dup))

d <- dup %>%
  group_by(fly,chrom) %>%
  summarize(n = sum(Freq), .groups = "drop")

ggplot(d,aes(x=fly,y=chrom,fill=n)) +
  geom_tile() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))




