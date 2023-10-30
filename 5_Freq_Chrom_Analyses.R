# Title: 5_Chrom_Analyses.R
# Author: Nathan Duda
# Date: 10/29/2023
# Purpose: Make a frequency distribution of the duplicates in the 52 flies.
#          Analyze the chromosomal locations of the duplicates.


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
all_annotations <- read.csv("C:/Users/17735/Downloads/Dmel_Duplicate_Genes/Annotations.tsv", sep="")
# keep only the full gene information 
all_annotations <- all_annotations[all_annotations$type=='gene',]

# get the number of genes on each chromosome for each fly 
chrom_gene_counts <- as.data.frame(table(all_annotations$chrom,all_annotations$fly))
colnames(chrom_gene_counts) <- c('chrom','fly','n_genes')

# merge number of genes with chromosome lengths 
chrom_info <- merge(chrom_lengths,chrom_gene_counts,by=c('chrom','fly'))

# plot the number of genes on a chromosome versus the chromosomes length
gg <- ggplot(chrom_info,aes(x=length,y=n_genes,color=chrom)) +
  geom_point() +
  theme_bw() +
  xlab('Length') +
  ylab('Number of genes') +
  geom_segment(y=mean(chrom_info[chrom_info$chrom=='2L',]$n_genes), 
                   yend=mean(chrom_info[chrom_info$chrom=='2L',]$n_genes), 
                   x=min(chrom_info[chrom_info$chrom=='2L',]$length), 
                   xend=max(chrom_info[chrom_info$chrom=='2L',]$length), color = "#F8766D") +
  geom_segment(y=mean(chrom_info[chrom_info$chrom=='2R',]$n_genes), 
                   yend=mean(chrom_info[chrom_info$chrom=='2R',]$n_genes), 
                   x=min(chrom_info[chrom_info$chrom=='2R',]$length), 
                   xend=max(chrom_info[chrom_info$chrom=='2R',]$length), color = "#A3A500") +
  geom_segment(y=mean(chrom_info[chrom_info$chrom=='3L',]$n_genes), 
                   yend=mean(chrom_info[chrom_info$chrom=='3L',]$n_genes), 
                   x=min(chrom_info[chrom_info$chrom=='3L',]$length), 
                   xend=max(chrom_info[chrom_info$chrom=='3L',]$length), color = "#00BF7D") +
  geom_segment(y=mean(chrom_info[chrom_info$chrom=='3R',]$n_genes), 
                   yend=mean(chrom_info[chrom_info$chrom=='3R',]$n_genes), 
                   x=min(chrom_info[chrom_info$chrom=='3R',]$length), 
                   xend=max(chrom_info[chrom_info$chrom=='3R',]$length), color = "#00B0F6") +
  geom_segment(y=mean(chrom_info[chrom_info$chrom=='X',]$n_genes), 
                   yend=mean(chrom_info[chrom_info$chrom=='X',]$n_genes), 
                   x=min(chrom_info[chrom_info$chrom=='X',]$length), 
                   xend=max(chrom_info[chrom_info$chrom=='X',]$length), color = "#E76BF3") +
    
    
  geom_segment(y=max(chrom_info[chrom_info$chrom=='2L',]$n_genes), 
                   yend=min(chrom_info[chrom_info$chrom=='2L',]$n_genes), 
                   x=mean(chrom_info[chrom_info$chrom=='2L',]$length), 
                   xend=mean(chrom_info[chrom_info$chrom=='2L',]$length), color = "#F8766D") +
  geom_segment(y=max(chrom_info[chrom_info$chrom=='2R',]$n_genes), 
                   yend=min(chrom_info[chrom_info$chrom=='2R',]$n_genes), 
                   x=mean(chrom_info[chrom_info$chrom=='2R',]$length), 
                   xend=mean(chrom_info[chrom_info$chrom=='2R',]$length), color = "#A3A500") +
  geom_segment(y=max(chrom_info[chrom_info$chrom=='3L',]$n_genes), 
                   yend=min(chrom_info[chrom_info$chrom=='3L',]$n_genes), 
                   x=mean(chrom_info[chrom_info$chrom=='3L',]$length), 
                   xend=mean(chrom_info[chrom_info$chrom=='3L',]$length), color = "#00BF7D") +
  geom_segment(y=max(chrom_info[chrom_info$chrom=='3R',]$n_genes), 
                   yend=min(chrom_info[chrom_info$chrom=='3R',]$n_genes), 
                   x=mean(chrom_info[chrom_info$chrom=='3R',]$length), 
                   xend=mean(chrom_info[chrom_info$chrom=='3R',]$length), color = "#00B0F6") +
  geom_segment(y=max(chrom_info[chrom_info$chrom=='X',]$n_genes), 
                   yend=min(chrom_info[chrom_info$chrom=='X',]$n_genes), 
                   x=mean(chrom_info[chrom_info$chrom=='X',]$length), 
                   xend=mean(chrom_info[chrom_info$chrom=='X',]$length), color = "#E76BF3") 
    
ggsave("./Plots/length_ngenes_scatter.jpg", plot = gg, width = 8, height = 6)


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

gg <- ggplot(one_sided_chrom_matrix,aes(x=chrom1,y=chrom2,fill=n)) +
  geom_tile() +
  theme_bw() +
  xlab('') +
  ylab('') +
  labs(fill = "Duplicates")

ggsave("./Plots/chrom_dup_heatmap.jpg", plot = gg, width = 7, height = 5)


# frequency plot 
dup1 <- dups[,c(3,2,21)]
dup2 <- dups[,c(1,2,32)]

colnames(dup1) <- c('dup','fly','chrom')
colnames(dup2) <- c('dup','fly','chrom')

dup <- rbind(dup1,dup2)


plot_freq_dist <- dup[c('dup','fly')]
plot_freq_dist <- unique(plot_freq_dist)

plot_freq_dist <- as.data.frame(table(plot_freq_dist$dup))
plot_freq_dist <- as.data.frame(table(plot_freq_dist$Freq))

gg <- ggplot(plot_freq_dist,aes(x=Var1,y=Freq)) +
  geom_bar(stat='identity') +
  scale_x_discrete(limits=factor(1:52)) +
  xlab('Number of flies duplicate is present in') +
  ylab('Freq') +
  theme_bw()

ggsave("./Plots/freq_distribution_barplot.jpg", plot = gg, width = 10, height = 3)


# distance when on same chrom
same_chrom <- dups[dups$dup1_chrom == dups$dup2_chrom,]
same_chrom$dist <- same_chrom$dup1_end_on_chrom - same_chrom$dup2_start_on_chrom
same_chrom$dist2 <- same_chrom$dup2_end_on_chrom - same_chrom$dup1_start_on_chrom

same_chrom$dist[same_chrom$dist < 0] <- NA
same_chrom$dist2[same_chrom$dist2 < 0] <- NA

same_chrom$dist <- coalesce(same_chrom$dist,same_chrom$dist2)
same_chrom <- same_chrom %>% select(-dist2)

gg <- ggplot(same_chrom,aes(y=dist,x=dup1_chrom,color=dup1_chrom)) +
  geom_boxplot() +
  xlab('Distance between duplicates') +
  ylab('') +
  guides(color = guide_legend(title = "Chromosome")) +
  theme_bw() +
  geom_hline(yintercept = mean(chrom_lengths[chrom_lengths$chrom=='2L',]$length), linewidth = 1.1, color = "#F8766D") +
  geom_hline(yintercept = mean(chrom_lengths[chrom_lengths$chrom=='2R',]$length), linewidth = 1.1, color = "#A3A500") +
  geom_hline(yintercept = mean(chrom_lengths[chrom_lengths$chrom=='3L',]$length), linewidth = 1.1, color = "#00BF7D") +
  geom_hline(yintercept = mean(chrom_lengths[chrom_lengths$chrom=='3R',]$length), linewidth = 1.1, color = "#00B0F6") +
  geom_hline(yintercept = mean(chrom_lengths[chrom_lengths$chrom=='X',]$length), linewidth = 1.1, color = "#E76BF3")

ggsave("./Plots/chrom_dist_between_dups_boxplot.jpg", plot = gg, width = 7, height = 5)


# number of duplicates per chromosome in each fly 
n_dups <- as.data.frame(table(dup))

n_dups <- n_dups %>%
  group_by(fly,chrom) %>%
  summarize(n = sum(Freq), .groups = "drop")

add <- n_dups %>%
  group_by(fly) %>%
  summarize(n = sum(n))

add$chrom <- 'Total'
n_dups <- bind_rows(n_dups, add)

gg <- ggplot(n_dups,aes(x=factor(chrom,levels=c('2L','2R','3L','3R','X','Total')),y=fly,fill=n,label=scales::comma(n))) +
  geom_tile() +
  geom_text(color = "white", size=2, hjust=1) +  
  theme_bw() +
  xlab('Chromosome') +
  ylab('Fly') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave("./Plots/dup_chrom_fly_heatmap.jpg", plot = gg, width = 6, height = 6)


