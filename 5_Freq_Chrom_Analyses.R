# Title: 5_Chrom_Analyses.R
# Author: Nathan Duda
# Date: 10/29/2023
# Purpose: Make a frequency distribution of the duplicates in the 52 flies.
#          Analyze the chromosomal locations of the duplicates.


source("startup.R")

# import my duplicate genes 
dups <- read.csv("./Duplicate_Proteins.tsv", sep="")

# remove the crappy assemblies
dups <- dups[!dups$fly %in% c('I23','ZH26','N25','T29A','B59'),]

# import chromosomes for every fly
genomes <- readDNAStringSet("./Genome_Assemblies/All_Genome_Assemblies.fasta")

# get chromosome lengths
chrom_lengths <- data.frame(chrom = names(genomes), length = width(genomes))
chrom_lengths <- chrom_lengths %>% 
  separate(chrom, into = c("fly", "chrom"), sep = ".fasta_")

# remove the whitespace in the fly column 
chrom_lengths$fly <- trimws(chrom_lengths$fly)


# import annotations for all flies 
all_annotations <- read.csv("./Annotations_Gene.tsv", sep="")

# remove the crappy assemblies
all_annotations <- all_annotations[!all_annotations$fly %in% c('I23','ZH26','N25','T29A','B59'),]

# get the number of genes on each chromosome for each fly 
chrom_gene_counts <- as.data.frame(table(all_annotations$chrom,all_annotations$fly))
colnames(chrom_gene_counts) <- c('chrom','fly','n_genes')

# merge number of genes with chromosome lengths 
chrom_info <- merge(chrom_lengths,chrom_gene_counts,by=c('chrom','fly'))

# plot the number of genes on a chromosome versus the chromosomes length

chrom_info <- chrom_info[!chrom_info$fly %in% c('I23','ZH26','N25','T29A','B59'),]

gg <- ggplot(chrom_info,aes(x=length,y=n_genes,color=chrom)) +
  geom_point() +
  theme_bw() +
  xlab('Length') +
  ylab('Number of genes') +
  guides(color = guide_legend(title = "Chromosome")) +
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

########################################
longest_chroms <- chrom_info %>%
  group_by(chrom) %>%
  filter(length == max(length))


outlier_threshold = 2.5

#chrom_info <- chrom_info_orig

chrom_info <- chrom_info %>%
  group_by(chrom) %>%
  mutate(length_outlier = ifelse(abs(length - median(length)) > outlier_threshold * IQR(length), "Outlier", "Not Outlier")) %>%
  mutate(n_genes_outlier = ifelse(abs(n_genes - median(n_genes)) > outlier_threshold * IQR(n_genes), "Outlier", "Not Outlier")) %>%
  mutate(outlier = case_when(n_genes_outlier == 'Outlier' | length_outlier == 'Outlier' ~ 'Outlier',
                             n_genes_outlier == 'Not Outlier' & length_outlier == 'Not Outlier' ~ 'Not Outlier'))

unique(chrom_info[chrom_info$outlier == 'Outlier',]$fly)

chrom_info <- chrom_info[!chrom_info$fly %in% (chrom_info[chrom_info$outlier == 'Outlier',]$fly),]



ggplot(chrom_info,aes(x=length,y=n_genes,color=chrom)) + 
  geom_point() +
  theme_bw() +
  xlab('Length') +
  ylab('Number of genes') +
  guides(color = guide_legend(title = "Chromosome")) +
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

########################################


# chromosomal location heatmap 
chrom_matrix <- as.data.frame(table(dups$dup_chrom,dups$dup_family,dups$fly))

colnames(chrom_matrix) <- c('chrom','dup_family','fly','Freq')


# number of unique duplicate families 
dup_fams <- chrom_matrix %>%
  mutate(Freq = case_when(Freq == 0 ~ NA,
         T ~ Freq)) %>%
  na.omit() %>%
  group_by(chrom,fly) %>%
  summarize(Freq=length(unique(dup_family)), .groups = 'keep')

gg <- ggplot(dup_fams,aes(x=fly,y=chrom,fill=Freq)) +
  geom_tile() +
  theme_bw() + 
  xlab('') +
  ylab('') +
  labs(fill = "Dup Families")
  
ggsave("./Plots/trash.jpg", plot = gg, width = 10, height = 8)


# total number of copies 
copies <- chrom_matrix %>%
  group_by(chrom,fly) %>%
  summarize(Freq=sum(Freq), .groups = 'keep')

gg <- ggplot(copies,aes(x=fly,y=chrom,fill=Freq)) +
  geom_tile() +
  theme_bw() + 
  xlab('') +
  ylab('') +
  labs(fill = "Copies    a")
  
ggsave("./Plots/trash2.jpg", plot = gg, width = 10, height = 8)


# average copies per family 
copies_per_fam <- copies
copies_per_fam$Freq <- copies_per_fam$Freq / dup_fams$Freq

gg <- ggplot(copies_per_fam,aes(x=fly,y=chrom,fill=Freq)) +
  geom_tile() +
  theme_bw() + 
  xlab('') +
  ylab('') +
  labs(fill = "Avg copies per fam")

ggsave("./Plots/trash3.jpg", plot = gg, width = 10, height = 8)





#





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
  xlab('Chromosome of duplicate 2') +
  ylab('Chromosome of duplicate 1') +
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
  #scale_x_discrete(limits=factor(1:52)) + # not 52 without crappy assemblies 
  xlab('Number of flies duplicate is present in') +
  ylab('Freq') +
  theme_bw()

ggsave("./Plots/freq_distribution_barplot.jpg", plot = gg, width = 10, height = 3)


# distance when on same chrom
dups <- dups %>%
  mutate(same_chrom = case_when(dup1_chrom == dup2_chrom ~ 'same_chrom',
                                dup1_chrom != dup2_chrom ~ 'diff_chrom')) %>%
  mutate(dist = case_when((same_chrom=='same_chrom') & (dup1_start_on_chrom > dup2_start_on_chrom) ~ (dup1_start_on_chrom - dup2_end_on_chrom),
                          (same_chrom=='same_chrom') & (dup2_start_on_chrom > dup1_start_on_chrom) ~ (dup2_start_on_chrom - dup1_end_on_chrom)))
write.table(dups,'Duplicate_Proteins_2.tsv')


same_chrom <- dups[dups$same_chrom == 'same_chrom',]

gg <- ggplot(same_chrom,aes(y=dist,x=dup1_chrom,color=dup1_chrom)) +
  geom_boxplot() +
  ylab('Distance between duplicates') +
  xlab('') +
  guides(color = guide_legend(title = "Chromosome")) +
  theme_bw() +
  ggsignif::geom_signif(comparisons = list(c('X','2L'),c('X','2R'),c('X','3L'),c('X','3R')),
                        map_signif_level = c('****'=0.0001,"***"=0.001, "**"=0.01, "*"=0.05),
                        textsize = 6, 
                        vjust = 0.5,
                        y_position = c(3.5e+07,3.75e+07,4e+07,4.25e+07),
                        test = 't.test') +
  geom_hline(yintercept = mean(chrom_lengths[chrom_lengths$chrom=='2L',]$length), linewidth = 1.1, color = "#F8766D") +
  geom_hline(yintercept = mean(chrom_lengths[chrom_lengths$chrom=='2R',]$length), linewidth = 1.1, color = "#A3A500") +
  geom_hline(yintercept = mean(chrom_lengths[chrom_lengths$chrom=='3L',]$length), linewidth = 1.1, color = "#00BF7D") +
  geom_hline(yintercept = mean(chrom_lengths[chrom_lengths$chrom=='3R',]$length), linewidth = 1.1, color = "#00B0F6") +
  geom_hline(yintercept = mean(chrom_lengths[chrom_lengths$chrom=='X',]$length), linewidth = 1.1, color = "#E76BF3")

ggsave("./Plots/chrom_dist_between_dups_boxplot.jpg", plot = gg, width = 7, height = 8)


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
  labs(fill = "Number of Duplicates") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave("./Plots/dup_chrom_fly_heatmap.jpg", plot = gg, width = 7, height = 6)



