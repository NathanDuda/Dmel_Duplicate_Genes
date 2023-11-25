# Title: 5_Chrom_Analyses.R
# Author: Nathan Duda
# Date: 10/29/2023
# Purpose: Make a frequency distribution of the duplicates in the 52 flies.
#          Analyze the chromosomal locations of the duplicates.


source("startup.R")

# import my duplicate genes 
dups <- read.csv("./Duplicate_Proteins.tsv", sep="")

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

######################################
# get chromosome statistics 
sd <- chrom_info %>%
  group_by(chrom) %>%
  summarize(sd_length = sd(length), mean_length = mean(length),
            sd_n_genes = sd(n_genes), mean_n_genes = mean(n_genes))

# calculate gc content
gc_content <- data.frame(
  Name = names(genomes),
  GC_Content = sapply(genomes, function(seq) {
    gc_count <- sum(strsplit(as.character(seq), "")[[1]] %in% c("G", "C"))
    gc_content <- gc_count / length(seq)
    return(gc_content)
  }), stringsAsFactors = FALSE)


gc_content <- gc_content %>% separate(Name,into = c('fly','chrom'),sep='.fasta_')
gc_content <- gc_content[!gc_content$fly %in% c('I23','ZH26','N25','T29A','B59'),]
gc_content <- gc_content %>%
  group_by(chrom) %>%
  summarize(sd_gc_content = sd(GC_Content), mean_gc_content = mean(GC_Content))

sd <- merge(sd,gc_content,by='chrom')

######################################


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


# chromosomal location heatmap 
chrom_matrix <- as.data.frame(table(dups$dup_chrom,dups$dup_family,dups$fly))

colnames(chrom_matrix) <- c('chrom','dup_family','fly','Freq')


# number of unique duplicate families per chrom
dup_fams_chrom <- chrom_matrix %>%
  mutate(Freq = case_when(Freq == 0 ~ NA,
         T ~ Freq)) %>%
  na.omit() %>%
  group_by(chrom,fly) %>%
  summarize(Freq=length(unique(dup_family)), .groups = 'keep')

# number of unique duplicate families 
dup_fams <- chrom_matrix %>%
  mutate(Freq = case_when(Freq == 0 ~ NA,
                          T ~ Freq)) %>%
  na.omit() %>%
  group_by(chrom,fly) %>%
  summarize(Freq=length(unique(dup_family)), .groups = 'keep') %>%
  ungroup()

  #group_by(fly) %>%
  #summarize(Freq=sum(Freq))


# total number of copies 
copies <- chrom_matrix %>%
  group_by(chrom,fly) %>%
  summarize(Freq=sum(Freq), .groups = 'keep') %>%
  ungroup()

total <- copies %>%
  group_by(fly) %>%
  summarise(Freq = sum(Freq)) %>%
  mutate(chrom='total')

copies <- rbind(total,copies)

gg <- ggplot(copies,aes(y=fly,x=factor(chrom,levels=c('2L','2R','3L','3R','X','total')),fill=Freq)) +
  geom_tile() +
  theme_bw() + 
  xlab('') +
  ylab('') +
  geom_text(aes(label = scales::comma(Freq)), color = "white",size=1.9) +
  labs(fill = "Duplicate Copies") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave("./Plots/copies_per_fly_chrom_heatmap.jpg", plot = gg, width = 5, height = 6)

# average copies per family 
copies_per_fam <- copies

colnames(copies_per_fam) <- c('fly','copies_per_fam','chrom')
colnames(dup_fams_chrom) <- c('chrom','fly','dup_fams_chrom')
copies_per_fam <- merge(copies_per_fam,dup_fams_chrom,by=c('fly','chrom'))
copies_per_fam$Freq <- copies_per_fam$copies_per_fam / copies_per_fam$dup_fams_chrom

gg <- ggplot(copies_per_fam,aes(x=fly,y=chrom,fill=Freq)) +
  geom_tile() +
  theme_bw() + 
  xlab('') +
  ylab('') +
  labs(fill = "Avg copies per fam") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave("./Plots/copies_per_dup_fam_heatmap.jpg", plot = gg, width = 10, height = 3)



#
#chrom_matrix <- chrom_matrix_orig

chrom_matrix$dup_family_fly <- paste0(chrom_matrix$dup_family, '_', chrom_matrix$fly)

chrom_matrix <- chrom_matrix %>% arrange(Freq)


#ggplot(chrom_matrix, aes(fill=chrom, y=Freq, x=fly)) + 
#  geom_bar(position="dodge", stat="identity") 


ggplot(chrom_matrix, aes(fill=chrom, y=Freq, x=fly)) + 
  geom_bar(position="stack", stat="identity") +
  theme_bw()
  


#

chrom_matrix <- chrom_matrix %>%
    pivot_wider(names_from = 'chrom',values_from = 'Freq') %>%
    mutate(zero_count = rowSums(. == 0, na.rm = TRUE)) %>%
    filter(zero_count < 5) %>%
    mutate(same_chrom_fam = case_when(zero_count == 4 ~ 'same_chrom',
                                      zero_count < 4 ~ 'diff_chrom(s)')) %>% 
    pivot_longer(cols = c('2L','2R','3L','3R','X'), names_to = 'chrom') %>%
    select(-zero_count)

# add same_chrom information to dups dataframe
same_chrom <- unique(chrom_matrix[c('same_chrom_fam','dup_family_fly')])
dups$dup_family_fly <- paste0(dups$dup_family, '_', dups$fly)
dups <- merge(dups,same_chrom,by='dup_family_fly')

# plot matrix 
table(chrom_matrix$same_chrom_fam,chrom_matrix$chrom)


t <- chrom_matrix[chrom_matrix$value > 0,]
table(t$same_chrom_fam,t$chrom)

t2 <- t %>%
  group_by(chrom,same_chrom_fam) %>%
  summarise(v = sum(value), .groups='keep')


gg <- ggplot(t2,aes(x=same_chrom_fam,y=chrom,fill=v)) +
  geom_tile() +
  theme_bw() +
  xlab('') +
  ylab('') +
  labs(fill='Copies')

ggsave("./Plots/chrom_dup_heatmap.jpg", plot = gg, width = 7, height = 5)


# frequency plot 
freq_dist <- dups[c('dup','fly')]

freq_dist <- as.data.frame(table(freq_dist$dup))
plot_freq_dist <- as.data.frame(table(freq_dist$Freq))

gg <- ggplot(plot_freq_dist,aes(x=Var1,y=Freq)) +
  geom_bar(stat='identity') +
  #scale_x_discrete(limits=factor(1:52)) + # not 52 without crappy assemblies 
  xlab('Number of flies duplicate is present in') +
  ylab('Freq') +
  theme_bw()

ggsave("./Plots/freq_distribution_barplot.jpg", plot = gg, width = 10, height = 3)


# frequency distribution with chromosomes 
plot_freq_dist_chrom <- dups[c('dup','fly','dup_chrom')]

colnames(freq_dist) <- c('dup','Freq')

plot_freq_dist_chrom <- as.data.frame(table(plot_freq_dist_chrom$dup,plot_freq_dist_chrom$dup_chrom))

plot_freq_dist_chrom <- as.data.frame(table(plot_freq_dist_chrom$Freq,plot_freq_dist_chrom$Var2))

plot_freq_dist_chrom <- plot_freq_dist_chrom[!plot_freq_dist_chrom$Var1==0,]

ggplot(plot_freq_dist_chrom, aes(fill=Var2, y=Freq, x=Var1)) + 
  geom_bar(position="stack", stat="identity") +
  theme_bw()


# distance when on same chrom
dups <- dups %>%
  group_by(dup_family_fly) %>%
  mutate(farthest_dist = case_when(same_chrom_fam == 'same_chrom' ~ abs(min(dup_end_on_chrom) - max(dup_start_on_chrom)),
                          same_chrom_fam == 'diff_chrom(s)' ~ NA))

write.table(dups,'Duplicate_Proteins_2.tsv')


same_chrom <- dups[dups$same_chrom_fam == 'same_chrom',]

#gg <- 
  ggplot(same_chrom,aes(y=farthest_dist,x=dup_chrom,color=dup_chrom)) +
  geom_boxplot() +
  ylab('Distance between duplicates') +
  xlab('') +
  guides(color = guide_legend(title = "Chromosome")) +
  theme_bw() +
  ylim(0,50000)





  scale_y_log10() +
  geom_hline(yintercept = mean(chrom_lengths[chrom_lengths$chrom=='2L',]$length), linewidth = 1.1, color = "#F8766D") +
  geom_hline(yintercept = mean(chrom_lengths[chrom_lengths$chrom=='2R',]$length), linewidth = 1.1, color = "#A3A500") +
  geom_hline(yintercept = mean(chrom_lengths[chrom_lengths$chrom=='3L',]$length), linewidth = 1.1, color = "#00BF7D") +
  geom_hline(yintercept = mean(chrom_lengths[chrom_lengths$chrom=='3R',]$length), linewidth = 1.1, color = "#00B0F6") +
  geom_hline(yintercept = mean(chrom_lengths[chrom_lengths$chrom=='X',]$length), linewidth = 1.1, color = "#E76BF3")
  #ggsignif::geom_signif(comparisons = list(c('X','2L'),c('X','2R'),c('X','3L'),c('X','3R')),
  #                      map_signif_level = c('****'=0.0001,"***"=0.001, "**"=0.01, "*"=0.05),
  #                      textsize = 6, 
  #                      vjust = 0.5,
  #                      y_position = c(3.5e+07,3.75e+07,4e+07,4.25e+07),
  #                      test = 't.test') 
  
ggsave("./Plots/chrom_dist_between_dups_boxplot.jpg", plot = gg, width = 7, height = 8)



