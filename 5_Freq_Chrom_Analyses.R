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

###################
annotations <- all_annotations
annotations$loc <- (annotations$len / 2) + annotations$start
annotations <- annotations[c('gene_group','chrom','loc')]



dups <- read.csv("./Duplicate_Proteins_3.tsv", sep="")


dups$loc <- (dups$dup_length / 2) + dups$dup_start_on_chrom

k <- dups %>%
  group_by(dup_family,dup_chrom) %>%
  summarize(loc = max(loc), .groups = 'keep')

k <- dups[dups$fly=='A2',]

dups$gene_group <- paste0(dups$fly, '_', dups$dup)
dups <- dups[c('gene_group','dup_chrom','loc')]
colnames(dups) <- c('gene_group','chrom','loc')
annotations$type <- 'all'
dups$type <- 'dup'
both <- rbind(dups,annotations)


gg <- ggplot(both,aes(x=loc,color=type)) +
  geom_density(bw=200000) + 
  facet_grid(.~chrom) +
  theme_bw() + theme(legend.position = "none")
ggsave("./trash.jpg", plot = gg, width = 20, height = 10)

###################


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



###########################################################################


equiv_genes <- read.csv("./Equivalent_Genes.tsv", sep="")
dups <- read.csv("./Duplicate_Proteins.tsv", sep="")

t <- equiv_genes %>% sample_n(10000, replace = FALSE)


t <- t %>%
  mutate(rowname=rownames(.)) %>%
  mutate(row_group = paste0('group_',rowname)) %>%
  select(-rowname) %>%
  pivot_longer(cols=c(A1:TOM.008)) %>%
  separate_rows(value, sep = ",\\s*")

t <- na.omit(t)

p <- as.data.frame(table(t$row_group,t$name))
k <- as.data.frame(table(p$Var1,p$Freq))


k <- k[k$Freq!=0 & k$Var2!=0,]

k <- as.data.frame(table(k$Freq))

ggplot(k,aes(x=Var1,y=Freq)) +
  geom_bar(stat='identity') +
  #scale_x_discrete(limits=factor(1:52)) + # not 52 without crappy assemblies 
  xlab('Number of flies duplicate is present in') +
  ylab('Freq') +
  theme_bw()



all_annotations <- read.csv("./Annotations_Gene.tsv", sep="")
all_annotations <- all_annotations[!all_annotations$fly %in% c('I23','ZH26','N25','T29A','B59'),]
all_annotations$loc <- (all_annotations$len / 2) + all_annotations$start
gn_chroms <- all_annotations[c('gene_group','loc','chrom')]
p <- merge(t,gn_chroms,by='gene_group')


k2 <- as.data.frame(table(p$row_group,p$name,p$chrom))
k2 <- as.data.frame(table(k2$Var1,k2$Var3,k2$Freq))

k2 <- k2[k2$Freq!=0 & k2$Var3!=0,]



#
k2 <- as.data.frame(table(k2$Freq,k2$Var2))


ggplot(k2, aes(fill=Var2, y=Freq, x=Var1)) + 
  geom_bar(position="fill", stat="identity") +
  theme_bw()



##############################################################



# get the number of duplicates each gene has in its family
dups <- dups %>%
  select(dup, dup_family) %>%
  group_by(dup_family) %>%
  mutate(n_copies_in_fam = n())

dups <- merge(equiv_genes, dups, by.x = 'gene_group', by.y = 'dup')


dups_equiv_counts <- dups %>%
  mutate_at(vars(2:48), ~ str_count(as.character(.), ",")) %>%
  mutate_at(vars(2:48), ~ . + 1)


# n duplicated same amount of copies 

n_same_copies <- dups_equiv_counts %>%
  mutate_at(vars(2:48), ~ coalesce(., n_copies_in_fam)) %>%
  rowwise() %>%
  mutate(n_same_n_duplicated = sum(across(2:48) == n_copies_in_fam, na.rm = TRUE)) 
  #%>%
  #mutate(n_same_n_duplicated = n_same_n_duplicated - 1)


n_same_copies <- as.data.frame(table(n_same_copies$n_same_n_duplicated))


plot_n_same_copies <- ggplot(n_same_copies,aes(x=Var1,y=Freq)) +
  geom_bar(stat='identity') +
  #scale_x_discrete(limits=factor(1:52)) + # not 52 without crappy assemblies 
  xlab('Number of flies duplicate is present in') +
  ylab('Freq') +
  theme_bw()

###

plot_n_same_copies



# n duplicated 

n_more_than_1_copy <- dups_equiv_counts %>%
  mutate_at(vars(2:48), ~ coalesce(., n_copies_in_fam)) %>%
  rowwise() %>%
  mutate(n_same_n_duplicated = sum(across(2:48) > 1, na.rm = TRUE))


n_more_than_1_copy <- as.data.frame(table(n_more_than_1_copy$n_same_n_duplicated))


plot_n_more_than_1_copy <- ggplot(n_more_than_1_copy,aes(x=Var1,y=Freq)) +
  geom_bar(stat='identity') +
  #scale_x_discrete(limits=factor(1:52)) + # not 52 without crappy assemblies 
  xlab('Number of flies duplicate is present in') +
  ylab('Freq') +
  theme_bw()



plot_n_same_copies # number of flies with the same amount of copies 
n_more_than_1_copy # number of flies with >1 copy (can be 2 copies while the rest have 300)
# bars add up to the 41,611 duplicate copies 
#########################################################################
# frequency plot 

ggsave("./Plots/freq_distribution_barplot.jpg", plot = plot_n_same_copies, width = 10, height = 3)



#t <- n_same_copies[c('dup_family','n_same_n_duplicated')]
#t <- t[!duplicated(t),]
#t_plot <- as.data.frame(table(t$n_same_n_duplicated))
#ggplot(t_plot,aes(x=Var1,y=Freq)) +
  #geom_bar(stat='identity') +
  #xlab('Number of flies duplicate is present in') +
  #ylab('Freq') +
  #theme_bw()








# frequency distribution with chromosomes 




# add chrom info

dup_chrom <- read.csv("./Duplicate_Proteins.tsv", sep="")
dup_chrom <- dup_chrom[c('dup','dup_chrom')]

colnames(dups_equiv_counts)[1] <- 'dup'
t <- merge(dup_chrom,dups_equiv_counts,by='dup')


t <- t %>%
  mutate_at(vars(3:49), ~ coalesce(., n_copies_in_fam))

sd <- t %>%
  rowwise() %>%
  mutate(sd = sd(c_across(3:49), na.rm = TRUE)) %>%
  mutate(identitcal_n_copies_across = case_when(all(c_across(3:49) == n_copies_in_fam) ~ T,
                                                T ~ F)) %>%
  mutate(duplicated_across = case_when(all(c_across(3:49) > 1) ~ T,
                                       T ~ F))

sd <- sd %>%
  mutate(mean = mean(c_across(3:49), na.rm = TRUE))












#


colnames(all_annotations)[1] <- 'id'
gene_chrom <- all_annotations[c('id','chrom')]


#equiv_genes_orig <- equiv_genes
equiv_genes <- equiv_genes_orig

#equiv_genes <- equiv_genes[c(12000:14000),]



result <- data.frame()

p = 0 

names <- c(
  'A1','A2','A3','A4','A5','A6','A7','AB8','AKA-017','AKA-018','B1','B2','B3','B4','B6','COR-014',
  'COR-018','COR-023','COR-025','GIM-012','GIM-024','ISO-1','JUT-008','JUT-011','KIE-094','LUN-004',
  'LUN-007','MUN-008','MUN-009','MUN-013','MUN-015','MUN-016','MUN-020','ORE','RAL-059','RAL-091',
  'RAL-176','RAL-177','RAL-375','RAL-426','RAL-737','RAL-855','SLA-001','STO-022','TEN-015','TOM-007',
  'TOM-008')

for (name in names){
  
  
  eq <- equiv_genes[c('gene_group',paste0(gsub('-','.',name)))]
  eq <- na.omit(eq)
  
  
  eq <- eq %>%
    mutate_at(2, ~strsplit(as.character(.), ",\\s*")) %>%
    pivot_longer(cols = 2, names_to = "column", values_to = "value") %>%
    unnest(value)
  
  colnames(eq) <- c('gene_group','column','id')
  
  
  
  trash <- merge(eq,gene_chrom,by='id')
  trash <- trash[c('gene_group','chrom')]
  trash <- as.data.frame(table(trash$gene_group,trash$chrom))
  colnames(trash) <- c('gene_group','chrom',paste0(name))
  

  
  if (name == names[1]){result <- trash}
  if (name != names[1]){result <- left_join(result,trash,by=c('gene_group','chrom'))}
  
  
  p = p + 1
  print(p)
  
}


#write.table(result,file='Equivalent_Genes_Chromosomes.tsv')
result <- read.csv("./Equivalent_Genes_Chromosomes.tsv", sep="")




t <- result %>%
  rowwise() %>%
  mutate(sd = sd(c_across(3:49), na.rm = TRUE),
         mean = mean(c_across(3:49), na.rm = TRUE))

 
# distibution of variance per chromosome 

ggplot(t,aes(x=sd,color=chrom)) +
  geom_density() +
  xlim(0,.5)







p <- t %>% 
  sample_n(100000, replace = FALSE)

gg <- ggplot(p,aes(x=gene_group,y=sd,color=chrom)) +
  geom_point() + 
  facet_grid(.~chrom)

ggsave("./Plots/chrom.jpg", plot = gg, width = 15, height = 10)


#

#write.table(t,file='Equivalent_Genes_Chromosomes_sd_mean.tsv')
t <- read.csv("./Equivalent_Genes_Chromosomes_sd_mean.tsv", sep="")

any(duplicated(t[c('gene_group','chrom')]))


#






t <- pivot_longer(result,cols=c(3:49))
#write.table(t,file='pivot.tsv')

#
t <- read.csv("./pivot.tsv",sep='')



c1 <- t[t$chrom=='2L',]



gg <- ggplot(c1,aes(x=name,y=gene_group,fill=value)) +
  geom_tile() 

ggsave("./Plots/chrom.jpg", plot = gg, width = 15, height = 10)




#


plot_freq_dist_chrom <- dups[c('dup','dup_chrom','')]

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



