

source("startup.R")


# import annotations to get intron information 
genes_w_introns <- read.csv("./Annotations.tsv", sep="")

# get a list of genes with introns 
genes_w_introns <- genes_w_introns[genes_w_introns$type=='intron',]
genes_w_introns <- genes_w_introns[c('gene_group','fly','type')]
genes_w_introns <- genes_w_introns[!duplicated(genes_w_introns),]

# import my duplicate genes 
dups <- read.csv("./Duplicate_Proteins.tsv", sep="")

# get intron presence information for my duplicates 
colnames(genes_w_introns) <- c('dup1','fly','dup1_intron')
dups <- left_join(dups,genes_w_introns,by=c('fly','dup1'))
dups$dup1_intron[is.na(dups$dup1_intron)] <- 'no_intron'

colnames(genes_w_introns) <- c('dup2','fly','dup2_intron')
dups <- left_join(dups,genes_w_introns,by=c('fly','dup2'))
dups$dup2_intron[is.na(dups$dup2_intron)] <- 'no_intron'


dups <- dups %>%
  mutate(mech = case_when((dup1_intron == 'intron') & (dup2_intron == 'intron') ~ 'DNA',
                          (dup1_intron == 'no_intron') & (dup2_intron == 'intron') ~ 'RNA',
                          (dup1_intron == 'intron') & (dup2_intron == 'no_intron') ~ 'RNA',
                          (dup1_intron == 'no_intron') & (dup2_intron == 'no_intron') ~ 'unk'))



mech <- as.data.frame(table(dups$mech))

ggplot(mech, aes(x="", y=Freq, group=Var1, fill=Var1)) +
  geom_bar(width = 1, stat = "identity", position = position_fill()) +
  geom_text(aes(label = Freq), position = position_fill(vjust = 0.5), colour = 'black', size = 3) +
  # facet_grid(.~factor(fly, levels=unique(mech$fly))) + 
  theme_void() +
  theme(legend.title = element_blank()) +
  coord_polar("y") 

mech_per_fly <- as.data.frame(table(dups$mech,dups$fly))
colnames(mech_per_fly) <- c('mech','fly','n')

t <- pivot_wider(mech_per_fly,names_from = 'mech',values_from = 'n')

t$dna_to_rna_ratio <- t$DNA / t$RNA








