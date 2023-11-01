

source("startup.R")


# import annotations to get intron information 
genes_w_introns <- read.csv("./Annotations.tsv", sep="")

# get only genes with introns between their start and stop codons 
genes_w_introns <- genes_w_introns %>%
  filter(type == 'start_codon' | type == 'stop_codon' | type == 'intron') %>%
  group_by(gene_group, fly) %>%
  filter((type == "intron" &
         ((strand == '+' & (lag(type == "start_codon") | lead(type == "stop_codon"))) |
          (strand == '-' & (lag(type == "stop_codon") | lead(type == "start_codon"))))))

# get a list of genes with introns 
genes_w_introns <- genes_w_introns[c('gene_group','fly','type','chrom')]
genes_w_introns <- genes_w_introns[!duplicated(genes_w_introns),]

# import my duplicate genes 
dups <- read.csv("./Duplicate_Proteins_2.tsv", sep="")

# get intron presence information for my duplicates 
colnames(genes_w_introns) <- c('dup1','fly','dup1_intron','dup1_chrom')
dups <- left_join(dups,genes_w_introns,by=c('fly','dup1','dup1_chrom'))
dups$dup1_intron[is.na(dups$dup1_intron)] <- 'no_intron'

colnames(genes_w_introns) <- c('dup2','fly','dup2_intron','dup2_chrom')
dups <- left_join(dups,genes_w_introns,by=c('fly','dup2','dup2_chrom'))
dups$dup2_intron[is.na(dups$dup2_intron)] <- 'no_intron'

# classify the duplicates into a mechanism using the introns
dups <- dups %>%
  mutate(mech = case_when((dup1_intron == 'intron') & (dup2_intron == 'intron') ~ 'DNA',
                          (dup1_intron == 'no_intron') & (dup2_intron == 'intron') ~ 'RNA',
                          (dup1_intron == 'intron') & (dup2_intron == 'no_intron') ~ 'RNA',
                          (dup1_intron == 'no_intron') & (dup2_intron == 'no_intron') ~ 'unk'))

# make a pie chart of all the duplication mechanisms 
mech <- as.data.frame(table(dups$mech))

ggplot(mech, aes(x="", y=Freq, group=Var1, fill=Var1)) +
  geom_bar(width = 1, stat = "identity", position = position_fill()) +
  geom_text(aes(label = Freq), position = position_fill(vjust = 0.5), colour = 'black', size = 3) +
  # facet_grid(.~factor(fly, levels=unique(mech$fly))) + 
  theme_void() +
  theme(legend.title = element_blank()) +
  coord_polar("y") 

# compare the duplication mechanisms between flies 
mech_per_fly <- as.data.frame(table(dups$mech,dups$fly))
colnames(mech_per_fly) <- c('mech','fly','n')


ggplot(mech_per_fly,aes(y=fly,x=mech,fill=n)) +
  geom_tile() +
  theme_bw()+
  xlab('Mechanism') +
  ylab('Fly') +
  labs(fill='Number of Duplicates')



# get paralog families

dups_ids <- dups[c('dup1','dup2')]


network <- igraph::graph.data.frame(dups_ids, directed = FALSE)


plot(network, layout = layout.fruchterman.reingold)

network_connections <- dups_ids %>%
  rownames_to_column(var = "Connection") %>%
  gather(key = "Direction", value = "Node", -Connection)

# Reorder the columns to have "Connection" first, then "Direction", and finally "Node"
network_connections <- network_connections[, c("Connection", "Direction", "Node")]

print(network_connections)



#






#################
t <- pivot_wider(mech_per_fly,names_from = 'mech',values_from = 'n')

t$dna_to_rna_ratio <- t$DNA / t$RNA


ratio <- t[c('fly','dna_to_rna_ratio')]
ratio$x <- 'dna_to_rna_ratio'

ggplot(ratio,aes(y=fly,x=x,fill=dna_to_rna_ratio)) +
  geom_tile()



same_chrom <- dups[dups$same_chrom == 'same_chrom',]

ggplot(same_chrom,aes(x=mech,y=dist)) +
  geom_boxplot()



table(dups$same_chrom,dups$mech)

#             DNA     RNA     unk
# diff_chrom  71,703  14,237  4,314
# same_chrom  29,896  5,691   6,501


trash <- as.data.frame(table(dups$same_chrom,dups$mech,dups$dup1_chrom,dups$dup2_chrom))




s <- genes_w_introns_orig


s <- s[s$type %in% c('start_codon','stop_codon','intron'),]
s2 <- s %>%
  group_by(gene_group,fly) %>%
  filter((type == "intron" & strand == '+' & (lag(type == "stop_codon") | lead(type == "start_codon"))) |
         (type == "intron" & strand == '-' & (lag(type == "start_codon") | lead(type == "stop_codon"))))







