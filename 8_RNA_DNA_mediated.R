

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
colnames(genes_w_introns) <- c('dup','fly','dup_intron','dup_chrom')
dups <- left_join(dups,genes_w_introns,by=c('fly','dup','dup_chrom'))
dups$dup_intron[is.na(dups$dup_intron)] <- 'no_intron'

# classify the duplicates into a mechanism using the introns
dups <- dups %>%
  group_by(dup_family_fly) %>%
  mutate(mech = case_when(
    all(dup_intron == "intron") ~ "dna",
    any(dup_intron == "intron") & any(dup_intron == "no_intron") ~ "rna",
    all(dup_intron == "no_intron") ~ "unknown",
    TRUE ~ NA)) %>%
  ungroup()


mech <- unique(dups[c('dup_family_fly','mech','fly')])


# make a pie chart of all the duplication mechanisms 
mech <- as.data.frame(table(mech$mech))

ggplot(mech, aes(x="", y=Freq, group=Var1, fill=Var1)) +
  geom_bar(width = 1, stat = "identity", position = position_fill()) +
  geom_text(aes(label = Freq), position = position_fill(vjust = 0.5), colour = 'black', size = 3) +
  # facet_grid(.~factor(fly, levels=unique(mech$fly))) + 
  theme_void() +
  theme(legend.title = element_blank()) +
  coord_polar("y") 


# compare the duplication mechanisms between flies 
mech_per_fly <- dups[c('mech','fly','dup_family')]
mech_per_fly <- mech_per_fly[!duplicated(mech_per_fly[c('mech','dup_family','fly')]),]
mech_per_fly <- as.data.frame(table(mech_per_fly$mech,mech_per_fly$fly))
colnames(mech_per_fly) <- c('mech','fly','n')


ggplot(mech_per_fly,aes(y=fly,x=mech,fill=n)) +
  geom_tile() +
  theme_bw()+
  xlab('Mechanism') +
  ylab('Fly') +
  labs(fill='Number of Duplicates')


# mech for families with all copies on same chromosome versus copies on different chromosomes
mech_chrom <- unique(dups[c('dup_family_fly','mech','same_chrom_fam')])

mech_chrom <- as.data.frame(table(mech_chrom$mech,mech_chrom$same_chrom_fam))

colnames(mech_chrom) <- c('mech','same_chrom','Freq')

ggplot(mech_chrom, aes(x="", y=Freq, group=same_chrom, fill=same_chrom)) +
  geom_bar(width = 1, stat = "identity", position = position_fill()) +
  geom_text(aes(label = Freq), position = position_fill(vjust = 0.5), colour = 'black', size = 3) +
  facet_grid(.~factor(mech)) + 
  theme_void() +
  theme(legend.title = element_blank()) +
  coord_polar("y") 


#


t <- as.data.frame(table(dups$dup_family_fly))
t <- as.data.frame(table(t$Freq))
max(t$Freq) / sum(t$Freq)


colnames(t) <- c('n_dups_in_fam','number_of_families_with')

ggplot(t, aes(x=n_dups_in_fam,y=number_of_families_with)) +
  geom_bar(stat='identity') +
  scale_x_discrete(limits=factor(2:40)) +
  theme_bw()




