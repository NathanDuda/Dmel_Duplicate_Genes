

# Author: Nathan Duda
# Purpose: Analyze the output of CNVSelectR

source("startup.R")


cnvselectr_output <- read.csv("C:/Users/17735/Downloads/Dmel_Duplicate_Genes/CNVSelectR_Output.tsv", sep="")

cnvselectr_output <- cnvselectr_output %>%
  group_by(group) %>%
  mutate(n_flies_in = n()) %>%
  mutate(signif = case_when(pval > 0.05 ~ 'No',
                            pval <= 0.05 ~ 'Yes'))


ggplot(cnvselectr_output, aes(x=ds, y=n_flies_in, color=signif)) +
  geom_jitter(size = 0.6) +
  scale_y_continuous(breaks = c(0, 10, 20, 30, 40, 50), labels = c(0, 10, 20, 30, 40, 50)) +
  xlim(0,1) +
  theme_bw() +
  labs(x='Ds', y='Frequency', color='Significant')


cnvselectr_plot <- ggplot(cnvselectr_output, aes(x=ds, y=n_flies_in, color=signif)) +
  geom_point() +
  scale_y_continuous(breaks = c(0, 10, 20, 30, 40, 50), labels = c(0, 10, 20, 30, 40, 50)) +
  xlim(0,1) +
  theme_bw() +
  labs(x='Ds', y='Frequency', color='Significant')

ggsave(cnvselectr_plot, file = './CNVSelectR_Output.jpeg', height = 6, width = 10)



# add chromosome info
genes <- read.csv("./Annotations_Gene.tsv", sep="")

cnvselectr_output <- cnvselectr_output %>%
  separate(name, into = c('dup_1','dup_2'), sep = '/')

genes <- genes[c('gene_group','start','end','chrom')]
colnames(genes) <- c('dup_1','dup_1_start','dup_1_end','dup_1_chrom')
cnvselectr_output <- merge(cnvselectr_output,genes,by='dup_1')

colnames(genes) <- c('dup_2','dup_2_start','dup_2_end','dup_2_chrom')
cnvselectr_output <- merge(cnvselectr_output,genes,by='dup_2')

chrom_counts <- as.data.frame(table(cnvselectr_output$signif,cnvselectr_output$dup_1_chrom, cnvselectr_output$dup_2_chrom))

# plot number of signif/not signif duplicate pairs per chrom 
ggplot(chrom_counts, aes(x=Var2, y=Var3, fill=Freq)) +
  geom_tile() +
  facet_wrap(~ Var1) +
  labs(x='',y='')

# keep only one unique signif value per equivalent gene group (each unique chrom combination is kept)
one_per_eq_group <- cnvselectr_output %>%
  distinct(group, signif, dup_1_chrom, dup_2_chrom) %>%
  select(-group)

one_per_eq_group <- as.data.frame(table(one_per_eq_group$signif, 
                                        one_per_eq_group$dup_1_chrom, 
                                        one_per_eq_group$dup_2_chrom))

# plot
ggplot(one_per_eq_group, aes(x=Var2, y=Var3, fill=Freq)) +
  geom_tile() +
  facet_wrap(~ Var1) +
  labs(x='',y='')

# check if most signif dup pairs are on same chromosome 
same_diff_chrom <- cnvselectr_output %>%
  mutate(same_diff_chrom = case_when(dup_1_chrom == dup_2_chrom ~ 'same',
                                     T ~ 'diff'))

table(same_diff_chrom$same_diff_chrom, same_diff_chrom$signif)


same_diff_chrom_one_per_eq_group <- one_per_eq_group %>%
  mutate(same_diff_chrom = case_when(Var2 == Var3 ~ 'same',
                                     T ~ 'diff')) %>%
  group_by(Var1,same_diff_chrom) %>%
  summarize(n = sum(Freq), .groups='keep')



# add fbgns
fbgn_fbpp <- read.delim("./fbgn_fbtr_fbpp_expanded_fb_2023_06.tsv", header=FALSE, comment.char="#")



fly_name_list <- read.table("./Individuals_List_47.txt", quote="\"", comment.char="")
my_FBgns <- as.data.frame(matrix(ncol=4,nrow=0))

# loop over all flies 
for (row in 1:nrow(fly_name_list)){
  
  # get name of current fly 
  name <- fly_name_list[row,1]
  
  # import blast result
  blastp_FBgn <- read.delim(paste0("./FBgn_Blast_Output/",name,"_blastp.tsv"), header=FALSE)
  colnames(blastp_FBgn) <- c('qseqid', 'sseqid', 'length', 'bitscore', 'pident', 'evalue')
  
  # keep only the gene with its FBgn and the evalue
  #blastp_FBgn <- blastp_FBgn[,c('qseqid','sseqid','evalue')]
  #colnames(blastp_FBgn) <- c('gn','FBpp','evalue')
  
  blastp_FBgn <- blastp_FBgn %>%
    group_by(qseqid) %>%
    filter(length >= 100 & pident >= 90) %>%
    filter(evalue == min(evalue)) %>%
    ungroup() %>%
    distinct(qseqid,evalue, .keep_all = T) %>% # keep only one if lowest evalue is same for multiple 
    select(qseqid,sseqid)
  # same FBpp can match to two different genes 
  
  blastp_FBgn$qseqid <- paste0(name,'_',blastp_FBgn$qseqid)
  my_FBgns <- rbind(my_FBgns,blastp_FBgn)
}



fbgn_fbpp <- fbgn_fbpp[c('V3','V8','V10')]
colnames(fbgn_fbpp) <- c('fbgn','fbtr','fbpp')

colnames(my_FBgns) <- c('gn','fbpp')
my_FBgns <- merge(my_FBgns,fbgn_fbpp,by='fbpp')

t <- as.data.frame(my_FBgns$fbgn[!duplicated(my_FBgns$fbgn)])

write.table(t, file = './my_FBgn_list.txt', row.names = FALSE, col.names = FALSE)

# separating FBgns by sig / not sig in cnvselectr output

dup_1 <- cnvselectr_output[c('dup_1','signif')]
dup_2 <- cnvselectr_output[c('dup_2','signif')]

colnames(dup_1) <- c('dup','signif')
colnames(dup_2) <- c('dup','signif')

dup <- rbind(dup_1,dup_2)
colnames(dup) <- c('gn','signif')

x <- merge(dup,my_FBgns,by='gn')
not_sig <- x[x$signif == 'No',]
sig <- x[x$signif == 'Yes',]

write.table(not_sig$fbgn, file = './NotSig_FBgn_list.txt', row.names = FALSE, col.names = FALSE)
write.table(sig$fbgn, file = './Sig_FBgn_list.txt', row.names = FALSE, col.names = FALSE)


#
#exp <- read.delim("./high-throughput_gene_expression_fb_2023_06.tsv", header=FALSE, comment.char="#")
#exp <- exp[c('V5','V6','V7','V9')]

#t <- pivot_wider(names_from = 'V5')

write.table(x, './gn_fbgn_under_selection_or_not.tsv')


YO_dmel_exp <- read.delim("C:/Users/17735/Downloads/Eight_Species/Raw_Data/YO_Expression/GSE99574_HiSAT2_dmel.nrc.YO.txt")


YO_dmel_exp <- YO_dmel_exp %>% 
  select(-jaccard, -YOgnID) %>%
  rename_all(~ ifelse(startsWith(., "w1118"), paste0(., "1"), .)) %>% 
  rename_all(~ gsub(paste0('w1118_'), "", .)) %>%
  rename_all(~ gsub(paste0('oreR_'), "", .)) %>% 
  rename_all(~ gsub("dmel_", "", .)) #%>%
 # mutate(across(2:ncol(YO_dmel_exp), as.numeric)) 


tissue_names <- c('f_ac','f_dg','f_go','f_hd','f_re','f_tx','f_wb',
                  'm_ac','m_dg','m_go','m_hd','m_re','m_tx','m_wb')

for (tissue in tissue_names) {
  YO_dmel_exp <- YO_dmel_exp %>%
    mutate("{tissue}" := rowMeans(select(., starts_with(tissue)), na.rm = TRUE))
}


YO_dmel_exp <- YO_dmel_exp %>%
  select(FBgnID,all_of(tissue_names)) %>%
  filter(rowSums(select(., 2:15)) > 0) %>%
  filter(!FBgnID=='-')

tau <- YO_dmel_exp

rownames(YO_dmel_exp) <- YO_dmel_exp$FBgnID
YO_dmel_exp <- YO_dmel_exp %>% select(-FBgnID)
YO_dmel_exp <- YO_dmel_exp / rowSums(YO_dmel_exp)
YO_dmel_exp$fbgn <- rownames(YO_dmel_exp)


gn_YOexp <- merge(my_FBgns,YO_dmel_exp,by='fbgn')



gn_YOexp <- gn_YOexp[c('gn','f_ac','f_dg','f_go','f_hd','f_tx','f_wb','m_ac','m_dg','m_go',
                       'm_hd','m_re','m_tx','m_wb')]
colnames(gn_YOexp) <- c('dup_1','dup1_f_ac','dup1_f_dg','dup1_f_go','dup1_f_hd','dup1_f_tx','dup1_f_wb','dup1_m_ac','dup1_m_dg','dup1_m_go',
                        'dup1_m_hd','dup1_m_re','dup1_m_tx','dup1_m_wb')
t <- merge(cnvselectr_output,gn_YOexp,by='dup_1')



colnames(gn_YOexp) <- c('dup_2','dup2_f_ac','dup2_f_dg','dup2_f_go','dup2_f_hd','dup2_f_tx','dup2_f_wb','dup2_m_ac','dup2_m_dg','dup2_m_go',
                        'dup2_m_hd','dup2_m_re','dup2_m_tx','dup2_m_wb')

t <- merge(t,gn_YOexp,by='dup_2')

p <- t[c(7,14:39)]

p <- pivot_longer(p, cols=c(2:27))

p$name <- gsub('dup1_','',p$name)
p$name <- gsub('dup2_','',p$name)

p <- p %>%
  group_by(name,signif) %>%
  summarize(value = mean(value), .groups='keep')


ggplot(p, aes(x=name, fill=value, y=signif)) +
  geom_tile() +
  theme_bw()

chrom_counts_d1 <- t[c('dup_1_chrom','signif')]
chrom_counts_d2 <- t[c('dup_2_chrom','signif')]
colnames(chrom_counts_d1) <- c('chrom','signif')
colnames(chrom_counts_d2) <- c('chrom','signif')
chrom_counts <- rbind(chrom_counts_d1,chrom_counts_d2)

ggplot(chrom_counts, aes(x=chrom, fill=signif)) +
  geom_bar() +
  theme_bw()
ggplot(chrom_counts, aes(x=chrom, fill=signif)) +
  geom_bar(position = 'fill') +
  theme_bw() +
  ylab('proportion')

# expression by chromosome
exp_sig_chrom <- t[c(7,10,13:39)]
exp_sig_chrom_1 <-  exp_sig_chrom %>% select(signif, dup_1_chrom, contains('dup1'))
exp_sig_chrom_2 <-  exp_sig_chrom %>% select(signif, dup_2_chrom, contains('dup2'))
colnames(exp_sig_chrom_2) <- colnames(exp_sig_chrom_1)

exp_sig_chrom <- rbind(exp_sig_chrom_1,exp_sig_chrom_2)

exp_sig_chrom <- pivot_longer(exp_sig_chrom, cols = c(3:15))
exp_sig_chrom <- exp_sig_chrom %>% mutate(name = gsub('dup1_','',name))


ggplot(exp_sig_chrom, aes(x = name, y = dup_1_chrom, fill = value)) +
  geom_tile() +
  facet_grid(.~signif)




# tau



rownames(tau) <- tau$FBgnID
tau <- tau %>% select(-FBgnID)

library(tispec)

tau <- calcTau(tau)
tau$fbgn <- rownames(tau)

tau <- tau[c('fbgn','tau')]

t <- merge(x,tau,by='fbgn')


ggplot(t, aes(x=signif, y=tau)) +
  geom_boxplot() +
  geom_signif(comparisons = list(c("No", "Yes")), test = "t.test", map_signif_level = TRUE)


# ppi



# import the three ppi data
ppi_1 <- read.delim("./PPI_Data/FlyBi_data_20210205.csv")
ppi_1 <- ppi_1[c('FBgn1','FBgn2')]

ppi_2 <- read.delim("./PPI_Data/fly_other_physical.txt")
ppi_2 <- ppi_2[c('FLY_GENE1','FLY_GENE2')]

ppi_3 <- read.delim("./PPI_Data/flybase_ppi.txt")
ppi_3 <- ppi_3[c('FLY_GENE1','FLY_GENE2')]

# combine all three datasets
ppi_2 <- rbind(ppi_2,ppi_3)
colnames(ppi_2) <- colnames(ppi_1)
ppi <- rbind(ppi_1,ppi_2)


ppi <- ppi[!duplicated(ppi),]

ppi_1 <- ppi[1]
ppi_2 <- ppi[2]

colnames(ppi_1) <- 'ppi'
colnames(ppi_2) <- 'ppi'

ppi <- rbind(ppi_1,ppi_2)

ppi <- as.data.frame(table(ppi$ppi))
colnames(ppi) <- c('fbgn','ppi')


p <-  merge(x,ppi,by='fbgn')


ggplot(p, aes(x=signif, y=ppi)) +
  geom_boxplot() +
  geom_signif(comparisons = list(c("No", "Yes")), test = "t.test", map_signif_level = TRUE)


# duplication mechanism


annotations <- read.csv("C:/Users/17735/Downloads/Dmel_Duplicate_Genes/Annotations.tsv", sep="")

exon_counts <- annotations %>%
  select(gene_group, type) %>%
  filter(type == 'exon') %>%
  select(-type) %>%
  group_by(gene_group) %>%
  mutate(n_exons = n()) %>%
  distinct() 

colnames(exon_counts) <- c('dup_1','dup_1_exons')
cnvselectr_output <- merge(exon_counts, cnvselectr_output, by = 'dup_1')

colnames(exon_counts) <- c('dup_2','dup_2_exons')
cnvselectr_output <- merge(exon_counts, cnvselectr_output, by = 'dup_2')


cnvselectr_output <- cnvselectr_output %>%
  mutate(mech = case_when(dup_1_exons > 1 & dup_2_exons == 1 ~ 'RNA_dup_2',
                          dup_2_exons > 1 & dup_1_exons == 1 ~ 'RNA_dup_1',
                          dup_1_exons > 1 & dup_2_exons > 1 ~ 'DNA',
                          dup_2_exons > 1 & dup_1_exons > 1 ~ 'DNA',
                          dup_2_exons == 1 & dup_1_exons == 1 ~ 'unknown'))
  
  
  
t <- as.data.frame(table(cnvselectr_output$mech, cnvselectr_output$n_flies_in))

ggplot(t, aes(x=Var2, y = Freq, group=Var1, color = Var1)) +
  geom_line()

ggplot(cnvselectr_output, aes(x=n_flies_in, fill=mech)) +
  geom_histogram(position = 'fill')

t <- as.data.frame(table(cnvselectr_output$mech, cnvselectr_output$dup_1_chrom, cnvselectr_output$dup_2_chrom))

ggplot(t, aes(x=Var2, y=Var3, fill=Freq)) +
  geom_tile() +
  facet_wrap(~ Var1) +
  labs(x='',y='')

t <- as.data.frame(table(cnvselectr_output$mech, cnvselectr_output$signif))

ggplot(t, aes(x=Var2, y=Var1, fill = Freq)) +
  geom_tile()



cnvselectr_output <- cnvselectr_output %>%
  group_by(group) %>%
  mutate(n_DNA_in_group = sum(mech == 'DNA'),
         n_unk_in_group = sum(mech == 'unknown'),
         n_RNA_in_group = sum(mech %in% c('RNA_dup_1', 'RNA_dup_2')))

# get percentage of each mechanism 
perc_mech <- cnvselectr_output %>%
  distinct(group, .keep_all = TRUE) %>%
  ungroup() %>%
  summarize(dna = sum(n_DNA_in_group),
            rna = sum(n_RNA_in_group),
            unk = sum(n_unk_in_group)) %>%
  t() %>%
  as.data.frame() %>%
  mutate(perc = V1 / sum(V1))

# check if DNA mediated duplicates are next to each other
same_chrom <- cnvselectr_output %>%
  filter(dup_1_chrom == dup_2_chrom) %>%
  mutate(distance = case_when(dup_1_start > dup_2_start ~ dup_1_start - dup_2_end,
                              dup_2_start > dup_1_start ~ dup_2_start - dup_1_end)) %>%
  mutate(combined_mech = case_when(str_detect(mech, 'RNA') ~ 'RNA', T ~ mech))

same_chrom %>% filter(combined_mech != 'unknown') %>%
  ggplot(., aes(x = dup_1_chrom, y = distance, color = combined_mech)) +
    geom_boxplot() +
    geom_point(position = position_dodge(width = 0.75)) +
    scale_y_log10() + 
    theme_bw()

diff_chrom <- cnvselectr_output %>%
  mutate(chrom_pair = paste0(pmin(dup_1_chrom, dup_2_chrom), "_", pmax(dup_1_chrom, dup_2_chrom))) %>%
  filter(dup_1_chrom != dup_2_chrom) %>%
  group_by(chrom_pair) %>%
  summarise(count = n())


# compare protein lengths
lengths <- cnvselectr_output %>%
  mutate(dup_1_len = dup_1_end - dup_1_start,
         dup_2_len = dup_2_end - dup_2_start)

ggplot(lengths, aes(x = signif, y = dup_1_len)) +
  geom_boxplot()



