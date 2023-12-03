


source("startup.R")





# import duplicates for dsim 
dsim_dups <- read.csv("./Dsim/Duplicate_Proteins.tsv", sep="")

# import duplicates for dmel 
dmel_dups <- read.csv("./Duplicate_Proteins_3.tsv", sep="")
dmel_dups_orig <- dmel_dups

# compare number of duplicate families between dmel flies and dsim 
length(unique(dsim_dups$dup_family))
dmel_dups %>%
  group_by(fly) %>%
  summarise(length(unique(dup_family)))


# read in list with the 47 fly names 
fly_name_list <- read.table("./Individuals_List_47.txt", quote="\"", comment.char="")

# create empty dataframes for output
any_ortholog <- as.data.frame(matrix(ncol=13,nrow=0))
all_ortholog <- as.data.frame(matrix(ncol=13,nrow=0))
dups_orthos <- as.data.frame(matrix(ncol=16,nrow=0))
all_orthologs <- as.data.frame(matrix(ncol=4,nrow=0)) 

# loop over all flies 
for (row in 1:nrow(fly_name_list)){
  
  # get name of current fly 
  name <- fly_name_list[row,1]
  
  orthologs <- read.delim2(paste0("./OrthoFinder_Outputs/",name,"_Output/Results_Nov10/Orthologues/Orthologues_",
                                  name,"_annotations/",name,"_annotations__v__Dsim_annotations.tsv"))
  
  # format orthologs into each pair
  orthologs <- orthologs %>% separate_rows(paste0(gsub('-','.',name),'_annotations'), sep = ",")
  orthologs <- orthologs %>% separate_rows(Dsim_annotations, sep = ",")
  
  colnames(orthologs) <- c('orthogroup','dup','dsim_ortholog')
  
  orthologs$dup <- trimws(orthologs$dup)
  orthologs$dsim_ortholog <- trimws(orthologs$dsim_ortholog)
  
  
  # merge with dups
  dmel_dups <- dmel_dups_orig[dmel_dups_orig$fly==name,]
  
  dmel_dups <- left_join(dmel_dups,orthologs,by='dup')
  
  dmel_dups <- dmel_dups %>%
    mutate(dsim_dup = case_when(dsim_ortholog %in% dsim_dups$dup ~ T,
                                T ~ F))
  
  # get duplicates where at least one copy in family has dsim ortholog
  any_ortholog_subset <- dmel_dups %>%
    group_by(dup_family) %>%
    filter(any(dsim_dup))
  
  # get duplicates where every copy in family has dsim ortholog
  all_ortholog_subset <- dmel_dups %>%
    group_by(dup_family) %>%
    filter(all(dsim_dup))
  
  # combine with all duplicates 
  any_ortholog <- rbind(any_ortholog,any_ortholog_subset)
  all_ortholog <- rbind(all_ortholog,all_ortholog_subset)
  dups_orthos <- rbind(dups_orthos,dmel_dups)
  
  orthologs$fly <- name
  all_orthologs <- rbind(all_orthologs,orthologs)
  
  print(row)

}

###########

# write to table
write.table(dups_orthos,file='./Dmel_Dsim_Ortho_Duplicates.tsv')
dups_orthos <- read.csv("./Dmel_Dsim_Ortho_Duplicates.tsv", sep="")

# number of orthologs 
t <- all_orthologs %>%
  group_by(fly) %>%
  summarise(n_orthologs = length(unique(dsim_ortholog)))



dups_orthos <- dups_orthos[!duplicated(dups_orthos[c('fly','dup')]),]




# frequency plot using any_orthologs 
freq_dist <- any_ortholog[c('dup','fly')]

freq_dist <- as.data.frame(table(freq_dist$dup))
plot_freq_dist <- as.data.frame(table(freq_dist$Freq))

ggplot(plot_freq_dist,aes(x=Var1,y=Freq)) +
  geom_bar(stat='identity') +
  #scale_x_discrete(limits=factor(1:52)) + # not 52 without crappy assemblies 
  xlab('Number of flies duplicate is present in') +
  ylab('Freq') +
  theme_bw()



# frequency plot using all_orthologs 
freq_dist <- all_ortholog[c('dup','fly')]

freq_dist <- as.data.frame(table(freq_dist$dup))
plot_freq_dist <- as.data.frame(table(freq_dist$Freq))

ggplot(plot_freq_dist,aes(x=Var1,y=Freq)) +
  geom_bar(stat='identity') +
  #scale_x_discrete(limits=factor(1:52)) + # not 52 without crappy assemblies 
  xlab('Number of flies duplicate is present in') +
  ylab('Freq') +
  theme_bw()


####

# classify dmel duplicates based on if they have dsim orthologs
dups_orthos <- dups_orthos %>%
  group_by(dup_family,fly) %>%
  mutate(any_ortholog = case_when(any(dsim_dup) ~ T, T ~ F)) %>%
  mutate(dsim_ortho = case_when(!is.na(dsim_ortholog) ~ T, is.na(dsim_ortholog) ~ F)) %>%
  mutate(all_ortholog =  case_when(all(dsim_dup) ~ T, T ~ F)) %>%
  mutate(dsim_dup = case_when(dsim_dup ~ 'Dsim_Ortholog_Duplicated',
                              !dsim_dup ~ 'No_Duplicated_Ortholog'))


# make a frequency distribution plot 
plot_freq_dist_colored <- dups_orthos[c('dup','fly','dsim_dup')]


#########
fly_counts <- plot_freq_dist_colored %>%
  group_by(dup) %>%
  mutate(x_axis = n_distinct(fly)) %>%
  ungroup() %>%
  distinct(dup,dsim_dup, .keep_all = T) %>%
  group_by(x_axis,dsim_dup) %>%
  mutate(y_axis = n()) %>%
  ungroup() %>%
  distinct(dsim_dup,x_axis,y_axis)


freq_dist <- ggplot(fly_counts, aes(fill = factor(dsim_dup,levels=c('No_Duplicated_Ortholog','Dsim_Ortholog_Duplicated')), 
                       x = x_axis, y = y_axis)) + 
  geom_bar(position = "stack", stat='identity') +
  scale_x_discrete(limits = factor(1:47)) +
  theme_bw() +
  ylab('Frequency') +
  xlab('') +
  scale_fill_manual(values = c("hotpink3", "cornflowerblue"), name = "") +
  guides(fill = guide_legend(title = "", override.aes = list(fill = c("hotpink3", "cornflowerblue"))))



plot_freq_dist_percent <- fly_counts %>%
  group_by(x_axis) %>%
  mutate(Percentage = y_axis / sum(y_axis) * 100)

perc <- ggplot(plot_freq_dist_percent, aes(x = x_axis, y = Percentage, fill = factor(dsim_dup,levels=c('No_Duplicated_Ortholog','Dsim_Ortholog_Duplicated')))) +
  geom_bar(stat = "identity", position = "stack") +
  scale_x_discrete(limits=factor(1:47)) +
  theme_bw() +
  xlab('') +
  ylab('Percentage of \n Duplicates ') +
  scale_fill_manual(values = c("hotpink3", "cornflowerblue"), name = "") +
  guides(fill = guide_legend(title = "", override.aes = list(fill = c("hotpink3", "cornflowerblue"))))

  
gg <- ggarrange(freq_dist,perc,nrow = 2, ncol = 1,
          labels = c('A','B'),
          align = 'v')

ggsave("./Plots/dsim_freq_dist_barplot.jpg", plot = gg, width = 12, height = 6)


############################




plot_freq_dist_colored <- dups_orthos[c('dup','fly','dsim_dup','dup_chrom')]


all_chrom <- plot_freq_dist_colored %>%
  select(-dsim_dup) %>%
  group_by(dup,fly) %>%
  mutate(x_axis = n())
  
  







table(plot_freq_dist_colored$fly,plot_freq_dist_colored$dup_chrom)



all_chrom <- 
  plot_freq_dist_colored %>%
  group_by(dup) %>%
  mutate(x_axis = n_distinct(fly)) %>%
  ungroup() %>%
  #distinct(dup,fly,dsim_dup,dup_chrom, .keep_all = T) %>%
  group_by(x_axis,dsim_dup,dup_chrom) %>%
  mutate(y_axis = n()) %>%
  ungroup() %>%
  distinct(dsim_dup,dup_chrom,x_axis,y_axis) %>%
  group_by(dup_chrom,x_axis) %>%
  mutate(y_axis=sum(y_axis)) %>%
  ungroup() %>%
  select(-dsim_dup) %>%
  distinct() 
  



plot_all_chrom <- 
  ggplot(all_chrom, aes(x = x_axis, y = y_axis,fill=dup_chrom)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_x_discrete(limits=factor(1:47)) +
  theme_bw() +
  ylab('Percentage of \n All Duplicates') +
  xlab('') +
  scale_y_continuous(breaks = c(0, .25, .50, .75, 1.00), labels = c("0", "25", "50", "75",'100')) +
  labs(fill = "")





dup_w_ortho_chrom <- plot_freq_dist_colored %>%
  group_by(dup) %>%
  mutate(x_axis = n_distinct(fly)) %>%
  ungroup() %>%
  distinct(dup,dsim_dup,dup_chrom, .keep_all = T) %>%
  group_by(x_axis,dsim_dup,dup_chrom) %>%
  mutate(y_axis = n()) %>%
  ungroup() %>%
  distinct(dsim_dup,dup_chrom,x_axis,y_axis) %>%
  filter(dsim_dup=='Dsim_Ortholog_Duplicated')

plot_dup_w_ortho_chrom <- 
  ggplot(dup_w_ortho_chrom, aes(x = x_axis, y = y_axis,fill=dup_chrom)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_x_discrete(limits=factor(1:47)) +
  theme_bw() +
  ylab('Percentage of \n Duplicates with Orthologs') +
  xlab('') +
  scale_y_continuous(breaks = c(0, .25, .50, .75, 1.00), labels = c("0", "25", "50", "75",'100')) +
  scale_fill_manual(values = c(colorRampPalette(c("#99CCFF", "royalblue4"))(5)), name = "")





dup_NO_ortho_chrom <- plot_freq_dist_colored %>%
  group_by(dup) %>%
  mutate(x_axis = n_distinct(fly)) %>%
  ungroup() %>%
  distinct(dup,dsim_dup,dup_chrom, .keep_all = T) %>%
  group_by(x_axis,dsim_dup,dup_chrom) %>%
  mutate(y_axis = n()) %>%
  ungroup() %>%
  distinct(dsim_dup,dup_chrom,x_axis,y_axis) %>%
  filter(dsim_dup=='No_Duplicated_Ortholog')


plot_dup_NO_ortho_chrom <- 
  ggplot(dup_NO_ortho_chrom, aes(x = x_axis, y = y_axis,fill=dup_chrom)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_x_discrete(limits=factor(1:47)) +
  theme_bw() +
  ylab('Percentage of \n Duplicates without Orthologs') +
  xlab('Number of Flies') +
  scale_y_continuous(breaks = c(0, .25, .50, .75, 1.00), labels = c("0", "25", "50", "75",'100')) +
  scale_fill_manual(values = c(colorRampPalette(c("plum1", "hotpink4"))(5)), name = "")


gg <- ggarrange(freq_dist,perc,plot_dup_w_ortho_chrom,plot_dup_NO_ortho_chrom,
                nrow=4,ncol=1,align='v',labels = c('A','B','C','D'))


ggsave("./Plots/dsim_freq_dist_barplots.jpg", plot = gg, width = 13, height = 12)


gg <- ggarrange(freq_dist,perc + xlab('Number of Flies'),nrow=2,ncol=1,align='v',labels = c('A','B'))
ggsave("./Plots/dsim_freq_perc_barplots.jpg", plot = gg, width = 13, height = 6)


##################################################

fly_counts <- fly_counts %>%
  select(-dsim_dup) %>%
  group_by(x_axis) %>%
  mutate(y_axis=sum(y_axis)) %>%
  ungroup() %>%
  distinct()


freq <- ggplot(fly_counts, aes(x = x_axis, y = y_axis)) + 
  geom_bar(position = "stack", stat='identity') +
  scale_x_discrete(limits = factor(1:47)) +
  theme_bw() +
  ylab('Frequency') +
  xlab('') 

gg <- ggarrange(freq,plot_all_chrom,nrow=2,ncol=1,align='v',labels=c('A','B'),
          common.legend = T, legend = "right")

ggsave("./Plots/freq_dist_chrom_barplots.jpg", plot = gg, width = 13, height = 6)




