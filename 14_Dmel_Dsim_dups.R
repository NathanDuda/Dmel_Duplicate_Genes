


source("startup.R")





# import duplicates for dsim 
dsim_dups <- read.csv("./Dsim/Duplicate_Proteins.tsv", sep="")

# import duplicates for dmel 
dmel_dups <- read.csv("./Duplicate_Proteins.tsv", sep="")
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
  
  print(row)

}

###########

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


# Creating the bar plot
ggplot(fly_counts, aes(fill = factor(dsim_dup,levels=c('No_Duplicated_Ortholog','Dsim_Ortholog_Duplicated')), 
                       x = x_axis, y = y_axis)) + 
  geom_bar(position = "stack", stat='identity') +
  scale_x_discrete(limits = factor(1:47)) +
  theme_bw() +
  ylab('Frequency') +
  xlab('') +
  scale_fill_manual(values = c("plum4", "steelblue3"), name = "") +
  guides(fill = guide_legend(title = "", override.aes = list(fill = c("plum4", "steelblue3"))))

##########

data <- data.frame(category = c("A", "B", "A", "B", "C", "C", "A"))

# Plot using ggplot2
ggplot(data, aes(x = category)) +
  geom_bar()


###########


freq_dist <- ggplot(plot_freq_dist_colored, aes(fill = factor(Var2,levels=c('No_Duplicated_Ortholog','Dsim_Ortholog_Duplicated')), y = Freq, x = Var1)) + 
  geom_bar(position = "stack", stat = "identity") +
  scale_x_discrete(limits = factor(1:47)) +
  theme_bw() +
  ylab('Frequency') +
  xlab('') +
  scale_fill_manual(values = c("plum4", "steelblue3"), name = "") +
  guides(fill = guide_legend(title = "", override.aes = list(fill = c("plum4", "steelblue3"))))

plot_freq_dist_percent <- plot_freq_dist_colored %>%
  group_by(Var1) %>%
  mutate(Percentage = Freq / sum(Freq) * 100)

perc <- ggplot(plot_freq_dist_percent, aes(x = Var1, y = Percentage, fill = factor(Var2,levels=c('No_Duplicated_Ortholog','Dsim_Ortholog_Duplicated')))) +
  geom_bar(stat = "identity", position = "stack") +
  scale_x_discrete(limits=factor(1:47)) +
  theme_bw() +
  xlab('') +
  ylab('Percentage of \n Duplicates ') +
  scale_fill_manual(values = c("plum4", "steelblue3"), name = "") +
  guides(fill = guide_legend(title = "", override.aes = list(fill = c("plum4", "steelblue3"))))

  
gg <- ggarrange(freq_dist,perc,nrow = 2, ncol = 1,
          labels = c('A','B'),
          align = 'v')

gg

ggsave("./Plots/dsim_freq_dist_barplot.jpg", plot = gg, width = 12, height = 6)


############################



plot_freq_dist_colored <- dups_orthos[c('dup','fly','dsim_dup','dup_chrom')]
#plot_freq_dist_colored$dsim_dup <- paste0(plot_freq_dist_colored$dsim_dup,'_',plot_freq_dist_colored$dup_chrom)
plot_freq_dist_colored <- as.data.frame(table(plot_freq_dist_colored$dup,
                                              plot_freq_dist_colored$dsim_dup,
                                              plot_freq_dist_colored$dup_chrom))
#plot_freq_dist_colored <- plot_freq_dist_colored[!plot_freq_dist_colored$Freq == 0,]
#plot_freq_dist_colored <- as.data.frame(table(plot_freq_dist_colored$Freq,plot_freq_dist_colored$Var2,plot_freq_dist_colored$Var3))





dup_w_ortho_chrom <- plot_freq_dist_colored[plot_freq_dist_colored$Var2 =='Dsim_Ortholog_Duplicated',]


dup_w_ortho_chrom <- dup_w_ortho_chrom %>%
  group_by(Var1) %>%
  mutate(Percentage = Freq / sum(Freq) * 100) %>%
  mutate(x_axis = sum(Freq)) %>%
  na.omit()
dup_w_ortho_chrom <- as.data.frame(table(dup_w_ortho_chrom$Freq,dup_w_ortho_chrom$Var3))
dup_w_ortho_chrom <- dup_w_ortho_chrom %>%
  group_by(Var1) %>%
  mutate(Percentage = Freq / sum(Freq) * 100) 

plot_dup_w_ortho_chrom <- 
  ggplot(t, aes(x = Var1, y = Percentage,fill=Var2)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_x_discrete(limits=factor(1:47)) +
  theme_bw() +
  ylab('Percentage of \n Duplicates with Orthologs') +
  xlab('') +
  scale_y_continuous(breaks = c(0, .25, .50, .75, 1.00), labels = c("0", "25", "50", "75",'100')) +
  scale_fill_manual(values = c('steelblue1','steelblue3','steelblue4','royalblue4','navy'), name = "")




plot_dup_NO_ortho_chrom <- plot_freq_dist_colored[plot_freq_dist_colored$Var2 == 'No_Duplicated_Ortholog',]

plot_dup_NO_ortho_chrom <- plot_dup_NO_ortho_chrom %>%
  group_by(Var1) %>%
  mutate(Percentage = Freq / sum(Freq) * 100) %>%
  mutate(x_axis = sum(Freq)) %>%
  ungroup()
plot_dup_NO_ortho_chrom <- as.data.frame(table(plot_dup_NO_ortho_chrom$Freq,plot_dup_NO_ortho_chrom$Var3))
plot_dup_NO_ortho_chrom <- plot_dup_NO_ortho_chrom %>%
  group_by(Var1) %>%
  mutate(Percentage = Freq / sum(Freq) * 100) 

plot_p <- ggplot(plot_dup_NO_ortho_chrom, aes(x = Var1, y = Percentage,fill=Var2)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_x_discrete(limits=factor(1:47)) +
  theme_bw() +
  ylab('Percentage of \n Duplicates without Orthologs') +
  xlab('Number of Flies') +
  scale_y_continuous(breaks = c(0, .25, .50, .75, 1.00), labels = c("0", "25", "50", "75",'100')) +
  scale_fill_manual(values = c('lightpink','hotpink','mediumorchid1','mediumorchid','mediumorchid4'), name = "")



gg <- ggarrange(freq_dist,perc,plot_dup_w_ortho_chrom,plot_p,nrow=4,ncol=1,align='v',
          labels = c('A','B','C','D'))


ggsave("./Plots/dsim_freq_dist_barplot.jpg", plot = gg, width = 12, height = 12)
ggsave("./Plots/dsim_freq_dist_barplot.jpg")
