









source("startup.R")

# import expression 
expr <- read.csv("C:/Users/17735/Downloads/Dmel_Duplicate_Genes/Expression_Values_for_Ranges.tsv", sep="")

# remove the two rows without expression values 
expr <- na.omit(expr)

# get the average of the expression values for the same tissue in the same fly 
expr <- expr %>%
  group_by(fly, chrom, start, end, tissue) %>%
  mutate(expression_value = mean(expression_value)) %>%
  distinct(fly, chrom, start, end, tissue, expression_value)

# import duplicates
dups <- read.csv("./mRNA_Duplicate_Sequences.tsv", sep="")
dups_orig <- dups

dups <- dups_orig
# merge hits with their expression values
colnames(expr) <- c('fly','hit_chrom','hit_start_on_chrom','hit_end_on_chrom','tissue','hit_expression_value')
dups <- left_join(dups,expr,by=c('fly','hit_chrom','hit_start_on_chrom','hit_end_on_chrom'))

# merge gns with their expression values
colnames(expr) <- c('fly','gn_chrom','gn_start','gn_end','tissue','gn_expression_value')
dups <- left_join(dups,expr,by=c('fly','gn_chrom','gn_start','gn_end','tissue'))

# set the expression values below 1 to 0
dups$gn_expression_value[dups$gn_expression_value < 1] <- 0
dups$hit_expression_value[dups$hit_expression_value < 1] <- 0

# keep rows that have an expression value for the gn or the hit or both 
dups <- subset(dups, !is.na(gn_expression_value) | !is.na(hit_expression_value))

# format dataframe so that there are columns for each of the different expression values 
dups <- dups %>%
  pivot_wider(names_from='tissue',values_from=c('gn_expression_value','hit_expression_value'))
  #select(-gn_expression_value_NA)

# calculate differences in expression between gns and their hits
dups$change_in_exp_ovary <- dups$hit_expression_value_Ovary - dups$gn_expression_value_Ovary
dups$change_in_exp_head <- dups$hit_expression_value_Head - dups$gn_expression_value_Head
dups$change_in_exp_gut <- dups$hit_expression_value_Gut - dups$gn_expression_value_Gut
dups$change_in_exp_wb <- dups$`hit_expression_value_Whole-body` - dups$`gn_expression_value_Whole-body`


# make boxplot of expression of duplicate versus its hit(s)
plot_expr_diff_boxplot <- dups[,c(1,16,23:34)]

plot_expr_diff_boxplot_ovary <- plot_expr_diff_boxplot[,c(1,2,4,8)]
plot_expr_diff_boxplot_ovary <- na.omit(plot_expr_diff_boxplot_ovary)
plot_expr_diff_boxplot_ovary <- reshape2::melt(plot_expr_diff_boxplot_ovary,id.vars=c(1,2))
plot_expr_diff_boxplot_ovary <- separate(plot_expr_diff_boxplot_ovary, variable, into = c("gn_or_hit", "tissue"), sep = "_expression_value_")

plot_expr_diff_boxplot_head <- plot_expr_diff_boxplot[,c(1,2,6,10)]
plot_expr_diff_boxplot_head <- na.omit(plot_expr_diff_boxplot_head)
plot_expr_diff_boxplot_head <- reshape2::melt(plot_expr_diff_boxplot_head,id.vars=c(1,2))
plot_expr_diff_boxplot_head <- separate(plot_expr_diff_boxplot_head, variable, into = c("gn_or_hit", "tissue"), sep = "_expression_value_")

plot_expr_diff_boxplot_gut <- plot_expr_diff_boxplot[,c(1,2,5,9)]
plot_expr_diff_boxplot_gut <- na.omit(plot_expr_diff_boxplot_gut)
plot_expr_diff_boxplot_gut <- reshape2::melt(plot_expr_diff_boxplot_gut,id.vars=c(1,2))
plot_expr_diff_boxplot_gut <- separate(plot_expr_diff_boxplot_gut, variable, into = c("gn_or_hit", "tissue"), sep = "_expression_value_")

plot_expr_diff_boxplot_wb <- plot_expr_diff_boxplot[,c(1,2,3,7)]
plot_expr_diff_boxplot_wb <- na.omit(plot_expr_diff_boxplot_wb)
plot_expr_diff_boxplot_wb <- reshape2::melt(plot_expr_diff_boxplot_wb,id.vars=c(1,2))
plot_expr_diff_boxplot_wb <- separate(plot_expr_diff_boxplot_wb, variable, into = c("gn_or_hit", "tissue"), sep = "_expression_value_")


#plot_expr_diff_boxplot_ovary <- na.omit(plot_expr_diff_boxplot[plot_expr_diff_boxplot$tissue == 'Ovary',])
#plot_expr_diff_boxplot_head <- na.omit(plot_expr_diff_boxplot[plot_expr_diff_boxplot$tissue == 'Head',])
#plot_expr_diff_boxplot_gut <- na.omit(plot_expr_diff_boxplot[plot_expr_diff_boxplot$tissue == 'Gut',])
#plot_expr_diff_boxplot_wb <- na.omit(plot_expr_diff_boxplot[plot_expr_diff_boxplot$tissue == 'Whole-body',])

ggarrange(
  ggpaired(plot_expr_diff_boxplot_ovary, x = "gn_or_hit", y = "value", color = "gn_or_hit", line.color = "gray", line.size = 0.1) + scale_y_log10(limits=c(0.01,100)) + ggtitle('ovary'),
  ggpaired(plot_expr_diff_boxplot_head, x = "gn_or_hit", y = "value", color = "gn_or_hit", line.color = "gray", line.size = 0.1) + scale_y_log10(limits=c(0.01,100)) + ggtitle('head'),
  ggpaired(plot_expr_diff_boxplot_gut, x = "gn_or_hit", y = "value", color = "gn_or_hit", line.color = "gray", line.size = 0.1) + scale_y_log10(limits=c(0.01,100)) + ggtitle('gut'),
  ggpaired(plot_expr_diff_boxplot_wb, x = "gn_or_hit", y = "value", color = "gn_or_hit", line.color = "gray", line.size = 0.1) + scale_y_log10(limits=c(0.01,100)) + ggtitle('whole body'),
  nrow = 1,
  ncol = 4
)


# make density plot of the differences in expression 
plot_expr_diff <- dups[,c(1,16,23:34)]

ggarrange(
  ggplot(plot_expr_diff, aes(x = change_in_exp_ovary)) + geom_density(bw = 0.1) + theme_bw() + xlim(-10,10),
  ggplot(plot_expr_diff, aes(x = change_in_exp_head)) + geom_density(bw = 0.1) + theme_bw() + xlim(-10,10),
  ggplot(plot_expr_diff, aes(x = change_in_exp_gut)) + geom_density(bw = 0.1) + theme_bw() + xlim(-10,10),
  ggplot(plot_expr_diff, aes(x = change_in_exp_wb)) + geom_density(bw = 0.1) + theme_bw() + xlim(-10,10),
  nrow = 1,
  ncol = 4
)


# make scatter plot of the change in expression compared to gn expression
ggarrange(
  ggplot(plot_expr_diff,aes(y=gn_expression_value_Ovary,x=change_in_exp_ovary)) + geom_point(size = 0.9) + theme_bw() + xlim(-20,30) + ylim(1,30),
  ggplot(plot_expr_diff,aes(y=gn_expression_value_Head,x=change_in_exp_head)) + geom_point(size = 0.9) + theme_bw() + xlim(-20,30) + ylim(1,30),
  ggplot(plot_expr_diff,aes(y=gn_expression_value_Gut,x=change_in_exp_gut)) + geom_point(size = 0.9) + theme_bw() + xlim(-20,30) + ylim(1,30),
  ggplot(plot_expr_diff,aes(y=`gn_expression_value_Whole-body`,x=change_in_exp_wb)) + geom_point(size = 0.9) + theme_bw() + xlim(-20,30) + ylim(1,30),
  nrow = 1,
  ncol = 4
)

# sequence divergence versus expression divergence 
blast_output <- read.csv("./Blast_Outputs/mRNA_raw_Duplicates_from_Blast.tsv", sep="")

blast_output <- blast_output[,-c(5:7,12)]
colnames(blast_output) <- c('gn_id','fly','hit_chrom','perc_identity',
                            'hit_start_on_chrom','hit_end_on_chrom',
                            'evalue','bitscore')

dups <- merge(dups,blast_output,by=c('gn_id','fly','hit_chrom','hit_start_on_chrom','hit_end_on_chrom'))


ggarrange(
  ggplot(dups,aes(y=bitscore,x=change_in_exp_ovary)) + geom_point(size = 0.9) + theme_bw() + xlim(-20,20) + ylim(0,22500),
  ggplot(dups,aes(y=bitscore,x=change_in_exp_head)) + geom_point(size = 0.9) + theme_bw() + xlim(-20,20) + ylim(0,22500),
  ggplot(dups,aes(y=bitscore,x=change_in_exp_gut)) + geom_point(size = 0.9) + theme_bw() + xlim(-20,20) + ylim(0,22500),
  ggplot(dups,aes(y=bitscore,x=change_in_exp_wb)) + geom_point(size = 0.9) + theme_bw() + xlim(-20,20) + ylim(0,22500),
  nrow = 1,
  ncol = 4
)

ggplot(dups,aes(y=`hit_expression_value_Whole-body`,x=`gn_expression_value_Whole-body`)) + 
  geom_point(size = 1.9) +
  theme_bw()

ggplot(dups,aes(y=perc_identity.x,x=change_in_exp_wb)) + 
  geom_point(size = 0.1) 

ggplot(dups,aes(y=evalue,x=change_in_exp_wb)) + 
  geom_point(size = 0.1) +
  scale_y_log10()

# base pair distance versus expression and sequence divergence

trash <- dups

trash <- dups %>%
  mutate(distance = case_when(hit_chrom == gn_chrom ~ abs(gn_start - hit_start_on_chrom),
                              hit_chrom != gn_chrom ~ NA))



# duplicates on different versus same chromosomes sequence and expression divergence 
dups <- dups %>%
  mutate(same_chrom = case_when(hit_chrom == gn_chrom ~ 'same_chrom',
                                hit_chrom != gn_chrom ~ 'diff_chrom')) %>%
  mutate(distance = case_when(hit_chrom == gn_chrom ~ abs(gn_start - hit_start_on_chrom),
                              hit_chrom != gn_chrom ~ NA))


ggplot(dups, aes(x=same_chrom, y=bitscore)) +
  geom_boxplot() +
  scale_y_log10()

ggplot(dups,aes(y=bitscore,x=distance)) + 
  geom_point(size = 0.9)

ggplot(dups, aes(x=same_chrom, y=change_in_exp_wb)) +
  geom_boxplot()
  #scale_y_log10() # gets rid of 0s and negatives

ggplot(dups, aes(x=bitscore, y=change_in_exp_wb, color = same_chrom)) +
  geom_point(size = 0.9) +
  scale_x_log10() +
  scale_y_log10()


#

same_chrom_dups <- dups[dups$same_chrom == 'same_chrom',]

# change in expression by chromosome (for duplicates on same chromosome)
ggplot(same_chrom_dups,aes(x=gn_chrom,y=change_in_exp_wb)) + 
  geom_boxplot() +
  geom_point() +
  ylim(-15,20) +
  theme_bw()

# heatmap of expression based on chromosome locations
heat_chrom_exp <- dups %>% 
  group_by(gn_chrom,hit_chrom) %>% 
  summarize(mean(change_in_exp_wb))

ggplot(dups,aes(x=hit_chrom,y=gn_chrom))



# frequency distribution
plot_freq_dist <- dups[,c(1,2)]

# keep only pairs 
# if gene is duplicated (2 or more copies) it is counted 
plot_freq_dist <- unique(plot_freq_dist)

plot_freq_dist <- as.data.frame(table(plot_freq_dist$gn_id))
plot_freq_dist <- as.data.frame(table(plot_freq_dist$Freq))

library(ggplot2)
ggplot(plot_freq_dist,aes(x=Var1,y=Freq)) +
  geom_bar(stat='identity') +
  scale_x_discrete(limits=c(1:50)) +
  xlab('number of flies duplicate is present in') +
  ylab('freq') +
  theme_bw()

# 400 times, a duplicated gene is present in only one fly





# hits to area with high expression 
# (surrounding expression not considered)




# investiage the few that have 0 expression in gn then have hit expressed 


# diff chrom


trash <- table(dups$gn_chrom,dups$hit_chrom)




