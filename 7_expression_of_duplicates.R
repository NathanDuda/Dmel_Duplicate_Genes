# Title: 7_expression_of_duplicates.R
# Author: Nathan Duda
# Date: 10/8/2023
# Purpose: Extract expression values for each gene and each hit for every tissue I have expression values for. 


source("startup.R")

# read in duplicates 
dups <- read.csv("./Duplicate_Sequences.tsv", sep="")
dups_orig <- dups

# import dataframe with fly names, expression file names, and tissues
individual_expression_file <- read.delim("./Supplementary_Tables_CopyPaste_RNA.txt", header=FALSE)
individual_expression_file <- individual_expression_file[individual_expression_file$V2 == 'Control'| individual_expression_file$V2 == '-', ]
individual_expression_file <- individual_expression_file[,c(1,4,3)]
colnames(individual_expression_file) <- c('fly','expression_file','tissue')

# get only the fly, chromosome, and coordinate, columns of both the gene and hit to get expression values for these ranges
dups_1 <- dups[,c(1,4,6,7)]
dups_2 <- dups[,c(1,4,11,12)]

colnames(dups_1) <- c('fly','chrom','start','end')
colnames(dups_2) <- c('fly','chrom','start','end')

dups <- rbind(dups_1,dups_2)


# read in all 50 expression files as objects 
read_and_assign <- function(file, name) {
  df <- import(paste0('./RNA_Seq/',file,'.bw',sep=''), format = "BigWig")
  assign(name, df, envir = .GlobalEnv)
}

lapply(1:nrow(individual_expression_file), function(i) {
  read_and_assign(individual_expression_file$expression_file[i], individual_expression_file$expression_file[i])
})


# keep only the flies that expression files exist for 
dups <- dups[dups$fly %in% individual_expression_file$fly, ]

# add expression file and tissue information columns to each range 
dups <- merge(dups,individual_expression_file,by=c('fly'))

# make empty column for expression values to be entered into 
dups$expression_value <- NA

# run for loop that fills in the expression values 
for (row_num in 1:nrow(dups)){
  
  row <- dups[row_num,]

  expr <- as.data.frame(get(row$expression_file))
  
  # end peaking into exp
  expression_subset_end <- expr[ (expr$start > dups$start[row_num]) & (expr$end > dups$end[row_num]) & (expr$start < dups$end[row_num]),]
  end_of_hit_peaking_into_exp <- expression_subset_end$score / (dups$end[row_num] - expression_subset_end$start)
  
  # expression inside gn
  expression_subset_mid <- expr[expr$start > dups$start[row_num] & expr$end < dups$end[row_num],]
  expression_subset_mid$score <- expression_subset_mid$score / expression_subset_mid$width
  expression_subset_mid <- expression_subset_mid$score
  
  # gn inside expression (most likely not expressed at all)
  expression_subset <- expr[expr$start < dups$start[row_num] & expr$end > dups$end[row_num],]
  expression_subset$score <- expression_subset$score / expression_subset$width
  expression_subset <- expression_subset$score
  
  # start of hit is in exp
  expression_subset_start <- expr[(expr$start < dups$start[row_num]) & (expr$end < dups$end[row_num]) & (expr$end > dups$start[row_num]),]
  start_of_hit_in_exp <- expression_subset_start$score / (expression_subset_start$end - dups$start[row_num])
  
  gene_exp_sections <- c(end_of_hit_peaking_into_exp,expression_subset_mid,expression_subset,start_of_hit_in_exp)
  dups[row_num,'expression_value'] <- mean(gene_exp_sections)


  print(paste('Progress: ', round(((row_num / nrow(dups)) * 100),1) , ' %'))
  
}

# write to file 
write.table(dups,'Expression_Values_for_Ranges.tsv')

