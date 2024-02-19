

# run on Temple's Compute cluster: 


# altering needleman function from ftrCOOL to speed up 
needleman <- function(seq1, seq2, max_length_diff, gap=-1, mismatch=-1, match=1){
    
  # Initialize col and rownames for matrices
  len1 = nchar(seq1) # Save number of chars in each sequence
  len2 = nchar(seq2) # Save number of chars in each sequence
  seq1 = unlist(strsplit(seq1, split = "")) # convert seq to character vector
  seq2 = unlist(strsplit(seq2, split = "")) # convert seq to character vector
  
  # Initialize matrix M (for scores)
  M = matrix(0, nrow = len1 + 1, ncol = len2 + 1) # Initialize matrix
  rownames(M) = c("-", seq1) # assign seq chars to matrix names
  colnames(M) = c("-", seq2) # assign seq chars to matrix names
  M[1, ] = cumsum(c(0, rep(gap, len2))) # Fill 1st row with gap penalites
  M[, 1] = cumsum(c(0, rep(gap, len1))) # Fill 1st col with gap penalites
  
  # Initialize matrix D (for directions)
  D = matrix(0, nrow = len1 + 1, ncol = len2 + 1) # Initialize matrix
  rownames(D) = c("-", seq1) # assign seq chars to matrix names
  colnames(D) = c("-", seq2) # assign seq chars to matrix names
  D[1, ] = rep("hor") # Fill 1st row with "hor" for horizontal moves
  D[, 1] = rep("ver") # Fill 1st col with "ver" for vertical moves
  type = c("dia", "hor", "ver") # Lookup vector
  
  n_diff = 0
  
  # Compute scores and save moves
  for (i in 2:(len1 + 1)){# for every (initially zero) row
    for (j in 2:(len2 + 1)){# for every (initially zero) col
      hor = M[i, j - 1] + gap # horizontal move = gap for seq1
      ver = M[i - 1, j] + gap # vertical move = gap for seq2
      dia = ifelse(rownames(M)[i] == colnames(M)[j], # diagonal = ifelse(chars equal, match, mismatch)
                   M[i - 1, j - 1] + match,
                   M[i - 1, j - 1] + mismatch)
      M[i, j] = max(dia, hor, ver) # Save current (best) score in M
      D[i, j] = type[which.max(c(dia, hor, ver))] # Save direction of move in D
      
      if (D[i, j] == 'hor' | D[i, j] == 'ver'){
        n_diff <- n_diff + 1
        if (n_diff > max_length_diff){
          return('diff')
        }
      }
      
    }
  }
  return(M[nrow(M), ncol(M)])
}



# run global alignment in parallel:
library(foreach)
library(doParallel)

cl <- makeCluster(detectCores())
registerDoParallel(cl)

fly_name_list <- read.table("/home/tun37257/workdir/Dmel_Pop/Global_Alignment/Individuals_List_47.txt", quote="\"", comment.char="")

global_align_fly <- function(fly) {  
  
  library(dplyr)
  
  all_annotations <- read.csv("/home/tun37257/workdir/Dmel_Pop/Global_Alignment/Annotations_Gene.tsv", sep="")
  all_prots <- all_annotations

  # get proteins for fly
  prots <- all_prots[all_prots$fly == fly, c('gene_group', 'chrom', 'nchar', 'prot')]
  colnames(prots) <- c('gene_group', 'chrom', 'prot_length', 'prot')

  # keep only proteins with at least 100 aa
  prots <- prots[prots$nchar >= 100,]
  
  # get all possible pairs 
  all_pairs <- merge(prots,prots, by = NULL)
  all_pairs <- all_pairs[all_pairs$gene_group.x != all_pairs$gene_group.y,]
  
  # keep only protein pairs that have a length difference smaller than 1% of the shorter protein length 
  all_pairs$length_diff <- abs(all_pairs$prot_length.x - all_pairs$prot_length.y)
  all_pairs <- all_pairs %>% mutate(shorter_length = case_when(prot_length.x > prot_length.y ~ prot_length.y,
                                                               prot_length.x <= prot_length.y ~ prot_length.x))
  all_pairs$max_length_diff <- floor(all_pairs$shorter_length * 0.01) 
  
  similar_length_pairs <- all_pairs %>% filter(length_diff <= max_length_diff)
  
  
  # keep only reciprocal hits 
  similar_length_pairs <- similar_length_pairs %>% 
    mutate(q_s = paste0(gene_group.x,chrom.x,'_',gene_group.y,chrom.y)) %>%
    mutate(s_q = paste0(gene_group.y,chrom.y,'_',gene_group.x,chrom.x)) %>%
    filter(s_q %in% q_s) 
  
  # keep only one of the reciprocal hits 
  similar_length_pairs <- similar_length_pairs %>%
    filter(!(q_s > s_q)) %>%
    select(-q_s,-s_q)
  
  # perform the global alignment 
  global_align_output <- similar_length_pairs %>% 
    rowwise() %>%
    mutate(score = needleman(prot.x, prot.y, max_length_diff, gap = -1, mismatch = -1, match = 1))

  global_align_output <- global_align_output[c('gene_group.x','gene_group.y','score')]
  
  write.table(global_align_output,file=paste0(fly,'_global_align.tsv'))
  
}


foreach(fly = fly_name_list[, 1]) %dopar% {
  global_align_fly(fly)
}

stopCluster(cl)

