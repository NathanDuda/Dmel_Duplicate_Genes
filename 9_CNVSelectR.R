

#install_github("peterbchi/CNVSelectR")
library(CNVSelectR)



CNVSelect_test_custom <- function(input_file, fasta_file){
  output <- list()
  
  ####
  ### inputs <- read_DAF(input_file)
  DAF <- read.csv(input_file, sep='')
  first.col <- dim(DAF)[1]
  output <- list()
  output$N <- as.numeric(DAF[1,2])   # N
  output$ploidy <- as.numeric(DAF[2,2])   # ploidy
  output$freqs <- as.numeric(DAF[5:first.col,1][1])  # frequencies
  inputs <- output
  ####
  
  ####
  ###dS <- get_dS(fasta_file)
  x <- read.alignment(fasta_file,format="fasta")
  n.seq <- length(x$seq)
  
  # Check if all sequences are of lengths that are multiples of 3
  # and also for internal stop codons
  for(i in 1:n.seq){
    tmp <- unlist(strsplit(x$seq[[i]], ""))
    if(length(tmp)/3 != round(length(tmp)/3, 0)){
      stop(paste("All sequences must be of lengths that are a multiple of 3; Sequence", i, "has", length(tmp), "sites"))
    }
    codon.starts <- seq(1, length(tmp), by=3)
    stop.codons <- c("taa", "tag", "tga", "TAA", "TAG", "TGA")
    for(j in codon.starts){
      codon <- paste(tmp[c(j, j+1, j+2)], collapse="")
      if(is.element(codon, stop.codons)){
        stop(paste("Sequences must not have internal stop codons; Sequence ", i, " has the following stop codon that starts at nucleotide position ", j, ": ", codon, sep=""))
      }
    }
  }
  
  dS <- NULL
  col_names <- NULL
  if(n.seq >= 2){ ################################################### only n.seq > 2
    seqs <- seq(1,n.seq, by=2)
    for(i in seqs){
      pair <- x
      pair$seq <- pair$seq[c(i, (i+1))]
      pair$nam <- pair$nam[c(i, (i+1))]
      pair$nb <- 2
      dS <- c(dS, kaks(pair)$ks)
      col_names <- c(col_names, paste(pair$nam[1], pair$nam[2], sep = "/"))
    }
  }
  dS <- matrix(dS, nrow=1)
  colnames(dS) <- col_names
  
  ####
  
  
  
  
  
  ##############
  
  #if(length(inputs$freqs) != length(dS)){
  #  stop("number of sequences provided must be twice the number of frequencies")
  #}
  
  output$freqs <- rep(inputs$freqs, length(dS))
  output$dS <- dS
  
  if(any(dS == 0)){return(0)} 
  
  for(i in 1:length(dS)){
    up <- 1e-3
    t <- (dS[i]*35)/up
    #if(!is.na(t)){
    out1 <- Genedupdip_neutralgenerator(inputs$N, up)
    out2 <- mexpv_dipneut(t=t, A=t(out1[[1]]), v=out1[[2]], N=inputs$N, Pos=out1[[3]])
    
    output$crit_lower[i] <- out2[[2]][length(out2[[2]])]
    output$crit_upper[i] <- out2[[3]][length(out2[[3]])]
    
    lower.tail <- out2[[5]][inputs$freqs[1]*2*inputs$N]
    if(lower.tail < 0.5){
      output$p_val[i] <- lower.tail * 2}
    else{
      output$p_val[i] <- (1-lower.tail) * 2}
  }
  return(output)

}




file_name <- 'group_2.fa'

nrow = 0 
output_df <- data.frame(
  group = character(),
  ds = character(),
  pval = character()
)


for (file_name in list.files('./CNVSelectR/Combined_Codon_Alignments/')) {
  nrow = nrow + 1
  
  fasta_file <- paste0("./CNVSelectR/Combined_Codon_Alignments/",file_name)
  freq <- sum(grepl("^>", readLines(fasta_file))) / 2
  
  nrep <- freq 
  
  add_frequencies <- rep(freq/47, nrep)
  add_NAs <- rep(NA, nrep + 1)
  
  input_file_table <- data.frame(
    A = c("Ne", "ploidy", "full/approximate", "frequency", add_frequencies),
    B = c(47, 2, 'full', add_NAs))
  write.table(input_file_table,file='./CNVSelectR/input_file.tsv')
  input_file <- './CNVSelectR/input_file.tsv'
  
  out <- CNVSelect_test_custom(input_file,fasta_file)
  

  if(any(out!=0) | typeof(out) == 'list'){
    output_df <- rbind(output_df, c(gsub("[^0-9]", "", file_name), c(out$dS), c(out$p_val)))
  } 

  colnames(output_df) <- c('group','ds','pval')
}


t <- data.frame()






