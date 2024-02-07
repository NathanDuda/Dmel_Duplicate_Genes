

#install_github("peterbchi/CNVSelectR")
library(CNVSelectR)

mexpv_dipneut_CUSTOM <- function(t, A, v,N,Pos, tol=1e-7, m=NULL){
  #TLS Have to modify number of arguments compared to standard expokit implementation to
  #account for extra variables N,Pos.
  n <- dim(A)[1]
  if(is.null(m)){
    m <- min(n,30)
  }
  
  anorm <- norm(A,'I')
  mxrej <- 10;  btol  <- 1.0e-7
  gamma <- 0.9; delta <- 1.2
  mb    <- m; t_out   <- abs(t)
  istep <- 0; t_new   <- 0
  t_now <- 0; s_error <- 0
  rndoff <-  anorm*2.2204e-16
  
  k1 <- 2; xm <- 1/m; normv <- sqrt(sum(v)^2); beta <- normv
  fact <- (((m+1)/exp(1))^(m+1))*sqrt(2*pi*(m+1))
  t_new <- (1/anorm)*((fact*tol)/(4*beta*anorm))^xm
  s <- 10^(floor(log10(t_new))-1); t_new <- ceiling(t_new/s)*s
  sgn <- sign(t); istep <- 0
  tlist <- NULL
  w <- v
  hump <- normv
  count1 <- 0
  
  # I think I need to initialize these since they are vectors
  dupgenomesavgMat <- NA
  PIlowerMat <- NA
  PIupperMat <- NA
  
  while (t_now < t_out){
    ####################################################################
    #TLS Main modifications come in here.
    tlist <- c(tlist, t_now)
    count1 <- count1+1
    # calculate proportion duplicated conditional on wt not fixing
    pn <- rep(0,2*N)#pn(i) stores conditional probably of haveing i haploid genomes with a duplicate at time t_now (each individual has two haploid genomes)
    if (w[1] > 0){
      w <- w/(1-w[1])
      w[1] <- 0
    }
    pnsummer <- 0
    pnflag <- 0
    for (i in 1:(N+1)){
      for (j in 1:(N+2-i)){
        if (i!=1 || j != 1){
          pnsummer <- pnsummer+w[Pos[i,j]]
          pn[2*(i-1)+(j-1)] <- pn[2*(i-1)+(j-1)] + w[Pos[i,j]]
          if (pnsummer == 1){
            pnflag <- 1
            break
          }
        }
        if (pnflag ==1){
          break
        }
      }
      if (pnflag == 1){
        break
      }
    }
    dupgenomesavg <- pn%*%(1:(2*N))/(2*N) #Stores the expected number of haploid genomes with a duplicate at time t_now
    
    #Now calculate the 95# prediction interval on the proportion of
    #duplicate haplotypes
    
    cumsumpn <- cumsum(pn)
    PIlower <- min(which(cumsumpn>=0.025,1))/(2*N)  # something is wrong here with larger N
    PIupper <- min(which(cumsumpn>=0.975,1))/(2*N)
    dupgenomesavgMat[count1] <- dupgenomesavg
    PIlowerMat[count1] <- PIlower
    PIupperMat[count1] <- PIupper
    
    
    #################################################
    #TLS - Now back to standard mexpv as implemented in expokit.
    
    istep <- istep + 1
    t_step <- min( t_out-t_now,t_new )
    V <- matrix(0, nrow=n,ncol=(m+1))
    H <- matrix(0, nrow=(m+2),ncol=(m+2))
    
    V[,1] <- as.double((1/beta)*as.matrix(w)) # this line seems to cause problems with bigger N
    for (j in 1:m){
      p <- A%*%V[,j]
      for (i in 1:j){
        H[i,j] <- as.double((V[,i]) %*% p)
        p <- p-H[i,j]*V[,i]
      }
      s <- norm(p, "f")   # Matlab defaults to Frobenius norm, make sure this is ok though since p is a vector
      if (s < btol){
        k1 <- 0
        mb <- j
        t_step <- t_out-t_now
        break
      }
      H[j+1,j] <- s
      V[,j+1] <- as.double((1/s)*p)
    }
    if (k1 != 0){
      H[m+2,m+1] <- 1
      avnorm <- norm(A%*%V[,m+1], "f")  # 4/29/21 I think this is right but check closer
    }
    ireject <- 0
    while (ireject <= mxrej){
      mx <- mb + k1
      F <- expm(sgn*t_step*H[1:mx,1:mx])
      if (k1 == 0){
        err_loc <- btol
        break
      } else {
        phi1 <- abs( beta*F[m+1,1] )
        phi2 <- abs( beta*F[m+2,1] * avnorm )
        if (phi1 > 10*phi2){
          err_loc <- phi2
          xm <- 1/m
        } else if (phi1 > phi2){
          err_loc <- (phi1*phi2)/(phi1-phi2)
          xm <- 1/m
        } else {
          err_loc <- phi1
          xm <- 1/(m-1)
        }
      }
      if (err_loc <= delta * t_step*tol){
        break
      } else {
        t_step <- gamma * t_step * (t_step*tol/err_loc)^xm
        s <- 10^(floor(log10(t_step))-1)
        t_step <- ceiling(t_step/s) * s
        if (ireject == mxrej){
          #error('The requested tolerance is too high.')
          return(NA)
        }
        ireject <- ireject + 1
      }
    }
    mx <- mb + max( c(0,(k1-1) ))
    w <- V[,(1:mx)]%*%(beta*F[1:mx,1])   # check this 4/29/21
    beta <- norm( w , "f")
    hump <- max(hump,beta)
    
    ineg <- 0
    for (i in 1:n){
      if (w[i] < 0){
        w[i] <- 0
        ineg <- ineg + 1
      }
    }
    wnorm <- norm(w)
    if (ineg > 0){
      w <- (1/wnorm)*w
    }
    roundoff <- abs(1-wnorm)/n   # this was 1.0d0, no idea what that means
    
    t_now <- t_now + t_step
    t_new <- gamma * t_step * (t_step*tol/err_loc)^xm
    s <- 10^(floor(log10(t_new))-1)
    t_new <- ceiling(t_new/s) * s
    
    err_loc <- max(err_loc,roundoff)
    err_loc <- max(err_loc,rndoff)
    s_error <- s_error + err_loc
  }
  err <- s_error
  hump <- hump / normv
  
  return(list(dupgenomesavgMat,PIlowerMat,PIupperMat,tlist, cumsumpn ,w))
}



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
  
  if(any(dS == 0)){return(NA)} 
  
  for(i in 1:length(dS)){
    up <- 1e-3
    t <- (dS[i]*35)/up
    #if(!is.na(t)){
    out1 <- Genedupdip_neutralgenerator(inputs$N, up)
    out2 <- mexpv_dipneut_CUSTOM(t=t, A=t(out1[[1]]), v=out1[[2]], N=inputs$N, Pos=out1[[3]])
    
    if (any(is.na(out2))){return(NA)}
    
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




file_name <- 'group_3664.fa'

group = 0 
output_df <- data.frame()

for (file_name in list.files('./CNVSelectR/Blastp_Dups_Combined_Codon_Alignments/')) {
  group = group + 1
  
  fasta_file <- paste0("./CNVSelectR/Blastp_Dups_Combined_Codon_Alignments/",file_name)
  freq <- sum(grepl("^>", readLines(fasta_file))) / 2
  
  nrep <- freq 
  
  add_frequencies <- rep(freq/47, nrep)
  add_NAs <- rep(NA, nrep + 1)
  
  input_file_table <- data.frame(
    A = c("Ne", "ploidy", "full/approximate", "frequency", add_frequencies),
    B = c(47, 2, 'full', add_NAs))
  write.table(input_file_table,file='./CNVSelectR/input_file.tsv')
  input_file <- './CNVSelectR/input_file.tsv'
  
  
  
  tryCatch({# error whenever internal stop codons 
    out <- CNVSelect_test_custom(input_file, fasta_file)
    
    if (!any(is.na(out))) {
      for (i in 1:length(out[['dS']])) {
        ds <- out[['dS']][i]
        name <- colnames(out[['dS']])[i]  
        pval <- out[['p_val']][i]
        
        output_df <- rbind(output_df, data.frame(name, ds, pval, group))
      }
    }
  }, error = function(e) {
    # If an error occurs, print a message and continue with the loop
    print(paste("Error occurred for file:", file_name))
  })
  
  print(group)
}


write.table(output_df, 'CNVSelectR_Output.tsv')

