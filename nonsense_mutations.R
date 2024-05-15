
source('./startup.R')


dups <- read.csv("./connected_dups.tsv", sep = "")
colnames(dups) <- c('gene_group', 'community')


seqs <- read.csv("./Annotations_Gene.tsv", sep = "")



dups <- merge(dups, seqs, by = 'gene_group')


# keep only duplicate pairs
t <- dups %>%
  mutate(community = paste0(community, '_', fly)) %>%
  group_by(community) %>%
  mutate(n_copies = n()) %>%
  filter(n_copies == 2)


# calculate dn, ds, and nonsense


calculate_mutations <- function(seq1, seq2) {
  # Function to calculate mutations
  nonsynonymous1 <- 0
  nonsynonymous2 <- 0
  synonymous <- 0
  nonsense1 <- 0
  nonsense2 <- 0
  
  nonsynonymous1_positions <- c()
  nonsynonymous2_positions <- c()
  synonymous_positions <- c()
  nonsense1_positions <- c()
  nonsense2_positions <- c()
  
  
  if (nchar(seq1) %% 3 != 0) {stop('length of seq1 is not a a multiple of 3')}
  if (nchar(seq2) %% 3 != 0) {stop('length of seq2 is not a a multiple of 3')}
  
  # Check length of sequences
  if (nchar(seq1) != nchar(seq2)) {
    stop("Sequences are not aligned properly")
  }
  
  for (i in 1:nchar(seq1)) {
    position = (i * 3) - 3
    codon1 <- substr(seq1, position + 1, position + 3)
    codon2 <- substr(seq2, position + 1, position + 3)
    
    # Check if codons are complete
    if (nchar(codon1) == 3 && nchar(codon2) == 3) {
      # Check for synonymous mutation
      if (codon1 != codon2) {
        if (translate_codon(codon1) == translate_codon(codon2)) {
          synonymous <- synonymous + 1
          synonymous_positions <- c(synonymous_positions, position)
        }
        
        stop_codons <- c("TAA", "TAG", "TGA")
        if (codon1 %in%  stop_codons) {
          nonsense1 <- nonsense1 + 1
          nonsense1_positions <- c(nonsense1_positions, position)
        }
        if (codon2 %in%  stop_codons) {
          nonsense2 <- nonsense2 + 1
          nonsense2_positions <- c(nonsense2_positions, position)
        }
        
        if ((translate_codon(codon1) != translate_codon(codon2)) & (!codon1 %in% stop_codons)) {
          nonsynonymous1 <- nonsynonymous1 + 1
          nonsynonymous1_positions <- c(nonsynonymous1_positions, position)
        }
        if ((translate_codon(codon1) != translate_codon(codon2)) & (!codon2 %in% stop_codons)) {
          nonsynonymous2 <- nonsynonymous2 + 1
          nonsynonymous2_positions <- c(nonsynonymous2_positions, position)
        }
        
          
        }
      }
    }
  
  
  return(list(nonsynonymous1 = nonsynonymous1, 
              nonsynonymous2 = nonsynonymous2,
              synonymous = synonymous, 
              nonsense1 = nonsense1,
              nonsense2 = nonsense2,
              
              nonsynonymous1_positions = nonsynonymous1_positions, 
              nonsynonymous2_positions = nonsynonymous2_positions,
              synonymous_positions = synonymous_positions, 
              nonsense1_positions = nonsense1_positions,
              nonsense2_positions = nonsense2_positions))
}


translate_codon <- function(codon) {
  # Function to translate codon to amino acid
  if (grepl("-", codon)) {return("-")}
  
  codon_table <- list(
    "TTT" = "F", "TTC" = "F", "TTA" = "L", "TTG" = "L",
    "TCT" = "S", "TCC" = "S", "TCA" = "S", "TCG" = "S",
    "TAT" = "Y", "TAC" = "Y", "TAA" = "*", "TAG" = "*",
    "TGT" = "C", "TGC" = "C", "TGA" = "*", "TGG" = "W",
    "CTT" = "L", "CTC" = "L", "CTA" = "L", "CTG" = "L",
    "CCT" = "P", "CCC" = "P", "CCA" = "P", "CCG" = "P",
    "CAT" = "H", "CAC" = "H", "CAA" = "Q", "CAG" = "Q",
    "CGT" = "R", "CGC" = "R", "CGA" = "R", "CGG" = "R",
    "ATT" = "I", "ATC" = "I", "ATA" = "I", "ATG" = "M",
    "ACT" = "T", "ACC" = "T", "ACA" = "T", "ACG" = "T",
    "AAT" = "N", "AAC" = "N", "AAA" = "K", "AAG" = "K",
    "AGT" = "S", "AGC" = "S", "AGA" = "R", "AGG" = "R",
    "GTT" = "V", "GTC" = "V", "GTA" = "V", "GTG" = "V",
    "GCT" = "A", "GCC" = "A", "GCA" = "A", "GCG" = "A",
    "GAT" = "D", "GAC" = "D", "GAA" = "E", "GAG" = "E",
    "GGT" = "G", "GGC" = "G", "GGA" = "G", "GGG" = "G"  
    )
  
  return(codon_table[[codon]])
}

# Example usage:
seq1 <- "ATTACATTTTGGAGCCAGTTTCTTGATGCCCCCTCGCTGGCTGTCCACATCCTGAGTGATCGCCTGCTGGAGTTCGTACTGGACTGGTATTTGGGCGAGCACAAGCGGGTATCGGGTCCACCGAGCCTGGAACAGCAGACGACGATTACGAAGAGCGCTGTGGAGTTGGCGCAACAGATACGGGAACGCAGGCAGCGCAGTTACGACATAGTCAAGGCCTACTGTGAGCGAATCGAGAGCGTTAATCGGGATCTCAATGCGGTGGTCGATGGCCCGTTTCCGGAGGCCTTGGACCAGGCACGCGAAATCGATCGCAAGTTGGACGAAAAGGAGTACAGCGATGAGGATCTGCGGCGCCTGCCCTTCCTAGGCGTACCCTTCTCCACCAAGGATAGCACCGCTGTGGCTGGGAAGCTACACACACTGGGCCTGCTAGCAAGGAAATCGGAACGTTCGACCACGGATGCGGAATGCGTGCGCCTGATGAAGGAAAGTGGCGCCATTATCATCGCCACCAGCAATGTGCCGGAGGTGAACAAGTGGATTGAGTCGAGGAACATGCTAATCGGCTGCACCAACAATCCCTACGATCTGCGGCGTTCCGTCGGCGGTTCTTCGGGTGGTGAAGCCGCTCTGATCGCTGCCTGCTGCACGGGTTTTGGTTTGGGAACAGACATTGGTGGTTCCATACGCATACCGGCCTTCAATTGCGGCATCTTCGGGCACAAACCCACTTCCGGTGCGGTCAATATGGCCGGATGTACGTTCAGAACGGGCAAGGAAAAGGATACGATGGTATGTGCCGGTCCAATGAGTCGTTCTGCGCGGGATCTACTGCCCATGATGCAGGTTCTGGTGGAGCCATCGTTGAAAGCCAAGCTAAAACTGGATCAGAAGGTGGATCTGAAGCGGCTGCGTTACTTCTATGTGTCCTCCAATGGCATGGCTCAGTGCAATCCTATTAACCGGGAAACCGAACGGGTTATGTACAAGATTCGCAAACATTTCGAGGCAGTGAGTGGAAAGGACGTGCGGCACGCTGACTTGCCCTATACCAAGCTGACCGGAAAGATGTGGCGCTATTGGATGACCCAGGAACCGGCCAACTTCAATCTGCTGTTGGGCAATGGCGCCGAGCTCAATCCATTTGTGGAGCTCTTCAAGAAGATACTGGGGCAAAGTGACTACAGCATGGCGGCCATTTATGGTCTGATCGACAGTGTGCTGCCCAAGGAAAAGGAGAAGCTGATGAGGGAAGCCACGGCCAAGTGCAAGAAGTCCGTCCAGGACTTGCTCGGCGACGATGGTGTGCTGTTCTTTCACAGTTCTCCCAGGACTGCACCATTCCACTACTATCCTCTGGTCAAGTTCAATGACTTCGCGTACTTCAGCCTCTTCAATGTGCTCCACTTGCCGGCCACCCAGGTGCCCATGGGTCTGGACTCCAAGGGCATGCCGCTGGGCATCCAAGTGGTGGCCAATCCGAACAACGATCGCCTTTGCCTGGCGGTTGCCGAGGAACTGGAGCGCACTTTTGGCGGCTGGGTGCCTCCATTTCCGCTGAAGAAC"
seq2 <- "ATGTCATTTTGGAGCCACGTACTTGATGCCCTGCTGGCACTGGTCCACATCCTGAGCGATCGCCTGCTGGAGTTCGTACTGGACTGGTATTTGGGCGAGCACAAGCGGGTATCGGGTCCACCGAGCCTGGAACAGCAGACGACGATTACGAAGAGCGCTGTGGAGTTGGCGCAACAGATACGGGAACGCAGGCAGCGCAGTTACGACATAGTCAAGGCCTACTGTGAGCGAATCGAGAGCGTTAATCGGGATCTCAATGCGGTGGTCGATGGCCCGTTTCCGGAGGCCTTGGACCAGGCACGCGAAATCGATCGCAAGTTGGACGAAAAGGAGTACAGCGATGAGGATCTGCGGCGCCTGCCCTTCCTAGGCGTACCCTTCTCCACCAAGGATAGCACCGCTGTGGCTGGGAAGCTACACACACTGGGCCTGCTAGCAAGGAAATCGGAACGTTCGACCACGGATGCGGAATGCGTGCGCCTGATGAAGGAAAGTGGCGCCATTATCATCGCCACCAGCAATGTGCCGGAGGTGAACAAGTGGATTGAGTCGAGGAACATGCTAATCGGCTGCACCAACAATCCCTACGATCTGCGGCGTTCCGTCGGCGGTTCTTCGGGTGGTGAAGCCGCTCTGATCGCTGCCTGCTGCACGGGTTTTGGTTTGGGAACAGACATTGGTGGTTCCATACGCATACCGGCCTTCAATTGCGGCATCTTCGGGCACAAACCCACTTCCGGTGCGGTCAATATGGCCGGATGTACGTTCAGAACGGGCAAGGAAAAGGATACGATGGTATGTGCCGGTCCAATGAGTCGTTCTGCGCGGGATCTACTGCCCATGATGCAGGTTCTGGTGGAGCCATCGTTGAAAGCCAAGCTAAAACTGGATCAGAAGGTGGATCTGAAGCGGCTGCGTTACTTCTATGTGTCCTCCAATGGCATGGCTCAGTGCAATCCTATTAACCGGGAAACCGAACGGGTTATGTACAAGATTCGCAAACATTTCGAGGCAGTGAGTGGAAAGGACGTGCGGCACGCTGACTTGCCCTATACCAAGCTGACCGGAAAGATGTGGCGCTATTGGATGACCCAGGAACCGGCCAACTTCAATCTGCTGTTGGGCAATGGCGCCGAGCTCAATCCATTTGTGGAGCTCTTCAAGAAGATACTGGGGCAAAGTGACTACAGCATGGCGGCCATTTATGGTCTGATCGACAGTGTGCTGCCCAAGGAAAAGGAGAAGCTGATGAGGGAAGCCACGGCCAAGTGCAAGAAGTCCGTCCAGGACTTGCTCGGCGACGATGGTGTGCTGTTCTTTCACAGTTCTCCCAGGACTGCACCATTCCACTACTATCCTCTGGTCAAGTTCAATGACTTCGCGTACTTCAGCCTCTTCAATGTGCTCCACTTGCCGGCCACCCAGGTGCCCATGGGTCTGGACTCCAAGGGCATGCCGCTGGGCATC-AAGTGGTGGCCAATCCGAACAACGATCGCCTTTGCCTGGCGGTTGCCGAGGAACTGGAGCGCACTTTTGGCGGCTGGGTGCCTCCATTTCCGC--------"
mutations <- calculate_mutations(seq1, seq2)
cat("Nonsynonymous mutations in seq1:", mutations$nonsynonymous1, "\n")
cat("Nonsynonymous mutations in seq2:", mutations$nonsynonymous1, "\n")
cat("Synonymous mutations:", mutations$synonymous, "\n")
cat("Nonsense mutations in seq1:", mutations$nonsense1, "\n")
cat("Nonsense mutations in seq2:", mutations$nonsense2, "\n")


df <- data.frame(
  Type = c("Nonsynonymous1", "Nonsynonymous2", "Synonymous", "Nonsense1", "Nonsense2"),
  Count = c(mutations$nonsynonymous1, mutations$nonsynonymous2, mutations$synonymous, mutations$nonsense1, mutations$nonsense2),
  Positions = I(list(mutations$nonsynonymous1_positions, mutations$nonsynonymous2_positions, mutations$synonymous_positions, mutations$nonsense1_positions, mutations$nonsense2_positions))
)

df$Positions <- ifelse(sapply(df$Positions, is.null), NA, df$Positions)

