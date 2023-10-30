





source("startup.R")


########
# get FBgns



#########



# import my duplicate genes 
dups <- read.csv("C:/Users/17735/Downloads/Dmel_Duplicate_Genes/Duplicate_Proteins.tsv", sep="")




fbgn_prots <- read.csv("./dmel-all-translation-r6.54.fasta", sep="", header = F)

# format the fasta into a dataframe 
colnames(fbgn_prots)[1] <- 'prot_seq'
fbgn_prots <- fbgn_prots %>%
  mutate(Group = cumsum(grepl("^>", prot_seq))) %>%
  group_by(Group) %>%
  mutate(prot = paste(prot_seq, collapse = "")) %>%
  ungroup() %>%
  mutate(prot_seq = prot) %>%
  select(-Group,-prot)

fbgn_prots$FBpp <- substr(fbgn_prots$prot_seq, start = 2, stop = 12)
fbgn_prots$prot_seq <- substring(fbgn_prots$prot_seq, first = 13)


# remove the unnecessary rows 
fbgn_prots <- fbgn_prots %>% mutate_all(~na_if(., "")) %>% na.omit()

# extract the FBgn id 
fbgn_prots$FBgn <- substr(fbgn_prots$V6, start = 8, stop = 18)


fbgn_prots <- fbgn_prots[,c(12,13,1,3,5:9)]


















