#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# loading necessary libraries
library(tidyverse)

# bringing in reference library of pathogens in FASTA format
fasta <- read.delim(args[1], header = F)

# retaining only the deflines/getting rid of sequences
deflines <- which(grepl(">", fasta$V1))
fasta <- fasta[deflines,] ; rownames(fasta) <- NULL

# splitting names into two columns
lookup <- data.frame("NCBI_RefSeq_Accession" = rep(NA, length(fasta)),
                     "Pathogen_Common_Name" = rep(NA, length(fasta)))
for (i in 1:length(fasta)){
  
  defline = str_remove(fasta[i], ">")
  accession <- unlist(strsplit(defline, " "))[1]
  lookup$NCBI_RefSeq_Accession[i] <- accession
  common_name <- str_remove(str_remove(defline, accession), " ")
  lookup$Pathogen_Common_Name[i] <- common_name
  
}

# export lookup as lightweight tab-delimited file
write.table(lookup, "pathogen_name_lookup.tsv", sep = "\t", 
            col.names = F, row.names = F, quote = F)
