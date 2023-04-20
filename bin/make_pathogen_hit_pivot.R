#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# loading necessary libraries
library(tidyverse)
library(readxl)
library(openxlsx)

# bring in virus name lookup table
pathogen_lookup <- read.delim(args[1], header = F, 
                           col.names = c("refseq_acc", "common"))

# record sequencing run name from samplesheet
run_name <- str_remove(str_remove(args[2], "samplesheet_"), "_samplesheet")

# create a list of all "hits" files and read each of them into one R data frame
# as their own column
hit_files <- list.files(path = ".", pattern = "*_hits.txt")
for (i in 1:length(hit_files)){
  
  if (i == 1){
    
    pivot <- read.delim(hit_files[i], header = F, col.names = c("refseq_acc", "reads"))
    pivot <- pivot[order(pivot$refseq_acc),] ; rownames(pivot) <- NULL
    colnames(pivot)[2] <- str_remove(hit_files[i], "_hits.txt")
    
  } else {
    
    hits <- read.delim(hit_files[i], header = F, col.names = c("refseq_acc", "reads"))
    hits <- hits[order(hits$refseq_acc),] ; rownames(hits) <- NULL
    sample_id <- str_remove(hit_files[i], "_hits.txt")
    colnames(hits)[2] <- sample_id
    
    pivot <- cbind(pivot, hits[,sample_id])
    colnames(pivot)[i+1] <- sample_id
    
  }
  
}

# make refseq accessions the rownames so they don't interfere with number
# parsing downstream
refseq_acc <- pivot$refseq_acc
rownames(pivot) <- refseq_acc
if (colnames(pivot)[1] == "refseq_acc"){
  pivot <- pivot[,2:ncol(pivot)]
}

# get rid of rows that only contain zeroes
pivot <- pivot[rowSums(pivot) > 1, ]

# add column for virus common names
pivot$`Virus Common Name` <- as.character(NA)
for (i in 1:nrow(pivot)){
  
  accession <- rownames(pivot)[i]
  common_name <- pathogen_lookup[pathogen_lookup$refseq_acc==accession,"common"]
  pivot$`Virus Common Name`[i] <- common_name
  
}

# writing and exporting pivot excel files
write.xlsx(pivot, paste(run_name, "_pathogen_hits.xlsx", sep = ""),
           rowNames = TRUE, colNames = TRUE, quote = FALSE, na = "")
