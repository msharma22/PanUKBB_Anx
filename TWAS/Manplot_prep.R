### Merge all 13 tissue files and write df to a csv 
library(dplyr)
library(readr)
AFR_total_genes.csv <- list.files(path="/home/maya/prediXcan/EditedprediXcanOUTPUT/ANX/AFR/ANX_AFR", full.names = TRUE) %>% 
  lapply(read_csv) %>% 
  bind_rows 
### Format p-value column in ascending order 
AFR_total_genes <-read.csv("~/AoU/sprediXcan/prediXcanOUTPUT/EDITEDoutputs/AFR/AFR_total_genes.csv")
### Compute how many genes meet p-value threshold 
filteredAFR_total_genes = filter(AFR_total_genes, AFR_total_genes$pvalue <= 8.288177e-08)
counttotalAFR <- sum(AFR_total_genes$pvalue <= 8.288177e-08, na.rm=TRUE)
