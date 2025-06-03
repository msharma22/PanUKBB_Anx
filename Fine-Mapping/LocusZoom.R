######Worrier/Anxious Feelings Phenotype(ID: 1980)#####
setwd("/home/maya/LocusZoom/")
library(data.table)
library(dplyr)
#read in the data file
example1980 <- read.csv("/home/maya/PANUKBB/Pan-UKBB_sumstats/hg37_varID_categorical-1980-both_sexes.tsv.gz", sep="\t")
#use dplyr to edit out unnecessary data
lz1980 <- example1980 %>%
  #select function chooses the columns in the data that you want to be included in the file
  select(chr, pos, ref, alt, beta_meta_hq, se_meta_hq, neglog10_pval_meta_hq)
#fwrite writes the data from the variable into a file
#the file needs to be tab-delimited for LocusZoom, so sep="\t"
#"exampleData" is the name of the file that will be created with this command
lz1980Filtered<-filter(lz1980, !is.na(neglog10_pval_meta_hq))
fwrite(lz1980Filtered, "1980METAData", append= FALSE, quote = "auto", sep="\t", row.names = FALSE, col.names=TRUE)
######################################Worrier/Anxious Feelings Phenotype - AFR (ID: 1980)##############################################
AFR1980 <- example1980 %>%
  select(chr, pos, ref, alt, beta_AFR, se_AFR, neglog10_pval_AFR)
AFR1980Filtered<-filter(AFR1980, !is.na(neglog10_pval_AFR))
fwrite(AFR1980Filtered, "1980AFR", append= FALSE, quote = "auto", sep="\t", row.names = FALSE, col.names=TRUE)
######################################Worrier/Anxious Feelings Phenotype - AMR (ID: 1980)##############################################
AMR1980 <- example1980 %>%
  select(chr, pos, ref, alt, beta_AMR, se_AMR, neglog10_pval_AMR)
AMR1980Filtered<-filter(AMR1980, !is.na(neglog10_pval_AMR))
fwrite(AMR1980Filtered, "1980AMR", append= FALSE, quote = "auto", sep="\t", row.names = FALSE, col.names=TRUE)
######################################Worrier/Anxious Feelings Phenotype - CSA (ID: 1980)##############################################
CSA1980 <- example1980 %>%
  select(chr, pos, ref, alt, beta_CSA, se_CSA, neglog10_pval_CSA)
CSA1980Filtered<-filter(CSA1980, !is.na(neglog10_pval_CSA))
fwrite(CSA1980Filtered, "1980CSA", append= FALSE, quote = "auto", sep="\t", row.names = FALSE, col.names=TRUE)
######################################Worrier/Anxious Feelings Phenotype - EAS (ID: 1980)##############################################
EAS1980 <- example1980 %>%
  select(chr, pos, ref, alt, beta_EAS, se_EAS, neglog10_pval_EAS)
EAS1980Filtered<-filter(EAS1980, !is.na(neglog10_pval_EAS))
fwrite(EAS1980Filtered, "1980EAS", append= FALSE, quote = "auto", sep="\t", row.names = FALSE, col.names=TRUE)
######################################Worrier/Anxious Feelings Phenotype - MID (ID: 1980)##############################################
MID1980 <- example1980 %>%
  select(chr, pos, ref, alt, beta_MID, se_MID, neglog10_pval_MID)
MID1980Filtered<-filter(MID1980, !is.na(neglog10_pval_MID))
fwrite(MID1980Filtered, "1980MID", append= FALSE, quote = "auto", sep="\t", row.names = FALSE, col.names=TRUE)
