library(data.table)
library(dplyr)
library(coloc)
library(hash)
###cis eQTL file: '/home/grace/coloc/eQTL_data/Hypothalamus/Brain_Spinal_cord_cervical_c-1.allpairs.txt.gz'
### Corrected cis eQTL file: '/home/grace/coloc/eQTL_data/Spinal_cord_cervical/Brain_SCC_eQTL.csv'

"%&%" = function(a,b) paste(a,b,sep="")

######### Convert GWAS file to correct format
#reformatting the GWAS data
data_GWAS <- fread("/home/maya/PANUKBB/Pan-UKBB_sumstats/hg37_varID_categorical-1980-both_sexes.tsv.gz") #GWAS file, build 38

data_GWAS <- data_GWAS %>%
  mutate(CHR = substring(chromosome,4), #get chromosome number
         BP = sapply(strsplit(panel_variant_id,"_"), `[`, 2)) %>% #get position 
  filter(nchar(effect_allele) == 1, #filter out indels
         nchar(non_effect_allele) == 1) %>% #filter out indels
  rename(FREQ = frequency, #rename to correct column names
         BETA = effect_size,
         SE = standard_error,
         P = pvalue)

data_GWAS <- data_GWAS %>%
  mutate(ID =(CHR %&% "_" %&% BP), 
         SNP_hg38 = BP) %>% #makes ID column
  rename(A1 = effect_allele, #rename to correct column names
         A2 = non_effect_allele)

GWASpath <- "/home/maya/COLOC/gwasFILES/META_ANX_coloc"
fwrite(data_GWAS, GWASpath, append=FALSE, quote ='auto', sep=',') #write new GWAS file
