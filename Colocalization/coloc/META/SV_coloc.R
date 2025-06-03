###1. Retrieve tissue specific cis-QTL data from GTEx Portal
### 2. Convert eQTL file to correct format - I used converted eQTLs files from grace's directory 
### 3. Convert GWAS file to correct format
library(data.table) 
library(dplyr)
library(coloc)
library(hash)

"%&%" = function(a,b) paste(a,b,sep="")

#reformatting the GWAS data
data_GWAS <- fread("/home/maya/prediXcan/harmonizedGWAS/Anx/hg38_PUKBB_META_Anx.tsv.gz") #GWAS file, build 38

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

GWASpath <- "/home/maya/COLOC/gwasFILES/META_Anx_coloc"
fwrite(data_GWAS, GWASpath, append=FALSE, quote ='auto', sep=',') #write new GWAS file
### 4. Create eQTL files for each significant gene. I pulled these files from grace's directory 

################################################## CADM2 - Hypothalamus ###########################################################################
### 5. Retrieve summary stat data from gene and phenotype of interest
library(data.table)
library(dplyr)
library(coloc)
library(hash)

"%&%" = function(a,b) paste(a,b,sep="")

GWASdir = "/home/maya/COLOC/gwasFILES/META_Anx_coloc"
eQTLdir = "/home/maya/COLOC/eQTL_data/Hypothalamus/META_Anx_colocs/CADM2_data.csv"

eQTLs = fread(eQTLdir) #read in significant gene eQTL
gwas = fread(GWASdir)

#get intersecting SNPs
joined = inner_join(eQTLs,gwas, by=c("ID"="ID"))

joined <- arrange(joined, P.y, by_group=FALSE)


#check for complement bases
#build hash table (like a Python dictionary)

bases = hash()
bases[["A"]] <- "T"
bases[["C"]] <- "G"
bases[["G"]] <- "C"
bases[["T"]] <- "A"


#pull SNPs that  match (assumes C|G and A|T SNPs aren't flipped)
match = joined[(joined$A1.x==joined$A1.y & joined$REF == joined$A2),]

#remove ambiguous strand SNPs (A/T or C/G)
a = joined[!(joined$REF=="A" & joined$ALT=="T") & !(joined$REF=="T" & joined$ALT=="A") & !(joined$REF=="C" & joined$ALT=="G") &
             !(joined$REF=="G" & joined$ALT=="C") ]

#of non-ambiguous, pull SNPs that are flipped (or the complement bases match or are flipped)

compmatch = a[(a$A1.x==values(bases,keys=a$A1.y) & a$REF == values(bases,keys=a$A2)),]

flipped = a[(a$A1.x==a$A2 & a$REF == a$A1.y),]
compflipped = a[(a$A1.x==values(bases,keys=a$A2) & a$REF == values(bases,keys=a$A1.y)),]

#if flipped, change sign of cancer beta, check to see if any SNPs are in flipped df's with if stmt
if(dim(flipped)[1] > 0){
  flipped = mutate(flipped,BETA.y = -1*BETA.y)
}
if(dim(compflipped)[1] > 0){
  compflipped = mutate(compflipped,BETA.y = -1*BETA.y)
}

#bind all and sort by position
matchsnps = rbind(match, compmatch, flipped, compflipped) %>% arrange(POS)

#matchsnps <- matchsnps %>%
#na.omit()

matchsnps <- matchsnps %>%
  distinct(BP, .keep_all = TRUE) %>%
  select(-n_cases)
matchsnps <- matchsnps %>%
  na.omit()
### 6. Format each set of summary stats for coloc
gwascoloc = list("beta" = matchsnps$BETA.y, "varbeta" = (matchsnps$SE.y)^2, "snp" = matchsnps$ID, "position" = matchsnps$BP,
                 "type" = "cc")
eqtlcoloc = list("beta" = matchsnps$BETA.x, "varbeta" = (matchsnps$SE.x)^2, "snp" = matchsnps$ID, "position" = matchsnps$POS,
                 "type" = "quant", "N" = matchsnps$OBS_CT[1], "MAF" = matchsnps$A1_FREQ, "sdY"=1 )


plot_dataset(gwascoloc)
plot_dataset(eqtlcoloc)

###7. Run coloc assuming a single causal variant
my.res = coloc.abf(dataset1=gwascoloc, dataset2=eqtlcoloc)
my.res
sensitivity(my.res,"H4 > 0.9")

################################################## CADM2 - Spinal Cord Cervical ###########################################################################
### 5. Retrieve summary stat data from gene and phenotype of interest
library(data.table)
library(dplyr)
library(coloc)
library(hash)

"%&%" = function(a,b) paste(a,b,sep="")

GWASdir = "/home/maya/COLOC/gwasFILES/META_Anx_coloc"
eQTLdir = "/home/maya/COLOC/eQTL_data/Spinal_cord_cervical/META_Anx_colocs/CADM2_data.csv"

eQTLs = fread(eQTLdir) #read in significant gene eQTL
gwas = fread(GWASdir)

#get intersecting SNPs
joined = inner_join(eQTLs,gwas, by=c("ID"="ID"))

joined <- arrange(joined, P.y, by_group=FALSE)


#check for complement bases
#build hash table (like a Python dictionary)

bases = hash()
bases[["A"]] <- "T"
bases[["C"]] <- "G"
bases[["G"]] <- "C"
bases[["T"]] <- "A"


#pull SNPs that  match (assumes C|G and A|T SNPs aren't flipped)
match = joined[(joined$A1.x==joined$A1.y & joined$REF == joined$A2),]

#remove ambiguous strand SNPs (A/T or C/G)
a = joined[!(joined$REF=="A" & joined$ALT=="T") & !(joined$REF=="T" & joined$ALT=="A") & !(joined$REF=="C" & joined$ALT=="G") &
             !(joined$REF=="G" & joined$ALT=="C") ]

#of non-ambiguous, pull SNPs that are flipped (or the complement bases match or are flipped)

compmatch = a[(a$A1.x==values(bases,keys=a$A1.y) & a$REF == values(bases,keys=a$A2)),]

flipped = a[(a$A1.x==a$A2 & a$REF == a$A1.y),]
compflipped = a[(a$A1.x==values(bases,keys=a$A2) & a$REF == values(bases,keys=a$A1.y)),]

#if flipped, change sign of cancer beta, check to see if any SNPs are in flipped df's with if stmt
if(dim(flipped)[1] > 0){
  flipped = mutate(flipped,BETA.y = -1*BETA.y)
}
if(dim(compflipped)[1] > 0){
  compflipped = mutate(compflipped,BETA.y = -1*BETA.y)
}

#bind all and sort by position
matchsnps = rbind(match, compmatch, flipped, compflipped) %>% arrange(POS)

#matchsnps <- matchsnps %>%
#na.omit()

matchsnps <- matchsnps %>%
  distinct(BP, .keep_all = TRUE) %>%
  select(-n_cases)
matchsnps <- matchsnps %>%
  na.omit()
### 6. Format each set of summary stats for coloc
gwascoloc = list("beta" = matchsnps$BETA.y, "varbeta" = (matchsnps$SE.y)^2, "snp" = matchsnps$ID, "position" = matchsnps$BP,
                 "type" = "cc")
eqtlcoloc = list("beta" = matchsnps$BETA.x, "varbeta" = (matchsnps$SE.x)^2, "snp" = matchsnps$ID, "position" = matchsnps$POS,
                 "type" = "quant", "N" = matchsnps$OBS_CT[1], "MAF" = matchsnps$A1_FREQ, "sdY"=1 )


plot_dataset(gwascoloc)
plot_dataset(eqtlcoloc)

###7. Run coloc assuming a single causal variant
my.res = coloc.abf(dataset1=gwascoloc, dataset2=eqtlcoloc)
my.res
sensitivity(my.res,"H4 > 0.9")

################################################## AF131215.8 - Hippocampus ###########################################################################
### 5. Retrieve summary stat data from gene and phenotype of interest
library(data.table)
library(dplyr)
library(coloc)
library(hash)

"%&%" = function(a,b) paste(a,b,sep="")

GWASdir = "/home/maya/COLOC/gwasFILES/META_Anx_coloc"
eQTLdir = "/home/maya/COLOC/eQTL_data/Hippocampus/META_Anx_colocs/AF131215.8_data.csv"

eQTLs = fread(eQTLdir) #read in significant gene eQTL
gwas = fread(GWASdir)

#get intersecting SNPs
joined = inner_join(eQTLs,gwas, by=c("ID"="ID"))

joined <- arrange(joined, P.y, by_group=FALSE)


#check for complement bases
#build hash table (like a Python dictionary)

bases = hash()
bases[["A"]] <- "T"
bases[["C"]] <- "G"
bases[["G"]] <- "C"
bases[["T"]] <- "A"


#pull SNPs that  match (assumes C|G and A|T SNPs aren't flipped)
match = joined[(joined$A1.x==joined$A1.y & joined$REF == joined$A2),]

#remove ambiguous strand SNPs (A/T or C/G)
a = joined[!(joined$REF=="A" & joined$ALT=="T") & !(joined$REF=="T" & joined$ALT=="A") & !(joined$REF=="C" & joined$ALT=="G") &
             !(joined$REF=="G" & joined$ALT=="C") ]

#of non-ambiguous, pull SNPs that are flipped (or the complement bases match or are flipped)

compmatch = a[(a$A1.x==values(bases,keys=a$A1.y) & a$REF == values(bases,keys=a$A2)),]

flipped = a[(a$A1.x==a$A2 & a$REF == a$A1.y),]
compflipped = a[(a$A1.x==values(bases,keys=a$A2) & a$REF == values(bases,keys=a$A1.y)),]

#if flipped, change sign of cancer beta, check to see if any SNPs are in flipped df's with if stmt
if(dim(flipped)[1] > 0){
  flipped = mutate(flipped,BETA.y = -1*BETA.y)
}
if(dim(compflipped)[1] > 0){
  compflipped = mutate(compflipped,BETA.y = -1*BETA.y)
}

#bind all and sort by position
matchsnps = rbind(match, compmatch, flipped, compflipped) %>% arrange(POS)

#matchsnps <- matchsnps %>%
#na.omit()

matchsnps <- matchsnps %>%
  distinct(BP, .keep_all = TRUE) %>%
  select(-n_cases)
matchsnps <- matchsnps %>%
  na.omit()

### 6. Format each set of summary stats for coloc
gwascoloc = list("beta" = matchsnps$BETA.y, "varbeta" = (matchsnps$SE.y)^2, "snp" = matchsnps$ID, "position" = matchsnps$BP,
                 "type" = "cc")
eqtlcoloc = list("beta" = matchsnps$BETA.x, "varbeta" = (matchsnps$SE.x)^2, "snp" = matchsnps$ID, "position" = matchsnps$POS,
                 "type" = "quant", "N" = matchsnps$OBS_CT[1], "MAF" = matchsnps$A1_FREQ, "sdY"=1 )


plot_dataset(gwascoloc)
plot_dataset(eqtlcoloc)

###7. Run coloc assuming a single causal variant
my.res = coloc.abf(dataset1=gwascoloc, dataset2=eqtlcoloc)
my.res
sensitivity(my.res,"H4 > 0.9")

################################################## AF131215.8 - Caudate Basal Ganglia ###########################################################################
### 5. Retrieve summary stat data from gene and phenotype of interest
library(data.table)
library(dplyr)
library(coloc)
library(hash)

"%&%" = function(a,b) paste(a,b,sep="")

GWASdir = "/home/maya/COLOC/gwasFILES/META_Anx_coloc"
eQTLdir = "/home/maya/COLOC/eQTL_data/Caudate_basal_ganglia/META_Anx_colocs/AF131215.8_data.csv"

eQTLs = fread(eQTLdir) #read in significant gene eQTL
gwas = fread(GWASdir)

#get intersecting SNPs
joined = inner_join(eQTLs,gwas, by=c("ID"="ID"))

joined <- arrange(joined, P.y, by_group=FALSE)


#check for complement bases
#build hash table (like a Python dictionary)

bases = hash()
bases[["A"]] <- "T"
bases[["C"]] <- "G"
bases[["G"]] <- "C"
bases[["T"]] <- "A"


#pull SNPs that  match (assumes C|G and A|T SNPs aren't flipped)
match = joined[(joined$A1.x==joined$A1.y & joined$REF == joined$A2),]

#remove ambiguous strand SNPs (A/T or C/G)
a = joined[!(joined$REF=="A" & joined$ALT=="T") & !(joined$REF=="T" & joined$ALT=="A") & !(joined$REF=="C" & joined$ALT=="G") &
             !(joined$REF=="G" & joined$ALT=="C") ]

#of non-ambiguous, pull SNPs that are flipped (or the complement bases match or are flipped)

compmatch = a[(a$A1.x==values(bases,keys=a$A1.y) & a$REF == values(bases,keys=a$A2)),]

flipped = a[(a$A1.x==a$A2 & a$REF == a$A1.y),]
compflipped = a[(a$A1.x==values(bases,keys=a$A2) & a$REF == values(bases,keys=a$A1.y)),]

#if flipped, change sign of cancer beta, check to see if any SNPs are in flipped df's with if stmt
if(dim(flipped)[1] > 0){
  flipped = mutate(flipped,BETA.y = -1*BETA.y)
}
if(dim(compflipped)[1] > 0){
  compflipped = mutate(compflipped,BETA.y = -1*BETA.y)
}

#bind all and sort by position
matchsnps = rbind(match, compmatch, flipped, compflipped) %>% arrange(POS)

#matchsnps <- matchsnps %>%
#na.omit()

matchsnps <- matchsnps %>%
  distinct(BP, .keep_all = TRUE) %>%
  select(-n_cases)
matchsnps <- matchsnps %>%
  na.omit()

### 6. Format each set of summary stats for coloc
gwascoloc = list("beta" = matchsnps$BETA.y, "varbeta" = (matchsnps$SE.y)^2, "snp" = matchsnps$ID, "position" = matchsnps$BP,
                 "type" = "cc")
eqtlcoloc = list("beta" = matchsnps$BETA.x, "varbeta" = (matchsnps$SE.x)^2, "snp" = matchsnps$ID, "position" = matchsnps$POS,
                 "type" = "quant", "N" = matchsnps$OBS_CT[1], "MAF" = matchsnps$A1_FREQ, "sdY"=1 )


plot_dataset(gwascoloc)
plot_dataset(eqtlcoloc)

###7. Run coloc assuming a single causal variant
my.res = coloc.abf(dataset1=gwascoloc, dataset2=eqtlcoloc)
my.res
sensitivity(my.res,"H4 > 0.9")
################################################## AF131215.8 - Hypothalamus ###########################################################################
### 5. Retrieve summary stat data from gene and phenotype of interest
library(data.table)
library(dplyr)
library(coloc)
library(hash)

"%&%" = function(a,b) paste(a,b,sep="")

GWASdir = "/home/maya/COLOC/gwasFILES/META_Anx_coloc"
eQTLdir = "/home/maya/COLOC/eQTL_data/Hypothalamus/META_Anx_colocs/AF131215.8_data.csv"

eQTLs = fread(eQTLdir) #read in significant gene eQTL
gwas = fread(GWASdir)

#get intersecting SNPs
joined = inner_join(eQTLs,gwas, by=c("ID"="ID"))

joined <- arrange(joined, P.y, by_group=FALSE)


#check for complement bases
#build hash table (like a Python dictionary)

bases = hash()
bases[["A"]] <- "T"
bases[["C"]] <- "G"
bases[["G"]] <- "C"
bases[["T"]] <- "A"


#pull SNPs that  match (assumes C|G and A|T SNPs aren't flipped)
match = joined[(joined$A1.x==joined$A1.y & joined$REF == joined$A2),]

#remove ambiguous strand SNPs (A/T or C/G)
a = joined[!(joined$REF=="A" & joined$ALT=="T") & !(joined$REF=="T" & joined$ALT=="A") & !(joined$REF=="C" & joined$ALT=="G") &
             !(joined$REF=="G" & joined$ALT=="C") ]

#of non-ambiguous, pull SNPs that are flipped (or the complement bases match or are flipped)

compmatch = a[(a$A1.x==values(bases,keys=a$A1.y) & a$REF == values(bases,keys=a$A2)),]

flipped = a[(a$A1.x==a$A2 & a$REF == a$A1.y),]
compflipped = a[(a$A1.x==values(bases,keys=a$A2) & a$REF == values(bases,keys=a$A1.y)),]

#if flipped, change sign of cancer beta, check to see if any SNPs are in flipped df's with if stmt
if(dim(flipped)[1] > 0){
  flipped = mutate(flipped,BETA.y = -1*BETA.y)
}
if(dim(compflipped)[1] > 0){
  compflipped = mutate(compflipped,BETA.y = -1*BETA.y)
}

#bind all and sort by position
matchsnps = rbind(match, compmatch, flipped, compflipped) %>% arrange(POS)

#matchsnps <- matchsnps %>%
#na.omit()

matchsnps <- matchsnps %>%
  distinct(BP, .keep_all = TRUE) %>%
  select(-n_cases)
matchsnps <- matchsnps %>%
  na.omit()

### 6. Format each set of summary stats for coloc
gwascoloc = list("beta" = matchsnps$BETA.y, "varbeta" = (matchsnps$SE.y)^2, "snp" = matchsnps$ID, "position" = matchsnps$BP,
                 "type" = "cc")
eqtlcoloc = list("beta" = matchsnps$BETA.x, "varbeta" = (matchsnps$SE.x)^2, "snp" = matchsnps$ID, "position" = matchsnps$POS,
                 "type" = "quant", "N" = matchsnps$OBS_CT[1], "MAF" = matchsnps$A1_FREQ, "sdY"=1 )


plot_dataset(gwascoloc)
plot_dataset(eqtlcoloc)

###7. Run coloc assuming a single causal variant
my.res = coloc.abf(dataset1=gwascoloc, dataset2=eqtlcoloc)
my.res
sensitivity(my.res,"H4 > 0.9")

################################################## AF131215.8 - Nucleus Accumbens Basal Ganglia ###########################################################################
### 5. Retrieve summary stat data from gene and phenotype of interest

GWASdir = "/home/maya/COLOC/gwasFILES/META_Anx_coloc"
eQTLdir = "/home/maya/COLOC/eQTL_data/Nucleus_accumbens_basal_ganglia/META_Anx_colocs/AF131215.8_data.csv"

eQTLs = fread(eQTLdir) #read in significant gene eQTL
gwas = fread(GWASdir)

#get intersecting SNPs
joined = inner_join(eQTLs,gwas, by=c("ID"="ID"))

joined <- arrange(joined, P.y, by_group=FALSE)


#check for complement bases
#build hash table (like a Python dictionary)

bases = hash()
bases[["A"]] <- "T"
bases[["C"]] <- "G"
bases[["G"]] <- "C"
bases[["T"]] <- "A"


#pull SNPs that  match (assumes C|G and A|T SNPs aren't flipped)
match = joined[(joined$A1.x==joined$A1.y & joined$REF == joined$A2),]

#remove ambiguous strand SNPs (A/T or C/G)
a = joined[!(joined$REF=="A" & joined$ALT=="T") & !(joined$REF=="T" & joined$ALT=="A") & !(joined$REF=="C" & joined$ALT=="G") &
             !(joined$REF=="G" & joined$ALT=="C") ]

#of non-ambiguous, pull SNPs that are flipped (or the complement bases match or are flipped)

compmatch = a[(a$A1.x==values(bases,keys=a$A1.y) & a$REF == values(bases,keys=a$A2)),]

flipped = a[(a$A1.x==a$A2 & a$REF == a$A1.y),]
compflipped = a[(a$A1.x==values(bases,keys=a$A2) & a$REF == values(bases,keys=a$A1.y)),]

#if flipped, change sign of cancer beta, check to see if any SNPs are in flipped df's with if stmt
if(dim(flipped)[1] > 0){
  flipped = mutate(flipped,BETA.y = -1*BETA.y)
}
if(dim(compflipped)[1] > 0){
  compflipped = mutate(compflipped,BETA.y = -1*BETA.y)
}

#bind all and sort by position
matchsnps = rbind(match, compmatch, flipped, compflipped) %>% arrange(POS)

#matchsnps <- matchsnps %>%
#na.omit()

matchsnps <- matchsnps %>%
  distinct(BP, .keep_all = TRUE) %>%
  select(-n_cases)
matchsnps <- matchsnps %>%
  na.omit()

### 6. Format each set of summary stats for coloc
gwascoloc = list("beta" = matchsnps$BETA.y, "varbeta" = (matchsnps$SE.y)^2, "snp" = matchsnps$ID, "position" = matchsnps$BP,
                 "type" = "cc")
eqtlcoloc = list("beta" = matchsnps$BETA.x, "varbeta" = (matchsnps$SE.x)^2, "snp" = matchsnps$ID, "position" = matchsnps$POS,
                 "type" = "quant", "N" = matchsnps$OBS_CT[1], "MAF" = matchsnps$A1_FREQ, "sdY"=1 )


plot_dataset(gwascoloc)
plot_dataset(eqtlcoloc)

###7. Run coloc assuming a single causal variant
my.res = coloc.abf(dataset1=gwascoloc, dataset2=eqtlcoloc)
my.res
sensitivity(my.res,"H4 > 0.9")
################################################## AF131215.8 - Anterior Cingulate Cortex ###########################################################################
### 5. Retrieve summary stat data from gene and phenotype of interest

GWASdir = "/home/maya/COLOC/gwasFILES/META_Anx_coloc"
eQTLdir = "/home/maya/COLOC/eQTL_data/Anterior_cingulate_cortex/META_Anx_colocs/AF131215.8_data.csv"

eQTLs = fread(eQTLdir) #read in significant gene eQTL
gwas = fread(GWASdir)

#get intersecting SNPs
joined = inner_join(eQTLs,gwas, by=c("ID"="ID"))

joined <- arrange(joined, P.y, by_group=FALSE)


#check for complement bases
#build hash table (like a Python dictionary)

bases = hash()
bases[["A"]] <- "T"
bases[["C"]] <- "G"
bases[["G"]] <- "C"
bases[["T"]] <- "A"


#pull SNPs that  match (assumes C|G and A|T SNPs aren't flipped)
match = joined[(joined$A1.x==joined$A1.y & joined$REF == joined$A2),]

#remove ambiguous strand SNPs (A/T or C/G)
a = joined[!(joined$REF=="A" & joined$ALT=="T") & !(joined$REF=="T" & joined$ALT=="A") & !(joined$REF=="C" & joined$ALT=="G") &
             !(joined$REF=="G" & joined$ALT=="C") ]

#of non-ambiguous, pull SNPs that are flipped (or the complement bases match or are flipped)

compmatch = a[(a$A1.x==values(bases,keys=a$A1.y) & a$REF == values(bases,keys=a$A2)),]

flipped = a[(a$A1.x==a$A2 & a$REF == a$A1.y),]
compflipped = a[(a$A1.x==values(bases,keys=a$A2) & a$REF == values(bases,keys=a$A1.y)),]

#if flipped, change sign of cancer beta, check to see if any SNPs are in flipped df's with if stmt
if(dim(flipped)[1] > 0){
  flipped = mutate(flipped,BETA.y = -1*BETA.y)
}
if(dim(compflipped)[1] > 0){
  compflipped = mutate(compflipped,BETA.y = -1*BETA.y)
}

#bind all and sort by position
matchsnps = rbind(match, compmatch, flipped, compflipped) %>% arrange(POS)

#matchsnps <- matchsnps %>%
#na.omit()

matchsnps <- matchsnps %>%
  distinct(BP, .keep_all = TRUE) %>%
  select(-n_cases)
matchsnps <- matchsnps %>%
  na.omit()

### 6. Format each set of summary stats for coloc
gwascoloc = list("beta" = matchsnps$BETA.y, "varbeta" = (matchsnps$SE.y)^2, "snp" = matchsnps$ID, "position" = matchsnps$BP,
                 "type" = "cc")
eqtlcoloc = list("beta" = matchsnps$BETA.x, "varbeta" = (matchsnps$SE.x)^2, "snp" = matchsnps$ID, "position" = matchsnps$POS,
                 "type" = "quant", "N" = matchsnps$OBS_CT[1], "MAF" = matchsnps$A1_FREQ, "sdY"=1 )


plot_dataset(gwascoloc)
plot_dataset(eqtlcoloc)

###7. Run coloc assuming a single causal variant
my.res = coloc.abf(dataset1=gwascoloc, dataset2=eqtlcoloc)
my.res
sensitivity(my.res,"H4 > 0.9")
################################################## AF131215.8 - Putamen Basal Ganglia ###########################################################################
### 5. Retrieve summary stat data from gene and phenotype of interest

GWASdir = "/home/maya/COLOC/gwasFILES/META_Anx_coloc"
eQTLdir = "/home/maya/COLOC/eQTL_data/Putamen_basal_ganglia/META_Anx_colocs/AF131215.8_data.csv"

eQTLs = fread(eQTLdir) #read in significant gene eQTL
gwas = fread(GWASdir)

#get intersecting SNPs
joined = inner_join(eQTLs,gwas, by=c("ID"="ID"))

joined <- arrange(joined, P.y, by_group=FALSE)


#check for complement bases
#build hash table (like a Python dictionary)

bases = hash()
bases[["A"]] <- "T"
bases[["C"]] <- "G"
bases[["G"]] <- "C"
bases[["T"]] <- "A"


#pull SNPs that  match (assumes C|G and A|T SNPs aren't flipped)
match = joined[(joined$A1.x==joined$A1.y & joined$REF == joined$A2),]

#remove ambiguous strand SNPs (A/T or C/G)
a = joined[!(joined$REF=="A" & joined$ALT=="T") & !(joined$REF=="T" & joined$ALT=="A") & !(joined$REF=="C" & joined$ALT=="G") &
             !(joined$REF=="G" & joined$ALT=="C") ]

#of non-ambiguous, pull SNPs that are flipped (or the complement bases match or are flipped)

compmatch = a[(a$A1.x==values(bases,keys=a$A1.y) & a$REF == values(bases,keys=a$A2)),]

flipped = a[(a$A1.x==a$A2 & a$REF == a$A1.y),]
compflipped = a[(a$A1.x==values(bases,keys=a$A2) & a$REF == values(bases,keys=a$A1.y)),]

#if flipped, change sign of cancer beta, check to see if any SNPs are in flipped df's with if stmt
if(dim(flipped)[1] > 0){
  flipped = mutate(flipped,BETA.y = -1*BETA.y)
}
if(dim(compflipped)[1] > 0){
  compflipped = mutate(compflipped,BETA.y = -1*BETA.y)
}

#bind all and sort by position
matchsnps = rbind(match, compmatch, flipped, compflipped) %>% arrange(POS)

#matchsnps <- matchsnps %>%
#na.omit()

matchsnps <- matchsnps %>%
  distinct(BP, .keep_all = TRUE) %>%
  select(-n_cases)
matchsnps <- matchsnps %>%
  na.omit()

### 6. Format each set of summary stats for coloc
gwascoloc = list("beta" = matchsnps$BETA.y, "varbeta" = (matchsnps$SE.y)^2, "snp" = matchsnps$ID, "position" = matchsnps$BP,
                 "type" = "cc")
eqtlcoloc = list("beta" = matchsnps$BETA.x, "varbeta" = (matchsnps$SE.x)^2, "snp" = matchsnps$ID, "position" = matchsnps$POS,
                 "type" = "quant", "N" = matchsnps$OBS_CT[1], "MAF" = matchsnps$A1_FREQ, "sdY"=1 )


plot_dataset(gwascoloc)
plot_dataset(eqtlcoloc)

###7. Run coloc assuming a single causal variant
my.res = coloc.abf(dataset1=gwascoloc, dataset2=eqtlcoloc)
my.res
sensitivity(my.res,"H4 > 0.9")
################################################## PRAG1 - Cortex ###########################################################################
### 5. Retrieve summary stat data from gene and phenotype of interest

GWASdir = "/home/maya/COLOC/gwasFILES/META_Anx_coloc"
eQTLdir = "/home/maya/COLOC/eQTL_data/Cortex/META_Anx_colocs/PRAG1_data.csv"

eQTLs = fread(eQTLdir) #read in significant gene eQTL
gwas = fread(GWASdir)

#get intersecting SNPs
joined = inner_join(eQTLs,gwas, by=c("ID"="ID"))

joined <- arrange(joined, P.y, by_group=FALSE)


#check for complement bases
#build hash table (like a Python dictionary)

bases = hash()
bases[["A"]] <- "T"
bases[["C"]] <- "G"
bases[["G"]] <- "C"
bases[["T"]] <- "A"


#pull SNPs that  match (assumes C|G and A|T SNPs aren't flipped)
match = joined[(joined$A1.x==joined$A1.y & joined$REF == joined$A2),]

#remove ambiguous strand SNPs (A/T or C/G)
a = joined[!(joined$REF=="A" & joined$ALT=="T") & !(joined$REF=="T" & joined$ALT=="A") & !(joined$REF=="C" & joined$ALT=="G") &
             !(joined$REF=="G" & joined$ALT=="C") ]

#of non-ambiguous, pull SNPs that are flipped (or the complement bases match or are flipped)

compmatch = a[(a$A1.x==values(bases,keys=a$A1.y) & a$REF == values(bases,keys=a$A2)),]

flipped = a[(a$A1.x==a$A2 & a$REF == a$A1.y),]
compflipped = a[(a$A1.x==values(bases,keys=a$A2) & a$REF == values(bases,keys=a$A1.y)),]

#if flipped, change sign of cancer beta, check to see if any SNPs are in flipped df's with if stmt
if(dim(flipped)[1] > 0){
  flipped = mutate(flipped,BETA.y = -1*BETA.y)
}
if(dim(compflipped)[1] > 0){
  compflipped = mutate(compflipped,BETA.y = -1*BETA.y)
}

#bind all and sort by position
matchsnps = rbind(match, compmatch, flipped, compflipped) %>% arrange(POS)

#matchsnps <- matchsnps %>%
#na.omit()

matchsnps <- matchsnps %>%
  distinct(BP, .keep_all = TRUE) %>%
  select(-n_cases)
matchsnps <- matchsnps %>%
  na.omit()

### 6. Format each set of summary stats for coloc
gwascoloc = list("beta" = matchsnps$BETA.y, "varbeta" = (matchsnps$SE.y)^2, "snp" = matchsnps$ID, "position" = matchsnps$BP,
                 "type" = "cc")
eqtlcoloc = list("beta" = matchsnps$BETA.x, "varbeta" = (matchsnps$SE.x)^2, "snp" = matchsnps$ID, "position" = matchsnps$POS,
                 "type" = "quant", "N" = matchsnps$OBS_CT[1], "MAF" = matchsnps$A1_FREQ, "sdY"=1 )


plot_dataset(gwascoloc)
plot_dataset(eqtlcoloc)

###7. Run coloc assuming a single causal variant
my.res = coloc.abf(dataset1=gwascoloc, dataset2=eqtlcoloc)
my.res
sensitivity(my.res,"H4 > 0.9")
################################################## PRAG1 - Substantia Nigra ###########################################################################
### 5. Retrieve summary stat data from gene and phenotype of interest

GWASdir = "/home/maya/COLOC/gwasFILES/META_Anx_coloc"
eQTLdir = "/home/maya/COLOC/eQTL_data/Substantia_nigra/META_Anx_colocs/PRAG1_data.csv"

eQTLs = fread(eQTLdir) #read in significant gene eQTL
gwas = fread(GWASdir)

#get intersecting SNPs
joined = inner_join(eQTLs,gwas, by=c("ID"="ID"))

joined <- arrange(joined, P.y, by_group=FALSE)


#check for complement bases
#build hash table (like a Python dictionary)

bases = hash()
bases[["A"]] <- "T"
bases[["C"]] <- "G"
bases[["G"]] <- "C"
bases[["T"]] <- "A"


#pull SNPs that  match (assumes C|G and A|T SNPs aren't flipped)
match = joined[(joined$A1.x==joined$A1.y & joined$REF == joined$A2),]

#remove ambiguous strand SNPs (A/T or C/G)
a = joined[!(joined$REF=="A" & joined$ALT=="T") & !(joined$REF=="T" & joined$ALT=="A") & !(joined$REF=="C" & joined$ALT=="G") &
             !(joined$REF=="G" & joined$ALT=="C") ]

#of non-ambiguous, pull SNPs that are flipped (or the complement bases match or are flipped)

compmatch = a[(a$A1.x==values(bases,keys=a$A1.y) & a$REF == values(bases,keys=a$A2)),]

flipped = a[(a$A1.x==a$A2 & a$REF == a$A1.y),]
compflipped = a[(a$A1.x==values(bases,keys=a$A2) & a$REF == values(bases,keys=a$A1.y)),]

#if flipped, change sign of cancer beta, check to see if any SNPs are in flipped df's with if stmt
if(dim(flipped)[1] > 0){
  flipped = mutate(flipped,BETA.y = -1*BETA.y)
}
if(dim(compflipped)[1] > 0){
  compflipped = mutate(compflipped,BETA.y = -1*BETA.y)
}

#bind all and sort by position
matchsnps = rbind(match, compmatch, flipped, compflipped) %>% arrange(POS)

#matchsnps <- matchsnps %>%
#na.omit()

matchsnps <- matchsnps %>%
  distinct(BP, .keep_all = TRUE) %>%
  select(-n_cases)
matchsnps <- matchsnps %>%
  na.omit()

### 6. Format each set of summary stats for coloc
gwascoloc = list("beta" = matchsnps$BETA.y, "varbeta" = (matchsnps$SE.y)^2, "snp" = matchsnps$ID, "position" = matchsnps$BP,
                 "type" = "cc")
eqtlcoloc = list("beta" = matchsnps$BETA.x, "varbeta" = (matchsnps$SE.x)^2, "snp" = matchsnps$ID, "position" = matchsnps$POS,
                 "type" = "quant", "N" = matchsnps$OBS_CT[1], "MAF" = matchsnps$A1_FREQ, "sdY"=1 )


plot_dataset(gwascoloc)
plot_dataset(eqtlcoloc)

###7. Run coloc assuming a single causal variant
my.res = coloc.abf(dataset1=gwascoloc, dataset2=eqtlcoloc)
my.res
sensitivity(my.res,"H4 > 0.9")
################################################## FAM167A - Hypothalamus ###########################################################################
### 5. Retrieve summary stat data from gene and phenotype of interest
library(data.table)
library(dplyr)
library(coloc)
library(hash)

"%&%" = function(a,b) paste(a,b,sep="")

GWASdir = "/home/maya/COLOC/gwasFILES/META_Anx_coloc"
eQTLdir = "/home/maya/COLOC/eQTL_data/Hypothalamus/META_Anx_colocs/FAM167A_data.csv"

eQTLs = fread(eQTLdir) #read in significant gene eQTL
gwas = fread(GWASdir)

#get intersecting SNPs
joined = inner_join(eQTLs,gwas, by=c("ID"="ID"))

joined <- arrange(joined, P.y, by_group=FALSE)


#check for complement bases
#build hash table (like a Python dictionary)

bases = hash()
bases[["A"]] <- "T"
bases[["C"]] <- "G"
bases[["G"]] <- "C"
bases[["T"]] <- "A"


#pull SNPs that  match (assumes C|G and A|T SNPs aren't flipped)
match = joined[(joined$A1.x==joined$A1.y & joined$REF == joined$A2),]

#remove ambiguous strand SNPs (A/T or C/G)
a = joined[!(joined$REF=="A" & joined$ALT=="T") & !(joined$REF=="T" & joined$ALT=="A") & !(joined$REF=="C" & joined$ALT=="G") &
             !(joined$REF=="G" & joined$ALT=="C") ]

#of non-ambiguous, pull SNPs that are flipped (or the complement bases match or are flipped)

compmatch = a[(a$A1.x==values(bases,keys=a$A1.y) & a$REF == values(bases,keys=a$A2)),]

flipped = a[(a$A1.x==a$A2 & a$REF == a$A1.y),]
compflipped = a[(a$A1.x==values(bases,keys=a$A2) & a$REF == values(bases,keys=a$A1.y)),]

#if flipped, change sign of cancer beta, check to see if any SNPs are in flipped df's with if stmt
if(dim(flipped)[1] > 0){
  flipped = mutate(flipped,BETA.y = -1*BETA.y)
}
if(dim(compflipped)[1] > 0){
  compflipped = mutate(compflipped,BETA.y = -1*BETA.y)
}

#bind all and sort by position
matchsnps = rbind(match, compmatch, flipped, compflipped) %>% arrange(POS)

#matchsnps <- matchsnps %>%
#na.omit()

matchsnps <- matchsnps %>%
  distinct(BP, .keep_all = TRUE) %>%
  select(-n_cases)
matchsnps <- matchsnps %>%
  na.omit()
### 6. Format each set of summary stats for coloc
gwascoloc = list("beta" = matchsnps$BETA.y, "varbeta" = (matchsnps$SE.y)^2, "snp" = matchsnps$ID, "position" = matchsnps$BP,
                 "type" = "cc")
eqtlcoloc = list("beta" = matchsnps$BETA.x, "varbeta" = (matchsnps$SE.x)^2, "snp" = matchsnps$ID, "position" = matchsnps$POS,
                 "type" = "quant", "N" = matchsnps$OBS_CT[1], "MAF" = matchsnps$A1_FREQ, "sdY"=1 )


plot_dataset(gwascoloc)
plot_dataset(eqtlcoloc)

###7. Run coloc assuming a single causal variant
my.res = coloc.abf(dataset1=gwascoloc, dataset2=eqtlcoloc)
my.res
sensitivity(my.res,"H4 > 0.9")

################################################## RP11-981G7.6 - Amygdala ###########################################################################
### 5. Retrieve summary stat data from gene and phenotype of interest
library(data.table)
library(dplyr)
library(coloc)
library(hash)

"%&%" = function(a,b) paste(a,b,sep="")

GWASdir = "/home/maya/COLOC/gwasFILES/META_Anx_coloc"
eQTLdir = "/home/maya/COLOC/eQTL_data/Amygdala/META_Anx_colocs/RP11-981G7.6_data.csv"

eQTLs = fread(eQTLdir) #read in significant gene eQTL
gwas = fread(GWASdir)

#get intersecting SNPs
joined = inner_join(eQTLs,gwas, by=c("ID"="ID"))

joined <- arrange(joined, P.y, by_group=FALSE)


#check for complement bases
#build hash table (like a Python dictionary)

bases = hash()
bases[["A"]] <- "T"
bases[["C"]] <- "G"
bases[["G"]] <- "C"
bases[["T"]] <- "A"


#pull SNPs that  match (assumes C|G and A|T SNPs aren't flipped)
match = joined[(joined$A1.x==joined$A1.y & joined$REF == joined$A2),]

#remove ambiguous strand SNPs (A/T or C/G)
a = joined[!(joined$REF=="A" & joined$ALT=="T") & !(joined$REF=="T" & joined$ALT=="A") & !(joined$REF=="C" & joined$ALT=="G") &
             !(joined$REF=="G" & joined$ALT=="C") ]

#of non-ambiguous, pull SNPs that are flipped (or the complement bases match or are flipped)

compmatch = a[(a$A1.x==values(bases,keys=a$A1.y) & a$REF == values(bases,keys=a$A2)),]

flipped = a[(a$A1.x==a$A2 & a$REF == a$A1.y),]
compflipped = a[(a$A1.x==values(bases,keys=a$A2) & a$REF == values(bases,keys=a$A1.y)),]

#if flipped, change sign of cancer beta, check to see if any SNPs are in flipped df's with if stmt
if(dim(flipped)[1] > 0){
  flipped = mutate(flipped,BETA.y = -1*BETA.y)
}
if(dim(compflipped)[1] > 0){
  compflipped = mutate(compflipped,BETA.y = -1*BETA.y)
}

#bind all and sort by position
matchsnps = rbind(match, compmatch, flipped, compflipped) %>% arrange(POS)

#matchsnps <- matchsnps %>%
#na.omit()

matchsnps <- matchsnps %>%
  distinct(BP, .keep_all = TRUE) %>%
  select(-n_cases)
matchsnps <- matchsnps %>%
  na.omit()
### 6. Format each set of summary stats for coloc
gwascoloc = list("beta" = matchsnps$BETA.y, "varbeta" = (matchsnps$SE.y)^2, "snp" = matchsnps$ID, "position" = matchsnps$BP,
                 "type" = "cc")
eqtlcoloc = list("beta" = matchsnps$BETA.x, "varbeta" = (matchsnps$SE.x)^2, "snp" = matchsnps$ID, "position" = matchsnps$POS,
                 "type" = "quant", "N" = matchsnps$OBS_CT[1], "MAF" = matchsnps$A1_FREQ, "sdY"=1 )


plot_dataset(gwascoloc)
plot_dataset(eqtlcoloc)

###7. Run coloc assuming a single causal variant
my.res = coloc.abf(dataset1=gwascoloc, dataset2=eqtlcoloc)
my.res
sensitivity(my.res,"H4 > 0.9")
