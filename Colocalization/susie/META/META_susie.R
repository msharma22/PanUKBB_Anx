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

#################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################
################################################## CADM2 - Hypothalamus ######################################################################################################################################################
#################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################
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

#bind all and sort by position
matchsnps = rbind(match, compmatch, flipped, compflipped) %>% arrange(POS)


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

### 8. Add LD info to data frames for coloc.susie (multiple causal SNPs)
refdir = "/home/grace/coloc/missingness_hwe_steps/"
#make a list of EUR individuals in 1000 Genomes
pops = fread(refdir %&% "02geno_0.001_maf_0.05_filtered.fam")
eurlist = mutate(pops, V1=0) %>% dplyr::select(V1, V2)
fwrite(eurlist,"/home/maya/COLOC/EUR_list",col.names=FALSE, sep="\t")

#need chromosome code, start of range, end of range, set id for plink call
chr <- c(3)
first <- c(84458981)
##84958981 - 500K = 84458981
last <- c(86574429)
##86074429 + 500K = 86574429
id <- c('ENSG00000175161.13')
## chr3:84958981-86074429:+
# chr3: 84959546-86066835

range <- data.frame(chr, first, last, id)
fwrite(range, "/home/maya/COLOC/range.txt", col.names = FALSE, sep='\t')

#system call to plink to retrieve desired SNPs and individuals in .raw format
#genotypes must vary to be useful, maf>0.01

system("plink --bfile " %&% refdir %&% "02geno_0.001_maf_0.05_filtered --extract range /home/maya/COLOC/range.txt --keep /home/maya/COLOC/EUR_list --maf 0.01 --recode A --make-just-bim --out /home/maya/COLOC/EUR_coloc-region --allow-extra-chr")

#read the EUR .raw file generated by plink into R.
geno = fread("/home/maya/COLOC/EUR_coloc-region.raw")

#we want to put genotypes in NxP matrix (N=people, P=snps) with no other columns
genomat = as.matrix(geno[,7:length(geno)])
#now filter the gwascoloc and eqtlcoloc df's to those SNPs that were in EUR (in genomat)
#get a list of SNPs in genomat
snplist = colnames(genomat)
#remove the last 2 characters (_N) from genomat colnames to match topsleep with substr()
snplist = substr(snplist, 1, nchar(snplist)-2)
#rename col names of R to match coloc df's
colnames(genomat) = snplist

#check eQTL and 1000G EUR alleles (ask, what is coded as 1?), if need to, flip BETA signs.
#col 5 in bim is A1 (assigned 1 in dosage raw file)
bim = fread("/home/maya/COLOC/EUR_coloc-region.bim")

bim <- bim %>%
  mutate(ID = V1 %&% "_" %&% V4)

susiesnpsbim = inner_join(matchsnps, bim, by=c("ID"="ID"))

#pull SNPs that  match (assumes C|G and A|T SNPs aren't flipped)
matchbim = susiesnpsbim[(susiesnpsbim$A1.x==susiesnpsbim$V5 & susiesnpsbim$REF == susiesnpsbim$V6),]

#remove ambiguous strand SNPs (A/T or C/G)
b = susiesnpsbim[!(susiesnpsbim$REF=="A" & susiesnpsbim$ALT=="T") & !(susiesnpsbim$REF=="T" & susiesnpsbim$ALT=="A") & !(susiesnpsbim$REF=="C" & susiesnpsbim$ALT=="G") & !(susiesnpsbim$REF=="G" & susiesnpsbim$ALT=="C") ]

#of non-ambiguous, pull SNPs that are flipped (or the complement bases match or are flipped)
compmatchbim = b[(b$A1.x==values(bases,keys=b$V5) & b$REF == values(bases,keys=b$V6)),]
flippedbim = b[(b$A1.x==b$V6 & b$REF == b$V5),]
compflippedbim = b[(b$A1.x==values(bases,keys=b$V6) & b$REF == values(bases,keys=b$V5)),]

#if flipped, change sign of cancer and eqtl beta, check to see if any SNPs are in flipped df's with if stmt
if(dim(flippedbim)[1] > 0){
  flippedbim = mutate(flippedbim,BETA.x = -1*BETA.x,BETA.y = -1*BETA.y)
}
if(dim(compflippedbim)[1] > 0){
  compflippedbim = mutate(compflippedbim,BETA.x = -1*BETA.x,BETA.y = -1*BETA.y)
}

#bind all and sort by position
matchbimsnps = rbind(matchbim, compmatchbim, flippedbim, compflippedbim) %>% arrange(POS)

#update snplist
snplist = matchbimsnps$V2
#filter genomat to snps in snplist
x = genomat[,colnames(genomat) %in% snplist]
#calculate the correlation matrix of the genotypes, this is needed for susie
R = cor(x)

#filter the matchbimsnps df to just the snps in snplist with dplyr
susiesnps = filter(matchbimsnps, V2 %in% snplist)
#add LD to filtered coloc df's
gwascolocsusie = list("beta" = susiesnps$BETA.y, "varbeta" = (susiesnps$SE.y)^2, "snp" = susiesnps$V2, "position" = susiesnps$POS,
                      "type" = "cc", "LD"=R, "N" = 88) #need to add N, check ref later
eqtlcolocsusie = list("beta" = setNames(susiesnps$BETA.x, susiesnps$V2), "varbeta" = setNames(susiesnps$SE.x^2, susiesnps$V2), "snp" = susiesnps$V2,
                      "position" = susiesnps$POS,"type" = "quant", "N" = susiesnps$OBS_CT[1], "MAF" = susiesnps$A1_FREQ, "sdY"=1, "LD"= R)

check_dataset(gwascolocsusie,req="LD")

check_dataset(eqtlcolocsusie,req="LD")

plot_dataset(gwascolocsusie)

plot_dataset(eqtlcolocsusie)
###9. Run susie on each df and then run coloc
#run susie
sgwas = runsusie(gwascolocsusie)
seqtl = runsusie(eqtlcolocsusie)
#up to 10 credible sets by default
summary(sgwas)
summary(seqtl)
#run coloc
susie.res=coloc.susie(sgwas, seqtl)
print(susie.res$summary)
if(!is.na(susie.res)[1]){
  #retrive row numbers with H4>0.5
  sigrows = which(susie.res$summary$PP.H4.abf > 0.5)
  sigrows
  
  #sensitivity plot for row 1
  sensitivity(susie.res,"H4 > 0.5",row=1,dataset1=gwascolocsusie,dataset2=eqtlcolocsusie)
  
  #sensitivity plots for rows with H4>0.5
  for(i in sigrows){
    sensitivity(susie.res,"H4 > 0.5",row=i,dataset1=gwascolocsusie,dataset2=eqtlcolocsusie)
  }
}


###10. Run susie with lower coverage on each df and then run coloc
#run susie
sgwas = runsusie(gwascolocsusie, coverage=0.1)
seqtl = runsusie(eqtlcolocsusie, coverage=0.1)
#up to 10 credible sets by default
summary(sgwas)
summary(seqtl)
#run coloc
susie.res=coloc.susie(sgwas, seqtl)
print(susie.res$summary)
if(!is.na(susie.res)[1]){
  #retrieve row numbers with H4>0.4
  sigrows = which(susie.res$summary$PP.H4.abf > 0.4)
  sigrows
  
  #sensitivity plot for row 1
  sensitivity(susie.res,"H4 > 0.5",row=1,dataset1=gwascolocsusie,dataset2=eqtlcolocsusie)
  
  #sensitivity plots for rows with H4>0.4
  for(i in sigrows){
    sensitivity(susie.res,"H4 > 0.4",row=i,dataset1=gwascolocsusie,dataset2=eqtlcolocsusie)
  }
}
#################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################
################################################## CADM2 - Spinal Cord Cervical ######################################################################################################################################################
#################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################
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

### 8. Add LD info to data frames for coloc.susie (multiple causal SNPs)
refdir = "/home/grace/coloc/missingness_hwe_steps/"
#make a list of EUR individuals in 1000 Genomes
pops = fread(refdir %&% "02geno_0.001_maf_0.05_filtered.fam")
eurlist = mutate(pops, V1=0) %>% dplyr::select(V1, V2)
fwrite(eurlist,"/home/maya/COLOC/EUR_list",col.names=FALSE, sep="\t")

#need chromosome code, start of range, end of range, set id for plink call
chr <- c(3)
first <- c(84458981)
last <- c(86574429)
id <- c('ENSG00000175161.13')
#chr3:84958981-86074429:+
# chr3:84959546-86066835
range <- data.frame(chr, first, last, id)
fwrite(range, "/home/maya/COLOC/range.txt", col.names = FALSE, sep='\t')

#system call to plink to retrieve desired SNPs and individuals in .raw format
#genotypes must vary to be useful, maf>0.01

system("plink --bfile " %&% refdir %&% "02geno_0.001_maf_0.05_filtered --extract range /home/maya/COLOC/range.txt --keep /home/maya/COLOC/EUR_list --maf 0.01 --recode A --make-just-bim --out /home/maya/COLOC/EUR_coloc-region --allow-extra-chr")

#read the EUR .raw file generated by plink into R.
geno = fread("/home/maya/COLOC/EUR_coloc-region.raw")

#we want to put genotypes in NxP matrix (N=people, P=snps) with no other columns
genomat = as.matrix(geno[,7:length(geno)])
#now filter the gwascoloc and eqtlcoloc df's to those SNPs that were in EUR (in genomat)
#get a list of SNPs in genomat
snplist = colnames(genomat)
#remove the last 2 characters (_N) from genomat colnames to match topsleep with substr()
snplist = substr(snplist, 1, nchar(snplist)-2)
#rename col names of R to match coloc df's
colnames(genomat) = snplist

#check eQTL and 1000G EUR alleles (ask, what is coded as 1?), if need to, flip BETA signs.
#col 5 in bim is A1 (assigned 1 in dosage raw file)
bim = fread("/home/maya/COLOC/EUR_coloc-region.bim")

bim <- bim %>%
  mutate(ID = V1 %&% "_" %&% V4)

susiesnpsbim = inner_join(matchsnps, bim, by=c("ID"="ID"))

#pull SNPs that  match (assumes C|G and A|T SNPs aren't flipped)
matchbim = susiesnpsbim[(susiesnpsbim$A1.x==susiesnpsbim$V5 & susiesnpsbim$REF == susiesnpsbim$V6),]

#remove ambiguous strand SNPs (A/T or C/G)
b = susiesnpsbim[!(susiesnpsbim$REF=="A" & susiesnpsbim$ALT=="T") & !(susiesnpsbim$REF=="T" & susiesnpsbim$ALT=="A") & !(susiesnpsbim$REF=="C" & susiesnpsbim$ALT=="G") & !(susiesnpsbim$REF=="G" & susiesnpsbim$ALT=="C") ]

#of non-ambiguous, pull SNPs that are flipped (or the complement bases match or are flipped)
compmatchbim = b[(b$A1.x==values(bases,keys=b$V5) & b$REF == values(bases,keys=b$V6)),]
flippedbim = b[(b$A1.x==b$V6 & b$REF == b$V5),]
compflippedbim = b[(b$A1.x==values(bases,keys=b$V6) & b$REF == values(bases,keys=b$V5)),]

#if flipped, change sign of cancer and eqtl beta, check to see if any SNPs are in flipped df's with if stmt
if(dim(flippedbim)[1] > 0){
  flippedbim = mutate(flippedbim,BETA.x = -1*BETA.x,BETA.y = -1*BETA.y)
}
if(dim(compflippedbim)[1] > 0){
  compflippedbim = mutate(compflippedbim,BETA.x = -1*BETA.x,BETA.y = -1*BETA.y)
}

#bind all and sort by position
matchbimsnps = rbind(matchbim, compmatchbim, flippedbim, compflippedbim) %>% arrange(POS)

#update snplist
snplist = matchbimsnps$V2
#filter genomat to snps in snplist
x = genomat[,colnames(genomat) %in% snplist]
#calculate the correlation matrix of the genotypes, this is needed for susie
R = cor(x)

#filter the matchbimsnps df to just the snps in snplist with dplyr
susiesnps = filter(matchbimsnps, V2 %in% snplist)
#add LD to filtered coloc df's
gwascolocsusie = list("beta" = susiesnps$BETA.y, "varbeta" = (susiesnps$SE.y)^2, "snp" = susiesnps$V2, "position" = susiesnps$POS,
                      "type" = "cc", "LD"=R, "N" = 42) #need to add N, check ref later
eqtlcolocsusie = list("beta" = setNames(susiesnps$BETA.x, susiesnps$V2), "varbeta" = setNames(susiesnps$SE.x^2, susiesnps$V2), "snp" = susiesnps$V2,
                      "position" = susiesnps$POS,"type" = "quant", "N" = susiesnps$OBS_CT[1], "MAF" = susiesnps$A1_FREQ, "sdY"=1, "LD"= R)

check_dataset(gwascolocsusie,req="LD")

check_dataset(eqtlcolocsusie,req="LD")

plot_dataset(gwascolocsusie)

plot_dataset(eqtlcolocsusie)
###9. Run susie on each df and then run coloc
#run susie
sgwas = runsusie(gwascolocsusie)
seqtl = runsusie(eqtlcolocsusie)
#up to 10 credible sets by default
summary(sgwas)
summary(seqtl)
#run coloc
susie.res=coloc.susie(sgwas, seqtl)
print(susie.res$summary)
if(!is.na(susie.res)[1]){
  #retrive row numbers with H4>0.5
  sigrows = which(susie.res$summary$PP.H4.abf > 0.5)
  sigrows
  
  #sensitivity plot for row 1
  sensitivity(susie.res,"H4 > 0.5",row=1,dataset1=gwascolocsusie,dataset2=eqtlcolocsusie)
  
  #sensitivity plots for rows with H4>0.5
  for(i in sigrows){
    sensitivity(susie.res,"H4 > 0.5",row=i,dataset1=gwascolocsusie,dataset2=eqtlcolocsusie)
  }
}


###10. Run susie with lower coverage on each df and then run coloc
#run susie
sgwas = runsusie(gwascolocsusie, coverage=0.1)
seqtl = runsusie(eqtlcolocsusie, coverage=0.1)
#up to 10 credible sets by default
summary(sgwas)
summary(seqtl)
#run coloc
susie.res=coloc.susie(sgwas, seqtl)
print(susie.res$summary)
if(!is.na(susie.res)[1]){
  #retrieve row numbers with H4>0.4
  sigrows = which(susie.res$summary$PP.H4.abf > 0.4)
  sigrows
  
  #sensitivity plot for row 1
  sensitivity(susie.res,"H4 > 0.5",row=1,dataset1=gwascolocsusie,dataset2=eqtlcolocsusie)
  
  #sensitivity plots for rows with H4>0.4
  for(i in sigrows){
    sensitivity(susie.res,"H4 > 0.4",row=i,dataset1=gwascolocsusie,dataset2=eqtlcolocsusie)
  }
}
