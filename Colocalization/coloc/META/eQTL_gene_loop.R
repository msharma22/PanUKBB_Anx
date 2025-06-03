#making eQTL files for each gene
library(data.table)
library(dplyr)


##########################################################################AMYGDALA ##########################################################################- do for 12 other brain tisssues -- need to have genes.csv file made in each tisse directory 

"%&%" = function(a,b) paste(a,b,sep="")
tissue = 'Amygdala'

eQTL_Bigdir = "/home/grace/coloc/eQTL_data/" %&% tissue %&% "/" %&% "Brain_Amy_eQTL.csv" 
eQTLs <- fread(eQTL_Bigdir) #read in eQTLs from tissue 

popList = c("META_Anx_colocs")

#read in gene file for each 
for (population in popList) {
  new_pop = gsub('_coloc','',population)
  geneFile_path <- "/home/maya/panUKBB-COLOC/eQTL_data/" %&% tissue %&% "/" %&% population %&% "/genes.csv"
  geneFile <- fread(geneFile_path)
  data <- fread('/home/maya/panUKBB-prediXcan/prediXcanOUTPUT/ANX/META/META_Anx_' %&% tissue %&% '.csv')
  if ("gene_name" %in% colnames(geneFile)){
    join <- right_join(data,geneFile, by = "gene_name")
    gene_ids <- as.list(join$gene)
    genes <- as.list(join$gene_name)
  } else{
    join <- right_join(data,geneFile, by = c("gene_name"="gene_name.x"))
    gene_ids <- as.list(join$gene)
    genes <- as.list(join$gene_name)
  }
  for (k in 1:length(genes)) {
    geneID <- toString(gene_ids[k])
    gene <- toString(genes[k])
    
    new_eQTLs <- eQTLs %>%
      filter(gene_id == geneID) %>% #only use snps with ensembl id of the significant gene
      arrange(P, by_group = FALSE) #arrange by p value
    name <- "/home/maya/panUKBB-COLOC/eQTL_data/" %&% tissue %&% "/" %&% population %&% "/" %&% "/" %&% gene %&% "_data.csv"
    fwrite(new_eQTLs, name, append=FALSE,quote = "auto", sep=",") #write new eQTL file
  }
}
########################################################################## Anterior_Cingulate_Cortex ##########################################################################- do for 12 other brain tisssues -- need to have genes.csv file made in each tisse directory 

"%&%" = function(a,b) paste(a,b,sep="")
tissue = 'Anterior_cingulate_cortex'

eQTL_Bigdir = "/home/grace/coloc/eQTL_data/" %&% tissue %&% "/" %&% "Brain_ACC_eQTL.csv" 
eQTLs <- fread(eQTL_Bigdir) #read in eQTLs from tissue 

popList = c("META_Anx_colocs")

#read in gene file for each 
for (population in popList) {
  new_pop = gsub('_coloc','',population)
  geneFile_path <- "/home/maya/COLOC/eQTL_data/" %&% tissue %&% "/" %&% population %&% "/genes.csv"
  geneFile <- fread(geneFile_path)
  data <- fread('/home/maya/prediXcan/prediXcanOUTPUT/ANX/META/META_Anx_Anterior_Cingulate_Cortex_BA24.csv')
  if ("gene_name" %in% colnames(geneFile)){
    join <- right_join(data,geneFile, by = "gene_name")
    gene_ids <- as.list(join$gene)
    genes <- as.list(join$gene_name)
  } else{
    join <- right_join(data,geneFile, by = c("gene_name"="gene_name.x"))
    gene_ids <- as.list(join$gene)
    genes <- as.list(join$gene_name)
  }
  for (k in 1:length(genes)) {
    geneID <- toString(gene_ids[k])
    gene <- toString(genes[k])
    
    new_eQTLs <- eQTLs %>%
      filter(gene_id == geneID) %>% #only use snps with ensembl id of the significant gene
      arrange(P, by_group = FALSE) #arrange by p value
    name <- "/home/maya/COLOC/eQTL_data/" %&% tissue %&% "/" %&% population %&% "/" %&% "/" %&% gene %&% "_data.csv"
    fwrite(new_eQTLs, name, append=FALSE,quote = "auto", sep=",") #write new eQTL file
  }
}

########################################################################## Caudate_basal_ganglia ##########################################################################- do for 12 other brain tisssues -- need to have genes.csv file made in each tisse directory 

"%&%" = function(a,b) paste(a,b,sep="")
tissue = 'Caudate_basal_ganglia'

eQTL_Bigdir = "/home/grace/coloc/eQTL_data/" %&% tissue %&% "/" %&% "Brain_CBC_eQTL.csv" 
eQTLs <- fread(eQTL_Bigdir) #read in eQTLs from tissue 

popList = c("META_Anx_colocs")

#read in gene file for each 
for (population in popList) {
  new_pop = gsub('_coloc','',population)
  geneFile_path <- "/home/maya/COLOC/eQTL_data/" %&% tissue %&% "/" %&% population %&% "/genes.csv"
  geneFile <- fread(geneFile_path)
  data <- fread('/home/maya/prediXcan/prediXcanOUTPUT/ANX/META/META_Anx_Caudate_Basal_Ganglia.csv')
  if ("gene_name" %in% colnames(geneFile)){
    join <- right_join(data,geneFile, by = "gene_name")
    gene_ids <- as.list(join$gene)
    genes <- as.list(join$gene_name)
  } else{
    join <- right_join(data,geneFile, by = c("gene_name"="gene_name.x"))
    gene_ids <- as.list(join$gene)
    genes <- as.list(join$gene_name)
  }
  for (k in 1:length(genes)) {
    geneID <- toString(gene_ids[k])
    gene <- toString(genes[k])
    
    new_eQTLs <- eQTLs %>%
      filter(gene_id == geneID) %>% #only use snps with ensembl id of the significant gene
      arrange(P, by_group = FALSE) #arrange by p value
    name <- "/home/maya/COLOC/eQTL_data/" %&% tissue %&% "/" %&% population %&% "/" %&% "/" %&% gene %&% "_data.csv"
    fwrite(new_eQTLs, name, append=FALSE,quote = "auto", sep=",") #write new eQTL file
  }
}

########################################################################## Caudate_basal_ganglia ##########################################################################- do for 12 other brain tisssues -- need to have genes.csv file made in each tisse directory 

"%&%" = function(a,b) paste(a,b,sep="")
tissue = 'Caudate_basal_ganglia'

eQTL_Bigdir = "/home/grace/coloc/eQTL_data/" %&% tissue %&% "/" %&% "Brain_CBC_eQTL.csv" 
eQTLs <- fread(eQTL_Bigdir) #read in eQTLs from tissue 

popList = c("META_Anx_colocs")

#read in gene file for each 
for (population in popList) {
  new_pop = gsub('_coloc','',population)
  geneFile_path <- "/home/maya/COLOC/eQTL_data/" %&% tissue %&% "/" %&% population %&% "/genes.csv"
  geneFile <- fread(geneFile_path)
  data <- fread('/home/maya/prediXcan/prediXcanOUTPUT/ANX/META/META_Anx_Caudate_Basal_Ganglia.csv')
  if ("gene_name" %in% colnames(geneFile)){
    join <- right_join(data,geneFile, by = "gene_name")
    gene_ids <- as.list(join$gene)
    genes <- as.list(join$gene_name)
  } else{
    join <- right_join(data,geneFile, by = c("gene_name"="gene_name.x"))
    gene_ids <- as.list(join$gene)
    genes <- as.list(join$gene_name)
  }
  for (k in 1:length(genes)) {
    geneID <- toString(gene_ids[k])
    gene <- toString(genes[k])
    
    new_eQTLs <- eQTLs %>%
      filter(gene_id == geneID) %>% #only use snps with ensembl id of the significant gene
      arrange(P, by_group = FALSE) #arrange by p value
    name <- "/home/maya/COLOC/eQTL_data/" %&% tissue %&% "/" %&% population %&% "/" %&% "/" %&% gene %&% "_data.csv"
    fwrite(new_eQTLs, name, append=FALSE,quote = "auto", sep=",") #write new eQTL file
  }
}
########################################################################## Cerebellar_Hemisphere ##########################################################################- do for 12 other brain tisssues -- need to have genes.csv file made in each tisse directory 

"%&%" = function(a,b) paste(a,b,sep="")
tissue = 'Cerebellar_Hemisphere'

eQTL_Bigdir = "/home/grace/coloc/eQTL_data/" %&% tissue %&% "/" %&% "Brain_CH_eQTL.csv" 
eQTLs <- fread(eQTL_Bigdir) #read in eQTLs from tissue 

popList = c("META_Anx_colocs")

#read in gene file for each 
for (population in popList) {
  new_pop = gsub('_coloc','',population)
  geneFile_path <- "/home/maya/COLOC/eQTL_data/" %&% tissue %&% "/" %&% population %&% "/genes.csv"
  geneFile <- fread(geneFile_path)
  data <- fread('/home/maya/prediXcan/prediXcanOUTPUT/ANX/META/META_Anx_Cerebellar_Hemisphere.csv')
  if ("gene_name" %in% colnames(geneFile)){
    join <- right_join(data,geneFile, by = "gene_name")
    gene_ids <- as.list(join$gene)
    genes <- as.list(join$gene_name)
  } else{
    join <- right_join(data,geneFile, by = c("gene_name"="gene_name.x"))
    gene_ids <- as.list(join$gene)
    genes <- as.list(join$gene_name)
  }
  for (k in 1:length(genes)) {
    geneID <- toString(gene_ids[k])
    gene <- toString(genes[k])
    
    new_eQTLs <- eQTLs %>%
      filter(gene_id == geneID) %>% #only use snps with ensembl id of the significant gene
      arrange(P, by_group = FALSE) #arrange by p value
    name <- "/home/maya/COLOC/eQTL_data/" %&% tissue %&% "/" %&% population %&% "/" %&% "/" %&% gene %&% "_data.csv"
    fwrite(new_eQTLs, name, append=FALSE,quote = "auto", sep=",") #write new eQTL file
  }
}
########################################################################## Cerebellar_Hemisphere ##########################################################################- do for 12 other brain tisssues -- need to have genes.csv file made in each tisse directory 
"%&%" = function(a,b) paste(a,b,sep="")
tissue = 'Cerebellar_Hemisphere'

eQTL_Bigdir = "/home/grace/coloc/eQTL_data/" %&% tissue %&% "/" %&% "Brain_CH_eQTL.csv" 
eQTLs <- fread(eQTL_Bigdir) #read in eQTLs from tissue 

popList = c("META_Anx_colocs")

#read in gene file for each 
for (population in popList) {
  new_pop = gsub('_coloc','',population)
  geneFile_path <- "/home/maya/COLOC/eQTL_data/" %&% tissue %&% "/" %&% population %&% "/genes.csv"
  geneFile <- fread(geneFile_path)
  data <- fread('/home/maya/prediXcan/prediXcanOUTPUT/ANX/META/META_Anx_Cerebellar_Hemisphere.csv')
  if ("gene_name" %in% colnames(geneFile)){
    join <- right_join(data,geneFile, by = "gene_name")
    gene_ids <- as.list(join$gene)
    genes <- as.list(join$gene_name)
  } else{
    join <- right_join(data,geneFile, by = c("gene_name"="gene_name.x"))
    gene_ids <- as.list(join$gene)
    genes <- as.list(join$gene_name)
  }
  for (k in 1:length(genes)) {
    geneID <- toString(gene_ids[k])
    gene <- toString(genes[k])
    
    new_eQTLs <- eQTLs %>%
      filter(gene_id == geneID) %>% #only use snps with ensembl id of the significant gene
      arrange(P, by_group = FALSE) #arrange by p value
    name <- "/home/maya/COLOC/eQTL_data/" %&% tissue %&% "/" %&% population %&% "/" %&% "/" %&% gene %&% "_data.csv"
    fwrite(new_eQTLs, name, append=FALSE,quote = "auto", sep=",") #write new eQTL file
  }
}
########################################################################## Cerebellar_Hemisphere ##########################################################################- do for 12 other brain tisssues -- need to have genes.csv file made in each tisse directory 

"%&%" = function(a,b) paste(a,b,sep="")
tissue = 'Cerebellar_Hemisphere'

eQTL_Bigdir = "/home/grace/coloc/eQTL_data/" %&% tissue %&% "/" %&% "Brain_CH_eQTL.csv" 
eQTLs <- fread(eQTL_Bigdir) #read in eQTLs from tissue 

popList = c("META_Anx_colocs")

#read in gene file for each 
for (population in popList) {
  new_pop = gsub('_coloc','',population)
  geneFile_path <- "/home/maya/COLOC/eQTL_data/" %&% tissue %&% "/" %&% population %&% "/genes.csv"
  geneFile <- fread(geneFile_path)
  data <- fread('/home/maya/prediXcan/prediXcanOUTPUT/ANX/META/META_Anx_Cerebellar_Hemisphere.csv')
  if ("gene_name" %in% colnames(geneFile)){
    join <- right_join(data,geneFile, by = "gene_name")
    gene_ids <- as.list(join$gene)
    genes <- as.list(join$gene_name)
  } else{
    join <- right_join(data,geneFile, by = c("gene_name"="gene_name.x"))
    gene_ids <- as.list(join$gene)
    genes <- as.list(join$gene_name)
  }
  for (k in 1:length(genes)) {
    geneID <- toString(gene_ids[k])
    gene <- toString(genes[k])
    
    new_eQTLs <- eQTLs %>%
      filter(gene_id == geneID) %>% #only use snps with ensembl id of the significant gene
      arrange(P, by_group = FALSE) #arrange by p value
    name <- "/home/maya/COLOC/eQTL_data/" %&% tissue %&% "/" %&% population %&% "/" %&% "/" %&% gene %&% "_data.csv"
    fwrite(new_eQTLs, name, append=FALSE,quote = "auto", sep=",") #write new eQTL file
  }
}
########################################################################## Cerebellar_Hemisphere ##########################################################################- do for 12 other brain tisssues -- need to have genes.csv file made in each tisse directory 
"%&%" = function(a,b) paste(a,b,sep="")
tissue = 'Cerebellar_Hemisphere'

eQTL_Bigdir = "/home/grace/coloc/eQTL_data/" %&% tissue %&% "/" %&% "Brain_CH_eQTL.csv" 
eQTLs <- fread(eQTL_Bigdir) #read in eQTLs from tissue 

popList = c("META_Anx_colocs")

#read in gene file for each 
for (population in popList) {
  new_pop = gsub('_coloc','',population)
  geneFile_path <- "/home/maya/COLOC/eQTL_data/" %&% tissue %&% "/" %&% population %&% "/genes.csv"
  geneFile <- fread(geneFile_path)
  data <- fread('/home/maya/prediXcan/prediXcanOUTPUT/ANX/META/META_Anx_Cerebellar_Hemisphere.csv')
  if ("gene_name" %in% colnames(geneFile)){
    join <- right_join(data,geneFile, by = "gene_name")
    gene_ids <- as.list(join$gene)
    genes <- as.list(join$gene_name)
  } else{
    join <- right_join(data,geneFile, by = c("gene_name"="gene_name.x"))
    gene_ids <- as.list(join$gene)
    genes <- as.list(join$gene_name)
  }
  for (k in 1:length(genes)) {
    geneID <- toString(gene_ids[k])
    gene <- toString(genes[k])
    
    new_eQTLs <- eQTLs %>%
      filter(gene_id == geneID) %>% #only use snps with ensembl id of the significant gene
      arrange(P, by_group = FALSE) #arrange by p value
    name <- "/home/maya/COLOC/eQTL_data/" %&% tissue %&% "/" %&% population %&% "/" %&% "/" %&% gene %&% "_data.csv"
    fwrite(new_eQTLs, name, append=FALSE,quote = "auto", sep=",") #write new eQTL file
  }
}
########################################################################## Cerebellum ##########################################################################- do for 12 other brain tisssues -- need to have genes.csv file made in each tisse directory 
"%&%" = function(a,b) paste(a,b,sep="")
tissue = 'Cerebellum'

eQTL_Bigdir = "/home/grace/coloc/eQTL_data/" %&% tissue %&% "/" %&% "Brain_Cerebellum_eQTL.csv" 
eQTLs <- fread(eQTL_Bigdir) #read in eQTLs from tissue 

popList = c("META_Anx_colocs")

#read in gene file for each 
for (population in popList) {
  new_pop = gsub('_coloc','',population)
  geneFile_path <- "/home/maya/COLOC/eQTL_data/" %&% tissue %&% "/" %&% population %&% "/genes.csv"
  geneFile <- fread(geneFile_path)
  data <- fread('/home/maya/prediXcan/prediXcanOUTPUT/ANX/META/META_Anx_Cerebellum.csv')
  if ("gene_name" %in% colnames(geneFile)){
    join <- right_join(data,geneFile, by = "gene_name")
    gene_ids <- as.list(join$gene)
    genes <- as.list(join$gene_name)
  } else{
    join <- right_join(data,geneFile, by = c("gene_name"="gene_name.x"))
    gene_ids <- as.list(join$gene)
    genes <- as.list(join$gene_name)
  }
  for (k in 1:length(genes)) {
    geneID <- toString(gene_ids[k])
    gene <- toString(genes[k])
    
    new_eQTLs <- eQTLs %>%
      filter(gene_id == geneID) %>% #only use snps with ensembl id of the significant gene
      arrange(P, by_group = FALSE) #arrange by p value
    name <- "/home/maya/COLOC/eQTL_data/" %&% tissue %&% "/" %&% population %&% "/" %&% "/" %&% gene %&% "_data.csv"
    fwrite(new_eQTLs, name, append=FALSE,quote = "auto", sep=",") #write new eQTL file
  }
}
##########################################################################  Cortex  ##########################################################################- do for 12 other brain tisssues -- need to have genes.csv file made in each tisse directory 
"%&%" = function(a,b) paste(a,b,sep="")
tissue = 'Cortex'

eQTL_Bigdir = "/home/grace/coloc/eQTL_data/" %&% tissue %&% "/" %&% "Brain_Cortex_eQTL.csv" 
eQTLs <- fread(eQTL_Bigdir) #read in eQTLs from tissue 

popList = c("META_Anx_colocs")

#read in gene file for each 
for (population in popList) {
  new_pop = gsub('_coloc','',population)
  geneFile_path <- "/home/maya/COLOC/eQTL_data/" %&% tissue %&% "/" %&% population %&% "/genes.csv"
  geneFile <- fread(geneFile_path)
  data <- fread('/home/maya/prediXcan/prediXcanOUTPUT/ANX/META/META_Anx_Cortex.csv')
  if ("gene_name" %in% colnames(geneFile)){
    join <- right_join(data,geneFile, by = "gene_name")
    gene_ids <- as.list(join$gene)
    genes <- as.list(join$gene_name)
  } else{
    join <- right_join(data,geneFile, by = c("gene_name"="gene_name.x"))
    gene_ids <- as.list(join$gene)
    genes <- as.list(join$gene_name)
  }
  for (k in 1:length(genes)) {
    geneID <- toString(gene_ids[k])
    gene <- toString(genes[k])
    
    new_eQTLs <- eQTLs %>%
      filter(gene_id == geneID) %>% #only use snps with ensembl id of the significant gene
      arrange(P, by_group = FALSE) #arrange by p value
    name <- "/home/maya/COLOC/eQTL_data/" %&% tissue %&% "/" %&% population %&% "/" %&% "/" %&% gene %&% "_data.csv"
    fwrite(new_eQTLs, name, append=FALSE,quote = "auto", sep=",") #write new eQTL file
  }
}
##########################################################################  Frontal_Cortex  ##########################################################################- do for 12 other brain tisssues -- need to have genes.csv file made in each tisse directory 
"%&%" = function(a,b) paste(a,b,sep="")
tissue = 'Frontal_Cortex'

eQTL_Bigdir = "/home/grace/coloc/eQTL_data/" %&% tissue %&% "/" %&% "Brain_Frontal_Cortex_eQTL.csv" 
eQTLs <- fread(eQTL_Bigdir) #read in eQTLs from tissue 

popList = c("META_Anx_colocs")

#read in gene file for each 
for (population in popList) {
  new_pop = gsub('_coloc','',population)
  geneFile_path <- "/home/maya/COLOC/eQTL_data/" %&% tissue %&% "/" %&% population %&% "/genes.csv"
  geneFile <- fread(geneFile_path)
  data <- fread('/home/maya/prediXcan/prediXcanOUTPUT/ANX/META/META_Anx_Frontal_Cortex_BA9.csv')
  if ("gene_name" %in% colnames(geneFile)){
    join <- right_join(data,geneFile, by = "gene_name")
    gene_ids <- as.list(join$gene)
    genes <- as.list(join$gene_name)
  } else{
    join <- right_join(data,geneFile, by = c("gene_name"="gene_name.x"))
    gene_ids <- as.list(join$gene)
    genes <- as.list(join$gene_name)
  }
  for (k in 1:length(genes)) {
    geneID <- toString(gene_ids[k])
    gene <- toString(genes[k])
    
    new_eQTLs <- eQTLs %>%
      filter(gene_id == geneID) %>% #only use snps with ensembl id of the significant gene
      arrange(P, by_group = FALSE) #arrange by p value
    name <- "/home/maya/COLOC/eQTL_data/" %&% tissue %&% "/" %&% population %&% "/" %&% "/" %&% gene %&% "_data.csv"
    fwrite(new_eQTLs, name, append=FALSE,quote = "auto", sep=",") #write new eQTL file
  }
}
##########################################################################  Hippocampus  ##########################################################################- do for 12 other brain tisssues -- need to have genes.csv file made in each tisse directory 
"%&%" = function(a,b) paste(a,b,sep="")
tissue = 'Hippocampus'

eQTL_Bigdir = "/home/grace/coloc/eQTL_data/" %&% tissue %&% "/" %&% "Brain_Hippocampus_eQTL.csv" 
eQTLs <- fread(eQTL_Bigdir) #read in eQTLs from tissue 

popList = c("META_Anx_colocs")

#read in gene file for each 
for (population in popList) {
  new_pop = gsub('_coloc','',population)
  geneFile_path <- "/home/maya/COLOC/eQTL_data/" %&% tissue %&% "/" %&% population %&% "/genes.csv"
  geneFile <- fread(geneFile_path)
  data <- fread('/home/maya/prediXcan/prediXcanOUTPUT/ANX/META/META_Anx_Hippocampus.csv')
  if ("gene_name" %in% colnames(geneFile)){
    join <- right_join(data,geneFile, by = "gene_name")
    gene_ids <- as.list(join$gene)
    genes <- as.list(join$gene_name)
  } else{
    join <- right_join(data,geneFile, by = c("gene_name"="gene_name.x"))
    gene_ids <- as.list(join$gene)
    genes <- as.list(join$gene_name)
  }
  for (k in 1:length(genes)) {
    geneID <- toString(gene_ids[k])
    gene <- toString(genes[k])
    
    new_eQTLs <- eQTLs %>%
      filter(gene_id == geneID) %>% #only use snps with ensembl id of the significant gene
      arrange(P, by_group = FALSE) #arrange by p value
    name <- "/home/maya/COLOC/eQTL_data/" %&% tissue %&% "/" %&% population %&% "/" %&% "/" %&% gene %&% "_data.csv"
    fwrite(new_eQTLs, name, append=FALSE,quote = "auto", sep=",") #write new eQTL file
  }
}
##########################################################################  Hypothalamus  ##########################################################################- do for 12 other brain tisssues -- need to have genes.csv file made in each tisse directory 
"%&%" = function(a,b) paste(a,b,sep="")
tissue = 'Hypothalamus'

eQTL_Bigdir = "/home/grace/coloc/eQTL_data/" %&% tissue %&% "/" %&% "Brain_Hypo_eQTL.csv" 
eQTLs <- fread(eQTL_Bigdir) #read in eQTLs from tissue 

popList = c("META_Anx_colocs")

#read in gene file for each 
for (population in popList) {
  new_pop = gsub('_coloc','',population)
  geneFile_path <- "/home/maya/COLOC/eQTL_data/" %&% tissue %&% "/" %&% population %&% "/genes.csv"
  geneFile <- fread(geneFile_path)
  data <- fread('/home/maya/prediXcan/prediXcanOUTPUT/ANX/META/META_Anx_Hippocampus.csv')
  if ("gene_name" %in% colnames(geneFile)){
    join <- right_join(data,geneFile, by = "gene_name")
    gene_ids <- as.list(join$gene)
    genes <- as.list(join$gene_name)
  } else{
    join <- right_join(data,geneFile, by = c("gene_name"="gene_name.x"))
    gene_ids <- as.list(join$gene)
    genes <- as.list(join$gene_name)
  }
  for (k in 1:length(genes)) {
    geneID <- toString(gene_ids[k])
    gene <- toString(genes[k])
    
    new_eQTLs <- eQTLs %>%
      filter(gene_id == geneID) %>% #only use snps with ensembl id of the significant gene
      arrange(P, by_group = FALSE) #arrange by p value
    name <- "/home/maya/COLOC/eQTL_data/" %&% tissue %&% "/" %&% population %&% "/" %&% "/" %&% gene %&% "_data.csv"
    fwrite(new_eQTLs, name, append=FALSE,quote = "auto", sep=",") #write new eQTL file
  }
}
##########################################################################  Nucleus_accumbens_basal_ganglia  ##########################################################################- do for 12 other brain tisssues -- need to have genes.csv file made in each tisse directory 
"%&%" = function(a,b) paste(a,b,sep="")
tissue = 'Nucleus_accumbens_basal_ganglia'

eQTL_Bigdir = "/home/grace/coloc/eQTL_data/" %&% tissue %&% "/" %&% "Brain_NABG_eQTL.csv" 
eQTLs <- fread(eQTL_Bigdir) #read in eQTLs from tissue 

popList = c("META_Anx_colocs")

#read in gene file for each 
for (population in popList) {
  new_pop = gsub('_coloc','',population)
  geneFile_path <- "/home/maya/COLOC/eQTL_data/" %&% tissue %&% "/" %&% population %&% "/genes.csv"
  geneFile <- fread(geneFile_path)
  data <- fread('/home/maya/prediXcan/prediXcanOUTPUT/ANX/META/META_Anx_Hippocampus.csv')
  if ("gene_name" %in% colnames(geneFile)){
    join <- right_join(data,geneFile, by = "gene_name")
    gene_ids <- as.list(join$gene)
    genes <- as.list(join$gene_name)
  } else{
    join <- right_join(data,geneFile, by = c("gene_name"="gene_name.x"))
    gene_ids <- as.list(join$gene)
    genes <- as.list(join$gene_name)
  }
  for (k in 1:length(genes)) {
    geneID <- toString(gene_ids[k])
    gene <- toString(genes[k])
    
    new_eQTLs <- eQTLs %>%
      filter(gene_id == geneID) %>% #only use snps with ensembl id of the significant gene
      arrange(P, by_group = FALSE) #arrange by p value
    name <- "/home/maya/COLOC/eQTL_data/" %&% tissue %&% "/" %&% population %&% "/" %&% "/" %&% gene %&% "_data.csv"
    fwrite(new_eQTLs, name, append=FALSE,quote = "auto", sep=",") #write new eQTL file
  }
}
########################################################################## Putamen_basal_ganglia  ##########################################################################- do for 12 other brain tisssues -- need to have genes.csv file made in each tisse directory 
"%&%" = function(a,b) paste(a,b,sep="")
tissue = 'Putamen_basal_ganglia'

eQTL_Bigdir = "/home/grace/coloc/eQTL_data/" %&% tissue %&% "/" %&% "Brain_PBG_eQTL.csv" 
eQTLs <- fread(eQTL_Bigdir) #read in eQTLs from tissue 

popList = c("META_Anx_colocs")

#read in gene file for each 
for (population in popList) {
  new_pop = gsub('_coloc','',population)
  geneFile_path <- "/home/maya/COLOC/eQTL_data/" %&% tissue %&% "/" %&% population %&% "/genes.csv"
  geneFile <- fread(geneFile_path)
  data <- fread('/home/maya/prediXcan/prediXcanOUTPUT/ANX/META/META_Anx_Putamen_Basal_Ganglia.csv')
  if ("gene_name" %in% colnames(geneFile)){
    join <- right_join(data,geneFile, by = "gene_name")
    gene_ids <- as.list(join$gene)
    genes <- as.list(join$gene_name)
  } else{
    join <- right_join(data,geneFile, by = c("gene_name"="gene_name.x"))
    gene_ids <- as.list(join$gene)
    genes <- as.list(join$gene_name)
  }
  for (k in 1:length(genes)) {
    geneID <- toString(gene_ids[k])
    gene <- toString(genes[k])
    
    new_eQTLs <- eQTLs %>%
      filter(gene_id == geneID) %>% #only use snps with ensembl id of the significant gene
      arrange(P, by_group = FALSE) #arrange by p value
    name <- "/home/maya/COLOC/eQTL_data/" %&% tissue %&% "/" %&% population %&% "/" %&% "/" %&% gene %&% "_data.csv"
    fwrite(new_eQTLs, name, append=FALSE,quote = "auto", sep=",") #write new eQTL file
  }
}
########################################################################## Spinal_cord_cervical  ##########################################################################- do for 12 other brain tisssues -- need to have genes.csv file made in each tisse directory 
"%&%" = function(a,b) paste(a,b,sep="")
tissue = 'Spinal_cord_cervical'

eQTL_Bigdir = "/home/grace/coloc/eQTL_data/" %&% tissue %&% "/" %&% "Brain_SCC_eQTL.csv" 
eQTLs <- fread(eQTL_Bigdir) #read in eQTLs from tissue 

popList = c("META_Anx_colocs")

population = c("META_Anx_colocs")

#read in gene file for each 
for (population in popList) {
  new_pop = gsub('_coloc','',population)
  geneFile_path <- "/home/maya/COLOC/eQTL_data/" %&% tissue %&% "/" %&% population %&% "/genes.csv"
  geneFile <- fread(geneFile_path)
  data <- fread('/home/maya/prediXcan/prediXcanOUTPUT/ANX/META/META_Anx_Spinal_Cord_Cervical_c-1.csv')
  if ("gene_name" %in% colnames(geneFile)){
    join <- right_join(data,geneFile, by = "gene_name")
    gene_ids <- as.list(join$gene)
    genes <- as.list(join$gene_name)
  } else{
    join <- right_join(data,geneFile, by = c("gene_name"="gene_name.x"))
    gene_ids <- as.list(join$gene)
    genes <- as.list(join$gene_name)
  }
  for (k in 1:length(genes)) {
    geneID <- toString(gene_ids[k])
    gene <- toString(genes[k])
    
    new_eQTLs <- eQTLs %>%
      filter(gene_id == geneID) %>% #only use snps with ensembl id of the significant gene
      arrange(P, by_group = FALSE) #arrange by p value
    name <- "/home/maya/COLOC/eQTL_data/" %&% tissue %&% "/" %&% population %&% "/" %&% "/" %&% gene %&% "_data.csv"
    fwrite(new_eQTLs, name, append=FALSE,quote = "auto", sep=",") #write new eQTL file
  }
}

########################################################################## Substantia_nigra  ##########################################################################- do for 12 other brain tisssues -- need to have genes.csv file made in each tisse directory 
"%&%" = function(a,b) paste(a,b,sep="")
tissue = 'Substantia_nigra'

eQTL_Bigdir = "/home/grace/coloc/eQTL_data/" %&% tissue %&% "/" %&% "Brain_SN_eQTL.csv" 
eQTLs <- fread(eQTL_Bigdir) #read in eQTLs from tissue 

popList = c("META_Anx_colocs")

#read in gene file for each 
for (population in popList) {
  new_pop = gsub('_coloc','',population)
  geneFile_path <- "/home/maya/COLOC/eQTL_data/" %&% tissue %&% "/" %&% population %&% "/genes.csv"
  geneFile <- fread(geneFile_path)
  data <- fread('/home/maya/prediXcan/prediXcanOUTPUT/ANX/META/META_Anx_Substantia_nigra.csv')
  if ("gene_name" %in% colnames(geneFile)){
    join <- right_join(data,geneFile, by = "gene_name")
    gene_ids <- as.list(join$gene)
    genes <- as.list(join$gene_name)
  } else{
    join <- right_join(data,geneFile, by = c("gene_name"="gene_name.x"))
    gene_ids <- as.list(join$gene)
    genes <- as.list(join$gene_name)
  }
  for (k in 1:length(genes)) {
    geneID <- toString(gene_ids[k])
    gene <- toString(genes[k])
    
    new_eQTLs <- eQTLs %>%
      filter(gene_id == geneID) %>% #only use snps with ensembl id of the significant gene
      arrange(P, by_group = FALSE) #arrange by p value
    name <- "/home/maya/COLOC/eQTL_data/" %&% tissue %&% "/" %&% population %&% "/" %&% "/" %&% gene %&% "_data.csv"
    fwrite(new_eQTLs, name, append=FALSE,quote = "auto", sep=",") #write new eQTL file
  }
}
