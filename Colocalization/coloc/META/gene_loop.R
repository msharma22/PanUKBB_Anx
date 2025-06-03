################################################################################ AMYGDALA ################################################################################IS DONE -- do for 12 other brain tisssues -- need to have genes.csv file made in each tisse directory 

library(data.table)
library(dplyr)

"%&%" = function(a,b) paste(a,b,sep="")

tissueList = c("Amygdala")
anno <- fread("/home/chris/gencode/gene_annotation_v38_strand.txt", header = T, sep=' ')

popList = c("META_Anx_colocs")

for (population in popList) {
  for (tissue in tissueList) {
    new_pop = gsub('_coloc','',population)
    data <- fread('/home/maya/prediXcan/prediXcanOUTPUT/ANX/META/META_Anx_' %&% tissue %&% '.csv')
    new_data <- data %>%
      mutate(new_gene = sapply(strsplit(gene, "[.]"), `[`, 1)) %>%
      mutate(pop = new_pop)
    
    genes <- new_data %>%
      select(new_gene, gene_name, pvalue, pop)
    
    genes <- genes %>%
      filter(pvalue <= 5e-8)
    join <- left_join(genes,anno, by = c("new_gene"="gene_id"))
    join <- join %>%
      rename(gene_name = gene_name.x) %>%
      select(-c(strand,gene_type,gene_name.y))
    path = "/home/maya/COLOC/eQTL_data/" %&% tissue %&% "/" %&% population %&% "/genes.csv"
    fwrite(join, path)
  }
  
}

#################################################### Anterior_cingulate_cortex ################################################################################IS DONE -- do for 12 other brain tisssues -- need to have genes.csv file made in each tisse directory 

tissueList = c("Anterior_cingulate_cortex")
anno <- fread("/home/chris/gencode/gene_annotation_v38_strand.txt", header = T, sep=' ')

popList = c("META_Anx_colocs")

for (population in popList) {
  for (tissue in tissueList) {
    new_pop = gsub('_coloc','',population)
    data <- fread('/home/maya/prediXcan/prediXcanOUTPUT/ANX/META/META_Anx_Anterior_Cingulate_Cortex_BA24.csv')
    new_data <- data %>%
      mutate(new_gene = sapply(strsplit(gene, "[.]"), `[`, 1)) %>%
      mutate(pop = new_pop)
    
    genes <- new_data %>%
      select(new_gene, gene_name, pvalue, pop)
    
    genes <- genes %>%
      filter(pvalue <= 5e-8)
    join <- left_join(genes,anno, by = c("new_gene"="gene_id"))
    join <- join %>%
      rename(gene_name = gene_name.x) %>%
      select(-c(strand,gene_type,gene_name.y))
    path = "/home/maya/COLOC/eQTL_data/" %&% tissue %&% "/" %&% population %&% "/genes.csv"
    fwrite(join, path)
  }
  
}

################################################################################ Caudate_basal_ganglia ################################################################################IS DONE -- do for 12 other brain tisssues -- need to have genes.csv file made in each tisse directory 

tissueList = c("Caudate_basal_ganglia")
anno <- fread("/home/chris/gencode/gene_annotation_v38_strand.txt", header = T, sep=' ')

popList = c("META_Anx_colocs")

for (population in popList) {
  for (tissue in tissueList) {
    new_pop = gsub('_coloc','',population)
    data <- fread('/home/maya/prediXcan/prediXcanOUTPUT/ANX/META/META_Anx_Caudate_Basal_Ganglia.csv')
    new_data <- data %>%
      mutate(new_gene = sapply(strsplit(gene, "[.]"), `[`, 1)) %>%
      mutate(pop = new_pop)
    
    genes <- new_data %>%
      select(new_gene, gene_name, pvalue, pop)
    
    genes <- genes %>%
      filter(pvalue == 5.30681059897932e-21)
    join <- left_join(genes,anno, by = c("new_gene"="gene_id"))
    join <- join %>%
      rename(gene_name = gene_name.x) %>%
      select(-c(strand,gene_type,gene_name.y))
    path = "/home/maya/COLOC/eQTL_data/" %&% tissue %&% "/" %&% population %&% "/genes.csv"
    fwrite(join, path)
  }
  
}

################################################################################ Cerebellar_Hemisphere ################################################################################IS DONE -- do for 12 other brain tisssues -- need to have genes.csv file made in each tisse directory 
tissueList = c("Cerebellar_Hemisphere")
anno <- fread("/home/chris/gencode/gene_annotation_v38_strand.txt", header = T, sep=' ')

popList = c("META_Anx_colocs")

for (population in popList) {
  for (tissue in tissueList) {
    new_pop = gsub('_coloc','',population)
    data <- fread('/home/maya/prediXcan/prediXcanOUTPUT/ANX/META/META_Anx_Cerebellar_Hemisphere.csv')
    new_data <- data %>%
      mutate(new_gene = sapply(strsplit(gene, "[.]"), `[`, 1)) %>%
      mutate(pop = new_pop)
    
    genes <- new_data %>%
      select(new_gene, gene_name, pvalue, pop)
    
    genes <- genes %>%
      filter(pvalue <= 5e-8)
    join <- left_join(genes,anno, by = c("new_gene"="gene_id"))
    join <- join %>%
      rename(gene_name = gene_name.x) %>%
      select(-c(strand,gene_type,gene_name.y))
    path = "/home/maya/COLOC/eQTL_data/" %&% tissue %&% "/" %&% population %&% "/genes.csv"
    fwrite(join, path)
  }
  
}
################################################################################ Cerebellum ################################################################################IS DONE -- do for 12 other brain tisssues -- need to have genes.csv file made in each tisse directory 

tissueList = c("Cerebellum")
anno <- fread("/home/chris/gencode/gene_annotation_v38_strand.txt", header = T, sep=' ')

popList = c("META_Anx_colocs")

for (population in popList) {
  for (tissue in tissueList) {
    new_pop = gsub('_coloc','',population)
    data <- fread('/home/maya/prediXcan/prediXcanOUTPUT/ANX/META/META_Anx_Cerebellum.csv')
    new_data <- data %>%
      mutate(new_gene = sapply(strsplit(gene, "[.]"), `[`, 1)) %>%
      mutate(pop = new_pop)
    
    genes <- new_data %>%
      select(new_gene, gene_name, pvalue, pop)
    
    genes <- genes %>%
      filter(pvalue <= 5e-8)
    join <- left_join(genes,anno, by = c("new_gene"="gene_id"))
    join <- join %>%
      rename(gene_name = gene_name.x) %>%
      select(-c(strand,gene_type,gene_name.y))
    path = "/home/maya/COLOC/eQTL_data/" %&% tissue %&% "/" %&% population %&% "/genes.csv"
    fwrite(join, path)
  }
  
}

################################################################################ Cortex ################################################################################IS DONE -- do for 12 other brain tisssues -- need to have genes.csv file made in each tisse directory 
tissueList = c("Cortex")
anno <- fread("/home/chris/gencode/gene_annotation_v38_strand.txt", header = T, sep=' ')

popList = c("META_Anx_colocs")

for (population in popList) {
  for (tissue in tissueList) {
    new_pop = gsub('_coloc','',population)
    data <- fread('/home/maya/prediXcan/prediXcanOUTPUT/ANX/META/META_Anx_Cortex.csv')
    new_data <- data %>%
      mutate(new_gene = sapply(strsplit(gene, "[.]"), `[`, 1)) %>%
      mutate(pop = new_pop)
    
    genes <- new_data %>%
      select(new_gene, gene_name, pvalue, pop)
    
    genes <- genes %>%
      filter(pvalue <= 5e-8)
    join <- left_join(genes,anno, by = c("new_gene"="gene_id"))
    join <- join %>%
      rename(gene_name = gene_name.x) %>%
      select(-c(strand,gene_type,gene_name.y))
    path = "/home/maya/COLOC/eQTL_data/" %&% tissue %&% "/" %&% population %&% "/genes.csv"
    fwrite(join, path)
  }
  
}

################################################################################ Frontal_Cortex ################################################################################
tissueList = c("Frontal_Cortex")
anno <- fread("/home/chris/gencode/gene_annotation_v38_strand.txt", header = T, sep=' ')

popList = c("META_Anx_colocs")

for (population in popList) {
  for (tissue in tissueList) {
    new_pop = gsub('_coloc','',population)
    data <- fread('/home/maya/prediXcan/prediXcanOUTPUT/ANX/META/META_Anx_Frontal_Cortex_BA9.csv')
    new_data <- data %>%
      mutate(new_gene = sapply(strsplit(gene, "[.]"), `[`, 1)) %>%
      mutate(pop = new_pop)
    
    genes <- new_data %>%
      select(new_gene, gene_name, pvalue, pop)
    
    genes <- genes %>%
      filter(pvalue <= 5e-8)
    join <- left_join(genes,anno, by = c("new_gene"="gene_id"))
    join <- join %>%
      rename(gene_name = gene_name.x) %>%
      select(-c(strand,gene_type,gene_name.y))
    path = "/home/maya/COLOC/eQTL_data/" %&% tissue %&% "/" %&% population %&% "/genes.csv"
    fwrite(join, path)
  }
  
}
################################################################################ Hippocampus ################################################################################IS DONE -- do for 12 other brain tisssues -- need to have genes.csv file made in each tisse directory 
tissueList = c("Hippocampus")
anno <- fread("/home/chris/gencode/gene_annotation_v38_strand.txt", header = T, sep=' ')

popList = c("META_Anx_colocs")

for (population in popList) {
  for (tissue in tissueList) {
    new_pop = gsub('_coloc','',population)
    data <- fread('/home/maya/prediXcan/prediXcanOUTPUT/ANX/META/META_Anx_Hippocampus.csv')
    new_data <- data %>%
      mutate(new_gene = sapply(strsplit(gene, "[.]"), `[`, 1)) %>%
      mutate(pop = new_pop)
    
    genes <- new_data %>%
      select(new_gene, gene_name, pvalue, pop)
    
    genes <- genes %>%
      filter(pvalue ==5.3068105989793174e-21)
    join <- left_join(genes,anno, by = c("new_gene"="gene_id"))
    join <- join %>%
      rename(gene_name = gene_name.x) %>%
      select(-c(strand,gene_type,gene_name.y))
    path = "/home/maya/COLOC/eQTL_data/" %&% tissue %&% "/" %&% population %&% "/genes.csv"
    fwrite(join, path)
  }
  
}
################################################################################ Hypothalamus ################################################################################IS DONE -- do for 12 other brain tisssues -- need to have genes.csv file made in each tisse directory 
tissueList = c("Hypothalamus")
anno <- fread("/home/chris/gencode/gene_annotation_v38_strand.txt", header = T, sep=' ')

popList = c("META_Anx_colocs")

for (population in popList) {
  for (tissue in tissueList) {
    new_pop = gsub('_coloc','',population)
    data <- fread('/home/maya/prediXcan/prediXcanOUTPUT/ANX/META/META_Anx_Hypothalamus.csv')
    new_data <- data %>%
      mutate(new_gene = sapply(strsplit(gene, "[.]"), `[`, 1)) %>%
      mutate(pop = new_pop)
    
    genes <- new_data %>%
      select(new_gene, gene_name, pvalue, pop)
    
    genes <- genes %>%
      filter(pvalue ==9.990130348046941e-05)
    join <- left_join(genes,anno, by = c("new_gene"="gene_id"))
    join <- join %>%
      rename(gene_name = gene_name.x) %>%
      select(-c(strand,gene_type,gene_name.y))
    path = "/home/maya/COLOC/eQTL_data/" %&% tissue %&% "/" %&% population %&% "/genes.csv"
    fwrite(join, path)
  }
  
}
################################################################################ Nucleus_accumbens_basal_ganglia ################################################################################IS DONE -- do for 12 other brain tisssues -- need to have genes.csv file made in each tisse directory 
tissueList = c("Nucleus_accumbens_basal_ganglia")
anno <- fread("/home/chris/gencode/gene_annotation_v38_strand.txt", header = T, sep=' ')

popList = c("META_Anx_colocs")

for (population in popList) {
  for (tissue in tissueList) {
    new_pop = gsub('_coloc','',population)
    data <- fread('/home/maya/prediXcan/prediXcanOUTPUT/ANX/META/META_Anx_Nucleus_Accumbens_Basal_Ganglia.csv')
    new_data <- data %>%
      mutate(new_gene = sapply(strsplit(gene, "[.]"), `[`, 1)) %>%
      mutate(pop = new_pop)
    
    genes <- new_data %>%
      select(new_gene, gene_name, pvalue, pop)
    
    genes <- genes %>%
      filter(pvalue <= 5e-8)
    join <- left_join(genes,anno, by = c("new_gene"="gene_id"))
    join <- join %>%
      rename(gene_name = gene_name.x) %>%
      select(-c(strand,gene_type,gene_name.y))
    path = "/home/maya/COLOC/eQTL_data/" %&% tissue %&% "/" %&% population %&% "/genes.csv"
    fwrite(join, path)
  }
  
}
################################################################################ Putamen_Basal_Ganglia ################################################################################IS DONE -- do for 12 other brain tisssues -- need to have genes.csv file made in each tisse directory 
tissueList = c("Putamen_basal_ganglia")
anno <- fread("/home/chris/gencode/gene_annotation_v38_strand.txt", header = T, sep=' ')

popList = c("META_Anx_colocs")

for (population in popList) {
  for (tissue in tissueList) {
    new_pop = gsub('_coloc','',population)
    data <- fread('/home/maya/prediXcan/prediXcanOUTPUT/ANX/META/META_Anx_Putamen_Basal_Ganglia.csv')
    new_data <- data %>%
      mutate(new_gene = sapply(strsplit(gene, "[.]"), `[`, 1)) %>%
      mutate(pop = new_pop)
    
    genes <- new_data %>%
      select(new_gene, gene_name, pvalue, pop)
    
    genes <- genes %>%
      filter(pvalue <= 5e-8)
    join <- left_join(genes,anno, by = c("new_gene"="gene_id"))
    join <- join %>%
      rename(gene_name = gene_name.x) %>%
      select(-c(strand,gene_type,gene_name.y))
    path = "/home/maya/COLOC/eQTL_data/" %&% tissue %&% "/" %&% population %&% "/genes.csv"
    fwrite(join, path)
  }
  
}
################################################################################ Spinal_cord_cervical ################################################################################IS DONE -- do for 12 other brain tisssues -- need to have genes.csv file made in each tisse directory 
tissueList = c("Spinal_cord_cervical")
anno <- fread("/home/chris/gencode/gene_annotation_v38_strand.txt", header = T, sep=' ')

popList = c("META_Anx_colocs")

for (population in popList) {
  for (tissue in tissueList) {
    new_pop = gsub('_coloc','',population)
    data <- fread('/home/maya/prediXcan/prediXcanOUTPUT/ANX/META/META_Anx_Spinal_Cord_Cervical_c-1.csv')
    new_data <- data %>%
      mutate(new_gene = sapply(strsplit(gene, "[.]"), `[`, 1)) %>%
      mutate(pop = new_pop)
    
    genes <- new_data %>%
      select(new_gene, gene_name, pvalue, pop)
    
    genes <- genes %>%
      filter(pvalue == 9.990130348046987e-05)
    join <- left_join(genes,anno, by = c("new_gene"="gene_id"))
    join <- join %>%
      rename(gene_name = gene_name.x) %>%
      select(-c(strand,gene_type,gene_name.y))
    path = "/home/maya/COLOC/eQTL_data/" %&% tissue %&% "/" %&% population %&% "/genes.csv"
    fwrite(join, path)
  }
  
}
################################################################################ Substantia_nigra ################################################################################IS DONE -- do for 12 other brain tisssues -- need to have genes.csv file made in each tisse directory 
tissueList = c("Substantia_nigra")
anno <- fread("/home/chris/gencode/gene_annotation_v38_strand.txt", header = T, sep=' ')

popList = c("META_Anx_colocs")

for (population in popList) {
  for (tissue in tissueList) {
    new_pop = gsub('_coloc','',population)
    data <- fread('/home/maya/prediXcan/prediXcanOUTPUT/ANX/META/META_Anx_Substantia_nigra.csv')
    new_data <- data %>%
      mutate(new_gene = sapply(strsplit(gene, "[.]"), `[`, 1)) %>%
      mutate(pop = new_pop)
    
    genes <- new_data %>%
      select(new_gene, gene_name, pvalue, pop)
    
    genes <- genes %>%
      filter(pvalue <= 5e-8)
    join <- left_join(genes,anno, by = c("new_gene"="gene_id"))
    join <- join %>%
      rename(gene_name = gene_name.x) %>%
      select(-c(strand,gene_type,gene_name.y))
    path = "/home/maya/COLOC/eQTL_data/" %&% tissue %&% "/" %&% population %&% "/genes.csv"
    fwrite(join, path)
  }
  
}


################################################################################ Bonferroni Correction ################################################################################IS DONE -- do for 12 other brain tisssues -- need to have genes.csv file made in each tisse directory 

#With Bonferroni correction
for (population in popList) {
  for (tissue in tissueList) {
    new_pop = gsub('_coloc','',population)
    data <- fread('/home/grace/mets-gwas/Brain_Tissues_PrediXcan/PUKBB_' %&% new_pop %&% '_Brain_' %&% tissue %&% '.csv')
    new_data <- data %>%
      mutate(new_gene = sapply(strsplit(gene, "[.]"), `[`, 1)) %>%
      mutate(pop = new_pop)
    
    genes <- new_data %>%
      select(new_gene, gene_name, pvalue, pop)
    
    bonf <- nrow(data)
    bonf <- 5e-8 / bonf
    genes <- genes %>%
      filter(pvalue <= bonf)
    join <- left_join(genes,anno, by = c("new_gene"="gene_id"))
    join <- join %>%
      rename(gene_name = gene_name.x) %>%
      select(-c(strand,gene_type,gene_name.y))
    path = "/home/grace/coloc/eQTL_data/" %&% tissue %&% "/" %&% population %&% "/genes.csv"
    fwrite(join, path)
  }
  
}
