################# Worrier/Anxious Feelings - META, EUR, AFR, CSA #################

#load libraries
library(data.table)
library(dplyr)
library(qqman)
install.packages("ggplot2")
library(ggplot2)
install.packages("MetBrewer")
library("MetBrewer")
### MetBrewer -- uses color palettes from paintings at the Metropolitan Museum of Art in NYC 

"%&%" = function(a,b) paste(a,b,sep="")
# reads in the gene_anno file and removes chrM and chrY cause they are not numeric 
#gene_anno <- fread("/home/chris/gencode/gene_annotation_v38_strand.txt", header = T, sep=' ') 
gene_anno <- read.table("/home/isabelle/snp_conversions/gencode38_geneid-to-name.txt", sep = ",", header = TRUE)

gene_anno <- gene_anno[gene_anno$chr != "chrM", ]
gene_anno <- gene_anno[gene_anno$chr != "chrY", ]
#gene_anno$chr <- as.numeric(gene_anno$chr)
gene_anno$chr_pos <- ifelse(gene_anno$chr %in% as.character(1:22), as.numeric(gene_anno$chr), NA)

#creates column (chr_pos) with chromosome number:start position
gene_anno <- gene_anno %>% mutate(chr_pos = ""%&% gene_anno$chr %&%":"%&% 
                                    gene_anno$start)
##META 
Total <- data.frame() #create starting data frame
#list of brain tissues
Tissues <- list("Amygdala","Anterior_Cingulate_Cortex_BA24", "Caudate_Basal_Ganglia",
                "Cerebellar_Hemisphere", "Cerebellum", "Cortex", 
                "Frontal_Cortex_BA9", "Hippocampus", "Hypothalamus", 
                "Nucleus_Accumbens_Basal_Ganglia", "Putamen_Basal_Ganglia", 
                "Spinal_Cord_Cervical_c-1", "Substantia_nigra")

#loop creating new data frame 
plots <- list('META')
phenos <- list('Anx')
for (plot in plots) 
  Total <- data.frame()
for (pheno in phenos) {
  Total <- data.frame()
  for (x in Tissues) {
    path <- paste("/home/maya/panUKBB-prediXcan/prediXcanOUTPUT/ANX/META/", 
                  plot,"_",pheno,"_",x,".csv", sep="") #finds file with specified tissue
    #opens file with the PrediXcan data from the specified tissue
    Tissue <- data.table::fread(path, header= T, sep=",") 
    Tissue <- mutate(Tissue, gene=strtrim(gene, 15)) %>% #gets rid of the decimal point
      mutate(Tissue, tissue= x) #adds tissue column
    Total <- rbind(Total,Tissue) #adds this data frame to the Total data frame
    Total <- arrange(Total,pvalue) #sorts by pval
  }
  newPath <- paste(pheno,'_',plot,'_Total.csv', sep='')
  data.table::fwrite(Total, newPath, append= FALSE, quote = 'auto', sep=',')
  
  gwas_file <- Total %>% 
  select(-c(pred_perf_pval,pred_perf_r2,pred_perf_qval)) 
  #edits the ensemble ids so there are no decimal points
  gwas_file <- mutate(gwas_file, gene = strtrim(gene, 15))
  #joins gwas with gene anno based on the ensemble ids
  gwas_file <- left_join(gwas_file,gene_anno, by=c("gene"="gene_id")) 
  #gets rid of NA data
  gwas_file <- na.omit(gwas_file)
  
  #adds column of -log10log(pvalue)
  gwas_file$logP <- as.numeric(-log10(gwas_file$pvalue))
  #adds column with transcription start position of each gene
  gwas_file <- gwas_file %>%  group_by(chr) %>% summarise(chr_len=max(start)) %>% 
    mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>% select(-chr_len) %>% 
    left_join(gwas_file, ., by=c("chr"="chr")) %>% 
    arrange(chr, start) %>% mutate(BPcum=start+tot)
  #creates data frame with the middle position of each chromosome in the "center" column
  axisdf <- gwas_file %>% group_by(chr) %>% summarize(center=(max(BPcum) + min(BPcum))/2)
  
  #title <- paste('Pan-UKB Phenotype: Worrier/Anxious Feelings Phenotype - META (n=501,478)', sep=' ')
  pathTitle <- paste("/home/maya/panUKBBpublication/plots/",pheno,'_',plot,'_Total.png',sep='')
  library("MetBrewer")
  ggplot(gwas_file,aes(x = BPcum, y = logP)) + 
    theme_bw() +
    geom_point(alpha=0.75, size = 2, aes(color = tissue, shape= tissue), show.legend=TRUE) + 
    ggtitle(title) + 
    scale_x_continuous(label = axisdf$chr, breaks = axisdf$center) + scale_color_manual(values = met.brewer("Nizami", 13, type = "continuous")) + 
    scale_shape_manual(values=c(1,2,3,4,5,6,7,8,9,10,11,12,13)) +
    labs(x = "Chromosome",y = "-log10(p)") + 
    geom_hline(yintercept = -log10(3.85e-06), color = "grey40", linetype = "dashed") + 
    theme(axis.text.x = element_text(angle = 60, size = 12, vjust = 0.5), legend.text= element_text(size=6))
  
  ggsave(pathTitle, plot=last_plot(), width = 10,height = 7)
  
}
}

##AFR
Total <- data.frame() #create starting data frame
#list of brain tissues
Tissues <- list("Amygdala","Anterior_Cingulate_Cortex_BA24", "Caudate_Basal_Ganglia",
                "Cerebellar_Hemisphere", "Cerebellum", "Cortex", 
                "Frontal_Cortex_BA9", "Hippocampus", "Hypothalamus", 
                "Nucleus_Accumbens_Basal_Ganglia", "Putamen_Basal_Ganglia", 
                "Spinal_Cord_Cervical_c-1", "Substantia_nigra")

#loop creating new data frame 
plots <- list('AFR')
phenos <- list('Anx')
for (plot in plots) 
  Total <- data.frame()
for (pheno in phenos) {
  Total <- data.frame()
  for (x in Tissues) {
    path <- paste("/home/maya/panUKBB-prediXcan/prediXcanOUTPUT/ANX/AFR/", 
                  plot,"_",pheno,"_",x,".csv", sep="") #finds file with specified tissue
    #opens file with the PrediXcan data from the specified tissue
    Tissue <- fread(path, header= T, sep=",") 
    Tissue <- mutate(Tissue, gene=strtrim(gene, 15)) %>% #gets rid of the decimal point
      mutate(Tissue, tissue= x) #adds tissue column
    Total <- rbind(Total,Tissue) #adds this data frame to the Total data frame
    Total <- arrange(Total,pvalue) #sorts by pval
  }
  newPath <- paste(pheno,'_',plot,'_Total.csv', sep='')
  fwrite(Total, newPath, append= FALSE, quote = 'auto', sep=',')
  
  gwas_file <- Total %>% 
    select(-c(pred_perf_pval,pred_perf_r2,pred_perf_qval)) 
  #edits the ensemble ids so there are no decimal points
  gwas_file <- mutate(gwas_file, gene = strtrim(gene, 15))
  #joins gwas with gene anno based on the ensemble ids
  gwas_file <- left_join(gwas_file,gene_anno, by=c("gene"="gene_id")) 
  #gets rid of NA data
  gwas_file <- na.omit(gwas_file)
  
  #adds column of -log10log(pvalue)
  gwas_file$logP <- as.numeric(-log10(gwas_file$pvalue))
  #adds column with transcription start position of each gene
  gwas_file <- gwas_file %>%  group_by(chr) %>% summarise(chr_len=max(start)) %>% 
    mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>% select(-chr_len) %>% 
    left_join(gwas_file, ., by=c("chr"="chr")) %>% 
    arrange(chr, start) %>% mutate(BPcum=start+tot)
  #creates data frame with the middle position of each chromosome in the "center" column
  axisdf <- gwas_file %>% group_by(chr) %>% summarize(center=(max(BPcum) + min(BPcum))/2)
  
  title <- paste('Pan-UKBB Phenotype: Worrier/Anxious Feelings Phenotype - AFR (n=6,233)', sep=' ')
  pathTitle <- paste("/home/maya/publication/panUKBBfiles/plots/",pheno,'_',plot,'_Total.png',sep='')
  ggplot(gwas_file,aes(x = BPcum, y = logP)) + 
    theme_bw() +
    geom_point(alpha=0.75, size = 2, aes(color = tissue, shape= tissue), show.legend=TRUE) + 
    ggtitle(title) + 
    scale_x_continuous(label = axisdf$chr, breaks = axisdf$center) + scale_color_manual(values = met.brewer("Nizami", 13, type = "continuous")) + 
    scale_shape_manual(values=c(1,2,3,4,5,6,7,8,9,10,11,12,13)) +
    labs(x = "Chromosome",y = "-log10(p)") + 
    geom_hline(yintercept = -log10(4.048596e-06), color = "grey40", linetype = "dashed") + 
    theme(axis.text.x = element_text(angle = 60, size = 6, vjust = 0.5), legend.text= element_text(size=6))
  
  ggsave(pathTitle, plot=last_plot(), width = 10,height = 7)
  
}
}

  
}
}
