#install required packages if not already installed
if (!require(data.table)) { install.packages("data.table") }
if (!require(dplyr)) { install.packages("dplyr") }
if (!require(coloc)) { install.packages("coloc") }
if (!require(hash)) { install.packages("hash") }
if (!require(optparse)) { install.packages("optparse") }
if (!require(R.utils)) { install.packages("R.utils") }
if (!require(ggplot2)) { install.packages("ggplot2") }
if (!require(httr)) { install.packages("httr") }
if (!require(jsonlite)) { install.packages("jsonlite") }
if (!require(RMySQL)) { install.packages("RMySQL") }

#install and load 'locuscomparer' from GitHub if not already installed
if (!require(locuscomparer)) {
  if (!require(devtools)) {
    install.packages("devtools")
  }
  devtools::install_github("boxiangliu/locuscomparer")
}
#load libraries with suppressed messages
suppressPackageStartupMessages(invisible(library(data.table)))
suppressPackageStartupMessages(invisible(library(dplyr)))
suppressPackageStartupMessages(invisible(library(coloc)))
suppressPackageStartupMessages(invisible(library(hash)))
suppressPackageStartupMessages(invisible(library(optparse)))
suppressPackageStartupMessages(invisible(library(R.utils)))
suppressPackageStartupMessages(invisible(library(ggplot2)))
suppressPackageStartupMessages(invisible(library(httr)))
suppressPackageStartupMessages(invisible(library(jsonlite)))
suppressPackageStartupMessages(invisible(library(RMySQL)))
suppressPackageStartupMessages(invisible(library(locuscomparer)))

#define a string concatenation operator
`%&%` <- function(a, b) paste0(a, b)

#define command-line options
option_list <- list(
  make_option("--process", type = "character", default = "pqtl", help = "If using eqtl data, add process flag, otherwise pqtl is the default."),
  make_option("--genes", type = "character", help = "List of genes of intrest, see READ.ME for specifications."),
  make_option("--seqIDdir", type = "character", help = "File path for ARIC pQTLs seqid.txt"),
  make_option("--GWASdir", type = "character", help = "File path for GWAS summary statistics"),
  make_option("--pQTLdir", type = "character", help = "Directory path for pQTL data"),
  make_option("--eQTLdir", type = "character", help = "Directory path for eQTL data"),
  make_option("--CHR_input", type = "character", help = "Column name for GWAS chromosome number"),
  make_option("--BP_input", type = "character", help = "Column name for GWAS bp postion"),
  make_option("--A1_input", type = "character", help = "Column name for GWAS affect allele (A1)"),
  make_option("--A2_input", type = "character", help = "Column name for GWAS alternate allele (A2)"),
  make_option("--BETA_input", type = "character", help = "Column name for GWAS BETA"),
  make_option("--SE_input", type = "character", help = "Column name for GWAS SE"),
  make_option("--ID_input", type = "character", help = "Column name for GWAS rsID"),
  make_option("--outputdir", type = "character", help = "Name of directory where you want output, this will be made for you")
  
)

#parse arguments
opt <- parse_args(OptionParser(option_list=option_list))

#gets working directory and the new dir the user wants to make for output
parent_dir <- file.path(getwd(), opt$outputdir) 
dir.create(parent_dir, showWarnings = FALSE, recursive = TRUE)  #create if it doesn't exist

cat('Reading in genes\n')
#read in the genes of intrest
genes = fread(opt$genes)


#loop through each gene and create a directory
for (gene in genes$genes) {
  dir_path <- file.path(parent_dir, gene)
  dir.create(dir_path) 
}


#if using eqtl data
if (opt$process == "eqtl"){
  #read in gwas and eqtl datasets
  gwas = data.table::fread(opt$GWASdir)
  eQTLs <- data.table::fread(opt$eQTLdir)
  cat('what about here\n')
  
  #format eqtl data, break up variant_id column and rename certain columns for better handling
  eQTLs <- eQTLs %>%
    mutate(`#CHROM` = sapply(strsplit(variant_id,"_"), `[`, 1), #chromosome number
           POS = sapply(strsplit(variant_id,"_"), `[`, 2), #position
           REF = sapply(strsplit(variant_id,"_"), `[`, 3), #reference allele, w/ indels
           A1 = sapply(strsplit(variant_id,"_"), `[`, 4)) %>% #alternate allele w/ indels
    rename(OBS_CT = OBS_CT, #number of allele observations
           A1_FREQ = A1_FREQ, #A1 allele frequency
           BETA = BETA, #regression coefficient for A1 allele
           SE = SE, #standard error
           P = P) #p-value
  
  #filter eqtl data
  eQTLs <- eQTLs %>%
    filter(nchar(A1)==1, #remove indels from alt alleles
           nchar(REF) == 1) %>% #remove indels from ref alleles
    mutate(CHR = gsub("chr", "", `X.CHROM`), #create chr column
           ALT = A1) #create ALT column
  eQTLs <- eQTLs %>%
    mutate(chr_pos = (CHR %&% ":" %&% POS)) #create chr_pos column, what will be used to match on
  cat('IS THIS WORKING\n')
  #remame GWAS input for better handling
  #gwas <- gwas %>%
    #rename_with(~ c("CHR", "POS", "A1", "A2", "BETA", "SE", "ID"), 
                #.cols = c(opt$CHR_input, opt$BP_input, opt$A1_input, opt$A2_input, opt$BETA_input, opt$SE_input, opt$ID_input))
    gwas <- gwas %>%
      mutate(chr_pos = paste0(CHR, ":", POS))
  
    #gwas <- gwas %>%
    #mutate(chr_pos = (CHR %&% ":" %&% POS)) #create chr_pos column, what will be used to match on
    cat('Wut about here\n')
    genes <- fread(opt$genes)
    if (!"genes" %in% colnames(genes)) {
      setnames(genes, 1, "genes")
    }
    
    for (target_gene in genes$genes) {
      gene_eQTL <- eQTLs %>% filter(grepl(paste0("^", target_gene), gene_id)) # filter eQTL data for specific gene
      cat("princess where are you\n")
      
      if (nrow(gene_eQTL) == 0) {
        cat("Skipping gene", target_gene, "as it is not found in eQTL data\n")
        next
      }
      
      # Continue with processing...
    }
    
  #for each gene in the target genes list, perform the following:
  for (target_gene in genes$genes) {
    gene_eQTL <- eQTLs %>% filter(grepl(paste0("^", target_gene), gene_id)) #filter eqtl data for specific gene
    cat('princess where are you&\n')
    if (nrow(gene_eQTL) == 0) {
      cat("Skipping gene", target_gene, "as it is not found in eQTL data\n") #if gene is not in data, skip it
      next
    }
    
    joined = inner_join(gene_eQTL,gwas, by=c("chr_pos"="chr_pos")) #join gwas and eqtl data based on chr_pos
    
    #check for complement bases
    #build hash table (like a Python dictionary)
    bases = hash()
    bases[["A"]] <- "T"
    bases[["C"]] <- "G"
    bases[["G"]] <- "C"
    bases[["T"]] <- "A"
    cat('princess where are you!!!!!!!!\n')
    #pull SNPs that  match (assumes C|G and A|T SNPs aren't flipped)
    match = joined[(joined$A1.x==joined$A1.y & joined$REF == joined$A2),]
    cat('help!!!!!!!!\n')
    #remove ambiguous strand SNPs (A/T or C/G)
    a = joined[!(joined$REF=="A" & joined$ALT=="T") & !(joined$REF=="T" & joined$ALT=="A") & !(joined$REF=="C" & joined$ALT=="G") &
                 !(joined$REF=="G" & joined$ALT=="C") ]
    cat('is this still working!\n')
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
    matchsnps = rbind(match, compmatch, flipped, compflipped) %>% arrange(chr_pos)
    
    matchsnps <- matchsnps %>%
      distinct(POS.y, .keep_all = TRUE) #keep only distinct
    matchsnps <- matchsnps %>%
      na.omit() #omit NAs
    cat('ugh!\n')
    #format gwas data for coloc analysis
    gwascoloc = list("beta" = matchsnps$BETA.y, "varbeta" = (matchsnps$SE.y)^2, "snp" = matchsnps$chr_pos, "position" = matchsnps$POS.y,
                     "type" = "cc")
    cat('si funciona?\n')
    #format eqtl data for coloc analysis
    eqtlcoloc = list("beta" = matchsnps$BETA.x, "varbeta" = (matchsnps$SE.x)^2, "snp" = matchsnps$chr_pos, "position" = matchsnps$POS.x,
                     "type" = "quant", "N" = 10000, "MAF" = matchsnps$A1_FREQ, "sdY"=1 )
    cat('e qua?\n')
    #save the gwas/eqtl coloc formatted data in a binary file for downstream analysis
    saveRDS(gwascoloc, file = file.path(parent_dir, target_gene, paste0(target_gene, "_gwascoloc")))
    saveRDS(eqtlcoloc, file = file.path(parent_dir, target_gene, paste0(target_gene, "_eqtlcoloc")))
    
    #save matchsnp data for LD downstream   
    saveRDS(matchsnps, file = file.path(parent_dir, target_gene, paste0(target_gene, "_matchsnps")))
    
    #run single variant coloc
    cat("SVA results for:", target_gene, "\n")
    # Check if SNPs match and are in same order
    if(!all(gwascoloc$snp == eqtlcoloc$snp)) {
      stop("Error: SNP IDs in GWAS and eQTL datasets do not match or are not in the same order.")
    }
    
    if(length(gwascoloc$snp) < 2) {
      stop("Error: Need at least 2 SNPs for coloc analysis.")
    }
    
     my.res = coloc.abf(dataset1=gwascoloc, dataset2=eqtlcoloc)
    
    #write out the results from the SVA coloc
    write.table(my.res$summary, file = file.path(parent_dir, target_gene, paste0(target_gene, "_sva_summary.tsv")), sep = "\t", row.names = FALSE, quote = FALSE)
    write.table(my.res$results, file = file.path(parent_dir, target_gene, paste0(target_gene, "_sva_results.tsv")), sep = "\t", row.names = FALSE, quote = FALSE)
    
    #make locus compare plots
    cat("Working on locuscompare results for:", target_gene, "\n")
    
    print(head(matchsnps))
    
    #format matchsnps data into file to be used as GWAS data input
    gwas_out <- matchsnps %>%
      dplyr::select(rsid.y, UK_pvalue) %>%
      rename(rsid = rsid.y , pval = UK_pvalue) %>%
      mutate(logp = -log10(pval))
    
    #format matchsnps data into file to be used as qtl data input
    eqtl_out <- matchsnps %>%
      select(rsid.x, P) %>%
      rename(rsid = rsid.x, pval = P) %>%
      mutate(logp = -log10(pval))
    
    #write to files
    gwas_locus <- file.path(parent_dir, target_gene, paste0(target_gene, "_gwas_locuscompare.txt"))
    eqtl_locus  <- file.path(parent_dir, target_gene, paste0(target_gene, "_eqtl_locuscompare.txt"))
    
    write.table(gwas_out, file = gwas_locus, sep = "\t", row.names = FALSE, quote = FALSE)
    write.table(eqtl_out,  file = eqtl_locus,  sep = "\t", row.names = FALSE, quote = FALSE)
    
    #plot and save it as a png
    plot_obj = locuscompare(in_fn1 = gwas_locus,
                            in_fn2 = eqtl_locus,
                            title1 = paste0("GWAS"),
                            title2 = paste0( "eQTL"),
                            population='AFR')  # Set population to AFR)
    
    ggsave(file.path(parent_dir, target_gene, paste0(target_gene, "_locuscompare.pdf")), plot_obj, width = 6, height = 5)
    
    
  }
  
}
