##############################################################Worrier/Anxious Feelings (1980)###################################################################
#pan-ukb sumstats in /home/maya/PANUKBB/Pan-UKBB_sumstats/hg37_varID_categorical-1980-both_sexes.tsv.gz
#!/bin/bash

#reformat pan-ukb for GWAS_TOOLS, you only need to do this once 
#awk statement creates a chr_pos_ref_alt_build SNP string and then pulls the rest of the columns
zcat /home/kayla/PrediXcan/continuous-1389-both_sexes.tsv.bgz |awk '{print "chr" $1 "_" $2 "_" $3 "_" $4 "_b37\t" $0}' > hg37_varID_continuous-1389-both_sexes.tsv 
gzip hg37_varID_continuous-1389-both_sexes.tsv #note: this command takes a few minutes


GWAS_TOOLS=/home/wheelerlab3/summary-gwas-imputation/src
METAXCAN=/usr/local/bin/MetaXcan_software
DATA=/home/wheelerlab3/mets-sleep/data

#After reformatting sumstats --> liftover and harmonize sumstats to hg38 (run once for each phenotype) - converting from hg37 --> hg38

#GWAS harmonization - AFR 
echo '#1'
/usr/bin/python $GWAS_TOOLS/gwas_parsing.py \
-gwas_file /home/maya/PANUKBB/Pan-UKBB_sumstats/hg37_varID_categorical-1980-both_sexes.tsv.gz \
-liftover $DATA/liftover/hg19ToHg38.over.chain.gz \
-snp_reference_metadata $DATA/reference_panel_1000G/variant_metadata.txt.gz METADATA \
-output_column_map chrchr_pos_ref_alt_b37 variant_id \
-output_column_map ref non_effect_allele \
-output_column_map alt effect_allele \
-output_column_map beta_AFR effect_size \
-output_column_map se_AFR standard_error \
-output_column_map chr chromosome \
--chromosome_format \
-output_column_map pos position \
-output_column_map af_AFR frequency \
--insert_value sample_size 6233  --insert_value n_cases NA \
-output_order variant_id panel_variant_id chromosome position effect_allele non_effect_allele frequency pvalue zscore effect_size standard_error sample_size n_cases \
-output /home/maya/prediXcan/harmonizedGWAS/Anx/hg38_PUKBB_AFR_Anx.tsv.gz


#GWAS harmonization - EUR
echo '#2'
/usr/bin/python $GWAS_TOOLS/gwas_parsing.py \
-gwas_file /home/maya/PANUKBB/Pan-UKBB_sumstats/hg37_varID_categorical-1980-both_sexes.tsv.gz \
-liftover $DATA/liftover/hg19ToHg38.over.chain.gz \
-snp_reference_metadata $DATA/reference_panel_1000G/variant_metadata.txt.gz METADATA \
-output_column_map chrchr_pos_ref_alt_b37 variant_id \
-output_column_map ref non_effect_allele \
-output_column_map alt effect_allele \
-output_column_map beta_EUR effect_size \
-output_column_map se_EUR standard_error \
-output_column_map chr chromosome \
--chromosome_format \
-output_column_map pos position \
-output_column_map af_EUR frequency \
--insert_value sample_size 409672  --insert_value n_cases NA \
-output_order variant_id panel_variant_id chromosome position effect_allele non_effect_allele frequency pvalue zscore effect_size standard_error sample_size n_cases \
-output /home/maya/prediXcan/harmonizedGWAS/Anx/hg38_PUKBB_EUR_Anx.tsv.gz

#GWAS harmonization - CSA
echo '#3'
/usr/bin/python $GWAS_TOOLS/gwas_parsing.py \
-gwas_file /home/maya/PANUKBB/Pan-UKBB_sumstats/hg37_varID_categorical-1980-both_sexes.tsv.gz \
-liftover $DATA/liftover/hg19ToHg38.over.chain.gz \
-snp_reference_metadata $DATA/reference_panel_1000G/variant_metadata.txt.gz METADATA \
-output_column_map chrchr_pos_ref_alt_b37 variant_id \
-output_column_map ref non_effect_allele \
-output_column_map alt effect_allele \
-output_column_map beta_CSA effect_size \
-output_column_map se_CSA standard_error \
-output_column_map chr chromosome \
--chromosome_format \
-output_column_map pos position \
-output_column_map af_CSA frequency \
--insert_value sample_size 8111  --insert_value n_cases NA \
-output_order variant_id panel_variant_id chromosome position effect_allele non_effect_allele frequency pvalue zscore effect_size standard_error sample_size n_cases \
-output /home/maya/prediXcan/harmonizedGWAS/Anx/hg38_PUKBB_CSA_Anx.tsv.gz

#GWAS harmonization - EAS
echo '#4'
/usr/bin/python $GWAS_TOOLS/gwas_parsing.py \
-gwas_file /home/maya/PANUKBB/Pan-UKBB_sumstats/hg37_varID_categorical-1980-both_sexes.tsv.gz \
-liftover $DATA/liftover/hg19ToHg38.over.chain.gz \
-snp_reference_metadata $DATA/reference_panel_1000G/variant_metadata.txt.gz METADATA \
-output_column_map chrchr_pos_ref_alt_b37 variant_id \
-output_column_map ref non_effect_allele \
-output_column_map alt effect_allele \
-output_column_map beta_EAS effect_size \
-output_column_map se_EAS standard_error \
-output_column_map chr chromosome \
--chromosome_format \
-output_column_map pos position \
-output_column_map af_EAS frequency \
--insert_value sample_size 2475  --insert_value n_cases NA \
-output_order variant_id panel_variant_id chromosome position effect_allele non_effect_allele frequency pvalue zscore effect_size standard_error sample_size n_cases \
-output /home/maya/prediXcan/harmonizedGWAS/Anx/hg38_PUKBB_EAS_Anx.tsv.gz

#GWAS harmonization - AMR
echo '#5'
/usr/bin/python $GWAS_TOOLS/gwas_parsing.py \
-gwas_file /home/maya/PANUKBB/Pan-UKBB_sumstats/hg37_varID_categorical-1980-both_sexes.tsv.gz \
-liftover $DATA/liftover/hg19ToHg38.over.chain.gz \
-snp_reference_metadata $DATA/reference_panel_1000G/variant_metadata.txt.gz METADATA \
-output_column_map chrchr_pos_ref_alt_b37 variant_id \
-output_column_map ref non_effect_allele \
-output_column_map alt effect_allele \
-output_column_map beta_AMR effect_size \
-output_column_map se_AMR standard_error \
-output_column_map chr chromosome \
--chromosome_format \
-output_column_map pos position \
-output_column_map af_AMR frequency \
--insert_value sample_size 929  --insert_value n_cases NA \
-output_order variant_id panel_variant_id chromosome position effect_allele non_effect_allele frequency pvalue zscore effect_size standard_error sample_size n_cases \
-output /home/maya/prediXcan/harmonizedGWAS/Anx/hg38_PUKBB_AMR_Anx.tsv.gz

#GWAS harmonization - MID
echo '#6'
/usr/bin/python $GWAS_TOOLS/gwas_parsing.py \
-gwas_file /home/maya/PANUKBB/Pan-UKBB_sumstats/hg37_varID_categorical-1980-both_sexes.tsv.gz \
-liftover $DATA/liftover/hg19ToHg38.over.chain.gz \
-snp_reference_metadata $DATA/reference_panel_1000G/variant_metadata.txt.gz METADATA \
-output_column_map chrchr_pos_ref_alt_b37 variant_id \
-output_column_map ref non_effect_allele \
-output_column_map alt effect_allele \
-output_column_map beta_MID effect_size \
-output_column_map se_MID standard_error \
-output_column_map chr chromosome \
--chromosome_format \
-output_column_map pos position \
-output_column_map af_MID frequency \
--insert_value sample_size 1474  --insert_value n_cases NA \
-output_order variant_id panel_variant_id chromosome position effect_allele non_effect_allele frequency pvalue zscore effect_size standard_error sample_size n_cases \
-output /home/maya/prediXcan/harmonizedGWAS/Anx/hg38_PUKBB_MID_Anx.tsv.gz


#GWAS harmonization - META
echo '#7'
/usr/bin/python $GWAS_TOOLS/gwas_parsing.py \
-gwas_file /home/maya/PANUKBB/Pan-UKBB_sumstats/hg37_varID_categorical-1980-both_sexes.tsv.gz \
-liftover $DATA/liftover/hg19ToHg38.over.chain.gz \
-snp_reference_metadata $DATA/reference_panel_1000G/variant_metadata.txt.gz METADATA \
-output_column_map chrchr_pos_ref_alt_b37 variant_id \
-output_column_map ref non_effect_allele \
-output_column_map alt effect_allele \
-output_column_map beta_meta effect_size \
-output_column_map se_meta standard_error \
-output_column_map chr chromosome \
--chromosome_format \
-output_column_map pos position \
-output_column_map af_meta frequency \
--insert_value sample_size 501485  --insert_value n_cases NA  \
-output_order variant_id panel_variant_id chromosome position effect_allele non_effect_allele frequency pvalue zscore effect_size standard_error sample_size n_cases \
-output /home/maya/prediXcan/harmonizedGWAS/Anx/hg38_PUKBB_META_Anx.tsv.gz


## Created separate ancestry specific directories for outputs! 
