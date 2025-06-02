#!/bin/bash
##############################################################Anxious/Worrier Feelings Phenotype (1980) -- AFR ##########################################################################
################################################################ SPrediXcan  - Brain Tissues ###################################################################
## Brain Tissues: Amygdala, Anterior Cingulate Cortex, Caudate Basal Ganglia, Cerebellar Hemisphere, Cerebellum, Cortex, Frontal Cortex, Hippocampus, Hypothalamus, Nucleus Accumbens Basal Ganglia, Putamen Basal Ganglia, Spinal Cord Cervical C-1, Substantia Nigra) 
### Limbic System tissues -  Amygdala + Hippocampus

METAXCAN=/usr/local/bin/MetaXcan_software


## AFR - Amygdala ### 
echo '#1'
/usr/bin/python $METAXCAN/SPrediXcan.py \
--gwas_file /home/maya/prediXcan/harmonizedGWAS/Anx/hg38_PUKBB_AFR_Anx.tsv.gz \
--snp_column panel_variant_id \
--effect_allele_column effect_allele \
--non_effect_allele_column non_effect_allele \
--beta_column effect_size \
--se_column standard_error \
--model_db_path /home/maya/prediXcan/mashrMODELS/mashr_Brain_Amygdala.db \
--covariance /home/maya/prediXcan/mashrMODELS/mashr_Brain_Amygdala.txt.gz \
--keep_non_rsid \
--additional_output \
--model_db_snp_key varID \
--throw \
--output_file /home/maya/prediXcan/prediXcanOUTPUT/ANX/AFR/AFR_Anx_Amygdala.csv

####  AFR - Anterior Cingulate Cortex ### 
echo '#2'
/usr/bin/python $METAXCAN/SPrediXcan.py \
--gwas_file /home/maya/prediXcan/harmonizedGWAS/Anx/hg38_PUKBB_AFR_Anx.tsv.gz \
--snp_column panel_variant_id \
--effect_allele_column effect_allele \
--non_effect_allele_column non_effect_allele \
--beta_column effect_size \
--se_column standard_error \
--model_db_path /home/maya/prediXcan/mashrMODELS/mashr_Brain_Anterior_cingulate_cortex_BA24.db \
--covariance /home/maya/prediXcan/mashrMODELS/mashr_Brain_Anterior_cingulate_cortex_BA24.txt.gz \
--keep_non_rsid \
--additional_output \
--model_db_snp_key varID \
--throw \
--output_file /home/maya/prediXcan/prediXcanOUTPUT/ANX/AFR/AFR_Anx_Anterior_Cingulate_Cortex.csv

###  AFR - Caudate Basal Ganglia ### 
echo '#3'
/usr/bin/python $METAXCAN/SPrediXcan.py \
--gwas_file /home/maya/prediXcan/harmonizedGWAS/Anx/hg38_PUKBB_AFR_Anx.tsv.gz \
--snp_column panel_variant_id \
--effect_allele_column effect_allele \
--non_effect_allele_column non_effect_allele \
--beta_column effect_size \
--se_column standard_error \
--model_db_path /home/maya/prediXcan/mashrMODELS/mashr_Brain_Caudate_basal_ganglia.db  \
--covariance /home/maya/prediXcan/mashrMODELS/mashr_Brain_Caudate_basal_ganglia.txt.gz \
--keep_non_rsid \
--additional_output \
--model_db_snp_key varID \
--throw \
--output_file /home/maya/prediXcan/prediXcanOUTPUT/ANX/AFR/AFR_Anx_Caudate_Basal_Ganglia.csv

#### AFR - Cerebellar Hemisphere ### 
echo '#4'
/usr/bin/python $METAXCAN/SPrediXcan.py \
--gwas_file /home/maya/prediXcan/harmonizedGWAS/Anx/hg38_PUKBB_AFR_Anx.tsv.gz \
--snp_column panel_variant_id \
--effect_allele_column effect_allele \
--non_effect_allele_column non_effect_allele \
--beta_column effect_size \
--se_column standard_error \
--model_db_path /home/maya/prediXcan/mashrMODELS/mashr_Brain_Cerebellar_Hemisphere.db  \
--covariance /home/maya/prediXcan/mashrMODELS/mashr_Brain_Cerebellar_Hemisphere.txt.gz \
--keep_non_rsid \
--additional_output \
--model_db_snp_key varID \
--throw \
--output_file /home/maya/prediXcan/prediXcanOUTPUT/ANX/AFR/AFR_Anx_Cerebellar_Hemisphere.csv

####AFR- Cerebellum ### 
echo '#5'
/usr/bin/python $METAXCAN/SPrediXcan.py \
--gwas_file /home/maya/prediXcan/harmonizedGWAS/Anx/hg38_PUKBB_AFR_Anx.tsv.gz \
--snp_column panel_variant_id \
--effect_allele_column effect_allele \
--non_effect_allele_column non_effect_allele \
--beta_column effect_size \
--se_column standard_error \
--model_db_path /home/maya/prediXcan/mashrMODELS/mashr_Brain_Cerebellum.db  \
--covariance /home/maya/prediXcan/mashrMODELS/mashr_Brain_Cerebellum.txt.gz \
--keep_non_rsid \
--additional_output \
--model_db_snp_key varID \
--throw \
--output_file /home/maya/prediXcan/prediXcanOUTPUT/ANX/AFR/AFR_Anx_Cerebellum.csv

#### AFR - Cortex ### 
echo '#6'
/usr/bin/python $METAXCAN/SPrediXcan.py \
--gwas_file /home/maya/prediXcan/harmonizedGWAS/Anx/hg38_PUKBB_AFR_Anx.tsv.gz \
--snp_column panel_variant_id \
--effect_allele_column effect_allele \
--non_effect_allele_column non_effect_allele \
--beta_column effect_size \
--se_column standard_error \
--model_db_path /home/maya/prediXcan/mashrMODELS/mashr_Brain_Cortex.db  \
--covariance /home/maya/prediXcan/mashrMODELS/mashr_Brain_Cortex.txt.gz \
--keep_non_rsid \
--additional_output \
--model_db_snp_key varID \
--throw \
--output_file /home/maya/prediXcan/prediXcanOUTPUT/ANX/AFR/AFR_Anx_Cortex.csv

# ####AFR - Frontal Cortex ### 
echo '#7'
/usr/bin/python $METAXCAN/SPrediXcan.py \
--gwas_file /home/maya/prediXcan/harmonizedGWAS/Anx/hg38_PUKBB_AFR_Anx.tsv.gz \
--snp_column panel_variant_id \
--effect_allele_column effect_allele \
--non_effect_allele_column non_effect_allele \
--beta_column effect_size \
--se_column standard_error \
--model_db_path /home/maya/prediXcan/mashrMODELS/mashr_Brain_Frontal_Cortex_BA9.db  \
--covariance /home/maya/prediXcan/mashrMODELS/mashr_Brain_Frontal_Cortex_BA9.txt.gz \
--keep_non_rsid \
--additional_output \
--model_db_snp_key varID \
--throw \
--output_file /home/maya/prediXcan/prediXcanOUTPUT/ANX/AFR/AFR_Anx_Frontal_Cortex.csv


### AFR - Hippocampus ### 
echo '#8'
/usr/bin/python $METAXCAN/SPrediXcan.py \
--gwas_file /home/maya/prediXcan/harmonizedGWAS/Anx/hg38_PUKBB_AFR_Anx.tsv.gz \
--snp_column panel_variant_id \
--effect_allele_column effect_allele \
--non_effect_allele_column non_effect_allele \
--beta_column effect_size \
--se_column standard_error \
--model_db_path /home/maya/prediXcan/mashrMODELS/mashr_Brain_Hippocampus.db  \
--covariance /home/maya/prediXcan/mashrMODELS/mashr_Brain_Hippocampus.txt.gz \
--keep_non_rsid \
--additional_output \
--model_db_snp_key varID \
--throw \
--output_file /home/maya/prediXcan/prediXcanOUTPUT/ANX/AFR/AFR_Anx_Hippocampus.csv

# #### AFR - Hypothalamus ### 
echo '#9'
/usr/bin/python $METAXCAN/SPrediXcan.py \
--gwas_file /home/maya/prediXcan/harmonizedGWAS/Anx/hg38_PUKBB_AFR_Anx.tsv.gz \
--snp_column panel_variant_id \
--effect_allele_column effect_allele \
--non_effect_allele_column non_effect_allele \
--beta_column effect_size \
--se_column standard_error \
--model_db_path /home/maya/prediXcan/mashrMODELS/mashr_Brain_Hypothalamus.db \
--covariance /home/maya/prediXcan/mashrMODELS/mashr_Brain_Hypothalamus.txt.gz \
--keep_non_rsid \
--additional_output \
--model_db_snp_key varID \
--throw \
--output_file /home/maya/prediXcan/prediXcanOUTPUT/ANX/AFR/AFR_Anx_Hypothalamus.csv


# #### AFR - Nucleus Accumbens Basal Ganglia ### 
echo '#10'
/usr/bin/python $METAXCAN/SPrediXcan.py \
--gwas_file /home/maya/prediXcan/harmonizedGWAS/Anx/hg38_PUKBB_AFR_Anx.tsv.gz \
--snp_column panel_variant_id \
--effect_allele_column effect_allele \
--non_effect_allele_column non_effect_allele \
--beta_column effect_size \
--se_column standard_error \
--model_db_path /home/maya/prediXcan/mashrMODELS/mashr_Brain_Nucleus_accumbens_basal_ganglia.db \
--covariance /home/maya/prediXcan/mashrMODELS/mashr_Brain_Nucleus_accumbens_basal_ganglia.txt.gz \
--keep_non_rsid \
--additional_output \
--model_db_snp_key varID \
--throw \
--output_file /home/maya/prediXcan/prediXcanOUTPUT/ANX/AFR/AFR_Anx_Nucleus_Accumbens_Basal_Ganglia.csv

# #### AFR - Putamen Basal Ganglia ### 
echo '#11'
/usr/bin/python $METAXCAN/SPrediXcan.py \
--gwas_file /home/maya/prediXcan/harmonizedGWAS/Anx/hg38_PUKBB_AFR_Anx.tsv.gz \
--snp_column panel_variant_id \
--effect_allele_column effect_allele \
--non_effect_allele_column non_effect_allele \
--beta_column effect_size \
--se_column standard_error \
--model_db_path /home/maya/prediXcan/mashrMODELS/mashr_Brain_Putamen_basal_ganglia.db \
--covariance /home/maya/prediXcan/mashrMODELS/mashr_Brain_Putamen_basal_ganglia.txt.gz \
--keep_non_rsid \
--additional_output \
--model_db_snp_key varID \
--throw \
--output_file /home/maya/prediXcan/prediXcanOUTPUT/ANX/AFR/AFR_Anx_Putamen_Basal_Ganglia.csv

#### AFR  - Spinal Cord Cervical C-1 ### 
echo '#12'
/usr/bin/python $METAXCAN/SPrediXcan.py \
--gwas_file /home/maya/prediXcan/harmonizedGWAS/Anx/hg38_PUKBB_AFR_Anx.tsv.gz \
--snp_column panel_variant_id \
--effect_allele_column effect_allele \
--non_effect_allele_column non_effect_allele \
--beta_column effect_size \
--se_column standard_error \
--model_db_path /home/maya/prediXcan/mashrMODELS/mashr_Brain_Spinal_cord_cervical_c-1.db \
--covariance /home/maya/prediXcan/mashrMODELS/mashr_Brain_Spinal_cord_cervical_c-1.txt.gz \
--keep_non_rsid \
--additional_output \
--model_db_snp_key varID \
--throw \
--output_file /home/maya/prediXcan/prediXcanOUTPUT/ANX/AFR/AFR_Anx_Spinal_Cord_Cervical_c-1.csv

#### AFR - Substantia Nigra ### 
echo '#13'
/usr/bin/python $METAXCAN/SPrediXcan.py \
--gwas_file /home/maya/prediXcan/harmonizedGWAS/Anx/hg38_PUKBB_AFR_Anx.tsv.gz \
--snp_column panel_variant_id \
--effect_allele_column effect_allele \
--non_effect_allele_column non_effect_allele \
--beta_column effect_size \
--se_column standard_error \
--model_db_path /home/maya/prediXcan/mashrMODELS/mashr_Brain_Substantia_nigra.db \
--covariance /home/maya/prediXcan/mashrMODELS/mashr_Brain_Substantia_nigra.txt.gz \
--keep_non_rsid \
--additional_output \
--model_db_snp_key varID \
--throw \
--output_file /home/maya/prediXcan/prediXcanOUTPUT/ANX/AFR/AFR_Anx_Substantia_nigra.csv 
