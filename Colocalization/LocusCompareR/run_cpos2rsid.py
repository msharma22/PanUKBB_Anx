#!/bin/bash

python3 cpos2rsid.py --input bim_AFR_gwas_sumstats.tsv.gz \
--chr chr \
--pos pos \
--output /home/maya/panUKBBpublication/LocusCompareR/rsid_GWAS.txt \
--build 38
