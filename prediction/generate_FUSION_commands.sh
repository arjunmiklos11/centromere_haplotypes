#!/bin/bash

# Paths to executables
GEMMA=/home/amm422/.conda/envs/R/bin/gemma
PLINK=/vast/palmer/apps/avx2/software/PLINK/1.9b_6.21-x86_64/plink
GCTA=/gpfs/gibbs/pi/ycgh/amm422/software/fusion_twas-master/gcta_nr_robust

# Creating tmp and output directories
mkdir -p /gpfs/gibbs/pi/ycgh/amm422/project/centromere/prediction/tmp \
    /gpfs/gibbs/pi/ycgh/amm422/project/centromere/prediction/weights/no_cent \
    /gpfs/gibbs/pi/ycgh/amm422/project/centromere/prediction/weights/cent

JOBLIST=/gpfs/gibbs/pi/ycgh/amm422/project/centromere/prediction/FUSION_jobs.txt
> $JOBLIST

for GENE in $(ls /gpfs/gibbs/pi/ycgh/amm422/project/centromere/prediction/inputs_PLINK); do

    # --- No centromere ---
    echo "module load miniconda; conda activate /home/amm422/.conda/envs/R; module load PLINK/1.9b_6.21-x86_64; cd /gpfs/gibbs/pi/ycgh/amm422/software/fusion_twas-master; Rscript FUSION.compute_weights.R \
        --bfile /gpfs/gibbs/pi/ycgh/amm422/project/centromere/prediction/inputs_PLINK/$GENE/$GENE \
        --pheno /gpfs/gibbs/pi/ycgh/amm422/project/centromere/prediction/inputs_PLINK/$GENE/$GENE.tsv \
        --covar /gpfs/gibbs/pi/ycgh/amm422/project/centromere/prediction/covariates_PLINK.tsv \
        --tmp ../../project/centromere/prediction/tmp/${GENE}_no_cent \
        --out /gpfs/gibbs/pi/ycgh/amm422/project/centromere/prediction/weights/no_cent/${GENE} \
        --models blup,bslmm,enet,lasso,top1 \
        --hsq_p 1 \
        --save_hsq \
        --PATH_gemma $GEMMA \
        --PATH_plink $PLINK \
        --PATH_gcta $GCTA" >> $JOBLIST

    # --- With centromere ---
    echo "module load miniconda; conda activate /home/amm422/.conda/envs/R; module load PLINK/1.9b_6.21-x86_64; cd /gpfs/gibbs/pi/ycgh/amm422/software/fusion_twas-master; Rscript FUSION.compute_weights.R \
        --bfile /gpfs/gibbs/pi/ycgh/amm422/project/centromere/prediction/inputs_PLINK/$GENE/${GENE}_centromeres \
        --pheno /gpfs/gibbs/pi/ycgh/amm422/project/centromere/prediction/inputs_PLINK/$GENE/$GENE.tsv \
        --covar /gpfs/gibbs/pi/ycgh/amm422/project/centromere/prediction/covariates_PLINK.tsv \
        --tmp ../../project/centromere/prediction/tmp/${GENE}_cent \
        --out /gpfs/gibbs/pi/ycgh/amm422/project/centromere/prediction/weights/cent/${GENE} \
        --models blup,bslmm,enet,lasso,top1 \
        --hsq_p 1 \
        --save_hsq \
        --PATH_gemma $GEMMA \
        --PATH_plink $PLINK \
        --PATH_gcta $GCTA" >> $JOBLIST

done
