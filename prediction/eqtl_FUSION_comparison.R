# Want to compare results from Shuangjia's eQTL analysis of centromere haplotypes to FUSION results
# She ran tests looking at the effect of centromere haplotypes on genes on the same chromosome (so no trans effects)
# My predicitive modeling runs included centromere haplotypes across all chromosomes--> will compare results for now but could re-run the eQTL analysis to match

# How to compare:
# eQTL gives p-values telling us if variants have a significant effect on expression of genes (one p-value for each variant-gene pair)
# FUSION gives the R2 of each model that considers all variants in order to predict the expression of a single gene (so not a one variant-gene pair)
# Can look to see what the best performing models were for the genes of interest (those w/ differing performance between models considering centromeres vs. those that were not)
# Depending on best performing models, can look at model weights to see what particular variants are driving model performance most
# Note that model type will determine range of weights (certain only focus on a few variants, others consider all)
# Can see if the weights of the models are similar to the significant eQTLs for the set of genes
# Can also just look at the top eQTLs for the genes of interest and see if any centromere haplotypes are among them

library(tidyverse)
library(ggplot2)
library(dplyr)
library(arrow)

# Reading in FUSION model performance (w/ centromeres) across all tested genes
cent_df = read.table("/gpfs/gibbs/pi/ycgh/amm422/project/centromere/prediction/weights/tables/cent_df.tsv", sep= "\t")
no_cent_df = read.table("/gpfs/gibbs/pi/ycgh/amm422/project/centromere/prediction/weights/tables/no_cent_df.tsv", sep= "\t")

# Reading in genes identified as having differing R2 w/ and w/out centromeres (R2 diff > 0.05, at least one test w/ R2 > 0.1)
better_w_df <- read.table("/gpfs/gibbs/pi/ycgh/amm422/project/centromere/prediction/weights/tables/better_w_cent.tsv", sep = "\t", header = TRUE)
better_wo_df <- read.table("/gpfs/gibbs/pi/ycgh/amm422/project/centromere/prediction/weights/tables/better_wo_cent.tsv", sep = "\t", header = TRUE)

# Directory for model weights (need to read in .RDat file for each gene's model to load the weights df in the R workspace)
cent_results_dir <- "/gpfs/gibbs/pi/ycgh/amm422/project/centromere/prediction/weights/cent"
cent_RDat_files <- list.files(cent_results_dir, pattern = "\\.wgt.RDat$", full.names = TRUE)

## eQTL analysis includes both GATK and 1KGSV variants--> indicated by suffix of the variant_id column--> only need GATK variants
# Directory for centromere haplotype eQTL analysis (just centromere haplotype:gene pairs)
# Files named: chr*.nominal.cis_qtl_pairs.chr*.parquet
cent_eqtls_dir <- "/gpfs/gibbs/pi/ycgh/lushjia/project/SV/AFGR/RNA/centromere/eqtl/nominal"

# Directory for nearby variant eQTL analysis (variants w/in +/- 1 Mb)
# Files named: chr*.5callset.nominal.cis_qtl_pairs.chr*.parquet
nearby_eqtls_dir <- "/gpfs/gibbs/pi/ycgh/lushjia/project/SV/AFGR/RNA/hprc_v2/edge/collapse_eqtl/separately/gatk_1kgsv"

# Reading in eQTL results for genes w/ better performance w/ centromeres
cent_results_dfs <- list()
nearby_results_dfs <- list()
for (i in 1:nrow(better_w_df)) {
    gene <- better_w_df[i, "gene"]
    chr <- better_w_df[i, "chr"]
    ## Centromere haplotypes
    # Opening connection to cent eqtl dataset
    cent_pointer <- open_dataset(paste0(cent_eqtls_dir, "/", chr, ".nominal.cis_qtl_pairs.", chr, ".parquet"))
    # Filtering to only read in pairs w/ the gene of interest
    subset <- cent_pointer %>%
        filter(phenotype_id == gene) %>%
        collect()
    # Adding to list of subset dfs
    cent_results_dfs[[gene]] <- subset
    ## Nearby variants
    nearby_pointer <- open_dataset(paste0(nearby_eqtls_dir, "/", chr, ".5callset.nominal.cis_qtl_pairs.", chr, ".parquet"))
    # Filtering to only read in pairs w/ the gene of interest
    subset <- nearby_pointer %>%
        filter(phenotype_id == gene & str_detect(variant_id, ":GATK")) %>%
        collect()
    # Adding to list of subset dfs
    nearby_results_dfs[[gene]] <- subset
}
# Concatenating all results
cent_eqtls_df <- bind_rows(cent_results_dfs)
nearby_eqtls_df <- bind_rows(nearby_results_dfs)
better_w_eqtls <- bind_rows(list(cent_eqtls_df, nearby_eqtls_df))
better_w_eqtls <- better_w_eqtls %>% arrange(phenotype_id)

# Reading in eQTL results for genes w/ better performance w/out centromeres
cent_results_dfs <- list()
nearby_results_dfs <- list()
for (i in 1:nrow(better_wo_df)) {
    gene <- better_wo_df[i, "gene"]
    chr <- better_wo_df[i, "chr"]
    ## Centromere haplotypes
    # Opening connection to cent eqtl dataset
    cent_pointer <- open_dataset(paste0(cent_eqtls_dir, "/", chr, ".nominal.cis_qtl_pairs.", chr, ".parquet"))
    # Filtering to only read in pairs w/ the gene of interest
    subset <- cent_pointer %>%
        filter(phenotype_id == gene) %>%
        collect()
    # Adding to list of subset dfs
    cent_results_dfs[[gene]] <- subset
    ## Nearby variants
    nearby_pointer <- open_dataset(paste0(nearby_eqtls_dir, "/", chr, ".5callset.nominal.cis_qtl_pairs.", chr, ".parquet"))
    # Filtering to only read in pairs w/ the gene of interest
    subset <- nearby_pointer %>%
        filter(phenotype_id == gene & str_detect(variant_id, ":GATK")) %>%
        collect()
    # Adding to list of subset dfs
    nearby_results_dfs[[gene]] <- subset
}
# Concatenating all results
cent_eqtls_df <- bind_rows(cent_results_dfs)
nearby_eqtls_df <- bind_rows(nearby_results_dfs)
better_wo_eqtls <- bind_rows(list(cent_eqtls_df, nearby_eqtls_df))
better_wo_eqtls <- better_wo_eqtls %>% arrange(phenotype_id)

### Looking to see if genes w/ better predictive modeling performance w/ centromeres have centromere haplotype eQTLs
# Looking at the lists of significant eQTLs for each gene (nominal p-value < 0.05)
# NOTE: these are uncorrected p-values! So should do correction going forward
for (id in unique(better_w_eqtls$phenotype_id)) {
    cur_df <- better_w_eqtls %>% 
        filter(phenotype_id == id & pval_nominal < 0.05 & str_detect(variant_id, "centro")) %>% 
        arrange(pval_nominal)
    View(cur_df)
}

# Loading weights for ENSG00000145247.11 (OCIAD2)
# load(paste0(cent_results_dir))