#!/usr/bin/env python3
"""
FROM SHUANGJIA

generate input files for gene expression prediction:
outputs 
    1. covariates table (sample * covariates)
    covariates.tsv
    2. per gene -> ./{phenotype_id}/
     - genotype: row: sample x column: variants (within 1M window of TSS of each gene)
     - variant info: row: variant_id chr pos (not used in prediction but helpful for knowing which variant in each column of the genotype matrix)
     - gene expression: one column, each row is one value of normalized gene expression for each sample

    # output example 
    covariates: /gpfs/gibbs/pi/ycgh/lushjia/project/SV/AFGR/RNA/hprc_v2/edge/susie/covariates.tsv
    genotype: /gpfs/gibbs/pi/ycgh/lushjia/project/SV/AFGR/RNA/hprc_v2/edge/collapse/var_1_to_1_rescue/separate_prediction/ENSG00000271913.5/ENSG00000271913.5.1kg_rescue.geno.tsv
    variant info: /gpfs/gibbs/pi/ycgh/lushjia/project/SV/AFGR/RNA/hprc_v2/edge/collapse/var_1_to_1_rescue/separate_prediction/ENSG00000271913.5/ENSG00000271913.5.1kg_rescue.variant_info.txt
    gene expression: /gpfs/gibbs/pi/ycgh/lushjia/project/SV/AFGR/RNA/hprc_v2/edge/susie/genes/ENSG00000227627.2/ENSG00000227627.2.txt

Usage: python input_generator.py chr1
"""

from __future__ import print_function
import tensorqtl
import pandas as pd
import numpy as np
from tensorqtl import genotypeio, cis, post
import qtl.plot

from datetime import datetime
import os
import re
import pickle
import argparse
from collections import defaultdict
import sys
import torch 
from xarray import DataArray
from pandas_plink import read_plink1_bin, write_plink1_bin

sys.path.insert(1, os.path.dirname(tensorqtl.__file__))
from core import *
from post import *
import genotypeio, cis, trans, susie

### Fix for writing files w/ pandas_plink
_old_to_csv = pd.DataFrame.to_csv
def _fixed_to_csv(*args, **kwargs):
    if 'line_terminator' in kwargs:
        kwargs['lineterminator'] = kwargs.pop('line_terminator')
    return _old_to_csv(*args, **kwargs)
pd.DataFrame.to_csv = _fixed_to_csv

input_chr = sys.argv[1] # don't run for sex chromosomes

# input path of phenotype + covariants
phenotypes = '/gpfs/gibbs/pi/ycgh/lushjia/project/SV/AFGR/RNA/sample_430/covariance_update/normalized.RNA.no_chrY.bed.gz' # 16572 genes w/ measured gene expression--> remove genes on sex chromosomes to get 16034 genes
mode = 'cis'
covariates = '/gpfs/gibbs/pi/ycgh/lushjia/project/SV/AFGR/RNA/sample_430/covariance_update/covariates.samp430.txt.gz'
seed = 1234

# read in phenotypes + covariants 
phenotype_df, phenotype_pos_df = read_phenotype_bed(phenotypes) # Phenotype df is gene expression
covariates_df = pd.read_csv(covariates, sep='\t', index_col=0).T # Covariate matrix has PCs and surrogate variables
assert phenotype_df.columns.equals(covariates_df.index)
paired_covariate_df = None
interaction_df = None
group_s = None

# genotype for each chr 
# read in genotype of GATK (with GT) using tensorqtl functions
genotype_path_GATK = f"/gpfs/gibbs/pi/ycgh/lushjia/project/SV/AFGR/RNA/sample_430/GATK/GATK_{input_chr}.samp430.add_resource"
genotype_df_GATK, variant_df_GATK = genotypeio.load_genotypes(genotype_path_GATK, dosages = False)
# (1070401, 430)

# apply minor allele count filter to GAKT
maf_threshold = 10 / (430 * 2) # allele count >=10; total sample = 430, total alleles = 430 * 2

# Using median genotype dosage (not mean) rounded to nearest whole number for missing genotype values
def impute_median(genotypes_df, missing=-9):
    """Impute missing genotypes to median"""
    m = genotypes_df == missing # df w/ True where value is missing (== -9), False otherwise
    ix = np.nonzero(m)[0] # Coordinates of missing values
    if len(ix) > 0:
        genotypes_df_copy = genotypes_df.copy()
        medians = round(genotypes_df_copy.replace(missing, np.nan).median(axis=1)) # Getting median value for each genotype w/ missing values and rounding to nearest whole number
        genotypes_df = genotypes_df.mask(m, medians, axis=0) # Replacing missing values w/ medians
        ### For mean instead:
        # a = genotypes_df.sum(1) # Sum of genotype dosages by variant
        # b = m.sum(1) #.float(), # Number of missing values by variant
        # mu = (a - missing*b) / (genotypes_df.shape[1] - b) # Mean dosage for samples w/ dosage values
        # genotypes_df = genotypes_df.mask(m, mu, axis=0) # Replace missing values w/ mean dosage
    return genotypes_df

# GATK 
genotype_df_GATK = impute_median(genotype_df_GATK)

af_t = genotype_df_GATK.sum(1) / (2 * genotype_df_GATK.shape[1]) # Calculating SNP allele frequency--> sum of allele dosage across all samples / number of total alleles (2x number of samples)
maf_t = np.where(af_t > 0.5, 1 - af_t, af_t) # Getting minor allele frequency --> returns 1 - SNP allele frequency if SNP allele frequency > 0.5
mask_t = maf_t >= maf_threshold # Filtering by minor allele frequency cutoff

genotype_df_GATK = genotype_df_GATK[mask_t]
variant_df_GATK = variant_df_GATK[mask_t]
# output GATK variant number after MAC filter
print(f"GATK variant number after MAC filter: {sum(mask_t)}")

# log
logger = SimpleLogger()

# phenotype of this chr --> selecting gene expression for only genes on the current chromosome
phenotype_pos_df_chr = phenotype_pos_df.loc[phenotype_pos_df['chr'] == input_chr]
phenotype_df_chr = phenotype_df.loc[phenotype_pos_df['chr'] == input_chr]

genotype_df_1kg = genotype_df_GATK
variant_df_1kg = variant_df_GATK
# sort --> by position on the chromosome
variant_df_1kg_sort = variant_df_1kg.sort_values(by="pos",ascending=True) 
genotype_df_1kg_sort = genotype_df_1kg.reindex(variant_df_1kg_sort.index)

# Reading in centromere haplotypes--> does it make sense to keep the haplotypes for the sex chromosomes given that we're excluding those genes from the analysis?
centromere_df = pd.read_csv("/gpfs/gibbs/pi/ycgh/lushjia/project/SV/AFGR/RNA/centromere/afgr.clade.sample.matrix.v2.tsv.gz", sep="\t")
centromere_df.set_index("clade", inplace = True)
centromere_pos_df = pd.read_csv("/gpfs/gibbs/pi/ycgh/lushjia/project/SV/AFGR/RNA/centromere/b38/centro.chr.pos.txt", sep="\t", header=None, names = ["chromosome", "start", "end"])
# Sorting columns (samples) to match inputs to InputGeneratorCis
centromere_df = centromere_df[list(genotype_df_1kg_sort.columns)]
# Making df that stores position of each centromere haplotype (matching format of variant_df_1kg_sort: columns: [snp, chrom, pos, end]--> end is extra column!)
haplotypes = list(centromere_df.index)
chromosomes = list(map(lambda x: x.split("_")[0], haplotypes))
haplotype_pos_df = pd.DataFrame(data={"snp": haplotypes,
                                      "chrom": chromosomes})
# Adding position info
haplotype_pos_df = pd.merge(haplotype_pos_df, centromere_pos_df, left_on="chrom", right_on="chromosome").drop(columns=["chromosome"])
# Renaming columns
haplotype_pos_df.rename(columns={'start': 'pos'}, inplace=True)
# Setting index to haplotype
haplotype_pos_df.set_index("snp", inplace=True)


igc = genotypeio.InputGeneratorCis(genotype_df_1kg_sort, variant_df_1kg_sort, phenotype_df_chr, phenotype_pos_df_chr, window=1000000) # Window defined as +/- 1 Mb around transcription start site
for k, (phenotype, genotypes, genotype_range, phenotype_id) in enumerate(igc.generate_data(verbose=True), 1):
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    genotypes_t = torch.tensor(genotypes, dtype=torch.float).to(device)
    variant_info = variant_df_1kg_sort.iloc[genotype_range[0]:genotype_range[-1]+1]
    # filter monomorphic variants
    mask_t = ~(genotypes_t == genotypes_t[:, [0]]).all(1)
    if mask_t.any():
        genotypes_t = genotypes_t[mask_t]
        mask = mask_t.cpu().numpy().astype(bool)
        variant_info = variant_info[mask]
        genotype_range = genotype_range[mask]
    # create a folder for each gene 
    os.makedirs(f"inputs_PLINK/{phenotype_id}", exist_ok=True)
    # os.makedirs(f"inputs/{phenotype_id}", exist_ok=True)
    # phenotype - gene expression to file 
    # np.savetxt(f"./inputs/{phenotype_id}/{phenotype_id}.txt", phenotype, fmt="%.6f")  
    ## genotype to file
    geno_matrix = genotypes_t.detach().cpu().numpy()
    # Saving version w/out centromere haplotypes
    geno_matrix_no_cent = geno_matrix
    # np.savetxt(f"./inputs/{phenotype_id}/{phenotype_id}.1kg_rescue.geno_no_centromeres.tsv", geno_matrix.T, delimiter="\t", fmt="%.6g")
    # Concatenating centromere haplotypes
    geno_matrix = np.concatenate((geno_matrix, centromere_df.to_numpy()), axis = 0)
    geno_matrix_cent = geno_matrix
    # np.savetxt(f"./inputs/{phenotype_id}/{phenotype_id}.1kg_rescue.geno.tsv", geno_matrix.T, delimiter="\t", fmt="%.6g")
    ## variant info to file
    # Saving version w/ out centromere haplotypes
    variant_info_no_cent = variant_info
    # variant_info.to_csv(f"./inputs/{phenotype_id}/{phenotype_id}.1kg_rescue.variant_info_no_centromeres.txt", sep='\t', index=True, header=False)
    # Adding centromere position info
    variant_info = pd.concat([variant_info, haplotype_pos_df], ignore_index=False)
    variant_info_cent = variant_info
    # variant_info.to_csv(f"./inputs/{phenotype_id}/{phenotype_id}.1kg_rescue.variant_info.txt", sep='\t', index=True, header=False)
    """
    To make plink file, need:  
        - sample_ids (in same order as columns in genotype matrix)
        - variant_ids (in same order as the subset that matches the +/- 1 Mb window around TSS) --> in format: 1:1022260:C:T:GATK (chrom:pos:a0:a1:GATK)
        - chromosomes
        - position
        - a0 (reference allele)
        - a1 (effect allele)
        - fid (family ID) --> for FAM file; same as sample id
        - iid (individual ID) --> for FAM file; same as sample id
        - father (array of 0s, length of number of samples)
        - mother (array of 0s, length of number of samples)
        - gender (array of 0s, length of number of samples)
        - trait (phenotype info--> not used by FUSION, so set to array of -9s, length of number of samples)
    """ 
    ### For no-centromere version
    sample_ids = list(genotype_df_1kg_sort.columns)
    num_samples = len(sample_ids)
    variant_ids_no_cent = list(variant_info_no_cent.index)
    chromosomes = list(map(lambda x: x.split(":")[0], variant_ids_no_cent))
    positions = list(map(lambda x: x.split(":")[1], variant_ids_no_cent))
    a0s = list(map(lambda x: x.split(":")[2], variant_ids_no_cent))
    a1s = list(map(lambda x: x.split(":")[3], variant_ids_no_cent))
    # Formatting xarray as input for pandas_plink to write to plink format
    G = DataArray(
        geno_matrix_no_cent.T,
        dims=["sample", "variant"],
        coords={
            "sample": sample_ids,
            "variant": variant_ids_no_cent,
            "snp": ("variant", variant_ids_no_cent),
            "chrom": ("variant", chromosomes),
            "pos": ("variant", positions),
            "a0": ("variant", a0s),
            "a1": ("variant", a1s),
            "fid": ("sample", sample_ids),
            "iid": ("sample", sample_ids),
            "father": ("sample", ["0"] * num_samples),
            "mother": ("sample", ["0"] * num_samples),
            "gender": ("sample", ["0"] * num_samples),
            "trait": ("sample", ["-9"] * num_samples)
        }
    )
    # Writing to plink file
    write_plink1_bin(G, f"./inputs_PLINK/{phenotype_id}/{phenotype_id}.bed")

    ### For centromere version
    # Just need to get additional info for chromosome haplotypes and append to existing lists
    num_haplotypes = len(haplotype_pos_df)
    variant_ids_cent = variant_ids_no_cent + list(haplotype_pos_df.index)
    chromosomes_cent = chromosomes + list(map(lambda x: x[3:], list(haplotype_pos_df["chrom"])))
    positions_cent = positions + list(haplotype_pos_df["pos"])
    a0s_cent = a0s + (["A"] * num_haplotypes) # Since I don't know the actual bps of the centromere haplotypes, just setting all centromere haplotype a0s to "A"
    a1s_cent = a1s + (["T"] * num_haplotypes) # Setting all centromere haplotype a1s to "T"
    # Formatting xarray as input for pandas_plink to write to plink format
    G_cent = DataArray(
        geno_matrix_cent.T,
        dims=["sample", "variant"],
        coords={
            "sample": sample_ids,
            "variant": variant_ids_cent,
            "snp": ("variant", variant_ids_cent),
            "chrom": ("variant", chromosomes_cent),
            "pos": ("variant", positions_cent),
            "a0": ("variant", a0s_cent),
            "a1": ("variant", a1s_cent),
            "fid": ("sample", sample_ids),
            "iid": ("sample", sample_ids),
            "father": ("sample", ["0"] * num_samples),
            "mother": ("sample", ["0"] * num_samples),
            "gender": ("sample", ["0"] * num_samples),
            "trait": ("sample", ["-9"] * num_samples)
        }
    )
    # Writing to plink file
    write_plink1_bin(G_cent, f"./inputs_PLINK/{phenotype_id}/{phenotype_id}_centromeres.bed")

    ### For PLINK phenotype file--> three columns: FID, IID, and phenotype--> rows are samples in same order as columns in genotype matrix; write w/ out headers
    plink_phenotype_df = pd.DataFrame(data = {"FID": sample_ids,
                                            "IID": sample_ids,
                                            "phenotype": phenotype})
    # Writing to file
    plink_phenotype_df.to_csv(f"./inputs_PLINK/{phenotype_id}/{phenotype_id}.tsv", header = False, sep="\t", index = False)

# output covariates table (sample * covariates) 
if input_chr == "chr1":
    covariates_df.to_csv("covariates.tsv", sep="\t", index=True, header=True, float_format="%.6f") # 431 * 41
    ### For PLINK covariate file--> need to add first two columns to existing covariate table to be FID and IID--> FID == IID == sample_ids
    new_cols = pd.DataFrame(data = {"FID": sample_ids,
                                    "IID": sample_ids})
    plink_covariates_df = pd.concat([new_cols, covariates_df.reset_index(drop=True)], axis="columns")
    # Writing to file
    plink_covariates_df.to_csv("covariates_PLINK.tsv", sep="\t", index=False)

