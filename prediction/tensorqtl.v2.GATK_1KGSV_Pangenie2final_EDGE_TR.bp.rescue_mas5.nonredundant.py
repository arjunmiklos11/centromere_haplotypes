#!/usr/bin/env python3
"""
FROM SHUANGJIA

Want to use to generate expanded variant call set (GATK + others) --> re-run predictive modeling experiment w/ and w/out centromeres included
for previously identified genes of interest

Unlike before, can't use FUSION to run predictive models b/c genotype data contains continuous dosages
Will instead use Shuangjia's script to run the 5 models (w/ susie swapped for bslmm) in R

## use V2 MC: GATK GT + 1KGSV + Pangenie(v2) + EDGE + TR genotype to run permutation
Usage: python tensorqtl.v2.GATK_1KGSV_Pangenie2final_EDGE_TR.bp.rescue_mas5.nonredundant.py chr1 gene

# rescue GATK, Pangenie, and edge, if in another callset MAS >= 10 and MAS>=5 in this callset 

# use final verison of pangenie v2 *
# GATK back based on vcfeval *

# apply genomeSTRiP monomophic filter to edge and node 
# update RNA-seq, covariates 
# run imputation missing to mean before apply maf filter to GATK, 1KGSV, Pangenie 

# add GATK back if GATK is represented by filtered edge - only through vcf coordinate 

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

sys.path.insert(1, os.path.dirname(tensorqtl.__file__))
from core import *
from post import *
import genotypeio, cis, trans, susie

# MODIFY TO TAKE BOTH CHROMOSOME AND GENE--> since we only want to run the models again for select genes--> usage: python tensorqtl.v2.GATK_1KGSV_Pangenie2final_EDGE_TR.bp.rescue_mas5.nonredundant.py CHR GENE
input_chr = sys.argv[1] 
input_gene = sys.argv[2]
# input_chr = "chr4"
# input_gene = "ENSG00000145247.11"

# input path of phenotype + covariants
phenotypes = '/gpfs/gibbs/pi/ycgh/lushjia/project/SV/AFGR/RNA/sample_430/covariance_update/normalized.RNA.no_chrY.bed.gz'
# prefix = 'method_5.ma_10'
mode = 'cis'
covariates = '/gpfs/gibbs/pi/ycgh/lushjia/project/SV/AFGR/RNA/sample_430/covariance_update/covariates.samp430.txt.gz'
# maf_threshold = 0 # filter by minor alelle count 10 later
seed = 1234

# read in phenotypes + covariants 
phenotype_df, phenotype_pos_df = read_phenotype_bed(phenotypes)
covariates_df = pd.read_csv(covariates, sep='\t', index_col=0).T
assert phenotype_df.columns.equals(covariates_df.index)
paired_covariate_df = None
interaction_df = None
group_s = None


# genotypee for each chr 
# read in genotype of GATK (with GT) using tensorqtl functions
genotype_path_GATK = f"/gpfs/gibbs/pi/ycgh/lushjia/project/SV/AFGR/RNA/sample_430/GATK/GATK_{input_chr}.samp430.add_resource"
genotype_df_GATK, variant_df_GATK = genotypeio.load_genotypes(genotype_path_GATK, select_samples=phenotype_df.columns,dosages = False)
# (1070401, 430)
# chr10 

# 1KGSV
genotype_path_1KGSV = f"/gpfs/gibbs/pi/ycgh/lushjia/project/SV/AFGR/RNA/sample_430/1KGSV/1KGSV.{input_chr}.samp430.add_resource"
genotype_df_1KGSV, variant_df_1KGSV = genotypeio.load_genotypes(genotype_path_1KGSV, select_samples=phenotype_df.columns,dosages = False)
# (2643, 430)

# Pangenie
genotype_path_Pangenie = f"/gpfs/gibbs/pi/ycgh/lushjia/project/SV/AFGR/RNA/hprc_v2/edge/pangenie_final/liftover/chr/{input_chr}.pangenie_all-samples_filtered_shapeit.samp430.lifted_over.concat.sorted.mod_id"
genotype_df_Pangenie, variant_df_Pangenie = genotypeio.load_genotypes(genotype_path_Pangenie, select_samples=phenotype_df.columns,dosages = False)
# (239725, 430)

# read in genotype of WW_edge using pandas csv
WW_edge_path = f"/gpfs/gibbs/pi/ycgh/lushjia/project/SV/AFGR/RNA/hprc_v2/edge/bubble_puncture/{input_chr}.independent.filter0.filter1000.filter_genomeSTRiP_rescue.AB.sorted.id.vcf.gz"
WW_genotype_edge_all = pd.read_csv(WW_edge_path,compression='gzip', sep="\t", comment='#',header = None)

# create genotype df for WW_edge
WW_genotype_edge_df = WW_genotype_edge_all.iloc[:, 9: ]
WW_genotype_edge_df.columns = phenotype_df.columns
WW_genotype_edge_df.index = WW_genotype_edge_all.iloc[:,2] # (1525771, 430) = (NUM VARIANTS ON CHR (DIFF. FOR EACH CHR), NUM SAMPLES (430))

# create variant df for WW_node 
WW_variant_edge_df = WW_genotype_edge_all.iloc[:,0:2]
WW_variant_edge_df.columns = ['chrom','pos']
WW_variant_edge_df.index =  WW_genotype_edge_all.iloc[:,2]
WW_variant_edge_df.index.name = 'id'
# (NUM VARIANTS, 2)

# no need to filter because already filter on vcf file 

# add rescued edge 
WW_edge_rescue = f"/gpfs/gibbs/pi/ycgh/lushjia/project/SV/AFGR/RNA/hprc_v2/edge/collapse/edge_back/rescue_edge/{input_chr}.rescue_gatk_pangenie_mas10.AB.sorted.id.vcf.gz"
WW_genotype_edge_all_rescue = pd.read_csv(WW_edge_rescue,compression='gzip', sep="\t", comment='#',header = None)

# create genotype df for WW_edge
WW_genotype_edge_df_rescue = WW_genotype_edge_all_rescue.iloc[:, 9: ]
WW_genotype_edge_df_rescue.columns = phenotype_df.columns
WW_genotype_edge_df_rescue.index = WW_genotype_edge_all_rescue.iloc[:,2] # (NUM VARIANTS, 430)

# create variant df for WW_node 
WW_variant_edge_df_rescue = WW_genotype_edge_all_rescue.iloc[:,0:2]
WW_variant_edge_df_rescue.columns = ['chrom','pos']
WW_variant_edge_df_rescue.index =  WW_genotype_edge_all_rescue.iloc[:,2]
WW_variant_edge_df_rescue.index.name = 'id'

# combine rescued edge with original edge
WW_genotype_edge_df_combined = pd.concat([WW_genotype_edge_df, WW_genotype_edge_df_rescue])
WW_variant_edge_df_combined = pd.concat([WW_variant_edge_df, WW_variant_edge_df_rescue])
print(f"EDGE variant number after filter: {WW_genotype_edge_df_combined.shape[0]}")


# read in genotype of Tandem repeat using pandas csv (Mark Chassion lab)
TR_path = f"/gpfs/gibbs/pi/ycgh/lushjia/project/SV/AFGR/RNA/hprc_v2/edge/tandem_repeat/{input_chr}.GTmat_bc.locusXsample.no_minus1.sorted.samp430.vcf.gz"
TR_genotype_all = pd.read_csv(TR_path,compression='gzip', sep="\t", comment='#',header = None)

# create genotype df for Tandem repeat
TR_genotype_df = TR_genotype_all.iloc[:, 9: ]
# missing "." replace by -9
TR_genotype_df = TR_genotype_df.replace('.', -9) # replace missing with -9
TR_genotype_df.columns = phenotype_df.columns
TR_genotype_df.index = TR_genotype_all.iloc[:,2] # (724, 430)
# make sure all elements are float 
TR_genotype_df = TR_genotype_df.astype(float)

# create variant df for Tandem repeat
TR_variant_df = TR_genotype_all.iloc[:,0:2]
TR_variant_df.columns = ['chrom','pos']
TR_variant_df.index =  TR_genotype_all.iloc[:,2]
TR_variant_df.index.name = 'id'

# filter rows with all 0 genotypes in TR_genotype_df
TR_variant_df = TR_variant_df.loc[(TR_genotype_df != 0).any(axis=1)] 
TR_genotype_df = TR_genotype_df.loc[(TR_genotype_df != 0).any(axis=1)] 
# filter not missing in <10 samples 
TR_variant_df = TR_variant_df.loc[(TR_genotype_df != -9).sum(axis=1) >= 10] 
TR_genotype_df = TR_genotype_df.loc[(TR_genotype_df != -9).sum(axis=1) >= 10] 
print(f"TR variant number after filter: {TR_genotype_df.shape[0]}")


# apply minor allele count filter to GAKT+1KGSV+Pangenie genotype and variant df 
# maf_threshold = 10 / (430 * 2) # allele count >=10 # total sample = 430
mas_threshold = 10

# NOTE: previously used median variant values for missing values
def impute_mean(genotypes_df, missing=-9):
    """Impute missing genotypes to mean"""
    m = genotypes_df == missing 
    ix = np.nonzero(m)[0]
    if len(ix) > 0:
        a = genotypes_df.sum(1) 
        b = m.sum(1)#.float()
        mu = (a - missing*b) / (genotypes_df.shape[1] - b)
        genotypes_df = genotypes_df.mask(m, mu, axis=0)
    return genotypes_df

# GATK 
genotype_df_GATK = impute_mean(genotype_df_GATK)

af_t = genotype_df_GATK.sum(1) / (2 * genotype_df_GATK.shape[1])
ix_t = af_t <= 0.5
a = (genotype_df_GATK > 0.5).sum(1).astype(int)
b = (genotype_df_GATK < 1.5).sum(1).astype(int)
ma_samples_t = np.where(ix_t, a, b)
mask_t = ma_samples_t >= mas_threshold

# output GATK variant number after MAC filter
print(f"GATK variant number after MAC filter: {sum(mask_t)}")

# readin GATK variant id in filtered edge - vcfeval 
gatk_edge_file = open(f"/gpfs/gibbs/pi/ycgh/lushjia/project/SV/AFGR/RNA/hprc_v2/edge/collapse/{input_chr}.gatk.rescue.mas5.txt",'r')
# 14:18224526:C:T:GATK
gatk_edge = set()
for line in gatk_edge_file:
    line = line.strip()
    gatk_edge.add(line)
# create boolean mask for GATK variant id in edge using index of genotype_df_GATK
mask_gatk_edge = variant_df_GATK.index.isin(gatk_edge)

# ouput final GATK variant number after both filter 
print(f"GATK variant number after MAC filter or edge filter: {sum(mask_t|mask_gatk_edge)}")

genotype_df_GATK = genotype_df_GATK[mask_t|mask_gatk_edge]
variant_df_GATK = variant_df_GATK[mask_t|mask_gatk_edge]


# 1KGSV 
genotype_df_1KGSV = impute_mean(genotype_df_1KGSV)

af_t = genotype_df_1KGSV.sum(1) / (2 * genotype_df_1KGSV.shape[1])
ix_t = af_t <= 0.5
a = (genotype_df_1KGSV > 0.5).sum(1).astype(int)
b = (genotype_df_1KGSV < 1.5).sum(1).astype(int)
ma_samples_t = np.where(ix_t, a, b)
mask_t = ma_samples_t >= mas_threshold

genotype_df_1KGSV = genotype_df_1KGSV[mask_t]
variant_df_1KGSV = variant_df_1KGSV[mask_t]
# output 1KGSV variant number after MAC filter
print(f"1KGSV variant number after MAC filter: {sum(mask_t)}")


# Pangenie
genotype_df_Pangenie = impute_mean(genotype_df_Pangenie)

af_t = genotype_df_Pangenie.sum(1) / (2 * genotype_df_Pangenie.shape[1])
ix_t = af_t <= 0.5
a = (genotype_df_Pangenie > 0.5).sum(1).astype(int)
b = (genotype_df_Pangenie < 1.5).sum(1).astype(int)
ma_samples_t = np.where(ix_t, a, b)
mask_t = ma_samples_t >= mas_threshold

# output Pangenie variant number after MAC filter
print(f"Pangenie variant number after MAC filter: {sum(mask_t)}")


# readin pangenie variant id in filtered edge - vcfeval 
pangenie_edge_file = open(f"/gpfs/gibbs/pi/ycgh/lushjia/project/SV/AFGR/RNA/hprc_v2/edge/collapse/{input_chr}.pangenie.rescue.mas5.txt",'r')
# 14:18224526:C:T:GATK
pangenie_edge = set()
for line in pangenie_edge_file:
    line = line.strip()
    pangenie_edge.add(line)
# create boolean mask for GATK variant id in edge using index of genotype_df_GATK
mask_pangenie_edge = variant_df_Pangenie.index.isin(pangenie_edge)

# ouput final Pangenie variant number after both filter 
print(f"Pangenie variant number after MAC filter or edge filter: {sum(mask_t|mask_pangenie_edge)}")

genotype_df_Pangenie = genotype_df_Pangenie[mask_t|mask_pangenie_edge]
variant_df_Pangenie = variant_df_Pangenie[mask_t|mask_pangenie_edge]


# pangenie if in redundant gatk - pangenie, remove 
pangenie_redundant_file = open(f"/gpfs/gibbs/pi/ycgh/lushjia/project/SV/AFGR/RNA/hprc_v2/edge/collapse/var_1_to_1_rescue/redundant/{input_chr}.gatk_pangenie_1_1.r2_0.5.txt",'r')
pangenie_redundant = set()
for line in pangenie_redundant_file:
    # var_chr1_3287:Pangenie  1:789640:T:C:GATK       1.0
    pangenie_var, gatk_var, r2 = line.strip().split()
    pangenie_redundant.add(pangenie_var)
# create boolean mask for Pangenie variant id in redundant gatk using index of genotype_df_Pangenie
mask_pangenie_redundant = variant_df_Pangenie.index.isin(pangenie_redundant)
genotype_df_Pangenie = genotype_df_Pangenie[~mask_pangenie_redundant]
variant_df_Pangenie = variant_df_Pangenie[~mask_pangenie_redundant]

# ouput final Pangenie variant number after both filter 
print(f"Pangenie variant number after MAC filter or edge filter or redundant with GATK: {variant_df_Pangenie.shape[0]}")


# edge if in redundant gatk - edge and pangenie - edge 
edge_redundant_file_gatk = open(f"/gpfs/gibbs/pi/ycgh/lushjia/project/SV/AFGR/RNA/hprc_v2/edge/collapse/var_1_to_1_rescue/redundant/{input_chr}.gatk_edge_1_1.r2_0.5.txt",'r')
edge_redundant = set()
for line in edge_redundant_file_gatk:
    # 10:103056510:G:A:GATK   >2852950>2852952:edge   0.9970856855031399
    gatk_var, edge_var, r2 = line.strip().split()
    edge_redundant.add(edge_var)
edge_redundant_file_pangenie = open(f"/gpfs/gibbs/pi/ycgh/lushjia/project/SV/AFGR/RNA/hprc_v2/edge/collapse/var_1_to_1_rescue/redundant/{input_chr}.pangenie_edge_1_1.r2_0.5.txt",'r')
for line in edge_redundant_file_pangenie:
    # var_chr20_545963:Pangenie       >98495986>98495987:edge 0.9969710825912574
    pangenie_var, edge_var, r2 = line.strip().split()
    edge_redundant.add(edge_var)
# create boolean mask for edge variant id in redundant gatk and pangenie using index of genotype_df_WW_edge
mask_edge_redundant = WW_variant_edge_df_combined.index.isin(edge_redundant)

WW_genotype_edge_df_combined = WW_genotype_edge_df_combined[~mask_edge_redundant]
WW_variant_edge_df_combined = WW_variant_edge_df_combined[~mask_edge_redundant]

print(f"EDGE variant number after removing redundant with GATK and Pangenie: {WW_genotype_edge_df_combined.shape[0]}")



# final genotyep df contain all variants from GATK+1KGSV+Pangennie+WW_node+WW_edge
genotype_df_4 = pd.concat([genotype_df_GATK, genotype_df_1KGSV, genotype_df_Pangenie, WW_genotype_edge_df_combined, TR_genotype_df]) 
# (2992306, 430)

variant_df_4 =  pd.concat([variant_df_GATK, variant_df_1KGSV, variant_df_Pangenie, WW_variant_edge_df_combined, TR_variant_df])
# (2992306, 2)


# sort genotype and variant dataframe 
# only input 1 chr so no need to sort chr 
variant_df_4_sort = variant_df_4.sort_values(by="pos",ascending=True) 
genotype_df_4_sort = genotype_df_4.reindex(variant_df_4_sort.index)

# Now want to output tables w/ and w/out centromere haplotypes added to run through predictive models
# phenotype of this gene --> selecting gene expression for input gene
phenotype_pos_df_gene = phenotype_pos_df.loc[input_gene].to_frame().transpose()
phenotype_df_gene = phenotype_df.loc[input_gene].to_frame().transpose()

# Reading in centromere haplotypes--> does it make sense to keep the haplotypes for the sex chromosomes given that we're excluding those genes from the analysis?
centromere_df = pd.read_csv("/gpfs/gibbs/pi/ycgh/lushjia/project/SV/AFGR/RNA/centromere/afgr.clade.sample.matrix.v2.tsv.gz", sep="\t")
centromere_df.set_index("clade", inplace = True)
centromere_pos_df = pd.read_csv("/gpfs/gibbs/pi/ycgh/lushjia/project/SV/AFGR/RNA/centromere/b38/centro.chr.pos.txt", sep="\t", header=None, names = ["chromosome", "start", "end"])
# Sorting columns (samples) to match inputs to InputGeneratorCis
centromere_df = centromere_df[list(genotype_df_4_sort.columns)]
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


# Generating input tables to run predictive models
logger = SimpleLogger()

igc = genotypeio.InputGeneratorCis(genotype_df_4_sort, variant_df_4_sort, phenotype_df_gene, phenotype_pos_df_gene, window=1000000) # Window defined as +/- 1 Mb around transcription start site
for k, (phenotype, genotypes, genotype_range, phenotype_id) in enumerate(igc.generate_data(verbose=True), 1):
    print(phenotype_id)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    genotypes_t = torch.tensor(genotypes, dtype=torch.float).to(device)
    variant_info = variant_df_4_sort.iloc[genotype_range[0]:genotype_range[-1]+1]
    # filter monomorphic variants
    mask_t = ~(genotypes_t == genotypes_t[:, [0]]).all(1)
    if mask_t.any():
        genotypes_t = genotypes_t[mask_t]
        mask = mask_t.cpu().numpy().astype(bool)
        variant_info = variant_info[mask]
        genotype_range = genotype_range[mask]
    # create a folder for each gene 
    os.makedirs(f"/gpfs/gibbs/pi/ycgh/amm422/project/centromere/prediction/full_variant_tests/inputs/{phenotype_id}", exist_ok=True)
    # phenotype - gene expression to file 
    np.savetxt(f"/gpfs/gibbs/pi/ycgh/amm422/project/centromere/prediction/full_variant_tests/inputs/{phenotype_id}/{phenotype_id}.txt", phenotype, fmt="%.6f")  
    ## genotype to file
    geno_matrix = genotypes_t.detach().cpu().numpy()
    # Saving version w/out centromere haplotypes
    np.savetxt(f"/gpfs/gibbs/pi/ycgh/amm422/project/centromere/prediction/full_variant_tests/inputs/{phenotype_id}/{phenotype_id}.full_variants.geno_no_centromeres.tsv", geno_matrix.T, delimiter="\t", fmt="%.6g")
    # Concatenating centromere haplotypes
    geno_matrix = np.concatenate((geno_matrix, centromere_df.to_numpy()), axis = 0)
    np.savetxt(f"/gpfs/gibbs/pi/ycgh/amm422/project/centromere/prediction/full_variant_tests/inputs/{phenotype_id}/{phenotype_id}.full_variants.geno.tsv", geno_matrix.T, delimiter="\t", fmt="%.6g")
    ## variant info to file
    # Saving version w/ out centromere haplotypes
    variant_info.to_csv(f"/gpfs/gibbs/pi/ycgh/amm422/project/centromere/prediction/full_variant_tests/inputs/{phenotype_id}/{phenotype_id}.full_variants.variant_info_no_centromeres.txt", sep='\t', index=True, header=False)
    # Adding centromere position info
    variant_info = pd.concat([variant_info, haplotype_pos_df], ignore_index=False)
    variant_info.to_csv(f"/gpfs/gibbs/pi/ycgh/amm422/project/centromere/prediction/full_variant_tests/inputs/{phenotype_id}/{phenotype_id}.full_variants.variant_info.txt", sep='\t', index=True, header=False)