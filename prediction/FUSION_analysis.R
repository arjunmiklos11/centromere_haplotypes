### File to read in FUSION results and perform analysis
# Note: when running FUSION, I set the --hsa_p option to 1 to try and run the program for all genes
# --> so some of the genes have very low heritablitiy estimates--> worth filterting these out after looking at the results
# --> also, FUSION has a check to skip genes w/ heritability < 0, so some genes were skipped (won't have a .wgt.RDat file

# Goal is to compare the R2 values from runs w/out centromere haplotypes to those w/ centromere haplotypes
# --> if R2 higher w/ centromere haplotypes, then they are contributing to variance explained by the model

library(tidyverse)
library(ggplot2)
library(dplyr)
library(matrixStats)
library(biomaRt)
library(conflicted)
library(karyoploteR)
library(GenomicRanges)
library(regioneR)
library(GenomeInfoDb)

### Heritability analysis
# Starting w/ looking at distributions of the heritability estimates for all genes w/ and w/out centromere haplotypes
cent_results_dir <- "/gpfs/gibbs/pi/ycgh/amm422/project/centromere/prediction/weights/cent"
no_cent_results_dir <- "/gpfs/gibbs/pi/ycgh/amm422/project/centromere/prediction/weights/no_cent"

# Looping over .hsq files to get hertitability estimates (columns: gene, V(g)/Vp, V(g)/Vp se, V(g)/Vp p-value; single row per file)
cent_hsq_files <- list.files(cent_results_dir, pattern = "\\.hsq$", full.names = TRUE)
no_cent_hsq_files <- list.files(no_cent_results_dir, pattern = "\\.hsq$", full.names = TRUE)

# W/ centromeres
tables <- list()
for (file in cent_hsq_files) {
    table <- read.table(file, header = FALSE, col.names = c("gene", "Vg_Vp", "se", "pval"))
    table$gene <- basename(table$gene)
    tables <- append(tables, list(table))
}
cent_hsq_df <- bind_rows(tables)

# W/out centromeres
tables <- list()
for (file in no_cent_hsq_files) {
    table <- read.table(file, header = FALSE, col.names = c("gene", "Vg_Vp", "se", "pval"))
    table$gene <- basename(table$gene)
    tables <- append(tables, list(table))
}
no_cent_hsq_df <- bind_rows(tables)

# Making histograms of heritabililty w/ and w/out centromeres
# W/ centromeres
ggplot(cent_hsq_df, aes(x = Vg_Vp)) + 
    geom_histogram(binwidth = 0.01, fill = "skyblue", color = "black", alpha = 0.7) +
    labs(title = "Heritability Estimates w/ Centromere Haplotypes", x = "V(G)/Vp", y = "Number of Genes") + 
    theme_minimal() +
    theme(
        axis.title = element_text(size = 16),   # Size for both X and Y axis labels
        axis.text = element_text(size = 14),    # Size for axis tick marks/numbers
        plot.title = element_text(size = 20)    # Size for the main plot title
    ) +
    coord_cartesian(
        xlim = c(-.3, 1.4),
        ylim = c(0, 575)
    )

# W/out centromeres
ggplot(no_cent_hsq_df, aes(x = Vg_Vp)) + 
    geom_histogram(binwidth = 0.01, fill = "skyblue", color = "black", alpha = 0.7) +
    labs(title = "Heritability Estimates w/out Centromere Haplotypes", x = "V(G)/Vp", y = "Number of Genes") + 
    theme_minimal() +
    theme(
        axis.title = element_text(size = 16),   # Size for both X and Y axis labels
        axis.text = element_text(size = 14),    # Size for axis tick marks/numbers
        plot.title = element_text(size = 20)    # Size for the main plot title
    ) +
    coord_cartesian(
        xlim = c(-.3, 1.4),
        ylim = c(0, 575)
    )

# Printing means/medians
print(paste0("Mean heritability w/ centromeres: ", mean(cent_hsq_df$Vg_Vp), "; Median heritability w/ centromeres: ", median(cent_hsq_df$Vg_Vp)))
print(paste0("Mean heritability w/out centromeres: ", mean(no_cent_hsq_df$Vg_Vp), "; Median heritability w/out centromeres: ", median(no_cent_hsq_df$Vg_Vp)))

# Printing number of genes excluded from FUSION modeling due to V(G)/Vp < 0
print(paste0("Number of genes w/ V(G)/Vp < 0 w/ centromeres: ", nrow(cent_hsq_df %>% filter(Vg_Vp < 0))))
print(paste0("Number of genes w/ V(G)/Vp < 0 w/out centromeres: ", nrow(no_cent_hsq_df %>% filter(Vg_Vp < 0))))

# Running KS-test to compare distributions of heritability estimates
ks_results <- ks.test(cent_hsq_df$Vg_Vp, no_cent_hsq_df$Vg_Vp)
print(paste0("KS-test p-value: ", ks_results$p.value)) # p-value ~ 5.6 * 10^-113

# Making scatter plot comparing heritability estimates w/ and w/out centromeres
# Making df for plot--> joining by gene name
merged_hsq_df <- inner_join(cent_hsq_df %>% select("gene", "Vg_Vp"), no_cent_hsq_df %>% select("gene", "Vg_Vp"), by = "gene")
# Renaming columns
colnames(merged_hsq_df) <- c("gene", "Vg_Vp_cent", "Vg_Vp_no_cent")
# Scatter plot
ggplot(merged_hsq_df, aes(x = Vg_Vp_no_cent, y = Vg_Vp_cent)) +
    geom_point(color = "skyblue", alpha = 0.3) +
    geom_smooth(method = "lm", se = FALSE, color = "red") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
    labs(title = "Heritability Estimates w/ Centromeres vs w/out Centromeres", x = "V(G)/Vp w/out Centromeres", y = "V(G)/Vp w/ Centromeres") +
    theme_minimal() +
    theme(
        axis.title = element_text(size = 16),   # Size for both X and Y axis labels
        axis.text = element_text(size = 14),    # Size for axis tick marks/numbers
        plot.title = element_text(size = 20)    # Size for the main plot title
    )

### Predictive modeling analysis (looking at R2 values from cross-validation)
# Could add option to filter genes based on heritability (some genes already excluded b/c V(G)/Vp < 0)
# --> could exclude via minimum V(G)/Vp or keep genes w/ V(G)/Vp +/- SE > 0 or filter by p-value

## Getting dfs--> can skip now and just read in generated .tsv files
cent_RDat_files <- list.files(cent_results_dir, pattern = "\\.wgt.RDat$", full.names = TRUE)
no_cent_RDat_files <- list.files(no_cent_results_dir, pattern = "\\.wgt.RDat$", full.names = TRUE)

# W/ centromeres
# Build df w/ rows as genes and columns for R2 and p-values for each model + V(G)/Vp, SE, and p-value
genes <- list()
blup_r2s <- list()
blup_pvals <- list()
bslmm_r2s <- list()
bslmm_pvals <- list()
enet_r2s <- list()
enet_pvals <- list()
lasso_r2s <- list()
lasso_pvals <- list()
top1_r2s <- list()
top1_pvals <- list()
hsqs <- list()
hsq_ses <- list()
hsq_pvals <- list()

for (file in cent_RDat_files) {
    load(paste0(file))
    genes <- append(genes, str_sub(basename(file), start = 1, end = -10))
    blup_r2s <- append(blup_r2s, cv.performance["rsq", "blup"])
    blup_pvals <- append(blup_pvals, cv.performance["pval", "blup"])
    bslmm_r2s <- append(bslmm_r2s, cv.performance["rsq", "bslmm"])
    bslmm_pvals <- append(bslmm_pvals, cv.performance["pval", "bslmm"])
    enet_r2s <- append(enet_r2s, cv.performance["rsq", "enet"])
    enet_pvals <- append(enet_pvals, cv.performance["pval", "enet"])
    lasso_r2s <- append(lasso_r2s, cv.performance["rsq", "lasso"])
    lasso_pvals <- append(lasso_pvals, cv.performance["pval", "lasso"])
    top1_r2s <- append(top1_r2s, cv.performance["rsq", "top1"])
    top1_pvals <- append(top1_pvals, cv.performance["pval", "top1"])
    hsqs <- append(hsqs, hsq[1])
    hsq_ses <- append(hsq_ses, hsq[2])
    hsq_pvals <- append(hsq_pvals, hsq.pv)
}
cent_df <- data.frame(
    gene = unlist(genes),
    blup_r2 = unlist(blup_r2s),
    blup_pval = unlist(blup_pvals),
    bslmm_r2 = unlist(bslmm_r2s),
    bslmm_pval = unlist(bslmm_pvals),
    enet_r2 = unlist(enet_r2s),
    enet_pval = unlist(enet_pvals),
    lasso_r2 = unlist(lasso_r2s),
    lasso_pval = unlist(lasso_pvals),
    top1_r2 = unlist(top1_r2s),
    top1_pval = unlist(top1_pvals),
    hsq = unlist(hsqs),
    hsq_se = unlist(hsq_ses),
    hsq_pval = unlist(hsq_pvals)
)
cent_df <- cent_df %>% column_to_rownames(var = "gene")

# W/out centromeres
# Build df w/ rows as genes and columns for R2 and p-values for each model + V(G)/Vp, SE, and p-value
genes <- list()
blup_r2s <- list()
blup_pvals <- list()
bslmm_r2s <- list()
bslmm_pvals <- list()
enet_r2s <- list()
enet_pvals <- list()
lasso_r2s <- list()
lasso_pvals <- list()
top1_r2s <- list()
top1_pvals <- list()
hsqs <- list()
hsq_ses <- list()
hsq_pvals <- list()

for (file in no_cent_RDat_files) {
    load(paste0(file))
    genes <- append(genes, str_sub(basename(file), start = 1, end = -10))
    blup_r2s <- append(blup_r2s, cv.performance["rsq", "blup"])
    blup_pvals <- append(blup_pvals, cv.performance["pval", "blup"])
    bslmm_r2s <- append(bslmm_r2s, cv.performance["rsq", "bslmm"])
    bslmm_pvals <- append(bslmm_pvals, cv.performance["pval", "bslmm"])
    enet_r2s <- append(enet_r2s, cv.performance["rsq", "enet"])
    enet_pvals <- append(enet_pvals, cv.performance["pval", "enet"])
    lasso_r2s <- append(lasso_r2s, cv.performance["rsq", "lasso"])
    lasso_pvals <- append(lasso_pvals, cv.performance["pval", "lasso"])
    top1_r2s <- append(top1_r2s, cv.performance["rsq", "top1"])
    top1_pvals <- append(top1_pvals, cv.performance["pval", "top1"])
    hsqs <- append(hsqs, hsq[1])
    hsq_ses <- append(hsq_ses, hsq[2])
    hsq_pvals <- append(hsq_pvals, hsq.pv)
}
no_cent_df <- data.frame(
    gene = unlist(genes),
    blup_r2 = unlist(blup_r2s),
    blup_pval = unlist(blup_pvals),
    bslmm_r2 = unlist(bslmm_r2s),
    bslmm_pval = unlist(bslmm_pvals),
    enet_r2 = unlist(enet_r2s),
    enet_pval = unlist(enet_pvals),
    lasso_r2 = unlist(lasso_r2s),
    lasso_pval = unlist(lasso_pvals),
    top1_r2 = unlist(top1_r2s),
    top1_pval = unlist(top1_pvals),
    hsq = unlist(hsqs),
    hsq_se = unlist(hsq_ses),
    hsq_pval = unlist(hsq_pvals)
)
no_cent_df <- no_cent_df %>% column_to_rownames(var = "gene")

# Replacing NAs w/ -Inf (cases where models didn't converge; don't want to select these as best models)
cent_df_copy <- cent_df
cent_df_copy[is.na(cent_df_copy)] <- -Inf
no_cent_df_copy <- no_cent_df
no_cent_df_copy[is.na(no_cent_df_copy)] <- -Inf

# Getting best model by R2 for each gene
models <- c("blup", "bslmm", "enet", "lasso", "top1")

cent_df$best_model <- models[max.col(cent_df_copy[, paste0(models, "_r2")])]
cent_df$best_r2 <- apply(cent_df[, paste0(models, "_r2")], 1, max, na.rm = TRUE)
cent_df$best_pval <- sapply(1:nrow(cent_df), function(i) {
    model_name <- cent_df$best_model[i]
    col_name <- paste0(model_name, "_pval")
    if (col_name %in% colnames(cent_df)) {
        return(cent_df[i, col_name])
    } else {
        return(NA)
    }
})
no_cent_df$best_model <- models[max.col(no_cent_df_copy[, paste0(models, "_r2")])]
no_cent_df$best_r2 <- apply(no_cent_df[, paste0(models, "_r2")], 1, max, na.rm = TRUE)
no_cent_df$best_pval <- sapply(1:nrow(no_cent_df), function(i) {
    model_name <- no_cent_df$best_model[i]
    col_name <- paste0(model_name, "_pval")
    if (col_name %in% colnames(no_cent_df)) {
        return(no_cent_df[i, col_name])
    } else {
        return(NA)
    }
})

# Writing to file
write.table(no_cent_df, file = "/gpfs/gibbs/pi/ycgh/amm422/project/centromere/prediction/weights/tables/no_cent_df.tsv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
write.table(cent_df, file = "/gpfs/gibbs/pi/ycgh/amm422/project/centromere/prediction/weights/tables/cent_df.tsv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

## Reading in dfs and performing analysis
cent_df = read.table("/gpfs/gibbs/pi/ycgh/amm422/project/centromere/prediction/weights/tables/cent_df.tsv", sep= "\t")
no_cent_df = read.table("/gpfs/gibbs/pi/ycgh/amm422/project/centromere/prediction/weights/tables/no_cent_df.tsv", sep="\t")

## Plotting histograms of best R2 values w/ and w/out centromeres
# W/ centromeres
ggplot(cent_df, aes(x = best_r2)) + 
    geom_histogram(binwidth = 0.01, fill = "skyblue", color = "black", alpha = 0.7) +
    labs(title = "Best R2 w/ Centromere Haplotypes", x = "R2", y = "Number of Genes") + 
    theme_minimal() +
    theme(
        axis.title = element_text(size = 16),   # Size for both X and Y axis labels
        axis.text = element_text(size = 14),    # Size for axis tick marks/numbers
        plot.title = element_text(size = 20)    # Size for the main plot title
    ) +
    coord_cartesian(
        xlim = c(0, .85),
        ylim = c(0, 3000)
    )

# W/out centromeres
ggplot(no_cent_df, aes(x = best_r2)) + 
    geom_histogram(binwidth = 0.01, fill = "skyblue", color = "black", alpha = 0.7) +
    labs(title = "Best R2 w/out Centromere Haplotypes", x = "R2", y = "Number of Genes") + 
    theme_minimal() +
    theme(
        axis.title = element_text(size = 16),   # Size for both X and Y axis labels
        axis.text = element_text(size = 14),    # Size for axis tick marks/numbers
        plot.title = element_text(size = 20)    # Size for the main plot title
    ) +
    coord_cartesian(
        xlim = c(0, .85),
        ylim = c(0, 3000)
    )

## Printing means/medians of best R2 values
print(paste0("Mean best R2 w/ centromeres: ", mean(cent_df$best_r2), "; Median best R2 w/ centromeres: ", median(cent_df$best_r2)))
print(paste0("Mean best R2 w/out centromeres: ", mean(no_cent_df$best_r2), "; Median best R2 w/out centromeres: ", median(no_cent_df$best_r2)))
# KS-test to compare distributions of best R2 values
ks_results <- ks.test(cent_df$best_r2, no_cent_df$best_r2)
print(paste0("KS-test p-value: ", ks_results$p.value)) # R2 distribution significantly different w/ and w/out centromeres (p ~ 0.003)

## Making scatter plot comparing best R2 values w/ and w/out centromeres
# Making df for plot--> joining by gene name
# Converting rownames to column
cent_df <- cent_df %>% rownames_to_column(var = "gene")
no_cent_df <- no_cent_df %>% rownames_to_column(var = "gene")
merged_r2_df <- inner_join(cent_df %>% select("gene", "best_r2"), no_cent_df %>% select("gene", "best_r2"), by = "gene")
# Renaming columns
colnames(merged_r2_df) <- c("gene", "best_r2_cent", "best_r2_no_cent")
# Scatter plot
ggplot(merged_r2_df, aes(x = best_r2_no_cent, y = best_r2_cent)) +
    geom_point(color = "skyblue", alpha = 0.3) +
    geom_smooth(method = "lm", se = FALSE, color = "red") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
    labs(title = "Best R2 w/ Centromeres vs w/out Centromeres", x = "R2 w/out Centromeres", y = "R2 w/ Centromeres") +
    theme_minimal() +
    theme(
        axis.title = element_text(size = 16),   # Size for both X and Y axis labels
        axis.text = element_text(size = 14),    # Size for axis tick marks/numbers
        plot.title = element_text(size = 20)    # Size for the main plot title
    )
# Printing number of included genes (genes w/ hsq > 0 both w/ and w/out centromeres)
print(paste0("Number of genes w/ hsq > 0 both w/ and w/out centromeres: ", nrow(merged_r2_df)))

## Looking for genes w/ biggest differences in best R2 values between w/ and w/out centromeres
# Looking in both directions; requiring that best R2 >0.1 in at least one of the tests
# Adding column for difference in best R2--> will do R2 w/ centromeres - R2 w/out centromeres (positive means better w/ centromeres and vice versa)
merged_r2_df$r2_diff <- merged_r2_df$best_r2_cent - merged_r2_df$best_r2_no_cent
# Searching for genes w/ better performance w/ centromeres--> >X difference and best R2 > 0.1 in at least one test
DIFF_THRESHOLD <- 0.05 # Choosing 0.05 for now (8 genes w/ better performance w/ centromeres, 22 w/ better performance w/out centromeres)
better_w_cent <- merged_r2_df %>% filter(r2_diff > DIFF_THRESHOLD & (best_r2_cent > 0.1 | best_r2_no_cent > 0.1)) %>% arrange(-r2_diff) %>% pull(gene)
# Searching for genes w/ better performance w/out centromeres--> <-X difference and best R2 > 0 .1 in at least one test
better_w_no_cent <- merged_r2_df %>% filter(r2_diff < -DIFF_THRESHOLD & (best_r2_cent > 0.1 | best_r2_no_cent > 0.1)) %>% arrange(r2_diff) %>% pull(gene)

## Remaking the scatter plot coloring the genes w/ large R2 differences
# Adding column to mark selected genes (better_w, better_wo, and normal)
merged_r2_df$category <- sapply(merged_r2_df$gene, function(x) {
    if (x %in% better_w_cent) {
        return("better_w")
    } 
    else if (x %in% better_w_no_cent) {
       return("better_wo")
    }
    else {
        return("normal")
    }
})
# Remaking the scatter plot w/ colored dots
ggplot(merged_r2_df, aes(x = best_r2_no_cent, y = best_r2_cent)) +
    geom_point(aes(color = category, alpha = category)) +
    geom_smooth(method = "lm", se = FALSE, color = "red") +
    scale_color_manual(values = c("better_w" = "green", "normal" = "skyblue", "better_wo" = "red")) +
    scale_alpha_manual(values = c("better_w" = 1, "normal" = 0.3, "better_wo" = 1)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
    geom_abline(slope = 1, intercept = DIFF_THRESHOLD, linetype = "dashed", color = "grey") +
    geom_abline(slope = 1, intercept = -DIFF_THRESHOLD, linetype = "dashed", color = "grey") +
    labs(title = "Best R2 w/ Centromeres vs w/out Centromeres", x = "R2 w/out Centromeres", y = "R2 w/ Centromeres") +
    theme_minimal() +
    theme(
        axis.title = element_text(size = 16),   # Size for both X and Y axis labels
        axis.text = element_text(size = 14),    # Size for axis tick marks/numbers
        plot.title = element_text(size = 20),    # Size for the main plot title
        legend.position = "none"
    )

## Making tables showing info of the identified genes (w/ high R2 differences)
# Include: gene coordinates (start, stop), distance to centromere, and gene name/description
# Reading in gene position info (for phenotype .bed file) + centromere position info
# Gene info
phenotype_df <- read_tsv("/gpfs/gibbs/pi/ycgh/lushjia/project/SV/AFGR/RNA/sample_430/covariance_update/normalized.RNA.no_chrY.bed.gz")
phenotype_df <- phenotype_df %>% rename(chr = `#chr`)
phenotype_pos_df <- phenotype_df %>% select(c(TargetID, chr, start, end))
phenotype_df <- phenotype_df %>% select(!c(chr, start, end))
# Centromere info
centromere_pos_df <- read_tsv("/gpfs/gibbs/pi/ycgh/lushjia/project/SV/AFGR/RNA/centromere/b38/centro.chr.pos.txt", col_names = c("chr", "start", "end"))

# Getting gene position + centromere distance
# Better w/ centromeres
better_w_df <- merged_r2_df %>% filter(category == "better_w") %>% arrange(-r2_diff) %>% select(c(gene, best_r2_cent, best_r2_no_cent, r2_diff))
better_w_df <- left_join(better_w_df, phenotype_pos_df, by = c("gene" = "TargetID"))
# Better w/out centromeres
better_wo_df <- merged_r2_df %>% filter(category == "better_wo") %>% arrange(r2_diff) %>% select(c(gene, best_r2_cent, best_r2_no_cent, r2_diff))
better_wo_df <- left_join(better_wo_df, phenotype_pos_df, by = c("gene" = "TargetID"))
# Getting distance to centromere on chromosome of gene
better_w_df <- better_w_df %>%
    rowwise() %>%
    mutate(dist_to_cent = (function(x, y, z) {
        cent_start <- centromere_pos_df %>% filter(z == chr) %>% pull(start)
        cent_end <- centromere_pos_df %>% filter(z == chr) %>% pull(end)
        if (cent_start > y){
            return(cent_start - y)
        }
        else{
            return(x - cent_end)
        }
    })(start, end, chr)) %>%
    ungroup()
better_wo_df <- better_wo_df %>%
    rowwise() %>%
    mutate(dist_to_cent = (function(x, y, z) {
        cent_start <- centromere_pos_df %>% filter(z == chr) %>% pull(start)
        cent_end <- centromere_pos_df %>% filter(z == chr) %>% pull(end)
        if (cent_start > y){
            return(cent_start - y)
        }
        else{
            return(x - cent_end)
        }
    })(start, end, chr)) %>%
    ungroup()
# Setting up Ensembl database (to query gene IDs against to get gene names/descriptions)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
# Getting gene info (symbol and description) from Ensembl
# Better w/ centromeres
better_w_df$ensembl_gene_id <- sapply(better_w_df$gene, function(x) {
    return(str_split_i(x, "\\.", 1))
})
gene_map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol", "description"),
  filters    = "ensembl_gene_id",
  values     = better_w_df$ensembl_gene_id,
  mart       = ensembl
)
better_w_df <- left_join(better_w_df, gene_map, by = "ensembl_gene_id")
better_w_df <- better_w_df %>% relocate(ensembl_gene_id, .after = gene)
# Better w/out centromeres
better_wo_df$ensembl_gene_id <- sapply(better_wo_df$gene, function(x) {
    return(str_split_i(x, "\\.", 1))
})
gene_map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol", "description"),
  filters    = "ensembl_gene_id",
  values     = better_wo_df$ensembl_gene_id,
  mart       = ensembl
)
better_wo_df <- left_join(better_wo_df, gene_map, by = "ensembl_gene_id")
better_wo_df <- better_wo_df %>% relocate(ensembl_gene_id, .after = gene)

# Writing tables w/ genes better or worse w/ centromeres to files
write.table(better_w_df, file = "/gpfs/gibbs/pi/ycgh/amm422/project/centromere/prediction/weights/tables/better_w_cent.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(better_wo_df, file = "/gpfs/gibbs/pi/ycgh/amm422/project/centromere/prediction/weights/tables/better_wo_cent.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

## Making plots for genes w/ better/worse performance w/ centromere haplotypes included
# Want to make box/dot plot showing distance to centromeres
# Want to make chromosome map plot showing locations of centromeres + genes w/ worse/better performance
# Boxplot:
# Making df (genes, distance to centromeres, category)
boxplot_df <- bind_rows(better_w_df, better_wo_df)
boxplot_df <- boxplot_df %>%
    rowwise() %>%
    mutate(category = (function(x) {
        if (x %in% (better_w_df %>% pull(gene))) {
            return("better_w")
        }
        else if (x %in% (better_wo_df %>% pull(gene))) {
           return("better_wo")
        }
    })(gene)) %>%
    ungroup()
boxplot_df <- boxplot_df %>% select(gene, dist_to_cent, category)
# Making plot
ggplot(boxplot_df, aes(x = factor(category), y = dist_to_cent, fill = factor(category))) +
    geom_boxplot(outlier.shape = NA) +
    # Keeps the scattered points aligned perfectly with the colored sub-boxes
    geom_jitter(width = 0.2, alpha = 1, color = "black") +
    labs(title = "Distance of genes w/ differing R2 to centromeres", x = "Performance (R2) with centromere haplotypes included", y = "Distance to centromere (on chromosome of considered gene)") +
    scale_x_discrete(labels = c("Better", "Worse")) +
    theme_minimal() +
    theme(
        axis.title = element_text(size = 16),   # Size for both X and Y axis labels
        axis.text = element_text(size = 14),    # Size for axis tick marks/numbers
        plot.title = element_text(size = 20),   # Size for the main plot title
        legend.position = "none",
        axis.text.x = element_text(
            size = 16
        )
    )
# Printing means/medians of distances to centromeres
print(paste0("Distance to centromeres for genes w/ better R2 w/ centromeres: mean: ", mean(better_w_df %>% pull(dist_to_cent)), ", median: ", median(better_w_df %>% pull(dist_to_cent))))
print(paste0("Distance to centromeres for genes w/ worse R2 w/ centromeres: mean: ", mean(better_wo_df %>% pull(dist_to_cent)), ", median: ", median(better_wo_df %>% pull(dist_to_cent))))
# Doing Welch's T-test to compare means--> p = 0.5001 (no significant difference)
t.test(better_w_df %>% pull(dist_to_cent), better_wo_df %>% pull(dist_to_cent))

# Chromosome map plot (using karyotypeR)
# Plotting against hg38 genome (just showing autosomes since we didn't consider genes on the sex chromosomes)
# Getting chromosome sizes of hg38 genome
chrom_info <- getChromInfoFromUCSC("hg38")
chrom_sizes <- chrom_info %>% 
    filter(chrom %in% paste0("chr", 1:22)) %>%
    arrange(as.numeric(str_sub(chrom, 4, -1)))
# Build full cytoband dataframe with arms + custom centromeres
custom_bands <- centromere_pos_df %>%
    filter(chr %in% paste0("chr", 1:22)) %>%
    left_join(chrom_sizes, by = c("chr" = "chrom")) %>%
    rename(cen_start = start, cen_end = end) %>%
    arrange(as.numeric(str_sub(chr, 4, -1)))
custom_bands_long <- bind_rows(
    # Left arm (p arm)
    custom_bands %>%
        transmute(chr      = chr,
                  start    = 1,
                  end      = cen_start - 1,
                  name     = paste0(str_sub(chr, 4, -1), "p"),
                  gieStain = "gneg"),
    # Centromere
    custom_bands %>%
        transmute(chr      = chr,
                  start    = cen_start,
                  end      = cen_end,
                  name     = paste0("cen", str_sub(chr, 4, -1)),
                  gieStain = "acen"),
    # Right arm (q arm)
    custom_bands %>%
        transmute(chr      = chr,
                  start    = cen_end + 1,
                  end      = size,
                  name     = paste0(str_sub(chr, 4, -1), "q"),
                  gieStain = "gneg")
) %>%
    arrange(as.numeric(str_sub(chr, 4, -1)), start)
# Convert to GRanges
custom_bands_gr <- GRanges(
    seqnames = custom_bands_long$chr,
    ranges   = IRanges(start = custom_bands_long$start, end = custom_bands_long$end),
    name     = custom_bands_long$name,
    gieStain = custom_bands_long$gieStain
)
# Plot
kp <- plotKaryotype(genome      = "hg38",
                    cytobands   = custom_bands_gr,
                    chromosomes = paste0("chr", 1:22))
# Adding markers for genes w/ better R2 w/ centromeres
better_w_gr <- GRanges(
    seqnames = better_w_df$chr,
    ranges = IRanges(start = better_w_df$start, end = better_w_df$end),
    labels = better_w_df$hgnc_symbol
)
better_wo_gr <- GRanges(
    seqnames = better_wo_df$chr,
    ranges = IRanges(start = better_wo_df$start, end = better_wo_df$end),
    labels = better_wo_df$hgnc_symbol
)
kpPlotMarkers(kp, better_w_gr,
              labels        = better_w_gr$labels,
              text.orientation = "horizontal",
              r1            = 0.5,
              label.color   = "green",
              line.color    = "green")

kpPlotMarkers(kp, better_wo_gr,
              labels        = better_wo_gr$labels,
              text.orientation = "horizontal",
              r1            = 0.5,
              label.color   = "red",
              line.color    = "red")
