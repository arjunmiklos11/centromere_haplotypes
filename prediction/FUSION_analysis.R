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
# Want to add option to filter genes based on heritability (some genes already excluded b/c V(G)/Vp < 0)
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
# Writing to file
write.table(cent_df, file = "/gpfs/gibbs/pi/ycgh/amm422/project/centromere/prediction/weights/tables/cent_df.tsv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

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
# Writing to file
write.table(no_cent_df, file = "/gpfs/gibbs/pi/ycgh/amm422/project/centromere/prediction/weights/tables/no_cent_df.tsv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

## Reading in dfs and performing analysis
cent_df = read.table("/gpfs/gibbs/pi/ycgh/amm422/project/centromere/prediction/weights/tables/cent_df.tsv", sep= "\t")
no_cent_df = read.table("/gpfs/gibbs/pi/ycgh/amm422/project/centromere/prediction/weights/tables/no_cent_df.tsv", sep="\t")

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
