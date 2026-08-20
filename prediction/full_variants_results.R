# Re-ran predictive modeling using the more complete variant call set
# - First trial: just GATK variants
# - Second trial: GATK variants + SV + pangenie + EdgeDepth + danbing-tk
# Want to see if we still maintain the difference in model performance w/ and w/out centromeres w/ the additional variants
# NOTE: used different Bayesian predictive models between the two trials:
# First was w/ FUSION and used bslmm, second couldn't use FUSION b/c used continuous genotypes so had to substitute susie
# Both use the same other models (top 1, lasso, ridge, e-net)

# Want to compare the R2 differences w/ and w/out centromeres between the two trials
# Want to examine the model weights (check which was the best performing model for each trial) to see which variants in particular were driving the models (are there centromeric drivers)
# Compare weights to Shaungjia's eQTL analysis--> are the variants w/ high weights also identified as eQTLs?

library(tidyverse)

# Just looking at model performance for now (R2 comparisons, not looking at weights yet)
# Loading results from trial 1 (including lists of genes w/ better/worse performance w/ centromeres included)
better_w_t1_df <- read.table("/gpfs/gibbs/pi/ycgh/amm422/project/centromere/prediction/weights/tables/better_w_cent.tsv", sep="\t", header=TRUE)
better_wo_t1_df <- read.table("/gpfs/gibbs/pi/ycgh/amm422/project/centromere/prediction/weights/tables/better_wo_cent.tsv", sep="\t", header=TRUE)

cent_t1_df <- read.table("/gpfs/gibbs/pi/ycgh/amm422/project/centromere/prediction/weights/tables/cent_df.tsv", sep="\t", header=TRUE)
no_cent_t1_df <- read.table("/gpfs/gibbs/pi/ycgh/amm422/project/centromere/prediction/weights/tables/no_cent_df.tsv", sep="\t", header=TRUE)

# Processing results from trial 2
# Results stored in gene-specific folders within overall directory; stored in .tsv files w/ names GENE.cv.performance.w_cent.txt and GENE.cv.performance.wo_cent.txt
trial2_dir <- "/gpfs/gibbs/pi/ycgh/amm422/project/centromere/prediction/full_variant_tests/predictive_modeling"
subdirs <- list.files(trial2_dir, full.names = TRUE)

# W/ cents
dfs <- list()
for (subdir in subdirs){
    gene <- basename(subdir)
    performance_df <- as.data.frame(t(read.table(paste0(trial2_dir, "/", gene, "/", gene, ".cv.performance.w_cent.txt"))))
    performance_df$gene <- gene
    dfs[[gene]] <- performance_df
}
t2_w_cent_results <- bind_rows(dfs)
rownames(t2_w_cent_results) <- t2_w_cent_results$gene
t2_w_cent_results <- t2_w_cent_results %>% select(!gene)

# W/out cents
dfs <- list()
for (subdir in subdirs){
    gene <- basename(subdir)
    performance_df <- as.data.frame(t(read.table(paste0(trial2_dir, "/", gene, "/", gene, ".cv.performance.wo_cent.txt"))))
    performance_df$gene <- gene
    dfs[[gene]] <- performance_df
}
t2_wo_cent_results <- bind_rows(dfs)
rownames(t2_wo_cent_results) <- t2_wo_cent_results$gene
t2_wo_cent_results <- t2_wo_cent_results %>% select(!gene)

# Getting best performing models (get model name, R2, and p-value)
# Getting best model by R2 for each gene
models <- c("lasso", "ridge", "enet", "top1", "susie") # NOTE: ridge = blup, susie used intead of bslmm

t2_w_cent_results$best_model <- models[max.col(t2_w_cent_results[, paste0(models, "_rsq")])]
t2_w_cent_results$best_rsq <- apply(t2_w_cent_results[, paste0(models, "_rsq")], 1, max, na.rm = TRUE)
t2_w_cent_results$best_pval <- sapply(1:nrow(t2_w_cent_results), function(i) {
    model_name <- t2_w_cent_results$best_model[i]
    col_name <- paste0(model_name, "_pval")
    if (col_name %in% colnames(t2_w_cent_results)) {
        return(t2_w_cent_results[i, col_name])
    } else {
        return(NA)
    }
})
t2_wo_cent_results$best_model <- models[max.col(t2_wo_cent_results[, paste0(models, "_rsq")])]
t2_wo_cent_results$best_rsq <- apply(t2_wo_cent_results[, paste0(models, "_rsq")], 1, max, na.rm = TRUE)
t2_wo_cent_results$best_pval <- sapply(1:nrow(t2_wo_cent_results), function(i) {
    model_name <- t2_wo_cent_results$best_model[i]
    col_name <- paste0(model_name, "_pval")
    if (col_name %in% colnames(t2_wo_cent_results)) {
        return(t2_wo_cent_results[i, col_name])
    } else {
        return(NA)
    }
})

# Want to make line + dot plot that shows change in R2 from GATK variant set to full variant set
# Make separate plots for those w/ better/worse performance w/ centromeres using the GATK set
# Color lines/dots depending on if the R2 improved or got worse w/ the full variant set
# Make plots using R2 values w/ centromeres included
# Better w/ centromeres:
# Making plot df: need gene id, best R2 using GATK set, best R2 using full variant set--> use to get R2 diff + category (improved/worsened R2 w/ full set)
# Also getting best models for t1 and t2
genes <- list()
best_r2_gatk <- list()
best_r2_full <- list()
best_model_gatk <- list()
best_model_full <- list()
for (gene in better_w_t1_df %>% pull(gene)){
    genes <- append(genes, gene)
    best_r2_gatk <- append(best_r2_gatk, cent_t1_df[gene, "best_r2"])
    best_model_gatk <- append(best_model_gatk, cent_t1_df[gene, "best_model"])
    best_r2_full <- append(best_r2_full, t2_w_cent_results[gene, "best_rsq"])
    best_model_full <- append(best_model_full, t2_w_cent_results[gene, "best_model"])
}
plot_df <- data.frame(
    gene = unlist(genes),
    gatk = unlist(best_r2_gatk),
    best_model_gatk = unlist(best_model_gatk),
    full = unlist(best_r2_full),
    best_model_full = unlist(best_model_full)
)
plot_df$r2_diff <- plot_df$full - plot_df$gatk
plot_df <- plot_df %>%
    rowwise() %>%
    mutate(category = (function(x) {
        if (x > 0){
            return("increase")
        }
        else if (x == 0) {
           return("no change")
        }
        else {
           return("decrease")
        }
    })(r2_diff)) %>%
    ungroup()
# Converting to long format
plot_df <- plot_df %>%
    pivot_longer(
        cols = c("gatk", "full"),
        names_to = "variant_set",
        values_to = "best_r2"
    )
# Making plot
order <- c("gatk", "full")
plot_df$variant_set <- factor(plot_df$variant_set, levels = order)
colors <- c(
    "increase" = "green",
    "no change" = "grey",
    "decrease" = "red"
)
ggplot(plot_df, aes(x = variant_set, y = best_r2, color = category, group = gene)) + 
    geom_line(linewidth = 1) +
    geom_point(size = 4) +
    scale_color_manual(values = colors, name = "Change") +
    labs(title = "Change in R2 from GATK to Full Variant Set", x = "Variant Set (w/ centromere haplotypes)", y = "Best R2") +
    theme_minimal() +
    theme(
        plot.title = element_text(size = 20), # Main title
        axis.title = element_text(size = 18),               # Both X and Y axis labels
        axis.text = element_text(size = 16)                 # Both X and Y tick values
    )

# Better w/out centromeres:
# Making plot df:
genes <- list()
best_r2_gatk <- list()
best_r2_full <- list()
best_model_gatk <- list()
best_model_full <- list()
for (gene in better_wo_t1_df %>% pull(gene)){
    genes <- append(genes, gene)
    best_r2_gatk <- append(best_r2_gatk, cent_t1_df[gene, "best_r2"])
    best_model_gatk <- append(best_model_gatk, cent_t1_df[gene, "best_model"])
    best_r2_full <- append(best_r2_full, t2_w_cent_results[gene, "best_rsq"])
    best_model_full <- append(best_model_full, t2_w_cent_results[gene, "best_model"])
}
plot_df <- data.frame(
    gene = unlist(genes),
    gatk = unlist(best_r2_gatk),
    best_model_gatk = unlist(best_model_gatk),
    full = unlist(best_r2_full),
    best_model_full = unlist(best_model_full)
)
plot_df$r2_diff <- plot_df$full - plot_df$gatk
plot_df <- plot_df %>%
    rowwise() %>%
    mutate(category = (function(x) {
        if (x > 0){
            return("increase")
        }
        else if (x == 0) {
           return("no change")
        }
        else {
           return("decrease")
        }
    })(r2_diff)) %>%
    ungroup()
# Converting to long format
plot_df <- plot_df %>%
    pivot_longer(
        cols = c("gatk", "full"),
        names_to = "variant_set",
        values_to = "best_r2"
    )
# Making plot
order <- c("gatk", "full")
plot_df$variant_set <- factor(plot_df$variant_set, levels = order)
colors <- c(
    "increase" = "green",
    "no change" = "grey",
    "decrease" = "red"
)
ggplot(plot_df, aes(x = variant_set, y = best_r2, color = category, group = gene)) + 
    geom_line(linewidth = 1) +
    geom_point(size = 4) +
    scale_color_manual(values = colors, name = "Change") +
    labs(title = "Change in R2 from GATK to Full Variant Set", x = "Variant Set (w/ centromere haplotypes)", y = "Best R2") +
    theme_minimal() +
    theme(
        plot.title = element_text(size = 20), # Main title
        axis.title = element_text(size = 18),               # Both X and Y axis labels
        axis.text = element_text(size = 16)                 # Both X and Y tick values
    )

# Want to make line + dot plot that shows change in R2 difference w/ and w/out centromeres considered between the GATK and full variant sets
# Make separate plots for those w/ better/worse performance w/ centromeres using the GATK set
# Color lines/dots depending on if the R2 difference increased or decreased w/ the full variant set
# Better w/ centromeres w/ the GATK call set:
# Making the plot df: Need:
# - gene id
# - Best R2 w/ centromeres using GATK (and best model)
# - Best R2 w/out centromeres using GATK (and best model)
# - Best R2 w/ centromeres using full (and best model)
# - Best R2 w/out centromeres using GATK (and best model)
# --> use to get the R2 diff using GATK and R2 diff using full + category (R2 diff increase, no change, or decrease)
genes <- list()
best_r2_gatk_cent <- list()
best_model_gatk_cent <- list()
best_r2_gatk_no_cent <- list()
best_model_gatk_no_cent <- list()
best_r2_full_cent <- list()
best_model_full_cent <- list()
best_r2_full_no_cent <- list()
best_model_full_no_cent <- list()
for (gene in better_w_t1_df %>% pull(gene)){
    genes <- append(genes, gene)
    best_r2_gatk_cent <- append(best_r2_gatk_cent, cent_t1_df[gene, "best_r2"])
    best_model_gatk_cent <- append(best_model_gatk_cent, cent_t1_df[gene, "best_model"])
    best_r2_gatk_no_cent <- append(best_r2_gatk_no_cent, no_cent_t1_df[gene, "best_r2"])
    best_model_gatk_no_cent <- append(best_model_gatk_no_cent, no_cent_t1_df[gene, "best_model"])
    best_r2_full_cent <- append(best_r2_full_cent, t2_w_cent_results[gene, "best_rsq"])
    best_model_full_cent <- append(best_model_full_cent, t2_w_cent_results[gene, "best_model"])
    best_r2_full_no_cent <- append(best_r2_full_no_cent, t2_wo_cent_results[gene, "best_rsq"])
    best_model_full_no_cent <- append(best_model_full_no_cent, t2_wo_cent_results[gene, "best_model"])
}
plot_df <- data.frame(
    gene = unlist(genes),
    best_r2_gatk_cent = unlist(best_r2_gatk_cent),
    best_model_gatk_cent = unlist(best_model_gatk_cent),
    best_r2_gatk_no_cent = unlist(best_r2_gatk_no_cent),
    best_model_gatk_no_cent = unlist(best_model_gatk_no_cent),
    best_r2_full_cent = unlist(best_r2_full_cent),
    best_model_full_cent = unlist(best_model_full_cent),
    best_r2_full_no_cent = unlist(best_r2_full_no_cent),
    best_model_full_no_cent = unlist(best_model_full_no_cent)
)
plot_df$gatk <- plot_df$best_r2_gatk_cent - plot_df$best_r2_gatk_no_cent
plot_df$full <- plot_df$best_r2_full_cent - plot_df$best_r2_full_no_cent
plot_df$difference_in_r2_diffs <- plot_df$full - plot_df$gatk
plot_df <- plot_df %>%
    rowwise() %>%
    mutate(category = (function(x) {
        if (x > 0){
            return("increase")
        }
        else if (x == 0) {
           return("no change")
        }
        else {
           return("decrease")
        }
    })(difference_in_r2_diffs)) %>%
    ungroup()
# Converting to long format
plot_df <- plot_df %>%
    pivot_longer(
        cols = c("gatk", "full"),
        names_to = "variant_set",
        values_to = "r2_diff"
    )
# Making plot
order <- c("gatk", "full")
plot_df$variant_set <- factor(plot_df$variant_set, levels = order)
colors <- c(
    "increase" = "green",
    "no change" = "grey",
    "decrease" = "red"
)
ggplot(plot_df, aes(x = variant_set, y = r2_diff, color = category, group = gene)) + 
    geom_line(linewidth = 1) +
    geom_point(size = 4) +
    geom_text(
        data = subset(plot_df, variant_set == "full"),
        aes(label = gene, color = "black"),
        hjust = -0.1,
        size = 4,
        show.legend = FALSE
    ) +
    scale_color_manual(values = colors, name = "Change") +
    labs(title = "Change in R2 difference w/ and w/out centromere haplotypes from GATK to full variant set", x = "Variant Set", y = "Difference in R2 w/ and w/out centromere haplotypes") +
    theme_minimal() +
    theme(
        plot.title = element_text(size = 20),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16)
    )

# Better w/out centromeres w/ the GATK call set:
# Making the plot df:
genes <- list()
best_r2_gatk_cent <- list()
best_model_gatk_cent <- list()
best_r2_gatk_no_cent <- list()
best_model_gatk_no_cent <- list()
best_r2_full_cent <- list()
best_model_full_cent <- list()
best_r2_full_no_cent <- list()
best_model_full_no_cent <- list()
for (gene in better_wo_t1_df %>% pull(gene)){
    genes <- append(genes, gene)
    best_r2_gatk_cent <- append(best_r2_gatk_cent, cent_t1_df[gene, "best_r2"])
    best_model_gatk_cent <- append(best_model_gatk_cent, cent_t1_df[gene, "best_model"])
    best_r2_gatk_no_cent <- append(best_r2_gatk_no_cent, no_cent_t1_df[gene, "best_r2"])
    best_model_gatk_no_cent <- append(best_model_gatk_no_cent, no_cent_t1_df[gene, "best_model"])
    best_r2_full_cent <- append(best_r2_full_cent, t2_w_cent_results[gene, "best_rsq"])
    best_model_full_cent <- append(best_model_full_cent, t2_w_cent_results[gene, "best_model"])
    best_r2_full_no_cent <- append(best_r2_full_no_cent, t2_wo_cent_results[gene, "best_rsq"])
    best_model_full_no_cent <- append(best_model_full_no_cent, t2_wo_cent_results[gene, "best_model"])
}
plot_df <- data.frame(
    gene = unlist(genes),
    best_r2_gatk_cent = unlist(best_r2_gatk_cent),
    best_model_gatk_cent = unlist(best_model_gatk_cent),
    best_r2_gatk_no_cent = unlist(best_r2_gatk_no_cent),
    best_model_gatk_no_cent = unlist(best_model_gatk_no_cent),
    best_r2_full_cent = unlist(best_r2_full_cent),
    best_model_full_cent = unlist(best_model_full_cent),
    best_r2_full_no_cent = unlist(best_r2_full_no_cent),
    best_model_full_no_cent = unlist(best_model_full_no_cent)
)
plot_df$gatk <- plot_df$best_r2_gatk_cent - plot_df$best_r2_gatk_no_cent
plot_df$full <- plot_df$best_r2_full_cent - plot_df$best_r2_full_no_cent
plot_df$difference_in_r2_diffs <- plot_df$full - plot_df$gatk
plot_df <- plot_df %>%
    rowwise() %>%
    mutate(category = (function(x) {
        if (x > 0){
            return("increase")
        }
        else if (x == 0) {
           return("no change")
        }
        else {
           return("decrease")
        }
    })(difference_in_r2_diffs)) %>%
    ungroup()
# Converting to long format
plot_df <- plot_df %>%
    pivot_longer(
        cols = c("gatk", "full"),
        names_to = "variant_set",
        values_to = "r2_diff"
    )
# Making plot
order <- c("gatk", "full")
plot_df$variant_set <- factor(plot_df$variant_set, levels = order)
colors <- c(
    "increase" = "green",
    "no change" = "grey",
    "decrease" = "red"
)
ggplot(plot_df, aes(x = variant_set, y = r2_diff, color = category, group = gene)) + 
    geom_line(linewidth = 1) +
    geom_point(size = 4) +
    # geom_text( # To add gene id labels
    #     data = subset(plot_df, variant_set == "full"),
    #     aes(label = gene, color = "black"),
    #     hjust = -0.1,
    #     size = 4,
    #     show.legend = FALSE
    # ) +
    scale_color_manual(values = colors, name = "Change") +
    labs(title = "Change in R2 difference w/ and w/out centromere haplotypes from GATK to full variant set", x = "Variant Set", y = "Difference in R2 w/ and w/out centromere haplotypes") +
    theme_minimal() +
    theme(
        plot.title = element_text(size = 20),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16)
    )

### Examining model weights
# Want to look at model weight breakdown for models w/ centromere haplotypes included--> how much of the signal is being driven by the centromere haplotypes
# Will use the full variant models for now (but can also repeat for the GATK variant models--> need to access weights by loading .RDat files for each gene)
# Think it makes sense to look at the weights for the same model type in both cases--> elastic net seems like a good choice (balance between lasso and ridge)
# --> NOTE: e-net is not the best performing model for each gene, but difficult to compare weights for different model types
# Maybe could make sense to look at the model weights for each model type, but will start w/ e-net
# Will categorize weights into non-centromeric, cis-centromeric (centromere haplotypes on same chromosome as gene), and trans-centromeric (centroemere haplotypes on different chromosome as gene)

full_variant_dir <- "/gpfs/gibbs/pi/ycgh/amm422/project/centromere/prediction/full_variant_tests/predictive_modeling"

# Genes w/ improved performance w/ centromere haplotypes (using GATK variants, t1)
weight_dfs <- list()
for (gene_id in better_w_t1_df |> pull(gene)){
    # Getting chr for gene
    chr <- better_w_t1_df |>
        filter(gene == gene_id) |>
        pull(chr)
    # Getting weights w/ and w/out centromeres
    weights_w_cent <- read.table(paste0(full_variant_dir, "/", gene_id, "/", gene_id, ".model.weights.w_cent.tsv"), sep="\t", header=TRUE, row.names = 1)
    weights_wo_cent <- read.table(paste0(full_variant_dir, "/", gene_id, "/", gene_id, ".model.weights.wo_cent.tsv"), sep="\t", header=TRUE, row.names = 1)
    # Adding column to indicate dataset of variant (GATK, 1KSV, Pangenie, edge --> edgeDepth, TR --> tandem repeat, danbing-tk, centromere (cis/trans)))
    datasets <- c()
    chrs <- c()
    for (id in rownames(weights_w_cent)) {
        if (str_detect(id, ":")){
            dataset <- str_split_i(id, ":", i = -1)
            datasets <- c(datasets, dataset)
            chrs <- c(chrs, chr)
        }
        else{
            cur_chr <- str_split_i(id, "_", i = 1)
            chrs <- c(chrs, cur_chr)
            if (cur_chr == chr){
                datasets <- c(datasets, "centromeric (cis)")
            }
            else{
                datasets <- c(datasets, "centromeric (trans)")
            }
        }
    }
    weights_w_cent$dataset <- datasets
    weights_w_cent$chr <- chrs
    merged_weights <- weights_w_cent |>
        rownames_to_column("variant_id") |>
        rename(enet_w_cent = enet) |>
        select(variant_id, enet_w_cent, dataset, chr) |>
        left_join(weights_wo_cent |> 
            rownames_to_column("variant_id") |>
            rename(enet_wo_cent = enet) |>
            select(variant_id, enet_wo_cent)) |>
        replace_na(list(enet_wo_cent = 0))
    # Adding gene id as column
    merged_weights$gene <- gene_id
    weight_dfs[[gene_id]] <- merged_weights
}
weight_df <- bind_rows(weight_dfs)

# Making stacked bar plot--> showing total proportion of weights from each dataset
# Pair w/ and w/out centromeres for each gene
# Getting proportions of weights for each variant id (handling w/ and w/out centromeres separately)
proportion_df <- weight_df |>
    group_by(gene) |>
    mutate(
        prop_weight_w_cent = (abs(enet_w_cent) / sum(abs(enet_w_cent))) * 100,
        prop_weight_wo_cent = (abs(enet_wo_cent) / sum(abs(enet_wo_cent))) * 100
    ) |>
    ungroup()
# Making df for plot
plot_df <- proportion_df |>
    group_by(gene, dataset) |>
    summarize(
        w_cent = sum(prop_weight_w_cent), 
        wo_cent = sum(prop_weight_wo_cent)) |>
    pivot_longer(
        cols = c("w_cent", "wo_cent"),
        names_to = "test",
        values_to = "prop_weight"
    )
# Ordering genes by decreasing difference in R2 between w/ and w/out centromeres using full variant sets
better_w_t2_df <- read.table("/gpfs/gibbs/pi/ycgh/amm422/project/centromere/prediction/full_variant_tests/tables/better_w_cent_results.tsv", sep="\t", header=TRUE)
order <- better_w_t2_df |>
    arrange(-full) |>
    pull(gene)
plot_df$gene <- factor(plot_df$gene, levels = order)
# Ordering datasets
order <- c("GATK", "1KGSV", "Pangenie", "edge", "TR", "centromeric (cis)", "centromeric (trans)")
plot_df$dataset <- factor(plot_df$dataset, levels = order)
# Ordering tests
order <- c("wo_cent", "w_cent")
plot_df$test <- factor(plot_df$test, levels = order)
# Plotting
ggplot(plot_df, aes(x = test, y = prop_weight, fill = dataset)) +
    geom_col(position = "stack", width = 0.7) +
    facet_grid(~ gene, switch = "x") +
    theme_minimal() +
    theme(
        strip.placement = "outside",        # Moves gene labels below the X-axis
        strip.background = element_blank(), # Removes the background boxes for labels
        panel.spacing = unit(0.5, "lines"), # Controls spacing between the main paired groups
        plot.title = element_text(face = "bold", size = 20),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 10)
    ) +
  labs(
    title = "Proportion of Model Weights by Variant Type",
    x = "Gene (tests w/ and w/out centromere haplotypes included)",
    y = "Proportion of Weights (using elastic net model)",
    fill = "Variant Type"
  )

# Genes w/ worse performance w/ centromere haplotypes (using GATK variants, t1)
weight_dfs <- list()
for (gene_id in better_wo_t1_df |> pull(gene)){
    # Getting chr for gene
    chr <- better_wo_t1_df |>
        filter(gene == gene_id) |>
        pull(chr)
    # Getting weights w/ and w/out centromeres
    weights_w_cent <- read.table(paste0(full_variant_dir, "/", gene_id, "/", gene_id, ".model.weights.w_cent.tsv"), sep="\t", header=TRUE, row.names = 1)
    weights_wo_cent <- read.table(paste0(full_variant_dir, "/", gene_id, "/", gene_id, ".model.weights.wo_cent.tsv"), sep="\t", header=TRUE, row.names = 1)
    # Adding column to indicate dataset of variant (GATK, 1KSV, Pangenie, edge --> edgeDepth, TR --> tandem repeat, danbing-tk, centromere (cis/trans)))
    datasets <- c()
    chrs <- c()
    for (id in rownames(weights_w_cent)) {
        if (str_detect(id, ":")){
            dataset <- str_split_i(id, ":", i = -1)
            datasets <- c(datasets, dataset)
            chrs <- c(chrs, chr)
        }
        else{
            cur_chr <- str_split_i(id, "_", i = 1)
            chrs <- c(chrs, cur_chr)
            if (cur_chr == chr){
                datasets <- c(datasets, "centromeric (cis)")
            }
            else{
                datasets <- c(datasets, "centromeric (trans)")
            }
        }
    }
    weights_w_cent$dataset <- datasets
    weights_w_cent$chr <- chrs
    merged_weights <- weights_w_cent |>
        rownames_to_column("variant_id") |>
        rename(enet_w_cent = enet) |>
        select(variant_id, enet_w_cent, dataset, chr) |>
        left_join(weights_wo_cent |> 
            rownames_to_column("variant_id") |>
            rename(enet_wo_cent = enet) |>
            select(variant_id, enet_wo_cent)) |>
        replace_na(list(enet_wo_cent = 0))
    # Adding gene id as column
    merged_weights$gene <- gene_id
    weight_dfs[[gene_id]] <- merged_weights
}
weight_df <- bind_rows(weight_dfs)

# Making stacked bar plot--> showing total proportion of weights from each dataset
# Pair w/ and w/out centromeres for each gene
# Getting proportions of weights for each variant id (handling w/ and w/out centromeres separately)
proportion_df <- weight_df |>
    group_by(gene) |>
    mutate(
        prop_weight_w_cent = (abs(enet_w_cent) / sum(abs(enet_w_cent))) * 100,
        prop_weight_wo_cent = (abs(enet_wo_cent) / sum(abs(enet_wo_cent))) * 100
    ) |>
    ungroup()
# Making df for plot
plot_df <- proportion_df |>
    group_by(gene, dataset) |>
    summarize(
        w_cent = sum(prop_weight_w_cent), 
        wo_cent = sum(prop_weight_wo_cent)) |>
    pivot_longer(
        cols = c("w_cent", "wo_cent"),
        names_to = "test",
        values_to = "prop_weight"
    )
# Ordering genes by decreasing difference in R2 between w/ and w/out centromeres using full variant sets
better_wo_t2_df <- read.table("/gpfs/gibbs/pi/ycgh/amm422/project/centromere/prediction/full_variant_tests/tables/better_wo_cent_results.tsv", sep="\t", header=TRUE)
order <- better_wo_t2_df |>
    arrange(full) |>
    pull(gene)
plot_df$gene <- factor(plot_df$gene, levels = order)
# Ordering datasets
order <- c("GATK", "1KGSV", "Pangenie", "edge", "TR", "centromeric (cis)", "centromeric (trans)")
plot_df$dataset <- factor(plot_df$dataset, levels = order)
# Ordering tests
order <- c("wo_cent", "w_cent")
plot_df$test <- factor(plot_df$test, levels = order)
# Plotting
ggplot(plot_df, aes(x = test, y = prop_weight, fill = dataset)) +
    geom_col(position = "stack", width = 0.7) +
    facet_grid(~ gene, switch = "x") +
    theme_minimal() +
    theme(
        strip.placement = "outside",        # Moves gene labels below the X-axis
        strip.background = element_blank(), # Removes the background boxes for labels
        panel.spacing = unit(0.5, "lines"), # Controls spacing between the main paired groups
        plot.title = element_text(face = "bold", size = 20),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 10)
    ) +
  labs(
    title = "Proportion of Model Weights by Variant Type",
    x = "Gene (tests w/ and w/out centromere haplotypes included)",
    y = "Proportion of Weights (using elastic net model)",
    fill = "Variant Type"
  )