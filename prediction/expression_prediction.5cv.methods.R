# expression_prediction.5cv.methods.R

### FROM SHUANGJIA

# output: 
# gene.cv.calls.txt - store the predicted values for all samples across folds plus the true phenotype values
# gene.cv.performance.txt - store the cross-validation performance metrics for each method, including R-squared and p-value
# gene.train.performance.txt - store the training performance metrics for each method, including R-squared for each fold
# gene.model.weights.tsv - store the coefficients for each variant for each model (rows: variants, columns: models)

# Usage: Rscript expression_prediction.5cv.methods.R <gene_name> # ENSG00000196476.11

library(susieR)
library(Matrix)
library(Rfast)
library(data.table)
library(glmnet)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

# gene name
gene = args[1]
# gene <- "ENSG00000042445.13"

# Reading in list of variants (w/out centromeres)
variant_df <- read.table(paste0("/gpfs/gibbs/pi/ycgh/amm422/project/centromere/prediction/full_variant_tests/inputs/", gene, "/", gene, ".full_variants.variant_info_no_centromeres.txt"),
                sep = "\t",
                header = FALSE,
                col.names = c("variant_id", "chr", "start")) # Add "end" w/ cents included
row.names(variant_df) <- variant_df$variant_id
variant_df <- variant_df %>% select(c(chr, start)) # Add "end" w/ cents included

# Reading in genotypes (w/out centromeres)
X = read.table(paste(c("/gpfs/gibbs/pi/ycgh/amm422/project/centromere/prediction/full_variant_tests/inputs/", gene, "/", gene, ".full_variants.geno_no_centromeres.tsv"), collapse = ""),
                sep = "\t",
                header = FALSE,
                check.names = FALSE,
                stringsAsFactors = FALSE) # 430 * number of SNPs in the window around the TSS + number of centromere haplotypes

# change inf (in TR callset) to 0 
X[X == Inf] <- 0

# read in covariates
Z = read.table("/gpfs/gibbs/pi/ycgh/amm422/project/centromere/prediction/covariates.tsv",
               sep = "\t",
               header = TRUE,
               row.names = 1,
               check.names = FALSE,
               stringsAsFactors = FALSE) # 430 * 40

# phenotype - gene expression
y = read.table(paste(c("/gpfs/gibbs/pi/ycgh/amm422/project/centromere/prediction/full_variant_tests/inputs/", gene, "/", gene, ".txt"), collapse = ""),
               sep = "\t",
               header = FALSE,
               check.names = FALSE,
               stringsAsFactors = FALSE) # 430 * 1

y <- as.matrix(y)
Z <- as.matrix(Z) # nolint: object_name_linter.
X <- as.matrix(X)

# remove covariates from both genotype and phenotype - residualization
remove.covariate.effects <- function (X, Z, y) {
  # include the intercept term
  if (any(Z[,1]!=1)) Z = cbind(1, Z)
  A   <- forceSymmetric(crossprod(Z))
  SZy <- as.vector(solve(A,c(y %*% Z)))
  SZX <- as.matrix(solve(A,t(Z) %*% X))
  y <- y - c(Z %*% SZy)
  X <- X - Z %*% SZX
  return(list(X = X,y = y,SZy = SZy,SZX = SZX))
}

out = remove.covariate.effects(X, Z, y[,1])

#  important : genotypes are standardized and scaled here:
scale_Xcov = scale(out$X) # samp * var 
# scale phenotype # should have already been scaled before in my analysis; but not after correcting for covariates
scale_ycov = scale(out$y)

# CROSSVALIDATION ANALYSES
set.seed(1)

# create a matrix to store the cross-validation performance metrics for each method
cv.performance = matrix(NA,nrow=10,ncol=1) # 1 col
rownames(cv.performance) = c("lasso_rsq","lasso_pval","ridge_rsq","ridge_pval","enet_rsq","enet_pval","top1_rsq","top1_pval","susie_rsq","susie_pval")
colnames(cv.performance) = c('all_variants')

# create a matrix to store the training performance metrics for each method
train.performance = matrix(NA,nrow=5,ncol=5) # 5 col: gatk * lasso, ridge, elastic net, top1, bayesian
rownames(train.performance) = c('lasso','ridge','enet', 'top1','susie')

# split data into 5 folds
N = 430 
sample = sample(430) # randomly shuffle the sample order 
folds = cut(seq(1,N), breaks=5, labels=FALSE) # assign fold labels to samples

scale_ycov_shuffle = scale_ycov[sample, , drop=FALSE] # shuffle the phenotype according to the sample order
scale_Xcov_shuffle = scale_Xcov[sample, , drop=FALSE] # shuffle the genotype according to the sample order

# create a matrix to store the predicted values for all samples across folds
cv.calls = matrix(NA,nrow=N,ncol=6) # store the predicted values for all samples across folds
colnames(cv.calls) = c("pred_lasso", "pred_ridge", "pred_enet", "pred_top1", "pred_susie", "true_y")

for (i in 1:5){
    test_index = which(folds == i) # get the test index for the current fold
    train_index = which(folds != i) # get the train index for the current fold
    
    X_train = scale_Xcov_shuffle[train_index, , drop=FALSE] # get the training genotype data
    y_train = scale_ycov_shuffle[train_index, , drop=FALSE] # get the training phenotype data

    # # run susie fit on the training data
    # fitted_susie = susie(X_train, y_train, L = 10, standardize = TRUE, compute_univariate_zscore = TRUE)
    # # record training performance rsq for the current fold
    # reg_train = summary(lm( y_train ~ predict(fitted_susie) ))
    # train.performance[ 1, i ] = reg_train$adj.r.sq

    X_test = scale_Xcov_shuffle[test_index, , drop=FALSE] # get the test genotype data
    y_test = scale_ycov_shuffle[test_index, , drop=FALSE] # get the test phenotype data
    # pred_y = predict(fitted_susie, newx = X_test)
    # cv.calls[sample[test_index], 1] = pred_y # store the predicted values for the test samples

    # run lasso 
    set.seed(123)
    lasso_fit <- cv.glmnet(X_train, y_train, alpha=1, family="gaussian")
    # make predictions on the test set using the lasso model
    # record training performance rsq for the current fold
    pred_train_lasso <- predict(lasso_fit, newx = X_train, s = "lambda.min")
    reg_train_lasso = summary(lm( y_train ~ pred_train_lasso ))
    train.performance[ 1, i ] = reg_train_lasso$adj.r.sq

    pred_y_lasso <- predict(lasso_fit, newx = X_test, s = "lambda.min")
    cv.calls[sample[test_index], 1] = pred_y_lasso

    # run ridge regression 
    set.seed(123)
    ridge_fit <- cv.glmnet(X_train, y_train, alpha=0, family="gaussian")
    # make predictions on the test set using the ridge model
    # record training performance rsq for the current fold
    pred_train_ridge <- predict(ridge_fit, newx = X_train, s = "lambda.min")
    reg_train_ridge = summary(lm( y_train ~ pred_train_ridge ))
    train.performance[ 2, i ] = reg_train_ridge$adj.r.sq

    pred_y_ridge <- predict(ridge_fit, newx = X_test, s = "lambda.min")
    cv.calls[sample[test_index], 2] = pred_y_ridge

    # run elastic net
    set.seed(123)
    enet_fit <- cv.glmnet(X_train, y_train, alpha=0.5, family="gaussian")
    # make predictions on the test set using the elastic net model
    # record training performance rsq for the current fold
    pred_train_enet <- predict(enet_fit, newx = X_train, s = "lambda.min")
    reg_train_enet = summary(lm( y_train ~ pred_train_enet ))
    train.performance[ 3, i ] = reg_train_enet$adj.r.sq

    pred_y_enet <- predict(enet_fit, newx = X_test, s = "lambda.min")
    cv.calls[sample[test_index], 3] = pred_y_enet

    # top1 
    eff.wgt = t( X_train ) %*% (y_train) / ( nrow(y_train) - 1)
	eff.wgt[ - which.max( eff.wgt^2 ) ] = 0
    train.performance[ 4, i ] = summary(lm( y_train ~ X_train %*% eff.wgt ))$adj.r.sq
    cv.calls[sample[test_index], 4] = X_test %*% eff.wgt

    # susie --> Different Bayesian model from FUSION
    fitted_susie = susie(X_train, y_train, L = 10, standardize = TRUE, compute_univariate_zscore = TRUE)
    train.performance[ 5, i ] = summary(lm( y_train ~ predict(fitted_susie, newx = X_train) ))$adj.r.sq
    pred_y_susie = predict(fitted_susie, newx = X_test)
    cv.calls[sample[test_index], 5] = pred_y_susie # store the predicted values for the test samples
}

# calculate metrics (e.g. mean squared error) for each fold, average across folds

# compute rsq + P-value
# lasso 
reg = summary(lm( scale_ycov ~ cv.calls[,1] ))
cv.performance[ 1, 1 ] = reg$adj.r.sq
cv.performance[ 2, 1 ] = reg$coef[2,4]
# ridge 
reg = summary(lm( scale_ycov ~ cv.calls[,2] ))
cv.performance[ 3, 1 ] = reg$adj.r.sq
cv.performance[ 4, 1 ] = reg$coef[2,4]
# elastic net
reg = summary(lm( scale_ycov ~ cv.calls[,3] ))
cv.performance[ 5, 1 ] = reg$adj.r.sq
cv.performance[ 6, 1 ] = reg$coef[2,4]
# top1 
reg = summary(lm( scale_ycov ~ cv.calls[,4] ))
cv.performance[ 7, 1 ] = reg$adj.r.sq
cv.performance[ 8, 1 ] = reg$coef[2,4]
# susie
reg = summary(lm( scale_ycov ~ cv.calls[,5] ))
cv.performance[ 9, 1 ] = reg$adj.r.sq
cv.performance[ 10, 1 ] = reg$coef[2,4]

# add y to cv.calls for downstream analyses
cv.calls[, 6] = scale_ycov

# Final refit on all samples, then store weights in a SNP x model matrix
# Rows = SNPs, columns = model weights
# Intercept is excluded from this table

# Make a dataframe with SNPs as row names
final_weights <- data.frame(
  lasso = rep(NA_real_, ncol(scale_Xcov)),
  ridge = rep(NA_real_, ncol(scale_Xcov)),
  enet = rep(NA_real_, ncol(scale_Xcov)),
  top1 = rep(NA_real_, ncol(scale_Xcov)),
  susie = rep(NA_real_, ncol(scale_Xcov)),
  row.names = rownames(variant_df)
)

# Lasso
lasso_final <- cv.glmnet(scale_Xcov, scale_ycov, alpha = 1, family = "gaussian")
lasso_coef <- as.vector(coef(lasso_final, s = "lambda.min"))
lasso_coef <- lasso_coef[-1]   # remove intercept
names(lasso_coef) <- rownames(variant_df)
final_weights$lasso[match(names(lasso_coef), rownames(final_weights))] <- as.numeric(lasso_coef)

# Ridge
ridge_final <- cv.glmnet(scale_Xcov, scale_ycov, alpha = 0, family = "gaussian")
ridge_coef <- as.vector(coef(ridge_final, s = "lambda.min"))
ridge_coef <- ridge_coef[-1]
names(ridge_coef) <- rownames(variant_df)
final_weights$ridge[match(names(ridge_coef), rownames(final_weights))] <- as.numeric(ridge_coef)

# Elastic net
enet_final <- cv.glmnet(scale_Xcov, scale_ycov, alpha = 0.5, family = "gaussian")
enet_coef <- as.vector(coef(enet_final, s = "lambda.min"))
enet_coef <- enet_coef[-1]
names(enet_coef) <- rownames(variant_df)
final_weights$enet[match(names(enet_coef), rownames(final_weights))] <- as.numeric(enet_coef)

# Top1
top1_w <- t(scale_Xcov) %*% scale_ycov / (nrow(scale_Xcov) - 1)
top1_w[-which.max(top1_w^2)] <- 0
names(top1_w) <- rownames(variant_df)
final_weights$top1[match(names(top1_w), rownames(final_weights))] <- as.numeric(top1_w)

# SuSiE
susie_final <- susie(scale_Xcov, scale_ycov, L = 10, standardize = TRUE, compute_univariate_zscore = TRUE)
susie_coef <- coef(susie_final)

# If coef() returns a list, pull out the coefficient vector
if (is.list(susie_coef)) {
  if ("coefficients" %in% names(susie_coef)) {
    susie_coef <- susie_coef$coefficients
  } else if ("beta" %in% names(susie_coef)) {
    susie_coef <- susie_coef$beta
  } else {
    susie_coef <- unlist(susie_coef, use.names = TRUE)
  }
}

# If the vector length is one longer than the number of SNPs, the first element is the intercept
if (length(susie_coef) == ncol(scale_Xcov) + 1L) {
  susie_coef <- susie_coef[-1]
}

# Remove any unnamed/empty-string entries
if (!is.null(names(susie_coef))) {
  keep <- !(names(susie_coef) %in% c("", NA_character_, "(Intercept)"))
  susie_coef <- susie_coef[keep]
}

# Now map to SNP names
names(susie_coef) <- rownames(variant_df)
final_weights$susie[match(names(susie_coef), rownames(final_weights))] <- as.numeric(susie_coef)

# Making output directory if it doesn't already exist
dir.create(paste0("/gpfs/gibbs/pi/ycgh/amm422/project/centromere/prediction/full_variant_tests/predictive_modeling/", gene), showWarnings = FALSE)

# save cv.performance and cv.calls for downstream analyses
# w/ centromeres --> suffix: *.w_cent.txt; w/out centromeres --> suffix: *.wo_cent.txt (run w/out centromeres later)
write.table(cv.calls,paste(c("/gpfs/gibbs/pi/ycgh/amm422/project/centromere/prediction/full_variant_tests/predictive_modeling/", gene, "/",gene, ".cv.calls.wo_cent.txt"), collapse = ""),quote=F,row.names =FALSE, col.names = TRUE, sep='\t')

write.table(cv.performance,paste(c("/gpfs/gibbs/pi/ycgh/amm422/project/centromere/prediction/full_variant_tests/predictive_modeling/", gene, "/",gene, ".cv.performance.wo_cent.txt"), collapse = ""),quote=F,row.names =TRUE,col.names = TRUE,sep='\t')

# save training performance for downstream analyses
write.table(train.performance,paste(c("/gpfs/gibbs/pi/ycgh/amm422/project/centromere/prediction/full_variant_tests/predictive_modeling/", gene, "/",gene, ".train.performance.wo_cent.txt"), collapse = ""),quote=F,row.names =TRUE,col.names = FALSE,sep='\t')

# write the SNP x model weight matrix
write.table(
  final_weights,
  paste0("/gpfs/gibbs/pi/ycgh/amm422/project/centromere/prediction/full_variant_tests/predictive_modeling/",
         gene, "/", gene, ".model.weights.wo_cent.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = NA
)
