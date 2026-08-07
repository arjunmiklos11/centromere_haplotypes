# expression_prediction.5cv.methods.R

### FROM SHUANGJIA

# output: 
# gene.cv.calls.1kg.methods.txt - store the predicted values for all samples across folds plus the true phenotype values
# gene.cv.performance.1kg.methods.txt - store the cross-validation performance metrics for each method, including R-squared and p-value
# gene.train.performance.1kg.methods.txt - store the training performance metrics for each method, including R-squared for each fold

# Usage: Rscript expression_prediction.5cv.methods.R <gene_name> # ENSG00000196476.11

### INSTEAD OF RUNNING THIS SCRIPT, CAN RUN FUSION INSTEAD (Runs the same models to predict expression, but uses different Bayesian approach than the one here (susie))

library(susieR)
library(Matrix)
library(Rfast)
library(data.table)
library(glmnet)

args <- commandArgs(trailingOnly = TRUE)

# gene name 
# gene = args[1] 
gene = "ENSG00000000938.12"

# 1. rescue gatk + 1kgsv prediction 
# Reading in genotypes--> should match the per-gene genotype files outputted by input_generator.py (all variants w/in 1 Mb of TSS--> need to add centromere haplotypes too)
# Path to match produced files from input_generator.py: paste0("/gpfs/gibbs/pi/ycgh/amm422/project/centromere/prediction/inputs/", gene, "/", gene, ".1kg_rescue.geno.tsv")
X = read.table(paste(c("/gpfs/gibbs/pi/ycgh/lushjia/project/SV/AFGR/RNA/hprc_v2/edge/collapse/var_1_to_1_rescue/separate_prediction/", gene, "/",gene, ".1kg_rescue.geno.tsv"), collapse = ""),
                sep = "\t",
                header = FALSE,
                check.names = FALSE,
                stringsAsFactors = FALSE) # 430 * number of SNPs in the window around the TSS + number of centromere haplotypes

# change inf (in TR callset) to 0 
X[X == Inf] <- 0

# read in covariates--> should match the covariate table outputted by input_generator.py (which was unchanged from what was read into input_generator.py)
# Path produced by input_generator.py: /gpfs/gibbs/pi/ycgh/amm422/project/centromere/prediction/covariates.tsv
Z = read.table("/gpfs/gibbs/pi/ycgh/lushjia/project/SV/AFGR/RNA/hprc_v2/edge/susie/covariates.tsv",
               sep = "\t",
               header = TRUE,
               row.names = 1,
               check.names = FALSE,
               stringsAsFactors = FALSE) # 430 * 40

# phenotype - gene expression--> should also match per-gene gene expression output from input_generator.py (just the gene expression of the given gene)
# Path to match produced files from input_generator.py: paste0("/gpfs/gibbs/pi/ycgh/amm422/project/centromere/prediction/inputs/", gene, "/", gene, ".txt")
y = read.table(paste(c("/gpfs/gibbs/pi/ycgh/lushjia/project/SV/AFGR/RNA/hprc_v2/edge/susie/genes/", gene, "/",gene, ".txt"), collapse = ""),
               sep = "\t",
               header = FALSE,
               check.names = FALSE,
               stringsAsFactors = FALSE) # 430 * 1

y <- as.matrix(y)
Z <- as.matrix(Z) # nolint: object_name_linter.
X <- as.matrix(X)
# any(is.infinite(X))

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
### Why are we running these analyses outside of FUSION? They appear to be the same tests?
set.seed(1)

# create a matrix to store the cross-validation performance metrics for each method
cv.performance = matrix(NA,nrow=10,ncol=1) # 1 col
rownames(cv.performance) = c("lasso_rsq","lasso_pval","ridge_rsq","ridge_pval","enet_rsq","enet_pval","top1_rsq","top1_pval","susie_rsq","susie_pval")
colnames(cv.performance) = c('gatk')

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

# save cv.performance and cv.calls for downstream analyses
write.table(cv.calls,paste(c("/gpfs/gibbs/pi/ycgh/lushjia/project/SV/AFGR/RNA/hprc_v2/edge/collapse/var_1_to_1_rescue/separate_prediction/", gene, "/",gene, ".cv.calls.1kg.methods.txt"), collapse = ""),quote=F,row.names =FALSE, col.names = TRUE, sep='\t')

write.table(cv.performance,paste(c("/gpfs/gibbs/pi/ycgh/lushjia/project/SV/AFGR/RNA/hprc_v2/edge/collapse/var_1_to_1_rescue/separate_prediction/", gene, "/",gene, ".cv.performance.1kg.methods.txt"), collapse = ""),quote=F,row.names =TRUE,col.names = TRUE,sep='\t')

# save training performance for downstream analyses
write.table(train.performance,paste(c("/gpfs/gibbs/pi/ycgh/lushjia/project/SV/AFGR/RNA/hprc_v2/edge/collapse/var_1_to_1_rescue/separate_prediction/", gene, "/",gene, ".train.performance.1kg.methods.txt"), collapse = ""),quote=F,row.names =TRUE,col.names = FALSE,sep='\t')



