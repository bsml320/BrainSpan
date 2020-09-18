setwd("F:/4.UThealth_study/23.Brainspan_WGCNA/Manuscript_Github_code/Tensor_imputation")

#Part IV  Tensor imputation evaluation

# Load library for tensor imputation
#install.packages("tensorBF")
library(tensorBF)

load("Brainspan_tensor_core_structure.rda")

# Here we only use the first 500 gene for test.
n = 500 
# If you want to use all genes to repeat our results, please set: n = dim(BrainSpan_tensor_core)[3]
Y = BrainSpan_tensor_core[, , 1:n]

#==================================================
# Please select either plan A or B for downstream imputation.
# Plan 1. Random remove all genes in one samples
random_index = sample(which(!is.na(as.numeric(Y[,,1]))), 1)		
# For leave one out, please try to use 1 to 560 instead of random_index.
missing.inds = dim(Y)[1]*dim(Y)[2]*(0:(dim(Y)[3]-1)) + random_index
Yobs = Y[missing.inds]
#Yobs == Y[(random_index %% dim(Y)[1]), (random_index %/% dim(Y)[1])+1, ]
Y[missing.inds] <- NA

#=================================================
# Plan 2. Random remove genes
# Replace the observed value by missing value
missing.inds = sample(prod(dim(Y)), n)
Yobs = Y[missing.inds]
Y[missing.inds] <- NA
dim(Y)
#==================================================

# Model construct: take long time  
# Retry higher K, if Error in eigen(covU) : infinite or missing values in 'x', or try to use all gene for imputation. 
print (Sys.time())
res <- tensorBF(Y = Y, method = "CP", K = NULL, opts = NULL, fiberCentering = NULL, slabScaling = NULL, noiseProp = c(0.5, 0.5))
print (Sys.time())

# Missing value prediction
Ypred = predictTensorBF(Y = Y, res = res)

# Replace the nagetive value by 0
Ypred[Ypred < 0] <- 0

# Replace the potential outliers introduce by imputation
data_obs <- BrainSpan_protein_coding_FPKM

for (i in 1:dim(Ypred)[3]){
	gene_max <- quantile(data_obs[i, ])[5]
	gene_Q3 <- quantile(data_obs[i,])[4]
	gene_Q1 <- quantile(data_obs[i,])[2]
	gene_threhold <- gene_max + (gene_Q3 - gene_Q1) * 1.0
	Ypred[,,i][Ypred[,,i] > as.numeric(gene_threhold)] <- as.numeric(gene_threhold)
	print (i)
}

# Load caret library for imputation evaluation
#install.packages("caret")
library(caret)

Ypredict <- round(Ypred[missing.inds], 3)
Ymean <- as.numeric(round(rowMeans(BrainSpan_protein_coding_FPKM), 3))[1:length(Ypredict)]

Ypred_obs_compare = data.frame(cbind(Yobs, Ypredict, Ymean))

Corr_Ypred = round(cor(Yobs, Ypredict), 3)
R2_Ypred = round(R2(Yobs, Ypredict), 3)
RMSE_Ypred = round(RMSE(Yobs, Ypredict), 2)
MAE_Ypred = round(MAE(Yobs, Ypredict), 2)

Corr_Ymean = round(cor(Yobs, Ymean), 3)
R2_Ymean = round(R2(Yobs, Ymean), 3)
RMSE_Ymean = round(RMSE(Yobs, Ymean), 2)
MAE_Ymean = round(MAE(Yobs, Ymean), 2)

Ypred_obs_eval = data.frame(cbind(Corr_Ypred, R2_Ypred, RMSE_Ypred, MAE_Ypred, Corr_Ymean, R2_Ymean, RMSE_Ymean, MAE_Ymean))

write.table(Ypred_obs_compare, "Tensor_imputed_comparison", sep = "\t", quote = F)
write.table(Ypred_obs_eval, "Tensor_imputed_evaluation", sep = "\t", quote = F)