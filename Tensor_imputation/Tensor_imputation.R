setwd("F:/4.UThealth_study/23.Brainspan_WGCNA/Manuscript_Github_code/Tensor_imputation")

#Part I  BrainSpan Tensor construction
#==================================================================================
#Raw BrainSpan RNAseq data is available from https://www.brainspan.org/api/v2/well_known_file_download/267666525

BrainSpan_FPKM = read.csv("expression_matrix.csv", head = F, row.names = 1)
BrainSpan_row_info = read.csv("rows_metadata.csv")
BrainSpan_col_info = read.csv("columns_metadata.csv")

load("BrainSpan_protein_coding_info.rda")

### Filter protein coding genes
BrainSpan_protein_coding_FPKM <- BrainSpan_FPKM[match(Protein_coding_info$ensembl_gene_id, BrainSpan_row_info$ensembl_gene_id),]

### Construct Meta information matrix 
BrainSpan_MetaInfo = matrix(0, nrow = length(unique(BrainSpan_col_info$donor_name)), ncol = length(unique(BrainSpan_col_info$structure_acronym)))

rownames(BrainSpan_MetaInfo) = unique(BrainSpan_col_info$donor_name)
colnames(BrainSpan_MetaInfo) = unique(BrainSpan_col_info$structure_acronym)

for (i in 1:nrow(BrainSpan_MetaInfo)){
	for (j in 1:ncol(BrainSpan_MetaInfo)){ 
		if (colnames(BrainSpan_MetaInfo)[j] %in% BrainSpan_col_info[which(BrainSpan_col_info$donor_name %in% rownames(BrainSpan_MetaInfo)[i]),]$structure_acronym) { 
			BrainSpan_MetaInfo[i, j] = 1
		}
	}
}

### Construct tensor structure
BrainSpan_tensor <- array(NA, dim = c(nrow(BrainSpan_MetaInfo), ncol(BrainSpan_MetaInfo), nrow(BrainSpan_protein_coding_FPKM)))
dim(BrainSpan_tensor)

### Import data to tensor
transform_coordinate = matrix(0, nrow = length(BrainSpan_col_info$donor_name), ncol = 2)

for (i in 1:length(BrainSpan_col_info$donor_name)){
		transform_coordinate[i,1] = match(BrainSpan_col_info$donor_name[i], rownames(BrainSpan_MetaInfo))
		transform_coordinate[i,2] = match(BrainSpan_col_info$structure_acronym[i], colnames(BrainSpan_MetaInfo))
}

for (i in 1:nrow(BrainSpan_protein_coding_FPKM)){
	array_temp = matrix(NA, nrow(BrainSpan_MetaInfo), ncol(BrainSpan_MetaInfo))
	for (j in 1:nrow(transform_coordinate)){
		array_temp[transform_coordinate[j, ][1], transform_coordinate[j, ][2]] = BrainSpan_protein_coding_FPKM[i, j]
	}
	print (i)
	BrainSpan_tensor[, , i] = array_temp
}

rownames(array_temp) = rownames(BrainSpan_MetaInfo)
colnames(array_temp) = colnames(BrainSpan_MetaInfo)

### Names of tensor 
tensor_names_dimA = rownames(BrainSpan_MetaInfo)
tensor_names_dimB = colnames(BrainSpan_MetaInfo)
tensor_names_dimC = Protein_coding_info$ensembl_gene_id

### Save BrainSpan data into tensor structure for downstream imputation
save(BrainSpan_tensor, tensor_names_dimA, tensor_names_dimB, tensor_names_dimC, BrainSpan_protein_coding_FPKM, BrainSpan_MetaInfo, file = "Brainspan_tensor_structure.rda")

#==================================================================================
#Part II  BrainSpan core Tensor construction 

### Filter out individual or region without enougth samples number
load("Brainspan_tensor_structure.rda")

#individual_id_low <- match(c("H376.IIA.50", "H376.IIA.51", names(which(rowSums(BrainSpan_MetaInfo) <= 5))), tensor_names_dimA)
individual_id_low <- c(match(c("H376.IIA.50", "H376.IIA.51"), tensor_names_dimA), as.numeric((which(rowSums(BrainSpan_MetaInfo) <= 5))))
tissue_id_low <- match(names(which(colSums(BrainSpan_MetaInfo) <= 5)), tensor_names_dimB)

BrainSpan_tensor_core <- BrainSpan_tensor[-individual_id_low, -tissue_id_low, ]
tensor_names_dimA_core <- tensor_names_dimA[-individual_id_low]
tensor_names_dimB_core <- tensor_names_dimB[-tissue_id_low]
tensor_names_dimC_core <- tensor_names_dimC

save(BrainSpan_tensor_core, tensor_names_dimA_core, tensor_names_dimB_core, tensor_names_dimC_core, BrainSpan_protein_coding_FPKM, file = "Brainspan_tensor_core_structure.rda")

#==================================================================================
#Part III  Tensor imputation

# Load library for tensor imputation
#install.packages("tensorBF")
library(tensorBF)

load("Brainspan_tensor_core_structure.rda")

# Here we only use the first 500 gene for test.
n = 500 
# If you want to use all genes to repeat our results, please set: n = dim(BrainSpan_tensor_core)[3]
Y = BrainSpan_tensor_core[, , 1:n]

# Model construct: take long time  
# Retry higher K, if Error in eigen(covU) : infinite or missing values in 'x'
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

# Save the tensor data as table file
BrainSpan_core_tensor_imputed = as.data.frame(matrix(nrow = dim(Y)[3], ncol = dim(Y)[1]*dim(Y)[2]))

for (i in 1:dim(Y)[3]){
        BrainSpan_core_tensor_imputed[i, ] = as.numeric(t(Ypred[, , i]))
}

rownames(BrainSpan_core_tensor_imputed) = tensor_names_dimC_core[1:dim(Y)[3]]

for (i in 1:dim(Y)[1]){
        for (j in 1:dim(Y)[2]){
                colnames(BrainSpan_core_tensor_imputed)[(i-1)*dim(Y)[2] + j] = paste0(tensor_names_dimA_core[i], "_", tensor_names_dimB_core[j])
        }
}

# Save imputation results
write.table(BrainSpan_core_tensor_imputed, file = "BrainSpan_imputated.txt", sep = "\t", quote = F)

