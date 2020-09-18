# BrainSpan gene expression imputation
&#8194;&#8194; As the most complex organ of the human body, the brain is composed of diverse regions, each consisting of distinct cell types and their respective cellular interactions. Human brain development involves a finely-tuned cascade of interactive events. These include spatiotemporal gene expression changes and dynamic alterations in cell type composition. However, our understanding of this process is still largely incomplete due to the difficulty of brain spatiotemporal transcriptome collection. In this study, we developed a tensor-based approach to impute gene expression on a transcriptome-wide level. After rigorous computational benchmarking, we applied our approach to infer missing data points in the widely used BrainSpan resource and completed the entire grid of spatiotemporal transcriptomics. Next, we conducted deconvolutional analyses to comprehensively characterize major cell type dynamics across the entire BrainSpan resource to estimate the cellular temporal changes and distinct neocortical areas across development. Moreover, integration of these results with GWAS summary statistics for 13 brain associated traits revealed multiple novel trait-cell type associations and trait-spatiotemporal relationships. In summary, our imputed BrainSpan transcriptomics provide a unique and valuable resource for the research community and our findings significantly extend the current understanding of transcriptional and cellular dynamics during spatiotemporal development of the human brain and it is responsible for cognition, behavior and neuropsychiatric disorders.

# 1. Raw BrainSpan data
&#8194;&#8194; Raw BrainSpan RNAseq data is available from https://www.brainspan.org/api/v2/well_known_file_download/267666525

# 2. Tensor construction and imputation
&#8194;&#8194; We transformed the BrainSpan transcriptome data with missing values into a 3D tensor X, with axes corresponding to the individual, region, and protein-coding gene, respectively. Since all individuals have their corresponding ages, the axis for individuals also represented the dimension of lifespan. Thus, the BrainSpan data imputation problem was transformed as a tensor decomposition and completion task. We applied the Bayesian tensor decomposition method using the trilinear CANDECOMP/PARAFAC (CP) factorization (Khan and Ammad-ud-din 2016) algorithm. The CP method factorized an input tensor into a low-dimensional component space U, V and W vector (rank-1 tensor), corresponding to the individual, region and gene tensor respectively. The number of components, denoted by R, was tested from 10 to 100. In our case, R was automatically determined and optimized by the package. The factorized tensors were optimized to approximate the measured data by minimizing reconstructed variance. R is a positive integer which represents the rank of the tensor and u_r,v_r,w_r denote rank-1 tensor with appropriate dimensions. Here, the notation ‘∘’ represents the outer product of tensors.

&#8194;&#8194; After model fitting, we used the factorized tensors to reconstruct the original tensor, which included not only the approximation of the measured samples, but also the newly generated values for the originally missing samples. 

&#8194;&#8194; After downloading the raw BrainSpan data, users can load them into R environment.  
&#8194;&#8194;&#8194;&#8194;&#8194;&#8194;  `>  BrainSpan_FPKM = read.csv("expression_matrix.csv", head = F, row.names = 1)`  
&#8194;&#8194;&#8194;&#8194;&#8194;&#8194;  `>  BrainSpan_row_info = read.csv("rows_metadata.csv")`  
&#8194;&#8194;&#8194;&#8194;&#8194;&#8194;  `>  BrainSpan_col_info = read.csv("columns_metadata.csv")`  
&#8194;&#8194; The prepocess and imputation Rscript is avaiable at Tensor_imputation/Tensor_imputation.R. In this study, only protein coding genes were selected for downstream analysis. Each round imputation will take ~3 days when using single thread on Intel(R) Xeon(R) Platinum 8276L CPU. Here, as an example, you can only use the first 500 genes for testing.  
&#8194;&#8194; The models' computational complexity is O (K^3), and given that K is much smaller than the number of dimensions, the models is practically useful for K values up to a few hundred. The tensor models are highly suitable for speedup using GPU and the tensorBF model can be scaled well with parallel computing as well. Moreover, the model is not memory intensive. For reasonable values of K, the number of parameters to be learned is much smaller than the size of tensor.

# 3. Evaluation
&#8194;&#8194; We performed Leave-One-Out (LOO) cross-validation to evaluate the results. For each sample with measured transcriptome data, we constructed a tensor excluding this sample (i.e., the holdout sample), applied CP factorization, and imputed the missing transcriptome. Imputation performance was evaluated using four measurements: Pearson correlation coefficient (PCC), R-squared (R2), Root Mean Squared Error (RMSE), and Mean Absolute Error (MAE). We provide R codes at folder Tensor_imputation/Tensor_evaluation.R for imputation performance evaluation.

# 4. Imputation result
&#8194;&#8194; To obtain robust results, we repeated the imputation procedure 100 times. For each gene, we used the median of 100 imputed values as final results. We provide the final completed BrainSpan data at folder Tensor_imputation/BrainSpan_imputation.zip. 

# 5. Cell type deconvolution analysis
&#8194;&#8194; To repeat our cell type deconvolution analysis results, user can download the cell type signatures from Cell_type_deconvolution folder, then apply CIBERSORT (https://cibersort.stanford.edu/).

### Requirements
The R script rely on necessary R package, tensorBF, caret.

## Citation
Pei G, Wang Y, Simon L, Dai Y, Zhao Z, Jia P. Gene expression imputation and cell type deconvolution in human brain with spatiotemporal precision and its implications for brain-related disorders. Genome Research. Under review
