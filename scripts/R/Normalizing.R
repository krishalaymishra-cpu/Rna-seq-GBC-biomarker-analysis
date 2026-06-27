### FOR NORMALIZING WE USED DESeq2 AND FOR FEATURE SLEECTION WE USED ELASTICNET REGRESSION MODEL

#NORMALIZING THE DATA
count_mat<-read.csv("final_filtered_counts.csv")
col_mat<-read.csv("col_data.csv")
library(DESeq2)
dds<-DESeqDataSetFromMatrix(countData = count_mat,
                            colData = col_mat,
                            design = ~gender + condition)
dds<-DESeq(dds)
vsd<-vst(dds, blind = FALSE)
norm_counts<-assay(vsd)
View(norm_counts)    # gene IDs

#FEATURE_SELECTION FOR DOWNREGULATED
downreg<-read.csv("common_downreg_genes.csv")
downreg$X=NULL
colnames(downreg)[1]="genes"
rownames(downreg)<-downreg[,1]
downreg$genes=NULL
length(intersect(rownames(norm_counts), rownames(downreg)))
# Subset properly and force matrix structure
X <- norm_counts[rownames(norm_counts) %in% rownames(downreg), , drop = FALSE]

# Check dimensions (should be genes × samples)
dim(X)

# Transpose for glmnet (samples × genes)
X <- t(X)

# Sanity check
dim(X)                # samples × genes
rownames(X)[1:5]      # sample IDs
colnames(X)[1:5]
#now i have the X dimension as norm_counts
#so now i need the Y dimension
BiocManager::install("glmnet")
library(glmnet)
Y<-factor(col_mat$condition)
#run glmnet
# Run Elastic Net (alpha between 0 and 1, e.g. 0.5)
cvfit <- cv.glmnet(
  x = as.matrix(X),
  y = Y,
  family = "binomial",   # classification
  alpha = 0.5            # Elastic Net balance
)

# Extract selected genes
coef_matrix <- as.matrix(coef(cvfit, s = "lambda.min"))
selected_genes_downreg <- rownames(coef_matrix)[coef_matrix[,1] != 0]
selected_genes_down <- setdiff(selected_genes_downreg, "(Intercept)")

print(selected_genes_down)
length(selected_genes_downreg)
#with confidance
# Coefficient matrix at lambda.min
coef_matrix <- as.matrix(coef(cvfit, s = "lambda.min"))

# Extract selected genes + their coefficients
selected <- coef_matrix[coef_matrix[,1] != 0, , drop = FALSE]
selected <- selected[rownames(selected) != "(Intercept)", , drop = FALSE]

# Gene names
selected_genes <- rownames(selected)

# Gene coefficients (effect sizes)
selected_coefs <- selected[,1]

# Combine into a data frame
selected_features_down <- data.frame(
  Gene = selected_genes,
  Coefficient = selected_coefs
)

# View top few
selected_features_down

#UPREGULATED FEATURE SELECTION
upreg<-read.csv("common_upreg_genes.csv")
upreg$X=NULL
rownames(upreg)=upreg[,1]
length(intersect(rownames(norm_counts),rownames(upreg)))
X<-norm_counts[rownames(norm_counts)%in% rownames(upreg),,drop=FALSE]
dim(X)
X<-t(X)
rownames(X)[1:15]
colnames(X)[1:4]
#as usual i have the X dimension
#for Y
library(glmnet)
Y<-as.factor(col_mat$condition)
#run glm
cvfitup<-cv.glmnet(
  x = as.matrix(X),
  y = Y,
  family = "binomial",
  alpha = 0.45
)

# Extract selected genes
coef_matrix <- as.matrix(coef(cvfitup, s = "lambda.min"))
selected_genes_upreg <- rownames(coef_matrix)[coef_matrix[,1] != 0]
selected_genes_up <- setdiff(selected_genes_upreg, "(Intercept)")

print(selected_genes_up)
length(selected_genes_upreg)
#with confidance
# Coefficient matrix at lambda.min
coef_matrix <- as.matrix(coef(cvfitup, s = "lambda.min"))

# Extract selected genes + their coefficients
selected <- coef_matrix[coef_matrix[,1] != 0, , drop = FALSE]
selected <- selected[rownames(selected) != "(Intercept)", , drop = FALSE]

# Gene names
selected_genes <- rownames(selected)

# Gene coefficients (effect sizes)
selected_coefs <- selected[,1]

# Combine into a data frame
selected_features_up <- data.frame(
  Gene = selected_genes,
  Coefficient = selected_coefs
)

# View top few
selected_features_up
