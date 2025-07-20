library("survival")
library("survminer")
library("openxlsx")

# Load data, including sample names, variable values, OS, and OS.time
survival <- read.xlsx('survival.xlsx')
rownames(survival) <- survival$Sample  # Set the sample names as row names

# Remove variables with too many NA values (keep those with less than 20% missing data)
survival <- survival[, which(colMeans(!is.na(survival)) > 0.8)]
covariates <- colnames(survival)
covariates = covariates[-c('Sample', 'OS', 'OS.time')]  # Remove 'Sample', 'OS', and 'OS.time' from covariates

# Perform univariate Cox analysis for each variable
surv_obj <- Surv(time = survival$OS.time, event = survival$OS)  # Define survival object
cox_model <- coxph(surv_obj ~ ., data = survival[, -c('Sample', 'OS', 'OS.time')])  # Fit Cox model

# Process the results
result <- cbind(coef(cox_model), summary(cox_model)$coefficients[, "Pr(>|z|)"])  # Extract coefficients and p-values
rownames(result) <- covariates  # Set covariates as row names
colnames(result) <- c('Coefficient', 'pval')  # Assign column names
result <- as.data.frame(result)  # Convert to data frame
result$FDR <- p.adjust(result$pval, method = 'fdr')  # Apply FDR correction
result$sig <- ifelse(result$FDR <= 0.05, 'sig', 'notsig')  # Mark significant variables
