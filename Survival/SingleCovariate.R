library("survival")
library("survminer")
library("openxlsx")

# Load data, including sample names, variable values, OS, and OS.time
survival <- read.xlsx('survival.xlsx')
rownames(survival) <- survival$Sample  # Set the sample names as row names

# Remove variables with too many NA values (keep those with less than 20% missing data)
survival <- survival[, which(colMeans(!is.na(survival)) > 0.8)]  # Filter out columns with more than 20% missing data
covariates <- colnames(survival)  # Extract column names (covariates)
covariates = covariates[-c('Sample', 'OS', 'OS.time')]  # Remove 'Sample', 'OS', and 'OS.time' from covariates

# Initialize results data frame
result <- data.frame(matrix(ncol = 2, nrow = length(covariates)))  # Create a data frame with 2 columns and rows equal to the number of covariates
rownames(result) = covariates  # Set row names as covariates
colnames(result) = c('Coefficient', 'pval')  # Set column names to 'Coefficient' and 'pval'

# Perform univariate Cox regression for each covariate
surv_obj <- Surv(time = survival$OS.time, event = survival$OS)  # Define the survival object
for (var in covariates) {  # Loop over each covariate
  cox_model <- coxph(surv_obj ~ survival[, var], data = survival)  # Fit Cox regression model
  result[var, 1] <- log(exp(coef(cox_model)))  # Store the log of the coefficient
  result[var, 2] <- summary(cox_model)$coefficients[, "Pr(>|z|)"]  # Store the p-value
}

# Multiple testing correction using FDR
result$FDR = p.adjust(result$pval, method = 'fdr')  # Apply FDR correction to p-values
result$sig = ifelse(result$FDR <= 0.05, 'sig', 'notsig')  # Mark significant variables (FDR ≤ 0.05)
