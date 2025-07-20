library("survival")
library("survminer")
library("openxlsx")

# Load data, results from MultiCovariate.R
result <- read.xlsx('MultiCovariate.xlsx')
result <- result[order(result$FDR), ]  # Sort the results by FDR

# Select variables
num_selected <- 2 # Input a number to specify how many variables to select
survival_selected <- survival[, rownames(result[1:num_selected,])]  # Select the variables based on row names in 'result'
survival_selected$OS <- survival$OS  # Add survival status to the data
survival_selected$OS.time <- survival$OS.time  # Add survival time to the data

# Build the prognostic model using Cox proportional hazards regression
coxphdata <- coxph(Surv(OS.time, OS) ~ ., data = survival_selected)  # Fit Cox model
coefficients <- as.matrix(coxphdata$coefficients)  # Extract coefficients
survival_selected <- survival_selected[, 1:num_selected]  # Select the first 'num_selected' columns for the model
survival_selected <- as.matrix(survival_selected)  # Convert to matrix
risk <- t(coefficients) %*% t(survival_selected)  # Calculate risk score

# Save the result to an Excel file
write.xlsx(risk, 'risk.xlsx')
