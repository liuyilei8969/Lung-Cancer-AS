library("survival")
library("survminer")
library("openxlsx")

# Load data from the 'survival.xlsx' file, which contains grouping information, OS, and OS.time
lung = read.xlsx('survival.xlsx')

# Fit survival curve based on PSI (possibly a variable indicating a certain condition)
fit <- survfit(Surv(OS.time, OS) ~ PSI, data = lung)  # Create a Kaplan-Meier survival fit for OS based on PSI
summary(fit)  # Print the summary of the survival fit

# Plot the survival curve
gg <- ggsurvplot(fit,
           pval = TRUE,  # Show p-value on the plot
           linetype = "strata")  # Different line types for strata (groups)
