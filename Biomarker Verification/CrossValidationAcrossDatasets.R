library(openxlsx)
library(dplyr)
library(stringr)
library(ROCR)
library(caret)

# Read sample information
annotation <- read.xlsx('/path/to/sample/information.xlsx') # Prepare an xlsx file with sample names and corresponding subtypes
annotation <- cbind(annotation$Filename, annotation$Subtype, str_sub(annotation$Filename, -1))
colnames(annotation) <- c('Sample_id', 'Subtype', 'Group')
annotation <- as.data.frame(annotation)
rownames(annotation) <- annotation$Sample_id

# Read the AS matrix file
AS <- read.table('splicing_matrix.txt', header = TRUE, sep = '\t')
rownames(AS) <- AS$`#event_id`
AS <- AS[, c(-1, -2)]
# AS <- AS[str_sub(colnames(AS), -1) == 'T', ] # Subtype-Specific Biomarkers

# Read Biomarker file
markers <- read.xlsx('markers.xlsx')
markers <- markers$markers

# Extract PSI for AS markers
ASmarkers <- AS[markers,]

ASmarkers <- t(apply(ASmarkers, 1, function(x) {
  # Replace NA values with row mean
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  return(x)
}))
ASmarkers <- as.data.frame(ASmarkers)

set.seed(123)

# Plot ROC using ASmarkers and perform 5 random samplings to plot ROC
roc_list <- list()
auc_values <- numeric(6)

# Perform 5-fold CV with ASmarkers
ASmarkers = as.data.frame(t(ASmarkers))
ASmarkers$label = str_sub(rownames(ASmarkers), -1) # Tumor-Specific Biomarkers
#ASmarkers$label = annotation[colnames(ASmarkers), ]$Subtype # Subtype-Specific Biomarkers
ASmarkers$label = as.factor(ASmarkers$label)
ASmarkers$label = ifelse(ASmarkers$label == 'N', 0, 1) # Tumor-Specific Biomarkers
#ASmarkers$label = ifelse(ASmarkers$label == 'A', 0, 1) # Subtype-Specific Biomarkers

# Similarly, read and process AS matrix and sample information for CCLE and TCGA, storing them as Validation_ASmarkers
# ......

set.seed(321) 

train_data <- ASmarkers
test_data <- Validation_ASmarkers

# Train KNN model
pre <- knn(train = train_data[, c(1:24)], test = test_data[, c(1:24)], cl = train_data[, c(25)], k = 7, prob = TRUE)   # Tumor-Specific Biomarkers
#pre <- knn(train = train_data[, c(1:8)], test = test_data[, c(1:8)], cl = train_data[, c(9)], k = 7, prob = TRUE)  # Subtype-Specific Biomarkers

# Prediction probabilities
prob <- attr(pre, "prob")
prob <- 2 * ifelse(pre == "0", 1 - prob, prob) - 1
all_predictions <- c(all_predictions, prob)
all_actuals <- c(all_actuals, test_data[, c(25)])
#all_actuals <- c(all_actuals, test_data[, c(8)])

# Calculate ROC and AUC
pred_knn <- prediction(all_predictions, all_actuals)
pred_knn <- performance(pred_knn, "tpr", "fpr")
roc_result <- pred_knn
result <- roc(response = all_actuals, predictor = all_predictions)
auc_values <- result$auc

# Visualization
par(mfrow =c(1,1), omi = c(0.2,0.2,0.2,0.2), bty = 'l', font.lab = 1, font.axis = 1, cex.axis = 1, cex.lab = 1.3)

plot(roc_result, col = 'firebrick', lwd = 3, xlab = 'False Positive Rate', ylab = 'True Positive Rate', xlim = c(0, 1), ylim = c(0, 1))

# Add diagonal line
abline(a = 0, b = 1, col = 'black', lty = 2)

# Add legend
legend('bottomright', legend = c(paste('ASmarkers', 'AUC:', round(auc_values, 3))),
       col = c('firebrick'), lwd = 2)
