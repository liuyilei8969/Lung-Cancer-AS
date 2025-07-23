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
AS <- read.table('AS_matrix.txt', header = TRUE, sep = '\t')
rownames(AS) <- AS$`#event_id`
AS <- AS[, c(-1, -2)]

# Calculate the NA ratio for each row
na_ratio <- apply(AS, 1, function(x) mean(is.na(x)))

# Calculate the standard deviation for each row
std_dev <- apply(AS, 1, sd, na.rm = TRUE)

# Filter the AS matrix based on NA ratio and standard deviation
AS_filtered = AS[na_ratio < 0.7 & std_dev > 0.15, ]

# Impute missing values with the mean for each row
AS_filtered <- t(apply(AS_filtered, 1, function(x) {
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  return(x)
}))

AS_filtered = as.data.frame(t(AS_filtered))

# Add group information
AS_filtered$Group = str_sub(rownames(AS_filtered), -1)
AS_filtered$Group = as.factor(AS_filtered$Group)
AS_filtered$Group = ifelse(AS_filtered$Group == 'N', 0, 1)

# Subtype-Specific 
# AS_filtered$Group = annotation$Subtype
# AS_filtered = AS_filtered[str_sub(rownames(AS_filtered), -1) == 'T', ]
# AS_filtered = AS_filtered[AS_filtered$Group %in% c('A', 'S'), ]
# AS_filtered$Group = as.factor(AS_filtered$Group)
# AS_filtered$Group = ifelse(AS_filtered$Group == 'A', 0, 1)

set.seed(123)
folds <- createFolds(1:356, k = 5)
# folds <- createFolds(1:148, k = 5) # Subtype-Specific                
i = 1
trainlst = list()
markerlst = list()
all_predictions <- c()
all_actuals <- c()

# Split data into folds
for (fold in folds) {
  train_data <- AS_filtered[-fold, ]
  test_data <- AS_filtered[fold, ]
  trainlst[[i]] = rownames(train_data)
  i = i + 1
}
========================================================================================================================================
# Run the script from Biomarker Identification to find features for each fold from the trainlst, saving the results in markerlst
========================================================================================================================================
i = 1
for (fold in folds) {
  train_data <- AS_filtered[-fold, ]
  test_data <- AS_filtered[fold, ]
  
  markers = markerlst[[i]]
  # Train KNN model
  pre <- knn(train = train_data[, c(markers)], test = test_data[, c(markers)], cl = train_data[, c('Group')], k = 7, prob = TRUE)
  # Get prediction probabilities
  prob <- attr(pre, "prob")
  prob <- 2 * ifelse(pre == "0", 1 - prob, prob) - 1
  all_predictions <- c(all_predictions, prob)
  all_actuals <- c(all_actuals, test_data[, c('Group')])
  i = i + 1
}
# Convert prediction probabilities to class predictions
threshold <- 0
predicted_classes <- ifelse(all_predictions > threshold, 1, 0)

# Convert actual labels to numeric (assuming actual labels are 0 or 1)
actual_classes <- as.numeric(all_actuals)
# Compute confusion matrix
conf_matrix <- table(Predicted = predicted_classes, Actual = actual_classes)

# Compute evaluation metrics
accuracy <- sum(diag(conf_matrix)) / sum(conf_matrix)
precision <- conf_matrix[2, 2] / sum(conf_matrix[, 2])
recall <- conf_matrix[2, 2] / sum(conf_matrix[2, ])
f1_score <- 2 * precision * recall / (precision + recall)

# Output the results
cat("Accuracy:", accuracy, "\n")
cat("Precision:", precision, "\n")
cat("Recall:", recall, "\n")
cat("F1-Score:", f1_score, "\n")

# ROC curve and AUC
pred_knn <- prediction(all_predictions, all_actuals)
pred_knn <- performance(pred_knn, "tpr", "fpr")
plot(pred_knn, col = 'firebrick', lwd = 1.5)
result <- roc(response = all_actuals, predictor = all_predictions, direction = "<")
auc_value <- auc(result)

# Precision-Recall curve and AUC
pr = pr.curve(scores.class0 = all_predictions, weights.class0 = all_actuals, curve = TRUE)
prauc_value <- pr$auc.integral
plot(pr, col = 'firebrick', lwd = 1.5)
