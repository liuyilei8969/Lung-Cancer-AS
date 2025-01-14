library(openxlsx)
library(dplyr)
library(stringr)
library(ROCR)
library(caret)


# 读取样本信息
annotation <- read.xlsx('/path/to/sample/information.xlsx') # 准备xlsx表格，有样本名称和对应亚型的信息
annotation <- cbind(annotation$Filename,annotation$Subtype,str_sub(annotation$Filename,-1))
colnames(annotation) <- c('Sample_id','Subtype','Group')
annotation <- as.data.frame(annotation)
rownames(annotation) <- annotation$Sample_id

# 读取 AS matrix 文件
AS <- read.table('splicing_matrix.txt', header = True, sep = '\t')
rownames(AS) <- AS$`#event_id`
AS <- AS[,c(-1,-2)]
# AS <- AS[str_sub(colnames(AS), -1) == 'T', ] # Subtype-Specific Biomarkers

# 读取 Biomarker 文件
markers <- read.xlsx('markers.xlsx')
markers <- markers$markers

# 取出 ASmarker 的 PSI
ASmarkers <- AS[markers,]

ASmarkers <- t(apply(ASmarkers, 1, function(x) {
  # 用行均值替换 NA 值
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  return(x)
}))
ASmarkers <- as.data.frame(ASmarkers)

# 去除意义不大的 AS 事件，增加随机取样的可信度
row_sd <- apply(AS, 1, sd, na.rm = TRUE)
na_ratio <- rowMeans(is.na(AS))
AS_filtered <- AS[row_sd > 0.1 & na_ratio < 0.9, ]
AS_filtered <- t(apply(AS_filtered, 1, function(x) {
  # 用行均值替换 NA 值
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  return(x)
}))

# 提取分组信息
Group <- substr(colnames(AS_filtered), -1) # Tumor-Specific Biomarkers
# Group <- annotation[colnames(AS_filtered), ]$Subtype # Subtype-Specific Biomarkers

# 转置数据，使每行成为样本，列成为特征
AS_filtered <- t(AS_filtered)

set.seed(123)

# 利用 ASmarker 绘制 ROC ，并执行5次随机抽样，绘制 ROC
roc_list <- list()
auc_values <- numeric(6)

# 用 ASmarker 做一个5 fold CV
ASmarkers = as.data.frame(t(ASmarkers))
ASmarkers$label = str_sub(rownames(ASmarkers),-1) # Tumor-Specific Biomarkers
#ASmarkers$label = annotation[colnames(ASmarkers), ]$Subtype # Subtype-Specific Biomarkers
ASmarkers$label = as.factor(ASmarkers$label)
ASmarkers$label =ifelse(ASmarkers$label=='N',0,1) # Tumor-Specific Biomarkers
#ASmarkers$label =ifelse(ASmarkers$label=='A',0,1) # Subtype-Specific Biomarkers

folds <- createFolds(1:356, k=5)

# 初始化变量以保存预测和实际值
all_predictions <- c()
all_actuals <- c()
# 在每一折交叉验证中进行训练和预测
for (fold in folds) {
  train_data <- ASmarkers[-fold, ]
  test_data <- ASmarkers[fold, ]
  
  # 训练 KNN 模型
  pre <- knn(train = train_data[, c(1:24)], test = test_data[, c(1:24)], cl = train_data[, c(25)], k = 7, prob = TRUE)   # Tumor-Specific Biomarkers
  #pre <- knn(train = train_data[, c(1:8)], test = test_data[, c(1:8)], cl = train_data[, c(9)], k = 7, prob = TRUE)  # Subtype-Specific Biomarkers
  
  # 预测概率
  prob <- attr(pre, "prob")
  prob <- 2 * ifelse(pre == "0", 1 - prob, prob) - 1
  all_predictions <- c(all_predictions, prob)
  all_actuals <- c(all_actuals, test_data[, c(25)])
  #all_actuals <- c(all_actuals, test_data[, c(8)])
}

# 计算ROC和AUC
pred_knn <- prediction(all_predictions, all_actuals)
pred_knn <- performance(pred_knn, "tpr", "fpr")
roc_list[[1]] <- pred_knn
result <- roc(response = all_actuals, predictor = all_predictions)
auc_values[1] <- result$auc

set.seed(7)

for (i in 1:5) {
  # 随机选择特征
  selected_features <- sample(ncol(AS_filtered), 24) # Tumor-Specific Biomarkers
  #selected_features <- sample(ncol(AS_filtered), 8) # Subtype-Specific Biomarkers
  AS_filtered_selected <- AS_filtered[, selected_features]
  
  # 将分组信息添加到数据框中
  AS_filtered_selected <- data.frame(AS_filtered_selected)
  AS_filtered_selected$label <- Group
  
  # 确保目标变量是因子
  AS_filtered_selected$label <- as.factor(AS_filtered_selected$label)
  AS_filtered_selected$label =ifelse(AS_filtered_selected$label=='N',0,1) # Tumor-Specific Biomarkers
  #AS_filtered_selected$label =ifelse(AS_filtered_selected$label=='A',0,1) # Subtype-Specific Biomarkers
  
  # 使用 K 折交叉验证
  folds <- createFolds(1:356, k=5)
  
  # 初始化变量以保存预测和实际值
  all_predictions <- c()
  all_actuals <- c()
  
  # 在每一折交叉验证中进行训练和预测
  for (fold in folds) {
    train_data <- AS_filtered_selected[-fold, ]
    test_data <- AS_filtered_selected[fold, ]
    
    # 训练 KNN 模型
    pre <- knn(train = train_data[, 1:24], test = test_data[, 1:24], cl = train_data[, 25], k = 7, prob = TRUE) # Tumor-Specific Biomarkers
    #pre <- knn(train = train_data[, c(1:8)], test = test_data[, c(1:8)], cl = train_data[, c(9)], k = 7, prob = TRUE)  # Subtype-Specific Biomarkers
    
    # 预测概率
    prob <- attr(pre, "prob")
    prob <- 2 * ifelse(pre == "0", 1 - prob, prob) - 1
    all_predictions <- c(all_predictions, prob)
    all_actuals <- c(all_actuals, test_data[, 25])
  }
  
  # 计算ROC和AUC
  pred_knn <- prediction(all_predictions, all_actuals)
  pred_knn <- performance(pred_knn, "tpr", "fpr")
  roc_list[[i+1]] <- pred_knn
  result <- roc(response = all_actuals, predictor = all_predictions)
  auc_values[i+1] <- result$auc
}

# 可视化
par(mfrow =c(1,1),omi = c(0.2,0.2,0.2,0.2),bty = 'l',font.lab = 1, font.axis = 1,cex.axis = 1,cex.lab = 1.3)

plot(roc_list[[1]], col = 'firebrick', lwd = 3, xlab = 'False Positive Rate', ylab = 'True Positive Rate', xlim = c(0, 1), ylim = c(0, 1))

for (i in 2:6) {
  plot(roc_list[[i]], add = TRUE, col = 'grey', lwd = 1.5)
}

# 添加斜对角线
abline(a = 0, b = 1, col = 'black', lty = 2)

# 添加图例
legend('bottomright', legend = c(paste('ASmarkers', 'AUC:', round(auc_values[1], 2)),
                                 paste('Random', 1:5, 'AUC:', round(auc_values[2:6], 2))),
       col = c('firebrick',rep('grey', 5) ), lwd = 2)
