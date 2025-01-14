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

# 同样方法读取并处理 CCLE 和 TCGA 的 ASmatrix 和样本信息，存储为 Validation_ASmarkers
# ......

set.seed(321) 

train_data <- ASmarkers
test_data <- Validation_ASmarkers
  
# 训练 KNN 模型
pre <- knn(train = train_data[, c(1:24)], test = test_data[, c(1:24)], cl = train_data[, c(25)], k = 7, prob = TRUE)   # Tumor-Specific Biomarkers
#pre <- knn(train = train_data[, c(1:8)], test = test_data[, c(1:8)], cl = train_data[, c(9)], k = 7, prob = TRUE)  # Subtype-Specific Biomarkers
  
# 预测概率
prob <- attr(pre, "prob")
prob <- 2 * ifelse(pre == "0", 1 - prob, prob) - 1
all_predictions <- c(all_predictions, prob)
all_actuals <- c(all_actuals, test_data[, c(25)])
#all_actuals <- c(all_actuals, test_data[, c(8)])

# 计算ROC和AUC
pred_knn <- prediction(all_predictions, all_actuals)
pred_knn <- performance(pred_knn, "tpr", "fpr")
roc_result <- pred_knn
result <- roc(response = all_actuals, predictor = all_predictions)
auc_values <- result$auc

# 可视化
par(mfrow =c(1,1),omi = c(0.2,0.2,0.2,0.2),bty = 'l',font.lab = 1, font.axis = 1,cex.axis = 1,cex.lab = 1.3)

plot(roc_result, col = 'firebrick', lwd = 3, xlab = 'False Positive Rate', ylab = 'True Positive Rate', xlim = c(0, 1), ylim = c(0, 1))

# 添加斜对角线
abline(a = 0, b = 1, col = 'black', lty = 2)

# 添加图例
legend('bottomright', legend = c(paste('ASmarkers', 'AUC:', round(auc_values, 3))),
       col = c('firebrick'), lwd = 2)
