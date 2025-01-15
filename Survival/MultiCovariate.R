library("survival")
library("survminer")
library("openxlsx")

# 加载数据，有样本名称，变量值，OS，OS.time
survival <- read.xlsx('survival.xlsx')
rownames(survival) <- survival$Sample

# 去除 NA 值过多的变量
survival <- survival[ ,which(colMeans(!is.na(survival)) > 0.8) ]
covariates <- colnames(survival)
covariates = covariates[,-c('Sample','OS','OS.time')]

# 遍历每个变量进行单变量 Cox 分析
surv_obj <- Surv(time = survival$OS.time, event = survival$OS)
cox_model <- coxph(surv_obj ~., data = survival[,-c('Sample','OS','OS.time')])

# 处理结果
result <- cbind(coef(cox_model),summary(cox_model)$coefficients[, "Pr(>|z|)"])
rownames(result) <- covariates
colnames(result) <- c('Coefficient','pval')
result <- as.data.frame(result)
result$FDR <- p.adjust(result$pval,method = 'fdr')
result$sig <- ifelse(result_df$FDR <= 0.05,'sig','notsig')
