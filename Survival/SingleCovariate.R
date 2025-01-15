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

# 初始化结果
result <- data.frame(matrix(ncol = 2,nrow = length(covariates)))
rownames(result) = covariates
colnames(result) = c('Coefficient','pval')

# 遍历每个变量进行单变量 Cox 分析
surv_obj <- Surv(time = survival$OS.time, event = survival$OS)
for (var in covariates){
  cox_model <- coxph(surv_obj ~ survival[,var], data = survival)
  result[var,1] <- log(exp(coef(cox_model)))
  result[var,2] <- summary(cox_model)$coefficients[, "Pr(>|z|)"]
}

# 多重检验校正
result$FDR = p.adjust(result_df$pval,method = 'fdr')
result$sig = ifelse(result_df$FDR <= 0.05,'sig','notsig')

