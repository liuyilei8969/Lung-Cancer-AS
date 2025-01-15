library("survival")
library("survminer")
library("openxlsx")

# 加载数据，MultiCovariate.R 的结果
result <- read.xlsx('MultiCovariate.xlsx')
result <- result[order(result$FDR), ]

# 挑选变量
num_selected <- 2 # 输入一个数字
survival_selected <- survival[,rownames(result[1:num_selected,])]
survival_selected$OS <- survival$OS
survival_selected$OS.time <- survival$OS.time

# 建立预后模型
coxphdata = coxph(Surv(OS.time,OS)~.,data=survival_selected) 
coefficients = as.matrix(coxphdata$coefficients)
survival_selected = survival_selected[,1:num_selected]
survival_selected = as.matrix(survival_selected)
risk = t(coefficients) %*% t(survival_selected)

# 保存结果
write.xlsx(risk,'risk.xlsx')
