library("survival")
library("survminer")
library("openxlsx")

lung = read.xlsx('survival.xlsx') # xlsx中有分组信息，OS，OS.time

fit <- survfit(Surv(OS.time, OS) ~ PSI, data = lung)
summary(fit)

gg <- ggsurvplot(fit,
           pval = TRUE, 
           linetype = "strata")
