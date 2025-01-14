library(Boruta)
library(dplyr)
library(stringr)
library(openxlsx)

# 读取分组信息
annotation <- read.xlsx('/path/to/sample/information.xlsx') # 准备xlsx表格，有样本名称和对应亚型的信息
annotation <- cbind(annotation$Filename,annotation$Subtype,str_sub(annotation$Filename,-1))
colnames(annotation) <- c('Sample_id','Subtype','Group')
annotation <- as.data.frame(annotation)
rownames(annotation) <- annotation$Sample_id
#annotation <- annotation[annotation$Group =='T',] # Subtype-Specific 需要运行此行

# 读取AS矩阵
AS <- read.table('/path/to/AS_matrix.txt', header = TRUE, sep = '\t')
rownames(AS) <- AS$`#event_id`
AS <- AS[,c(-1,-2)]
#AS <- AS[,str_sub(colnames(AS),-1)=='T'] # Subtype-Specific 需要运行此行

# 筛选AS事件
na_ratio <- apply(AS, 1, function(x) mean(is.na(x)))
std_dev <- apply(AS, 1, sd, na.rm = TRUE)
                  
AS_filtered <- AS[na_ratio < 0.7 & std_dev >0.15, ]

# 填充NA
AS_filtered <- t(apply(AS, 1, function(x) {
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  return(x)
}))

AS_filtered <- as.data.frame(t(AS_filtered))

AS_filtered$Group <- annotation$Group
#AS_filtered$Group <- annotation$Subtype # Subtype-Specific 需要运行此行
AS_filtered$Group <- as.factor(AS_filtered$Group)
#AS_filtered = AS_filtered[AS_filtered$Group %in% c('A','S','AS'),] # Subtype-Specific 需要运行此行
 
set.seed(101)
boruta_result <- Boruta(Group ~ ., data = AS_filtered,  mcAdj = TRUE, doTrace = 3)


# 获取所有特征的状态
feature_status <- attStats(boruta_result)

#print(feature_status)
# 提取被确认的重要特征
confirmed_features <- feature_status %>% filter(decision == "Confirmed") %>% filter(normHits==1) %>%rownames()

# 保存
write.xlsx(confirmed_features,'Boruta_features.xlsx')
