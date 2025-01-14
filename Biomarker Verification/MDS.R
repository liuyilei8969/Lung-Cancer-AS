library(openxlsx)
library(dplyr)
library(stringr)
library(ggplot2)

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


# MDS 降维
dist <- as.matrix(dist(t(ASmarkers),method='euclidean'))
mds <- cmdscale(dist)
mds.data <- data.frame(Sample=rownames(mds),
                       X=mds[,1],
                       Y=mds[,2])

mds.data <- mds.data%>%
  mutate(Group=str_sub(Sample,-1)) # Tumor-Specific Biomarkers
mds.data$Group = ifelse(mds.data$Group == 'T','Tumor','Normal') # Tumor-Specific Biomarkers

mds.data$Group <- annotation[rownames(mds.data),]$Subtype # Subtype-Specific Biomarkers

# 可视化
ggplot(data=mds.data, aes(x=X, y=Y, fill=Group, color=Group)) +
  geom_point(size=2,stroke = 1, shape = 21)


