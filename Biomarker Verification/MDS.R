library(openxlsx)
library(dplyr)
library(stringr)
library(ggplot2)

# 读取 AS matrix 文件
AS <- read.table('splicing_matrix.txt', header = True, sep = '\t')
rownames(AS) <- AS$`#event_id`
AS <- AS[,c(-1,-2)]

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

dist <- as.matrix(dist(t(ASmarkers),method='euclidean'))
mds <- cmdscale(dist)
mds.data <- data.frame(Sample=rownames(mds),
                       X=mds[,1],
                       Y=mds[,2])

mds.data <- mds.data%>%
  mutate(Group=str_sub(Sample,-1)) # Tumor-Specific Biomarkers
mds.data$Group = ifelse(mds.data$Group == 'T','Tumor','Normal') # Tumor-Specific Biomarkers



# 可视化
ggplot(data=mds.data, aes(x=X, y=Y, fill=Group, color=Group)) +
  geom_point(size=2,stroke = 1, shape = 21)


