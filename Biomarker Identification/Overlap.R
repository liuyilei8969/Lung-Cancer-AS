library(openxlsx)

# 读取DS事件
DS <- read.xlsx("/path/to/DS/events.xlsx",sheet =2)
DS <- c(DS$up,DS$dn)
DS <- na.omit(DS)

# 读取 Boruta 结果
Boruata_features <- read.xlsx('/path/to/Boruta_features.xlsx')
Boruata_features <- Boruata_features[1]

# 取交集
markers <- intersect(Boruata_features, DS)
markers <- as.data.frame(markers)

write.xlsx(markers,'markers.xlsx')
