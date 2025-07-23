library(openxlsx)
library(dplyr)
library(stringr)
library(ggplot2)

# Read sample information
annotation <- read.xlsx('/path/to/sample/information.xlsx') # Prepare an xlsx file with sample names and corresponding subtypes
annotation <- cbind(annotation$Filename, annotation$Subtype, str_sub(annotation$Filename, -1))
colnames(annotation) <- c('Sample_id', 'Subtype', 'Group')
annotation <- as.data.frame(annotation)
rownames(annotation) <- annotation$Sample_id

# Read the AS matrix file
AS <- read.table('AS_matrix.txt', header = TRUE, sep = '\t')
rownames(AS) <- AS$`#event_id`
AS <- AS[, c(-1, -2)]
# AS <- AS[str_sub(colnames(AS), -1) == 'T', ] # Subtype-Specific Biomarkers

# Read Biomarker file
markers <- read.xlsx('markers.xlsx')
markers <- markers$markers

# Extract PSI values for AS markers
ASmarkers <- AS[markers,]

ASmarkers <- t(apply(ASmarkers, 1, function(x) {
  # Replace NA values with the row mean
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  return(x)
}))
ASmarkers <- as.data.frame(ASmarkers)

# MDS dimensionality reduction
dist <- as.matrix(dist(t(ASmarkers), method = 'euclidean'))
mds <- cmdscale(dist)
mds.data <- data.frame(Sample = rownames(mds),
                       X = mds[,1],
                       Y = mds[,2])

mds.data <- mds.data %>%
  mutate(Group = str_sub(Sample, -1)) # Tumor-Specific Biomarkers
mds.data$Group = ifelse(mds.data$Group == 'T', 'Tumor', 'Normal') # Tumor-Specific Biomarkers

mds.data$Group <- annotation[rownames(mds.data),]$Subtype # Subtype-Specific Biomarkers

# Visualization
ggplot(data = mds.data, aes(x = X, y = Y, fill = Group, color = Group)) +
  geom_point(size = 2, stroke = 1, shape = 21)
