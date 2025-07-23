library(openxlsx)

# Read DS events
DS <- read.xlsx("/path/to/DS/events.xlsx")
DS <- c(DS$up, DS$dn)
DS <- na.omit(DS)

# Read Boruta results
Boruata_features <- read.xlsx('/path/to/Boruta_features.xlsx')
Boruata_features <- Boruata_features[1]

# Find the intersection
markers <- intersect(Boruata_features, DS)
markers <- as.data.frame(markers)

# Save the results
write.xlsx(markers, 'markers.xlsx')
