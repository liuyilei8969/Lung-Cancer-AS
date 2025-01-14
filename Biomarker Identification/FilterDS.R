library(openxlsx)

# How we get Tumor-Specific 127 DS events

AASdiff <- read.table('A_NvsT-splicing_diff.txt',sep = '\t',header = T)
SASdiff <- read.table('S_NvsT-splicing_diff.txt',sep = '\t',header = T)
ASASdiff <- read.table('AS_NvsT-splicing_diff.txt',sep = '\t',header = T)
LASdiff <- read.table('L_NvsT-splicing_diff.txt',sep = '\t',header = T)

AASsigup <- AASdiff[AASdiff$dI_g1_vs_g2 > 0.1 & AASdiff$FDR < 0.01 & AASdiff$coverage > 20,]
ASASsigup <- ASASdiff[ASASdiff$dI_g1_vs_g2 > 0.1 & ASASdiff$FDR < 0.01 & ASASdiff$coverage > 20,]
SASsigup <- SASdiff[SASdiff$dI_g1_vs_g2 > 0.1 & SASdiff$FDR < 0.01 & SASdiff$coverage > 20,]
LASsigup <- LASdiff[LASdiff$dI_g1_vs_g2 > 0.1 & LASdiff$FDR < 0.01 & LASdiff$coverage > 20,]

ASsigup <- Reduce(intersect,list(AASsigup$name,ASASsigup$name,SASsigup$name,LASsigup$name))

AASsigdn <- AASdiff[AASdiff$dI_g1_vs_g2 < -0.1 & AASdiff$FDR <0.01 & AASdiff$coverage >20,]
ASASsigdn <- ASASdiff[ASASdiff$dI_g1_vs_g2 < -0.1 & ASASdiff$FDR <0.01 & ASASdiff$coverage >20,]
SASsigdn <- SASdiff[SASdiff$dI_g1_vs_g2 < -0.1 & SASdiff$FDR <0.01 & SASdiff$coverage >20,]
LASsigdn <- LASdiff[LASdiff$dI_g1_vs_g2 < -0.1 & LASdiff$FDR <0.01 & LASdiff$coverage >20,]

ASsigdn <- Reduce(intersect,list(AASsigdn$name,ASASsigdn$name,SASsigdn$name,LASsigdn$name))

ASsig <- cbind(ASsigup, ASsigdn)
write.xlsx(ASsig,'DS.xlsx')

# How we get Subtype-Specific 178 DS events

AASdiff <- read.table('A-splicing_diff.txt',sep = '\t',header = T)
SASdiff <- read.table('S-splicing_diff.txt',sep = '\t',header = T)

AASsigup <- AASdiff[AASdiff$dI_g1_vs_g2 > 0.1 & AASdiff$FDR < 0.01 & AASdiff$coverage > 20,]
SASsigup <- SASdiff[SASdiff$dI_g1_vs_g2 > 0.1 & SASdiff$FDR < 0.01 & SASdiff$coverage > 20,]

AASsigdn <- AASdiff[AASdiff$dI_g1_vs_g2 < -0.1 & AASdiff$FDR <0.01 & AASdiff$coverage >20,]
SASsigdn <- SASdiff[SASdiff$dI_g1_vs_g2 < -0.1 & SASdiff$FDR <0.01 & SASdiff$coverage >20,]

Asigup <- Reduce(intersect,list(AASsigup$name,SASsigdn$name))
Ssigup <- Reduce(intersect,list(SASsigup$name,AASsigdn$name))

ASsig <- cbind(Asigup, Ssigup)
write.xlsx(ASsig,'DS.xlsx')
