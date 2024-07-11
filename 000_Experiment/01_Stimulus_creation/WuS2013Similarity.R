setwd("M:/Arbeit/Programmierung")

#### similarity based on Rating (Wong & Szuecs, 2013) ####
WuS = read.table("Wong&Szuecs2013AppendixA.txt", header = TRUE, sep = "")

WuS$Mean = (WuS$C + WuS$A1+ WuS$A2)/3

plot(WuS$D1, WuS$Mean)

WuS$pairs = factor(interaction(WuS$D1, WuS$D2, sep = ""), levels =interaction(WuS$D1, WuS$D2, sep = ""), ordered = TRUE)

plot(WuS$pairs, WuS$Mean, "p")

mean(subset(WuS$Mean, WuS$D1 >= 4 & WuS$D2 >= 4))
mean(subset(WuS$Mean, WuS$D1 < 8 & WuS$D2 <8 & WuS$D1 != 4 & WuS$D2 != 4))

#### similarity based on font pixel overlay like in Wong & Szuecs (2013) ####
library(readxl)
font_comp = read_xlsx("font_comp_table.xlsx")
font_comp_rel = read_xlsx("font_comp_table_rel.xlsx")

comp_split = strsplit(font_comp$Comparison, split = "", fixed = TRUE)
comp_split_rel = strsplit(font_comp$Comparison, split = "", fixed = TRUE)

font_comp$C1 = NA
font_comp$C2 = NA
for (i in 1:45) {
  font_comp$C1[i] = as.double(comp_split[[i]][1])
  font_comp$C2[i] = as.double(comp_split[[i]][2])
}

font_comp_rel$C1 = NA
font_comp_rel$C2 = NA
for (i in 1:45) {
  font_comp_rel$C1[i] = as.double(comp_split_rel[[i]][1])
  font_comp_rel$C2[i] = as.double(comp_split_rel[[i]][2])
}


comp_123 = colMeans(subset(font_comp[,2:163], font_comp$C1 %in% c(0,4,5,6,7,8,9) & font_comp$C2 %in% c(0,4,5,6,7,8,9)))
comp_456 = colMeans(subset(font_comp[,2:163], font_comp$C1 %in% c(1,2,3,4,5,6,7) & font_comp$C2 %in% c(1,2,3,4,5,6,7)))

comp_rel_123 = colMeans(subset(font_comp_rel[,2:163], font_comp_rel$C1 %in% c(0,4,5,6,7,8,9) & font_comp_rel$C2 %in% c(0,4,5,6,7,8,9)))
comp_rel_456 = colMeans(subset(font_comp_rel[,2:163], font_comp_rel$C1 %in% c(1,2,3,4,5,6,7) & font_comp_rel$C2 %in% c(1,2,3,4,5,6,7)))


mean_comp = data.frame(names(font_comp)[2:163], comp_123, comp_456, comp_rel_123, comp_rel_456)
names(mean_comp)[1] = "fonts"
rownames(mean_comp) = 1:162
