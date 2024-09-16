# Matrix plot spin test

setwd("...")

library(corrplot)

varNames = c("Resting-state timescales", "RC", "DC","SFC", "Ji atlas", "|x|-coordinate", "y-coordinate", 
             "Degree centrality", "Betweenness centrality ", "Closeness centrality", 
             "Participation coefficient", "Within-module degree z-score", "Cortical thickness", 
             "ISC (Single)","ISC (Triplet)", "ISC (Nonet)", "ISC (Complete)")

# load matrices
effectSizes = as.matrix(read.delim("EffectSizes_MeanVariables.txt", sep = ',', header = FALSE))
rownames(effectSizes) = colnames(effectSizes) = varNames
rawValues = as.matrix(read.delim("rawValues_MeanVariables.txt", sep = ',', header = FALSE))
rownames(rawValues) = colnames(rawValues) = varNames
spinP_FDR = as.matrix(read.delim("spin_pFDRtrue_MeanVariables.txt", sep = ',', header = FALSE))
rownames(spinP_FDR) = colnames(spinP_FDR) = varNames
spinP_sig = as.matrix(read.delim("spin_pFDRtrue_sig_MeanVariables.txt", sep = ',', header = FALSE))
rownames(spinP_sig) = colnames(spinP_sig) = varNames
spinP_FDRlow = as.matrix(read.delim("spin_pFDRlower_MeanVariables.txt", sep = ',', header = FALSE))
rownames(spinP_FDRlow) = colnames(spinP_FDRlow) = varNames
spinP_siglow = as.matrix(read.delim("spin_pFDRlower_sig_MeanVariables.txt", sep = ',', header = FALSE))
rownames(spinP_siglow) = colnames(spinP_siglow) = varNames
spinP_FDRup = as.matrix(read.delim("spin_pFDRupper_MeanVariables.txt", sep = ',', header = FALSE))
rownames(spinP_FDRup) = colnames(spinP_FDRup) = varNames
spinP_sigup = as.matrix(read.delim("spin_pFDRupper_sig_MeanVariables.txt", sep = ',', header = FALSE))
rownames(spinP_sigup) = colnames(spinP_sigup) = varNames

trace(corrplot, edit=TRUE)
# EDDITED 445: 
'place_points = function(sig.locs, point) {
text(pos.pNew[, 1][sig.locs], pos.pNew[, 2][sig.locs],
     labels = point, col = pch.col, cex = pch.cex,
     lwd = 2)'
# TO
'place_points = function(sig.locs, point) {
text((pos.pNew[, 1][sig.locs])+0.4, (pos.pNew[, 2][sig.locs]) + 0.3,
     labels = point, col = pch.col, cex = pch.cex,
     lwd = 2)'
# AND SAVED --> moves stars a little up an a little to the right
library(scico)
scico_palette_data("imola",categorical = TRUE)
library(scales)
viridisHexa = viridis_pal()(20)
colorPaletteChoice = colorRampPalette(viridisHexa)
spinPdisp = spinP_FDR
spinPdisp[spinP_FDR >= 0.05] = NA
rawDisp = rawValues
rawDisp[spinP_FDR >= 0.05] = NA


# reorder
effectSizes_cut = effectSizes[c(1:6,8:18),c(1:6,8:18)]
spinP_FDR_cut = spinP_FDR[c(1:6,8:18),c(1:6,8:18)]
corrplot(effectSizes_cut, method = 'color', type = 'upper', diag = FALSE, 
         p.mat = spinP_FDR_cut,
         insig = "label_sig",
         pch.cex = 1.2,
         pch.col = "black",
         tl.col="black",
         tl.cex=0.8, tl.srt = 45,
         addCoef.col = "black",
         tl.pos="td", number.font = 2,
         outline=FALSE,
         col=colorRampPalette(viridisHexa)(200))
spinPdisp_cut = spinP_FDR_cut
spinPdisp_cut[spinPdisp_cut >= 0.05] = 1
