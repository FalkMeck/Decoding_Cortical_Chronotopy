library(brms)
library(tidyr)
library(ggplot2)
library(ggdist)
library(scico)
library(writexl)

scico_palette_show()


analysis_dir = ".../"
outDir = paste0(analysis_dir, "_Bayesian_Model_comp_final")
models = intersect(list.files(outDir, pattern = ".RData"), list.files(outDir, pattern = "^TS"))
modelNames = sapply(strsplit(models, ".", fixed = TRUE), "[[",1)
setwd(outDir)

modelsRes = c("TSRaut_aXxY", "TSRaut_DC", "TSRaut_JiItoPerc", "TSRaut_RC", "TSRaut_SFCindi")
modelResExplore = c("TSRaut_GM",
 "TSRaut_GMDeg","TSRaut_GMBtw","TSRaut_GMCls","TSRaut_GMPtC","TSRaut_GMWMz",
 "TSRaut_SBMt")

load(paste0(analysis_dir, "TS_Raut_Data_final_lausanne250_2024_09_05.RData"))

# RESULTS ####
# RC ####
lazyLoad("TSRaut_RC")
model = brmodel_TSRaut_RC
summary(model)
h <- c("RC_Class_SingRC_Feeder > 0",
       "RC_Class_SingRC_Club  > 0",
       "RC_Class_SingRC_Club - RC_Class_SingRC_Feeder > 0")
hRC = brms::hypothesis(model, h)
write_xlsx(hRC$hypothesis, path = paste0(outDir, "/TSRaut_RC/TSRaut_RC_h.xlsx"))

# DC ####
lazyLoad("TSRaut_DC")
model = brmodel_TSRaut_DC
summary(brmodel_TSRaut_DC)

h <- c("DC_Class_SingDC_Feeder > 0",
       "DC_Class_SingDC_Club  > 0",
       "DC_Class_SingDC_Club - DC_Class_SingDC_Feeder > 0")
hDC  = brms::hypothesis(model, h)
write_xlsx(hDC$hypothesis, path = paste0(outDir, "/TSRaut_DC/TSRaut_DC_h.xlsx"))

# Ji Ito ####
yellowMulti = c("#cf9a30", "#7b5c1e") 
lazyLoad("TSRaut_JiItoPerc")
model = brmodel_TSRaut_JiItoPerc
summary(model)
h <- c( "z_Ji_Ito_multi_perc > 0")
hJi_Ito = brms::hypothesis(model, h)
write_xlsx(hJi_Ito$hypothesis,  path = paste0(outDir, "/TSRaut_JiItoPerc/TSRaut_JiItoPerc_h.xlsx"))

# SFC indi ####
lazyLoad("TSRaut_SFCindi")
model = brmodel_TSRaut_SFCindi
summary(model)
h <- c( "z_SFC_S7_Multi_seed > 0")
hSFC = brms::hypothesis(model, h)
write_xlsx(hSFC$hypothesis,  path = paste0(outDir, "/TSRaut_SFCindi/TSRaut_SFCindi_h.xlsx"))
 
# XYZ ####
lazyLoad("TSRaut_aXxY")
model = brmodel_TSRaut_aXxY
summary(model)
h <- c( "z_absCoord_x <  0",
         "z_Coord_y > 0",
        "z_absCoord_x:z_Coord_y < 0")
hXY = brms::hypothesis(model, h)
write_xlsx(hXY$hypothesis, path = paste0(outDir, "/TSRaut_aXxY/TSRaut_aXxY_h.xlsx"))


# Graphmeasures all ####
lazyLoad("TSRaut_GM")
model = brmodel_TSRaut_GM
summary(model)
h <- c( "z_Degree >  0",
        "z_Closeness > 0",
        "z_Betweenness > 0", 
        "z_PartCoef_DC > 0",
        "z_WithinMod_DC > 0")
hGM = brms::hypothesis(model, h)
write_xlsx(hGM$hypothesis, path = paste0(outDir, "/TSRaut_GM/TSRaut_GM_h.xlsx"))

# Surface based morphometry thickness ####
lazyLoad("TSRaut_SBMt")
model = brmodel_TSRaut_SBMt
summary(model)
h <- c( "z_Stats_ThickAvg > 0")
hSBMt = brms::hypothesis(model, h)
write_xlsx(hSBMt$hypothesis,  path = paste0(outDir, "/TSRaut_SBMt/TSRaut_SBMt_h.xlsx"))
#
