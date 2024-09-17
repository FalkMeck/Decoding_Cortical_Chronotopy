# ISC Making Plots
setwd(".../_Baysian_Models_control_z")

library(brms)
library(tidyr)
library(ggplot2)
library(ggdist)
library(writexl)
library(extrafont)
library(purrr)
library(ggh4x)

colCond = c("#39568C", "#440154","#FDE725", "#55C667")
colTextCond = c("black", "black", "white", "white")
colsCond = c("Single" = "#FDE725", "Triplet" = "#55C667", "Nonet" = "#39568C", "Complete" = "#440154")
strip = strip_themed(background_x = elem_list_rect(fill = colsCond),   text_x = elem_list_text(color = colTextCond))

# Rich CLub ####
load("brmISCz_CxRC.RData")

draws = as_draws_df(brmSNV_ISCz_CxRC, variable = "^b_", regex = TRUE)
names(draws)
# 30000 = 10 chains * 3000 (4000 - 1000 Warm up)
drawRC = data.frame(matrix(NA, length(draws$b_Intercept)), 0)
drawRC$Single_Local = draws$b_Intercept
drawRC$Single_Feeder = draws$b_Intercept + draws$b_RC_Class_SingRC_Feeder
drawRC$Single_RC = draws$b_Intercept + draws$b_RC_Class_SingRC_Club

drawRC$Triplet_Local = draws$b_Intercept + draws$b_ConditionTriplet
drawRC$Triplet_Feeder = draws$b_Intercept +  draws$b_ConditionTriplet + draws$b_RC_Class_SingRC_Feeder +  draws$`b_ConditionTriplet:RC_Class_SingRC_Feeder`
drawRC$Triplet_RC = draws$b_Intercept +  draws$b_ConditionTriplet + draws$b_RC_Class_SingRC_Club +   draws$`b_ConditionTriplet:RC_Class_SingRC_Club`

drawRC$Nonet_Local = draws$b_Intercept + draws$b_ConditionNonet
drawRC$Nonet_Feeder = draws$b_Intercept +  draws$b_ConditionNonet + draws$b_RC_Class_SingRC_Feeder +  draws$`b_ConditionNonet:RC_Class_SingRC_Feeder`
drawRC$Nonet_RC = draws$b_Intercept +  draws$b_ConditionNonet + draws$b_RC_Class_SingRC_Club +   draws$`b_ConditionNonet:RC_Class_SingRC_Club`

drawRC$Complete_Local = draws$b_Intercept + draws$b_ConditionComplete
drawRC$Complete_Feeder = draws$b_Intercept +  draws$b_ConditionComplete + draws$b_RC_Class_SingRC_Feeder +  draws$`b_ConditionComplete:RC_Class_SingRC_Feeder`
drawRC$Complete_RC = draws$b_Intercept +  draws$b_ConditionComplete + draws$b_RC_Class_SingRC_Club +   draws$`b_ConditionComplete:RC_Class_SingRC_Club`

drawRC = drawRC[,3:length(drawRC)]

drawLong = pivot_longer(drawRC, col = names(drawRC))
condNames = strsplit(drawLong$name, split = "_")
drawLong$Condition = unlist(map(condNames, 1))
drawLong$RC_Class = unlist(map(condNames, 2))
names(drawLong)[2] = c("ISC")

drawLong$RC_Class = factor(drawLong$RC_Class, levels = c("Local", "Feeder", "RC"), 
                           labels = c("Local node", "Feeder node", "RC hub"))
drawLong$Condition = factor(drawLong$Condition, levels = c("Single", "Triplet", "Nonet", "Complete"))

# #440154
colRC = c("#3d014c","#8e02b1","#d01cfd")
colsRC =c("Local node" = "#3d014c","Feeder node" = "#8e02b1","RC hub" = "#d01cfd")
RCplot= ggplot(drawLong, aes(y = ISC, x = RC_Class, fill = RC_Class)) +
  stat_eye(position = "dodge",
           aes(fill_ramp = after_stat(level)), .width = c(0.05, 0.5, 0.95)) +
  theme_bw() +   theme(legend.position = "none") + guides(fill = "none") + 
  labs(x = "Rich Club classification",
       y = "Inter-subject correaltion (Fisher z-transformed)",
       title = " ")+   facet_wrap2(vars(drawLong$Condition), ncol = 2, strip = strip) +
  theme_bw() + scale_fill_manual(name = "", values = colsRC) +
  theme(legend.position="none", strip.background=element_rect(fill=colsCond), strip.text = element_text(colour = "white"))
RCplot

save(RCplot, drawLong, brmSNV_ISCz_CxRC, colRC, colsRC, file = paste0(getwd(), "/PlotISCz_RC.RData"))
load(paste0(getwd(), "/PlotISCz_RC.RData"))

RCplot_paper = RCplot + theme(legend.position = "none") +
  theme(axis.text=element_text(size=11,  family="Calibri", color = "black")) +
  theme(text=element_text(size=11,  family="Calibri", color = "black")) 
 
ggsave("RCPlot_ISC_PaperPlot.png", plot = RCplot_paper,
       path = ".../Results_ISC", 
       width = 82.3*2, height = 50*2, units = "mm", dpi = 300)

# Diverse Club ####
load("brmISCz_CxDC.RData")

draws = as_draws_df(brmSNV_ISCz_CxDC, variable = "^b_", regex = TRUE)
names(draws)
# 30000 = 10 chains * 3000 (4000 - 1000 Warm up)
drawDC = data.frame(matrix(NA, length(draws$b_Intercept)), 0)
drawDC$Single_Local = draws$b_Intercept
drawDC$Single_Feeder = draws$b_Intercept + draws$b_DC_Class_SingDC_Feeder
drawDC$Single_DC = draws$b_Intercept + draws$b_DC_Class_SingDC_Club

drawDC$Triplet_Local = draws$b_Intercept + draws$b_ConditionTriplet
drawDC$Triplet_Feeder = draws$b_Intercept +  draws$b_ConditionTriplet + draws$b_DC_Class_SingDC_Feeder +  draws$`b_ConditionTriplet:DC_Class_SingDC_Feeder`
drawDC$Triplet_DC = draws$b_Intercept +  draws$b_ConditionTriplet + draws$b_DC_Class_SingDC_Club +   draws$`b_ConditionTriplet:DC_Class_SingDC_Club`

drawDC$Nonet_Local = draws$b_Intercept + draws$b_ConditionNonet
drawDC$Nonet_Feeder = draws$b_Intercept +  draws$b_ConditionNonet + draws$b_DC_Class_SingDC_Feeder +  draws$`b_ConditionNonet:DC_Class_SingDC_Feeder`
drawDC$Nonet_DC = draws$b_Intercept +  draws$b_ConditionNonet + draws$b_DC_Class_SingDC_Club +   draws$`b_ConditionNonet:DC_Class_SingDC_Club`

drawDC$Complete_Local = draws$b_Intercept + draws$b_ConditionComplete
drawDC$Complete_Feeder = draws$b_Intercept +  draws$b_ConditionComplete + draws$b_DC_Class_SingDC_Feeder +  draws$`b_ConditionComplete:DC_Class_SingDC_Feeder`
drawDC$Complete_DC = draws$b_Intercept +  draws$b_ConditionComplete + draws$b_DC_Class_SingDC_Club +   draws$`b_ConditionComplete:DC_Class_SingDC_Club`

drawDC = drawDC[,3:length(drawDC)]

drawLong = pivot_longer(drawDC, col = names(drawDC))
condNames = strsplit(drawLong$name, split = "_")
drawLong$Condition = unlist(map(condNames, 1))
drawLong$DC_Class = unlist(map(condNames, 2))
names(drawLong)[2] = c("ISC")

drawLong$DC_Class = factor(drawLong$DC_Class, levels = c("Local", "Feeder", "DC"), 
                           labels = c("Local node", "Feeder node", "DC hub"))
drawLong$Condition = factor(drawLong$Condition, levels = c("Single", "Triplet", "Nonet", "Complete"))

colDC = c("#212445","#49529c","#979cce")
colsDC = c("Local node" ="#212445" , "Feeder node" = "#49529c", "DC hub" = "#979cce")

DCplot= ggplot(drawLong, aes(y = ISC, x = DC_Class, fill = DC_Class)) +
  stat_eye(position = "dodge",
           aes(fill_ramp = after_stat(level)), .width = c(0.05, 0.5, 0.95)) +
  theme_bw() +   theme(legend.position = "none") + guides(fill = "none") + 
  labs(x = "Diverse Club classification",
       y = "Inter-subject correaltion (Fisher z-transformed)",
       title = " ")+   facet_wrap2(vars(drawLong$Condition), ncol = 2, strip = strip) +
  theme_bw() + scale_fill_manual(name = "", values = colsDC) +
  theme(legend.position="none", strip.background=element_rect(fill=colsCond), strip.text = element_text(colour = "white"))
DCplot

save(DCplot, drawLong, brmSNV_ISCz_CxDC, colDC, colsDC, file = paste0(getwd(), "/PlotISCz_TSRaut_DC.RData"))
load(paste0(getwd(), "/PlotISCz_TSRaut_DC.RData"))

DCplot_paper = DCplot + theme(legend.position = "none") +
  theme(axis.text=element_text(size=11,  family="Calibri", color = "black")) +
  theme(text=element_text(size=11,  family="Calibri", color = "black")) 

ggsave("DCPlot_ISC_PaperPlot.png", plot = DCplot_paper,
       path = ".../Results_ISC", 
       width = 82.3*2, height = 50*2, units = "mm", dpi = 300)

# Multimodality ####
load("brmISCz_CxSFC.RData")
load("brmISCz_CxJiIto.RData")

dataSFC = data.frame(id = brmSNV_ISCz_CxSFC$data$Subject, Condition = brmSNV_ISCz_CxSFC$data$Condition,
                     z_SFC_S7_Multi_seed = brmSNV_ISCz_CxSFC$data$z_SFC_S7_Multi_seed, z_Stats_Surf = 0, z_ICV_l = 0, z_Coord_z = 0)
dataJiIto = data.frame(id = brmSNV_ISCz_CxJiIto$data$Subject, Condition = brmSNV_ISCz_CxJiIto$data$Condition, 
                       z_Ji_Ito_multi_perc = brmSNV_ISCz_CxJiIto$data$z_Ji_Ito_multi_perc, z_Stats_Surf = 0, z_ICV_l = 0, z_Coord_z = 0)

dataSFC = cbind(dataSFC, fitted(brmSNV_ISCz_CxSFC, dataSFC, re_formula = NA))
names(dataSFC) = c("id", "Condition", "SFC","ICV", "Surf","z","estimate_SFC", "error_SFC", "lwr_SFC", "upr_SFC")
dataJiIto = cbind(dataJiIto, fitted(brmSNV_ISCz_CxJiIto, dataJiIto, re_formula = NA))
names(dataJiIto) = c("id", "Condition", "JI","ICV", "Surf","z","estimate_JI", "error_JI", "lwr_JI", "upr_JI")

dataMULT = data.frame(dataSFC$id, dataSFC$Condition, dataSFC$SFC, dataJiIto$JI,
                      dataSFC$estimate_SFC, dataSFC$lwr_SFC, dataSFC$upr_SFC,
                      dataJiIto$estimate_JI, dataJiIto$lwr_JI, dataJiIto$upr_JI)
names(dataMULT) = c("id", "Condition", "SFC", "JI",
                    "estimate_SFC", "lwr_SFC", "upr_SFC", 
                    "estimate_JI", "lwr_JI", "upr_JI")

colMulti = c("#1c5863", "#38afc7")
colsMulti = c("Stepwise functional connectivity"="#1c5863","Ji atlas"="#38afc7")

MULTplot <- ggplot() +
  # SFC
  geom_point(data = brmSNV_ISCz_CxSFC$data, aes(z_SFC_S7_Multi_seed, ISC_FishZ),
             alpha = .2, color = colMulti[1], size = .05) +
  # Ji atlas
  geom_point(data = brmSNV_ISCz_CxJiIto$data, aes(z_Ji_Ito_multi_perc, ISC_FishZ),
             alpha = .2, color = colMulti[2], size = .05) +
  # SFC
  geom_ribbon(data = dataMULT, aes(x = SFC, y = estimate_SFC, ymin = lwr_SFC, ymax = upr_SFC), alpha = .3, fill = colMulti[1]) +
  geom_line(data = dataMULT, aes(SFC, estimate_SFC, color =  "Stepwise functional connectivity"), size = .5) +
  # Ji atlas
  geom_ribbon(data = dataMULT, aes(x = JI, y = estimate_JI, ymin = lwr_JI, ymax = upr_JI), alpha = .3, fill = colMulti[2]) +
  geom_line(data = dataMULT, aes(JI, estimate_JI, color =  "Ji atlas"), size = .5) +
  
  facet_wrap2(dataMULT$Condition, ncol = 2, strip = strip) +
  labs(x ="Multimodality (z-standardized)",
       y = "Inter-subject correaltion (Fisher z-transformed)",
       title = " ")+ 
  theme_bw() + scale_color_manual(name = "", values = colsMulti, breaks = labels(colsMulti))
MULTplot

save(MULTplot, dataSFC, dataJiIto, dataMULT, colsMulti, colMulti, brmSNV_ISCz_CxJiIto, brmSNV_ISCz_CxSFC, file = paste0(getwd(), "/PlotISCz_MULT.RData"))
load(paste0(getwd(), "/PlotISCz_MULT.RData"))

MULTplot_paper = MULTplot + theme(legend.position = "top") +
  theme(axis.text=element_text(size=11,  family="Calibri", color = "black")) +
  theme(text=element_text(size=11,  family="Calibri", color = "black")) 


ggsave("MULTplotz_ISCz_PaperPlot.png", plot = MULTplot_paper,
       path =".../Results_ISC", 
       width = 82.3*2, height = 50*2, units = "mm", dpi = 300)

# Coordiantes ####

load("brmISCz_CxXY.RData")
model = brmSNV_ISCz_CxXY

# https://kzee.github.io/PlotFixef_Demo.html
dataX = data.frame(id = model$data$Subject, Condition = model$data$Condition,
                          z_absCoord_x = model$data$z_absCoord_x, z_Coord_y = 0, z_Coord_z = 0, 
                          z_Stats_Surf = 0, z_ICV_l = 0)
dataY = data.frame(id = model$data$Subject, Condition = model$data$Condition,
                   z_absCoord_x = 0, z_Coord_y = model$data$z_Coord_y, z_Coord_z = 0, 
                   z_Stats_Surf = 0, z_ICV_l = 0)
dataXY = data.frame(id = model$data$Subject, Condition = model$data$Condition, z_absCoord_x = model$data$z_absCoord_x, 
                     z_Coord_y = model$data$z_Coord_y, z_ICV_l = 0, z_Stats_Surf = 0, z_Coord_z = 0)

dataX = cbind(dataX, fitted(model, dataX, re_formula = NA))
names(dataX) = c("id", "Condition", "x", "y", "z", "ICV", "Surf","estimate_x", "error_x", "lwr_x", "upr_x")
dataY = cbind(dataY, fitted(model, dataY, re_formula = NA))
names(dataY) = c("id", "Condition", "x", "y", "z", "ICV", "Surf", "estimate_y", "error_y", "lwr_y", "upr_y")
dataXY = cbind(dataXY, fitted(model, dataXY, re_formula = NA))
names(dataXY) = c("id", "x", "y", "ICV", "Surf", "z", "estimate", "error", "lwr", "upr")

dataXY = data.frame(dataX$id, dataX$Condition, dataX$x, dataY$y, # dataZ$z, 
                dataX$estimate_x, dataX$lwr_x, dataX$upr_x,
                dataY$estimate_y, dataY$lwr_y, dataY$upr_y,
                dataXY$estimate, dataXY$lwr, dataXY$upr) # , 
names(dataXY) = c("id", "Condition", "x", "y",
                   "estimate_x", "lwr_x", "upr_x", 
                   "estimate_y", "lwr_y", "upr_y", 
                "estiamte_XY", "lwr_XY", "upr_XY")

colXY = c("#13533c", "#5ad8ac")
colsXY = c("|x|"="#13533c","y"="#5ad8ac")

XYplot <- ggplot() +
  # X
  geom_point(data = model$data, aes(z_absCoord_x, ISC_FishZ),
             alpha = .2, color = colXY[1], size = .05) +
  # Y
  geom_point(data = model$data, aes(z_Coord_y, ISC_FishZ),
             alpha = .2, color = colXY[2], size = .05) +
  # X
  geom_ribbon(data = dataXY, aes(x = x, y = estimate_x, ymin = lwr_x, ymax = upr_x), alpha = .3, fill = colXY[1]) +
  geom_line(data = dataXY, aes(x, estimate_x, color =  "|x|"), size = .5) +
  # Y
  geom_ribbon(data = dataXY, aes(x = y, y = estimate_y, ymin = lwr_y, ymax = upr_y), alpha = .3, fill = colXY[2]) +
  geom_line(data = dataXY, aes(y, estimate_y, color =  "y"), size = .5) +

  facet_wrap2(model$data$Condition, ncol = 2, strip = strip) +
  labs(x = "Coordinate (z-standardized)",
       y = "Inter-subject correaltion (Fisher z-transformed)",
       title = " ")+ 
  theme_bw() + scale_color_manual(name = "", values = colsXY)
XYplot

save(XYplot, dataX, dataY,dataXY ,
     colsXY, colXY, model, file = paste0(getwd(), "/PlotISCz_aXxY.RData"))
load(paste0(getwd(), "/PlotISCz_aXxY.RData"))

XYplot_paper = XYplot + theme(legend.position = "top") +
  theme(axis.text=element_text(size=11,  family="Calibri", color = "black")) +
  theme(text=element_text(size=11,  family="Calibri", color = "black")) 

ggsave("XYPlot_ISCz_PaperPlot.png", plot = XYplot_paper,
       path = ".../Results_ISC", 
       width = 82.3*2, height = 50*2, units = "mm", dpi = 300)


# Graph Measures # IS NOT GREAT YET... ####
load("brmISCz_CxGM.RData")
model = brmSNV_ISCz_CxGM

# https://kzee.github.io/PlotFixef_Demo.html
dataDeg = data.frame(id = model$data$Subject, Condition = model$data$Condition,
                     z_Degree = model$data$z_Degree, z_Closeness = 0, z_Betweenness = 0, 
                     z_PartCoef_DC = 0, z_WithinMod_DC = 0,
                   z_Stats_Surf = 0, z_ICV_l = 0, z_Coord_z = 0)
dataCls = data.frame(id = model$data$Subject, Condition = model$data$Condition,
                     z_Degree = 0, z_Closeness =  model$data$z_Closeness, z_Betweenness = 0, 
                     z_PartCoef_DC = 0, z_WithinMod_DC = 0,
                     z_Stats_Surf = 0, z_ICV_l = 0, z_Coord_z = 0)
dataBtw = data.frame(id = model$data$Subject, Condition = model$data$Condition,
                     z_Degree = 0, z_Closeness = 0, z_Betweenness =  model$data$z_Betweenness, 
                     z_PartCoef_DC = 0, z_WithinMod_DC = 0,
                     z_Stats_Surf = 0, z_ICV_l = 0, z_Coord_z = 0)
dataPtC = data.frame(id = model$data$Subject, Condition = model$data$Condition,
                     z_Degree =0, z_Closeness = 0, z_Betweenness = 0, 
                     z_PartCoef_DC =  model$data$z_PartCoef_DC, z_WithinMod_DC = 0,
                     z_Stats_Surf = 0, z_ICV_l = 0, z_Coord_z = 0)
dataWMz = data.frame(id = model$data$Subject, Condition = model$data$Condition,
                     z_Degree = 0, z_Closeness = 0, z_Betweenness = 0, 
                     z_PartCoef_DC = 0, z_WithinMod_DC = model$data$z_WithinMod_DC,
                     z_Stats_Surf = 0, z_ICV_l = 0, z_Coord_z = 0)

dataDeg = cbind(dataDeg, fitted(model, dataDeg, re_formula = NA))
names(dataDeg) = c("id", "Condition", "Deg", "Cls", "Btw", "PtC", "WMz", "ICV", "Surf","z","estimate_Deg", "error_Deg", "lwr_Deg", "upr_Deg")
dataCls = cbind(dataCls, fitted(model, dataCls, re_formula = NA))
names(dataCls) = c("id", "Condition", "Deg", "Cls", "Btw", "PtC", "WMz", "ICV", "Surf","z","estimate_Cls", "error_Cls", "lwr_Cls", "upr_Cls")
dataBtw = cbind(dataBtw, fitted(model, dataBtw, re_formula = NA))
names(dataBtw) = c("id", "Condition", "Deg", "Cls", "Btw", "PtC", "WMz", "ICV", "Surf","z","estimate_Btw", "error_Btw", "lwr_Btw", "upr_Btw")
dataPtC = cbind(dataPtC, fitted(model, dataPtC, re_formula = NA))
names(dataPtC) = c("id", "Condition", "Deg", "Cls", "Btw", "PtC", "WMz", "ICV", "Surf","z","estimate_PtC", "error_PtC", "lwr_PtC", "upr_PtC")
dataWMz = cbind(dataWMz, fitted(model, dataWMz, re_formula = NA))
names(dataWMz) = c("id", "Condition", "Deg", "Cls", "Btw", "PtC", "WMz", "ICV", "Surf","z","estimate_WMz", "error_WMz", "lwr_WMz", "upr_WMz")

dataGM = data.frame(dataDeg$id, dataDeg$Condition, dataDeg$Deg, dataCls$Cls, dataBtw$Btw,dataPtC$PtC, dataWMz$WMz, 
                    dataDeg$estimate_Deg, dataDeg$lwr_Deg, dataDeg$upr_Deg,
                    dataCls$estimate_Cls, dataCls$lwr_Cls, dataCls$upr_Cls, 
                    dataBtw$estimate_Btw, dataBtw$lwr_Btw, dataBtw$upr_Btw, 
                    dataPtC$estimate_PtC, dataPtC$lwr_PtC, dataPtC$upr_PtC, 
                    dataWMz$estimate_WMz, dataWMz$lwr_WMz, dataWMz$upr_WMz)
names(dataGM) = c("id", "Condition", "Deg", "Cls", "Btw", "PtC", "WMz",
                   "estimate_Deg", "lwr_Deg", "upr_Deg", 
                   "estimate_Cls", "lwr_Cls", "upr_Cls", 
                   "estimate_Btw", "lwr_Btw", "upr_Btw",
                   "estimate_PtC", "lwr_PtC", "upr_PtC",
                   "estimate_WMz", "lwr_WMz", "upr_WMz")

# #95d840
colGM = c("#456a16","#6fa923","#95d840", "#b9e580", "#dcf2c0")
colsGM = c("Degree centrality" ="#456a16", "Closeness centrality" = "#6fa923",
           "Betweenness centrality" = "#95d840", 
           "Participation coefficient"="#b9e580", "Within-module degree z-score"="#dcf2c0")

GMplot <- ggplot() +
  # Deg
  geom_point(data = model$data, aes(z_Degree, ISC_FishZ),
             alpha = .2, color = colGM[1], size = .05) +
  # Cls
  geom_point(data = model$data, aes(z_Closeness, ISC_FishZ),
             alpha = .2, color = colGM[2], size = .05) +
  # Btw
  geom_point(data = model$data, aes(z_Betweenness, ISC_FishZ),
             alpha = .2, color = colGM[3], size = .05) +
  # PtC
  geom_point(data = model$data, aes(z_PartCoef_DC, ISC_FishZ),
             alpha = .2, color = colGM[4], size = .05) +
  # WMz
  geom_point(data = model$data, aes(z_WithinMod_DC, ISC_FishZ),
             alpha = .2, color = colGM[5], size = .05) +
  # Deg
  geom_ribbon(data = dataGM, aes(x = Deg, y = estimate_Deg, ymin = lwr_Deg, ymax = upr_Deg), alpha = .3, fill = colGM[1]) +
  geom_line(data = dataGM, aes(Deg, estimate_Deg, color =  "Degree centrality"), size = .5) +
  # Cls
  geom_ribbon(data = dataGM, aes(x = Cls, y = estimate_Cls, ymin = lwr_Cls, ymax = upr_Cls), alpha = .3, fill = colGM[2]) +
  geom_line(data = dataGM, aes(Cls, estimate_Cls, color =  "Closeness centrality"), size = .5) +
  # Btw
  geom_ribbon(data = dataGM, aes(x = Btw, y = estimate_Btw, ymin = lwr_Btw, ymax = upr_Btw), alpha = .3, fill = colGM[3]) +
  geom_line(data = dataGM, aes(Btw, estimate_Btw, color =  "Betweenness centrality"), size = .5) +
  # PtC
  geom_ribbon(data = dataGM, aes(x = PtC, y = estimate_PtC, ymin = lwr_PtC, ymax = upr_PtC), alpha = .3, fill = colGM[4]) +
  geom_line(data = dataGM, aes(PtC, estimate_PtC, color =  "Participation coefficient"), size = .5) +
  # WMz
  geom_ribbon(data = dataGM, aes(x = WMz, y = estimate_WMz, ymin = lwr_WMz, ymax = upr_WMz), alpha = .3, fill = colGM[5]) +
  geom_line(data = dataGM, aes(WMz, estimate_WMz, color =  "Within-module degree z-score"), size = .5) +
  
  facet_wrap2(dataGM$Condition, ncol = 2, strip = strip) +
  labs(x = "Centrality measure (z-standardized)",
       y =  "Inter-subject correaltion (Fisher z-transformed)",
       title = " ")+ 
  theme_bw() + scale_color_manual(name = "", values = colsGM, breaks = labels(colsGM))
GMplot

save(GMplot, dataDeg, dataCls, dataBtw, dataPtC, dataWMz, dataGM, colsGM, colGM, model, file = paste0(getwd(), "/PlotISCz_GM.RData"))
load(paste0(getwd(), "/PlotISCz_GM.RData"))

GMplot_paper = GMplot + theme(legend.position = "top", legend.text=element_text(size= 5.5)) +
  theme(axis.text=element_text(size=11,  family="Calibri", color = "black")) +
  theme(text=element_text(size=11,  family="Calibri", color = "black")) 

ggsave("GMPlot_ISCz_PaperPlot.png", plot = GMplot_paper,
       path = ".../Results_ISC", 
       width = 82.3*2, height = 50*2, units = "mm", dpi = 300)

# Cortical Thickness ####
load("brmISCz_CxSBMt.RData")
model = brmSNV_ISCz_CxSBMt

# https://kzee.github.io/PlotFixef_Demo.html
dataSBMt = data.frame(id = model$data$Subject, Condition = model$data$Condition,
                      z_Stats_Thick = model$data$z_Stats_Thick, z_Stats_Surf = 0, z_ICV_l = 0, z_Coord_z = 0)

dataSBMt = cbind(dataSBMt, fitted(model, dataSBMt, re_formula = NA))
names(dataSBMt) = c("id", "Condition", "SBMt", "ICV", "Surf","z","estimate_SBMt", "error_SBMt", "lwr_SBMt", "upr_SBMt")

colThick = "#fce303"
colsThick = c("Cortical thickness" = "#fce303")

Thickplot <- ggplot() +
  geom_point(data = model$data, aes(z_Stats_Thick, ISC_FishZ),
             alpha = .2, color = colThick, size = .05) +
  geom_ribbon(data = dataSBMt, aes(x = SBMt, y = estimate_SBMt, ymin = lwr_SBMt, ymax = upr_SBMt), alpha = .3, fill = colThick) +
  geom_line(data = dataSBMt, aes(SBMt, estimate_SBMt, color =  "Cortical thickness"), size = .5) +
 
  facet_wrap2(dataSBMt$Condition, ncol = 2, strip = strip) +
  labs(x = "Cortical thickness (z-standardized)",
       y = "Inter-subject correaltion (Fisher z-transformed)",
       title = " ")+ 
  theme_bw() + scale_color_manual(name = "", values = colsThick)
Thickplot

save(Thickplot, dataSBMt, colThick, colsThick, model, file = paste0(getwd(), "/PlotISCz_SBMt.RData"))
load(paste0(getwd(), "/PlotISCz_SBMt.RData"))

Thickplot_paper = Thickplot + theme(legend.position = "top") +
  theme(axis.text=element_text(size=11,  family="Calibri", color = "black")) +
  theme(text=element_text(size=11,  family="Calibri", color = "black")) 


ggsave("CT_Plot_ISCz_PaperPlot.png", plot = Thickplot_paper,
       path = ".../Results_ISC", 
       width = 82.3*2, height = 50*2, units = "mm", dpi = 300)
