### TS RESULTS (wihtout z control) 

# LIBRARIES
library(brms)
library(tidyr)
library(ggplot2)
library(ggdist)
library(scico)
library(writexl)

library(extrafont)

# Directories
analysis_dir = ".../"
outDir = paste0(analysis_dir, "_Bayesian_Model_comp_finalZ")
models = intersect(list.files(outDir, pattern = ".RData"), list.files(outDir, pattern = "^TS"))
modelNames = sapply(strsplit(models, ".", fixed = TRUE), "[[",1)
setwd(outDir)

load(paste0(analysis_dir, "TS_Raut_Data_final_lausanne250_2024_09_05.RData"))

# RESULTS ####
# RC ####
lazyLoad("TSRaut_RC")
model = brmodel_TSRaut_RC
summary(model)

draws = as_draws_df(model, variable = "^b_", regex = TRUE)
# 30000 = 10 chains * 3000 (4000 - 1000 Warm up)
drawRC = data.frame(matrix(NA, length(draws$b_Intercept)), 0)
drawRC$Local = draws$b_Intercept
drawRC$Feeder = draws$b_Intercept + draws$b_RC_Class_SingRC_Feeder
drawRC$RC = draws$b_Intercept + draws$b_RC_Class_SingRC_Club
drawRC = drawRC[,3:length(drawRC)]

drawLong = pivot_longer(drawRC, col = names(drawRC))
names(drawLong) = c("RC_Class", "z_logTSRaut")

drawLong$RC_Class = factor(drawLong$RC_Class, levels = c("Local", "Feeder", "RC"), 
                           labels = c("Local node", "Feeder node", "RC hub"))

# colors
# #440154
colRC = c("#3d014c","#8e02b1","#d01cfd")
colsRC =c("Local node" = "#3d014c","Feeder node" = "#8e02b1","RC hub" = "#d01cfd") 
RCplot= ggplot(drawLong, aes(y = z_logTSRaut, x = RC_Class, fill = RC_Class)) +
 stat_eye(position = "dodge",
         aes(fill_ramp = after_stat(level)), .width = c(0.05, 0.5, 0.95)) +
  theme_bw() +   theme(legend.position = "none") + guides(fill = "none") + 
  labs(x = "Rich Club classification",
       y = "Resting-state timescale estimate  \n (logarithmized & z-standardized)",
       title = " ")+ 
  theme_bw() + scale_fill_manual(name = "", values = colsRC)

RCplot
save(RCplot, drawLong,model, colRC, colsRC, file = paste0(outDir, "/TSRaut_RC/Plot_TSRaut_RC.RData"))
rm(brmodel_TSRaut_RC)

RCplot_paper = RCplot + theme(legend.position = "none") +
  theme(axis.text=element_text(size=11,  family="Calibri", color = "black")) +
  theme(text=element_text(size=11,  family="Calibri", color = "black")) 

ggsave(paste0("RCplot_","PaperPlot.png"), plot = RCplot_paper,
       path = ".../Results_TS", 
       width = 82.3*2, height = 50*2, units = "mm", dpi = 300)

# DC ####
lazyLoad("TSRaut_DC")
model = brmodel_TSRaut_DC
summary(brmodel_TSRaut_DC)

draws = as_draws_df(model, variable = "^b_", regex = TRUE)
# 30000 = 10 chains * 3000 (4000 - 1000 Warm up)
drawDC = data.frame(matrix(NA, length(draws$b_Intercept)), 0)
drawDC$Local = draws$b_Intercept
drawDC$Feeder = draws$b_Intercept + draws$b_DC_Class_SingDC_Feeder
drawDC$DC = draws$b_Intercept + draws$b_DC_Class_SingDC_Club
drawDC = drawDC[,3:length(drawDC)]

drawLong = pivot_longer(drawDC, col = names(drawDC))
names(drawLong) = c("DC_Class", "z_logTSRaut")

drawLong$DC_Class = factor(drawLong$DC_Class, levels = c("Local", "Feeder", "DC"), 
                           labels = c("Local node", "Feeder node", "DC hub"))

#404788
colDC = c("#212445","#49529c","#979cce")
colsDC = c("Local node" ="#212445" , "Feeder node" = "#49529c", "DC hub" = "#979cce")

DCplot= ggplot(drawLong, aes(y = z_logTSRaut, x = DC_Class, fill = DC_Class)) +
  stat_eye(position = "dodge",
           aes(fill_ramp = after_stat(level)), .width = c(0.05, 0.5, 0.95)) +
  theme_bw() +   theme(legend.position = "none") + guides(fill = "none") + 
  labs(x = "Diverse Club classification",
       y = "Resting-state timescale estimate  \n (logarithmized & z-standardized)",
       title = " ")+ 
  theme_bw() + scale_fill_manual(name = "", values = colsDC)
DCplot
save(DCplot, drawLong, model, colDC, colsDC, file = paste0(outDir, "/TSRaut_DC/Plot_TSRaut_DC.RData"))


DCplot_paper = DCplot + theme(legend.position = "none") +
  theme(axis.text=element_text(size=11,  family="Calibri", color = "black")) +
  theme(text=element_text(size=11,  family="Calibri", color = "black")) 

ggsave(paste0("DCPlot_","PaperPlot_20240610.png"), plot = DCplot_paper,
       path = ".../Results_TS", 
       width = 82.3*2, height = 50*2, units = "mm", dpi = 300)

# Multi Plot #####
load( paste0(outDir, "/TSRaut_SFCindi/Plot_TSRaut_MULTplot.RData"))

colMulti = c("#1c5863", "#38afc7")
colsMulti = c("Step-wise functional connectivity"="#1c5863","Ji atlas"="#38afc7")


MULTplot = ggplot () + 
  geom_point(data = brmodel_TSRaut_SFCindi$data, aes(z_SFC_S7_Multi_seed, z_logTS_Raut),
             alpha = .2, color = colMulti[1], size = .5) +
  geom_point(data = brmodel_TSRaut_JiItoPerc$data, aes(z_Ji_Ito_multi_perc, z_logTS_Raut),
             alpha = .2, color = colMulti[2], size = .5) +
  geom_ribbon(data = dataSFC, aes(x = SFC, y = estimate, ymin = lwr, ymax = upr), alpha = .3, fill = colMulti[1]) +
  geom_line(data = dataSFC, aes(SFC, estimate, color = "Step-wise functional connectivity"), size = 1) +
  geom_ribbon(data = dataJiIto, aes(x = JiIto, y = estimate, ymin = lwr, ymax = upr), alpha = .3,fill = colMulti[2]) +
  geom_line(data = dataJiIto, aes(JiIto, estimate, color = "Ji atlas"), size = 1) +
  labs(x = "Multimodality (z-standardized)",
       y = "Resting-state timescale estimate  \n (logarithmized & z-standardized)",
       title = " ") +
  theme_bw() + scale_color_manual(name = "", values = colsMulti, breaks = labels(colsMulti))
MULTplot
save(MULTplot, colMulti, colsMulti, brmodel_TSRaut_JiItoPerc, brmodel_TSRaut_SFCindi, dataJiIto, dataSFC,
     file = paste0(outDir, "/TSRaut_SFCindi/Plot_TSRaut_MULTplot.RData"))
rm(brmodel_TSRaut_SFCindi,brmodel_TSRaut_JiItoPerc)


MULTplot_paper = MULTplot + theme(legend.position = c (0.8,0.15)) +
  theme(axis.text=element_text(size=11,  family="Calibri", color = "black")) +
  theme(text=element_text(size=11,  family="Calibri", color = "black")) 

ggsave(paste0("MultiPlot_","PaperPlot.png"), plot = MULTplot_paper,
       path = "...Graphics/Results_TS/", 
       width = 82.3*2, height = 50*2, units = "mm", dpi = 300)

# XY ####
load(paste0(outDir, "/TSRaut_aXxY/Plot_TSRaut_aXxY.RData"))

colXY = c("#13533c", "#5ad8ac")
colsXY = c("|x|"="#13533c","y"="#5ad8ac")

XYplot <- ggplot() +
  # X
  geom_point(data = model$data, aes(z_absCoord_x, z_logTS_Raut),
    alpha = .2, color = colXY[1], size = .5) +
  # Y
  geom_point(data = model$data, aes(z_Coord_y, z_logTS_Raut),
             alpha = .2, color = colXY[2], size = .5) +
  # X
  geom_ribbon(data = dataX, aes(x = x, y = estimate, ymin = lwr, ymax = upr), alpha = .3, fill = colXY[1]) +
  geom_line(data = dataX, aes(x, estimate, color =  "|x|"), size = 1) +
  # Y
  geom_ribbon(data = dataY, aes(x = y, y = estimate, ymin = lwr, ymax = upr), alpha = .3, fill = colXY[2]) +
  geom_line(data = dataY, aes(y, estimate, color =  "y"), size = 1) +

    labs(x = "Coordinate (z-standardized)",
       y = "Resting-state timescale estimate  \n (logarithmized & z-standardized)",
       title = " ")+ 
  theme_bw() + scale_color_manual(name = "", values = colsXY)
XYplot

save(XYplot, dataX, dataY, dataXY, colsXY, colXY, model, file = paste0(outDir, "/TSRaut_aXxY/Plot_TSRaut_aXxY.RData"))

XYplot_paper = XYplot + theme(legend.position = c (0.93,0.2)) +
  theme(axis.text=element_text(size=11,  family="Calibri", color = "black")) +
  theme(text=element_text(size=11,  family="Calibri", color = "black")) 


ggsave(paste0("XYPlot_","PaperPlot.png"), plot = XYplot_paper,
       path = ".../Results_TS", 
       width = 82.3*2, height = 50*2, units = "mm", dpi = 300)


# Graphmeasures all ####
load(paste0(outDir, "/TSRaut_GM/Plot_TSRaut_GM.RData"))

# #95d840
colGM = c("#456a16","#6fa923","#95d840", "#b9e580", "#dcf2c0")
colsGM = c("Degree centrality" ="#456a16", "Closeness centrality" = "#6fa923",
           "Betweenness centrality" = "#95d840", 
           "Participation coefficient"="#b9e580", "Within-module degree z-score"="#dcf2c0")

GMplot <- ggplot() +
  # Deg
  geom_point(data = model$data, aes(z_Degree, z_logTS_Raut),
             alpha = .2, color = colGM[1], size = .5) +
  # Cls
  geom_point(data = model$data, aes(z_Closeness, z_logTS_Raut),
             alpha = .2, color = colGM[2], size = .5) +
  # Btw
  geom_point(data = model$data, aes(z_Betweenness, z_logTS_Raut), 
             alpha = .2, color = colGM[3], size = .5) +
  # PtC
  geom_point(data = model$data, aes(z_PartCoef_DC, z_logTS_Raut),
             alpha = .2, color = colGM[4], size = .5) +
  #WMz
  geom_point(data = model$data, aes(z_WithinMod_DC, z_logTS_Raut),
             alpha = .2, color = colGM[5], size = .5) +
  # Deg
  geom_ribbon(data = dataDeg, aes(x = Deg, y = estimate, ymin = lwr, ymax = upr), alpha = .3, fill = colGM[1]) +
  geom_line(data = dataDeg, aes(Deg, estimate, color = "Degree centrality"), size = 1) +
  # Cls
  geom_ribbon(data = dataCls, aes(x = Cls, y = estimate, ymin = lwr, ymax = upr), alpha = .3, fill= colGM[2]) +
  geom_line(data = dataCls, aes(Cls, estimate, color = "Closeness centrality"), size = 1) +
  # Btw
  geom_ribbon(data = dataBtw, aes(x = Btw, y = estimate, ymin = lwr, ymax = upr), alpha = .3, fill = colGM[3]) +
  geom_line(data = dataBtw, aes(Btw, estimate, color =  "Betweenness centrality"), size = 1) +
  # PtC
  geom_ribbon(data = dataPtC, aes(x = PtC, y = estimate, ymin = lwr, ymax = upr), alpha = .3, fill = colGM[4]) +
  geom_line(data = dataPtC, aes(PtC, estimate, color = "Participation coefficient"), size = 1) +
  #WMz
  geom_ribbon(data = dataWMz, aes(x = WMz, y = estimate, ymin = lwr, ymax = upr), alpha = .3, fill = colGM[5]) +
  geom_line(data = dataWMz, aes(WMz, estimate, color = "Within-module degree z-score"), size = 1) +

 # xlim(-3, 3) +
  labs(x = "Centrality measure  (z-standardized)",
       y = "Resting-state timescale estimate  \n (logarithmized & z-standardized)",
       title = " ") +
  theme_bw() + scale_color_manual(name = "", values = colsGM, breaks = labels (colsGM))
GMplot

save(GMplot, dataDeg, dataCls, dataBtw, dataPtC, dataWMz, model, colGM,colsGM,
     file = paste0(outDir, "/TSRaut_GM/Plot_TSRaut.RData"))


GMplot_paper = GMplot + theme(legend.position = c(0.8,0.25)) +
  theme(axis.text=element_text(size=11,  family="Calibri", color = "black")) +
  theme(text=element_text(size=11,  family="Calibri", color = "black")) 

ggsave(paste0("GMPlot_","PaperPlot.png"), plot = GMplot_paper,
       path = ".../Results_TS", 
       width = 82.3*2, height = 50*2, units = "mm", dpi = 300)

# Surface based morphometry thickness ####
load(paste0(outDir, "/TSRaut_SBMt/Plot_TSRaut_SBMt.RData"))

colThick = "#fce303" # #fde725 ##fce303
colsThick = c("Cortical thickness" = "#fce303")

SBMt_plot <- ggplot() +
  # Thickness
  geom_point(data = model$data, aes(z_Stats_ThickAvg, z_logTS_Raut),
             alpha = .2, color = colThick, size = .5) +
  geom_ribbon(data = dataSBMt, aes(x = Thick, y = estimate, ymin = lwr, ymax = upr), alpha = .3, fill = colThick) +
  geom_line(data = dataSBMt, aes(Thick, estimate, color = "Cortical thickness"), size = 1) +
  labs(x = "Cortical thickness (z-standardized)",
       y = "Resting-state timescale estimate  \n (logarithmized & z-standardized)",
       title = " ") +
  theme_bw() + scale_color_manual(name = "", values = colsThick)
SBMt_plot

save(SBMt_plot, dataSBMt, colThick, colsThick, model, file = paste0(outDir, "/TSRaut_SBMt/Plot_TSRaut_SBMt.RData"))


SBMt_plot_paper = SBMt_plot + theme(legend.position = c (0.85,0.1)) +
  theme(axis.text=element_text(size=11,  family="Calibri", color = "black")) +
  theme(text=element_text(size=11,  family="Calibri", color = "black")) 


ggsave(paste0("SBMt_plot_","PaperPlot.png"), plot = SBMt_plot_paper,
       path =".../Results_TS", 
       width = 82.3*2, height = 50*2, units = "mm", dpi = 300)

rm(brmodel_TSRaut_SBMt)


#### Model Comparison ####
analysis_dir = ".../2_RestingState_Timescales/"
outDir = paste0(analysis_dir, "_Bayesian_Model_comp_finalZ")
models = list.files(outDir, pattern = ".RData")
models = models[c(8,10,11,13,15,17,18,20)]
modelNames = sapply(strsplit(models, ".", fixed = TRUE), "[[",1)
setwd(outDir)

library(brms)

loo_command = "loo_compare("

loo_results = data.frame(matrix(NA,ncol = 3+8, nrow = length(models)))
names(loo_results)[1:3] = c("Model", "elpd_loo", "se_elpd_loo")

for (i in 1:length(models)){  
  lazyLoad(modelNames[i])
  loo_results$Model[i] = modelNames[i]
  if (i == length(models)) {
    loo_command = paste0(loo_command, "loo_", modelNames[i], ")")
  } else {
    loo_command = paste0(loo_command, "loo_", modelNames[i], ",")}
  elpd_loo_i = eval(parse(text = paste0("loo_",modelNames[i],"$estimates")))
  loo_results$elpd_loo[i] = elpd_loo_i[1,1]
  loo_results$se_elpd_loo[i] = elpd_loo_i[1,2]
}

# Model Comparision
looComp = eval(parse(text = loo_command))

names(loo_results)[4:11] = colnames(looComp)
rows_brms = paste0("TSRaut_", sapply(strsplit(rownames(looComp), "_", fixed = TRUE), "[[",3))

for (i in 1:length(models)){
  loo_results[i,4:11] = looComp[which(rows_brms == loo_results$Model[i]),1:8]
}

save(loo_results, file = paste0(outDir, "/_loo_comparison.RData"))
write.csv(loo_results, file = paste0(outDir, "/_loo_comparison.csv"))
