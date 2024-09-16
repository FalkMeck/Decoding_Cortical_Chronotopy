### BAYESIAN ANALYSIS OF Cortical TIMESCALES by RAUT

# Libraries
library(tidyverse)
library(tidybayes)
library(brms)
library(bayesplot)
library(bayestestR)
library(ggplot2)
library(RcppEigen)

### Data
analysis_dir = ".../"
outDir = paste0(analysis_dir, "_Bayesian_Model_comp_finalZ")
ifelse(!dir.exists(file.path(outDir)), dir.create(file.path(outDir)), FALSE)
setwd(outDir)

load(paste0(analysis_dir, "TS_Raut_Data_final_lausanne250_2024_09_05.RData"))



### BASE MODEL ####
brmodel_TSRaut_BazeRC = brm(formula = z_logTS_Raut ~ 1 +
                              + z_ICV_l + z_Stats_Surf + z_Coord_z + # control variables: cortical volume, surface area/ICV, z-Coodrdinate
                              (1|Subject), # multilevel
                            family = gaussian,  data = R2data_RCRautFinal, 
                            warmup = 1000, iter = 4000, chains = 10, cores = 5,
                            control = list(adapt_delta = 0.99), seed = 15)
summary(brmodel_TSRaut_BazeRC)
loo_TSRaut_BazeRC = loo(brmodel_TSRaut_BazeRC)
save(brmodel_TSRaut_BazeRC, loo_TSRaut_BazeRC, file = paste0(outDir, "/TSRaut_BazeRC.RData"))

### RICH CLUB ####
# Rich Club Single # Warning: 
brmodel_TSRaut_RC =  brm(formula = z_logTS_Raut ~ RC_Class_Sing + 
                           + z_ICV_l + z_Stats_Surf +  z_Coord_z + # control variables: cortical volume, surface area/ICV, z-Coodrdinate
                           (1|Subject), # multilevel 
                         family = gaussian,  data = R2data_RCRautFinal, 
                         warmup = 1000, iter = 4000, chains = 10, cores = 5,
                         control = list(adapt_delta = 0.99), seed = 15)
loo_TSRaut_RC = loo(brmodel_TSRaut_RC)
save(brmodel_TSRaut_RC, loo_TSRaut_RC, file = paste0(outDir, "/TSRaut_RC.RData"))
summary(brmodel_TSRaut_RC)
rm(list = ls())
rstudioapi::restartSession()

### DIVERSE CLUB ####
# Diverse Club Single
brmodel_TSRaut_DC =  brm(formula = z_logTS_Raut ~ DC_Class_Sing + 
                           + z_ICV_l + z_Stats_Surf +  z_Coord_z + # control variables: cortical volume, surface area/ICV, z-Coodrdinate
                           (1|Subject), # multilevel  
                         family = gaussian,  data = R2data_RCRautFinal, 
                         warmup = 1000, iter = 4000, chains = 10, cores = 5,
                         control = list(adapt_delta = 0.99), seed = 15)
loo_TSRaut_DC = loo(brmodel_TSRaut_DC)
save(brmodel_TSRaut_DC, loo_TSRaut_DC, file = paste0(outDir, "/TSRaut_DC.RData"))
summary(brmodel_TSRaut_DC)
rm(list = ls())
rstudioapi::executeCommand("restartR")

### MULTIMODAL ####
## SFC
#  Seven Step Single SFC
brmodel_TSRaut_SFCindi =  brm(formula = z_logTS_Raut ~ z_SFC_S7_Multi_seed + 
                                + z_ICV_l + z_Stats_Surf + z_Coord_z +  z_Coord_z + # control variables: cortical volume, surface area/ICV, z-Coodrdinate
                                (1|Subject), # multilevel 
                              family = gaussian,  data = R2data_RCRautFinal, 
                              warmup = 1000, iter = 4000, chains = 10, cores = 5,
                              control = list(adapt_delta = 0.99), seed = 15)
loo_TSRaut_SFCindi = loo(brmodel_TSRaut_SFCindi)
save(brmodel_TSRaut_SFCindi, loo_TSRaut_SFCindi, file = paste0(outDir, "/TSRaut_SFCindi.RData"))
summary(brmodel_TSRaut_SFCindi)
rm(list = ls())
rstudioapi::executeCommand("restartR")

## Ji-Ito-Atlas
# Ji-Ito multimodal percentage # Warning: 
brmodel_TSRaut_JiItoPerc=  brm(formula = z_logTS_Raut ~ z_Ji_Ito_multi_perc + 
                                 + z_ICV_l + z_Stats_Surf + z_Coord_z + # control variables: cortical volume, surface area/ICV, z-Coodrdinate
                                 (1|Subject), # multilevel 
                               family = gaussian,  data = R2data_RCRautFinal, 
                               warmup = 1000, iter = 4000, chains = 10, cores = 5,
                               control = list(adapt_delta = 0.99), seed = 15)
loo_TSRaut_JiItoPerc = loo(brmodel_TSRaut_JiItoPerc)
save(brmodel_TSRaut_JiItoPerc, loo_TSRaut_JiItoPerc, file = paste0(outDir, "/TSRaut_JiItoPerc.RData"))
summary(brmodel_TSRaut_JiItoPerc)
rm(list = ls())
rstudioapi::executeCommand("restartR")

### GRAPH MEASURES ####
# all
brmodel_TSRaut_GM=  brm(formula = z_logTS_Raut ~ z_Degree + z_Closeness + z_Betweenness + z_PartCoef_DC + z_WithinMod_DC + 
                          + z_ICV_l + z_Stats_Surf +  z_Coord_z + # control variables: cortical volume, surface area/ICV, z-Coodrdinate
                          (1|Subject), # multilevel 
                        family = gaussian,  data = R2data_RCRautFinal, 
                        warmup = 1000, iter = 4000, chains = 10, cores = 5,
                        control = list(adapt_delta = 0.99), seed = 15)
loo_TSRaut_GM = loo(brmodel_TSRaut_GM)
save(brmodel_TSRaut_GM, loo_TSRaut_GM, file = paste0(outDir, "/TSRaut_GM.RData"))
summary(brmodel_TSRaut_GM)
rm(list = ls())
rstudioapi::executeCommand("restartR")



### Coordinates ####
# interaction of absX and Y 
brmodel_TSRaut_aXxY=  brm(formula = z_logTS_Raut ~ z_absCoord_x*z_Coord_y +
                             + z_ICV_l + z_Stats_Surf +  z_Coord_z + # control variables: cortical volume, surface area/ICV, z-Coodrdinate
                             (1|Subject), # multilevel  
                           family = gaussian,  data = R2data_RCRautFinal, 
                           warmup = 1000, iter = 4000, chains = 10, cores = 5,
                           control = list(adapt_delta = 0.99), seed = 15)
loo_TSRaut_aXxY = loo(brmodel_TSRaut_aXxY)
save(brmodel_TSRaut_aXxY, loo_TSRaut_aXxY, file = paste0(outDir, "/TSRaut_aXxY.RData"))
summary(brmodel_TSRaut_aXxY)
rm(list = ls())
rstudioapi::executeCommand("restartR")


### SBM: SURFACE and THICKNESS ####
# only Thickness
brmodel_TSRaut_SBMt=  brm(formula = z_logTS_Raut ~ z_Stats_ThickAvg + 
                            + z_ICV_l + z_Stats_Surf +  z_Coord_z + # control variables: cortical volume, surface area/ICV, z-Coodrdinate
                            (1|Subject), # multilevel 
                          family = gaussian,  data = R2data_RCRautFinal, 
                          warmup = 1000, iter = 4000, chains = 10, cores = 5,
                          control = list(adapt_delta = 0.99), seed = 15)
loo_TSRaut_SBMt = loo(brmodel_TSRaut_SBMt)
save(brmodel_TSRaut_SBMt, loo_TSRaut_SBMt, file = paste0(outDir, "/TSRaut_SBMt.RData"))
summary(brmodel_TSRaut_SBMt)
rm(list = ls())
rstudioapi::executeCommand("restartR")

###### Test ASSUMPTIONS/DISTRIBUTION #############
analysis_dir = ".../"
outDir = paste0(analysis_dir, "_Bayesian_Model_comp_finalZ")
models = list.files(outDir, pattern = ".RData")
modelNames = sapply(strsplit(models, ".", fixed = TRUE), "[[",1)
setwd(outDir)

# adjust the global plotting theme
theme_set(theme_gray(base_size = 13, base_family = "Arial") + theme(panel.grid = element_blank()))
# define own functions for ppc_stat/pp_check
dispersion <- function(x) {var(x)/mean(x)}


for (i in 2:length(modelNames)) {


# convert .RData -> .rdb/.rdx
e = local({load(models[i]); environment()})
tools:::makeLazyLoadDB(e, modelNames[i])
lazyLoad(modelNames[i])
# ls()
eval(parse(text = paste0("fit = brmodel_", modelNames[i])))

assumDir = paste0(outDir, "/", modelNames[i])
ifelse(!dir.exists(file.path(assumDir)), dir.create(file.path(assumDir)), FALSE)

# general density overlay
png(file = paste0(assumDir, "/Dens_", modelNames[i],".png"),width = 1200,height = 800)
densPlot = pp_check(fit,
         type = 'dens_overlay',
         ndraws = 2000) + theme_bw(16)
print(densPlot)
dev.off()

# mean
png(file = paste0(assumDir, "/Mean_", modelNames[i],".png"),width = 1200,height = 800)
meanPlot = pp_check(fit,
         ndraws = 2000,
         type ='stat',
         stat = "mean",
         binwidth = .001)
print(meanPlot)
dev.off()

# mean & SD
png(file = paste0(assumDir, "/M_SD_", modelNames[i],".png"),width = 1200,height = 800)
meanSDplot = pp_check(fit,
         type = "stat_2d",
         ndraws = 2000)+ # can also be 5000
  theme_bw(16)+
  xlab("Mean")+
  ylab("Standard Deviation (SD)")
print(meanSDplot)
dev.off()

# Predict from Posterior
# based on samples and prior
post_p <- posterior_predict(fit, ndraws = 2000) 

# dispersion
png(file = paste0(assumDir, "/Disp_", modelNames[i],".png"),width = 1200,height = 800)
dispPlot = ppc_stat(y = fit$data$z_logTS_Raut, 
         yrep = post_p,
         stat="dispersion", binwidth = 1)+
  theme_bw(16)
print(dispPlot)
dev.off()

png(file = paste0(assumDir, "/performance_", modelNames[i],".png"),width = 1200,height = 800)
perfomPlot = performance::check_model(fit)
print(perfomPlot)
dev.off()

eval(parse(text = paste0("rm(brmodel_", modelNames[i],",loo_",modelNames[i],",e)")))
}

rm(list = ls())
rstudioapi::executeCommand("restartR")

#### Model Comparison ####
analysis_dir = ".../"
outDir = paste0(analysis_dir, "_Bayesian_Model_comp_finalZ")
models = list.files(outDir, pattern = ".RData")
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
save.image(".../_Bayesian_Model_comp_finalZ/TSRaut_loo_comp_final.RData")
