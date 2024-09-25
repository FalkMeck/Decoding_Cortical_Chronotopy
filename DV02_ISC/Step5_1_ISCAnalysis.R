###################################
#####                         #####
#####       ISC ANALYSIS      #####
#####                         #####
###################################

## BAYESIAN ####
library(parallel)
library(rstan)
library(brms)
library(bayesplot)
library(rstanarm)

library(car)

library(tidyverse)
library(tidybayes)
library(bayestestR)


iscDir = "/..."
curDate = "_2024_09_05"
load(paste0(iscDir, "/iscR2data_Final", curDate , ".RData")) # use the provided data set, if you want to replicate
outDir = paste0(iscDir, "/_Bayesian_Models_final_Z")
ifelse(!dir.exists(file.path(outDir)), dir.create(file.path(outDir)), FALSE)
setwd(outDir)

# MODELS ####

# RICH CLUB ####
brmSNV_ISCz_CxRC = brm(formula = ISC_FishZ  ~ 1 + Condition * RC_Class_Sing + 
                       z_Stats_Surf + z_ICV_l + z_Coord_z + 
                       (1|Subject), 
                       data = iscR2data_Final,
                       prior = c(set_prior("normal(0,5)", class = "b"),
                               set_prior("student_t(3,0,2.5)", class = "sd")), 
                      family = skew_normal(link = "identity", link_sigma = "log", link_alpha = "identity"), 
                      warmup = 1000, iter = 4000, chains = 10, cores = 5,
                      control = list(adapt_delta = 0.8, max_treedepth = 15), seed = 15, save_pars = save_pars(all=TRUE), sample_prior = TRUE)
save(brmSNV_ISCz_CxRC, file = paste0(outDir, "/brmISCz_CxRC.RData"))
summary(brmSNV_ISCz_CxRC)
rm(list = ls())
rstudioapi::restartSession()

# Diverse CLUB ####
brmSNV_ISCz_CxDC = brm(formula = ISC_FishZ  ~ 1 + Condition * DC_Class_Sing + 
                         z_Stats_Surf + z_ICV_l + z_Coord_z + 
                         (1|Subject), 
                       data = iscR2data_Final,
                       prior = c(set_prior("normal(0,5)", class = "b"),
                                 set_prior("student_t(3,0,2.5)", class = "sd")), 
                       family = skew_normal(link = "identity", link_sigma = "log", link_alpha = "identity"), 
                       warmup = 1000, iter = 4000, chains = 10, cores = 5,
                       control = list(adapt_delta = 0.8, max_treedepth = 15), seed = 15, save_pars = save_pars(all=TRUE), sample_prior = TRUE)
save(brmSNV_ISCz_CxDC, file = paste0(outDir, "/brmISCz_CxDC.RData"))
summary(brmSNV_ISCz_CxDC)
rm(list = ls())
rstudioapi::restartSession()

# Uni-multi ####
# SFC
# brmSNV_ISCz_CxSFC = brm(formula = ISC_FishZ  ~ 1 + Condition * z_SFC_S7_Multi_seed + 
#                           z_Stats_Surf + z_ICV_l +  z_Coord_z + 
#                           (1|Subject), 
#                         data = iscR2data_Final,
#                         prior = c(set_prior("normal(0,5)", class = "b"),
#                                   set_prior("student_t(3,0,2.5)", class = "sd")), 
#                         family = skew_normal(link = "identity", link_sigma = "log", link_alpha = "identity"), 
#                         warmup = 1000, iter = 4000, chains = 10, cores = 5,
#                         control = list(adapt_delta = 0.8, max_treedepth = 15), seed = 15, save_pars = save_pars(all=TRUE), sample_prior = TRUE)
# save(brmSNV_ISCz_CxSFC, file = paste0(outDir, "/brmISCz_CxSFC.RData"))
# summary(brmSNV_ISCz_CxSFC)
# rm(list = ls())
rstudioapi::restartSession()

# EDIT: SFC Diff 7 1
brmSNV_ISCz_CxSFCDiff = brm(formula = ISC_FishZ  ~ 1 + Condition * z_SFC_S7_1_Diff_Multi_seed +
                              z_Stats_Surf + z_ICV_l +  z_Coord_z +
                              (1|Subject),
                            data = iscR2data_Final,
                            prior = c(set_prior("normal(0,5)", class = "b"),
                                      set_prior("student_t(3,0,2.5)", class = "sd")),
                            family = skew_normal(link = "identity", link_sigma = "log", link_alpha = "identity"),
                            warmup = 1000, iter = 4000, chains = 10, cores = 5,
                            control = list(adapt_delta = 0.8, max_treedepth = 15), seed = 15, save_pars = save_pars(all=TRUE), sample_prior = TRUE)
save(brmSNV_ISCz_CxSFCDiff, file = paste0(outDir, "/brmISCz_CxSFCDiff.RData"))
summary(brmSNV_ISCz_CxSFCDiff)
rm(list = ls())
rstudioapi::restartSession()


# Ji Ito
brmSNV_ISCz_CxJiIto = brm(formula = ISC_FishZ  ~ 1 + Condition * z_Ji_Ito_multi_perc + 
                            z_Stats_Surf + z_ICV_l +  z_Coord_z + 
                            (1|Subject), 
                          data = iscR2data_Final,
                          prior = c(set_prior("normal(0,5)", class = "b"),
                                    set_prior("student_t(3,0,2.5)", class = "sd")), 
                          family = skew_normal(link = "identity", link_sigma = "log", link_alpha = "identity"), 
                          warmup = 1000, iter = 4000, chains = 10, cores = 5,
                          control = list(adapt_delta = 0.8, max_treedepth = 15), seed = 15, save_pars = save_pars(all=TRUE), sample_prior = TRUE)
save(brmSNV_ISCz_CxJiIto, file = paste0(outDir, "/brmISCz_CxJiIto.RData"))
summary(brmSNV_ISCz_CxJiIto)
rm(list = ls())
rstudioapi::restartSession()

# X, Y ####
brmSNV_ISCz_CxXY = brm(formula = ISC_FishZ  ~ 1 + Condition * (z_absCoord_x *z_Coord_y) + z_Coord_z + 
                        z_Stats_Surf + z_ICV_l +  
                        (1|Subject), 
                      data = iscR2data_Final,
                      prior = c(set_prior("normal(0,5)", class = "b"),
                                set_prior("student_t(3,0,2.5)", class = "sd")), 
                      family = skew_normal(link = "identity", link_sigma = "log", link_alpha = "identity"), 
                      warmup = 1000, iter = 4000, chains = 10, cores = 5,
                      control = list(adapt_delta = 0.8, max_treedepth = 15), seed = 15, save_pars = save_pars(all=TRUE), sample_prior = TRUE)
save(brmSNV_ISCz_CxXY, file = paste0(outDir, "/brmISCz_CxXY.RData"))
summary(brmSNV_ISCz_CxXY)
rm(list = ls())
rstudioapi::restartSession()

# Graph Measures ####
brmSNV_ISCz_CxGM = brm(formula = ISC_FishZ  ~ 1 + Condition * (z_Degree+ z_Closeness + z_Betweenness +z_PartCoef_DC +z_WithinMod_DC) + 
                           z_Stats_Surf + z_ICV_l +  z_Coord_z + 
                           (1|Subject), 
                         data = iscR2data_Final,
                         prior = c(set_prior("normal(0,5)", class = "b"),
                                   set_prior("student_t(3,0,2.5)", class = "sd")), 
                         family = skew_normal(link = "identity", link_sigma = "log", link_alpha = "identity"), 
                         warmup = 1000, iter = 4000, chains = 10, cores = 5,
                         control = list(adapt_delta = 0.8, max_treedepth = 15), seed = 15, save_pars = save_pars(all=TRUE), sample_prior = TRUE)
save(brmSNV_ISCz_CxGM, file = paste0(outDir, "/brmISCz_CxGM.RData"))
summary(brmSNV_ISCz_CxGM)
rm(list = ls())
rstudioapi::restartSession()

# Cortical Thickness ####
brmSNV_ISCz_CxSBMt = brm(formula = ISC_FishZ  ~ 1 + Condition * z_Stats_Thick + 
                           z_Stats_Surf + z_ICV_l +  z_Coord_z + 
                           (1|Subject), 
                         data = iscR2data_Final,
                         prior = c(set_prior("normal(0,5)", class = "b"),
                                   set_prior("student_t(3,0,2.5)", class = "sd")), 
                         family = skew_normal(link = "identity", link_sigma = "log", link_alpha = "identity"), 
                         warmup = 1000, iter = 4000, chains = 10, cores = 5,
                         control = list(adapt_delta = 0.8, max_treedepth = 15), seed = 15, save_pars = save_pars(all=TRUE), sample_prior = TRUE)
save(brmSNV_ISCz_CxSBMt, file = paste0(outDir, "/brmISCz_CxSBMt.RData"))
summary(brmSNV_ISCz_CxSBMt)
rm(list = ls())
rstudioapi::restartSession()

# Null model ####
brmSNV_ISCz_C_Null = brm(formula = ISC_FishZ  ~ 1 + Condition + 
                        z_Stats_Surf + z_ICV_l + z_Coord_z + 
                        (1|Subject), 
                      data = iscR2data_Final,
                      prior = c(set_prior("normal(0,5)", class = "b"),
                                set_prior("student_t(3,0,2.5)", class = "sd")), 
                      family = skew_normal(link = "identity", link_sigma = "log", link_alpha = "identity"), 
                      warmup = 1000, iter = 4000, chains = 10, cores = 5,
                      control = list(adapt_delta = 0.8, max_treedepth = 15), seed = 15, save_pars = save_pars(all=TRUE), sample_prior = TRUE)
save(brmSNV_ISCz_C_Null, file = paste0(outDir, "/brmISCz_C_Null.RData"))
summary(brmSNV_ISCz_C_Null)

rm(list = ls())
rstudioapi::restartSession()

# Model Comparison ####

library(brms)

iscDir = "/..."
outDir = paste0(iscDir, "/_Bayesian_Models_final_Z")
setwd(outDir)

looModels = list.files(path = outDir, pattern = glob2rx("loo*RData"))
looNames= stringr::str_split(looModels, pattern = "[.]")
looNames = sapply(looNames,"[[",1)

loo_command = "loo_compare("

loo_results = data.frame(matrix(NA,ncol = 3+8, nrow = length(looModels)))
names(loo_results)[1:3] = c("Model", "elpd_loo", "se_elpd_loo")

for (i in 1:length(looModels)){  
  load(looModels[i])
  loo_results$Model[i] = looNames[i]
  if (i == length(looNames)) {
    loo_command = paste0(loo_command, looNames[i], ")")
  } else {
    loo_command = paste0(loo_command, looNames[i], ",")}
  elpd_loo_i = eval(parse(text = paste0(looNames[i],"$estimates")))
  loo_results$elpd_loo[i] = elpd_loo_i[1,1]
  loo_results$se_elpd_loo[i] = elpd_loo_i[1,2]
}

# Model Comparision
looComp = eval(parse(text = loo_command))

names(loo_results)[4:11] = colnames(looComp)
rows_brms = paste0("loo_brm", substring(row.names(looComp), 8))

for (i in 1:length(looNames)){
  loo_results[i,4:11] = looComp[which(rows_brms == loo_results$Model[i]),1:8]
}

save(loo_results, file = paste0(outDir, "/_loo_comparison_20240506.RData"))
write.csv(loo_results, file = paste0(outDir, "/_loo_comparison_20240506.csv"))

