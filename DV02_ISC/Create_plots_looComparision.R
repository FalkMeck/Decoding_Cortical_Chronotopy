# Bayesian Model Comparision ####

rm(list = ls())
rstudioapi::restartSession()

library(ggplot2)
library(ggdist)
library(extrafont)
library(ggh4x)

tsDir = ".../DV01_RS_Timescales/_Bayesian_Model_comp_finalZ/"
iscDir = ".../DV02_ISC/_Bayesian_Models_final_Z/"

looRS = read.csv(paste0(tsDir, "_loo_comparison.csv"))
looRS$org = c("Posterior/lateral-to-anterior/medial Gradient", "Null model",
              "Diverse Club Architecture", "Graph Measures of Centrality",
              "Uni-to-multi-modal Gradient (Ji atlas)", "Rich Club Architecture",
               "Cytoarchitecture (Cortical Thickness)", "Uni-to-multi-modal Gradient (SFC)")


looISC = read.csv(paste0(iscDir, "_loo_comparison_20240907.csv"))
looISC$org = c("Null model", "Diverse Club Architecture", "Graph Measures of Centrality",
               "Uni-to-multi-modal Gradient (Ji atlas)", "Rich Club Architecture",
               "Cytoarchitecture (Cortical Thickness)", "Uni-to-multi-modal Gradient (SFC)",  
               "Posterior/lateral-to-anterior/medial Gradient")
looRS$org = factor(looRS$org, levels = c("Null model", "Rich Club Architecture", "Diverse Club Architecture",
                                         "Uni-to-multi-modal Gradient (SFC)", "Uni-to-multi-modal Gradient (Ji atlas)", 
                                         "Posterior/lateral-to-anterior/medial Gradient","Graph Measures of Centrality",
                                         "Cytoarchitecture (Cortical Thickness)"))
looISC$org = factor(looISC$org, levels = c("Null model", "Rich Club Architecture", "Diverse Club Architecture",
                                           "Uni-to-multi-modal Gradient (SFC)",  "Uni-to-multi-modal Gradient (Ji atlas)",
                                           "Posterior/lateral-to-anterior/medial Gradient", "Graph Measures of Centrality",
                                           "Cytoarchitecture (Cortical Thickness)"))

colsmodels = c("Null model" = "#a6a6a6", "Rich Club Architecture" = "#440154", "Diverse Club Architecture" = "#404788", 
               "Uni-to-multi-modal Gradient (SFC)"  = "#1c5863", "Uni-to-multi-modal Gradient (Ji atlas)" = "#38afc7", 
               "Posterior/lateral-to-anterior/medial Gradient" = "#29AF7F",
               "Graph Measures of Centrality" = "#95d840", "Cytoarchitecture (Cortical Thickness)" = "#fce303")
shortModelsNames = c("Null model" = "Null", "Rich Club Architecture" = "RC", "Diverse Club Architecture" = "DC",
                     "Uni-to-multi-modal Gradient (SFC)"  = "SFC","Uni-to-multi-modal Gradient (Ji atlas)" = "Ji", 
                     "Posterior/lateral-to-anterior/medial Gradient" = "PLAM",
                     "Graph Measures of Centrality" = "GM",
                     "Cytoarchitecture (Cortical Thickness)" = "CT")

dataBC = data.frame(looRS$org, as.numeric(looRS$elpd_loo), as.numeric(looRS$se_elpd_loo))
names(dataBC) = c("Organizational_Scheme", "elpd_rs", "se_elpd_rs")

dataBC$elpd_isc = NA
dataBC$se_elpd_isc = NA
for (i in 1:length(dataBC$Organizational_Scheme)) {
  dataBC$elpd_isc[i] = looISC$elpd_loo[looISC$org == dataBC$Organizational_Scheme[i]] 
  dataBC$se_elpd_isc[i] = looISC$se_elpd_loo[looISC$org == dataBC$Organizational_Scheme[i]] 
}
dataBC$elpd_rs_max = dataBC$elpd_rs + dataBC$se_elpd_rs
dataBC$elpd_rs_min = dataBC$elpd_rs - dataBC$se_elpd_rs

dataBC$elpd_isc_max = dataBC$elpd_isc + dataBC$se_elpd_isc
dataBC$elpd_isc_min = dataBC$elpd_isc - dataBC$se_elpd_isc

dataLong = data.frame(c(dataBC$Organizational_Scheme, dataBC$Organizational_Scheme), 
                      c(dataBC$elpd_rs, dataBC$elpd_isc), c(dataBC$se_elpd_rs, dataBC$se_elpd_isc),
                      c(dataBC$elpd_rs_min, dataBC$elpd_isc_min), c(dataBC$elpd_rs_max, dataBC$elpd_isc_max))
names(dataLong) = c("OrgaScheme", "elpd", "se_elpd", "min_eldp", "max_elpd")
dataLong$analysis = c(rep("RS-Timescales", 8), rep("ISC", 8))
dataLong$analysis = factor(dataLong$analysis, levels = c("RS-Timescales", "ISC"))


BC_rs = ggplot(aes(x= Organizational_Scheme, y = elpd_rs+ 15000, fill = Organizational_Scheme), data = dataBC) +
  geom_bar(stat = "identity") + scale_fill_manual(name = "", values = colsmodels) +
  geom_errorbar(aes(ymin=elpd_rs_min+ 15000, ymax=elpd_rs_max+ 15000), width=.5, linewidth = .75) + coord_cartesian(ylim=c(3000,3300))+
  scale_y_continuous(labels = c(3000-15000, 3100-15000, 3200-15000, 3300-15000)) + ylab("elpd") + xlab("Organizational scheme") + # Expected log Pointwise Predictive Density
  scale_x_discrete(labels = shortModelsNames) +
  theme_classic() + theme(legend.position = "none") +
  theme(axis.text=element_text(size=22,  family="Calibri", color = "black")) +
  theme(text=element_text(size=22,  family="Calibri", color = "black"))

BC_isc = ggplot(aes(x= Organizational_Scheme, y = elpd_isc, fill = Organizational_Scheme), data = dataBC) +
  geom_bar(stat = "identity") + scale_fill_manual(name = "", values = colsmodels) +
  geom_errorbar(aes(ymin=elpd_isc_min, ymax=elpd_isc_max), width=.5, linewidth = .75) + coord_cartesian(ylim=c(90300,91300))+
  scale_y_continuous(position = "left") + ylab("elpd") + xlab("Organizational scheme") + # Expected log Pointwise Predictive Density
  scale_x_discrete(labels = shortModelsNames) +
  theme_classic() + theme(legend.position = "none") +
  theme(axis.text=element_text(size=22,  family="Calibri", color = "black")) +
  theme(text=element_text(size=22,  family="Calibri", color = "black"))

library(ggpubr)
ggarrange(BC_rs, BC_isc, labels = c("Resting-state timescales", "Inter-subject correlation"), nrow = 1, ncol = 2,
          font.label = list(size=22,  family="Calibri", color = "black"))

ggsave("BC_Plot_TimeScales_order.png", plot = BC_rs,
       path = "Q:/fmecklenbrauck/07_WS23_24/RICHIE2/Review/Repeat/RICHIE2_Review_Repeat", 
       width = 82.3*2, height = 50*2, units = "mm", dpi = 300)
ggsave("BC_Plot_ISC_order.png", plot = BC_isc,
       path = "Q:/fmecklenbrauck/07_WS23_24/RICHIE2/Review/Repeat/RICHIE2_Review_Repeat", 
       width = 82.3*2, height = 50*2, units = "mm", dpi = 300)
