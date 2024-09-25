### ISC hypothesis tests 
setwd(".../_Bayesian_Models_final_Z")

# libraries
library(brms)
library(writexl)

# Rich Club ####
load("brmISCz_CxRC.RData")
summary(brmSNV_ISCz_CxRC)

# Main effect RC in different conditions
h3 <- c( "RC_Class_SingRC_Feeder > 0", "RC_Class_SingRC_Club > 0", "RC_Class_SingRC_Club - RC_Class_SingRC_Feeder > 0",

         "RC_Class_SingRC_Feeder + ConditionTriplet:RC_Class_SingRC_Feeder  > 0",
         "RC_Class_SingRC_Club  + ConditionTriplet:RC_Class_SingRC_Club > 0",
         "RC_Class_SingRC_Club + ConditionTriplet:RC_Class_SingRC_Club - RC_Class_SingRC_Feeder - ConditionTriplet:RC_Class_SingRC_Feeder > 0",

         "RC_Class_SingRC_Feeder + ConditionNonet:RC_Class_SingRC_Feeder  > 0",
         "RC_Class_SingRC_Club  + ConditionNonet:RC_Class_SingRC_Club > 0",
         "RC_Class_SingRC_Club + ConditionNonet:RC_Class_SingRC_Club - RC_Class_SingRC_Feeder - ConditionNonet:RC_Class_SingRC_Feeder > 0",

         "RC_Class_SingRC_Feeder + ConditionComplete:RC_Class_SingRC_Feeder  > 0",
         "RC_Class_SingRC_Club  + ConditionComplete:RC_Class_SingRC_Club > 0",
         "RC_Class_SingRC_Club + ConditionComplete:RC_Class_SingRC_Club - RC_Class_SingRC_Feeder - ConditionComplete:RC_Class_SingRC_Feeder > 0")
h3RCeffect_inCondtions = brms::hypothesis(brmSNV_ISCz_CxRC, h3)
write_xlsx(h3RCeffect_inCondtions$hypothesis, path = paste0(getwd(),"/h3RCeffect_inCondtions","_ISC_RC",".xlsx"))

# Main effect RC in different conditions
h3 <- c( "RC_Class_SingRC_Feeder = 0", "RC_Class_SingRC_Club = 0", "RC_Class_SingRC_Club - RC_Class_SingRC_Feeder = 0",
         
         "RC_Class_SingRC_Feeder + ConditionTriplet:RC_Class_SingRC_Feeder  = 0",
         "RC_Class_SingRC_Club  + ConditionTriplet:RC_Class_SingRC_Club = 0",
         "RC_Class_SingRC_Club + ConditionTriplet:RC_Class_SingRC_Club - RC_Class_SingRC_Feeder - ConditionTriplet:RC_Class_SingRC_Feeder = 0",
         
         "RC_Class_SingRC_Feeder + ConditionNonet:RC_Class_SingRC_Feeder  = 0",
         "RC_Class_SingRC_Club  + ConditionNonet:RC_Class_SingRC_Club = 0",
         "RC_Class_SingRC_Club + ConditionNonet:RC_Class_SingRC_Club - RC_Class_SingRC_Feeder - ConditionNonet:RC_Class_SingRC_Feeder = 0",
         
         "RC_Class_SingRC_Feeder + ConditionComplete:RC_Class_SingRC_Feeder  = 0",
         "RC_Class_SingRC_Club  + ConditionComplete:RC_Class_SingRC_Club = 0",
         "RC_Class_SingRC_Club + ConditionComplete:RC_Class_SingRC_Club - RC_Class_SingRC_Feeder - ConditionComplete:RC_Class_SingRC_Feeder = 0")
h3RCeffect_inCondtions = brms::hypothesis(brmSNV_ISCz_CxRC, h3)
write_xlsx(h3RCeffect_inCondtions$hypothesis, path = paste0(getwd(),"/h3RCeffect_inCondtions","_ISC_RC","undirected",".xlsx"))


# Main effect of condtions in different levels
h4 = c( # in local
      "ConditionTriplet =  0", # a
      "ConditionTriplet >  0", # b
      "ConditionNonet - ConditionTriplet = 0",
      "ConditionComplete - ConditionNonet = 0",
        # #in feeder
        "ConditionTriplet + ConditionTriplet:RC_Class_SingRC_Feeder >  0", # Triplet > Single
        "ConditionNonet + ConditionNonet:RC_Class_SingRC_Feeder > 0", # Nonet > Single
        "ConditionComplete + ConditionComplete:RC_Class_SingRC_Feeder > 0", # Complete > Single
        "ConditionNonet + ConditionNonet:RC_Class_SingRC_Feeder - ConditionTriplet - ConditionTriplet:RC_Class_SingRC_Feeder =  0", # a1 Nonet >= Triplet
      "ConditionNonet + ConditionNonet:RC_Class_SingRC_Feeder - ConditionTriplet - ConditionTriplet:RC_Class_SingRC_Feeder >  0", #b1
        "ConditionComplete + ConditionComplete:RC_Class_SingRC_Feeder - ConditionNonet - ConditionNonet:RC_Class_SingRC_Feeder = 0", # a2 Complete >= Nonet
      "ConditionComplete + ConditionComplete:RC_Class_SingRC_Feeder - ConditionNonet - ConditionNonet:RC_Class_SingRC_Feeder > 0", # b2
        # # in RC
        "ConditionTriplet + ConditionTriplet:RC_Class_SingRC_Club =  0", #a1 Triplet (>)= Single
      "ConditionTriplet + ConditionTriplet:RC_Class_SingRC_Club >  0", #(b1)
        "ConditionNonet + ConditionNonet:RC_Class_SingRC_Club - ConditionTriplet - ConditionTriplet:RC_Class_SingRC_Club = 0", # a2 Nonet >= Triplet
      "ConditionNonet + ConditionNonet:RC_Class_SingRC_Club - ConditionTriplet - ConditionTriplet:RC_Class_SingRC_Club > 0", #b2
      "ConditionNonet + ConditionNonet:RC_Class_SingRC_Club >  0", # Nonet > Single
      
      "ConditionComplete + ConditionComplete:RC_Class_SingRC_Club > 0",  # Complete > Single
      "ConditionComplete + ConditionComplete:RC_Class_SingRC_Club - ConditionTriplet - ConditionTriplet:RC_Class_SingRC_Club = 0", # Complete >(=) Triplet
        "ConditionComplete + ConditionComplete:RC_Class_SingRC_Club - ConditionNonet - ConditionNonet:RC_Class_SingRC_Club = 0", #a3 Complete >= Nonet
        "ConditionComplete + ConditionComplete:RC_Class_SingRC_Club - ConditionNonet - ConditionNonet:RC_Class_SingRC_Club > 0" #b3
      )
h4RCondtionsEffect_inRC = brms::hypothesis(brmSNV_ISCz_CxRC, h4)
write_xlsx(h4RCondtionsEffect_inRC$hypothesis, path = paste0(getwd(),"/h4RCondtionsEffect_inRC","_ISC_RC",".xlsx"))


# Diverse Club ####
load("brmISCz_CxDC.RData")
summary(brmSNV_ISCz_CxDC)

# Main effect CC in different conditions
h3 <- c( "DC_Class_SingDC_Feeder > 0", "DC_Class_SingDC_Club > 0", "DC_Class_SingDC_Club - DC_Class_SingDC_Feeder > 0",
         
         "DC_Class_SingDC_Feeder + ConditionTriplet:DC_Class_SingDC_Feeder  > 0",
         "DC_Class_SingDC_Club  + ConditionTriplet:DC_Class_SingDC_Club > 0",
         "DC_Class_SingDC_Club + ConditionTriplet:DC_Class_SingDC_Club - DC_Class_SingDC_Feeder - ConditionTriplet:DC_Class_SingDC_Feeder > 0",
         
         "DC_Class_SingDC_Feeder + ConditionNonet:DC_Class_SingDC_Feeder  > 0",
         "DC_Class_SingDC_Club  + ConditionNonet:DC_Class_SingDC_Club > 0",
         "DC_Class_SingDC_Club + ConditionNonet:DC_Class_SingDC_Club - DC_Class_SingDC_Feeder - ConditionNonet:DC_Class_SingDC_Feeder > 0",
         
         "DC_Class_SingDC_Feeder + ConditionComplete:DC_Class_SingDC_Feeder  > 0",
         "DC_Class_SingDC_Club  + ConditionComplete:DC_Class_SingDC_Club > 0",
         "DC_Class_SingDC_Club + ConditionComplete:DC_Class_SingDC_Club - DC_Class_SingDC_Feeder - ConditionComplete:DC_Class_SingDC_Feeder > 0")
h3DCeffect_inCondtions = brms::hypothesis(brmSNV_ISCz_CxDC, h3)
write_xlsx(h3DCeffect_inCondtions$hypothesis, path = paste0(getwd(),"/h3DCeffect_inCondtions","_ISC_DC",".xlsx"))

# Diverse Club ####
load("brmISCz_CxDC.RData")
summary(brmSNV_ISCz_CxDC)

# Main effect CC in different conditions
h3 <- c( "DC_Class_SingDC_Feeder = 0", "DC_Class_SingDC_Club = 0", "DC_Class_SingDC_Club - DC_Class_SingDC_Feeder = 0",
         
         "DC_Class_SingDC_Feeder + ConditionTriplet:DC_Class_SingDC_Feeder  = 0",
         "DC_Class_SingDC_Club  + ConditionTriplet:DC_Class_SingDC_Club = 0",
         "DC_Class_SingDC_Club + ConditionTriplet:DC_Class_SingDC_Club - DC_Class_SingDC_Feeder - ConditionTriplet:DC_Class_SingDC_Feeder = 0",
         
         "DC_Class_SingDC_Feeder + ConditionNonet:DC_Class_SingDC_Feeder  = 0",
         "DC_Class_SingDC_Club  + ConditionNonet:DC_Class_SingDC_Club = 0",
         "DC_Class_SingDC_Club + ConditionNonet:DC_Class_SingDC_Club - DC_Class_SingDC_Feeder - ConditionNonet:DC_Class_SingDC_Feeder = 0",
         
         "DC_Class_SingDC_Feeder + ConditionComplete:DC_Class_SingDC_Feeder  = 0",
         "DC_Class_SingDC_Club  + ConditionComplete:DC_Class_SingDC_Club = 0",
         "DC_Class_SingDC_Club + ConditionComplete:DC_Class_SingDC_Club - DC_Class_SingDC_Feeder - ConditionComplete:DC_Class_SingDC_Feeder = 0")
h3DCeffect_inCondtions = brms::hypothesis(brmSNV_ISCz_CxDC, h3)
write_xlsx(h3DCeffect_inCondtions$hypothesis, path = paste0(getwd(),"/h3DCeffect_inCondtions","_ISC_DC","undirected",".xlsx"))

# Main effect of condtions in different levels
h4 = c( # in local
  "ConditionTriplet =  0", # a
  "ConditionTriplet >  0", # b
  "ConditionNonet - ConditionTriplet = 0",
  "ConditionComplete - ConditionNonet = 0",
  # #in feeder
  "ConditionTriplet + ConditionTriplet:DC_Class_SingDC_Feeder >  0", # Triplet > Single
  "ConditionNonet + ConditionNonet:DC_Class_SingDC_Feeder > 0", # Nonet > Single
  "ConditionComplete + ConditionComplete:DC_Class_SingDC_Feeder > 0", # Complete > Single
  "ConditionNonet + ConditionNonet:DC_Class_SingDC_Feeder - ConditionTriplet - ConditionTriplet:DC_Class_SingDC_Feeder =  0", # a1 Nonet >= Triplet
  "ConditionNonet + ConditionNonet:DC_Class_SingDC_Feeder - ConditionTriplet - ConditionTriplet:DC_Class_SingDC_Feeder >  0", #b1
  "ConditionComplete + ConditionComplete:DC_Class_SingDC_Feeder - ConditionNonet - ConditionNonet:DC_Class_SingDC_Feeder = 0", # a2 Complete >= Nonet
  "ConditionComplete + ConditionComplete:DC_Class_SingDC_Feeder - ConditionNonet - ConditionNonet:DC_Class_SingDC_Feeder > 0", # b2
  # # in DC
  "ConditionTriplet + ConditionTriplet:DC_Class_SingDC_Club =  0", #a1 Triplet (>)= Single
  "ConditionTriplet + ConditionTriplet:DC_Class_SingDC_Club >  0", #(b1)
  "ConditionNonet + ConditionNonet:DC_Class_SingDC_Club - ConditionTriplet - ConditionTriplet:DC_Class_SingDC_Club = 0", # a2 Nonet >= Triplet
  "ConditionNonet + ConditionNonet:DC_Class_SingDC_Club - ConditionTriplet - ConditionTriplet:DC_Class_SingDC_Club > 0", #b2
  "ConditionNonet + ConditionNonet:DC_Class_SingDC_Club >  0", # Nonet > Single
  "ConditionComplete + ConditionComplete:DC_Class_SingDC_Club > 0",  # Complete > Single
  "ConditionComplete + ConditionComplete:DC_Class_SingDC_Club - ConditionTriplet - ConditionTriplet:DC_Class_SingDC_Club = 0", # Complete >(=) Triplet
  "ConditionComplete + ConditionComplete:DC_Class_SingDC_Club - ConditionNonet - ConditionNonet:DC_Class_SingDC_Club = 0", #a3 Complete >= Nonet
  "ConditionComplete + ConditionComplete:DC_Class_SingDC_Club - ConditionNonet - ConditionNonet:DC_Class_SingDC_Club > 0" #b3
)
h4DCondtionsEffect_inDC = brms::hypothesis(brmSNV_ISCz_CxDC, h4)
write_xlsx(h4DCondtionsEffect_inDC$hypothesis, path = paste0(getwd(),"/h4DCondtionsEffect_inDC","_ISC_DC",".xlsx"))

# SFC Multimodal ####
load("brmISCz_CxSFC.RData")
summary(brmSNV_ISCz_CxSFC)

# Interaction effect, I think I just have to test the interactions...
h3 = c("ConditionTriplet:z_SFC_S7_Multi_seed  > 0",
       "ConditionNonet:z_SFC_S7_Multi_seed - ConditionTriplet:z_SFC_S7_Multi_seed > 0",
       "ConditionComplete:z_SFC_S7_Multi_seed - ConditionNonet:z_SFC_S7_Multi_seed > 0")
h3SFCInteraction = brms::hypothesis(brmSNV_ISCz_CxSFC, h3)
write_xlsx(h3SFCInteraction$hypothesis, path = paste0(getwd(),"/h3SFCInteraction","_ISC_SFC",".xlsx"))

# Interaction effect, I think I just have to test the interactions...
h4 = c("z_SFC_S7_Multi_seed  > 0", # Single
       "z_SFC_S7_Multi_seed + ConditionTriplet:z_SFC_S7_Multi_seed > 0", # Triplet
       "z_SFC_S7_Multi_seed + ConditionNonet:z_SFC_S7_Multi_seed > 0", # Nonet
       "z_SFC_S7_Multi_seed + ConditionComplete:z_SFC_S7_Multi_seed > 0")
h4SFCinCond = brms::hypothesis(brmSNV_ISCz_CxSFC, h4)
write_xlsx(h4SFCinCond$hypothesis, path = paste0(getwd(),"/h4SFCinCond","_ISC_SFC",".xlsx"))

h4 = c("z_SFC_S7_Multi_seed  = 0", # Single
       "z_SFC_S7_Multi_seed + ConditionTriplet:z_SFC_S7_Multi_seed = 0", # Triplet
       "z_SFC_S7_Multi_seed + ConditionNonet:z_SFC_S7_Multi_seed = 0", # Nonet
       "z_SFC_S7_Multi_seed + ConditionComplete:z_SFC_S7_Multi_seed = 0")
h4SFCinCond = brms::hypothesis(brmSNV_ISCz_CxSFC, h4)
write_xlsx(h4SFCinCond$hypothesis, path = paste0(getwd(),"/h4SFCinCond","_ISC_SFC","undirected",".xlsx"))

# # SFC Multimodal EDIT: Difference Map 7 - Map 1 ####
load("brmISCz_CxSFCDiff.RData")
summary(brmSNV_ISCz_CxSFCDiff)

# Interaction effect, I think I just have to test the interactions...
h3 = c("ConditionTriplet:z_SFC_S7_1_Diff_Multi_seed   > 0",
       "ConditionNonet:z_SFC_S7_1_Diff_Multi_seed  - ConditionTriplet:z_SFC_S7_1_Diff_Multi_seed  > 0",
       "ConditionComplete:z_SFC_S7_1_Diff_Multi_seed  - ConditionNonet:z_SFC_S7_1_Diff_Multi_seed  > 0")
h3SFCInteraction = brms::hypothesis(brmSNV_ISCz_CxSFCDiff, h3)
write_xlsx(h3SFCInteraction$hypothesis, path = paste0(getwd(),"/h3SFCInteraction","_ISC_SFCDiff",".xlsx"))

# Interaction effect, I think I just have to test the interactions...
h4 = c("z_SFC_S7_1_Diff_Multi_seed   > 0", # Single
       "z_SFC_S7_1_Diff_Multi_seed  + ConditionTriplet:z_SFC_S7_1_Diff_Multi_seed  > 0", # Triplet
       "z_SFC_S7_1_Diff_Multi_seed  + ConditionNonet:z_SFC_S7_1_Diff_Multi_seed  > 0", # Nonet
       "z_SFC_S7_1_Diff_Multi_seed  + ConditionComplete:z_SFC_S7_1_Diff_Multi_seed  > 0")
h4SFCinCond = brms::hypothesis(brmSNV_ISCz_CxSFCDiff, h4)
write_xlsx(h4SFCinCond$hypothesis, path = paste0(getwd(),"/h4SFCinCond","_ISC_SFCDiff",".xlsx"))

h4 = c("z_SFC_S7_1_Diff_Multi_seed   = 0", # Single
       "z_SFC_S7_1_Diff_Multi_seed  + ConditionTriplet:z_SFC_S7_1_Diff_Multi_seed  = 0", # Triplet
       "z_SFC_S7_1_Diff_Multi_seed  + ConditionNonet:z_SFC_S7_1_Diff_Multi_seed  = 0", # Nonet
       "z_SFC_S7_1_Diff_Multi_seed  + ConditionComplete:z_SFC_S7_1_Diff_Multi_seed  = 0")
h4SFCinCond = brms::hypothesis(brmSNV_ISCz_CxSFCDiff, h4)
write_xlsx(h4SFCinCond$hypothesis, path = paste0(getwd(),"/h4SFCinCond","_ISC_SFCDiff","undirected",".xlsx"))

# Ji Ito Multimodal ####
load("brmISCz_CxJiIto.RData")
summary(brmSNV_ISCz_CxJiIto)

# Interaction effect, I think I just have to test the interactions...
h3 = c("ConditionTriplet:z_Ji_Ito_multi_perc   > 0",
       "ConditionNonet:z_Ji_Ito_multi_perc  - ConditionTriplet:z_Ji_Ito_multi_perc  > 0",
       "ConditionComplete:z_Ji_Ito_multi_perc  - ConditionNonet:z_Ji_Ito_multi_perc  > 0")
h3JiItoInteraction = brms::hypothesis(brmSNV_ISCz_CxJiIto, h3)
write_xlsx(h3JiItoInteraction$hypothesis, path = paste0(getwd(),"/h3JiItoInteraction","_ISC_JiIto",".xlsx"))

# Interaction effect, I think I just have to test the interactions...
h4 = c("z_Ji_Ito_multi_perc  > 0", # Single
       "z_Ji_Ito_multi_perc + ConditionTriplet:z_Ji_Ito_multi_perc > 0", # Triplet
       "z_Ji_Ito_multi_perc + ConditionNonet:z_Ji_Ito_multi_perc > 0", # Nonet
       "z_Ji_Ito_multi_perc + ConditionComplete:z_Ji_Ito_multi_perc > 0") # complete
h4JiItoInCond = brms::hypothesis(brmSNV_ISCz_CxJiIto, h4)
write_xlsx(h4JiItoInCond$hypothesis, path = paste0(getwd(),"/h4JiItoInCond","_ISC_JiIto",".xlsx"))

# Interaction effect, I think I just have to test the interactions...
h4 = c("z_Ji_Ito_multi_perc  = 0", # Single
       "z_Ji_Ito_multi_perc + ConditionTriplet:z_Ji_Ito_multi_perc = 0", # Triplet
       "z_Ji_Ito_multi_perc + ConditionNonet:z_Ji_Ito_multi_perc = 0", # Nonet
       "z_Ji_Ito_multi_perc + ConditionComplete:z_Ji_Ito_multi_perc = 0") # complete
h4JiItoInCond = brms::hypothesis(brmSNV_ISCz_CxJiIto, h4)
write_xlsx(h4JiItoInCond$hypothesis, path = paste0(getwd(),"/h4JiItoInCond","_ISC_JiIto","undirected",".xlsx"))


# XYZ Coordinates ####
load("brmISCz_CxXYZ.RData")
summary(brmSNV_ISCz_CxXYZ)

# Interaction effect ???
h3 = c("ConditionTriplet:z_absCoord_x < 0",
       "ConditionTriplet:z_Coord_y  > 0",
       "ConditionTriplet:z_absCoord_x + ConditionTriplet:z_Coord_y + ConditionTriplet:z_absCoord_x:z_Coord_y < 0", 
       
       "ConditionNonet:z_absCoord_x - ConditionTriplet:z_absCoord_x< 0",
       "ConditionNonet:z_Coord_y- ConditionTriplet:z_Coord_y  > 0",
       "ConditionNonet:z_absCoord_x + ConditionNonet:z_Coord_y + ConditionNonet:z_absCoord_x:z_Coord_y - ConditionTriplet:z_absCoord_x - ConditionTriplet:z_Coord_y - ConditionTriplet:z_absCoord_x:z_Coord_y < 0",
       
       "ConditionComplete:z_absCoord_x - ConditionNonet:z_absCoord_x< 0",
       "ConditionComplete:z_Coord_y- ConditionNonet:z_Coord_y  > 0",
       "ConditionComplete:z_absCoord_x + ConditionComplete:z_Coord_y + ConditionComplete:z_absCoord_x:z_Coord_y - ConditionNonet:z_absCoord_x - ConditionNonet:z_Coord_y - ConditionNonet:z_absCoord_x:z_Coord_y < 0")
h3CoordInteraction = brms::hypothesis(brmSNV_ISCz_CxXYZ, h3)
write_xlsx(h3CoordInteraction$hypothesis, path = paste0(getwd(),"/h3CoordInteraction","_ISC_XYZ",".xlsx"))

# increases of ISC in different condtions
h4 =  c(# Single
        "z_absCoord_x <  0",
        "z_Coord_y > 0",
        "z_absCoord_x + z_Coord_y + z_absCoord_x:z_Coord_y < 0",
        # Triplet
        "z_absCoord_x + ConditionTriplet:z_absCoord_x < 0",
        "z_Coord_y  + ConditionTriplet:z_Coord_y  > 0",
        "z_absCoord_x + z_Coord_y + z_absCoord_x:z_Coord_y + ConditionTriplet:z_absCoord_x + ConditionTriplet:z_Coord_y + ConditionTriplet:z_absCoord_x:z_Coord_y < 0", 
        # Noent
        "z_absCoord_x + ConditionNonet:z_absCoord_x < 0",
        "z_Coord_y  + ConditionNonet:z_Coord_y  > 0",
        "z_absCoord_x + z_Coord_y + z_absCoord_x:z_Coord_y + ConditionNonet:z_absCoord_x + ConditionNonet:z_Coord_y + ConditionNonet:z_absCoord_x:z_Coord_y < 0", 
        # Complete
        "z_absCoord_x + ConditionComplete:z_absCoord_x < 0",
        "z_Coord_y  + ConditionComplete:z_Coord_y  > 0",
        "z_absCoord_x + z_Coord_y + z_absCoord_x:z_Coord_y + ConditionComplete:z_absCoord_x + ConditionComplete:z_Coord_y + ConditionComplete:z_absCoord_x:z_Coord_y < 0")
h4CoordInCond = brms::hypothesis(brmSNV_ISCz_CxXYZ, h4)
write_xlsx(h4CoordInCond$hypothesis, path = paste0(getwd(),"/h4CoordInCond","_ISC_XYZ",".xlsx"))

# increases of ISC in different condtions
h4 =  c(# Single
  "z_absCoord_x =  0",
  "z_Coord_y = 0",
  #      "z_Coord_z > 0", 
  "z_absCoord_x + z_Coord_y + z_absCoord_x:z_Coord_y = 0",
  # Triplet
  "z_absCoord_x + ConditionTriplet:z_absCoord_x = 0",
  "z_Coord_y  + ConditionTriplet:z_Coord_y  = 0",
  #      "z_Coord_z + ConditionTriplet:z_Coord_z > 0", 
  "z_absCoord_x + z_Coord_y + z_absCoord_x:z_Coord_y + ConditionTriplet:z_absCoord_x + ConditionTriplet:z_Coord_y + ConditionTriplet:z_absCoord_x:z_Coord_y = 0", 
  # Noent
  "z_absCoord_x + ConditionNonet:z_absCoord_x = 0",
  "z_Coord_y  + ConditionNonet:z_Coord_y  = 0",
  #      "z_Coord_z + ConditionNonet:z_Coord_z > 0", 
  "z_absCoord_x + z_Coord_y + z_absCoord_x:z_Coord_y + ConditionNonet:z_absCoord_x + ConditionNonet:z_Coord_y + ConditionNonet:z_absCoord_x:z_Coord_y = 0", 
  # Complete
  "z_absCoord_x + ConditionComplete:z_absCoord_x = 0",
  "z_Coord_y  + ConditionComplete:z_Coord_y  = 0",
  #       "z_Coord_z + ConditionComplete:z_Coord_z > 0", 
  "z_absCoord_x + z_Coord_y + z_absCoord_x:z_Coord_y + ConditionComplete:z_absCoord_x + ConditionComplete:z_Coord_y + ConditionComplete:z_absCoord_x:z_Coord_y = 0")
h4CoordInCond = brms::hypothesis(brmSNV_ISCz_CxXY, h4)
write_xlsx(h4CoordInCond$hypothesis, path = paste0(getwd(),"/h4CoordInCond","_ISCz_XY","undirected",".xlsx"))


# Graph Measures #####
load("brmISCz_CxGM.RData")
summary(brmSNV_ISCz_CxGM)
# Interaction effect, I think I just have to test the interactions...
h3 = c( # Degree
  "ConditionTriplet:z_Degree   > 0",
  "ConditionNonet:z_Degree  - ConditionTriplet:z_Degree  > 0",
  "ConditionComplete:z_Degree  - ConditionNonet:z_Degree  > 0",
  # Cls
  "ConditionTriplet:z_Closeness   > 0",
  "ConditionNonet:z_Closeness  - ConditionTriplet:z_Closeness  > 0",
  "ConditionComplete:z_Closeness  - ConditionNonet:z_Closeness  > 0",
  # Btw
  "ConditionTriplet:z_Betweenness   > 0",
  "ConditionNonet:z_Betweenness  - ConditionTriplet:z_Betweenness  > 0",
  "ConditionComplete:z_Betweenness  - ConditionNonet:z_Betweenness  > 0",
  #PtC
  "ConditionTriplet:z_PartCoef_DC   > 0",
  "ConditionNonet:z_PartCoef_DC  - ConditionTriplet:z_PartCoef_DC  > 0",
  "ConditionComplete:z_PartCoef_DC  - ConditionNonet:z_PartCoef_DC  > 0",
  #WMz
  "ConditionTriplet:z_WithinMod_DC   > 0",
  "ConditionNonet:z_WithinMod_DC  - ConditionTriplet:z_WithinMod_DC  > 0",
  "ConditionComplete:z_WithinMod_DC  - ConditionNonet:z_WithinMod_DC  > 0")
h3GMInteraction = brms::hypothesis(brmSNV_ISCz_CxGM, h3)
write_xlsx(h3GMInteraction$hypothesis, path = paste0(getwd(),"/h3GMInteraction","_ISC_GM",".xlsx"))

# Interaction effect, I think I just have to test the interactions...
h4 = c(# Single
  "z_Degree  > 0", 
  "z_Closeness  > 0", 
  "z_Betweenness  > 0", 
  "z_PartCoef_DC  > 0", 
  "z_WithinMod_DC  > 0",
  # Triplet
  "z_Degree + ConditionTriplet:z_Degree > 0",
  "z_Closeness + ConditionTriplet:z_Closeness > 0",
  "z_Betweenness + ConditionTriplet:z_Betweenness > 0",
  "z_PartCoef_DC + ConditionTriplet:z_PartCoef_DC > 0",
  "z_WithinMod_DC + ConditionTriplet:z_WithinMod_DC > 0",
  # Nonet
  "z_Degree + ConditionNonet:z_Degree > 0",
  "z_Closeness + ConditionNonet:z_Closeness > 0",
  "z_Betweenness + ConditionNonet:z_Betweenness > 0",
  "z_PartCoef_DC + ConditionNonet:z_PartCoef_DC > 0",
  "z_WithinMod_DC + ConditionNonet:z_WithinMod_DC > 0",
  # Complete
  "z_Degree + ConditionComplete:z_Degree > 0",
  "z_Closeness + ConditionComplete:z_Closeness > 0",
  "z_Betweenness + ConditionComplete:z_Betweenness > 0",
  "z_PartCoef_DC + ConditionComplete:z_PartCoef_DC > 0",
  "z_WithinMod_DC + ConditionComplete:z_WithinMod_DC > 0")
h4GMInCond = brms::hypothesis(brmSNV_ISCz_CxGM, h4)
write_xlsx(h4GMInCond$hypothesis, path = paste0(getwd(),"/h4GMInCond","_ISC_GM",".xlsx"))

# Interaction effect, I think I just have to test the interactions...
h4 = c(# Single
  "z_Degree  = 0", 
  "z_Closeness  = 0", 
  "z_Betweenness  = 0", 
  "z_PartCoef_DC  = 0", 
  "z_WithinMod_DC  = 0",
  # Triplet
  "z_Degree + ConditionTriplet:z_Degree = 0",
  "z_Closeness + ConditionTriplet:z_Closeness = 0",
  "z_Betweenness + ConditionTriplet:z_Betweenness = 0",
  "z_PartCoef_DC + ConditionTriplet:z_PartCoef_DC = 0",
  "z_WithinMod_DC + ConditionTriplet:z_WithinMod_DC = 0",
  # Nonet
  "z_Degree + ConditionNonet:z_Degree = 0",
  "z_Closeness + ConditionNonet:z_Closeness = 0",
  "z_Betweenness + ConditionNonet:z_Betweenness = 0",
  "z_PartCoef_DC + ConditionNonet:z_PartCoef_DC = 0",
  "z_WithinMod_DC + ConditionNonet:z_WithinMod_DC = 0",
  # Complete
  "z_Degree + ConditionComplete:z_Degree = 0",
  "z_Closeness + ConditionComplete:z_Closeness = 0",
  "z_Betweenness + ConditionComplete:z_Betweenness = 0",
  "z_PartCoef_DC + ConditionComplete:z_PartCoef_DC = 0",
  "z_WithinMod_DC + ConditionComplete:z_WithinMod_DC = 0")
h4GMInCond = brms::hypothesis(brmSNV_ISCz_CxGM, h4)
write_xlsx(h4GMInCond$hypothesis, path = paste0(getwd(),"/h4GMInCond","_ISC_GM","undirected",".xlsx"))


# Suface Based Morphometry: Cortical Thickness #### 
load("brmISCz_CxSBMt.RData")
summary(brmSNV_ISCz_CxSBMt)

# Interaction effect, I think I just have to test the interactions...
h3 = c("ConditionTriplet:z_Stats_Thick  > 0",
       "ConditionNonet:z_Stats_Thick - ConditionTriplet:z_Stats_Thick > 0",
       "ConditionComplete:z_Stats_Thick - ConditionNonet:z_Stats_Thick > 0")
h3SBMtInteraction = brms::hypothesis(brmSNV_ISCz_CxSBMt, h3)
write_xlsx(h3SBMtInteraction$hypothesis, path = paste0(getwd(),"/h3SBMtInteraction","_ISCz_SBMt",".xlsx"))

# Interaction effect, I think I just have to test the interactions...
h4 = c("z_Stats_Thick  > 0", # Single
       "z_Stats_Thick + ConditionTriplet:z_Stats_Thick > 0", # Triplet
       "z_Stats_Thick + ConditionNonet:z_Stats_Thick > 0", # Nonet
       "z_Stats_Thick + ConditionComplete:z_Stats_Thick > 0")
h4SBMtinCond = brms::hypothesis(brmSNV_ISCz_CxSBMt, h4)
write_xlsx(h4SBMtinCond$hypothesis, path = paste0(getwd(),"/h4SBMtinCond","_ISCz_SBMt",".xlsx"))

# Interaction effect, I think I just have to test the interactions...
h4 = c("z_Stats_Thick  = 0", # Single
       "z_Stats_Thick + ConditionTriplet:z_Stats_Thick = 0", # Triplet
       "z_Stats_Thick + ConditionNonet:z_Stats_Thick = 0", # Nonet
       "z_Stats_Thick + ConditionComplete:z_Stats_Thick = 0")
h4SBMtinCond = brms::hypothesis(brmSNV_ISCz_CxSBMt, h4)
write_xlsx(h4SBMtinCond$hypothesis, path = paste0(getwd(),"/h4SBMtinCond","_ISCz_SBMt","undirected",".xlsx"))

###### Test ASSUMPTIONS/DISTRIBUTION #############
# Libraries
library(brms)
library(bayesplot)
library(ggplot2)
library(RcppEigen)

outDir = ".../_Bayesian_Models_final_Z"
models = list.files(outDir, pattern=glob2rx("brmISCz*RData"))
modelNames = sapply(strsplit(models, ".", fixed = TRUE), "[[",1)
modellabels = substring(modelNames, 4)
setwd(outDir)

# adjust the global plotting theme
theme_set(theme_gray(base_size = 13, base_family = "Arial") + theme(panel.grid = element_blank()))
# define own functions for ppc_stat/pp_check
dispersion <- function(x) {var(x)/mean(x)}
modelNamesTrue = 1:5


for (i in modelNamesTrue) {
  # convert .RData -> .rdb/.rdx
  e = local({load(models[i]); environment()})
  tools:::makeLazyLoadDB(e, modelNames[i])
  lazyLoad(modelNames[i])
  # ls()
  eval(parse(text = paste0("fit = brmSNV_", modellabels[i])))
  
  assumDir = paste0(outDir, "/", modellabels[i])
  ifelse(!dir.exists(file.path(assumDir)), dir.create(file.path(assumDir)), FALSE)

  # general density overlay
  png(file = paste0(assumDir, "/Dens_", modellabels[i],".png"),width = 1200,height = 800)
  densPlot = pp_check(fit,
                      type = 'dens_overlay',
                      ndraws = 2000) + theme_bw(16)
  print(densPlot)
  dev.off()

  # grouped density overlay
  png(file = paste0(assumDir, "/DensGrouped_", modellabels[i],".png"),width = 1200,height = 800)
  densGrPlot = pp_check(fit,
                      type = 'dens_overlay_grouped', group = "Condition",
                      ndraws = 2000) + theme_bw(16)
  print(densGrPlot)
  dev.off()

  # mean
  png(file = paste0(assumDir, "/Mean_", modellabels[i],".png"),width = 1200,height = 800)
  meanPlot = pp_check(fit,
                      ndraws = 2000,
                      type ='stat',
                      stat = "mean",
                      binwidth = .00005)
  print(meanPlot)
  dev.off()

  # grouped mean
  png(file = paste0(assumDir, "/MeanGrouped_", modellabels[i],".png"),width = 1200,height = 800)
  meanGrPlot = pp_check(fit,
                      ndraws = 2000,
                      type ='stat_grouped', group = c("Condition"),
                      stat = "mean",
                      binwidth = .00005)
  print(meanGrPlot)
  dev.off()

  # mean & SD
  png(file = paste0(assumDir, "/M_SD_", modellabels[i],".png"),width = 1200,height = 800)
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
  png(file = paste0(assumDir, "/Disp_", modellabels[i],".png"),width = 1200,height = 800)
  dispPlot = ppc_stat(y = fit$data$ISC_FishZ, 
                      yrep = post_p,
                      stat="dispersion", binwidth = 0.0001)+
    theme_bw(16)
  print(dispPlot)
  dev.off()
  

  eval(parse(text = paste0("rm(brmSNV_", modellabels[i],",e)")))
}