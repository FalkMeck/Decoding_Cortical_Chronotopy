# libraries
library(stringr)
library(dplyr)

subjectNames = c("HFG_121")


templates= c("lausanne250") # focus on one first
tempSizes = c(219)

conditions = c("Single","Triplet","Nonet","Complete")

analysis_dir = ".../Analysis/"
iscDir = ".../ISC_Analysis/"
setwd(iscDir)

varNames = c("Subject","Template", "Region", # general Info
             "ISC_FishZ", "NumVox" , "Condition", # ISC Info
             "RC_Class_Sing","RC_Class_GroupPrev", "RC_Class_GroupFDR",  # Rich Club
             "DC_Class_Sing", "DC_Class_GroupPrev", "DC_Class_GroupFDR", # Diverse Club
             "SFC_S7_LeftAuditory", "SFC_S7_LeftSomatosensory", "SFC_S7_LeftVisual","SFC_S7_Multi_seed","SFC_S7_RightAuditory", "SFC_S7_RightSomatosensory", "SFC_S7_RightVisual", # SFC Step 7 individual
             "Lvl2_S7_LeftAuditory", "Lvl2_S7_LeftSomatosensory", "Lvl2_S7_LeftVisual","Lvl2_S7_Multi_seed","Lvl2_S7_RightAuditory", "Lvl2_S7_RightSomatosensory", "Lvl2_S7_RightVisual",  # SFC second Level Step 7 t-values
             "SFC_S7_1_Diff_Multi_seed", # EDIT (2024_09_21)
             "Network", "Network_percentage", "Ji_Ito", "Ji_Ito_adj", "Ji_Ito_multi_perc", # Overlap with Glasser Atlas
             "Coord_x","Coord_y","Coord_z", # position
             "Degree", "Closeness", "Betweenness","PartCoef_RC", "WithinMod_RC", # Grpah measures of centrality
             "PartCoef_DC", "WithinMod_DC", # Graph measures used for DC calculation
             "Stats_NumVert","Stats_SurfArea","Stats_GrayVol","Stats_ThickAvg","Stats_ThickStdv","Stats_MeanCurv","Stats_GausCurv","Stats_FoldInd","Stats_CurvInd", # Surface based morphometry
             # also cleared for outliers # DO I USE THAT????
             "Age","Sex","ICV","Beamer") # more general information
SubInfo = readxl::read_excel(paste0(analysis_dir,"VP_Info.xlsx"))
brainInfo = read.csv(paste0(analysis_dir, 'CorticalMeasuresENIGMA_SurfAvg.csv'))

beamerInfo = read.csv2(paste0(iscDir,"SubjectCodes_and_Beamer.csv"), header = FALSE)
names(beamerInfo) = c("SubNum","Subject","Beamer","StartBlock","List")
beamerInfo$Beamer = factor(beamerInfo$Beamer, levels = c("1","0"), labels = c("new","old"))

for (t in 1:length(templates)) {
  print(paste0("Reading for ", templates[t]))
  iscR2data = data.frame(matrix(data = NA, nrow = tempSizes[t]*length(subjectNames)*length(conditions),length(varNames)))
  names(iscR2data) = varNames
  
  # add subject names
  iscR2data$Subject = sort(rep(subjectNames, tempSizes[t]*length(conditions)))

  #add Conditions
  iscR2data$Condition =  rep(sort(rep(conditions, tempSizes[t])), length(subjectNames))
    
  # add template name
  iscR2data$Template = templates[t]
  
  if (file.exists(paste0(analysis_dir,"2_RCAnalysis/prev/",templates[t],"/RC_Class_",templates[t],"_prev.csv"))) {
    RCprev =  read.csv(paste0(analysis_dir,"2_RCAnalysis/prev/",templates[t],"/RC_Class_",templates[t],"_prev.csv"))
  }
  if (file.exists(paste0(analysis_dir,"2_RCAnalysis/tTestFDR/",templates[t],"/RC_Class_",templates[t],"_tTestFDR.csv"))) {
    RCFDR = read.csv(paste0(analysis_dir,"2_RCAnalysis/tTestFDR/",templates[t],"/RC_Class_",templates[t],"_tTestFDR.csv"))
  }
  
  DCprev = read.csv(paste0(analysis_dir,"2_RCAnalysis/prev/",templates[t],"/DC_Class_",templates[t],"_prev.csv"))
  DCFDR = read.csv(paste0(analysis_dir,"2_RCAnalysis/tTestFDR/",templates[t],"/DC_Class_",templates[t],"_tTestFDR.csv"))
  
  regions = RCprev$Region
  
  # Load 2nd Level SFC maps
  SFC_Lvl2 = read.csv(paste0(analysis_dir, "_SFC_Second_Level/Second_Level_tValues_Step_7_", templates[t], ".csv"))
  # Load atlas/Network overlap/classification based on Glasser atlas
  Glasser_overlap = readxl::read_xlsx(paste0(analysis_dir, "overview_", templates[t], "_Glasser.xlsx"))
  
  for (i in 1:length(subjectNames)) {
    print(paste0("Processing for ", subjectNames[i]))
    if (file.exists(paste0(analysis_dir,"2_RCAnalysis/",subjectNames[i],"/",templates[t],"/RC_Class_",templates[t],"_",subjectNames[i],".csv"))) {
      RCSing = read.csv(paste0(analysis_dir,"2_RCAnalysis/",subjectNames[i],"/",templates[t],"/RC_Class_",templates[t],"_",subjectNames[i],".csv"))
    }   
    
    DCSing = read.csv(paste0(analysis_dir,"2_RCAnalysis/",subjectNames[i],"/",templates[t],"/DC_Class_",templates[t],"_",subjectNames[i],".csv"))
    
    # Load evel SFC maps Differce (2024_09_21)
    SFC_Diff = read.csv(paste0(analysis_dir,subjectNames[i],'/SFC_analysis/', subjectNames[i],'_indiSFC_diff_7_',templates[t],'.csv'))
    
    Coords = read.csv(paste0(analysis_dir, subjectNames[i],"/Coord_", subjectNames[i],"_",templates[t],".csv"))
    Coords$Region = ifelse(Coords$hemi == "left", paste0("ctx-lh-",Coords$name),paste0("ctx-rh-",Coords$name))
    Coords = Coords[!(Coords$name == "unknown" | Coords$name == "corpuscallosum"),]
    
    CentRC = read.csv(paste0(analysis_dir,"2_RCAnalysis/",subjectNames[i],"/",templates[t],"/Ci_",templates[t],"_",subjectNames[i],".csv"))
    # also contains the information for DC analysis
    #CentDC = read.csv(paste0(analysis_dir,"2_RCAnalysis/",subjectNames[i],"/",templates[t],"/DCi_",templates[t],"_",subjectNames[i],".csv"))  
    
    StatsMorph = read.csv(paste0(analysis_dir, subjectNames[i],"/Stats_", subjectNames[i],"_",templates[t],".csv"))
    StatsMorph$Region = ifelse(StatsMorph$hemi == "left", paste0("ctx-lh-",StatsMorph$StructName),paste0("ctx-rh-",StatsMorph$StructName))
    
    for (c in 1:length(conditions)) {
        # Regions
        iscR2data$Region[iscR2data$Subject == subjectNames[i] & iscR2data$Condition == conditions[c]] = regions
          
        iscFile = read.csv(paste0(iscDir, subjectNames[i],"/", subjectNames[i], "_ISC_",conditions[c],".csv"))
    
    for (r in 1:length(regions)) {
      
      # ISC
      iscR2data$ISC_FishZ[iscR2data$Subject == subjectNames[i] & iscR2data$Condition == conditions[c] & iscR2data$Region == regions[r]] =
        iscFile$ISC[iscFile$Region == regions[r]]
      iscR2data$NumVox[iscR2data$Subject == subjectNames[i] & iscR2data$Condition == conditions[c] & iscR2data$Region == regions[r]] =
        iscFile$VoxelNum[iscFile$Region == regions[r]]
      
      
      if (file.exists(paste0(analysis_dir,"2_RCAnalysis/",subjectNames[i],"/",templates[t],"/RC_Class_",templates[t],"_",subjectNames[i],".csv"))) {
        # RC single subject
        iscR2data$RC_Class_Sing[iscR2data$Subject == subjectNames[i] & iscR2data$Condition == conditions[c] & iscR2data$Region == regions[r]] = 
          RCSing$RC_class[RCSing$Region == regions[r]]
      }
      
      # RC group methods
      if (file.exists(paste0(analysis_dir,"2_RCAnalysis/prev/",templates[t],"/RC_Class_",templates[t],"_prev.csv"))) {
        iscR2data$RC_Class_GroupPrev[iscR2data$Subject == subjectNames[i] & iscR2data$Condition == conditions[c] & iscR2data$Region == regions[r]] = 
          RCprev$RC_class[RCprev$Region == regions[r]]
      }
      if (file.exists(paste0(analysis_dir,"2_RCAnalysis/tTestFDR/",templates[t],"/RC_Class_",templates[t],"_tTestFDR.csv"))) {
        iscR2data$RC_Class_GroupFDR[iscR2data$Subject == subjectNames[i] & iscR2data$Condition == conditions[c] & iscR2data$Region == regions[r]] = 
          RCFDR$RC_class[RCFDR$Region == regions[r]]
      }
      
      # DC single subject
      iscR2data$DC_Class_Sing[iscR2data$Subject == subjectNames[i] & iscR2data$Condition == conditions[c] & iscR2data$Region == regions[r]] = 
        DCSing$DC_class[DCSing$Region == regions[r]]
      # DC group methods
      iscR2data$DC_Class_GroupPrev[iscR2data$Subject == subjectNames[i] & iscR2data$Condition == conditions[c] & iscR2data$Region == regions[r]] = 
        DCprev$DC_class[DCprev$Region == regions[r]]
      iscR2data$DC_Class_GroupFDR[iscR2data$Subject == subjectNames[i] & iscR2data$Condition == conditions[c] & iscR2data$Region == regions[r]] = 
        DCFDR$DC_class[DCFDR$Region == regions[r]]
      
      # SFC Individual
      if (file.exists(paste0(analysis_dir, subjectNames[i], "/SFC_analysis/", subjectNames[i], "_indiSFC_7_", templates[t],".csv"))){
          SFCIndiv = read.csv(paste0(analysis_dir, subjectNames[i], "/SFC_analysis/", subjectNames[i], "_indiSFC_7_", templates[t],".csv"))
          SFCpos = c(1:length(varNames))[varNames == "SFC_S7_LeftAuditory"]
          iscR2data[iscR2data$Subject == subjectNames[i] & iscR2data$Condition == conditions[c] & iscR2data$Region == regions[r], SFCpos:(SFCpos+6)] = 
            SFCIndiv[SFCIndiv$Region == regions[r], 4:10]
      }
      # SFC sencond Level t-Values
      SFCpos2 = c(1:length(varNames))[varNames == "Lvl2_S7_LeftAuditory"]
      iscR2data[iscR2data$Subject == subjectNames[i] & iscR2data$Condition == conditions[c] & iscR2data$Region == regions[r], SFCpos2:(SFCpos2+6)] = 
        SFC_Lvl2[SFC_Lvl2$Region == regions[r], 4:10]
      
      # EDIT: SFC Differences
      SFCposDiff = c(1:length(varNames))[varNames == "SFC_S7_1_Diff_Multi_seed"]
      iscR2data[iscR2data$Subject == subjectNames[i] & iscR2data$Condition == conditions[c] & iscR2data$Region == regions[r], SFCposDiff] = 
        SFC_Diff[SFC_Diff$Region == regions[r], 4]
      
      # Glasser Atlas
      ATLpos = c(1:length(varNames))[varNames == "Network"]
      iscR2data[iscR2data$Subject == subjectNames[i] & iscR2data$Condition == conditions[c] & iscR2data$Region == regions[r], ATLpos:(ATLpos+4)] = 
        Glasser_overlap[Glasser_overlap$Region == regions[r], 2:6]
      
      # Coordinates
      iscR2data$Coord_x[iscR2data$Subject == subjectNames[i] & iscR2data$Condition == conditions[c] & iscR2data$Region == regions[r]] = Coords$x[Coords$Region == regions[r]]
      iscR2data$Coord_y[iscR2data$Subject == subjectNames[i] & iscR2data$Condition == conditions[c] & iscR2data$Region == regions[r]] = Coords$y[Coords$Region == regions[r]]
      iscR2data$Coord_z[iscR2data$Subject == subjectNames[i] & iscR2data$Condition == conditions[c] & iscR2data$Region == regions[r]] = Coords$z[Coords$Region == regions[r]]
      
      # Centrality measures
      degPos = c(1:length(varNames))[varNames == "Degree"]
      iscR2data[iscR2data$Subject == subjectNames[i] & iscR2data$Condition == conditions[c] & iscR2data$Region == regions[r],degPos:(degPos+4)] = CentRC[Coords$Region == regions[r], c(3,4,5,8,6)]
      iscR2data[iscR2data$Subject == subjectNames[i] & iscR2data$Condition == conditions[c] & iscR2data$Region == regions[r],(degPos+5):(degPos+6)] = CentRC[Coords$Region == regions[r], c(9,7)]
      
      # Surface Morpholometry
      surfPos = c(1:length(varNames))[varNames == "Stats_NumVert"]
      iscR2data[iscR2data$Subject == subjectNames[i] & iscR2data$Condition == conditions[c] & iscR2data$Region == regions[r],surfPos:(surfPos+8)] = StatsMorph[StatsMorph$Region == regions[r], 3:11]
    }}
    
    # subject Data
    iscR2data$Age[iscR2data$Subject == subjectNames[i]] = SubInfo$Alter[SubInfo$`VP-Code` == subjectNames[i]]
    iscR2data$Sex[iscR2data$Subject == subjectNames[i]] = SubInfo$Geschlecht[SubInfo$`VP-Code` == subjectNames[i]]
    
    # ICV
    iscR2data$ICV[iscR2data$Subject == subjectNames[i]] = brainInfo$ICV[brainInfo$SubjID == subjectNames[i]]
    
    # Beamer
    iscR2data$Beamer[iscR2data$Subject == subjectNames[i]] = beamerInfo$Beamer[beamerInfo$Subject == subjectNames[i]]
    
    print(paste0("Done ",as.character(i),"/",as.character(length(subjectNames))))
  }
  
  # make factor what needs to be factor
  iscR2data$Subject = factor(iscR2data$Subject)
  iscR2data$Region = factor(iscR2data$Region)
  iscR2data$RC_Class_Sing = factor(iscR2data$RC_Class_Sing, levels = c("RC_Local","RC_Feeder" ,"RC_Club"))
  iscR2data$RC_Class_GroupPrev = factor(iscR2data$RC_Class_GroupPrev, levels = c("RC_Local","RC_Feeder" ,"RC_Club"))
  iscR2data$RC_Class_GroupFDR = factor(iscR2data$RC_Class_GroupFDR, levels = c("RC_Local","RC_Feeder" ,"RC_Club"))
  iscR2data$DC_Class_Sing = factor(iscR2data$DC_Class_Sing, levels = c("DC_Local","DC_Feeder" ,"DC_Club"))
  iscR2data$DC_Class_GroupPrev = factor(iscR2data$DC_Class_GroupPrev, levels = c("DC_Local","DC_Feeder" ,"DC_Club"))
  iscR2data$DC_Class_GroupFDR = factor(iscR2data$DC_Class_GroupFDR, levels = c("DC_Local","DC_Feeder" ,"DC_Club"))
  iscR2data$Network = factor(iscR2data$Network)
  iscR2data$Ji_Ito = factor(iscR2data$Ji_Ito, levels = c("Unimodal", "Transmodal"))
  iscR2data$Ji_Ito_adj = factor(iscR2data$Ji_Ito_adj, levels = c("Periphery","NoLabel" ,"Core"))
  iscR2data$Sex = factor(iscR2data$Sex, levels = c("m","w"))
  iscR2data$Beamer = factor(iscR2data$Beamer, levels = c("1","0"), labels = c("new", "old"))
  
  # save data frame for this parcellation
  curDate = format(Sys.time(), "_%Y_%m_%d")
  save(iscR2data, file = paste0(iscDir, "iscR2data_",templates[t],curDate,".RData")) # }

# t = 1
curDate = format(Sys.time(), "_%Y_%m_%d")
load(paste0(iscDir, "iscR2data_",templates[t],curDate,".RData"))

  # convert Volume to liter
  iscR2data$ICV_l = iscR2data$ICV/100^3
  
  iscR2data$absCoord_x = abs(iscR2data$Coord_x)
  iscR2data$absCoord_y = abs(iscR2data$Coord_y)
  
  # Scale continous variables 
  iscR2data_scaled = iscR2data
  nums = unlist(lapply(iscR2data, is.numeric), use.names = FALSE)  
  toScale = (c(13:length(iscR2data))[nums[13:length(iscR2data)]])
  iscR2data_scaled = cbind(iscR2data_scaled, data.frame(matrix(NA,nrow = length(iscR2data_scaled$Subject), ncol = length(toScale))))
  names(iscR2data_scaled)[(length(iscR2data)+1):length(iscR2data_scaled)] = paste0("z_",names(iscR2data[toScale]))
  iscR2data_scaled[(length(iscR2data)+1):length(iscR2data_scaled)] = as.data.frame(lapply(iscR2data_scaled[,toScale], scale))
  
  iscR2data_scaled$Stats_Surf = iscR2data_scaled$Stats_SurfArea/iscR2data_scaled$ICV
  iscR2data_scaled$Stats_Thick = iscR2data_scaled$Stats_ThickAvg/iscR2data_scaled$ICV
  
  iscR2data_scaled$z_Stats_Surf = as.numeric(scale(iscR2data_scaled$Stats_Surf))
  iscR2data_scaled$z_Stats_Thick = as.numeric(scale(iscR2data_scaled$Stats_Thick))

  iscR2data_scaled$Hemi = "R"
  iscR2data_scaled$Hemi[grepl("-lh-", iscR2data_scaled$Region,fixed = TRUE)] = "L"
  iscR2data_scaled$Hemi = factor(iscR2data_scaled$Hemi, levels = c("L", "R"))
  
  
  save(iscR2data_scaled, file = paste0(iscDir, "iscR2data_scaled_",templates[t],curDate,".RData"))

  

}

# libraries
library(plot3D)
library(ggplot2)
library(lattice)
library(ggpubr)
library(moments)

setwd(".../DV02_ISC/")

template = 'lausanne250'
curDate = format(Sys.time(), "_%Y_%m_%d")
#curDate = "_2024_03_13"
load(paste0(getwd(), "/iscR2data_scaled_",template,curDate,".RData"))

iscR2data_scaled$Condition = factor(iscR2data_scaled$Condition, levels = c("Single", "Triplet", "Nonet", "Complete"))
iscR2data_scaled$z_NumVox = as.numeric(scale(iscR2data_scaled$NumVox))

ggplot(iscR2data_scaled, aes(x = ISC_FishZ)) +
  geom_histogram(fill = "white", colour = "black", bins = 100) +
  facet_grid(Condition ~ .)

moments::skewness(iscR2data_scaled$ISC_FishZ)
moments::kurtosis(iscR2data_scaled$ISC_FishZ)
moments::jarque.test(iscR2data_scaled$ISC_FishZ)
# significant, not NV
outliers::outlier(iscR2data_scaled$ISC_FishZ)

qqnorm(iscR2data_scaled$ISC_FishZ)

# check fro other distribtions and maybe glmer
fitdistrplus::descdist(iscR2data_scaled$ISC_FishZ, boot = 100L)
rcompanion::plotNormalHistogram(iscR2data_scaled$ISC_FishZ, breaks = 100)


# Kurtosis is to large, outliers!
boxISC = boxplot(iscR2data_scaled$ISC_FishZ)
iscR2data_Final = iscR2data_scaled[!(iscR2data_scaled$ISC_FishZ > (boxISC$stats[4] + 1.5*IQR(iscR2data_scaled$ISC_FishZ)) | 
                                       iscR2data_scaled$ISC_FishZ < (boxISC$stats[2] - 1.5*IQR(iscR2data_scaled$ISC_FishZ))),]
table(iscR2data_Final$Subject)/876
rcompanion::plotNormalHistogram(iscR2data_Final$ISC_FishZ, breaks = 100)
fitdistrplus::descdist(iscR2data_Final$ISC_FishZ, boot = 100L)

ggdensity(iscR2data_Final, x = "ISC_FishZ", fill = "Subject", title = "ISC_FishZ") +
  scale_x_continuous(limits = c(-0.2, 1.6)) +
  stat_overlay_normal_density(color = "red", linetype = "dashed")

table(iscR2data_Final$Subject, iscR2data_Final$Condition)/219

## 
hist(iscR2data_Final$ISC_FishZ, breaks = 100)
boxplot(iscR2data_Final$ISC_FishZ)

moments::skewness(iscR2data_Final$ISC_FishZ) # 0.2343574
moments::kurtosis(iscR2data_Final$ISC_FishZ) #2.899748
moments::jarque.test(iscR2data_Final$ISC_FishZ)
# Hair et al. (2010) and Bryne (2010) argued that data is considered to be normal if skewness is between ‐2 to +2 and kurtosis is between ‐7 to +7.

iscR2data_Final$Beamer[is.na(iscR2data_Final$Beamer)] = "old"
save('iscR2data_Final', file = paste0(getwd(), "/iscR2data_Final", curDate , ".RData"))
load(paste0(getwd(), "/iscR2data_Final", curDate , ".RData"))