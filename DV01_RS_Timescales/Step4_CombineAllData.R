# libraries
library(stringr)
library(dplyr)

subjectNames= c("HFG_121")
templates= c("lausanne250")
tempSizes = c(219)

analysis_dir = ".../Analysis/"

varNames = c("Subject","Template", "Region", # general Info
             "Timescales_Ito", "Timescales_Raut" , # Time scales
             "RC_Class_Sing","RC_Class_GroupPrev", "RC_Class_GroupFDR",  # Rich Club
             "DC_Class_Sing", "DC_Class_GroupPrev", "DC_Class_GroupFDR", # Diverse Club
             "SFC_S7_LeftAuditory", "SFC_S7_LeftSomatosensory", "SFC_S7_LeftVisual","SFC_S7_Multi_seed","SFC_S7_RightAuditory", "SFC_S7_RightSomatosensory", "SFC_S7_RightVisual", # SFC Step 7 individual
             "Lvl2_S7_LeftAuditory", "Lvl2_S7_LeftSomatosensory", "Lvl2_S7_LeftVisual","Lvl2_S7_Multi_seed","Lvl2_S7_RightAuditory", "Lvl2_S7_RightSomatosensory", "Lvl2_S7_RightVisual",  # SFC second Level Step 7 t-values
             "Network", "Network_percentage", "Ji_Ito", "Ji_Ito_adj", "Ji_Ito_multi_perc", # Overlap with Glasser Atlas
             "Coord_x","Coord_y","Coord_z", # position
             "Degree", "Closeness", "Betweenness","PartCoef_RC", "WithinMod_RC", # Grpah measures of centrality
             "PartCoef_DC", "WithinMod_DC", # Graph measures used for DC calculation
             "Stats_NumVert","Stats_SurfArea","Stats_GrayVol","Stats_ThickAvg","Stats_ThickStdv","Stats_MeanCurv","Stats_GausCurv","Stats_FoldInd","Stats_CurvInd", # Surface based morphometry
             # also cleared for outliers # DO I USE THAT????
             "Age","Sex","ICV") # more general information
SubInfo = readxl::read_excel(paste0(analysis_dir,"VP_Info.xlsx"))
brainInfo = read.csv(paste0(analysis_dir, 'CorticalMeasuresENIGMA_SurfAvg.csv'))

for (t in 1:length(templates)) {
  print(paste0("Reading for ", templates[t]))
  R2data = data.frame(matrix(data = NA, nrow = tempSizes[t]*length(subjectNames),length(varNames)))
  names(R2data) = varNames
  
  # add subject names
  R2data$Subject = sort(rep(subjectNames, tempSizes[t]))
  
  # add template name
  R2data$Template = templates[t]
  
  Timescales_Ito = read.delim(paste0(analysis_dir,"TimeScale_Ito/Murray_taus_regions_autocorrelation_rest100_",templates[t],".txt"),
                              header = FALSE, sep = ",")
    
  Timescales_Raut =read.delim(paste0(analysis_dir,"TimeScale_Raut/SubjectWise_hwhm2_custom_0.01-0.1_",templates[t],".txt"),
                              header = FALSE, sep = ",")
  
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
  
   Coords = read.csv(paste0(analysis_dir, subjectNames[i],"/Coord_", subjectNames[i],"_",templates[t],".csv"))
   Coords$Region = ifelse(Coords$hemi == "left", paste0("ctx-lh-",Coords$name),paste0("ctx-rh-",Coords$name))
   Coords = Coords[!(Coords$name == "unknown" | Coords$name == "corpuscallosum"),]
   
   CentRC = read.csv(paste0(analysis_dir,"2_RCAnalysis/",subjectNames[i],"/",templates[t],"/Ci_",templates[t],"_",subjectNames[i],".csv"))
   # also contains the information for DC analysis
   #CentDC = read.csv(paste0(analysis_dir,"2_RCAnalysis/",subjectNames[i],"/",templates[t],"/DCi_",templates[t],"_",subjectNames[i],".csv"))  
   
   StatsMorph = read.csv(paste0(analysis_dir, subjectNames[i],"/Stats_", subjectNames[i],"_",templates[t],".csv"))
   StatsMorph$Region = ifelse(StatsMorph$hemi == "left", paste0("ctx-lh-",StatsMorph$StructName),paste0("ctx-rh-",StatsMorph$StructName))
   
   # Regions
   R2data$Region[R2data$Subject == subjectNames[i]] = regions
   
   # Timescale Ito
   R2data$Timescales_Ito[R2data$Subject == subjectNames[i]] = Timescales_Ito[,i]
   
   # Timescale Raut
   R2data$Timescales_Raut[R2data$Subject == subjectNames[i]] = Timescales_Raut[,i]
   
   for (r in 1:length(regions)) {
     if (file.exists(paste0(analysis_dir,"2_RCAnalysis/",subjectNames[i],"/",templates[t],"/RC_Class_",templates[t],"_",subjectNames[i],".csv"))) {
      # RC single subject
      R2data$RC_Class_Sing[R2data$Subject == subjectNames[i] & R2data$Region == regions[r]] = RCSing$RC_class[RCSing$Region == regions[r]]
     }
     
     # RC group methods
     if (file.exists(paste0(analysis_dir,"2_RCAnalysis/prev/",templates[t],"/RC_Class_",templates[t],"_prev.csv"))) {
        R2data$RC_Class_GroupPrev[R2data$Subject == subjectNames[i] & R2data$Region == regions[r]] = RCprev$RC_class[RCprev$Region == regions[r]]
     }
    if (file.exists(paste0(analysis_dir,"2_RCAnalysis/tTestFDR/",templates[t],"/RC_Class_",templates[t],"_tTestFDR.csv"))) {
       R2data$RC_Class_GroupFDR[R2data$Subject == subjectNames[i] & R2data$Region == regions[r]] = RCFDR$RC_class[RCFDR$Region == regions[r]]
     }
     
     # DC single subject
     R2data$DC_Class_Sing[R2data$Subject == subjectNames[i] & R2data$Region == regions[r]] = DCSing$DC_class[DCSing$Region == regions[r]]
     # DC group methods
     R2data$DC_Class_GroupPrev[R2data$Subject == subjectNames[i] & R2data$Region == regions[r]] = DCprev$DC_class[DCprev$Region == regions[r]]
     R2data$DC_Class_GroupFDR[R2data$Subject == subjectNames[i] & R2data$Region == regions[r]] = DCFDR$DC_class[DCFDR$Region == regions[r]]
     
     # SFC Individual
     # SFCIndiv = read.csv(paste0(analysis_dir, subjectNames[i], "/SFC_analysis/", subjectNames[i], "_SFC_7_", templates[t],".csv"))
     SFCIndiv = read.csv(paste0(analysis_dir, subjectNames[i], "/SFC_analysis/", subjectNames[i], "_indiSFC_7_", templates[t],".csv"))
     SFCpos = c(1:length(varNames))[varNames == "SFC_S7_LeftAuditory"]
     R2data[R2data$Subject == subjectNames[i] & R2data$Region == regions[r], SFCpos:(SFCpos+6)] = SFCIndiv[SFCIndiv$Region == regions[r], 4:10]
     
     # SFC second Level t-Values
     SFCpos2 = c(1:length(varNames))[varNames == "Lvl2_S7_LeftAuditory"]
     R2data[R2data$Subject == subjectNames[i] & R2data$Region == regions[r], SFCpos2:(SFCpos2+6)] = SFC_Lvl2[SFC_Lvl2$Region == regions[r], 4:10]
     
     # Glasser Atlas
     ATLpos = c(1:length(varNames))[varNames == "Network"]
     R2data[R2data$Subject == subjectNames[i] & R2data$Region == regions[r], ATLpos:(ATLpos+4)] = Glasser_overlap[Glasser_overlap$Region == regions[r], 2:6]
     
     # Coordinates
     R2data$Coord_x[R2data$Subject == subjectNames[i] & R2data$Region == regions[r]] = Coords$x[Coords$Region == regions[r]]
     R2data$Coord_y[R2data$Subject == subjectNames[i] & R2data$Region == regions[r]] = Coords$y[Coords$Region == regions[r]]
     R2data$Coord_z[R2data$Subject == subjectNames[i] & R2data$Region == regions[r]] = Coords$z[Coords$Region == regions[r]]
     
     # Centrality measures
     degPos = c(1:length(varNames))[varNames == "Degree"]
     R2data[R2data$Subject == subjectNames[i] & R2data$Region == regions[r],degPos:(degPos+4)] = CentRC[Coords$Region == regions[r], c(3,4,5,8,6)]
     R2data[R2data$Subject == subjectNames[i] & R2data$Region == regions[r],(degPos+5):(degPos+6)] = CentRC[Coords$Region == regions[r], c(9,7)]
     
     # Surface Morpholometry
     surfPos = c(1:length(varNames))[varNames == "Stats_NumVert"]
     R2data[R2data$Subject == subjectNames[i] & R2data$Region == regions[r],surfPos:(surfPos+8)] = StatsMorph[StatsMorph$Region == regions[r], 3:11]
   }
   
   # subject Data
   R2data$Age[R2data$Subject == subjectNames[i]] = SubInfo$Alter[SubInfo$`VP-Code` == subjectNames[i]]
   R2data$Sex[R2data$Subject == subjectNames[i]] = SubInfo$Geschlecht[SubInfo$`VP-Code` == subjectNames[i]]
   
   # ICV
   R2data$ICV[R2data$Subject == subjectNames[i]] = brainInfo$ICV[brainInfo$SubjID == subjectNames[i]]
   
   print(paste0("Done ",as.character(i),"/",as.character(length(subjectNames))))
  }
  
  # make factor what needs to be factor
  R2data$Subject = factor(R2data$Subject)
  R2data$Region = factor(R2data$Region)
  R2data$RC_Class_Sing = factor(R2data$RC_Class_Sing, levels = c("RC_Local","RC_Feeder" ,"RC_Club"))
  R2data$RC_Class_GroupPrev = factor(R2data$RC_Class_GroupPrev, levels = c("RC_Local","RC_Feeder" ,"RC_Club"))
  R2data$RC_Class_GroupFDR = factor(R2data$RC_Class_GroupFDR, levels = c("RC_Local","RC_Feeder" ,"RC_Club"))
  R2data$DC_Class_Sing = factor(R2data$DC_Class_Sing, levels = c("DC_Local","DC_Feeder" ,"DC_Club"))
  R2data$DC_Class_GroupPrev = factor(R2data$DC_Class_GroupPrev, levels = c("DC_Local","DC_Feeder" ,"DC_Club"))
  R2data$DC_Class_GroupFDR = factor(R2data$DC_Class_GroupFDR, levels = c("DC_Local","DC_Feeder" ,"DC_Club"))
  R2data$Network = factor(R2data$Network)
  R2data$Ji_Ito = factor(R2data$Ji_Ito, levels = c("Unimodal", "Transmodal"))
  R2data$Ji_Ito_adj = factor(R2data$Ji_Ito_adj, levels = c("Periphery","NoLabel" ,"Core"))
  R2data$Sex = factor(R2data$Sex, levels = c("m","w"))
  
  # save data frame for this parcellation
  curDate = format(Sys.time(), "_%Y_%m_%d")
  save(R2data, file = paste0(analysis_dir, "R2data_",templates[t],curDate,".RData"))
}





