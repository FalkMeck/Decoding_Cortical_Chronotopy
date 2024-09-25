# PREP DATA FOR MEAN SPINTEST #

analysis_dir = ".../"
setwd(analysis_dir)
load(paste0("...", "/DV01_RS_Timescales/TS_Raut_Data_final_lausanne250_2024_09_21.RData"))

# Timescales raw
meanTSRaut = tapply(R2data_RCRautFinal$Timescales_Raut,R2data_RCRautFinal$Region, mean)
data4spinTest = data.frame(levels(R2data_RCRautFinal$Region), meanTSRaut)
names(data4spinTest) = c("Region", "Timescale_Raut")
row.names(data4spinTest) = 1:219

# RC
data4spinTest$RC_perc = unlist(as.data.frame(tapply(R2data_RCRautFinal$RC_Class_Sing,R2data_RCRautFinal$Region, table))[ ,1])[(1:219)*3]/
  as.numeric(table(R2data_RCRautFinal$Region))

data4spinTest$RC_mode = NA 
for (i  in 1:219) {
  tmpData = R2data_RCRautFinal$RC_Class_Sing[R2data_RCRautFinal$Region == data4spinTest$Region[i]]
  tmpTable = table(tmpData)
  # if (length(names(tmpTable[tmpTable == max(tmpTable)]))>1) {
  #   print(i); print( data4spinTest$Region[i])
  # }
  data4spinTest$RC_mode[i] = names(tmpTable[tmpTable == max(tmpTable)])[1] # take lowest level
}

# uni-multi-modal: SFC, Ji-Ito-perc
data4spinTest$SFC = tapply(R2data_RCRautFinal$SFC_S7_1_Diff_Multi_seed,R2data_RCRautFinal$Region, mean)
data4spinTest$JiItoPerc = tapply(R2data_RCRautFinal$Ji_Ito_multi_perc,R2data_RCRautFinal$Region, mean)
  
# anterior-posterior and medial-lateral trend (+ z Coordiante): absX, Y, Z
data4spinTest$absX = tapply(R2data_RCRautFinal$absCoord_x,R2data_RCRautFinal$Region, mean)
data4spinTest$Y = tapply(R2data_RCRautFinal$Coord_y,R2data_RCRautFinal$Region, mean)
data4spinTest$Z = tapply(R2data_RCRautFinal$Coord_z,R2data_RCRautFinal$Region, mean)

# DC
data4spinTest$DC_perc = unlist(as.data.frame(tapply(R2data_RCRautFinal$DC_Class_Sing,R2data_RCRautFinal$Region, table))[ ,1])[(1:219)*3]/
  as.numeric(table(R2data_RCRautFinal$Region))

data4spinTest$DC_mode = NA 
for (i  in 1:219) {
  tmpData = R2data_RCRautFinal$DC_Class_Sing[R2data_RCRautFinal$Region == data4spinTest$Region[i]]
  tmpTable = table(tmpData)
  # if (length(names(tmpTable[tmpTable == max(tmpTable)]))>1) {
  #   print(i); print( data4spinTest$Region[i])
  # }
  data4spinTest$DC_mode[i] = names(tmpTable[tmpTable == max(tmpTable)])[1] # take lowest level
}

# Graph Measures
data4spinTest$Degree = tapply(R2data_RCRautFinal$Degree,R2data_RCRautFinal$Region, mean)
data4spinTest$Closeness = tapply(R2data_RCRautFinal$Closeness,R2data_RCRautFinal$Region, mean)
data4spinTest$Betweenness = tapply(R2data_RCRautFinal$Betweenness,R2data_RCRautFinal$Region, mean)
data4spinTest$PartCoef = tapply(R2data_RCRautFinal$PartCoef_DC,R2data_RCRautFinal$Region, mean)
data4spinTest$WithinModZ = tapply(R2data_RCRautFinal$WithinMod_DC,R2data_RCRautFinal$Region, mean)

# Thickness
data4spinTest$Thickness = tapply(R2data_RCRautFinal$Stats_ThickAvg,R2data_RCRautFinal$Region, mean)


load(paste0("...","/DV02_ISC/iscR2data_Final_2024_09_21.RData"))
sub_isc = iscR2data_Final[,c(1,3,4,6)]
wideISC = tidyr::pivot_wider(sub_isc, names_from = "Condition", values_from = "ISC_FishZ")
ISC_export = aggregate(Single ~ Region, data = wideISC, FUN = mean, na.action = na.omit)
ISC_export$ISC_FishZ.Triplet = aggregate(Triplet ~ Region, data = wideISC, FUN = mean, na.action = na.omit)$Triplet
ISC_export$ISC_FishZ.Nonet = aggregate(Nonet ~ Region, data = wideISC, FUN = mean, na.action = na.omit)$Nonet
ISC_export$ISC_FishZ.Complete = aggregate(Complete ~ Region, data = wideISC, FUN = mean, na.action = na.omit)$Complete
names(ISC_export)[2] = "ISC_FishZ.Single"

data4spinTest = cbind(data4spinTest, ISC_export[,2:5])

write.csv(x = data4spinTest, file = "data4spinTest.csv")


      