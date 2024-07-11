library(data.table)

analysis_dir = ".../"

subjectNames= c("HFG_121")
templates = c("lausanne250")

for (s in 1:length(subjectNames)) {
  outpath = paste0(analysis_dir, subjectNames[s])
  for (t in 1:length(templates)) {
    setwd(paste0(analysis_dir, subjectNames[s],'/FreeSurfer/',templates[t]))
    
    LabelLH = read.csv(paste0(analysis_dir,'LabelsLH_',templates[t],'.csv')) 
    LabelRH = read.csv(paste0(analysis_dir,'LabelsRH_',templates[t],'.csv')) 
    
    atlasPat_lh = dir(pattern = paste0(templates[t],"-lh"))
    atlasPat_rh = dir(pattern = paste0(templates[t],"-rh"))
    
    if (length(LabelLH$Label) != length(atlasPat_lh)) {
      LabelLH_true = LabelLH$Label[!(LabelLH$Label == "unknown" | LabelLH$Label == "corpuscallosum")]
    } else{LabelLH_true = LabelLH$Label}
    
    if (length(LabelRH$Label) != length(atlasPat_rh)) {
      LabelRH_true = LabelRH$Label[!(LabelRH$Label == "unknown" | LabelRH$Label == "corpuscallosum")]
    }else{LabelRH_true = LabelRH$Label}
    
    if (t == 5) {LabelLH_true = c("unknown", LabelLH_true); LabelRH_true = c("unknown", LabelRH_true)}
    

        # lh
       atlaslh_labels <- vector('list', length=length(atlasPat_lh))
        for (i in seq_along(atlasPat_lh)) {
          atlaslh_labels[[i]] =  fread(atlasPat_lh[i])
          atlaslh_labels[[i]][, name := LabelLH_true[i]]
        }
        atlaslh_labels <- rbindlist(atlaslh_labels, fill = TRUE)
        atlaslh_mean = atlaslh_labels[, .(x=mean(V2), y=mean(V3), z=mean(V4)), by=name]
        atlaslh_mean$hemi = "left"
        atlaslh_mean$subj = subjectNames[s]
        atlaslh_mean$template = templates[t]
        
        # rh
        atlasrh_labels <- vector('list', length=length(atlasPat_rh))
        for (i in seq_along(atlasPat_rh)) {
          atlasrh_labels[[i]] =  fread(atlasPat_rh[i])
          atlasrh_labels[[i]][, name := LabelRH_true[i]]
        }
        atlasrh_labels <- rbindlist(atlasrh_labels, fill = TRUE)
        atlasrh_mean = atlasrh_labels[, .(x=mean(V2), y=mean(V3), z=mean(V4)), by=name]
        atlasrh_mean$hemi = "right"
        atlasrh_mean$subj = subjectNames[s]
        atlasrh_mean$template = templates[t]
        
        atlas_coor = rbind(atlaslh_mean, atlasrh_mean)
        write.csv(atlas_coor, file = paste0(outpath,"/Coord_",subjectNames[s],"_",templates[t],".csv"))
  }}  
