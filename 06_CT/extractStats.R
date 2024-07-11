analysis_dir = "..."

subjectNames= c("HFG_121")
templates = c("aparc", "lausanne120", "lausanne250","lausanne500","economo","BB50human")

firstLine = c(61,60,60,60,60,60)

for (s in 1:length(subjectNames)) {
  fsStatsDir = paste0(analysis_dir, subjectNames[s], '/FreeSurfer/stats')
  outpath = paste0(analysis_dir,subjectNames[s])
  setwd(fsStatsDir)
  for (t in 3)) {
    LHstats = read.delim(paste0('lh.',templates[t],'.stats'),header = TRUE, skip = firstLine[t]-1, sep = "")
    RHstats = read.delim(paste0('rh.',templates[t],'.stats'),header = TRUE, skip = firstLine[t]-1, sep = "")
    ColNames = names(LHstats)[3:length(LHstats)]
    
    LHstats = LHstats[,1:length(ColNames)]; names(LHstats) = ColNames
    LHstats$hemi = "left"
    LHstats$subj = subjectNames[s]
    LHstats$template = templates[t]
    
    RHstats = RHstats[,1:length(ColNames)]; names(RHstats) = ColNames    
    RHstats$hemi = "right"
    RHstats$subj = subjectNames[s]
    RHstats$template = templates[t]
    
    stats = rbind(LHstats,RHstats)
    
    write.csv(stats, file = paste0(outpath,"/Stats_",subjectNames[s],'_',templates[t],'.csv'))
  }}
