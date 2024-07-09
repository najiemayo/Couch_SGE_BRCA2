#!/usr/local/biotools/r/R-4.2.2/bin/Rscript
## 
rm(list = ls())

library(rmarkdown)
RSTUDIO_PANDOC="/usr/local/biotools/pandoc/3.1.2/bin"
Sys.setenv(RSTUDIO_PANDOC="/usr/local/biotools/pandoc/3.1.2/bin")

args = commandArgs(T)

exons <- args[1]
loessSubset <- args[2]
EventCount_Filter <- args[3]
Ratio <- args[4] 
vset <- args[5]
fc <- args[6]
rmdnm <- args[7]
params <- args[8]
funR <- args[9]

source(params)
fn <- fns$fn[grep(exons, fns$exons)]
render_report2 = function(rmdnm, Gene, Exons, Input_file , Output_dir, loessSubset, EventCount_Filter, vset, Score_method, Lab_curation, Ratio, MM, Outlierremove, vexclude, clinvarfile, targetfile, HDRfile, extraFilter = "No", fn, funR) {
  #define main directory
  main_dir <- Output_dir
  
  #create directory if it doesn't exist
  if (!dir.exists(main_dir)) {
    dir.create(main_dir)
  }
  
  #define sub directory
  sub_dir <- "reports"
  
  #define full directory
  full_dir <- file.path(main_dir, sub_dir)
  
  #create directory if it doesn't exist
  if (!dir.exists(full_dir)) {
    dir.create(full_dir)
  }

  
  rmarkdown::render(
    rmdnm, 
    params = list(
      inputfile  = Input_file ,
      Output_dir = main_dir,
      loessSubset = loessSubset, 
      EventCount_Filter = EventCount_Filter,
      vset = vset,
      Score_method = Score_method,
      Lab_curation = Lab_curation,
      Exons = Exons,
      Ratio = Ratio,
      Gene = Gene,
      MM = MM,
      Outlierremove = Outlierremove,
      clinvarfile = clinvarfile,
      vexclude = vexclude ,
      targetfile = targetfile,
      extraFilter = extraFilter,
      HDRfile = HDRfile,
      funR = funR
    ),
    output_file = paste0( full_dir, "/", 
                          Gene, "_", vset, "_", loessSubset, "_", EventCount_Filter, "_FS", Score_method, 
                          "Lab_curation_", Lab_curation, "_E", fn,"_", Ratio,
                          "_MM_", MM,
                          "_Outlierremove_", Outlierremove,"_extraFilter_", extraFilter,".html")
  )

}

## function to generate inputfile 
createInput <- function(Gene, exonset, dtdir, Ratio,loessSubsetFiter,EventCount_Filter,  vset, fc, outdir){
  outfilename <- paste0(Gene, "_", vset, "_", gsub(", ", "", exonset), "_", Ratio, "_",loessSubsetFiter, "_",  EventCount_Filter,"_", fc)
  ration <- ifelse(Ratio == "D14vLib", "D14vLib", ifelse(Ratio == "D14vD5R", "ladjt_No_D14vD5", "ladjt_Yes_D14vD5"))
  filenames <- list.files(dtdir, pattern = glob2rx(paste0(Gene, "*",vset, "*","_", loessSubsetFiter, "_", 
                                                          EventCount_Filter, "_", ration, ".csv")))
  ##write(filenames, stderr())
  if(length(grep(",", exonset))>0){
    exonslist <- paste0("E", unlist(strsplit(exonset, ",")))
  }else{
    exonslist <- paste0("E", exonset) 
  }
  
  
  df <- data.frame(Gene = rep(Gene, length(exonslist)),
                   exon = exonslist,
                   dtdir = rep(dtdir, length(exonslist)),
                   filenames = rep(NA, length(exonslist)))
  ## fill in filenames in order of exonset to exonssset
  for(i in 1:length(exonslist)){
    #write(i, stderr())
    #write(filenames[grep(paste0(exonslist[i], "_SNV"), filenames)], stderr())
    df$filenames[i] <- filenames[grep(paste0(exonslist[i], "_SNV"), filenames)]
  }
  
  ## output file
  #create directory if it doesn't exist
  if (!dir.exists(outdir)) {
    dir.create(outdir)
  }
  write.csv(df, paste0(outdir, outfilename, ".csv"), row.names = FALSE)
  return(paste0(outdir, outfilename, ".csv"))
}


inputfile <- createInput(Gene, exons, dtdir, Ratio,loessSubset,EventCount_Filter,  vset, fc, outdir)
for(Score_method in Score_methods){
  for(Lab_curation in Lab_curations){
    for(MM in MMs){
      for(Outlierremove in Outlierremoves){
        for(exf in exfs){
          
          render_report2(
            rmdnm = rmdnm, 
            Input_file = inputfile,
            Output_dir = paste0(res_path, "/xexon_full/"),
            loessSubset = loessSubset,
            EventCount_Filter = EventCount_Filter,
            vset = vset,
            Score_method = Score_method,
            Lab_curation = Lab_curation,
            Exons = exons,
            Ratio  = Ratio,
            Gene = Gene,
            MM = MM,
            Outlierremove = Outlierremove,
            vexclude = vexclude,
            clinvarfile = clinvarfile,
            targetfile = targetfile,
            HDRfile = HDRfile,
            extraFilter = exf,
            fn=fn,
            funR = funR)
        }
      }
    }
  }
}
