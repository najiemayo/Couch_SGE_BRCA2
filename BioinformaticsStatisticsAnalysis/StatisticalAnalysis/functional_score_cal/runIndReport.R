#!/usr/local/biotools/r/R-4.2.2/bin/Rscript

rm(list = ls())

library(rmarkdown)
RSTUDIO_PANDOC="/usr/local/biotools/pandoc/3.1.2/bin"
Sys.setenv(RSTUDIO_PANDOC="/usr/local/biotools/pandoc/3.1.2/bin")

args <- commandArgs(trailingOnly = TRUE)


## this R script need to be run with argument (which is the parameter settings)

loessSubset <- args[1]
EventCount_Filter <- args[2]
ladjt <- args[3] 
vset <- args[4]
rmdnm <- args[5]
params <- args[6]
funR <- args[7]
rmdpath <- args[8]

source(params)

Inds <- read.table(datalist, sep = "\t", header = T)
Es <- Inds$Exon[which(Inds$Gene == Gene)]
Es <- substring(Es, 2)
Es <- unique(Es)

## where the rmd file reside

setwd(rmdpath) ## where the rmd file reside
render_report = function(rmdnm, Gene, Exon, ETfile, datalist, targetfile, HDRfile, vexclude, Output_dir, loessSubset, EventCount_Filter, vset, C_fc, ladjt ) {

  #create directory if it doesn't exist
  if (!dir.exists(Output_dir)) {
    dir.create(Output_dir)
  }
  
  rmarkdown::render(
    rmdnm, params = list(
      Exon = Exon, 
      ETfile = ETfile,
      datalist = datalist,
      targetfile = targetfile,
      HDRfile = HDRfile,
      vexclude = vexclude,
      Output_dir = Output_dir, 
      loessSubset = loessSubset, 
      EventCount_Filter = EventCount_Filter,
      vset = vset,
      Gene = Gene,
      C_fc = C_fc,
      ladjt = ladjt
    ),
    output_file = paste0(Output_dir, Gene, "_E", Exon, "_",  vset, "_", loessSubset, "_", EventCount_Filter,"_LadjtD14vD5_",  ladjt,".html")
  )
}


for(Exon in Es){
  rmdnmi <- paste0(Exon,"_", rmdnm)
  system(paste0("cp " , rmdpath, "/", rmdnm, " ", rmdpath,"/", rmdnmi))
  
  render_report(rmdnm = paste0(rmdpath, "/", rmdnmi),
                Exon= Exon, 
                ETfile = ETfile,
                datalist = datalist,
                targetfile = targetfile,
                HDRfile = HDRfile,
                vexclude = vexclude,
                Output_dir = dtdir, 
                loessSubset = loessSubset, 
                EventCount_Filter = EventCount_Filter,
                vset = vset,
                Gene = Gene,
                C_fc = 1,
                ladjt = ladjt) 
  
}


  
