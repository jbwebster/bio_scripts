


#if (!requireNamespace("BiocManager", quietly=TRUE))
#  install.packages("BiocManager")
#BiocManager::install("Rsubread")

library(Rsubread)

simReads(
  #transcript database and relative expression
  "simulationdata/inputs/TMPRSS2.ERG.all.fasta", #transcript file
  c(1, 1, 1, 1), # expression
  
  "simulationdata/outputs/TMPRSS2.ERG.all.outs/simulation1.v1", #output
  
  library.size = 10000,
  read.length = 100,
  
  paired.end = TRUE,
  
  # simulating sequencing errors
  simulate.sequencing.error = TRUE
  
  
)

