## Run this script from an R session with the current working directory set
## to here
library(EnrichmentBrowser)
odir <- 'pways-2015-07-20'
dir.create(odir)
pwys <- download.kegg.pathways("hsa", odir)
