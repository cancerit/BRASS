instLib = commandArgs(T)[1]

r = getOption("repos") # hard code the UK repo for CRAN
r["CRAN"] = "http://cran.uk.r-project.org"
options(repos = r)
rm(r)
source("http://bioconductor.org/biocLite.R")
biocLite("data.table", ask=FALSE, lib=instLib)
biocLite("gam", ask=FALSE, lib=instLib)
biocLite("VGAM", ask=FALSE, lib=instLib)
biocLite("stringr", ask=FALSE, lib=instLib)
biocLite("poweRlaw", ask=FALSE, lib=instLib)
biocLite("zlibbioc", ask=FALSE, lib=instLib)
biocLite("RColorBrewer", ask=FALSE, lib=instLib)
biocLite("copynumber", ask=FALSE, lib=instLib)
