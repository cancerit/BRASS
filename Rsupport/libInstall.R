instLib = commandArgs(T)[1]

r = getOption("repos") # hard code the UK repo for CRAN
r["CRAN"] = "http://cran.uk.r-project.org"
options(repos = r)
rm(r)
source("http://bioconductor.org/biocLite.R")

ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    biocLite(new.pkg, ask=FALSE, lib=instLib)
  sapply(pkg, library, character.only = TRUE)
}

biocPackages <- c("data.table", "gam")
ipak(biocPackages)

install.packages("VGAM_1.0-3.tar.gz", type="source", lib=instLib)

biocPackages <- c("stringr", "poweRlaw", "zlibbioc", "RColorBrewer")
ipak(biocPackages)

install.packages("devtools", lib=instLib)
library(devtools)
options(download.file.method = "auto")
install_github("sb43/copynumber", ref="f1688edc154f1a0e3aacf7781090afe02882f623")
