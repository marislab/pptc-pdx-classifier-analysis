###### Calculate TMB and create plots for Figures 4 and S4########
#
#     Authors: Jo Lynne Rokita, Gregory P. Way
#     Updated 2019-09-18
################################################################

# working directory (created with git clone)
mainDir <- "~/pptc-pdx-classifier-analysis/"
dataDir <- paste0(mainDir,"data/")
# set path to your git cloned repo
script.folder <- paste0(mainDir, "R/") 
##create directory for output files
ifelse(!dir.exists(file.path(paste0(mainDir, "results/"))), dir.create(file.path(paste0(mainDir, "results/"))), 
       "Directory exists!")

results.folder <- paste0(mainDir, "results/") 

##set wd
setwd(mainDir)

####Install package Dependencies if needed
source(paste0(script.folder, "install.packages.R"))
####Source plot theme script
source(paste0(script.folder, "theme.R"))
####Load histology color function
source(paste0(script.folder, "mutation-color-function.R"))

###run analyses
source(paste0(script.folder, "classifier-plots-revision.R"))
source(paste0(script.folder, "corplots.R"))