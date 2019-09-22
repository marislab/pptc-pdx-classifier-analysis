if (!require("dplyr")){
  install.packages("dplyr", repos='http://cran.us.r-project.org', dependencies = TRUE)
  library(dplyr)
}
if (!require("ggplot2")){
  install.packages("ggplot2", repos='http://cran.us.r-project.org', dependencies = TRUE)
  library(ggplot2)
}
if (!require("tidyr")){
  install.packages("tidyr", repos='http://cran.us.r-project.org', dependencies = TRUE)
  library(ggplot2)
}
if (!require("EnvStats")){
  install.packages("EnvStats", repos='http://cran.us.r-project.org', dependencies = TRUE)
  library(EnvStats)
}
if (!require("reshape2")){
  install.packages("reshape2", repos='http://cran.us.r-project.org', dependencies = TRUE)
  library(psych)
}
if (!require("ggbeeswarm")){
  install.packages("ggbeeswarm", repos='http://cran.us.r-project.org', dependencies = TRUE)
  library(gplots)
}
if (!require("devtools")){
  install.packages("https://cran.r-project.org/src/contrib/Archive/devtools/devtools_2.0.2.tar.gz", repos=NULL, dependencies = TRUE, type = "source")
  library(devtools)
}
if (!require("ggpubr")){
  devtools::install_github("kassambara/ggpubr")
  library(ggpubr)
}
