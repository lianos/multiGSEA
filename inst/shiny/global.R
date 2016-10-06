## use the imultiGSEA package for now
# on.rescomp <- dir.exists('/gne')
# if (on.rescomp) {
#   library(switchr)
#   switchrBaseDir("/gne/research/web/prd/shiny/ptd/shinyapps/.switchr")
#   switchrNoUnload(TRUE)
#   switchTo("facile")
# }

## Loading Custom Libraries ----------------------------------------------------
## Let's make sure to load custom-deployed packages first so that we don't
## accidentally load older ones that live deeper-down our .libPaths()
devtools::load_all('~/workspace/Rpkgs/GNE/multiGSEA')
library(DT)

## Loading "standard" Libraries ------------------------------------------------
# devtools::load_all('~/workspace/Rpkgs/GNE/multiGSEA.shiny')
# library(multiGSEA.shiny)

library(shiny)
library(shinydashboard)
library(rbokeh)
library(data.table)
library(dplyr)
library(dtplyr)

## By default shiny limits upload size to 5 MB, let's change this to 30MB
## (which is kind of big, no?)
options(shiny.maxRequestSize=30*1024^2)
options(multiGSEA.df.return='data.table')
