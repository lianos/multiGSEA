## use the imultiGSEA package for now
# on.rescomp <- dir.exists('/gne')
# if (on.rescomp) {
#   library(switchr)
#   switchrBaseDir("/gne/research/web/prd/shiny/ptd/shinyapps/.switchr")
#   switchrNoUnload(TRUE)
#   switchTo("facile")
# }

## Load multiGSEA? You shouldn't have to because this should be invoked by
## `multiGSEA::explore()`, but who knows how the user got here
# library(rprojroot)
# root <- find_root(is_r_package)
# devtools::load_all(root)
# library(multiGSEA)
devtools::load_all('~/workspace/Rpkgs/GNE/multiGSEA')
library(DT)

## Loading "standard" Libraries ------------------------------------------------
# devtools::load_all('~/workspace/Rpkgs/GNE/multiGSEA.shiny')
# library(multiGSEA.shiny)

library(shiny)
library(shinydashboard)
library(shinyjs)
# library(rbokeh)
library(plotly)
library(data.table)
library(dplyr)
library(dtplyr)

theme_set(theme_bw())

## By default shiny limits upload size to 5 MB, let's change this to 30MB
## (which is kind of big, no?)
options(shiny.maxRequestSize=30*1024^2)
options(multiGSEA.df.return='data.table')
