## library(multiGSEA)
devtools::load_all('~/workspace/Rpkgs/GNE/multiGSEA')
library(shiny)
library(shinydashboard)
library(DT)
library(ggplot2)
library(plotly)
library(data.table)
library(dplyr)
library(dtplyr)
theme_set(theme_bw())
## By default shiny limits upload size to 5 MB, let's change this to 20MB
options(shiny.maxRequestSize=29*1024^2)
options(multiGSEA.df.return='data.table')

source('utils.R')

xmg <- readRDS('~/workspace/projects/rutz/BAP1/reports/johnnycache/multiGSEA-NGS429-joint.rds')
xcol <- 'h'
xname <- 'HALLMARK_E2F_TARGETS'
