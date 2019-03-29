#' @import edgeR
#' @import limma
#' @import reshape2
#' @import data.table
#' @import methods
#' @import ggplot2
#' @import plotly
#' @import magrittr
#' @import checkmate
#' @importFrom circlize colorRamp2
#' @importFrom utils packageVersion getFromNamespace head tail write.csv
#' @importFrom stats setNames p.adjust model.matrix density filter phyper quantile runif
#' @importFrom stats as.dendrogram as.dist cor hclust order.dendrogram
#' @importFrom grDevices col2rgb rgb
#' @importFrom graphics abline axis legend lines rect segments hist pairs par
#' @importFrom graphics points strwidth
NULL

# A workaround to avoid R CMD check NOTEs for 'no visible binding for global variable'
# The following variables often appear within data.table[magic := stuff] when
# working with GSEA results.
utils::globalVariables(
  c(".",
    # appears in data.table[manip := ulations] of core multiGSEA tables
    "active", "collection", "featureId", "N", "n", "name",
    "mean.logFC.trim",
    "pval", "padj", "padj.by.collection", "padj.up", "padj.down",
    "value",  "significant", "x.idx",
    # random
    "group", "finalId", "ID", "id"
    ))
# https://github.com/jennybc/googlesheets/blob/master/R/googlesheets.R
