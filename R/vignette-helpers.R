## This functions are used in the vignette and placed here to keep the vignette
## clean (enough)

##' Returns a DGEList from the airway dataset for downstream analysis.
##'
##' The airway package provides a dataset that can be used for differential
##' expression analysis. The object is a SummarizedExperiment with Ensembl
##' gene IDs as its row identifiers. This function will convert this object
##' to a DGEList with entrez id's as gene identifiers. Genes with dupliate
##' entrez ids (or NA ids) will be removed.
##'
##' @export
createAirwayDGEList <- function() {
  if (!exists('airway')) {
    # If this gets called multiple times, I get random
    # 'Error: C stack usage 7970760 is too close to the limit' errors
    data("airway", package="airway")
  }
  stopifnot(require('org.Hs.eg.db'))

  airway$dex <- relevel(airway$dex, 'untrt')
  y.all <- DGEList(
    assay(airway),
    group=airway$dex,
    samples=as.data.frame(colData(airway))[, c('cell', 'dex')])
  y.all <- calcNormFactors(y.all)

  ## Annotate genes
  y.all$genes <- data.frame(
    ensg_id=rownames(y.all),
    entrez_id=mapIds(org.Hs.eg.db, rownames(y.all), 'ENTREZID', 'ENSEMBL'),
    symbol=mapIds(org.Hs.eg.db, rownames(y.all), 'SYMBOL', 'ENSEMBL'),
    length=sum(width(rowRanges(airway))),
    stringsAsFactors=FALSE)

  ## remove NA entrez_ids
  out <- y.all[!is.na(y.all$genes$entrez_id),]
  dup.entrez <- out$genes$entrez_id[duplicated(out$genes$entrez_id)]
  out <- out[!out$genes$entrez_id %in% dup.entrez,,keep.lib.sizes=FALSE]
  calcNormFactors(out)
  rownames(out) <- out$genes$entrez_id
  out
}
