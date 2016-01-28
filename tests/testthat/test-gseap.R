## For some reason, the results from GSEA are not getting stored in
## y@result(!!)

if (FALSE) {
  context("clusterProfiler::GSEA")


  test_that("clusterProfiler::GSEA works as expedted", {
    ## See: http://guangchuangyu.github.io/2015/05/use-clusterprofiler-as-an-universal-enrichment-analysis-tool/
    library(DOSE)
    library(clusterProfiler)
    data(geneList)

    vm <- exampleExpressionSet(do.voom=TRUE)
    gsi <- exampleGeneSets(vm)
    gsl <- exampleGeneSets()
    gsd <- conform(GeneSetDb(gsl), vm)
    my <- multiGSEA(gsd, vm, vm$design, ncol(vm$design), NULL)
    lfc <- logFC(my)
    df <- as.data.frame(gsd)
    term2gene <- data.frame(term=paste(df$collection, df$name, sep=';;'),
                            gene=df$featureId, stringsAsFactors=FALSE)
    term2name <- term2gene[, 'term', drop=FALSE]
    term2name$name <- sub('.*?;;', '', term2name$term)

    ##anno <- clusterProfiler:::build_Anno(term2gene, term2name)
    geneList <- setNames(lfc$t, lfc$entrez_id)
    gsea <- clusterProfiler::GSEA(geneList, TERM2GENE=term2gene)
  })
}

if (FALSE) {
  ## See: http://guangchuangyu.github.io/2015/05/use-clusterprofiler-as-an-universal-enrichment-analysis-tool/
  require(DOSE)
  require(clusterProfiler)
  data(geneList)
  deg <- names(geneList)[abs(geneList)>2]
  ## downloaded from http://www.disgenet.org/ds/DisGeNET/results/all_gene_disease_associations.tar.gz
  gda <- read.delim(system.file('testdata', "all_gene_disease_associations.txt", package='multiGSEA'))

  disease2gene=gda[, c("diseaseId", "geneId")]
  disease2name=gda[, c("diseaseId", "diseaseName")]
  x = enricher(deg, TERM2GENE=disease2gene, TERM2NAME=disease2name)
  head(summary(x))

  ## For some reason, the results from GSEA are not getting stored in
  ## y@result(!!)
  y = GSEA(geneList, TERM2GENE=disease2gene, TERM2NAME=disease2name)
  head(summary(y))

  gseaplot(y, "umls:C0003872")
}
