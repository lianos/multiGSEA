
<!-- README.md is generated from README.Rmd. Please edit that file -->

# multiGSEA

[![Build
Status](https://travis-ci.org/lianos/multiGSEA.svg?branch=develop)](https://travis-ci.org/lianos/multiGSEA)
[![Coverage
status](https://codecov.io/gh/lianos/multiGSEA/branch/master/graph/badge.svg)](https://codecov.io/github/lianos/multiGSEA?branch=master)

The `multiGSEA` package was built to facilitate the use of gene sets in
the analysis of high throughput genomics data (primarily RNA-seq).
Analysts can orchestrate any number of GSEA methods across a specific
contrast using the unified interface provided by the `multiGSEA`
function, and a shiny application is provided that facilitates the
exploration and interpretation of GSEA results.

  - The `multiGSEA` function is a wrapper that orchestrates the
    execution of any number of user-specified gene set enrichment
    analyses (GSEA) over a particular experimental contrast of interest.
    This will create a `MultiGSEAResult` object which stores the results
    of each GSEA method internally, allowing for easy query and
    retrieval.
  - A sister
    [`multiGSEA.shiny`](https://github.com/lianos/multiGSEA.shiny)
    package provides an `explore` function, which is invoked on
    `MultiGSEAREsult` objects returned from a call to `multiGSEA`. The
    shiny application facilitates interactive exploration of these GSEA
    results. This application can also be deployed to a shiny server and
    can be initialized by uploading a serialized `MultiGSEAResult`
    `*.rds` file.

Full details that outline the use of this software package is provided
in the [package’s vignette](vignettes/multiGSEA.Rmd), however a brief
description is outlined below.

# Usage

A [thorough vignette](vignettes/multiGSEA.Rmd) is provided with this
package that illustrates the multitudes of its gene-set-based analysis
functionality. For convenience, a small taster of the package’s
functionality is included here.

A subset of the RNA-seq data tumor/normal samples in the BRCA indication
from the TCGA are provided in this package. We will use that data to
perform a “camera” and “fry” gene set enrichment analysis of tumor vs
normal samples using the MSigDB hallmark and c2 gene set collections
with `multiGSEA`.

``` r
library(multiGSEA)
library(dplyr)
gdb <- getMSigGeneSetDb(c('h', 'c2'), 'human')
vm <- exampleExpressionSet(dataset='tumor-vs-normal', do.voom=TRUE)
mg <- multiGSEA(gdb, vm, vm$design, "tumor", methods=c("camera", "fry"))
```

We can view the top “camera” results with the smallest pvalues like so:

``` r
results(mg, "camera") %>% 
  arrange(pval) %>% 
  select(collection, name, padj) %>% 
  head
#>   collection                                        name         padj
#> 1         c2      SOTIRIOU_BREAST_CANCER_GRADE_1_VS_3_UP 2.604378e-38
#> 2         c2 ROSTY_CERVICAL_CANCER_PROLIFERATION_CLUSTER 3.601729e-37
#> 3         c2         NAKAYAMA_SOFT_TISSUE_TUMORS_PCA2_DN 1.004916e-28
#> 4         c2              KANG_DOXORUBICIN_RESISTANCE_UP 1.478337e-22
#> 5         c2                     BENPORATH_PROLIFERATION 2.437669e-22
#> 6         c2               CROONQUIST_IL6_DEPRIVATION_DN 1.934138e-21
```

The shift in expression of the genes within the top gene set can be
visualized with the `iplot` function below. This plot produces
interactive graphics, but rasterized versions are saved for use with
this `README`
file:

``` r
iplot(mg, 'c2', 'SOTIRIOU_BREAST_CANCER_GRADE_1_VS_3_UP', type="density")
```

<img src="vignettes/images/README_iplot_density.png" />

``` r
iplot(mg, 'c2', 'SOTIRIOU_BREAST_CANCER_GRADE_1_VS_3_UP', type="boxplot")
```

<img src="vignettes/images/README_iplot_boxplot.png" />

When these plots are rendered in your workspace or an Rmarkdown
document, the user can hover of the genes (dots) to see their name and
differential expression statistics.

# Installation

The multiGSEA suite of packages will soon be submitted to bioconductor
and installable via the recommended `biocLite` mechanism. In the
meantime, these packages can be installed like so:

``` r
# install.packages("devtools")
devtools::install_github("lianos/multiGSEA")
devtools::install_github("lianos/multiGSEA.shiny")
devtools::install_github("lianos/GeneSetDb.MSigDB.Hsapiens.v61")
devtools::install_github("lianos/GeneSetDb.MSigDB.Mmusculus.v61")
```
