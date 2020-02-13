# multiGSEA 1.0.0 (Bioconductor release)

## New Features

* Adds preliminary support for data-frame based inputs to `multiGSEA()` to more
  easily support pre-ranked and ORA-like tests. When `x` is a data.frame, use
  the `rank_by` parameter to specify the column name that holds the rank metric,
  and `rank_order` to specify how to order the ranks (`"descending"`, 
  `"ascending"`, or `"ordered"`).
  
  The implementation is still tied to the original design of package which only
  supported a numeric ranking vector, and hijacks the whole `xmeta.` workaround.
  Consider using `FacileAnalysis::ffsea(data.frame, gdb)` if you want a more
  coherent interface.


## Breaking Changes

* `featureId` used all-over-the-place and `sample` columns from single sample
  scoring are renamed to `feature_id` and `sample_id` to play better with
  facile.bio

* `enrichtest` renamed to `ora` (over representation analysis) to be more inline
  with nomenclature in the field
  
# multiGSEA 0.10

## User Visible Changes

* MSigGeneSetDb updated to v6.1 and data refactored out into GeneSetDb.* packages
* `geneSet,MultiGSEAResult` method does not remove duplicated columns present
  in both `geneSet(gdb)` and `logFC`, but rather appends a `*.gs` suffix to
  the duplicated columns from `geneSet`. This happened because I defined
  a geneset with a `logFC` column (from the experiment the set was derived)
  from, but this collided with (and superseded) the `logFC` value for the gene
  returned from `logFC(mgresult)`, and messed up interpretation of a geneset's
  move across a contrast.
* Introduces idea of gene- and geneset-level metadata in a GeneSetDb. The
  former (gene-level metadata) was always possible by creating a GeneSetDb
  with a data.frame which has more columns than necessary, ie. instead of
  just having a collection,name,featureId data.frame, you could add a symbol
  column. This would simply be part of the GeneSetDb@db data.table. If there
  are columns in your input data.frame that are all the same value per
  collection,name subset, then these are marked as "geneset level annotations"
  and tacked onto the GeneSetDb@table data.table. I put this in to pave the
  way for simple support of "tags" per geneset, but there are other uses
  as well.

# multiGSEA 0.6.x

# User Visible Changes

* MSigDB version updated to 5.2 in 0.6.5
* Interactive plots can be generated via iplot
* single sample scoring methods broken out into more modular functions
* svdScore added. Executes Jason's SVD scoring trick and also calculates
  "all the usual" prcomp stuff.
* Adds a `[.GeneSetDb` method so we can subset the geneSets out of a GeneSetDb
  more easily. Would still like to have a fully operational subset(gdb, ...).

# New Features

* Adds `svdGeneSetTest` method. A gene x sample expression matrix is
  transformed into a geneset x sampel matrix, and "normal" differential
  expression expression testing is run over the geneset scores.

# multiGESA 0.4.x

## User Visible Changes

* scoreSingleSamples now defualts to returning a melted data.frame of results
  (previously we returned a (list of) matrix of scores). `melted=FALSE`
  parameter was changed to `as.matrix=FALSE`

* Internal MSigDB collections updated to v5.1 (as of v0.4.30).
  The major difference in v5.1 vs v5.0 seems to be an updated c7 (immunologic)
  collection, which comes on the heels of this paper:
    Compendium of Immune Signatures Identifies Conserved and Species-Specific
    Biology in Response to Inflammation
    http://www.cell.com/immunity/abstract/S1074-7613(15)00532-4

* All methods now default to returning data.frame(s) instead of data.table(s)
  data.table is still used internally, however I understand that some people
  don't want to have them returned into their workspace.
  A global option (multiGSEA.df.return) is available for you to tweak this
  behavior, ie. set `options(multiGSEA.df.return='data.tabe')` to have
  `geneSets(gdb)` return a data.table instead of a data.frame.

* Changed geneSetFeatureStatistics to geneSetsStats

## New Features

* `use.treat` parameter has been added to relevant places to run the internal
  differential testing via the "treat" framework made available in limma and
  edgeR via the `treat` and `glmTreat` "pipelines," respectively. By default,
  all pipelines do not use the treat framework (but perhaps they should).

* goseq, fry, and romer testing methods have been added.

* Support for DGEList expression input fully baked in. Previously, a DGEList
  was shot through voom to be used internally. Internal differential gene
  expression statistics on DGELists use edgeR's quasi-likelihood framework.
  See `?edgeR::glmQLFit` for more information, and particulary reference
  Aaron Lun's tutorial on using this approach here:

      It's DE-licious: a recipe for differential expression analyses of RNA-seq
      experiments using quasi-likelihood methods in edgeR
      http://www.statsci.org/smyth/pubs/QLedgeRPreprint.pdf

* Adds support for PANTHER pathway (and GOSLIM) genesets from PANTHER.db
  package via the `getPantherGeneSetDb` function.

* Added geneSetFeatureStats to return the logFC (and membership) information
  for a geneset after a multiGSEA run.

* Added constructors and `as()` coercion functions to create GeneSetDb's
  from GeneSetCollection(s) and data.frames (and vice versa).

* `scoreSingleSamples` is another wrapper function which calls several single
  sample geneset scoring algorithms, including the ones provided by the
  GSVA package (such as plage and ssGSEA)

* Adds expression-utils.R for CPM/RPKM and meltx. Also has methods to
  facilitate "universal access" to pData, fData, and expression data from
  the different often-used expression containers
  (ExpressionSet, EList, DGEList, SummarizedExperiment)

* MSigDB definitions were updated to v5.0, which adds a new "hallmark (h)"
  gene set. A call to `getMSigDBset` without a version parameter will return
  the v5.0 data.

  Note that MSigDB provides the entrez IDs for each geneset as human genes.
  The mouse version of these gene sets was constructing by mapping the
  human entrez ID's to their mouse orthologs via `igis::orthologs`. When
  this mapping returns multiple mouse IDs for a single human ID, all of
  the mouse IDs are kept in the mouse geneset.

  To get the previous v4.0 genesets (provided by WEHI), you would simply
  specify version='v4.0' in the getMSigDBSet like so:
  
      ```r
      gdb5 <- getMSigDBsetc(c('c2', 'c7'), species='mouse')
      gdb4 <- getMSigDBsetc(c('c2', 'c7'), species='mouse', version='v4.0')
      ````

## Bug Fixes

* do.hyperGeometricTest was failing when the expression object `x` was a
  matrix.

# multiGSEA 0.1.0

## New Features

* `multiGSEA` returns a `MultiGSEAResult` object which holds all of the things
  one would need to analyze/plot/poke at the results of the function call.

* Plotting functions added so user can plot results of arbitrary genesets.

# User Visible Changes

* This version is a reboots the internals of this package so that its design
  is cleaner. All of the reporting bits were removed (rmd.plugins handles that
  now), as well as significant additions to the API of this packages.  The
  changes are too numerous to report here. The previous version before the
  refactor have been put here for posterity:
  
   * http://resscm.gene.com/bioinfo/projects/R/tags/multiGSEA_pre-refactor
   * http://resscm.gene.com/bioinfo/projects/R/tags/rmd.plugins_pre-refactor

