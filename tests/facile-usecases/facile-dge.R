library(FacileData)
library(FacileAnalysis)
library(dplyr)

tfds <- FacileDataSet("~/workspace/data/FacileData/dockerlink/FacileTcgaDataSet")

sstats <- samples(tfds) %>%
  with_sample_covariates() %>%
  group_by(indication) %>%
  summarize(n = n(),
            normal = sum(sample_type == "normal"),
            tumor = sum(sample_type == "tumor"))

tsamples <- tfds %>%
  filter_samples(indication == "KICH") %>%
  with_sample_covariates()

dge.tn <- tsamples %>%
  fdge_model_def("sample_type", "tumor", "normal", "sex") %>%
  fdge(method = "voom")

# GSEA Stuff ===================================================================
library(multiGSEA)
# devtools::load_all(".")
gdb <- getMSigGeneSetDb(c("h"), "human", "ensembl")
x <- dge.tn

# Prep the ffsea.data.frame ....................................................

xdf <- local({
  rank_by <- "logFC"
  signed <- TRUE
  max_padj <- 0.10
  min_logFC <- 1

  ranks. <- tidy(ranks(x, signed = signed))
  ranks. <- mutate(ranks.,
                   selected = padj <= max_padj, abs(logFC) >= min_logFC,
                   direction = ifelse(logFC > 0, "up", "down"))
  checkmate::assert_choice(rank_by, colnames(ranks.))

  take.cols <- c(
    rank_by, "selected", "direction",
    "symbol", "meta", "logFC", "t", "B",
    "AveExpr", "pval", "padj", "CI.L", "CI.R", "effective_length")
  take.cols <- intersect(take.cols, colnames(ranks.))

  select(ranks., feature_id, {{take.cols}})
})

# Exercise ffsea ...............................................................
dfres <- ffsea(xdf, gdb, methods = "cameraPR",
               rank_by = "logFC", rank_order = "ranked",
               select_by = "significant")

ttres <- ffsea(x, gdb, methods = "cameraPR")

dfenr <- ffsea(xdf, gdb, methods = c("cameraPR", "goseq"),
               rank_by = "logFC", rank_order = "ranked",
               select_by = "selected", split.updown = FALSE)
