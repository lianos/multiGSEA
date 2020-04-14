# Generate mouse and human entrez <-> ensembl gene conversion tables, primarily
# for us in `remap_identifiers`. These are stored in inst/extdata/identifiers.
#
# To help with the multimapping problem, we prefer mapping to protein_coding
# genes first, then the rest if a protein_coding doesn't exist

library(biomaRt)
library(dplyr)

factorize_biotype <- function(biotype) {
  bt.order <- c(
    'protein_coding', 'lincRNA',
    unique(biotype[grepl("lncrna", biotype, ignore.case = TRUE)]),
    'miRNA', 'rRNA', 'snoRNA',
    'scRNA', 'scaRNA', 'sRNA' ,"antisense", "ERCC")
  bt.order <- c(bt.order, setdiff(biotype, bt.order))
  bt.order <- intersect(bt.order, biotype)
  factor(biotype, bt.order)
}

# human ------------------------------------------------------------------------
hmart <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
hxref <- getBM(
  attributes = c("chromosome_name", "ensembl_gene_id", "entrezgene_id", "hgnc_symbol", "gene_biotype"),
  filters = "chromosome_name",
  values = c(1:22, "X", "Y", "MT"),
  mart = hmart)
table(hxref$chromosome_name)

hclean <- hxref %>%
  as.tbl() %>%
  dplyr::mutate(entrezgene_id = as.character(entrezgene_id)) %>%
  dplyr::filter(!is.na(ensembl_gene_id), nchar(ensembl_gene_id) > 0,
                !is.na(entrezgene_id), nchar(entrezgene_id) > 0) %>%
  dplyr::select(-chromosome_name) %>%
  group_by(ensembl_gene_id) %>%
  mutate(n_ensembl = n()) %>%
  group_by(entrezgene_id) %>%
  mutate(n_entrez = n()) %>%
  ungroup() %>%
  mutate(gene_biotype = factorize_biotype(gene_biotype))

hout.1 <- dplyr::filter(hclean, n_ensembl == 1, n_entrez == 1)
# for multi mappers, take the hit(s) that match the highest priority
# gene_biotype for the match
hout.multi <- hclean %>%
  dplyr::anti_join(hout.1, by = "ensembl_gene_id") %>%
  dplyr::arrange(ensembl_gene_id, gene_biotype) %>%
  dplyr::group_by(ensembl_gene_id) %>%
  dplyr::filter(gene_biotype == gene_biotype[1L]) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(entrezgene_id, gene_biotype) %>%
  dplyr::filter(gene_biotype == gene_biotype[1L]) %>%
  dplyr::ungroup()

hout <- bind_rows(hout.1, hout.multi)
# missing symbols
hmiss <- setdiff(hout$hgnc_symbol, hclean$hgnc_symbol)

hout <- rename(hout, symbol = "hgnc_symbol")
hfn <- "inst/extdata/identifiers/human-entrez-ensembl.csv"
write.csv(hout, hfn, row.names = FALSE)
system(paste("gzip", hfn))

# mouse ------------------------------------------------------------------------
mmart <- biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
mxref <- getBM(
  attributes = c("chromosome_name", "ensembl_gene_id", "entrezgene_id", "mgi_symbol", "gene_biotype"),
  filters = "chromosome_name",
  values = c(1:22, "X", "Y", "MT"),
  mart = mmart)
table(mxref$chromosome_name)

mclean <- mxref %>%
  as.tbl() %>%
  dplyr::mutate(entrezgene_id = as.character(entrezgene_id)) %>%
  dplyr::filter(!is.na(ensembl_gene_id), nchar(ensembl_gene_id) > 0,
                !is.na(entrezgene_id), nchar(entrezgene_id) > 0) %>%
  dplyr::select(-chromosome_name) %>%
  group_by(ensembl_gene_id) %>%
  mutate(n_ensembl = n()) %>%
  group_by(entrezgene_id) %>%
  mutate(n_entrez = n()) %>%
  ungroup() %>%
  mutate(gene_biotype = factorize_biotype(gene_biotype))

mout.1 <- dplyr::filter(mclean, n_ensembl == 1, n_entrez == 1)
# for multi mappers, take the hit(s) that match the highest priority
# gene_biotype for the match
mout.multi <- mclean %>%
  dplyr::anti_join(mout.1, by = "ensembl_gene_id") %>%
  dplyr::arrange(ensembl_gene_id, gene_biotype) %>%
  dplyr::group_by(ensembl_gene_id) %>%
  dplyr::filter(gene_biotype == gene_biotype[1L]) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(entrezgene_id, gene_biotype) %>%
  dplyr::filter(gene_biotype == gene_biotype[1L]) %>%
  dplyr::ungroup()

mout <- bind_rows(mout.1, mout.multi)

# missing symbols
mmiss <- setdiff(mout$mgi_symbol, mclean$mgi_symbol)

mout <- rename(mout, symbol = "mgi_symbol")
mfn <- "inst/extdata/identifiers/mouse-entrez-ensembl.csv"
write.csv(mout, mfn, row.names = FALSE)
system(paste("gzip", mfn))

