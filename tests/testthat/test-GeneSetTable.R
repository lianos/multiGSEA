context("GeneSetTable")

.c2 <- getMSigDBset('c2')
.c2$c2 <- head(.c2$c2[1:20])

.entrez.all <- unique(unlist(.c2))

test.x <- local({
  nrow <- length(.entrez.all) - 1000
  matrix(rnorm(nrow*4), nrow, 4, dimnames=list(sample(.entrez.all, nrow), NULL))
})

gst <- GeneSetTable(test.x, head(.c2, 20))
gst.50 <- GeneSetTable(test.x, head(.c2, 20), min.gs.size=50)
