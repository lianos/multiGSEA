library("multiGSEA")
library("testthat")
library("magrittr")
library("data.table")
library("dtplyr")
library("dplyr")

test_check("multiGSEA")

## Remove temporary files that were generated
test.dir <- system.file('tests', package='multiGSEA')
pdfs <- dir(test.dir, '\\.pdf$', full.names=TRUE)
if (length(pdfs)) {
  unlink(pdfs)
}
