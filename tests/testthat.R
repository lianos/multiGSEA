library("multiGSEA")
library("testthat")

test_check("multiGSEA")

## Remove temporary files that were generated
test.dir <- system.file('tests', package='multiGSEA')
pdfs <- dir(test.dir, '\\.pdf$', full.names=TRUE)
if (length(pdfs)) {
  unlink(pdfs)
}
