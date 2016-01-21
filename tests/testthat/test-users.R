context("Errors submitted by users")


if (FALSE) {
  test_that("Russ's error is his fault", {
    ## This error was due to a design matrix that was not full rank
    datnames <- load(system.file('testdata', 'messForSteve.Rdata', package='multiGSEA'))
    gsd <- conform(gsd.all, vm)
    mg <- multiGSEA(gsd, vm, vm$design[, -12], methods='camera')
  })

}
