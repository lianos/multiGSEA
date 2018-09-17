context("rename_rows")

test_that("rename_rows works from an embedded metadata gene frame ($genes)", {
  vm <- exampleExpressionSet()
  vm$genes$symbol[3] <- NA
  vms <- rename_rows(vm, "symbol")

  # Guarnatee same order of remapped object
  expect_equal(vms$E[,1], vm$E[,1], check.attributes = FALSE)

  is.symbol <- seq(nrow(vms)) != 3 & !duplicated(vm$genes$symbol)
  expect_equal(rownames(vms), ifelse(is.symbol, vm$genes$symbol, rownames(vm)))
})

test_that("rename_rows works from custom data.frame", {
  vm <- exampleExpressionSet()

  # A small remap data.frame
  remap <- subset(vm$genes, symbol %in% sample(symbol, 5))
  remap <- remap[, c("entrez_id", "symbol")]
  vms <- rename_rows(vm, remap)

  # Guarnatee same order of remapped object
  expect_equal(vms$E[,1], vm$E[,1], check.attributes = FALSE)

  # check that rows from remap were remapped OK
  remapped <- vms$genes$symbol %in% remap$symbol
  expect_equal(rownames(vms), ifelse(remapped, vm$genes$symbol, rownames(vm)))

  # Big remap data.frame
  vm.nona <- vm[!is.na(vm$genes$symbol),]
  vmsmall <- vm.nona[sample(nrow(vm.nona), 5),]
  remap2 <- vm.nona$genes[, c("entrez_id", "symbol")]
  vmss <- rename_rows(vmsmall, remap2)

  # Guarantee same order
  expect_equal(vmsmall$E[,1], vmss$E[,1], check.attributes = FALSE)
  expect_equal(rownames(vmss), vmsmall$genes$symbol)
})
