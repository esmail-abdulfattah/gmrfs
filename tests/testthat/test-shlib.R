test_that("INLA's library is looked up by the model's name, not by 'name'", {
  skip_if_not_installed("INLA")
  ## inla.external.lib() takes its argument's name (as library()); the old
  ## call inla.external.lib(name) looked for .../external/name/libname.so
  p <- do.call(INLA::inla.external.lib, list("fbesag"))
  expect_match(p, "external/fbesag/libfbesag\\.so$")
  expect_null(.gmrfs_inla_lib("no_such_model"))
})

test_that("INLA's library is used only if it has the model's function", {
  ## a copy of the package's own library stands in for one shipped by INLA:
  ## it has inla_cgeneric_gmrfs_fbesag, but not INLA's inla_cgeneric_fbesag
  so <- file.path(tempdir(), paste0("fake_inla_lib", .Platform$dynlib.ext))
  expect_true(file.copy(gmrfs_own_shlib(), so, overwrite = TRUE))
  on.exit(unlink(so))
  expect_true(.gmrfs_has_symbol(so, "inla_cgeneric_gmrfs_fbesag"))
  expect_false(.gmrfs_has_symbol(so, "inla_cgeneric_fbesag"))
  expect_false(.gmrfs_has_symbol(file.path(tempdir(), "missing.so"), "x"))

  own <- gmrfs_own_shlib()
  local_mocked_bindings(.gmrfs_inla_lib = function(name) so)
  expect_equal(gmrfs_resolve_shlib("fbesag", symbol = "inla_cgeneric_gmrfs_fbesag"), so)
  expect_equal(gmrfs_resolve_shlib("fbesag", symbol = "inla_cgeneric_fbesag"), own)
  expect_equal(gmrfs_resolve_shlib("fbesag", prefer_inla = FALSE,
                                   symbol = "inla_cgeneric_gmrfs_fbesag"), own)
  ## loading and unloading the copy leaves the package's own code usable
  cm <- cblocks(3, 1, list(diag(3)))
  expect_equal(as.matrix(gmrf_Q(cm, 0)), diag(3), ignore_attr = TRUE)
})

test_that("without a library from INLA, the models use the package's own", {
  skip_if_not_installed("INLA")
  W <- matrix(0, 4, 4)
  W[cbind(1:3, 2:4)] <- 1
  W <- W + t(W)
  expect_equal(fbesag(W, id = c(1, 1, 2, 2))$f$cgeneric$shlib, gmrfs_own_shlib())
})
