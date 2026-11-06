native_fixture <- function() {
  simulate_spatial_multiomics(heterospa_config(n_cells = 40, n_domains = 3,
    n_features = c(RNA = 12, ATAC = 15, Protein = 5), shared_dim = 2,
    specific_dim = 2, n_spatial_basis = 4, seed = 421))
}

native_runtime <- function() {
  skip_if_not_installed("reticulate")
  skip_if_not(reticulate::py_module_available("anndata"), "Optional Python anndata unavailable")
  reticulate::import("anndata", convert = FALSE)
}

expect_native_assay <- function(a, view, design) {
  has <- function(mapping, key) isTRUE(reticulate::py_to_r(mapping$`__contains__`(key)))
  expect_true(has(a$obsm, "spatial"))
  expect_true(has(a$layers, "counts"))
  expect_identical(unlist(reticulate::py_to_r(a$obs_names$to_list())), view$ids)
  expect_identical(unlist(reticulate::py_to_r(a$var_names$to_list())), colnames(view$counts))
  expect_equal(as.matrix(reticulate::py_to_r(a$X)), as.matrix(view$counts), ignore_attr = TRUE)
  expect_equal(as.matrix(reticulate::py_to_r(a$layers$`__getitem__`("counts"))), as.matrix(view$counts), ignore_attr = TRUE)
  expect_equal(reticulate::py_to_r(a$obsm$`__getitem__`("spatial")), unname(view$coordinates), ignore_attr = TRUE)
  expect_false(any(c("domain", "cell_type", "original_row") %in% unlist(reticulate::py_to_r(a$obs$columns$to_list()))))
  expect_equal(reticulate::py_to_r(a$uns$`__getitem__`("heterospa_design")), design)
  expect_false(reticulate::py_to_r(a$uns$`__getitem__`("include_truth")))
}

test_that("native Python H5AD actually persists spatial, counts and isolated metadata", {
  ad <- native_runtime(); x <- native_fixture(); before <- serialize(x, NULL)
  for (design in c("V1", "V2")) {
    dest <- tempfile("native_h5ad_"); paths <- export_h5ad(x, dest, design = design, backend = "python")
    expect_length(paths, 3)
    for (m in names(paths)) expect_native_assay(ad$read_h5ad(paths[[m]]), HeteroSpaSim:::.public_view(x, m, design), design)
    expect_true(file.exists(file.path(dest, "evaluator_truth", "RNA_unit_map.csv")))
    expect_false(file.exists(file.path(dest, "public_inputs", "RNA_unit_map.csv")))
  }
  # Python initialization can change the serialization header's native locale.
  # Compare the deserialized scientific object rather than header bytes.
  expect_identical(x, unserialize(before))
})

test_that("native H5MU retains every modality field and private maps stay in sidecar", {
  native_runtime()
  skip_if_not(reticulate::py_module_available("mudata"), "Optional Python mudata unavailable")
  md <- reticulate::import("mudata", convert = FALSE); x <- native_fixture()
  f <- tempfile(fileext = ".h5mu"); export_h5mu(x, f, design = "V2", backend = "python")
  obj <- md$read_h5mu(f)
  for (m in names(x$assays)) expect_native_assay(obj$mod$`__getitem__`(m), HeteroSpaSim:::.public_view(x, m, "V2"), "V2")
  expect_true(file.exists(file.path(paste0(f, "_evaluator_truth"), "RNA_unit_map.csv")))
})
