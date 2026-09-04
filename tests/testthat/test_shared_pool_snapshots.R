source(file.path(core_root, "R", "shared_pool_snapshots.R"), local = FALSE)

test_that("snapshot rows preserve joint values and requested column order", {
  path <- tempfile(fileext = ".csv")
  on.exit(unlink(path), add = TRUE)
  writeLines(c("# synthetic CmdStan header", "", "a,b,c", "1,2,3",
               "# adaptation ended", "4,5,6", "", "# elapsed time"), path)
  expect_identical(
    .sab_shared_pool_read_joint_row(path, c("c", "a"), 2L, 2L),
    c(c = 6, a = 4)
  )
})

test_that("snapshot parser rejects absent and duplicated CSV columns", {
  path <- tempfile(fileext = ".csv")
  on.exit(unlink(path), add = TRUE)
  writeLines(c("a,b", "1,2"), path)
  expect_error(.sab_shared_pool_read_joint_row(path, "c", 1L, 1L),
               "header is incompatible")
  writeLines(c("a,a", "1,2"), path)
  expect_error(.sab_shared_pool_read_joint_row(path, "a", 1L, 1L),
               "header is incompatible")
})

test_that("snapshot parser verifies total rows and selected row existence", {
  path <- tempfile(fileext = ".csv")
  on.exit(unlink(path), add = TRUE)
  writeLines(c("a,b", "1,2", "3,4"), path)
  expect_error(.sab_shared_pool_read_joint_row(path, "a", 1L, 3L),
               "retained row/count is invalid")
  expect_error(.sab_shared_pool_read_joint_row(path, "a", 3L, 2L),
               "retained row/count is invalid")
  writeLines(c("# comments only", ""), path)
  expect_error(.sab_shared_pool_read_joint_row(path, "a", 1L, 0L),
               "retained row/count is invalid")
})

test_that("snapshot parser checks malformed rows after the selected row", {
  path <- tempfile(fileext = ".csv")
  on.exit(unlink(path), add = TRUE)
  writeLines(c("a,b", "1,2", "3,4,5"), path)
  expect_error(.sab_shared_pool_read_joint_row(path, "a", 1L, 2L),
               "row width changed")
  writeLines(c("a,b", "1,2", "3"), path)
  expect_error(.sab_shared_pool_read_joint_row(path, "a", 1L, 2L),
               "row width changed")
})

test_that("snapshot parser fails closed on nonfinite selected values", {
  path <- tempfile(fileext = ".csv")
  on.exit(unlink(path), add = TRUE)
  for (value in c("NA", "NaN", "Inf", "-Inf")) {
    writeLines(c("a,b", paste0("1,", value)), path)
    expect_error(.sab_shared_pool_read_joint_row(path, c("a", "b"), 1L, 1L),
                 "retained row/count is invalid")
  }
})

test_that("snapshot indexing survives a streaming chunk boundary", {
  path <- tempfile(fileext = ".csv")
  on.exit(unlink(path), add = TRUE)
  writeLines(c("# preamble", "a,b",
               paste(seq_len(260L), -seq_len(260L), sep = ",")), path)
  expect_identical(.sab_shared_pool_read_joint_row(path, "b", 255L, 260L),
                   c(b = -255))
})
