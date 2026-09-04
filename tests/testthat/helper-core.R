core_roots <- c(
  normalizePath(".", mustWork = TRUE),
  normalizePath(file.path("..", ".."), mustWork = TRUE)
)
core_roots <- unique(core_roots)
core_root <- core_roots[
  vapply(
    core_roots,
    function(path) file.exists(file.path(path, "R", "system_a_adapter.R")),
    logical(1)
  )
][1]

if (is.na(core_root)) {
  stop("Could not locate the project R directory.")
}

source(file.path(core_root, "R", "system_a_adapter.R"), local = FALSE)
