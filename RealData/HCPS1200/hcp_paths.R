# Shared path handling for the HCP command-line entrypoints.
`%||%` <- function(x, y) if (!is.null(x)) x else y
hcp_parse_args <- function(args) {
  out <- list()
  for (a in args) {
    if (!grepl("=", a, fixed = TRUE)) next
    out[[sub("=.*$", "", a)]] <- sub("^[^=]*=", "", a)
  }
  out
}
hcp_resolve_path <- function(path, root) {
  if (grepl("^(/|[A-Za-z]:[/\\\\]|~)", path)) path.expand(path) else file.path(root, path)
}
