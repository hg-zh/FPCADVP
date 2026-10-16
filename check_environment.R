#!/usr/bin/env Rscript

# Run from the distributed repository root:
#   Rscript check_environment.R
#   Rscript check_environment.R full=TRUE
# Checks dependencies only; never installs packages or changes the environment.

local({
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) > 1L ||
      (length(args) == 1L && !args %in% c("full=TRUE", "full=FALSE"))) {
    stop("Usage: Rscript check_environment.R [full=TRUE|full=FALSE]", call. = FALSE)
  }
  full <- identical(args, "full=TRUE")

  # Prefer the active source() file when called from another R script.
  source_paths <- Filter(function(path) is.character(path) && length(path) == 1L,
                         lapply(sys.frames(), function(frame) frame$ofile))
  if (length(source_paths)) {
    this_file <- tail(source_paths, 1L)[[1L]]
  } else {
    file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    if (length(file_arg) != 1L) {
      stop("Run this file with Rscript or source its saved path.", call. = FALSE)
    }
    this_file <- sub("^--file=", "", file_arg)
    # Rscript encodes spaces in --file on some platforms.
    this_file <- gsub("~+~", " ", this_file, fixed = TRUE)
  }
  dependency_file <- file.path(dirname(normalizePath(this_file, mustWork = TRUE)),
                               "dependencies.tsv")
  dependencies <- read.delim(dependency_file, stringsAsFactors = FALSE,
                             check.names = FALSE)
  required_columns <- c("component", "reference_version", "required_for",
                        "version_policy", "notes")
  if (!identical(names(dependencies), required_columns) ||
      anyDuplicated(dependencies$component)) {
    stop("Invalid dependencies.tsv schema.", call. = FALSE)
  }

  selected <- dependencies[dependencies$required_for == "core" |
                             (full & dependencies$required_for == "html"), ]
  actual <- vapply(selected$component, function(package) {
    if (package == "R") return(as.character(getRversion()))
    if (!requireNamespace(package, quietly = TRUE)) return("MISSING")
    as.character(utils::packageVersion(package))
  }, character(1))
  status <- ifelse(actual == selected$reference_version, "reference", "different version")
  status[actual == "MISSING"] <- "missing"
  r_index <- which(selected$component == "R")
  if (getRversion() >= numeric_version(selected$reference_version[r_index])) {
    status[r_index] <- "meets recommended minimum"
  } else {
    status[r_index] <- "below recommended minimum"
  }
  report <- data.frame(component = selected$component,
                       reference = selected$reference_version,
                       installed = unname(actual), status = unname(status),
                       stringsAsFactors = FALSE)
  print(report, row.names = FALSE, right = FALSE)
  cat("\nR >= 4.5.0 is recommended; reference computations used R 4.5.0 or 4.5.1.\n")
  cat("Fresh sensitivity simulations enforce exactly R 4.5.0 and fdapace 0.6.0.\n")
  cat("Other listed package versions are reference versions; a mismatch is reported, not silently corrected.\n")
  if (!full) cat("HTML dependencies were not checked; use full=TRUE for Figure 2 HTML and optional mortality HTML.\n")
  cat("HCP preprocessing also requires Python 3; its scripts use only the standard library.\n")
  cat("Full cluster runs additionally require the documented Bash/Slurm setup.\n")

  missing <- selected$component[actual == "MISSING"]
  exact_mismatch <- selected$component[selected$version_policy == "exact" &
                                        actual != selected$reference_version]
  problems <- character(0)
  if (length(missing)) {
    problems <- c(problems, paste("Missing required packages:", paste(missing, collapse = ", ")))
  }
  if (length(exact_mismatch)) {
    problems <- c(problems, paste("Exact version required:", paste(exact_mismatch, collapse = ", "),
                                 "(see reference column)"))
  }
  if (length(problems)) {
    stop(paste(c(problems, "Install the documented dependencies, then rerun this check."),
               collapse = "\n"), call. = FALSE)
  }
  if (status[r_index] == "below recommended minimum") {
    warning("This R version is older than the recommended environment.", call. = FALSE)
  }
  cat("PASS: all selected dependencies are available and fdapace matches the reference.\n")
})
