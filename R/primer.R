#!/usr/bin/env Rscript
# =============================================================================
# R environment diagnostics script
# =============================================================================
# purpose: comprehensive environment and dependency check for R scripts,
#          particularly useful for computational biology and bioinformatics
#          workflows using renv for reproducibility.
#
# what it checks:
#   - library paths and renv configuration
#   - R session information (version, platform, working directory)
#   - system memory status (Linux/macOS/Windows compatible)
#   - critical package availability and versions
#   - parallel computing configuration
#   - temporary directory accessibility
#
# usage: source at the beginning of your analysis scripts or run standalone
#        source("diagnostics.R")
#        or
#        Rscript diagnostics.R
#
# platform: cross-platform (Linux, macOS, Windows)
# author: for computational biology workflows
# =============================================================================

cat("\n================ R / renv Diagnostics ================\n")

# ----- environment variables
# key renv and R-related environment variables that control package locations.
# useful for debugging renv configuration and package installation issues.
cat("\n--- Environment Variables ---\n")
env_vars <- c(
  RENV_PATHS_LIBRARY = Sys.getenv("RENV_PATHS_LIBRARY"),
  RENV_PATHS_CACHE   = Sys.getenv("RENV_PATHS_CACHE"),
  RENV_PATHS_SANDBOX = Sys.getenv("RENV_PATHS_SANDBOX"),
  R_LIBS_USER        = Sys.getenv("R_LIBS_USER")
)

# format output nicely
for (i in seq_along(env_vars)) {
  value <- if (env_vars[i] == "") "<not set>" else env_vars[i]
  cat(sprintf("  %-20s %s\n", names(env_vars)[i], value))
}

# ----- memory information
# shows system memory status and current R session memory usage.
# important for large omics datasets to catch memory issues early.
cat("\n--- Memory Info ---\n")
if (.Platform$OS.type == "unix") {
  if (Sys.info()["sysname"] == "Darwin") {
    # macOS: parse vm_stat and convert to human-readable format
    vm_info <- system("vm_stat", intern = TRUE)
    page_size_line <- grep("page size", vm_info, value = TRUE)
    page_size <- as.numeric(gsub(".*?(\\d+).*", "\\1", page_size_line))
    
    # extract key metrics
    parse_pages <- function(line) {
      as.numeric(gsub(".*?([0-9]+)\\.", "\\1", line))
    }
    
    free_pages <- parse_pages(grep("Pages free", vm_info, value = TRUE))
    active_pages <- parse_pages(grep("Pages active", vm_info, value = TRUE))
    inactive_pages <- parse_pages(grep("Pages inactive", vm_info, value = TRUE))
    
    # convert to GB
    free_gb <- round(free_pages * page_size / 1024^3, 2)
    active_gb <- round(active_pages * page_size / 1024^3, 2)
    inactive_gb <- round(inactive_pages * page_size / 1024^3, 2)
    total_gb <- free_gb + active_gb + inactive_gb
    
    cat(sprintf("  %-20s %.2f GB\n", "Memory free:", free_gb))
    cat(sprintf("  %-20s %.2f GB\n", "Memory active:", active_gb))
    cat(sprintf("  %-20s %.2f GB\n", "Memory inactive:", inactive_gb))
    cat(sprintf("  %-20s %.2f GB\n", "Total:", total_gb))
    
  } else {
    # Linux: use free -h
    mem_info <- system("free -h", intern = TRUE)
    cat(paste(mem_info, collapse = "\n"), "\n")
  }
} else {
  # Windows-specific memory limit
  mem_limit <- tryCatch(memory.limit(), error = function(e) "N/A")
  cat("  Memory limit:", mem_limit, "\n")
}

# R session memory in MB or GB
session_mem_mb <- round(sum(gc()[,2]))
if (session_mem_mb > 1024) {
  cat(sprintf("  %-20s %.2f GB\n", "R session memory:", session_mem_mb / 1024))
} else {
  cat(sprintf("  %-20s %d MB\n", "R session memory:", session_mem_mb))
}

# # ----- critical package availability
# # checks if commonly used bioinformatics and data analysis packages are
# # installed and shows their versions. Customize this list for your needs.
# cat("\n--- Critical Bioconductor/Omics Packages ---\n")
# critical_pkgs <- c("BiocManager", "DESeq2", "edgeR", "Seurat", 
#                    "data.table", "ggplot2", "dplyr", "tidyr")
# for (pkg in critical_pkgs) {
#   status <- if (requireNamespace(pkg, quietly = TRUE)) {
#     paste0("✓ (v", packageVersion(pkg), ")")
#   } else {
#     "✗ NOT INSTALLED"
#   }
#   cat(sprintf("  %-15s %s\n", pkg, status))
# }

# ----- temporary directory
# verifies the temp directory is accessible and writable.
# many packages write intermediate files here during analysis.
cat("\n--- Temporary Directory ---\n")
cat(sprintf("  %-20s %s\n", "tempdir():", tempdir()))
cat(sprintf("  %-20s %s\n", "Writable:", file.access(tempdir(), 2) == 0))

cat("\n======================================================\n")

# ----- renv status
# checks if renv is active and shows project-specific library paths.
# critical for ensuring reproducible package environments.
if (requireNamespace("renv", quietly = TRUE)) {
  cat("\n--- renv Status ---\n")
  cat(sprintf("  %-20s %s\n", "renv active:", "renv" %in% loadedNamespaces()))
  
  cat("\n--- renv Paths ---\n")
  renv_lib <- tryCatch(renv::paths$library(), error = function(e) "<error>")
  renv_cache <- tryCatch(renv::paths$cache(), error = function(e) "<error>")
  renv_sandbox <- tryCatch(renv::paths$sandbox(), error = function(e) "<error>")
  
  cat(sprintf("  %-20s %s\n", "library:", renv_lib))
  cat(sprintf("  %-20s %s\n", "cache:", renv_cache))
  cat(sprintf("  %-20s %s\n", "sandbox:", renv_sandbox))
} else {
  cat("\nrenv package not available\n")
}

# ----- parallel computing configuration
# shows available CPU cores and current threading settings.
# important for optimizing performance of multi-threaded operations
# (common in DESeq2, data.table, etc.)
cat("\n--- Parallel Backend ---\n")
cat(sprintf("  %-20s %d\n", "Cores available:", parallel::detectCores()))
dt_threads <- tryCatch(data.table::getDTthreads(), error = function(e) "N/A")
cat(sprintf("  %-20s %s\n", "data.table threads:", dt_threads))

# ----- R session information
# basic session info including R version, platform, and working directory.
# essential for reproducibility and debugging platform-specific issues.
cat("\n--- R Session Info ---\n")
cat(sprintf("  %-20s %s\n", "R version:", R.version.string))
cat(sprintf("  %-20s %s\n", "Platform:", R.version$platform))
cat(sprintf("  %-20s %s\n", "Working directory:", getwd()))

# ----- library paths
# shows where R looks for packages. Important for diagnosing package loading
# issues and confirming renv isolation is working correctly.
cat("\n--- Library Paths (.libPaths) ---\n")
lib_paths <- .libPaths()
for (i in seq_along(lib_paths)) {
  cat(sprintf("  [%d] %s\n", i, lib_paths[i]))
}

cat("\n======================================================\n")

# optional: return diagnostic info as a list (useful if sourced)
invisible(list(
  libpaths = .libPaths(),
  renv_active = "renv" %in% loadedNamespaces(),
  r_version = R.version.string,
  wd = getwd(),
  cores = parallel::detectCores(),
  tempdir = tempdir()
))