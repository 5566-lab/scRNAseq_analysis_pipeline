#!/usr/bin/env Rscript
repo <- normalizePath(getwd(), mustWork = TRUE)
files <- list.files(repo, pattern = "\\.[Rr]$", recursive = TRUE, full.names = TRUE)
files <- files[!grepl("archive/original_scripts", files)]
for (file in files) {
  parse(file)
  message("OK ", file)
}

required <- c(
  "README.md",
  "configs/config.yaml",
  "metadata/sample_manifest.csv",
  "docs/code_function_archive.md"
)
missing <- required[!file.exists(file.path(repo, required))]
if (length(missing) > 0) {
  stop("Missing required files: ", paste(missing, collapse = ", "), call. = FALSE)
}

message("Repository validation passed")
