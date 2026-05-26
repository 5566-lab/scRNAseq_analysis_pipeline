suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
})

parse_config_arg <- function(default_config = "configs/config.yaml") {
  option_list <- list(
    make_option(c("-c", "--config"), type = "character", default = default_config,
                help = "Path to YAML config file")
  )
  opt <- parse_args(OptionParser(option_list = option_list))
  normalizePath(opt$config, mustWork = TRUE)
}

read_pipeline_config <- function(path) {
  cfg <- yaml::read_yaml(path)
  cfg$.config_dir <- dirname(normalizePath(path, mustWork = TRUE))
  cfg$.repo_root <- normalizePath(file.path(cfg$.config_dir, ".."), mustWork = TRUE)
  cfg
}

project_path <- function(cfg, ...) {
  path <- file.path(...)
  if (grepl("^/", path)) {
    return(path)
  }
  file.path(cfg$.repo_root, path)
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  invisible(path)
}

ensure_parent_dir <- function(path) {
  ensure_dir(dirname(path))
  invisible(path)
}

message_step <- function(...) {
  message(sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...)))
}

load_config <- function(default_config = "configs/config.yaml") {
  read_pipeline_config(parse_config_arg(default_config))
}

