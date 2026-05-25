################################################
####### 04_1_hv_common.R (基础设施层)###########
################################################

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

log_info <- function(...) cat(sprintf("[%s] [INFO] %s\n", Sys.time(), sprintf(...)))
log_warn <- function(...) cat(sprintf("[%s] [WARN] %s\n", Sys.time(), sprintf(...)))
log_stop <- function(...) stop(sprintf("[%s] [ERROR] %s", Sys.time(), sprintf(...)), call. = FALSE)

read_config <- function() {
  cfg_path <- "config/config.yaml"
  if (!file.exists(cfg_path)) log_stop("未找到配置文件: %s", cfg_path)
  yaml::read_yaml(cfg_path)
}

save_session_info <- function(out_file) {
  sink(out_file)
  print(sessionInfo())
  sink()
}

safe_numeric_matrix <- function(mat) {
  mat <- as.matrix(mat)
  storage.mode(mat) <- "double"
  mat
}

strip_version <- function(x) {
  sub("\\..*$", "", x)
}

detect_id_type <- function(ids) {
  ids <- ids[!is.na(ids) & nzchar(ids)]
  if (length(ids) == 0) return("unknown")
  head_ids <- head(ids, 100)

  if (mean(grepl("^ENSMUSG", head_ids, ignore.case = TRUE)) > 0.5) return("mouse_ensembl")
  if (mean(grepl("^ENSG", head_ids, ignore.case = TRUE)) > 0.5) return("human_ensembl")
  return("symbol")
}

normalize_sample_id <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub('^"|"$', "", x)
  x <- gsub("\\.gz$", "", x, ignore.case = TRUE)
  x <- gsub("\\.txt$", "", x, ignore.case = TRUE)
  x <- gsub("\\.bam$", "", x, ignore.case = TRUE)
  x <- gsub("[[:space:]]+", "", x)
  toupper(x)
}

dedup_by_mean <- function(mat) {
  if (!anyDuplicated(rownames(mat))) return(mat)
  ord <- order(rownames(mat), -rowMeans(mat, na.rm = TRUE))
  mat2 <- mat[ord, , drop = FALSE]
  mat2[!duplicated(rownames(mat2)), , drop = FALSE]
}

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}
