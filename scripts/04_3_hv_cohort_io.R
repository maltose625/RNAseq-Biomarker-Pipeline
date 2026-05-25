################################################
##### 04_3_hv_cohort_io.R (队列数据清洗模块)#####
################################################

read_geo_pheno_manual <- function(series_matrix_path,
                                  biopsy_col = NULL,
                                  biopsy_keep_pattern = NULL,
                                  risk_col,
                                  high_label = "high-risk",
                                  title_col = "title") {
  log_info("手动解析 GEO series_matrix 文件: %s", series_matrix_path)

  if (!file.exists(series_matrix_path)) {
    log_stop("series_matrix 文件不存在: %s", series_matrix_path)
  }

  con <- if (grepl("\\.gz$", series_matrix_path, ignore.case = TRUE)) {
    gzfile(series_matrix_path, open = "rt")
  } else {
    file(series_matrix_path, open = "rt")
  }
  on.exit(close(con), add = TRUE)

  lines <- readLines(con, warn = FALSE)

  sample_lines <- grep("^!Sample_", lines, value = TRUE)
  if (length(sample_lines) == 0) {
    log_stop("未在 series_matrix 中找到 !Sample_ 开头的样本注释行: %s", series_matrix_path)
  }

  split_tab <- function(x) strsplit(x, "\t", fixed = TRUE)[[1]]

  parsed <- lapply(sample_lines, split_tab)
  max_len <- max(vapply(parsed, length, integer(1)))
  parsed <- lapply(parsed, function(x) {
    length(x) <- max_len
    x[is.na(x)] <- ""
    x
  })

  mat <- do.call(rbind, parsed)

  # 第一列是字段名，其余列是各个样本的值
  row_keys <- sub("^!Sample_", "", mat[, 1])
  value_mat <- mat[, -1, drop = FALSE]

  # 关键修复：先设置行名，再转置
  rownames(value_mat) <- row_keys
  df <- as.data.frame(t(value_mat), stringsAsFactors = FALSE, check.names = FALSE)

  clean_name <- function(x) {
    x <- trimws(x)
    x <- gsub('"', "", x, fixed = TRUE)
    x <- tolower(x)
    x <- gsub("[[:space:]]+", " ", x)
    x
  }

  make_unique_names <- function(x) {
    make.unique(x, sep = "__dup")
  }

  # 清洗列名
  colnames(df) <- make_unique_names(vapply(colnames(df), clean_name, character(1)))

  # 清洗值
  clean_value <- function(x) {
    x <- trimws(x)
    x <- gsub('^"|"$', "", x)
    x
  }
  df[] <- lapply(df, clean_value)

  # 展开 characteristics_ch1 / characteristics_ch1__dup*
  char_cols <- grep("^characteristics_ch1", colnames(df), value = TRUE)

  if (length(char_cols) > 0) {
    for (cc in char_cols) {
      vals <- trimws(df[[cc]])
      has_kv <- grepl(":", vals, fixed = TRUE)

      if (sum(has_kv, na.rm = TRUE) == 0) next

      key_part <- trimws(sub(":.*$", "", vals[has_kv]))
      val_part <- trimws(sub("^[^:]+:\\s*", "", vals))

      uniq_keys <- unique(key_part[nzchar(key_part)])

      # 只有当这一列大多数值都共享同一个 key 时，才展开成新列
      if (length(uniq_keys) == 1) {
        new_name <- clean_name(uniq_keys[1])

        # 避免重名覆盖
        if (new_name %in% colnames(df)) {
          new_name <- make.unique(c(colnames(df), new_name), sep = "__dup")[length(colnames(df)) + 1]
        }

        df[[new_name]] <- val_part
      }
    }
  }

  normalize_name <- function(x) {
    x <- tolower(trimws(x))
    x <- gsub("[[:space:]]+", " ", x)
    x
  }

  find_col <- function(target, available_cols) {
    if (is.null(target) || is.na(target) || !nzchar(target)) return(NULL)

    key <- normalize_name(target)
    norm_cols <- normalize_name(available_cols)

    idx <- which(norm_cols == key)
    if (length(idx) >= 1) return(available_cols[idx[1]])

    NULL
  }

  title_col_real  <- find_col(title_col, colnames(df))
  risk_col_real   <- find_col(risk_col, colnames(df))
  biopsy_col_real <- find_col(biopsy_col, colnames(df))

  missing_targets <- c()
  if (is.null(title_col_real)) missing_targets <- c(missing_targets, title_col)
  if (is.null(risk_col_real))  missing_targets <- c(missing_targets, risk_col)

  if (length(missing_targets) > 0) {
    log_stop(
      "解析后未找到 title/risk 对应列，请检查列名配置。当前列名: %s",
      paste(colnames(df), collapse = ", ")
    )
  }

  out <- data.frame(
    title = df[[title_col_real]],
    risk_raw = df[[risk_col_real]],
    stringsAsFactors = FALSE
  )

  if (!is.null(biopsy_col_real)) {
    out$biopsy <- df[[biopsy_col_real]]
  }

  out$title <- trimws(out$title)
  out$risk_raw <- trimws(out$risk_raw)

  # biopsy 过滤可选
  if ("biopsy" %in% colnames(out) &&
      !is.null(biopsy_keep_pattern) &&
      !is.na(biopsy_keep_pattern) &&
      nzchar(biopsy_keep_pattern)) {
    keep_idx <- grepl(biopsy_keep_pattern, out$biopsy, ignore.case = TRUE)
    out <- out[keep_idx, , drop = FALSE]
    log_info("按 biopsy_keep_pattern = '%s' 过滤后保留 %d 个样本",
             biopsy_keep_pattern, nrow(out))
  } else {
    log_warn("未提供可用的 biopsy 列或过滤条件，跳过 biopsy 过滤")
  }

  out <- out[!is.na(out$risk_raw) & nzchar(out$risk_raw), , drop = FALSE]

  high_label_norm <- tolower(trimws(high_label))
  risk_norm <- tolower(trimws(out$risk_raw))

  out$RiskLabel <- ifelse(risk_norm == high_label_norm, 1L, 0L)

  log_info(
    "表型解析完成：%d samples；high-risk = %d, low-risk = %d",
    nrow(out),
    sum(out$RiskLabel == 1, na.rm = TRUE),
    sum(out$RiskLabel == 0, na.rm = TRUE)
  )

  out
}



read_gct_matrix <- function(path) {
  if (!file.exists(path)) log_stop("未找到 GCT 文件: %s", path)
  log_info("读取人类表达矩阵 GCT: %s", path)

  df <- data.table::fread(
    path,
    skip = 2,
    header = TRUE,
    data.table = FALSE,
    check.names = FALSE,
    quote = ""
  )

  if (ncol(df) < 4) {
    log_stop("GCT 文件列数异常，当前仅 %d 列", ncol(df))
  }

  row_ids <- trimws(as.character(df[[1]]))
  row_ids <- strip_version(row_ids)

  expr_df <- df[, -(1:2), drop = FALSE]
  expr_mat <- as.matrix(expr_df)
  storage.mode(expr_mat) <- "double"

  rownames(expr_mat) <- make.unique(row_ids)
  expr_mat <- expr_mat[complete.cases(expr_mat), , drop = FALSE]
  expr_mat <- dedup_by_mean(expr_mat)

  log_info("GCT 读取完成：%d genes × %d samples", nrow(expr_mat), ncol(expr_mat))
  log_info("表达矩阵ID类型判断: %s", detect_id_type(rownames(expr_mat)))
  expr_mat
}

align_pheno_expr <- function(pheno, expr, title_col, cohort_name = "cohort", min_common = 20) {
  if (!title_col %in% colnames(pheno)) {
    log_stop("%s 临床表缺少样本标题列: %s；当前列名: %s",
             cohort_name, title_col, paste(colnames(pheno), collapse = ", "))
  }

  expr_samples_raw <- colnames(expr)
  pheno_samples_raw <- pheno[[title_col]]

  expr_samples_norm <- normalize_sample_id(expr_samples_raw)
  pheno_samples_norm <- normalize_sample_id(pheno_samples_raw)

  log_info("%s 表达矩阵样本数: %d", cohort_name, length(expr_samples_raw))
  log_info("%s 临床有效样本数: %d", cohort_name, nrow(pheno))
  log_info("%s 表达矩阵前6个样本: %s", cohort_name, paste(head(expr_samples_raw), collapse = ", "))
  log_info("%s 临床title前6个样本: %s", cohort_name, paste(head(pheno_samples_raw), collapse = ", "))

  common_samples_norm <- intersect(expr_samples_norm, pheno_samples_norm)
  log_info("%s 共同样本数: %d", cohort_name, length(common_samples_norm))

  if (length(common_samples_norm) < min_common) {
    log_stop("%s 共同样本数过少: %d，无法进行稳健分析", cohort_name, length(common_samples_norm))
  }

  expr_idx <- match(common_samples_norm, expr_samples_norm)
  pheno_idx <- match(common_samples_norm, pheno_samples_norm)

  expr2 <- expr[, expr_idx, drop = FALSE]
  pheno2 <- pheno[pheno_idx, , drop = FALSE]

  aligned_ids <- expr_samples_raw[expr_idx]
  colnames(expr2) <- aligned_ids
  rownames(pheno2) <- aligned_ids
  pheno2$MatchedSampleID <- aligned_ids

  log_info("%s 样本对齐完成：%d 个共同样本", cohort_name, ncol(expr2))

  list(pheno = pheno2, expr = expr2)
}
