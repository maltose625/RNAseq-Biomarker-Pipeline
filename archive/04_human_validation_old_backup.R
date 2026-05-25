#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(yaml)
  library(data.table)
  library(glmnet)
  library(pROC)
  library(ggplot2)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

log_info <- function(...) cat(sprintf("[%s] [INFO] %s\n", Sys.time(), sprintf(...)))
log_warn <- function(...) cat(sprintf("[%s] [WARN] %s\n", Sys.time(), sprintf(...)))
log_stop <- function(...) stop(sprintf("[%s] [ERROR] %s", Sys.time(), sprintf(...)))

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

resolve_mouse_gene_file <- function(config) {
  env_file <- Sys.getenv("MOUSE_GENE_FILE", unset = "")
  if (nzchar(env_file)) {
    if (!file.exists(env_file)) {
      log_stop("环境变量 MOUSE_GENE_FILE 指向的文件不存在: %s", env_file)
    }
    log_info("优先使用环境变量指定的小鼠输入基因文件: %s", env_file)
    return(env_file)
  }

  core_file <- config$paths$mouse_core_biomarkers %||% ""
  cand_file <- config$paths$mouse_candidate_genes %||% ""
  skip_file <- file.path(config$paths$results_ml %||% "results/03_MachineLearning",
                         "mouse_feature_selection", "ML_SKIP_REASON.txt")

  if (file.exists(skip_file) && nzchar(cand_file) && file.exists(cand_file)) {
    log_warn("检测到 ML_SKIP_REASON.txt，优先使用候选基因池: %s", cand_file)
    return(cand_file)
  }

  if (nzchar(core_file) && file.exists(core_file)) {
    log_info("使用核心标志物文件: %s", core_file)
    return(core_file)
  }

  if (nzchar(cand_file) && file.exists(cand_file)) {
    log_warn("未找到核心标志物文件，回退到候选基因池: %s", cand_file)
    return(cand_file)
  }

  log_stop("未找到可用的小鼠输入基因文件：请检查 MOUSE_GENE_FILE、mouse_core_biomarkers 或 mouse_candidate_genes")
}

read_candidate_genes <- function(path) {
  if (!file.exists(path)) log_stop("未找到候选基因文件: %s", path)

  log_info("读取小鼠候选基因文件: %s", path)

  if (grepl("\\.csv$", path, ignore.case = TRUE)) {
    df <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    df <- data.table::fread(path, data.table = FALSE, check.names = FALSE)
  }

  if (ncol(df) < 1) log_stop("候选基因文件为空: %s", path)

  gene_col <- intersect(
    c("MouseGene", "Gene", "gene", "Symbol", "symbol", "Feature"),
    colnames(df)
  )

  if (length(gene_col) > 0) {
    genes <- unique(trimws(as.character(df[[gene_col[1]]])))
  } else {
    genes <- unique(trimws(as.character(df[[1]])))
  }

  genes <- genes[!is.na(genes) & nzchar(genes)]

  if (length(genes) == 0) {
    log_stop("候选基因文件中没有有效基因")
  }

  log_info("读取到小鼠候选基因数: %d", length(genes))
  genes
}

read_mouse_deg_map <- function(path) {
  if (!file.exists(path)) {
    log_warn("未找到 DESeq2 结果文件用于 Ensembl->Symbol 辅助映射: %s", path)
    return(NULL)
  }

  deg <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)

  if (!"Ensembl" %in% colnames(deg) || !"Symbol" %in% colnames(deg)) {
    log_warn("DESeq2 结果缺少 Ensembl 或 Symbol 列，无法构建辅助映射")
    return(NULL)
  }

  map <- data.frame(
    MouseEnsembl = strip_version(as.character(deg$Ensembl)),
    MouseSymbol  = trimws(as.character(deg$Symbol)),
    stringsAsFactors = FALSE
  )

  map <- map[!is.na(map$MouseEnsembl) & nzchar(map$MouseEnsembl), , drop = FALSE]
  map <- map[!is.na(map$MouseSymbol)  & nzchar(map$MouseSymbol),  , drop = FALSE]
  unique(map)
}

read_ortholog_map <- function(path) {
  if (is.null(path) || !file.exists(path)) {
    log_stop("未找到严格同源映射文件 ortholog_map: %s", path)
  }

  log_info("读取离线同源映射表: %s", path)
  map <- data.table::fread(path, data.table = FALSE, check.names = FALSE)

  if (nrow(map) == 0 || ncol(map) < 2) {
    log_stop("ortholog_map 文件为空或列数不足: %s", path)
  }

  cn_raw <- colnames(map)
  cn <- tolower(gsub("[ .-]+", "_", cn_raw))
  colnames(map) <- cn

  required_cols <- c("mouse_ensembl", "mouse_symbol", "human_ensembl", "human_symbol")
  miss <- setdiff(required_cols, colnames(map))
  if (length(miss) > 0) {
    log_stop("ortholog_map 缺少必要列: %s", paste(miss, collapse = ", "))
  }

  if ("orthology_type" %in% colnames(map)) {
    map <- map[map$orthology_type == "ortholog_one2one", , drop = FALSE]
  }

  if ("orthology_confidence" %in% colnames(map)) {
    suppressWarnings({
      conf_num <- as.numeric(map$orthology_confidence)
    })
    keep_conf <- !is.na(conf_num) & conf_num == 1
    if (any(keep_conf)) {
      map <- map[keep_conf, , drop = FALSE]
    }
  }

  map <- unique(map)

  if (nrow(map) == 0) {
    log_stop("严格 ortholog 过滤后无可用映射，请检查文件内容")
  }

  list(data = map, cols = colnames(map))
}

build_mouse_to_human_mapping <- function(mouse_genes, expr_row_ids, ortholog_file = NULL, mouse_deg_file = NULL) {
  expr_id_type <- detect_id_type(expr_row_ids)
  mouse_id_type <- detect_id_type(mouse_genes)

  log_info("候选基因ID类型: %s", mouse_id_type)
  log_info("人类表达矩阵行名ID类型: %s", expr_id_type)

  deg_map <- read_mouse_deg_map(mouse_deg_file)
  orth <- read_ortholog_map(ortholog_file)

  candidate_df <- data.frame(
    MouseGene = unique(mouse_genes),
    stringsAsFactors = FALSE
  )

  if (!is.null(deg_map)) {
    candidate_df <- merge(
      candidate_df,
      deg_map,
      by.x = "MouseGene",
      by.y = "MouseEnsembl",
      all.x = TRUE
    )
  } else {
    candidate_df$MouseSymbol <- NA_character_
  }

  is_symbol_input <- !grepl("^ENSMUSG", candidate_df$MouseGene, ignore.case = TRUE)
  candidate_df$MouseSymbol[is_symbol_input & is.na(candidate_df$MouseSymbol)] <- candidate_df$MouseGene[is_symbol_input]

  map <- orth$data
  cols <- orth$cols

  mouse_symbol_cols  <- intersect(c("mouse_symbol", "mousesymbol", "mouse_gene", "mousegene", "mgi_symbol"), cols)
  mouse_ensembl_cols <- intersect(c("mouse_ensembl", "mouseensembl", "mouse_gene_stable_id", "mouse_ensembl_gene_id"), cols)
  human_symbol_cols  <- intersect(c("human_symbol", "humansymbol", "human_gene", "humangene", "hgnc_symbol"), cols)
  human_ensembl_cols <- intersect(c("human_ensembl", "humanensembl", "human_gene_stable_id", "human_ensembl_gene_id"), cols)

  target_human_col <- NULL
  if (expr_id_type == "human_ensembl" && length(human_ensembl_cols) > 0) {
    target_human_col <- human_ensembl_cols[1]
  } else if (length(human_symbol_cols) > 0) {
    target_human_col <- human_symbol_cols[1]
  } else if (length(human_ensembl_cols) > 0) {
    target_human_col <- human_ensembl_cols[1]
  }

  if (is.null(target_human_col)) {
    log_stop("ortholog_map 中无法找到适合当前表达矩阵ID类型的人类基因列")
  }

  if (mouse_id_type == "mouse_ensembl" && length(mouse_ensembl_cols) > 0) {
    mm_col <- mouse_ensembl_cols[1]
    tmp <- merge(
      candidate_df,
      map[, c(mm_col, target_human_col), drop = FALSE],
      by.x = "MouseGene",
      by.y = mm_col,
      all.x = TRUE
    )
    colnames(tmp)[colnames(tmp) == target_human_col] <- "HumanGene"
    tmp$Method <- "offline_ortholog_map_by_mouse_ensembl"
    tmp <- tmp[, c("MouseGene", "MouseSymbol", "HumanGene", "Method"), drop = FALSE]
    tmp <- tmp[!is.na(tmp$HumanGene) & nzchar(tmp$HumanGene), , drop = FALSE]
    tmp <- unique(tmp)

    if (nrow(tmp) > 0) {
      log_info("使用离线 ortholog_map（mouse_ensembl）成功映射 %d 条", nrow(tmp))
      return(tmp)
    }
  }

  if (length(mouse_symbol_cols) > 0 && any(!is.na(candidate_df$MouseSymbol))) {
    ms_col <- mouse_symbol_cols[1]
    candidate_df2 <- candidate_df[!is.na(candidate_df$MouseSymbol) & nzchar(candidate_df$MouseSymbol), , drop = FALSE]

    tmp <- merge(
      candidate_df2,
      map[, c(ms_col, target_human_col), drop = FALSE],
      by.x = "MouseSymbol",
      by.y = ms_col,
      all.x = TRUE
    )
    colnames(tmp)[colnames(tmp) == target_human_col] <- "HumanGene"
    tmp$Method <- "offline_ortholog_map_by_mouse_symbol"
    tmp <- tmp[, c("MouseGene", "MouseSymbol", "HumanGene", "Method"), drop = FALSE]
    tmp <- tmp[!is.na(tmp$HumanGene) & nzchar(tmp$HumanGene), , drop = FALSE]
    tmp <- unique(tmp)

    if (nrow(tmp) > 0) {
      log_info("使用离线 ortholog_map（mouse_symbol）成功映射 %d 条", nrow(tmp))
      return(tmp)
    }
  }

  log_stop("未能使用离线 ortholog_map 完成严格映射，请检查 ortholog_map 文件、基因ID类型及列名")
}

read_geo_pheno_manual <- function(series_matrix_path,
                                  biopsy_col = "biopsy",
                                  biopsy_keep_pattern = "1st",
                                  risk_col = "pls-nafld-based risk prediction at 1st biopsy",
                                  high_label = "high-risk",
                                  title_col = "title") {
  if (!file.exists(series_matrix_path)) {
    log_stop("未找到 series_matrix 文件: %s", series_matrix_path)
  }

  log_info("手动解析 GEO series_matrix 文件: %s", series_matrix_path)

  lines <- readLines(gzfile(series_matrix_path), warn = FALSE, encoding = "UTF-8")
  sample_lines <- lines[grepl("^!Sample_", lines)]

  if (length(sample_lines) == 0) {
    log_stop("series_matrix 中未找到 !Sample_ 行")
  }

  split_lines <- strsplit(sample_lines, "\t", fixed = TRUE)

  field_names <- sub("^!Sample_", "", vapply(split_lines, `[`, character(1), 1))
  value_list <- lapply(split_lines, function(x) x[-1])

  n_samples <- unique(vapply(value_list, length, integer(1)))
  if (length(n_samples) != 1) {
    log_stop("series_matrix 各 !Sample_ 行的样本数不一致: %s",
             paste(n_samples, collapse = ", "))
  }
  n_samples <- n_samples[1]

  pheno_raw <- data.frame(matrix(NA_character_, nrow = n_samples, ncol = length(field_names)),
                          stringsAsFactors = FALSE)
  colnames(pheno_raw) <- make.unique(field_names, sep = "__dup")
  for (i in seq_along(value_list)) {
    pheno_raw[[i]] <- gsub('^"|"$', "", value_list[[i]])
  }

  pheno <- pheno_raw

  char_idx <- grepl("^characteristics_ch1", colnames(pheno_raw))
  char_df <- pheno_raw[, char_idx, drop = FALSE]

  if (ncol(char_df) == 0) {
    log_stop("未找到 characteristics_ch1 列，无法解析临床特征")
  }

  for (j in seq_len(ncol(char_df))) {
    vals <- trimws(char_df[[j]])
    vals[is.na(vals)] <- ""

    has_colon <- grepl(":", vals)
    if (!any(has_colon)) next

    keys <- vals
    keys[has_colon] <- sub(":.*$", "", vals[has_colon])
    keys <- trimws(tolower(keys))

    values <- vals
    values[has_colon] <- sub("^[^:]+:\\s*", "", vals[has_colon])
    values <- trimws(values)

    uniq_keys <- unique(keys[has_colon])

    for (k in uniq_keys) {
      idx <- which(keys == k)
      if (length(idx) == 0) next

      colname <- k
      if (!colname %in% colnames(pheno)) {
        pheno[[colname]] <- NA_character_
      }

      oldv <- pheno[[colname]]
      newv <- values[idx]

      fill_idx <- is.na(oldv[idx]) | oldv[idx] == ""
      oldv[idx][fill_idx] <- newv[fill_idx]
      pheno[[colname]] <- oldv
    }
  }

  clean_names <- colnames(pheno)
  clean_names <- trimws(clean_names)
  clean_names <- gsub("\\s+", " ", clean_names)
  colnames(pheno) <- clean_names

  required_cols <- c(title_col, biopsy_col, risk_col)
  miss <- setdiff(required_cols, colnames(pheno))
  if (length(miss) > 0) {
    log_stop("手动解析后缺少必要列: %s；当前列名: %s",
             paste(miss, collapse = ", "),
             paste(colnames(pheno), collapse = ", "))
  }

  keep_idx <- grepl(biopsy_keep_pattern, pheno[[biopsy_col]], ignore.case = TRUE)
  pheno <- pheno[keep_idx, , drop = FALSE]

  risk_vec <- trimws(tolower(pheno[[risk_col]]))
  risk_vec[risk_vec %in% c("na", "", "n/a")] <- NA
  pheno <- pheno[!is.na(risk_vec), , drop = FALSE]
  risk_vec <- trimws(tolower(pheno[[risk_col]]))

  pheno$RiskLabel <- ifelse(risk_vec == tolower(high_label), 1, 0)

  log_info("series_matrix 解析完成：保留 1st biopsy 后样本数 = %d", nrow(pheno))
  log_info("高风险样本数 = %d，低风险样本数 = %d",
           sum(pheno$RiskLabel == 1, na.rm = TRUE),
           sum(pheno$RiskLabel == 0, na.rm = TRUE))

  pheno
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

check_binary_outcome <- function(y, min_count = 4, label = "y") {
  tab <- table(y)
  if (length(tab) != 2) {
    log_stop("%s 不是二分类。当前分布: %s", label, paste(names(tab), tab, collapse = ", "))
  }
  if (min(tab) < min_count) {
    log_stop("%s 最小类别样本数不足 %d。当前分布: %s",
             label, min_count, paste(names(tab), tab, collapse = ", "))
  }
  invisible(tab)
}

make_stratified_foldid <- function(y, k, seed = 123) {
  set.seed(seed)
  y <- as.factor(y)
  foldid <- integer(length(y))
  for (lv in levels(y)) {
    idx <- which(y == lv)
    foldid[idx] <- sample(rep(seq_len(k), length.out = length(idx)))
  }
  foldid
}

build_feature_matrix <- function(expr_mat, feature_names) {
  common_features <- intersect(feature_names, rownames(expr_mat))
  x <- matrix(0, nrow = ncol(expr_mat), ncol = length(feature_names))
  colnames(x) <- feature_names
  rownames(x) <- colnames(expr_mat)

  if (length(common_features) > 0) {
    x[, common_features] <- t(expr_mat[common_features, , drop = FALSE])
  }

  list(
    x = safe_numeric_matrix(x),
    common_features = common_features,
    missing_features = setdiff(feature_names, common_features)
  )
}

main <- function() {
  log_info("🚀 启动人类验证模块...")

  config <- read_config()

  out_root <- config$paths$results_ml %||% "results/03_MachineLearning"
  out_dir  <- file.path(out_root, "human_validation")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  series_matrix_file <- config$paths$human_series_matrix
  gct_file <- config$paths$human_gct_matrix

  ext_series_matrix_file <- config$paths$human_series_matrix_external %||% ""
  ext_gct_file <- config$paths$human_gct_matrix_external %||% ""

  candidate_file <- resolve_mouse_gene_file(config)
  ortholog_file <- config$paths$ortholog_map %||% NULL
  mouse_deg_file <- config$paths$deseq2_result_file

  seed <- as.integer(config$params$ml$seed %||% 123)
  alpha <- as.numeric(config$params$ml$alpha %||% 1)
  lambda_choice <- config$params$ml$lambda_choice %||% "lambda.1se"
  standardize <- as.logical(config$params$ml$standardize %||% TRUE)
  min_class_human <- as.integer(config$params$ml$min_class_count_human_model %||% 8)

  hv_cfg <- config$params$human_validation
  ext_cfg <- config$params$human_validation_external %||% hv_cfg

  biopsy_col <- hv_cfg$biopsy_col
  biopsy_keep_pattern <- hv_cfg$biopsy_keep_pattern %||% "1st"
  risk_col <- hv_cfg$risk_col
  high_label <- hv_cfg$high_label %||% "high-risk"
  title_col <- hv_cfg$title_col %||% "title"

  mouse_genes <- read_candidate_genes(candidate_file)
  if (length(mouse_genes) > 100) {
    log_warn("当前输入基因数 = %d，推测为候选基因池模式（candidate mode），而非最终核心 biomarkers 模式", length(mouse_genes))
  } else {
    log_info("当前输入基因数 = %d，推测为核心 biomarkers 模式", length(mouse_genes))
  }

  # =========================
  # Development cohort
  # =========================
  pheno_dev <- read_geo_pheno_manual(
    series_matrix_path = series_matrix_file,
    biopsy_col = biopsy_col,
    biopsy_keep_pattern = biopsy_keep_pattern,
    risk_col = risk_col,
    high_label = high_label,
    title_col = title_col
  )
  expr_dev <- read_gct_matrix(gct_file)

  aligned_dev <- align_pheno_expr(
    pheno = pheno_dev,
    expr = expr_dev,
    title_col = title_col,
    cohort_name = "development cohort",
    min_common = 20
  )
  pheno_dev <- aligned_dev$pheno
  expr_dev <- aligned_dev$expr

  y_dev <- factor(
    ifelse(pheno_dev$RiskLabel == 1, "high-risk", "low-risk"),
    levels = c("low-risk", "high-risk")
  )
  dev_tab <- check_binary_outcome(y_dev, min_count = min_class_human, label = "development cohort")

  mapping_df <- build_mouse_to_human_mapping(
    mouse_genes = mouse_genes,
    expr_row_ids = rownames(expr_dev),
    ortholog_file = ortholog_file,
    mouse_deg_file = mouse_deg_file
  )

  human_genes <- unique(mapping_df$HumanGene)
  available_genes <- intersect(human_genes, rownames(expr_dev))
  missing_genes   <- setdiff(human_genes, available_genes)

  log_info("映射后人类候选基因数: %d", length(human_genes))
  log_info("在人类开发队列表达矩阵中成功匹配基因数: %d", length(available_genes))
  log_info("未匹配基因数: %d", length(missing_genes))

  write.csv(mapping_df, file.path(out_dir, "mouse_to_human_mapping.csv"), row.names = FALSE)
  write.csv(data.frame(AvailableHumanGene = available_genes), file.path(out_dir, "available_human_genes.csv"), row.names = FALSE)
  write.csv(data.frame(MissingHumanGene = missing_genes), file.path(out_dir, "missing_human_genes.csv"), row.names = FALSE)

  if (length(available_genes) < 2) {
    log_stop("候选基因在人类开发队列表达矩阵中匹配不足 2 个，无法建模")
  }

  target_expr_dev <- expr_dev[available_genes, , drop = FALSE]
  x_dev_raw <- t(target_expr_dev)
  x_dev_raw <- safe_numeric_matrix(x_dev_raw)

  dat_dev <- data.frame(
    Sample = rownames(pheno_dev),
    y = y_dev,
    x_dev_raw,
    check.names = FALSE
  )

  x_train <- as.matrix(dat_dev[, setdiff(colnames(dat_dev), c("Sample", "y")), drop = FALSE])
  y_train <- dat_dev$y
  feature_names <- colnames(x_train)

  train_tab <- table(y_train)

  log_info("开发队列大小: %d；类别分布: %s",
           nrow(dat_dev), paste(names(train_tab), train_tab, collapse = ", "))

  if (length(train_tab) != 2 || min(train_tab) < 3) {
    log_stop("开发队列最小类别样本数 < 3，无法稳定执行 cv.glmnet")
  }

  center <- rep(0, ncol(x_train))
  scalev <- rep(1, ncol(x_train))
  names(center) <- colnames(x_train)
  names(scalev) <- colnames(x_train)

  if (standardize) {
    center <- colMeans(x_train)
    scalev <- apply(x_train, 2, sd)
    scalev[is.na(scalev) | scalev == 0] <- 1
    x_train <- scale(x_train, center = center, scale = scalev)
  }

  write.csv(
    data.frame(Feature = feature_names, Center = center, Scale = scalev),
    file.path(out_dir, "feature_scaling_params.csv"),
    row.names = FALSE
  )

  write.csv(
    data.frame(Sample = dat_dev$Sample, Truth = y_train),
    file.path(out_dir, "development_samples.csv"),
    row.names = FALSE
  )

  k <- min(10, min(train_tab))
  foldid <- make_stratified_foldid(y_train, k = k, seed = seed)

  log_info("正在训练开发队列模型：cv.glmnet (nfolds=%d, low-risk=%d, high-risk=%d)",
    k, train_tab["low-risk"], train_tab["high-risk"])
  set.seed(seed)
  cvfit <- cv.glmnet(
    x = x_train,
    y = y_train,
    family = "binomial",
    alpha = alpha,
    standardize = FALSE,
    nfolds = k,
    foldid = foldid,
    type.measure = "auc"
  )

  lam <- if (lambda_choice == "lambda.min") cvfit$lambda.min else cvfit$lambda.1se
  coef_mat <- as.matrix(coef(cvfit, s = lam))
  selected_features <- setdiff(rownames(coef_mat)[coef_mat[, 1] != 0], "(Intercept)")

  write.csv(
    data.frame(Feature = rownames(coef_mat), Coef = coef_mat[, 1]),
    file.path(out_dir, "human_model_coefficients.csv"),
    row.names = FALSE
  )

  write.csv(
    data.frame(HumanGene = selected_features),
    file.path(out_dir, "selected_human_features.csv"),
    row.names = FALSE
  )

  pdf(file.path(out_dir, "01_Human_LASSO_CV_Curve.pdf"), width = 7, height = 5)
  plot(cvfit, main = sprintf("Development Cohort LASSO CV (%s)", lambda_choice))
  dev.off()

  dev_pred_prob <- as.numeric(predict(cvfit, newx = x_train, s = lam, type = "response"))
  dev_pred_class <- ifelse(dev_pred_prob >= 0.5, "high-risk", "low-risk")
  dev_pred_class <- factor(dev_pred_class, levels = c("low-risk", "high-risk"))

  write.csv(
    data.frame(
      Sample = dat_dev$Sample,
      Truth = y_train,
      PredProb = dev_pred_prob,
      PredClass = dev_pred_class
    ),
    file.path(out_dir, "development_apparent_predictions.csv"),
    row.names = FALSE
  )

  dev_auc <- NA_real_
  if (length(train_tab) == 2) {
    roc_dev <- pROC::roc(
      response = y_train,
      predictor = dev_pred_prob,
      levels = c("low-risk", "high-risk"),
      quiet = TRUE
    )
    dev_auc <- as.numeric(pROC::auc(roc_dev))

    pdf(file.path(out_dir, "02_Development_Apparent_ROC_Curve.pdf"), width = 6, height = 6)
    plot(roc_dev, main = sprintf("Development Apparent ROC (AUC = %.3f)", dev_auc), col = "#2ca02c", lwd = 2)
    abline(a = 0, b = 1, lty = 2, col = "gray50")
    dev.off()

    log_info("开发队列表观 AUC = %.4f（注意：这不是外部验证 AUC）", dev_auc)
  }

  # =========================
  # External validation cohort
  # =========================
  external_auc <- NA_real_
  external_n <- 0L
  external_common_feature_n <- 0L

  if (nzchar(ext_series_matrix_file) && nzchar(ext_gct_file) &&
      file.exists(ext_series_matrix_file) && file.exists(ext_gct_file)) {

    log_info("检测到外部独立验证队列，开始执行 external validation ...")

  # 🌟 核心修正：单独指定外部验证集 GSE193080 的独有临床字段名
    pheno_ext <- read_geo_pheno_manual(
      series_matrix_path = config$paths$human_series_matrix_external,
      biopsy_col = "tissue",                         # 外部队列没有 biopsy 列，用 tissue 列占位
      biopsy_keep_pattern = "",                      # 空字符串 "" 在正则中会匹配所有，从而保留全部样本
      risk_col = "pls-nafld-based risk prediction",  # 外部队列的高低风险列名没有 "at 1st biopsy" 后缀
      high_label = high_label, 
      title_col = title_col
    )

    expr_ext <- read_gct_matrix(ext_gct_file)

    aligned_ext <- align_pheno_expr(
      pheno = pheno_ext,
      expr = expr_ext,
      title_col = ext_cfg$title_col %||% title_col,
      cohort_name = "external cohort",
      min_common = 10
    )
    pheno_ext <- aligned_ext$pheno
    expr_ext <- aligned_ext$expr

    y_ext <- factor(
      ifelse(pheno_ext$RiskLabel == 1, "high-risk", "low-risk"),
      levels = c("low-risk", "high-risk")
    )
    ext_tab <- check_binary_outcome(y_ext, min_count = 1, label = "external cohort")

    ext_mat_res <- build_feature_matrix(expr_ext, feature_names = feature_names)
    x_ext <- ext_mat_res$x
    external_common_feature_n <- length(ext_mat_res$common_features)
    external_missing_feature_n <- length(ext_mat_res$missing_features)
    external_n <- nrow(x_ext)

    log_info("外部队列成功匹配特征数: %d；缺失特征数: %d",
             external_common_feature_n, external_missing_feature_n)

    write.csv(
      data.frame(ExternalAvailableFeature = ext_mat_res$common_features),
      file.path(out_dir, "external_available_features.csv"),
      row.names = FALSE
    )
    write.csv(
      data.frame(ExternalMissingFeature = ext_mat_res$missing_features),
      file.path(out_dir, "external_missing_features.csv"),
      row.names = FALSE
    )
    write.csv(
      data.frame(Sample = rownames(pheno_ext), Truth = y_ext),
      file.path(out_dir, "external_validation_samples.csv"),
      row.names = FALSE
    )

    if (standardize) {
      x_ext <- scale(x_ext, center = center, scale = scalev)
    }

    pred_prob_ext <- as.numeric(predict(cvfit, newx = x_ext, s = lam, type = "response"))
    pred_class_ext <- ifelse(pred_prob_ext >= 0.5, "high-risk", "low-risk")
    pred_class_ext <- factor(pred_class_ext, levels = c("low-risk", "high-risk"))

    write.csv(
      data.frame(
        Sample = rownames(pheno_ext),
        Truth = y_ext,
        PredProb = pred_prob_ext,
        PredClass = pred_class_ext
      ),
      file.path(out_dir, "external_predictions.csv"),
      row.names = FALSE
    )

    cm_ext <- as.data.frame.matrix(table(Truth = y_ext, Pred = pred_class_ext))
    write.csv(cm_ext, file.path(out_dir, "external_confusion_matrix.csv"))

    if (length(ext_tab) == 2) {
      roc_ext <- pROC::roc(
        response = y_ext,
        predictor = pred_prob_ext,
        levels = c("low-risk", "high-risk"),
        quiet = TRUE
      )
      external_auc <- as.numeric(pROC::auc(roc_ext))

      pdf(file.path(out_dir, "03_External_ROC_Curve.pdf"), width = 6, height = 6)
      plot(roc_ext, main = sprintf("External Validation ROC (AUC = %.3f)", external_auc), col = "#1f77b4", lwd = 2)
      abline(a = 0, b = 1, lty = 2, col = "gray50")
      dev.off()

      log_info("外部独立验证 AUC = %.4f", external_auc)
    } else {
      log_warn("外部验证队列仅单一类别，跳过 ROC/AUC 计算")
      writeLines("外部验证队列仅单一类别，已跳过 ROC/AUC 计算。",
                 file.path(out_dir, "EXTERNAL_ROC_SKIPPED_REASON.txt"))
    }

  } else {
    log_warn("未检测到可用的外部验证队列文件，当前仅完成 development cohort 建模")
  }
  ## =========================================================
  ## Calibration + Threshold Optimization
  ## 建议放在 04_human_validation.R 的 ROC/AUC 计算之后
  ## =========================================================

  ## ---------- 0. 安全兜底：日志函数 ----------
  if (!exists("log_info")) {
    log_info <- function(fmt, ...) message(sprintf(fmt, ...))
  }
  if (!exists("log_warn")) {
    log_warn <- function(fmt, ...) warning(sprintf(fmt, ...), call. = FALSE)
  }
  if (!exists("log_stop")) {
    log_stop <- function(fmt, ...) stop(sprintf(fmt, ...), call. = FALSE)
  }

  ## ---------- 1. 安全兜底：输出目录 ----------
  if (!exists("outdir")) {
    if (exists("output_dir")) {
      outdir <- output_dir
    } else {
      outdir <- "results/03_MachineLearning/human_validation"
    }
  }
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  ## ---------- 2. 工具函数 ----------
  get_first_existing <- function(candidates, env = parent.frame()) {
    for (nm in candidates) {
      if (exists(nm, envir = env, inherits = TRUE)) {
        return(get(nm, envir = env, inherits = TRUE))
      }
    }
    NULL
  }

  as_binary01 <- function(y) {
    # 约定 high-risk = 1, low-risk = 0
    if (is.factor(y)) y <- as.character(y)
    if (is.logical(y)) return(as.integer(y))
    if (is.numeric(y)) {
      uy <- sort(unique(stats::na.omit(y)))
      if (all(uy %in% c(0, 1))) return(as.integer(y))
      if (all(uy %in% c(1, 2))) return(as.integer(y == 2))
    }
    y <- trimws(as.character(y))
    y_lower <- tolower(y)

    out <- ifelse(
      y_lower %in% c("high-risk", "high_risk", "high risk", "high", "1", "case", "positive"),
      1L,
      ifelse(
        y_lower %in% c("low-risk", "low_risk", "low risk", "low", "0", "control", "negative"),
        0L,
        NA_integer_
      )
    )

    if (any(is.na(out))) {
      bad <- unique(y[is.na(out)])
      log_stop("❌ 无法将标签转换为 0/1。未识别标签: %s", paste(bad, collapse = ", "))
    }
    out
  }

  brier_score <- function(truth01, prob) {
    mean((prob - truth01)^2, na.rm = TRUE)
  }

  make_confusion_df <- function(truth01, pred01) {
    tn <- sum(truth01 == 0 & pred01 == 0, na.rm = TRUE)
    fp <- sum(truth01 == 0 & pred01 == 1, na.rm = TRUE)
    fn <- sum(truth01 == 1 & pred01 == 0, na.rm = TRUE)
    tp <- sum(truth01 == 1 & pred01 == 1, na.rm = TRUE)

    data.frame(
      Truth = c("low-risk", "low-risk", "high-risk", "high-risk"),
      Predicted = c("low-risk", "high-risk", "low-risk", "high-risk"),
      N = c(tn, fp, fn, tp),
      stringsAsFactors = FALSE
    )
  }

  calc_metrics_at_threshold <- function(truth01, prob, threshold) {
    pred01 <- ifelse(prob >= threshold, 1L, 0L)

    tn <- sum(truth01 == 0 & pred01 == 0, na.rm = TRUE)
    fp <- sum(truth01 == 0 & pred01 == 1, na.rm = TRUE)
    fn <- sum(truth01 == 1 & pred01 == 0, na.rm = TRUE)
    tp <- sum(truth01 == 1 & pred01 == 1, na.rm = TRUE)

    sensitivity <- ifelse((tp + fn) == 0, NA_real_, tp / (tp + fn))
    specificity <- ifelse((tn + fp) == 0, NA_real_, tn / (tn + fp))
    ppv <- ifelse((tp + fp) == 0, NA_real_, tp / (tp + fp))
    npv <- ifelse((tn + fn) == 0, NA_real_, tn / (tn + fn))
    accuracy <- (tp + tn) / (tp + tn + fp + fn)
    balanced_accuracy <- mean(c(sensitivity, specificity), na.rm = TRUE)
    f1 <- ifelse((2 * tp + fp + fn) == 0, NA_real_, 2 * tp / (2 * tp + fp + fn))
    brier <- brier_score(truth01, prob)

    auc_val <- tryCatch({
      as.numeric(pROC::auc(truth01, prob))
    }, error = function(e) NA_real_)

    data.frame(
      threshold = threshold,
      AUC = auc_val,
      Brier = brier,
      Accuracy = accuracy,
      Sensitivity = sensitivity,
      Specificity = specificity,
      PPV = ppv,
      NPV = npv,
      BalancedAccuracy = balanced_accuracy,
      F1 = f1,
      TN = tn,
      FP = fp,
      FN = fn,
      TP = tp,
      stringsAsFactors = FALSE
    )
  }

  find_best_threshold_youden <- function(truth01, prob) {
    roc_obj <- pROC::roc(truth01, prob, quiet = TRUE, direction = "<")
    coords <- pROC::coords(
      roc_obj,
      x = "best",
      best.method = "youden",
      transpose = FALSE,
      ret = c("threshold", "sensitivity", "specificity", "accuracy")
    )
    as.data.frame(coords, stringsAsFactors = FALSE)
  }

  build_calibration_table <- function(truth01, prob, n_bins = 10) {
    n <- length(prob)
    if (n < 40) {
      n_bins <- min(5, n)
    } else if (n < 80) {
      n_bins <- min(8, n)
    }

    # 分位数分箱，尽量避免空箱
    qs <- unique(stats::quantile(prob, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE))
    if (length(qs) <= 2) {
      # 概率值过于集中时，退化为等宽分箱
      qs <- seq(min(prob, na.rm = TRUE), max(prob, na.rm = TRUE), length.out = n_bins + 1)
      qs <- unique(qs)
    }
    if (length(qs) <= 2) {
      # 再兜底：全部一个箱
      bin <- rep(1L, length(prob))
    } else {
      bin <- cut(prob, breaks = qs, include.lowest = TRUE, labels = FALSE)
    }

    df <- data.frame(truth01 = truth01, prob = prob, bin = bin)
    df <- df[!is.na(df$bin), , drop = FALSE]

    out <- do.call(
      rbind,
      lapply(split(df, df$bin), function(z) {
        n_bin <- nrow(z)
        obs <- mean(z$truth01)
        pred <- mean(z$prob)
        se <- sqrt(obs * (1 - obs) / max(n_bin, 1))
        lower <- max(0, obs - 1.96 * se)
        upper <- min(1, obs + 1.96 * se)

        data.frame(
          bin = unique(z$bin),
          n = n_bin,
          mean_pred = pred,
          obs_rate = obs,
          obs_lower = lower,
          obs_upper = upper
        )
      })
    )

    rownames(out) <- NULL
    out
  }

  plot_calibration_pdf <- function(calib_df, title, file) {
    grDevices::pdf(file, width = 6.5, height = 6)
    op <- par(no.readonly = TRUE)
    on.exit({
      par(op)
      grDevices::dev.off()
    }, add = TRUE)

    plot(
      calib_df$mean_pred, calib_df$obs_rate,
      xlim = c(0, 1), ylim = c(0, 1),
      xlab = "Mean predicted probability",
      ylab = "Observed event rate",
      main = title,
      pch = 19, cex = 1.2, col = "#1f77b4"
    )
    abline(0, 1, lty = 2, lwd = 2, col = "gray50")

    # 误差线
    segments(
      x0 = calib_df$mean_pred, y0 = calib_df$obs_lower,
      x1 = calib_df$mean_pred, y1 = calib_df$obs_upper,
      col = "#1f77b4", lwd = 1.5
    )

    # 连线
    ord <- order(calib_df$mean_pred)
    lines(
      calib_df$mean_pred[ord],
      calib_df$obs_rate[ord],
      col = "#1f77b4", lwd = 2
    )

    # 标注每箱样本数
    text(
      calib_df$mean_pred,
      calib_df$obs_rate,
      labels = paste0("n=", calib_df$n),
      pos = 3, cex = 0.8, col = "black"
    )
  }

  save_metrics_bundle <- function(prefix, truth01, prob, threshold, outdir) {
    metrics_df <- calc_metrics_at_threshold(truth01, prob, threshold)
    conf_df <- make_confusion_df(truth01, ifelse(prob >= threshold, 1L, 0L))

    utils::write.csv(
      metrics_df,
      file = file.path(outdir, paste0(prefix, "_metrics.csv")),
      row.names = FALSE
    )
    utils::write.csv(
      conf_df,
      file = file.path(outdir, paste0(prefix, "_confusion_matrix.csv")),
      row.names = FALSE
    )
  }

  ## ---------- 3. 自动识别 development / external 变量 ----------
  # 你也可以手动指定：
  # dev_truth_raw <- y_train
  # dev_prob <- train_prob
  # ext_truth_raw <- y_ext
  # ext_prob <- ext_prob

  dev_truth_raw <- get_first_existing(c("y_train", "y_dev", "development_y", "dev_y", "train_y"))
  dev_prob <- get_first_existing(c("dev_prob", "train_prob", "development_prob", "pred_train_prob"))

  ext_truth_raw <- get_first_existing(c("y_ext", "y_external", "external_y", "ext_y"))
  ext_prob <- get_first_existing(c("ext_prob", "external_prob", "pred_ext_prob"))

  # 如果 development prob 不存在，尝试用 cvfit + x_train 现算
  if (is.null(dev_prob)) {
    cvfit_obj <- get_first_existing(c("cvfit"))
    x_train_obj <- get_first_existing(c("x_train"))
    lambda_to_use <- get_first_existing(c("lambda_selected", "lambda_value", "lambda_use", "lambda_chosen"))

    if (is.null(lambda_to_use) && !is.null(cvfit_obj)) {
      lambda_to_use <- cvfit_obj$lambda.1se
    }

    if (!is.null(cvfit_obj) && !is.null(x_train_obj) && !is.null(lambda_to_use)) {
      dev_prob <- as.numeric(stats::predict(cvfit_obj, newx = x_train_obj, s = lambda_to_use, type = "response"))
    }
  }

  # 如果 external prob 不存在，尝试用 cvfit + x_ext 现算
  if (is.null(ext_prob)) {
    cvfit_obj <- get_first_existing(c("cvfit"))
    x_ext_obj <- get_first_existing(c("x_ext", "x_external", "external_x", "ext_x"))
    lambda_to_use <- get_first_existing(c("lambda_selected", "lambda_value", "lambda_use", "lambda_chosen"))

    if (is.null(lambda_to_use) && !is.null(cvfit_obj)) {
      lambda_to_use <- cvfit_obj$lambda.1se
    }

    if (!is.null(cvfit_obj) && !is.null(x_ext_obj) && !is.null(lambda_to_use)) {
      ext_prob <- as.numeric(stats::predict(cvfit_obj, newx = x_ext_obj, s = lambda_to_use, type = "response"))
    }
  }

  if (is.null(dev_truth_raw) || is.null(dev_prob)) {
    log_stop("❌ 无法识别 development cohort 的 truth/prob 变量，请检查 y_train / train_prob 等对象名。")
  }
  dev_truth01 <- as_binary01(dev_truth_raw)
  dev_prob <- as.numeric(dev_prob)

  has_external <- !(is.null(ext_truth_raw) || is.null(ext_prob))
  if (has_external) {
    ext_truth01 <- as_binary01(ext_truth_raw)
    ext_prob <- as.numeric(ext_prob)
  } else {
    log_warn("⚠️ 未识别到 external cohort 的 truth/prob 变量，将只输出 development 的 calibration 和 threshold 结果。")
  }

  ## ---------- 4. Brier score ----------
  dev_brier <- brier_score(dev_truth01, dev_prob)
  log_info("Development Brier score = %.4f", dev_brier)

  brier_df <- data.frame(
    Cohort = "development",
    Brier = dev_brier,
    stringsAsFactors = FALSE
  )

  if (has_external) {
    ext_brier <- brier_score(ext_truth01, ext_prob)
    log_info("External Brier score = %.4f", ext_brier)

    brier_df <- rbind(
      brier_df,
      data.frame(Cohort = "external", Brier = ext_brier, stringsAsFactors = FALSE)
    )
  }

  utils::write.csv(
    brier_df,
    file = file.path(outdir, "04_Brier_Scores.csv"),
    row.names = FALSE
  )

  ## ---------- 5. Calibration ----------
  dev_calib <- build_calibration_table(dev_truth01, dev_prob, n_bins = 10)
  utils::write.csv(
    dev_calib,
    file = file.path(outdir, "05_Development_Calibration_Table.csv"),
    row.names = FALSE
  )
  plot_calibration_pdf(
    dev_calib,
    title = "Development Calibration Plot",
    file = file.path(outdir, "05_Development_Calibration_Plot.pdf")
  )

  if (has_external) {
    ext_calib <- build_calibration_table(ext_truth01, ext_prob, n_bins = 10)
    utils::write.csv(
      ext_calib,
      file = file.path(outdir, "06_External_Calibration_Table.csv"),
      row.names = FALSE
    )
    plot_calibration_pdf(
      ext_calib,
      title = "External Calibration Plot",
      file = file.path(outdir, "06_External_Calibration_Plot.pdf")
    )
  }

  ## ---------- 6. 在 development 上优化阈值 ----------
  best_thr_df <- find_best_threshold_youden(dev_truth01, dev_prob)
  best_threshold <- as.numeric(best_thr_df$threshold[1])

  log_info("Development cohort 最佳阈值（Youden） = %.4f", best_threshold)

  utils::write.csv(
    best_thr_df,
    file = file.path(outdir, "07_Threshold_Optimization_Development_Youden.csv"),
    row.names = FALSE
  )

  ## ---------- 7. 默认阈值 0.5 vs 锁定阈值 best_threshold ----------
  # Development
  save_metrics_bundle(
    prefix = "08_Development_DefaultThreshold_0.5",
    truth01 = dev_truth01,
    prob = dev_prob,
    threshold = 0.5,
    outdir = outdir
  )

  save_metrics_bundle(
    prefix = "09_Development_LockedThreshold_Youden",
    truth01 = dev_truth01,
    prob = dev_prob,
    threshold = best_threshold,
    outdir = outdir
  )

  # External（注意：这里用 development 锁定下来的 threshold，不要在 external 上重新调阈值）
  if (has_external) {
    save_metrics_bundle(
      prefix = "10_External_DefaultThreshold_0.5",
      truth01 = ext_truth01,
      prob = ext_prob,
      threshold = 0.5,
      outdir = outdir
    )

    save_metrics_bundle(
      prefix = "11_External_LockedThreshold_YoudenFromDevelopment",
      truth01 = ext_truth01,
      prob = ext_prob,
      threshold = best_threshold,
      outdir = outdir
    )
  }

  ## ---------- 8. 保存样本级预测（附加阈值分类） ----------
  dev_pred_df <- data.frame(
    truth = dev_truth01,
    prob = dev_prob,
    pred_0.5 = ifelse(dev_prob >= 0.5, 1L, 0L),
    pred_locked = ifelse(dev_prob >= best_threshold, 1L, 0L),
    stringsAsFactors = FALSE
  )
  utils::write.csv(
    dev_pred_df,
    file = file.path(outdir, "12_Development_Predictions_With_Thresholds.csv"),
    row.names = FALSE
  )

  if (has_external) {
    ext_pred_df <- data.frame(
      truth = ext_truth01,
      prob = ext_prob,
      pred_0.5 = ifelse(ext_prob >= 0.5, 1L, 0L),
      pred_locked = ifelse(ext_prob >= best_threshold, 1L, 0L),
      stringsAsFactors = FALSE
    )
    utils::write.csv(
      ext_pred_df,
      file = file.path(outdir, "13_External_Predictions_With_Thresholds.csv"),
      row.names = FALSE
    )
  }

  ## ---------- 9. 汇总报告 ----------
  summary_rows <- list(
    data.frame(
      Cohort = "development",
      AUC = tryCatch(as.numeric(pROC::auc(dev_truth01, dev_prob)), error = function(e) NA_real_),
      Brier = dev_brier,
      DefaultThreshold = 0.5,
      LockedThreshold = best_threshold,
      stringsAsFactors = FALSE
    )
  )

  if (has_external) {
    summary_rows[[2]] <- data.frame(
      Cohort = "external",
      AUC = tryCatch(as.numeric(pROC::auc(ext_truth01, ext_prob)), error = function(e) NA_real_),
      Brier = ext_brier,
      DefaultThreshold = 0.5,
      LockedThreshold = best_threshold,
      stringsAsFactors = FALSE
    )
  }

  summary_df <- do.call(rbind, summary_rows)
  utils::write.csv(
    summary_df,
    file = file.path(outdir, "14_Calibration_Threshold_Summary.csv"),
    row.names = FALSE
  )

  log_info("✅ Calibration + threshold optimization 已完成。结果已保存到: %s", outdir)

  metrics_df <- data.frame(
    Metric = c(
      "DevelopmentN",
      "DevelopmentLowRiskN",
      "DevelopmentHighRiskN",
      "InputMouseGeneN",
      "MappedHumanGeneN",
      "AvailableHumanGeneN",
      "SelectedFeatureN",
      "LambdaChoice",
      "LambdaValue",
      "DevelopmentApparentAUC",
      "ExternalN",
      "ExternalMatchedFeatureN",
      "ExternalAUC"
    ),
    Value = as.character(c(
      nrow(dat_dev),
      unname(dev_tab["low-risk"]),
      unname(dev_tab["high-risk"]),
      length(mouse_genes),
      length(human_genes),
      length(available_genes),
      length(selected_features),
      lambda_choice,
      lam,
      dev_auc,
      external_n,
      external_common_feature_n,
      external_auc
    )),
    stringsAsFactors = FALSE
  )
  write.csv(metrics_df, file.path(out_dir, "model_metrics.csv"), row.names = FALSE)
  
  save_session_info(file.path(out_dir, "session_info.txt"))

  ## =========================================================
  ## Logistic Recalibration
  ## 说明：
  ## - 在 development cohort 上拟合 logistic recalibration:
  ##     recalibrated_logit = a + b * original_logit
  ## - 再应用到 development / external
  ## - 输出重校准前后对比图、指标表、阈值对比
  ## =========================================================
  ## ---------- 0. 前提检查 ----------
  ## ---------- 0. 前提检查 ----------
  current_env <- environment()
  required_objs <- c("dev_truth01", "dev_prob", "outdir")
  missing_objs <- required_objs[!vapply(required_objs, function(x) exists(x, envir = current_env), logical(1))]
  if (length(missing_objs) > 0) {
    stop(sprintf(
      "❌ Logistic recalibration 缺少必要对象：%s。请先运行前面的 calibration + threshold 代码块。",
      paste(missing_objs, collapse = ", ")
    ), call. = FALSE)
  }

  has_external_recal <- exists("ext_truth01", envir = current_env) && exists("ext_prob", envir = current_env)

  ## ---------- 1. 工具函数 ----------
  clip_prob <- function(p, eps = 1e-6) {
    p <- as.numeric(p)
    p[p < eps] <- eps
    p[p > 1 - eps] <- 1 - eps
    p
  }

  fit_logistic_recalibration <- function(truth01, prob) {
    prob <- clip_prob(prob)
    lp <- qlogis(prob)

    df <- data.frame(
      y = as.integer(truth01),
      lp = lp
    )

    fit <- glm(y ~ lp, data = df, family = binomial())

    list(
      model = fit,
      intercept = unname(coef(fit)[1]),
      slope = unname(coef(fit)[2])
    )
  }

  apply_logistic_recalibration <- function(prob, intercept, slope) {
    prob <- clip_prob(prob)
    lp <- qlogis(prob)
    plogis(intercept + slope * lp)
  }

  if (!exists("calc_metrics_at_threshold", inherits = TRUE)) {
    calc_metrics_at_threshold <- function(truth01, prob, threshold) {
      pred01 <- ifelse(prob >= threshold, 1L, 0L)

      tn <- sum(truth01 == 0 & pred01 == 0, na.rm = TRUE)
      fp <- sum(truth01 == 0 & pred01 == 1, na.rm = TRUE)
      fn <- sum(truth01 == 1 & pred01 == 0, na.rm = TRUE)
      tp <- sum(truth01 == 1 & pred01 == 1, na.rm = TRUE)

      sensitivity <- ifelse((tp + fn) == 0, NA_real_, tp / (tp + fn))
      specificity <- ifelse((tn + fp) == 0, NA_real_, tn / (tn + fp))
      ppv <- ifelse((tp + fp) == 0, NA_real_, tp / (tp + fp))
      npv <- ifelse((tn + fn) == 0, NA_real_, tn / (tn + fn))
      accuracy <- (tp + tn) / (tp + tn + fp + fn)
      balanced_accuracy <- mean(c(sensitivity, specificity), na.rm = TRUE)
      f1 <- ifelse((2 * tp + fp + fn) == 0, NA_real_, 2 * tp / (2 * tp + fp + fn))
      brier <- mean((prob - truth01)^2, na.rm = TRUE)

      auc_val <- tryCatch({
        as.numeric(pROC::auc(pROC::roc(truth01, prob, quiet = TRUE, direction = "<")))
      }, error = function(e) NA_real_)

      data.frame(
        threshold = threshold,
        AUC = auc_val,
        Brier = brier,
        Accuracy = accuracy,
        Sensitivity = sensitivity,
        Specificity = specificity,
        PPV = ppv,
        NPV = npv,
        BalancedAccuracy = balanced_accuracy,
        F1 = f1,
        TN = tn,
        FP = fp,
        FN = fn,
        TP = tp,
        stringsAsFactors = FALSE
      )
    }
  }

  if (!exists("build_calibration_table", inherits = TRUE)) {
    build_calibration_table <- function(truth01, prob, n_bins = 10) {
      n <- length(prob)
      if (n < 40) {
        n_bins <- min(5, n)
      } else if (n < 80) {
        n_bins <- min(8, n)
      }

      qs <- unique(stats::quantile(prob, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE))
      if (length(qs) <= 2) {
        qs <- seq(min(prob, na.rm = TRUE), max(prob, na.rm = TRUE), length.out = n_bins + 1)
        qs <- unique(qs)
      }

      if (length(qs) <= 2) {
        bin <- rep(1L, length(prob))
      } else {
        bin <- cut(prob, breaks = qs, include.lowest = TRUE, labels = FALSE)
      }

      df <- data.frame(truth01 = truth01, prob = prob, bin = bin)
      df <- df[!is.na(df$bin), , drop = FALSE]

      out <- do.call(
        rbind,
        lapply(split(df, df$bin), function(z) {
          n_bin <- nrow(z)
          obs <- mean(z$truth01)
          pred <- mean(z$prob)
          se <- sqrt(obs * (1 - obs) / max(n_bin, 1))
          lower <- max(0, obs - 1.96 * se)
          upper <- min(1, obs + 1.96 * se)

          data.frame(
            bin = unique(z$bin),
            n = n_bin,
            mean_pred = pred,
            obs_rate = obs,
            obs_lower = lower,
            obs_upper = upper
          )
        })
      )
      rownames(out) <- NULL
      out
    }
  }

  find_best_threshold_youden_quiet <- function(truth01, prob) {
    roc_obj <- pROC::roc(truth01, prob, quiet = TRUE, direction = "<")
    coords <- pROC::coords(
      roc_obj,
      x = "best",
      best.method = "youden",
      transpose = FALSE,
      ret = c("threshold", "sensitivity", "specificity", "accuracy")
    )
    as.data.frame(coords, stringsAsFactors = FALSE)
  }

  plot_calibration_compare_pdf <- function(calib_before, calib_after, title_before, title_after, file) {
    grDevices::pdf(file, width = 12, height = 5.8)
    op <- par(no.readonly = TRUE)
    on.exit({
      par(op)
      grDevices::dev.off()
    }, add = TRUE)

    par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))

    # before
    plot(
      calib_before$mean_pred, calib_before$obs_rate,
      xlim = c(0, 1), ylim = c(0, 1),
      xlab = "Mean predicted probability",
      ylab = "Observed event rate",
      main = title_before,
      pch = 19, cex = 1.1, col = "#d62728"
    )
    abline(0, 1, lty = 2, lwd = 2, col = "gray50")
    segments(
      x0 = calib_before$mean_pred, y0 = calib_before$obs_lower,
      x1 = calib_before$mean_pred, y1 = calib_before$obs_upper,
      col = "#d62728", lwd = 1.4
    )
    ord1 <- order(calib_before$mean_pred)
    lines(calib_before$mean_pred[ord1], calib_before$obs_rate[ord1], col = "#d62728", lwd = 2)
    text(calib_before$mean_pred, calib_before$obs_rate,
        labels = paste0("n=", calib_before$n), pos = 3, cex = 0.75)

    # after
    plot(
      calib_after$mean_pred, calib_after$obs_rate,
      xlim = c(0, 1), ylim = c(0, 1),
      xlab = "Mean predicted probability",
      ylab = "Observed event rate",
      main = title_after,
      pch = 19, cex = 1.1, col = "#1f77b4"
    )
    abline(0, 1, lty = 2, lwd = 2, col = "gray50")
    segments(
      x0 = calib_after$mean_pred, y0 = calib_after$obs_lower,
      x1 = calib_after$mean_pred, y1 = calib_after$obs_upper,
      col = "#1f77b4", lwd = 1.4
    )
    ord2 <- order(calib_after$mean_pred)
    lines(calib_after$mean_pred[ord2], calib_after$obs_rate[ord2], col = "#1f77b4", lwd = 2)
    text(calib_after$mean_pred, calib_after$obs_rate,
        labels = paste0("n=", calib_after$n), pos = 3, cex = 0.75)
  }

  ## ---------- 2. 在 development 上拟合重校准模型 ----------
  recal_fit <- fit_logistic_recalibration(dev_truth01, dev_prob)
  recal_intercept <- recal_fit$intercept
  recal_slope <- recal_fit$slope

  if (exists("log_info", inherits = TRUE)) {
    log_info("Logistic recalibration 拟合完成：intercept = %.6f, slope = %.6f",
            recal_intercept, recal_slope)
  } else {
    message(sprintf(
      "Logistic recalibration 拟合完成：intercept = %.6f, slope = %.6f",
      recal_intercept, recal_slope
    ))
  }

  recal_coef_df <- data.frame(
    term = c("intercept", "slope"),
    estimate = c(recal_intercept, recal_slope),
    stringsAsFactors = FALSE
  )
  utils::write.csv(
    recal_coef_df,
    file = file.path(outdir, "15_Logistic_Recalibration_Model.csv"),
    row.names = FALSE
  )

  ## ---------- 3. 应用重校准 ----------
  dev_prob_recal <- apply_logistic_recalibration(dev_prob, recal_intercept, recal_slope)

  if (has_external_recal) {
    ext_prob_recal <- apply_logistic_recalibration(ext_prob, recal_intercept, recal_slope)
  }

  ## ---------- 4. 重校准前后 calibration table ----------
  dev_calib_before_recal <- build_calibration_table(dev_truth01, dev_prob, n_bins = 10)
  dev_calib_after_recal  <- build_calibration_table(dev_truth01, dev_prob_recal, n_bins = 10)

  dev_calib_compare_df <- rbind(
    cbind(Stage = "Before_Recalibration", dev_calib_before_recal),
    cbind(Stage = "After_Recalibration",  dev_calib_after_recal)
  )
  utils::write.csv(
    dev_calib_compare_df,
    file = file.path(outdir, "16_Development_Calibration_BeforeAfter_Table.csv"),
    row.names = FALSE
  )

  plot_calibration_compare_pdf(
    calib_before = dev_calib_before_recal,
    calib_after  = dev_calib_after_recal,
    title_before = "Development Calibration (Before)",
    title_after  = "Development Calibration (After)",
    file = file.path(outdir, "16_Development_Calibration_BeforeAfter_Plot.pdf")
  )

  if (has_external_recal) {
    ext_calib_before_recal <- build_calibration_table(ext_truth01, ext_prob, n_bins = 10)
    ext_calib_after_recal  <- build_calibration_table(ext_truth01, ext_prob_recal, n_bins = 10)

    ext_calib_compare_df <- rbind(
      cbind(Stage = "Before_Recalibration", ext_calib_before_recal),
      cbind(Stage = "After_Recalibration",  ext_calib_after_recal)
    )
    utils::write.csv(
      ext_calib_compare_df,
      file = file.path(outdir, "17_External_Calibration_BeforeAfter_Table.csv"),
      row.names = FALSE
    )

    plot_calibration_compare_pdf(
      calib_before = ext_calib_before_recal,
      calib_after  = ext_calib_after_recal,
      title_before = "External Calibration (Before)",
      title_after  = "External Calibration (After)",
      file = file.path(outdir, "17_External_Calibration_BeforeAfter_Plot.pdf")
    )
  }

  ## ---------- 5. 重校准前后阈值优化（仍然只在 development 上做） ----------
  orig_thr_df <- find_best_threshold_youden_quiet(dev_truth01, dev_prob)
  orig_thr <- as.numeric(orig_thr_df$threshold[1])

  recal_thr_df <- find_best_threshold_youden_quiet(dev_truth01, dev_prob_recal)
  recal_thr <- as.numeric(recal_thr_df$threshold[1])

  thr_compare_df <- rbind(
    data.frame(
      Model = "Original",
      threshold = orig_thr,
      sensitivity = as.numeric(orig_thr_df$sensitivity[1]),
      specificity = as.numeric(orig_thr_df$specificity[1]),
      accuracy = as.numeric(orig_thr_df$accuracy[1]),
      stringsAsFactors = FALSE
    ),
    data.frame(
      Model = "Recalibrated",
      threshold = recal_thr,
      sensitivity = as.numeric(recal_thr_df$sensitivity[1]),
      specificity = as.numeric(recal_thr_df$specificity[1]),
      accuracy = as.numeric(recal_thr_df$accuracy[1]),
      stringsAsFactors = FALSE
    )
  )

  utils::write.csv(
    thr_compare_df,
    file = file.path(outdir, "18_Recalibration_Threshold_Comparison.csv"),
    row.names = FALSE
  )

  ## ---------- 6. 指标表：重校准前后对比 ----------
  make_metrics_compare <- function(truth01, prob_before, prob_after,
                                  thr_before, thr_after, cohort_name) {
    m1 <- calc_metrics_at_threshold(truth01, prob_before, 0.5)
    m1$Cohort <- cohort_name
    m1$Model <- "Original"
    m1$ThresholdType <- "Default_0.5"

    m2 <- calc_metrics_at_threshold(truth01, prob_before, thr_before)
    m2$Cohort <- cohort_name
    m2$Model <- "Original"
    m2$ThresholdType <- "Locked_Youden_From_Development"

    m3 <- calc_metrics_at_threshold(truth01, prob_after, 0.5)
    m3$Cohort <- cohort_name
    m3$Model <- "Recalibrated"
    m3$ThresholdType <- "Default_0.5"

    m4 <- calc_metrics_at_threshold(truth01, prob_after, thr_after)
    m4$Cohort <- cohort_name
    m4$Model <- "Recalibrated"
    m4$ThresholdType <- "Locked_Youden_From_Development"

    out <- rbind(m1, m2, m3, m4)

    out <- out[, c(
      "Cohort", "Model", "ThresholdType", "threshold",
      "AUC", "Brier", "Accuracy", "Sensitivity", "Specificity",
      "PPV", "NPV", "BalancedAccuracy", "F1",
      "TN", "FP", "FN", "TP"
    )]

    rownames(out) <- NULL
    out
  }

  dev_metrics_compare <- make_metrics_compare(
    truth01 = dev_truth01,
    prob_before = dev_prob,
    prob_after = dev_prob_recal,
    thr_before = orig_thr,
    thr_after = recal_thr,
    cohort_name = "development"
  )

  utils::write.csv(
    dev_metrics_compare,
    file = file.path(outdir, "19_Development_Metrics_BeforeAfter_Recalibration.csv"),
    row.names = FALSE
  )

  if (has_external_recal) {
    ext_metrics_compare <- make_metrics_compare(
      truth01 = ext_truth01,
      prob_before = ext_prob,
      prob_after = ext_prob_recal,
      thr_before = orig_thr,
      thr_after = recal_thr,
      cohort_name = "external"
    )

    utils::write.csv(
      ext_metrics_compare,
      file = file.path(outdir, "20_External_Metrics_BeforeAfter_Recalibration.csv"),
      row.names = FALSE
    )
  }

  ## ---------- 7. 保存样本级概率 ----------
  dev_pred_recal_df <- data.frame(
    truth = dev_truth01,
    prob_original = dev_prob,
    prob_recalibrated = dev_prob_recal,
    pred_original_0.5 = ifelse(dev_prob >= 0.5, 1L, 0L),
    pred_original_locked = ifelse(dev_prob >= orig_thr, 1L, 0L),
    pred_recal_0.5 = ifelse(dev_prob_recal >= 0.5, 1L, 0L),
    pred_recal_locked = ifelse(dev_prob_recal >= recal_thr, 1L, 0L),
    stringsAsFactors = FALSE
  )
  utils::write.csv(
    dev_pred_recal_df,
    file = file.path(outdir, "21_Development_Predictions_Original_vs_Recalibrated.csv"),
    row.names = FALSE
  )

  if (has_external_recal) {
    ext_pred_recal_df <- data.frame(
      truth = ext_truth01,
      prob_original = ext_prob,
      prob_recalibrated = ext_prob_recal,
      pred_original_0.5 = ifelse(ext_prob >= 0.5, 1L, 0L),
      pred_original_locked = ifelse(ext_prob >= orig_thr, 1L, 0L),
      pred_recal_0.5 = ifelse(ext_prob_recal >= 0.5, 1L, 0L),
      pred_recal_locked = ifelse(ext_prob_recal >= recal_thr, 1L, 0L),
      stringsAsFactors = FALSE
    )
    utils::write.csv(
      ext_pred_recal_df,
      file = file.path(outdir, "22_External_Predictions_Original_vs_Recalibrated.csv"),
      row.names = FALSE
    )
  }

  ## ---------- 8. 汇总表 ----------
  summary_list_recal <- list(
    data.frame(
      Cohort = "development",
      Original_AUC = tryCatch(as.numeric(pROC::auc(pROC::roc(dev_truth01, dev_prob, quiet = TRUE, direction = "<"))), error = function(e) NA_real_),
      Recalibrated_AUC = tryCatch(as.numeric(pROC::auc(pROC::roc(dev_truth01, dev_prob_recal, quiet = TRUE, direction = "<"))), error = function(e) NA_real_),
      Original_Brier = mean((dev_prob - dev_truth01)^2, na.rm = TRUE),
      Recalibrated_Brier = mean((dev_prob_recal - dev_truth01)^2, na.rm = TRUE),
      Original_Youden_Threshold = orig_thr,
      Recalibrated_Youden_Threshold = recal_thr,
      stringsAsFactors = FALSE
    )
  )

  if (has_external_recal) {
    summary_list_recal[[2]] <- data.frame(
      Cohort = "external",
      Original_AUC = tryCatch(as.numeric(pROC::auc(pROC::roc(ext_truth01, ext_prob, quiet = TRUE, direction = "<"))), error = function(e) NA_real_),
      Recalibrated_AUC = tryCatch(as.numeric(pROC::auc(pROC::roc(ext_truth01, ext_prob_recal, quiet = TRUE, direction = "<"))), error = function(e) NA_real_),
      Original_Brier = mean((ext_prob - ext_truth01)^2, na.rm = TRUE),
      Recalibrated_Brier = mean((ext_prob_recal - ext_truth01)^2, na.rm = TRUE),
      Original_Youden_Threshold = orig_thr,
      Recalibrated_Youden_Threshold = recal_thr,
      stringsAsFactors = FALSE
    )
  }

  summary_recal_df <- do.call(rbind, summary_list_recal)
  utils::write.csv(
    summary_recal_df,
    file = file.path(outdir, "23_Logistic_Recalibration_Summary.csv"),
    row.names = FALSE
  )

  ## ---------- 9. 日志 ----------
  if (exists("log_info", inherits = TRUE)) {
    log_info("✅ Logistic recalibration 已完成。输出文件：15~23 前缀结果已保存到: %s", outdir)
    log_info("重校准参数：intercept = %.6f, slope = %.6f", recal_intercept, recal_slope)
  } else {
    message(sprintf("✅ Logistic recalibration 已完成，输出目录: %s", outdir))
    message(sprintf("重校准参数：intercept = %.6f, slope = %.6f", recal_intercept, recal_slope))
  }

  log_info("✅ 人类验证模块完成。输出目录: %s", out_dir)
}
  
tryCatch(
  main(),
  error = function(e) {
    cat(sprintf("[%s] [ERROR] %s\n", Sys.time(), e$message))
    quit(save = "no", status = 1)
  }
)
