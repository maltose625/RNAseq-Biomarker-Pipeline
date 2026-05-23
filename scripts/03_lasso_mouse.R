#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(yaml)
  library(DESeq2)
  library(glmnet)
  library(data.table)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
}

log_info <- function(...) cat(sprintf("[%s] [INFO] %s\n", Sys.time(), sprintf(...)))
log_warn <- function(...) cat(sprintf("[%s] [WARN] %s\n", Sys.time(), sprintf(...)))
log_stop <- function(...) stop(sprintf("[%s] [ERROR] %s", Sys.time(), sprintf(...)))

read_config <- function() {
  cfg_path <- "config/config.yaml"
  if (!file.exists(cfg_path)) log_stop("未找到配置文件: %s", cfg_path)
  yaml::read_yaml(cfg_path)
}

read_count_matrix <- function(path, strip_ensembl_version = TRUE) {
  if (!file.exists(path)) log_stop("未找到 counts_matrix: %s", path)
  log_info("读取 counts_matrix: %s", path)

  # featureCounts 输出第一行为注释，真正表头从第二行开始
  df <- data.table::fread(
    path,
    sep = "\t",
    header = TRUE,
    skip = 1,
    data.table = FALSE,
    check.names = FALSE,
    quote = "",
    fill = TRUE
  )

  colnames(df) <- sub("^\ufeff", "", colnames(df))

  required_cols <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")
  if (!all(required_cols %in% colnames(df))) {
    log_stop(
      "counts_matrix 不是预期的 featureCounts 格式。当前列名: %s",
      paste(colnames(df), collapse = ", ")
    )
  }

  if (ncol(df) < 7) {
    log_stop("counts_matrix 列数异常，当前仅 %d 列，无法识别样本计数列", ncol(df))
  }

  gene_ids <- as.character(df$Geneid)
  if (strip_ensembl_version) {
    gene_ids <- sub("\\..*$", "", gene_ids)
  }

  rownames(df) <- make.unique(gene_ids)

  # 仅保留样本计数列（第 7 列开始）
  sample_df <- df[, -(1:6), drop = FALSE]

  # 清洗样本名：完整 bam 路径 -> SRR8581315
  clean_names <- basename(colnames(sample_df))
  clean_names <- sub("_sorted\\.bam$", "", clean_names)
  clean_names <- sub("\\.bam$", "", clean_names)
  colnames(sample_df) <- clean_names

  mat <- as.matrix(sample_df)
  storage.mode(mat) <- "double"

  if (anyNA(mat)) {
    log_warn("counts_matrix 中存在 NA，已保留；后续若出错请检查原始文件")
  }

  log_info("counts_matrix 读取完成：%d genes × %d samples", nrow(mat), ncol(mat))
  log_info("样本名: %s", paste(colnames(mat), collapse = ", "))

  mat
}

read_sample_info <- function(path, ref_group, trt_group) {
  if (!file.exists(path)) log_stop("未找到 sample_info: %s", path)
  log_info("读取 sample_info: %s", path)

  meta <- data.table::fread(
    path,
    sep = "\t",
    header = TRUE,
    data.table = FALSE,
    check.names = FALSE,
    quote = "",
    fill = TRUE
  )

  colnames(meta) <- sub("^\ufeff", "", colnames(meta))

  sample_candidates <- c("sample_id", "sample", "Sample", "SampleID", "run", "Run", "id", "ID")
  group_candidates  <- c("Condition", "condition", "group", "Group", "treatment", "Treatment")

  sample_col <- sample_candidates[sample_candidates %in% colnames(meta)][1]
  group_col  <- group_candidates[group_candidates %in% colnames(meta)][1]

  if (is.na(sample_col)) {
    log_stop("sample_info 缺少样本列。当前列名: %s", paste(colnames(meta), collapse = ", "))
  }
  if (is.na(group_col)) {
    log_stop("sample_info 缺少分组列。当前列名: %s", paste(colnames(meta), collapse = ", "))
  }

  meta <- meta[, c(sample_col, group_col), drop = FALSE]
  colnames(meta) <- c("Sample", "Condition")

  meta$Sample <- trimws(as.character(meta$Sample))
  meta$Condition <- trimws(as.character(meta$Condition))

  meta <- meta[meta$Condition %in% c(ref_group, trt_group), , drop = FALSE]
  if (nrow(meta) == 0) {
    log_stop("sample_info 中未找到 ref_group=%s 或 trt_group=%s 对应样本", ref_group, trt_group)
  }

  rownames(meta) <- meta$Sample
  meta$Condition <- factor(meta$Condition, levels = c(ref_group, trt_group))

  log_info("sample_info 读取完成：%d samples", nrow(meta))
  log_info("分组统计: %s", paste(names(table(meta$Condition)), table(meta$Condition), collapse = ", "))

  meta
}

read_deg_results <- function(path, strip_ensembl_version = TRUE) {
  if (!file.exists(path)) log_stop("未找到 DESeq2 结果文件: %s", path)

  deg <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)

  if (!"padj" %in% colnames(deg)) log_stop("DESeq2 结果缺少 padj 列")
  if (!"log2FoldChange" %in% colnames(deg)) log_stop("DESeq2 结果缺少 log2FoldChange 列")

  if ("Ensembl" %in% colnames(deg)) {
    deg$GeneKey <- as.character(deg$Ensembl)
  } else if ("Symbol" %in% colnames(deg)) {
    deg$GeneKey <- as.character(deg$Symbol)
  } else {
    deg$GeneKey <- as.character(deg[[1]])
  }

  if (strip_ensembl_version) {
    deg$GeneKey <- sub("\\..*$", "", deg$GeneKey)
  }

  deg
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

save_session_info <- function(out_file) {
  sink(out_file)
  print(sessionInfo())
  sink()
}

main <- function() {
  log_info("🚀 启动小鼠特征筛选模块...")

  config <- read_config()

  out_root <- config$paths$results_ml %||% "results/03_MachineLearning"
  out_dir  <- file.path(out_root, "mouse_feature_selection")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  ref_group <- config$params$deseq2$ref_group %||% "Control"
  trt_group <- config$params$deseq2$trt_group %||% "Treatment"

  strip_ensembl_version <- as.logical(config$params$deseq2$strip_ensembl_version %||% TRUE)
  min_count   <- as.integer(config$params$deseq2$min_count %||% 10)
  min_samples <- as.integer(config$params$deseq2$min_samples %||% 2)

  p_cutoff  <- as.numeric(config$thresholds$padj_cutoff %||% 0.05)
  fc_cutoff <- as.numeric(config$thresholds$log2fc_cutoff %||% 0.58)

  seed <- as.integer(config$params$ml$seed %||% 123)
  alpha <- as.numeric(config$params$ml$alpha %||% 1)
  lambda_choice <- config$params$ml$lambda_choice %||% "lambda.1se"
  standardize <- as.logical(config$params$ml$standardize %||% TRUE)

  # 这里我建议阈值设得更保守一点；小样本时不硬做 CV
  min_mouse_cv <- as.integer(config$params$ml$min_class_count_mouse_cv %||% 4)

  counts_file <- config$paths$counts_matrix
  sample_file <- config$paths$sample_info
  deg_file    <- config$paths$deseq2_result_file

  counts <- read_count_matrix(counts_file, strip_ensembl_version = strip_ensembl_version)
  meta   <- read_sample_info(sample_file, ref_group, trt_group)

  common_samples <- intersect(colnames(counts), meta$Sample)
  log_info("共同样本数: %d", length(common_samples))

  if (length(common_samples) < 4) {
    log_stop("共同样本数过少: %d。请检查 counts 列名与 sample_info 样本名是否一致", length(common_samples))
  }

  counts <- counts[, common_samples, drop = FALSE]
  meta   <- meta[common_samples, , drop = FALSE]

  # 低表达过滤
  keep <- rowSums(counts >= min_count) >= min_samples
  counts <- counts[keep, , drop = FALSE]
  log_info("低表达过滤后保留基因数: %d", nrow(counts))

  dds <- DESeqDataSetFromMatrix(
    countData = round(counts),
    colData = meta,
    design = ~ Condition
  )

  log_info("正在进行 VST 转换...")
  vst_mat <- assay(vst(dds, blind = TRUE))

  if (strip_ensembl_version) {
    rownames(vst_mat) <- sub("\\..*$", "", rownames(vst_mat))
  }

  deg <- read_deg_results(deg_file, strip_ensembl_version = strip_ensembl_version)

  candidate_genes <- unique(deg$GeneKey[
    !is.na(deg$padj) &
      !is.na(deg$log2FoldChange) &
      deg$padj < p_cutoff &
      abs(deg$log2FoldChange) > fc_cutoff
  ])

  candidate_genes <- intersect(candidate_genes, rownames(vst_mat))
  log_info("候选差异基因数: %d", length(candidate_genes))

  if (length(candidate_genes) < 2) {
    log_stop("候选差异基因不足 2 个，无法进行特征筛选")
  }

  x <- t(vst_mat[candidate_genes, , drop = FALSE])
  x <- as.matrix(x)
  storage.mode(x) <- "double"

  y <- factor(meta$Condition, levels = c(ref_group, trt_group))
  class_tab <- table(y)
  log_info("类别分布: %s", paste(names(class_tab), class_tab, collapse = ", "))

  if (length(class_tab) != 2) log_stop("当前不是二分类任务")
  if (min(class_tab) < 2) log_stop("某一类样本数 < 2，无法运行 binomial glmnet")

  # 小样本自动降级：只导出候选基因，不硬做 CV
  if (min(class_tab) < min_mouse_cv) {
    log_warn("最小类别样本数 < %d，跳过 LASSO CV，仅导出候选差异基因", min_mouse_cv)

    old_core_files <- c(
      file.path(out_dir, "Core_Biomarkers_1se.csv"),
      file.path(out_dir, "Core_Biomarkers_min.csv"),
      file.path(out_dir, "lasso_coefficients.csv"),
      file.path(out_dir, "lasso_coefficients_lambda_1se.csv"),
      file.path(out_dir, "lasso_coefficients_lambda_min.csv"),
      file.path(out_dir, "01_LASSO_CV_Curve.pdf"),
      file.path(out_dir, "NO_SELECTED_FEATURES.txt")
     )

    old_core_files <- old_core_files[file.exists(old_core_files)]
    if (length(old_core_files) > 0) {
    file.remove(old_core_files)
    log_info("已清理旧的 LASSO 核心结果文件，避免与本次降级输出混淆")
    }

    write.csv(
      data.frame(MouseGene = candidate_genes),
      file.path(out_dir, "Candidate_Genes_for_ML.csv"),
      row.names = FALSE
    )

    writeLines(
    c(
      "LASSO was skipped due to insufficient sample size or class imbalance.",
      sprintf("MinClassN = %d", min(class_tab)),
      sprintf("CandidateGeneN = %d", length(candidate_genes)),
      "WARNING: Candidate_Genes_for_ML.csv is NOT equivalent to final core biomarkers."
    ),
    file.path(out_dir, "ML_SKIP_REASON.txt")
    )

    write.csv(
      data.frame(Sample = rownames(meta), Condition = meta$Condition),
      file.path(out_dir, "samples_used.csv"),
      row.names = FALSE
    )

    write.csv(vst_mat, file.path(out_dir, "mouse_vst_matrix.csv"))
    save_session_info(file.path(out_dir, "session_info.txt"))

    log_info("✅ 已完成降级输出。输出目录: %s", out_dir)
    return(invisible(NULL))
  }

  # 仅在样本稍微够用时做 LASSO CV
  k <- min(3, min(class_tab))
  foldid <- make_stratified_foldid(y, k = k, seed = seed)

  log_info("正在执行 cv.glmnet (k=%d)...", k)
  set.seed(seed)
  cvfit <- cv.glmnet(
    x = x,
    y = y,
    family = "binomial",
    alpha = alpha,
    standardize = standardize,
    nfolds = k,
    foldid = foldid,
    type.measure = "class"
  )

  lam <- if (lambda_choice == "lambda.min") cvfit$lambda.min else cvfit$lambda.1se
  coef_mat <- as.matrix(coef(cvfit, s = lam))
  selected <- rownames(coef_mat)[coef_mat[, 1] != 0]
  selected <- setdiff(selected, "(Intercept)")

  log_info("在 %s 下筛到 %d 个候选标志物。", lambda_choice, length(selected))

  core_file <- if (lambda_choice == "lambda.min") {
    "Core_Biomarkers_min.csv"
  } else {
    "Core_Biomarkers_1se.csv"
  }

  coef_file <- if (lambda_choice == "lambda.min") {
    "lasso_coefficients_lambda_min.csv"
  } else {
    "lasso_coefficients_lambda_1se.csv"
  }

  write.csv(
    data.frame(MouseGene = candidate_genes),
    file.path(out_dir, "Candidate_Genes_for_ML.csv"),
    row.names = FALSE
  )

  write.csv(
    data.frame(MouseGene = selected),
    file.path(out_dir, core_file),
    row.names = FALSE
  )

  write.csv(
    data.frame(Feature = rownames(coef_mat), Coef = coef_mat[, 1]),
    file.path(out_dir, coef_file),
    row.names = FALSE
  )

  if (length(selected) == 0) {
    writeLines(
      c(
        sprintf("No non-zero features selected under %s.", lambda_choice),
        sprintf("CandidateGeneN = %d", length(candidate_genes))
      ),
      file.path(out_dir, "NO_SELECTED_FEATURES.txt")
    )
  }

  write.csv(
    data.frame(Sample = rownames(meta), Condition = meta$Condition),
    file.path(out_dir, "samples_used.csv"),
    row.names = FALSE
  )

  write.csv(vst_mat, file.path(out_dir, "mouse_vst_matrix.csv"))

  pdf(file.path(out_dir, "01_LASSO_CV_Curve.pdf"), width = 7, height = 5)
  plot(cvfit, main = sprintf("Mouse LASSO CV (%s)", lambda_choice))
  dev.off()

  save_session_info(file.path(out_dir, "session_info.txt"))
  log_info("✅ 小鼠特征筛选完成。输出目录: %s", out_dir)
}

tryCatch(
  main(),
  error = function(e) {
    cat(sprintf("[%s] [ERROR] %s\n", Sys.time(), e$message))
    quit(save = "no", status = 1)
  }
)
