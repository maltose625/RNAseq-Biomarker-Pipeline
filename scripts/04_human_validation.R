#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(yaml)
  library(data.table)
  library(glmnet)
  library(pROC)
  library(ggplot2)
})

source("scripts/04_1_hv_common.R")
source("scripts/04_2_hv_mapping_io.R")
source("scripts/04_3_hv_cohort_io.R")
source("scripts/04_4_hv_modeling.R")
source("scripts/04_5_hv_calibration_threshold.R")
source("scripts/04_6_hv_recalibration.R")

# =========================
# Local helper functions
# =========================

get_cvfit_from_result <- function(model_res) {
  if (inherits(model_res, "cv.glmnet")) return(model_res)
  if (is.list(model_res) && !is.null(model_res$cvfit)) return(model_res$cvfit)
  if (is.list(model_res) && !is.null(model_res$model)) return(model_res$model)
  log_stop("fit_development_model() 返回对象中未找到 cv.glmnet 模型（cvfit/model）")
}

safe_auc <- function(y01, prob, label = "unknown") {
  y01 <- as.integer(y01)
  prob <- as.numeric(prob)

  keep <- is.finite(y01) & is.finite(prob)
  y01 <- y01[keep]
  prob <- prob[keep]

  n_pos <- sum(y01 == 1)
  n_neg <- sum(y01 == 0)

  log_info("%s AUC 输入检查: n=%d, pos=%d, neg=%d", label, length(y01), n_pos, n_neg)

  if (length(unique(y01)) < 2 || n_pos == 0 || n_neg == 0) {
    log_warn("%s AUC 跳过：仅存在单一类别，返回 NA", label)
    return(NA_real_)
  }

  if (length(unique(round(prob, 12))) < 2) {
    log_warn("%s AUC 提示：预测概率近似常数，AUC 可能退化到 0.5", label)
  }

  roc_obj <- pROC::roc(response = y01, predictor = prob, quiet = TRUE, direction = "<")
  as.numeric(pROC::auc(roc_obj))
}

save_roc_curve_pdf <- function(y01, prob, file, title_text) {
  y01 <- as.integer(y01)
  prob <- as.numeric(prob)
  keep <- is.finite(y01) & is.finite(prob)
  y01 <- y01[keep]
  prob <- prob[keep]

  if (length(unique(y01)) < 2) {
    log_warn("ROC 曲线跳过：%s 仅存在单一类别", title_text)
    return(invisible(NULL))
  }

  roc_obj <- pROC::roc(response = y01, predictor = prob, quiet = TRUE, direction = "<")
  auc_val <- as.numeric(pROC::auc(roc_obj))

  pdf(file, width = 6, height = 6)
  plot(
    roc_obj,
    col = "#2C7FB8",
    lwd = 2,
    main = sprintf("%s\nAUC = %.4f", title_text, auc_val),
    legacy.axes = TRUE
  )
  abline(a = 0, b = 1, lty = 2, col = "grey70")
  dev.off()
}

extract_coef_df <- function(cvfit, lambda_selected) {
  coef_mat <- as.matrix(stats::coef(cvfit, s = lambda_selected))
  data.frame(
    Feature = rownames(coef_mat),
    Coefficient = as.numeric(coef_mat[, 1]),
    stringsAsFactors = FALSE
  )
}

extract_selected_features <- function(cvfit, lambda_selected) {
  coef_df <- extract_coef_df(cvfit, lambda_selected = lambda_selected)
  coef_df <- coef_df[coef_df$Feature != "(Intercept)", , drop = FALSE]
  coef_df <- coef_df[coef_df$Coefficient != 0, , drop = FALSE]
  coef_df[order(abs(coef_df$Coefficient), decreasing = TRUE), , drop = FALSE]
}

clip_prob_local <- function(p, eps = 1e-6) {
  pmin(pmax(as.numeric(p), eps), 1 - eps)
}

brier_score_local <- function(prob, y01) {
  prob <- clip_prob_local(prob)
  mean((prob - y01) ^ 2)
}

make_confusion_df_local <- function(truth01, pred01) {
  truth_lab <- ifelse(truth01 == 1, "high-risk", "low-risk")
  pred_lab  <- ifelse(pred01 == 1, "high-risk", "low-risk")
  tb <- table(
    Truth = factor(truth_lab, levels = c("low-risk", "high-risk")),
    Predicted = factor(pred_lab, levels = c("low-risk", "high-risk"))
  )
  as.data.frame.matrix(tb)
}

calc_metrics_at_threshold_local <- function(prob, y01, threshold) {
  prob <- clip_prob_local(prob)
  pred01 <- ifelse(prob >= threshold, 1L, 0L)

  tp <- sum(pred01 == 1 & y01 == 1)
  tn <- sum(pred01 == 0 & y01 == 0)
  fp <- sum(pred01 == 1 & y01 == 0)
  fn <- sum(pred01 == 0 & y01 == 1)

  sensitivity <- ifelse((tp + fn) == 0, NA, tp / (tp + fn))
  specificity <- ifelse((tn + fp) == 0, NA, tn / (tn + fp))
  accuracy <- (tp + tn) / length(y01)
  ppv <- ifelse((tp + fp) == 0, NA, tp / (tp + fp))
  npv <- ifelse((tn + fn) == 0, NA, tn / (tn + fn))
  bal_acc <- mean(c(sensitivity, specificity), na.rm = TRUE)
  f1 <- ifelse(is.na(ppv) || is.na(sensitivity) || (ppv + sensitivity) == 0,
               NA,
               2 * ppv * sensitivity / (ppv + sensitivity))

  list(
    threshold = as.numeric(threshold),
    Accuracy = accuracy,
    Sensitivity = sensitivity,
    Specificity = specificity,
    PPV = ppv,
    NPV = npv,
    BalancedAccuracy = bal_acc,
    F1 = f1,
    TN = tn,
    FP = fp,
    FN = fn,
    TP = tp,
    pred01 = pred01
  )
}

find_best_threshold_youden_local <- function(prob, y01) {
  prob <- as.numeric(prob)
  y01 <- as.integer(y01)

  if (length(unique(round(prob, 12))) < 2) {
    log_warn("预测概率近似常数，Youden 阈值优化退化，回退为默认阈值 0.5")
    metrics <- calc_metrics_at_threshold_local(prob, y01, 0.5)
    return(data.frame(
      threshold = 0.5,
      sensitivity = metrics$Sensitivity,
      specificity = metrics$Specificity,
      accuracy = metrics$Accuracy,
      stringsAsFactors = FALSE
    ))
  }

  roc_obj <- pROC::roc(response = y01, predictor = prob, quiet = TRUE, direction = "<")
  co <- pROC::coords(
    roc_obj,
    x = "best",
    best.method = "youden",
    ret = c("threshold", "sensitivity", "specificity"),
    transpose = FALSE
  )

  get_scalar <- function(obj, name) {
    if (is.data.frame(obj)) return(as.numeric(obj[[name]][1]))
    if (is.list(obj) && !is.null(obj[[name]])) return(as.numeric(obj[[name]][1]))
    if (!is.null(names(obj)) && name %in% names(obj)) return(as.numeric(obj[name][[1]]))
    stop(sprintf("无法从 coords() 结果中提取字段: %s", name))
  }

  threshold   <- get_scalar(co, "threshold")
  sensitivity <- get_scalar(co, "sensitivity")
  specificity <- get_scalar(co, "specificity")

  if (!is.finite(threshold)) {
    log_warn("Youden 阈值为非有限值（%s），回退为 0.5", as.character(threshold))
    threshold <- 0.5
  }

  metrics <- calc_metrics_at_threshold_local(prob, y01, threshold)

  data.frame(
    threshold = threshold,
    sensitivity = sensitivity,
    specificity = specificity,
    accuracy = metrics$Accuracy,
    stringsAsFactors = FALSE
  )
}

build_calibration_table_local <- function(prob, y01, n_bins = 8) {
  prob <- clip_prob_local(prob)
  qs <- unique(as.numeric(stats::quantile(
    prob,
    probs = seq(0, 1, length.out = n_bins + 1),
    na.rm = TRUE
  )))

  if (length(qs) < 3) {
    qs <- seq(0, 1, length.out = n_bins + 1)
  }

  bin <- cut(prob, breaks = qs, include.lowest = TRUE, labels = FALSE)
  df <- data.frame(prob = prob, y01 = y01, bin = bin)

  agg_n <- aggregate(y01 ~ bin, data = df, FUN = length)
  agg_pred <- aggregate(prob ~ bin, data = df, FUN = mean)
  agg_obs <- aggregate(y01 ~ bin, data = df, FUN = mean)

  out <- merge(merge(agg_n, agg_pred, by = "bin"), agg_obs, by = "bin")
  colnames(out) <- c("bin", "n", "mean_pred", "obs_rate")
  out <- out[order(out$bin), , drop = FALSE]
  rownames(out) <- NULL
  out
}

plot_calibration_pdf_local <- function(calib_df, file, title_text) {
  p <- ggplot(calib_df, aes(x = mean_pred, y = obs_rate, group = 1)) +
    geom_point(size = 2.5, color = "#2C7FB8") +
    geom_abline(intercept = 0, slope = 1, linetype = 2, color = "grey50") +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(
      title = title_text,
      x = "Mean predicted probability",
      y = "Observed event rate"
    ) +
    theme_bw(base_size = 12)

  if (nrow(calib_df) > 1) {
    p <- p + geom_line(color = "#2C7FB8", linewidth = 0.8)
  }

  ggsave(filename = file, plot = p, width = 6, height = 6)
}

fit_logistic_recalibration_local <- function(prob, y01) {
  prob <- clip_prob_local(prob)
  lp <- qlogis(prob)
  fit <- glm(y01 ~ lp, family = binomial())
  list(
    fit = fit,
    intercept = unname(stats::coef(fit)[1]),
    slope = unname(stats::coef(fit)[2])
  )
}

apply_logistic_recalibration_local <- function(prob, recal_model) {
  prob <- clip_prob_local(prob)
  lp <- qlogis(prob)
  plogis(recal_model$intercept + recal_model$slope * lp)
}

save_metrics_csv_local <- function(file, cohort, model_name, threshold_label, auc, brier, metrics) {
  df <- data.frame(
    Cohort = cohort,
    Model = model_name,
    ThresholdLabel = threshold_label,
    threshold = metrics$threshold,
    AUC = auc,
    Brier = brier,
    Accuracy = metrics$Accuracy,
    Sensitivity = metrics$Sensitivity,
    Specificity = metrics$Specificity,
    PPV = metrics$PPV,
    NPV = metrics$NPV,
    BalancedAccuracy = metrics$BalancedAccuracy,
    F1 = metrics$F1,
    TN = metrics$TN,
    FP = metrics$FP,
    FN = metrics$FN,
    TP = metrics$TP,
    stringsAsFactors = FALSE
  )
  utils::write.csv(df, file = file, row.names = FALSE)
  invisible(df)
}

save_confusion_csv_local <- function(file, truth01, pred01) {
  utils::write.csv(make_confusion_df_local(truth01, pred01), file = file, row.names = TRUE)
}

save_prediction_csv_local <- function(file, sample_ids, truth01, prob, threshold, pred01) {
  df <- data.frame(
    SampleID = sample_ids,
    Truth = ifelse(truth01 == 1, "high-risk", "low-risk"),
    Truth01 = truth01,
    Probability = prob,
    Threshold = threshold,
    Predicted = ifelse(pred01 == 1, "high-risk", "low-risk"),
    Predicted01 = pred01,
    stringsAsFactors = FALSE
  )
  utils::write.csv(df, file = file, row.names = FALSE)
}

# =========================
# Main
# =========================

main <- function() {
  log_info("🚀 启动人类验证模块...")

  config <- read_config()

  out_root <- config$paths$results_ml %||% "results/03_MachineLearning"
  out_dir  <- file.path(out_root, "human_validation")
  ensure_dir(out_dir)

  if (exists("save_session_info", mode = "function")) {
    try(save_session_info(file.path(out_dir, "00_session_info.txt")), silent = TRUE)
  }

  hv_cfg <- config$human_validation %||%
    config$params$human_validation %||%
    list()

  ext_cfg <- config$human_validation_external %||%
    config$params$human_validation_external %||%
    hv_cfg

  model_cfg <- config$model %||%
    config$params$ml %||%
    list()

  log_info("development risk_col = %s", hv_cfg$risk_col %||% "NULL")
  log_info("external risk_col = %s", ext_cfg$risk_col %||% "NULL")

  mouse_gene_file <- resolve_mouse_gene_file(config)
  mouse_genes <- read_candidate_genes(mouse_gene_file)

  if (length(mouse_genes) <= 100) {
    log_info("当前输入基因数 = %d，推测为核心 biomarkers 模式", length(mouse_genes))
  } else {
    log_warn("当前输入基因数 = %d，推测为候选基因池模式（candidate mode），而非最终核心 biomarkers 模式", length(mouse_genes))
  }

  # -------------------------
  # 1) Development cohort
  # -------------------------
  dev_pheno <- read_geo_pheno_manual(
    series_matrix_path = config$paths$human_series_matrix,
    biopsy_col = hv_cfg$biopsy_col %||% "biopsy",
    biopsy_keep_pattern = hv_cfg$biopsy_keep_pattern %||% "1st",
    risk_col = hv_cfg$risk_col %||% "pls-nafld-based risk prediction at 1st biopsy",
    high_label = hv_cfg$high_label %||% "high-risk",
    title_col = hv_cfg$title_col %||% "title"
  )

  dev_expr <- read_gct_matrix(config$paths$human_gct_matrix)

  dev_aligned <- align_pheno_expr(
    pheno = dev_pheno,
    expr = dev_expr,
    title_col = "title",
    cohort_name = "development",
    min_common = hv_cfg$min_common_samples %||% 20
  )

  # -------------------------
  # 2) Mouse -> Human mapping
  # -------------------------
  mapping_res <- build_mouse_to_human_mapping(
    mouse_genes = mouse_genes,
    expr_row_ids = rownames(dev_aligned$expr),
    ortholog_file = config$paths$ortholog_map,
    mouse_deg_file = config$paths$mouse_deg_annotated %||%
      file.path(config$paths$results_deseq2 %||% "results/01_DESeq2", "DESeq2_Annotated_Results.csv")
  )

  utils::write.csv(
    mapping_res$mapping_table,
    file = file.path(out_dir, "mouse_to_human_mapping.csv"),
    row.names = FALSE
  )

  utils::write.csv(
    data.frame(HumanGene = mapping_res$available_human_genes, stringsAsFactors = FALSE),
    file = file.path(out_dir, "available_human_genes.csv"),
    row.names = FALSE
  )

  utils::write.csv(
    data.frame(HumanGene = mapping_res$missing_human_genes, stringsAsFactors = FALSE),
    file = file.path(out_dir, "missing_human_genes_in_expr.csv"),
    row.names = FALSE
  )

  utils::write.csv(
    data.frame(MouseGene = mapping_res$unmapped_mouse_genes, stringsAsFactors = FALSE),
    file = file.path(out_dir, "unmapped_mouse_genes.csv"),
    row.names = FALSE
  )

  feature_names <- unique(mapping_res$available_human_genes)
  if (length(feature_names) == 0) {
    log_stop("❌ 严格直系同源映射未匹配到任何基因！请检查离线同源表或基因 ID 是否对应。")
  }

  # -------------------------
  # 3) Build development matrix
  # -------------------------
  dev_x_obj <- build_feature_matrix(dev_aligned$expr, feature_names)
  x_train <- dev_x_obj$x
  y_train <- factor(
    ifelse(dev_aligned$pheno$RiskLabel == 1, "high-risk", "low-risk"),
    levels = c("low-risk", "high-risk")
  )
  dev_truth01 <- as.integer(y_train == "high-risk")

  log_info(
    "开发队列大小: %d；类别分布: low-risk %d, high-risk %d",
    length(y_train),
    sum(y_train == "low-risk"),
    sum(y_train == "high-risk")
  )

  # -------------------------
  # 4) Fit development model
  # -------------------------
  model_res <- fit_development_model(
    x_train = x_train,
    y_train = y_train,
    alpha = model_cfg$alpha %||% 1,
    seed = model_cfg$seed %||% 123
  )

  cvfit <- get_cvfit_from_result(model_res)
  lambda_selected <- model_res$lambda_selected
  dev_prob <- model_res$dev_prob
  dev_auc <- safe_auc(dev_truth01, dev_prob, "development_original")
  dev_brier <- brier_score_local(dev_prob, dev_truth01)

  log_info("当前实际使用的 lambda_selected = %s", as.character(lambda_selected))
  log_info(
    "development prob summary: min=%.6f, q1=%.6f, median=%.6f, mean=%.6f, q3=%.6f, max=%.6f",
    min(dev_prob), quantile(dev_prob, 0.25), median(dev_prob), mean(dev_prob), quantile(dev_prob, 0.75), max(dev_prob)
  )

  plot_cv_curve_pdf(
    cvfit,
    file = file.path(out_dir, "01_Human_LASSO_CV_Curve.pdf")
  )

  save_scaling_params(
    center = model_res$center,
    scalev = model_res$scalev,
    file = file.path(out_dir, "scaling_params.csv")
  )

  save_roc_curve_pdf(
    y01 = dev_truth01,
    prob = dev_prob,
    file = file.path(out_dir, "02_Development_ROC_Curve.pdf"),
    title_text = "Development ROC Curve"
  )

  coef_df <- extract_coef_df(cvfit, lambda_selected = lambda_selected)
  selected_df <- extract_selected_features(cvfit, lambda_selected = lambda_selected)

  utils::write.csv(
    coef_df,
    file = file.path(out_dir, "human_model_coefficients.csv"),
    row.names = FALSE
  )

  utils::write.csv(
    selected_df,
    file = file.path(out_dir, "selected_human_features.csv"),
    row.names = FALSE
  )

  log_info("开发队列表观 AUC = %.4f（注意：这不是外部验证 AUC）", dev_auc)

  # -------------------------
  # 5) External cohort (optional)
  # -------------------------
  has_external <- !is.null(config$paths$human_series_matrix_external) &&
    !is.null(config$paths$human_gct_matrix_external) &&
    file.exists(config$paths$human_series_matrix_external) &&
    file.exists(config$paths$human_gct_matrix_external)

  ext_truth01 <- NULL
  ext_prob <- NULL
  ext_auc <- NA_real_
  ext_brier <- NA_real_
  x_ext <- NULL
  ext_aligned <- NULL

  if (has_external) {
    log_info("检测到外部验证队列，开始执行 external validation...")

    ext_pheno <- read_geo_pheno_manual(
      series_matrix_path = config$paths$human_series_matrix_external,
      biopsy_col = ext_cfg$biopsy_col,
      biopsy_keep_pattern = ext_cfg$biopsy_keep_pattern,
      risk_col = ext_cfg$risk_col %||% "pls-nafld-based risk prediction",
      high_label = ext_cfg$high_label %||% "high-risk",
      title_col = ext_cfg$title_col %||% "title"
    )

    ext_expr <- read_gct_matrix(config$paths$human_gct_matrix_external)

    ext_aligned <- align_pheno_expr(
      pheno = ext_pheno,
      expr = ext_expr,
      title_col = "title",
      cohort_name = "external",
      min_common = ext_cfg$min_common_samples %||% 20
    )

    ext_x_obj <- build_feature_matrix(ext_aligned$expr, feature_names)
    x_ext <- ext_x_obj$x
    ext_truth01 <- as.integer(ext_aligned$pheno$RiskLabel == 1)

    pred_ext <- predict_external_cohort(
      cvfit = cvfit,
      x_ext = x_ext,
      center = model_res$center,
      scalev = model_res$scalev,
      lambda_selected = lambda_selected
    )

    ext_prob <- pred_ext$ext_prob
    ext_auc <- safe_auc(ext_truth01, ext_prob, "external_original")
    ext_brier <- brier_score_local(ext_prob, ext_truth01)

    log_info(
      "external prob summary: min=%.6f, q1=%.6f, median=%.6f, mean=%.6f, q3=%.6f, max=%.6f",
      min(ext_prob), quantile(ext_prob, 0.25), median(ext_prob), mean(ext_prob), quantile(ext_prob, 0.75), max(ext_prob)
    )

    save_roc_curve_pdf(
      y01 = ext_truth01,
      prob = ext_prob,
      file = file.path(out_dir, "03_External_ROC_Curve.pdf"),
      title_text = "External ROC Curve"
    )

    log_info("外部验证完成：AUC = %.4f", ext_auc)
  } else {
    log_warn("未检测到完整的外部验证输入文件，将跳过 external validation。")
  }

  # -------------------------
  # 6) Brier scores
  # -------------------------
  brier_df <- data.frame(
    Cohort = c("Development", if (has_external) "External"),
    Brier = c(dev_brier, if (has_external) ext_brier),
    stringsAsFactors = FALSE
  )
  utils::write.csv(
    brier_df,
    file = file.path(out_dir, "04_Brier_Scores.csv"),
    row.names = FALSE
  )

  # -------------------------
  # 7) Calibration
  # -------------------------
  dev_calib <- build_calibration_table_local(dev_prob, dev_truth01, n_bins = 8)
  utils::write.csv(
    dev_calib,
    file = file.path(out_dir, "05_Development_Calibration_Table.csv"),
    row.names = FALSE
  )
  plot_calibration_pdf_local(
    dev_calib,
    file = file.path(out_dir, "05_Development_Calibration_Plot.pdf"),
    title_text = "Development Calibration Plot"
  )

  if (has_external) {
    ext_calib <- build_calibration_table_local(ext_prob, ext_truth01, n_bins = 8)
    utils::write.csv(
      ext_calib,
      file = file.path(out_dir, "06_External_Calibration_Table.csv"),
      row.names = FALSE
    )
    plot_calibration_pdf_local(
      ext_calib,
      file = file.path(out_dir, "06_External_Calibration_Plot.pdf"),
      title_text = "External Calibration Plot"
    )
  }

  # -------------------------
  # 8) Threshold optimization on development
  # -------------------------
  youden_df <- find_best_threshold_youden_local(dev_prob, dev_truth01)
  utils::write.csv(
    youden_df,
    file = file.path(out_dir, "07_Threshold_Optimization_Development_Youden.csv"),
    row.names = FALSE
  )

  th_default <- 0.5
  th_youden  <- youden_df$threshold[1]

  dev_m_default <- calc_metrics_at_threshold_local(dev_prob, dev_truth01, th_default)
  dev_m_locked  <- calc_metrics_at_threshold_local(dev_prob, dev_truth01, th_youden)

  save_metrics_csv_local(
    file = file.path(out_dir, "08_Development_DefaultThreshold_0.5_metrics.csv"),
    cohort = "Development",
    model_name = "Original",
    threshold_label = "Default_0.5",
    auc = dev_auc,
    brier = dev_brier,
    metrics = dev_m_default
  )

  save_confusion_csv_local(
    file = file.path(out_dir, "08_Development_DefaultThreshold_0.5_confusion_matrix.csv"),
    truth01 = dev_truth01,
    pred01 = dev_m_default$pred01
  )

  save_metrics_csv_local(
    file = file.path(out_dir, "09_Development_LockedThreshold_Youden_metrics.csv"),
    cohort = "Development",
    model_name = "Original",
    threshold_label = "Locked_Youden",
    auc = dev_auc,
    brier = dev_brier,
    metrics = dev_m_locked
  )

  save_confusion_csv_local(
    file = file.path(out_dir, "09_Development_LockedThreshold_Youden_confusion_matrix.csv"),
    truth01 = dev_truth01,
    pred01 = dev_m_locked$pred01
  )

  save_prediction_csv_local(
    file = file.path(out_dir, "12_Development_Predictions_With_Thresholds.csv"),
    sample_ids = dev_aligned$pheno$MatchedSampleID,
    truth01 = dev_truth01,
    prob = dev_prob,
    threshold = th_youden,
    pred01 = dev_m_locked$pred01
  )

  if (has_external) {
    ext_m_default <- calc_metrics_at_threshold_local(ext_prob, ext_truth01, th_default)
    ext_m_locked  <- calc_metrics_at_threshold_local(ext_prob, ext_truth01, th_youden)

    save_metrics_csv_local(
      file = file.path(out_dir, "10_External_DefaultThreshold_0.5_metrics.csv"),
      cohort = "External",
      model_name = "Original",
      threshold_label = "Default_0.5",
      auc = ext_auc,
      brier = ext_brier,
      metrics = ext_m_default
    )

    save_confusion_csv_local(
      file = file.path(out_dir, "10_External_DefaultThreshold_0.5_confusion_matrix.csv"),
      truth01 = ext_truth01,
      pred01 = ext_m_default$pred01
    )

    save_metrics_csv_local(
      file = file.path(out_dir, "11_External_LockedThreshold_YoudenFromDevelopment_metrics.csv"),
      cohort = "External",
      model_name = "Original",
      threshold_label = "Locked_YoudenFromDevelopment",
      auc = ext_auc,
      brier = ext_brier,
      metrics = ext_m_locked
    )

    save_confusion_csv_local(
      file = file.path(out_dir, "11_External_LockedThreshold_YoudenFromDevelopment_confusion_matrix.csv"),
      truth01 = ext_truth01,
      pred01 = ext_m_locked$pred01
    )

    save_prediction_csv_local(
      file = file.path(out_dir, "13_External_Predictions_With_Thresholds.csv"),
      sample_ids = ext_aligned$pheno$MatchedSampleID,
      truth01 = ext_truth01,
      prob = ext_prob,
      threshold = th_youden,
      pred01 = ext_m_locked$pred01
    )
  }

  summary_df <- data.frame(
    Cohort = c("Development", if (has_external) "External"),
    AUC = c(dev_auc, if (has_external) ext_auc),
    Brier = c(dev_brier, if (has_external) ext_brier),
    LockedThreshold_Youden = c(th_youden, if (has_external) th_youden),
    stringsAsFactors = FALSE
  )
  utils::write.csv(
    summary_df,
    file = file.path(out_dir, "14_Calibration_Threshold_Summary.csv"),
    row.names = FALSE
  )

  # -------------------------
  # 9) Logistic recalibration (exploratory)
  # -------------------------
  tryCatch({
    if (length(unique(round(dev_prob, 12))) < 2) {
      stop("development 预测概率近似常数，跳过 logistic recalibration")
    }

    recal_model <- fit_logistic_recalibration_local(dev_prob, dev_truth01)

    recal_model_df <- data.frame(
      intercept = recal_model$intercept,
      slope = recal_model$slope,
      stringsAsFactors = FALSE
    )
    utils::write.csv(
      recal_model_df,
      file = file.path(out_dir, "15_Logistic_Recalibration_Model.csv"),
      row.names = FALSE
    )

    dev_prob_recal <- apply_logistic_recalibration_local(dev_prob, recal_model)
    dev_auc_recal <- safe_auc(dev_truth01, dev_prob_recal, "development_recal")
    dev_brier_recal <- brier_score_local(dev_prob_recal, dev_truth01)

    dev_calib_recal <- build_calibration_table_local(dev_prob_recal, dev_truth01, n_bins = 8)
    utils::write.csv(
      dev_calib_recal,
      file = file.path(out_dir, "16_Development_Calibration_BeforeAfter_Table.csv"),
      row.names = FALSE
    )
    plot_calibration_pdf_local(
      dev_calib_recal,
      file = file.path(out_dir, "16_Development_Calibration_BeforeAfter_Plot.pdf"),
      title_text = "Development Calibration After Logistic Recalibration"
    )

    th_youden_recal <- find_best_threshold_youden_local(dev_prob_recal, dev_truth01)$threshold[1]

    threshold_cmp_df <- data.frame(
      Model = c("Original", "Recalibrated"),
      threshold = c(th_youden, th_youden_recal),
      stringsAsFactors = FALSE
    )
    utils::write.csv(
      threshold_cmp_df,
      file = file.path(out_dir, "18_Recalibration_Threshold_Comparison.csv"),
      row.names = FALSE
    )

    dev_m_orig_default <- calc_metrics_at_threshold_local(dev_prob, dev_truth01, 0.5)
    dev_m_orig_locked  <- calc_metrics_at_threshold_local(dev_prob, dev_truth01, th_youden)
    dev_m_recal_default <- calc_metrics_at_threshold_local(dev_prob_recal, dev_truth01, 0.5)
    dev_m_recal_locked  <- calc_metrics_at_threshold_local(dev_prob_recal, dev_truth01, th_youden_recal)

    dev_metrics_compare <- rbind(
      data.frame(
        Cohort = "Development", Model = "Original", ThresholdLabel = "Default_0.5",
        threshold = dev_m_orig_default$threshold, AUC = dev_auc, Brier = dev_brier,
        Accuracy = dev_m_orig_default$Accuracy, Sensitivity = dev_m_orig_default$Sensitivity,
        Specificity = dev_m_orig_default$Specificity, PPV = dev_m_orig_default$PPV,
        NPV = dev_m_orig_default$NPV, BalancedAccuracy = dev_m_orig_default$BalancedAccuracy,
        F1 = dev_m_orig_default$F1, TN = dev_m_orig_default$TN, FP = dev_m_orig_default$FP,
        FN = dev_m_orig_default$FN, TP = dev_m_orig_default$TP,
        stringsAsFactors = FALSE
      ),
      data.frame(
        Cohort = "Development", Model = "Original", ThresholdLabel = "Locked_Youden",
        threshold = dev_m_orig_locked$threshold, AUC = dev_auc, Brier = dev_brier,
        Accuracy = dev_m_orig_locked$Accuracy, Sensitivity = dev_m_orig_locked$Sensitivity,
        Specificity = dev_m_orig_locked$Specificity, PPV = dev_m_orig_locked$PPV,
        NPV = dev_m_orig_locked$NPV, BalancedAccuracy = dev_m_orig_locked$BalancedAccuracy,
        F1 = dev_m_orig_locked$F1, TN = dev_m_orig_locked$TN, FP = dev_m_orig_locked$FP,
        FN = dev_m_orig_locked$FN, TP = dev_m_orig_locked$TP,
        stringsAsFactors = FALSE
      ),
      data.frame(
        Cohort = "Development", Model = "Recalibrated", ThresholdLabel = "Default_0.5",
        threshold = dev_m_recal_default$threshold, AUC = dev_auc_recal, Brier = dev_brier_recal,
        Accuracy = dev_m_recal_default$Accuracy, Sensitivity = dev_m_recal_default$Sensitivity,
        Specificity = dev_m_recal_default$Specificity, PPV = dev_m_recal_default$PPV,
        NPV = dev_m_recal_default$NPV, BalancedAccuracy = dev_m_recal_default$BalancedAccuracy,
        F1 = dev_m_recal_default$F1, TN = dev_m_recal_default$TN, FP = dev_m_recal_default$FP,
        FN = dev_m_recal_default$FN, TP = dev_m_recal_default$TP,
        stringsAsFactors = FALSE
      ),
      data.frame(
        Cohort = "Development", Model = "Recalibrated", ThresholdLabel = "Locked_Youden",
        threshold = dev_m_recal_locked$threshold, AUC = dev_auc_recal, Brier = dev_brier_recal,
        Accuracy = dev_m_recal_locked$Accuracy, Sensitivity = dev_m_recal_locked$Sensitivity,
        Specificity = dev_m_recal_locked$Specificity, PPV = dev_m_recal_locked$PPV,
        NPV = dev_m_recal_locked$NPV, BalancedAccuracy = dev_m_recal_locked$BalancedAccuracy,
        F1 = dev_m_recal_locked$F1, TN = dev_m_recal_locked$TN, FP = dev_m_recal_locked$FP,
        FN = dev_m_recal_locked$FN, TP = dev_m_recal_locked$TP,
        stringsAsFactors = FALSE
      )
    )

    utils::write.csv(
      dev_metrics_compare,
      file = file.path(out_dir, "19_Development_Metrics_BeforeAfter_Recalibration.csv"),
      row.names = FALSE
    )

    save_prediction_csv_local(
      file = file.path(out_dir, "21_Development_Predictions_Original_vs_Recalibrated.csv"),
      sample_ids = dev_aligned$pheno$MatchedSampleID,
      truth01 = dev_truth01,
      prob = dev_prob_recal,
      threshold = th_youden_recal,
      pred01 = dev_m_recal_locked$pred01
    )

    if (has_external) {
      ext_prob_recal <- apply_logistic_recalibration_local(ext_prob, recal_model)
      ext_auc_recal <- safe_auc(ext_truth01, ext_prob_recal, "external_recal")
      ext_brier_recal <- brier_score_local(ext_prob_recal, ext_truth01)

      ext_calib_recal <- build_calibration_table_local(ext_prob_recal, ext_truth01, n_bins = 8)
      utils::write.csv(
        ext_calib_recal,
        file = file.path(out_dir, "17_External_Calibration_BeforeAfter_Table.csv"),
        row.names = FALSE
      )
      plot_calibration_pdf_local(
        ext_calib_recal,
        file = file.path(out_dir, "17_External_Calibration_BeforeAfter_Plot.pdf"),
        title_text = "External Calibration After Logistic Recalibration"
      )

      ext_m_orig_default <- calc_metrics_at_threshold_local(ext_prob, ext_truth01, 0.5)
      ext_m_orig_locked  <- calc_metrics_at_threshold_local(ext_prob, ext_truth01, th_youden)
      ext_m_recal_default <- calc_metrics_at_threshold_local(ext_prob_recal, ext_truth01, 0.5)
      ext_m_recal_locked  <- calc_metrics_at_threshold_local(ext_prob_recal, ext_truth01, th_youden_recal)

      ext_metrics_compare <- rbind(
        data.frame(
          Cohort = "External", Model = "Original", ThresholdLabel = "Default_0.5",
          threshold = ext_m_orig_default$threshold, AUC = ext_auc, Brier = ext_brier,
          Accuracy = ext_m_orig_default$Accuracy, Sensitivity = ext_m_orig_default$Sensitivity,
          Specificity = ext_m_orig_default$Specificity, PPV = ext_m_orig_default$PPV,
          NPV = ext_m_orig_default$NPV, BalancedAccuracy = ext_m_orig_default$BalancedAccuracy,
          F1 = ext_m_orig_default$F1, TN = ext_m_orig_default$TN, FP = ext_m_orig_default$FP,
          FN = ext_m_orig_default$FN, TP = ext_m_orig_default$TP,
          stringsAsFactors = FALSE
        ),
        data.frame(
          Cohort = "External", Model = "Original", ThresholdLabel = "Locked_Youden",
          threshold = ext_m_orig_locked$threshold, AUC = ext_auc, Brier = ext_brier,
          Accuracy = ext_m_orig_locked$Accuracy, Sensitivity = ext_m_orig_locked$Sensitivity,
          Specificity = ext_m_orig_locked$Specificity, PPV = ext_m_orig_locked$PPV,
          NPV = ext_m_orig_locked$NPV, BalancedAccuracy = ext_m_orig_locked$BalancedAccuracy,
          F1 = ext_m_orig_locked$F1, TN = ext_m_orig_locked$TN, FP = ext_m_orig_locked$FP,
          FN = ext_m_orig_locked$FN, TP = ext_m_orig_locked$TP,
          stringsAsFactors = FALSE
        ),
        data.frame(
          Cohort = "External", Model = "Recalibrated", ThresholdLabel = "Default_0.5",
          threshold = ext_m_recal_default$threshold, AUC = ext_auc_recal, Brier = ext_brier_recal,
          Accuracy = ext_m_recal_default$Accuracy, Sensitivity = ext_m_recal_default$Sensitivity,
          Specificity = ext_m_recal_default$Specificity, PPV = ext_m_recal_default$PPV,
          NPV = ext_m_recal_default$NPV, BalancedAccuracy = ext_m_recal_default$BalancedAccuracy,
          F1 = ext_m_recal_default$F1, TN = ext_m_recal_default$TN, FP = ext_m_recal_default$FP,
          FN = ext_m_recal_default$FN, TP = ext_m_recal_default$TP,
          stringsAsFactors = FALSE
        ),
        data.frame(
          Cohort = "External", Model = "Recalibrated", ThresholdLabel = "Locked_Youden",
          threshold = ext_m_recal_locked$threshold, AUC = ext_auc_recal, Brier = ext_brier_recal,
          Accuracy = ext_m_recal_locked$Accuracy, Sensitivity = ext_m_recal_locked$Sensitivity,
          Specificity = ext_m_recal_locked$Specificity, PPV = ext_m_recal_locked$PPV,
          NPV = ext_m_recal_locked$NPV, BalancedAccuracy = ext_m_recal_locked$BalancedAccuracy,
          F1 = ext_m_recal_locked$F1, TN = ext_m_recal_locked$TN, FP = ext_m_recal_locked$FP,
          FN = ext_m_recal_locked$FN, TP = ext_m_recal_locked$TP,
          stringsAsFactors = FALSE
        )
      )

      utils::write.csv(
        ext_metrics_compare,
        file = file.path(out_dir, "20_External_Metrics_BeforeAfter_Recalibration.csv"),
        row.names = FALSE
      )

      save_prediction_csv_local(
        file = file.path(out_dir, "22_External_Predictions_Original_vs_Recalibrated.csv"),
        sample_ids = ext_aligned$pheno$MatchedSampleID,
        truth01 = ext_truth01,
        prob = ext_prob_recal,
        threshold = th_youden_recal,
        pred01 = ext_m_recal_locked$pred01
      )

      recal_summary <- data.frame(
        Cohort = c("Development", "External"),
        AUC = c(dev_auc_recal, ext_auc_recal),
        Brier = c(dev_brier_recal, ext_brier_recal),
        YoudenThreshold = c(th_youden_recal, th_youden_recal),
        stringsAsFactors = FALSE
      )
    } else {
      recal_summary <- data.frame(
        Cohort = "Development",
        AUC = dev_auc_recal,
        Brier = dev_brier_recal,
        YoudenThreshold = th_youden_recal,
        stringsAsFactors = FALSE
      )
    }

    utils::write.csv(
      recal_summary,
      file = file.path(out_dir, "23_Logistic_Recalibration_Summary.csv"),
      row.names = FALSE
    )

    log_info("✅ Logistic recalibration（探索性）已完成。结果已保存到: %s", out_dir)

  }, error = function(e) {
    log_warn("logistic recalibration 阶段报错，已跳过：%s", conditionMessage(e))
  })

  model_metrics <- data.frame(
    Cohort = c("Development", if (has_external) "External"),
    N = c(length(dev_truth01), if (has_external) length(ext_truth01)),
    AUC = c(dev_auc, if (has_external) ext_auc),
    Brier = c(dev_brier, if (has_external) ext_brier),
    CandidateMouseGenes = length(mouse_genes),
    HumanMappedGenes = length(feature_names),
    SelectedHumanFeatures = nrow(selected_df),
    stringsAsFactors = FALSE
  )
  utils::write.csv(
    model_metrics,
    file = file.path(out_dir, "model_metrics.csv"),
    row.names = FALSE
  )

  log_info("✅ Calibration + threshold optimization 已完成。结果已保存到: %s", out_dir)
  log_info("✅ 人类验证模块完成。输出目录: %s", out_dir)
}

main()
