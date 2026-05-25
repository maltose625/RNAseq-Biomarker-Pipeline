#####################################################
##### 04_6_hv_recalibration.R (后处理重校准模块) ######
#####################################################

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

run_logistic_recalibration <- function(dev_truth01, dev_prob, ext_truth01 = NULL, ext_prob = NULL, outdir) {
  ensure_dir(outdir)

  dev_truth01 <- as_binary01(dev_truth01)
  dev_prob <- as.numeric(dev_prob)

  has_external <- !(is.null(ext_truth01) || is.null(ext_prob))
  if (has_external) {
    ext_truth01 <- as_binary01(ext_truth01)
    ext_prob <- as.numeric(ext_prob)
  }

  recal_fit <- fit_logistic_recalibration(dev_truth01, dev_prob)
  recal_intercept <- recal_fit$intercept
  recal_slope <- recal_fit$slope

  log_info("Logistic recalibration 拟合完成：intercept = %.6f, slope = %.6f",
           recal_intercept, recal_slope)

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

  dev_prob_recal <- apply_logistic_recalibration(dev_prob, recal_intercept, recal_slope)
  if (has_external) {
    ext_prob_recal <- apply_logistic_recalibration(ext_prob, recal_intercept, recal_slope)
  }

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

  if (has_external) {
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

  if (has_external) {
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

  if (has_external) {
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

  if (has_external) {
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

  log_info("✅ Logistic recalibration 已完成。输出文件：15~23 前缀结果已保存到: %s", outdir)
  log_info("重校准参数：intercept = %.6f, slope = %.6f", recal_intercept, recal_slope)

  list(
    intercept = recal_intercept,
    slope = recal_slope,
    summary = summary_recal_df
  )
}
