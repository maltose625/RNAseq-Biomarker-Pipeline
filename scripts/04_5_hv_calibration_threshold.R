##################################################################
###### 04_5_hv_calibration_threshold.R (模型评估与阈值模块)  ######
##################################################################

as_binary01 <- function(y) {
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

  segments(
    x0 = calib_df$mean_pred, y0 = calib_df$obs_lower,
    x1 = calib_df$mean_pred, y1 = calib_df$obs_upper,
    col = "#1f77b4", lwd = 1.5
  )

  ord <- order(calib_df$mean_pred)
  lines(
    calib_df$mean_pred[ord],
    calib_df$obs_rate[ord],
    col = "#1f77b4", lwd = 2
  )

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

run_calibration_threshold <- function(dev_truth01, dev_prob, ext_truth01 = NULL, ext_prob = NULL, outdir) {
  ensure_dir(outdir)

  dev_truth01 <- as_binary01(dev_truth01)
  dev_prob <- as.numeric(dev_prob)

  has_external <- !(is.null(ext_truth01) || is.null(ext_prob))
  if (has_external) {
    ext_truth01 <- as_binary01(ext_truth01)
    ext_prob <- as.numeric(ext_prob)
  }

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

  best_thr_df <- find_best_threshold_youden(dev_truth01, dev_prob)
  best_threshold <- as.numeric(best_thr_df$threshold[1])

  log_info("Development cohort 最佳阈值（Youden） = %.4f", best_threshold)

  utils::write.csv(
    best_thr_df,
    file = file.path(outdir, "07_Threshold_Optimization_Development_Youden.csv"),
    row.names = FALSE
  )

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

  summary_rows <- list(
    data.frame(
      Cohort = "development",
      AUC = tryCatch(as.numeric(pROC::auc(pROC::roc(dev_truth01, dev_prob, quiet = TRUE, direction = "<"))), error = function(e) NA_real_),
      Brier = dev_brier,
      DefaultThreshold = 0.5,
      LockedThreshold = best_threshold,
      stringsAsFactors = FALSE
    )
  )

  if (has_external) {
    summary_rows[[2]] <- data.frame(
      Cohort = "external",
      AUC = tryCatch(as.numeric(pROC::auc(pROC::roc(ext_truth01, ext_prob, quiet = TRUE, direction = "<"))), error = function(e) NA_real_),
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

  list(
    best_threshold = best_threshold,
    summary = summary_df
  )
}
