################################################
##### 04_4_hv_modeling.R (核心机器学习模块) ######
################################################

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

scale_train_matrix <- function(x) {
  center <- colMeans(x, na.rm = TRUE)
  scalev <- apply(x, 2, stats::sd, na.rm = TRUE)
  scalev[is.na(scalev) | scalev == 0] <- 1

  x_scaled <- scale(x, center = center, scale = scalev)
  x_scaled <- safe_numeric_matrix(x_scaled)

  list(
    x_scaled = x_scaled,
    center = center,
    scalev = scalev
  )
}

scale_new_matrix <- function(x, center, scalev) {
  all_cols <- colnames(x)
  center <- center[all_cols]
  scalev <- scalev[all_cols]
  scalev[is.na(scalev) | scalev == 0] <- 1
  x_scaled <- scale(x, center = center, scale = scalev)
  safe_numeric_matrix(x_scaled)
}

fit_development_model <- function(x_train, y_train, alpha = 1, seed = 123) {
  check_binary_outcome(y_train, min_count = 4, label = "development_y")
  train_tab <- table(y_train)

  scaled <- scale_train_matrix(x_train)
  x_train_scaled <- scaled$x_scaled
  center <- scaled$center
  scalev <- scaled$scalev

  k_target <- 10
  k <- min(k_target, min(train_tab))

  if (min(train_tab) < 2) {
    log_stop("开发队列某一类别样本数过少，无法进行交叉验证。")
  }

  foldid <- make_stratified_foldid(y_train, k = k, seed = seed)

  log_info(
    "正在训练开发队列模型：cv.glmnet (nfolds=%d, low-risk=%d, high-risk=%d)",
    k,
    unname(train_tab["low-risk"]),
    unname(train_tab["high-risk"])
  )

  y01 <- as.integer(y_train == "high-risk")

  set.seed(seed)
  cvfit <- glmnet::cv.glmnet(
    x = x_train_scaled,
    y = y01,
    family = "binomial",
    alpha = alpha,
    standardize = FALSE,
    nfolds = k,
    foldid = foldid,
    type.measure = "auc"
  )

  lambda_selected <- cvfit$lambda.1se
  coef_mat <- as.matrix(stats::coef(cvfit, s = lambda_selected))
  coef_df <- data.frame(
    Feature = rownames(coef_mat),
    Coefficient = as.numeric(coef_mat[, 1]),
    stringsAsFactors = FALSE
  )

  selected_features <- coef_df$Feature[coef_df$Feature != "(Intercept)" & coef_df$Coefficient != 0]

  if (length(selected_features) == 0) {
    log_warn("lambda.1se 下没有非零特征，回退到 lambda.min")
    lambda_selected <- cvfit$lambda.min
    coef_mat <- as.matrix(stats::coef(cvfit, s = lambda_selected))
    coef_df <- data.frame(
      Feature = rownames(coef_mat),
      Coefficient = as.numeric(coef_mat[, 1]),
      stringsAsFactors = FALSE
    )
    selected_features <- coef_df$Feature[coef_df$Feature != "(Intercept)" & coef_df$Coefficient != 0]
  }

  dev_prob <- as.numeric(stats::predict(cvfit, newx = x_train_scaled, s = lambda_selected, type = "response"))
  dev_auc <- as.numeric(pROC::auc(pROC::roc(y01, dev_prob, quiet = TRUE, direction = "<")))

  log_info("开发队列表观 AUC = %.4f（注意：这不是外部验证 AUC）", dev_auc)

  list(
    cvfit = cvfit,
    lambda_selected = lambda_selected,
    coef_df = coef_df,
    selected_features = selected_features,
    dev_prob = dev_prob,
    dev_auc = dev_auc,
    x_train_scaled = x_train_scaled,
    center = center,
    scalev = scalev,
    k = k
  )
}

predict_external_cohort <- function(cvfit, x_ext, center, scalev, lambda_selected) {
  x_ext_scaled <- scale_new_matrix(x_ext, center, scalev)
  ext_prob <- as.numeric(stats::predict(cvfit, newx = x_ext_scaled, s = lambda_selected, type = "response"))
  list(
    x_ext_scaled = x_ext_scaled,
    ext_prob = ext_prob
  )
}

plot_cv_curve_pdf <- function(cvfit, file) {
  grDevices::pdf(file, width = 7, height = 6)
  plot(cvfit)
  grDevices::dev.off()
}

plot_roc_pdf <- function(truth01, prob, title, file, color = "#1f77b4") {
  roc_obj <- pROC::roc(truth01, prob, quiet = TRUE, direction = "<")
  auc_val <- as.numeric(pROC::auc(roc_obj))

  grDevices::pdf(file, width = 6.5, height = 6)
  plot(
    roc_obj,
    col = color,
    lwd = 3,
    main = sprintf("%s\nAUC = %.4f", title, auc_val),
    legacy.axes = TRUE
  )
  abline(a = 0, b = 1, lty = 2, col = "gray60")
  grDevices::dev.off()
}

save_scaling_params <- function(center, scalev, file) {
  df <- data.frame(
    Feature = names(center),
    Center = as.numeric(center),
    Scale = as.numeric(scalev),
    stringsAsFactors = FALSE
  )
  utils::write.csv(df, file = file, row.names = FALSE)
}
