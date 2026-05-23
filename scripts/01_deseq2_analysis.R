#!/usr/bin/env Rscript
# =========================================================================
# 脚本名称: 01_deseq2_analysis.R
# 模块功能: DESeq2 差异分析 + 注释 + PCA/Volcano 可视化
# =========================================================================

suppressPackageStartupMessages({
  library(yaml)
  library(DESeq2)
  library(ggplot2)
  library(org.Mm.eg.db)
  library(ggrepel)
})

config <- read_yaml("config/config.yaml")
if (!dir.exists(config$paths$log_dir)) dir.create(config$paths$log_dir, recursive = TRUE)
log_file_path <- file.path(config$paths$log_dir, "01_deseq2.log")

log_con <- file(log_file_path, open = "wt")
sink(log_con)
sink(log_con, type = "message")

cat(sprintf("[%s] [INFO] 🚀 启动 DESeq2 自动化分析模块...\n", Sys.time()))

tryCatch({

  counts_file <- config$paths$counts_matrix
  if (is.null(counts_file) || !file.exists(counts_file)) {
  stop(sprintf("❌ 缺失 counts_matrix 文件：%s", counts_file))
  }
  out_dir     <- config$paths$results_deseq2
  ref_group   <- config$params$deseq2$ref_group
  trt_group   <- config$params$deseq2$trt_group
  p_cutoff <- NULL
  if (!is.null(config$thresholds$padj_cutoff)) {
  p_cutoff <- as.numeric(config$thresholds$padj_cutoff)
  }
  else if (!is.null(config$params$deseq2$pvalue_cutoff)) {
  p_cutoff <- as.numeric(config$params$deseq2$pvalue_cutoff)
  }
  
  fc_cutoff <- NULL
  if (!is.null(config$thresholds$log2fc_cutoff)) {
  fc_cutoff <- as.numeric(config$thresholds$log2fc_cutoff)
  }
  else if (!is.null(config$params$deseq2$log2fc_cutoff)) {
  fc_cutoff <- as.numeric(config$params$deseq2$log2fc_cutoff)
  }

  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  meta_file <- config$paths$sample_info
  if (is.null(meta_file) || !file.exists(meta_file)) {
  if (file.exists("config/sample_info.tsv")) {
    meta_file <- "config/sample_info.tsv"
  } else {
    stop(sprintf("❌ 缺失 metadata 文件：%s", meta_file))
  }}

  meta_data <- read.table(meta_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  raw_counts <- read.table(counts_file, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)

  bam_cols <- grep("\\.bam$", colnames(raw_counts), value = TRUE)
  if (length(bam_cols) == 0) stop("❌ counts matrix 中未找到 .bam 结尾的样本列")

  count_matrix <- raw_counts[, bam_cols, drop = FALSE]
  colnames(count_matrix) <- gsub(".*\\/|\\_sorted\\.bam", "", colnames(count_matrix))

  common_samples <- intersect(colnames(count_matrix), meta_data$SampleID)
  if (length(common_samples) < 3) stop("❌ 共同样本数过少，无法分析")

  count_matrix <- count_matrix[, common_samples, drop = FALSE]
  rownames(meta_data) <- meta_data$SampleID
  group_info <- meta_data[common_samples, , drop = FALSE]
  group_info$Condition <- factor(group_info$Group, levels = c(ref_group, trt_group))

  cat(sprintf("[%s] [INFO] 正在拟合 DESeq2 模型...\n", Sys.time()))
  dds <- DESeqDataSetFromMatrix(
    countData = round(as.matrix(count_matrix)),
    colData   = group_info,
    design    = ~ Condition
  )

  keep <- rowSums(counts(dds) >= 10) >= 3
  dds <- dds[keep, ]
  dds <- DESeq(dds)

  res <- results(dds, contrast = c("Condition", trt_group, ref_group), alpha = p_cutoff)
  res_df <- as.data.frame(res)

  # 注释
  res_df$Ensembl <- gsub("\\..*", "", rownames(res_df))
  res_df$Symbol <- mapIds(
    org.Mm.eg.db,
    keys      = res_df$Ensembl,
    column    = "SYMBOL",
    keytype   = "ENSEMBL",
    multiVals = "first"
  )
  res_df$Symbol[is.na(res_df$Symbol)] <- rownames(res_df)[is.na(res_df$Symbol)]

  # 显著性分类
  res_df$Significance <- "Not Significant"
  res_df$Significance[!is.na(res_df$padj) & res_df$padj < p_cutoff & res_df$log2FoldChange >  fc_cutoff] <- "Up-regulated"
  res_df$Significance[!is.na(res_df$padj) & res_df$padj < p_cutoff & res_df$log2FoldChange < -fc_cutoff] <- "Down-regulated"

  # 同时输出两个文件，避免下游读错
  write.csv(res_df, file = file.path(out_dir, "DESeq2_Annotated_Results.csv"), row.names = TRUE)
  write.csv(res_df, file = file.path(out_dir, "DESeq2_All_Results.csv"),       row.names = TRUE)

  cat(sprintf("[%s] [INFO] 已输出差异分析结果表。\n", Sys.time()))

  # PCA
  cat(sprintf("[%s] [INFO] 正在生成 PCA 图...\n", Sys.time()))
  vsd <- vst(dds, blind = FALSE)
  pca_data <- plotPCA(vsd, intgroup = c("Condition"), returnData = TRUE)
  percentVar <- round(100 * attr(pca_data, "percentVar"))
  pca_data$SampleID <- rownames(pca_data)

  pca_plot <- ggplot(pca_data, aes(PC1, PC2, color = Condition, label = SampleID)) +
    geom_point(size = 4, alpha = 0.85) +
    geom_text_repel(size = 3.5, max.overlaps = Inf) +
    theme_bw(base_size = 14) +
    labs(
      title = "PCA of RNA-seq Samples",
      x = paste0("PC1: ", percentVar[1], "% variance"),
      y = paste0("PC2: ", percentVar[2], "% variance")
    ) +
    scale_color_manual(values = c("#377eb8", "#e41a1c"))

  ggsave(file.path(out_dir, "01_PCA_Plot.pdf"), plot = pca_plot, width = 7, height = 5.5)

  # Volcano：上调前5 + 下调前5
  cat(sprintf("[%s] [INFO] 正在生成火山图...\n", Sys.time()))

  sig_df <- subset(
    res_df,
    !is.na(padj) & is.finite(log2FoldChange) &
      Significance %in% c("Up-regulated", "Down-regulated")
  )

  label_up <- head(sig_df[sig_df$Significance == "Up-regulated", ][order(sig_df[sig_df$Significance == "Up-regulated", ]$padj), ], 5)
  label_dn <- head(sig_df[sig_df$Significance == "Down-regulated", ][order(sig_df[sig_df$Significance == "Down-regulated", ]$padj), ], 5)

  label_ids <- unique(c(rownames(label_up), rownames(label_dn)))
  res_df$Label <- ifelse(rownames(res_df) %in% label_ids, res_df$Symbol, "")

  volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
    geom_point(alpha = 0.65, size = 1.5, na.rm = TRUE) +
    scale_color_manual(values = c(
      "Up-regulated"   = "#e41a1c",
      "Down-regulated" = "#377eb8",
      "Not Significant" = "grey80"
    )) +
    geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed", color = "black") +
    geom_text_repel(aes(label = Label), size = 3.8, box.padding = 0.5, color = "black", max.overlaps = 20, na.rm = TRUE) +
    theme_minimal(base_size = 14) +
    labs(
      title = "Volcano Plot of DEGs",
      x = expression(Log[2] * " Fold Change"),
      y = expression(-Log[10] * " Adjusted P-value")
    )

  ggsave(file.path(out_dir, "02_Volcano_Plot.pdf"), plot = volcano_plot, width = 8, height = 6)

  cat(sprintf("[%s] [INFO] ✅ DESeq2 模块运行完成。\n", Sys.time()))

}, error = function(e) {
  cat(sprintf("[%s] [ERROR] ❌ 程序崩溃: %s\n", Sys.time(), e$message))
  quit(status = 1)
}, finally = {
  cat(sprintf("[%s] [INFO] 模块运行结束。\n", Sys.time()))
  sink(type = "message")
  sink()
  close(log_con)
})
