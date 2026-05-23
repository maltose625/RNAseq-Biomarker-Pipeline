#!/usr/bin/env Rscript
# ==============================================================================
# 脚本名称: 02_enrichment.R
# 模块功能: GO / KEGG 富集分析（统一输入、统一阈值、只展示前10条）
# ==============================================================================

suppressPackageStartupMessages({
  library(yaml)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(enrichplot)
  library(AnnotationDbi)
  library(ggplot2)
  library(dplyr)
})

message("⏳ [1/5] 读取 config.yaml ...")

config_path <- "config/config.yaml"
if (!file.exists(config_path)) {
  stop("❌ 未找到配置文件：config/config.yaml")
}
config <- yaml::read_yaml(config_path)

# ---------- 读取阈值 ----------
p_cutoff <- NULL
fc_cutoff <- NULL

if (!is.null(config$thresholds$padj_cutoff)) {
  p_cutoff <- as.numeric(config$thresholds$padj_cutoff)
} else if (!is.null(config$params$enrichment$padj_cutoff)) {
  p_cutoff <- as.numeric(config$params$enrichment$padj_cutoff)
} else if (!is.null(config$params$deseq2$pvalue_cutoff)) {
  p_cutoff <- as.numeric(config$params$deseq2$pvalue_cutoff)
}

if (!is.null(config$thresholds$log2fc_cutoff)) {
  fc_cutoff <- as.numeric(config$thresholds$log2fc_cutoff)
} else if (!is.null(config$params$enrichment$log2fc_cutoff)) {
  fc_cutoff <- as.numeric(config$params$enrichment$log2fc_cutoff)
} else if (!is.null(config$params$deseq2$log2fc_cutoff)) {
  fc_cutoff <- as.numeric(config$params$deseq2$log2fc_cutoff)
}

if (length(p_cutoff) != 1 || is.na(p_cutoff)) {
  stop("❌ padj_cutoff 读取失败，请检查 config.yaml")
}
if (length(fc_cutoff) != 1 || is.na(fc_cutoff)) {
  stop("❌ log2fc_cutoff 读取失败，请检查 config.yaml")
}

message(sprintf("   [配置激活] padj < %.3f 且 |log2FC| > %.3f", p_cutoff, fc_cutoff))
# ---------- 输出目录 ----------
out_dir <- NULL
if (!is.null(config$paths$results_enrichment)) {
  out_dir <- as.character(config$paths$results_enrichment)
}
if (length(out_dir) != 1 || is.na(out_dir) || !nzchar(out_dir)) {
  out_dir <- "results/02_Enrichment"
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
message(sprintf("   [输出目录] %s", normalizePath(out_dir)))

# ---------- 读取 DESeq2 结果 ----------
message("⏳ [2/5] 加载 DESeq2 结果 ...")

res_file <- config$paths$deseq2_result_file
if (is.null(res_file) || !file.exists(res_file)) {
  stop(sprintf("❌ 未找到结果文件：%s", res_file))
}

deg_data <- read.csv(res_file, check.names = FALSE, stringsAsFactors = FALSE)

message(sprintf("   [文件路径] %s", normalizePath(res_file)))
message(sprintf("   [数据维度] %d 行 × %d 列", nrow(deg_data), ncol(deg_data)))
message(sprintf("   [列名] %s", paste(colnames(deg_data), collapse = ", ")))

if (!"padj" %in% colnames(deg_data)) stop("❌ DESeq2 结果缺少 padj 列")
if (!"log2FoldChange" %in% colnames(deg_data)) stop("❌ DESeq2 结果缺少 log2FoldChange 列")

deg_data$padj <- suppressWarnings(as.numeric(deg_data$padj))
deg_data$log2FoldChange <- suppressWarnings(as.numeric(deg_data$log2FoldChange))

manual_deg_n <- sum(
  !is.na(deg_data$padj) &
  !is.na(deg_data$log2FoldChange) &
  deg_data$padj < p_cutoff &
  abs(deg_data$log2FoldChange) > fc_cutoff
)

message(sprintf("   [手工复算] 满足 padj < %.3f 且 |log2FC| > %.3f 的基因数: %d",
                p_cutoff, fc_cutoff, manual_deg_n))

# ---------- 显著基因筛选 ----------
signif_degs <- subset(
  deg_data,
  !is.na(padj) &
    !is.na(log2FoldChange) &
    padj < p_cutoff &
    abs(log2FoldChange) > fc_cutoff
)

up_n <- sum(!is.na(deg_data$padj) & !is.na(deg_data$log2FoldChange) &
              deg_data$padj < p_cutoff & deg_data$log2FoldChange > fc_cutoff)
down_n <- sum(!is.na(deg_data$padj) & !is.na(deg_data$log2FoldChange) &
                deg_data$padj < p_cutoff & deg_data$log2FoldChange < -fc_cutoff)

message(sprintf("   [筛选结果] Up = %d, Down = %d, Total = %d",
                up_n, down_n, nrow(signif_degs)))

if (nrow(signif_degs) == 0) {
  write.csv(deg_data, "logs/02_enrichment_debug_deg_data.csv", row.names = FALSE)
  stop("❌ 未筛到显著差异基因，无法进行富集分析")
}

message(sprintf("   [DEGs] 共筛到 %d 个显著差异基因", nrow(signif_degs)))

message("⏳ [3/5] 转换 ENSEMBL -> ENTREZID/SYMBOL ...")
gene_ids <- bitr(
  unique(signif_degs$Ensembl),
  fromType = "ENSEMBL",
  toType   = c("ENTREZID", "SYMBOL"),
  OrgDb    = org.Mm.eg.db
)

bg_ids <- bitr(
  unique(deg_data$Ensembl),
  fromType = "ENSEMBL",
  toType   = "ENTREZID",
  OrgDb    = org.Mm.eg.db
)

if (is.null(gene_ids) || nrow(gene_ids) == 0) stop("❌ 差异基因无法成功映射到 ENTREZID")

parse_ratio <- function(x) {
  sapply(strsplit(x, "/"), function(z) as.numeric(z[1]) / as.numeric(z[2]))
}

wrap_text <- function(x, width = 45) {
  vapply(x, function(y) paste(strwrap(y, width = width), collapse = "\n"), character(1))
}

message("⏳ [4/5] GO 富集分析 ...")
go_res <- enrichGO(
  gene          = unique(gene_ids$ENTREZID),
  universe      = unique(bg_ids$ENTREZID),
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

if (!is.null(go_res) && nrow(as.data.frame(go_res)) > 0) {
  go_df <- as.data.frame(go_res) %>%
    arrange(p.adjust, pvalue) %>%
    mutate(GeneRatioNum = parse_ratio(GeneRatio))

  write.csv(go_df, file.path(out_dir, "GO_Enrichment_Results.csv"), row.names = FALSE)

  go_top <- go_df %>%
    slice_head(n = 10) %>%
    mutate(Description = wrap_text(Description, 45))

  p_go <- ggplot(go_top, aes(x = GeneRatioNum, y = reorder(Description, GeneRatioNum),
                             size = Count, color = ONTOLOGY)) +
    geom_point(alpha = 0.9) +
    theme_bw(base_size = 13) +
    labs(
      title = "Top 10 GO Terms",
      x = "Gene Ratio",
      y = NULL
    ) +
    scale_color_manual(values = c("BP" = "#d73027", "CC" = "#4575b4", "MF" = "#66a61e"))

  ggsave(file.path(out_dir, "02_GO_Dotplot.pdf"), plot = p_go, width = 10, height = 6.5)
} else {
  message("   [GO] 未获得显著富集结果")
}

message("⏳ [5/5] KEGG 富集分析 ...")
kegg_res <- enrichKEGG(
  gene          = unique(gene_ids$ENTREZID),
  universe      = unique(bg_ids$ENTREZID),
  organism      = "mmu",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05
)

if (!is.null(kegg_res) && nrow(as.data.frame(kegg_res)) > 0) {
  kegg_df <- as.data.frame(kegg_res) %>%
    arrange(p.adjust, pvalue) %>%
    mutate(GeneRatioNum = parse_ratio(GeneRatio))

  write.csv(kegg_df, file.path(out_dir, "KEGG_Enrichment_Results.csv"), row.names = FALSE)

  kegg_top <- kegg_df %>%
    slice_head(n = 10) %>%
    mutate(Description = wrap_text(Description, 45))

  p_kegg <- ggplot(kegg_top, aes(x = GeneRatioNum, y = reorder(Description, GeneRatioNum),
                                 size = Count, color = -log10(p.adjust))) +
    geom_point(alpha = 0.9) +
    theme_bw(base_size = 13) +
    labs(
      title = "Top 10 KEGG Pathways",
      x = "Gene Ratio",
      y = NULL,
      color = expression(-log[10](adj.P))
    )

  ggsave(file.path(out_dir, "02_KEGG_Dotplot.pdf"), plot = p_kegg, width = 10, height = 6.5)
} else {
  message("   [KEGG] 未获得显著富集结果")
}

message("🎉======================================================================")
message("🎉 [SUCCESS] GO 与 KEGG 富集分析完成！")
message("🎉======================================================================")
