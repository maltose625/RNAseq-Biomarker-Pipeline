#!/usr/bin/env Rscript
# =========================================================================
# 脚本名称: mouse-human_change.R
# 模块功能: 解析本地 MGI 同源表，离线构建严格 1:1 鼠人同源密码本
# =========================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
})

cat("[INFO] 🚀 开始离线解析 MGI 同源数据库...\n")

mgi_file <- "data/reference/HOM_MouseHumanSequence.rpt"
outfile <- "data/reference/mouse_human_ortholog_strict.tsv"

if (!file.exists(mgi_file)) {
  stop(sprintf("❌ 未找到 MGI 原始文件：%s，请检查路径。", mgi_file))
}

# 1. 读取 MGI 原始长表
mgi <- read_tsv(mgi_file, show_col_types = FALSE)

# 2. 规范化列名
colnames(mgi) <- c("class_key", "organism", "taxon_id", "symbol", "entrez_id", 
                   "mgi_id", "hgnc_id", "omim_id", "location", "coordinates", 
                   "refseq_nucl", "refseq_prot", "swiss_prot")

cat("[INFO] 正在提取小鼠与人类直系同源配对...\n")

# 3. 强制使用 dplyr:: 命名空间，防止 AnnotationDbi 冲突覆盖
mouse_df <- mgi %>%
  dplyr::filter(organism == "mouse, laboratory") %>%
  dplyr::select(class_key, mouse_symbol = symbol)

human_df <- mgi %>%
  dplyr::filter(organism == "human") %>%
  dplyr::select(class_key, human_symbol = symbol)

# 4. 内连接
ortholog_pairs <- dplyr::inner_join(mouse_df, human_df, by = "class_key")

cat("[INFO] 正在执行严格的跨物种 1:1 双向唯一性约束筛选...\n")

# 5. 核心过滤
strict_pairs <- ortholog_pairs %>%
  dplyr::add_count(mouse_symbol, name = "n_mouse") %>%
  dplyr::add_count(human_symbol, name = "n_human") %>%
  dplyr::filter(n_mouse == 1, n_human == 1) %>%
  dplyr::select(mouse_symbol, human_symbol)

cat("[INFO] 🌟 正在调用本地 org.Mm.eg.db 数据库为小鼠 Symbol 补全 Ensembl ID...\n")

# 6. 核心对齐
strict_pairs$mouse_ensembl <- mapIds(
  org.Mm.eg.db,
  keys      = strict_pairs$mouse_symbol,
  column    = "ENSEMBL",
  keytype   = "SYMBOL",
  multiVals = "first"
)

# 7. 过滤并构建最终矩阵
final_ortholog_matrix <- strict_pairs %>%
  dplyr::filter(!is.na(mouse_ensembl)) %>%
  dplyr::distinct(mouse_ensembl, .keep_all = TRUE) %>%
  dplyr::mutate(
    human_ensembl = NA_character_,
    orthology_type = "ortholog_one2one",
    orthology_confidence = 1
  ) %>%
  dplyr::select(mouse_ensembl, mouse_symbol, human_ensembl, human_symbol, 
         orthology_type, orthology_confidence)

# 8. 输出
dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
readr::write_tsv(final_ortholog_matrix, outfile)

cat(sprintf("🎉 [SUCCESS] 严格 1:1 同源密码本离线生成成功！\n"))
cat(sprintf("📂 输出路径: %s\n", outfile))
cat(sprintf("📊 共有 %d 个高置信度鼠-人基因对成功锁入密码本。\n", nrow(final_ortholog_matrix)))