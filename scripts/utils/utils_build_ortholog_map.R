suppressPackageStartupMessages({
  library(biomaRt)
  library(dplyr)
  library(readr)
})

build_strict_ortholog_map <- function(
  outfile = "data/reference/mouse_human_ortholog_strict.tsv"
) {
  dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)

  mouse <- useEnsembl("genes", dataset = "mmusculus_gene_ensembl")
  
  attrs <- c(
    "ensembl_gene_id",
    "mgi_symbol",
    "external_gene_name",
    "gene_biotype",
    "hsapiens_homolog_ensembl_gene",
    "hsapiens_homolog_associated_gene_name",
    "hsapiens_homolog_orthology_type",
    "hsapiens_homolog_orthology_confidence",
    "hsapiens_homolog_perc_id",
    "hsapiens_homolog_perc_id_r1"
  )

  raw_map <- getBM(
    attributes = attrs,
    mart = mouse
  )

  map1 <- raw_map %>%
    rename(
      mouse_ensembl = ensembl_gene_id,
      mouse_symbol = mgi_symbol,
      mouse_name = external_gene_name,
      mouse_biotype = gene_biotype,
      human_ensembl = hsapiens_homolog_ensembl_gene,
      human_symbol = hsapiens_homolog_associated_gene_name,
      orthology_type = hsapiens_homolog_orthology_type,
      orthology_confidence = hsapiens_homolog_orthology_confidence,
      mouse_to_human_pct = hsapiens_homolog_perc_id,
      human_to_mouse_pct = hsapiens_homolog_perc_id_r1
    ) %>%
    mutate(across(everything(), ~ifelse(.x == "", NA, .x))) %>%
    filter(
      !is.na(mouse_ensembl),
      !is.na(mouse_symbol),
      !is.na(human_ensembl),
      !is.na(human_symbol)
    ) %>%
    filter(
      orthology_type == "ortholog_one2one",
      orthology_confidence == 1
    ) %>%
    distinct()

  # 再做一次双向唯一性约束，避免隐性一对多
  map2 <- map1 %>%
    add_count(mouse_ensembl, name = "n_mouse") %>%
    add_count(human_ensembl, name = "n_human") %>%
    filter(n_mouse == 1, n_human == 1) %>%
    select(-n_mouse, -n_human)

  write_tsv(map2, outfile)
  message("Saved strict ortholog map: ", outfile)
  message("Rows: ", nrow(map2))
  invisible(map2)
}
# ===================================================
# 🌟 核心：配置清华镜像源并运行函数
# ===================================================
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor")

# 执行函数
build_strict_ortholog_map()