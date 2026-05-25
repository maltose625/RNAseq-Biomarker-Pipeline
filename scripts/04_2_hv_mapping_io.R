################################################
##### 04_2_hv_mapping_io.R (跨物种映射模块)######
################################################

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
  skip_file <- file.path(
    config$paths$results_ml %||% "results/03_MachineLearning",
    "mouse_feature_selection",
    "ML_SKIP_REASON.txt"
  )

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
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    log_warn("未找到 DESeq2 结果文件用于 Ensembl->Symbol 辅助映射: %s", path %||% "NULL")
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

  map
}

build_mouse_to_human_mapping <- function(mouse_genes, expr_row_ids, ortholog_file = NULL, mouse_deg_file = NULL) {
  expr_id_type <- detect_id_type(expr_row_ids)
  mouse_id_type <- detect_id_type(mouse_genes)

  log_info("候选基因ID类型: %s", mouse_id_type)
  log_info("人类表达矩阵行名ID类型: %s", expr_id_type)

  deg_map <- read_mouse_deg_map(mouse_deg_file)
  orth <- read_ortholog_map(ortholog_file)

  orth$mouse_ensembl <- strip_version(as.character(orth$mouse_ensembl))
  orth$human_ensembl <- strip_version(as.character(orth$human_ensembl))
  orth$mouse_symbol  <- trimws(as.character(orth$mouse_symbol))
  orth$human_symbol  <- trimws(as.character(orth$human_symbol))

  candidate_df <- data.frame(
    MouseGene = unique(trimws(as.character(mouse_genes))),
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

  if (mouse_id_type == "mouse_ensembl") {
    candidate_df$JoinKey <- strip_version(candidate_df$MouseGene)
    orth$JoinKey <- orth$mouse_ensembl
  } else {
    candidate_df$JoinKey <- toupper(trimws(candidate_df$MouseSymbol))
    orth$JoinKey <- toupper(trimws(orth$mouse_symbol))
  }

  target_human_col <- if (expr_id_type == "human_ensembl") "human_ensembl" else "human_symbol"

  mapped <- merge(
    candidate_df,
    orth[, c("JoinKey", "mouse_ensembl", "mouse_symbol", "human_ensembl", "human_symbol"), drop = FALSE],
    by = "JoinKey",
    all.x = FALSE,
    all.y = FALSE
  )

  if (nrow(mapped) == 0) {
    log_stop("❌ 严格直系同源映射未匹配到任何基因！请检查离线同源表或基因 ID 是否对应。")
  }

  mapped$HumanGene <- mapped[[target_human_col]]
  mapped <- mapped[!is.na(mapped$HumanGene) & nzchar(mapped$HumanGene), , drop = FALSE]
  mapped <- unique(mapped)

  expr_ids_norm <- if (expr_id_type == "human_ensembl") {
    strip_version(expr_row_ids)
  } else {
    trimws(expr_row_ids)
  }

  available_human_genes <- unique(mapped$HumanGene[mapped$HumanGene %in% expr_ids_norm])
  missing_human_genes   <- setdiff(unique(mapped$HumanGene), available_human_genes)
  unmapped_mouse_genes  <- setdiff(unique(mouse_genes), unique(mapped$MouseGene))

  log_info("使用离线 ortholog_map（%s）成功映射 %d 条", mouse_id_type, nrow(mapped))
  log_info("映射后人类候选基因数: %d", length(unique(mapped$HumanGene)))
  log_info("在人类开发队列表达矩阵中成功匹配基因数: %d", length(available_human_genes))
  log_info("未匹配基因数: %d", length(missing_human_genes))

  list(
    mouse_id_type = mouse_id_type,
    expr_id_type = expr_id_type,
    mapping_table = mapped,
    mapped_human_genes = unique(mapped$HumanGene),
    available_human_genes = available_human_genes,
    missing_human_genes = missing_human_genes,
    unmapped_mouse_genes = unmapped_mouse_genes
  )
}
