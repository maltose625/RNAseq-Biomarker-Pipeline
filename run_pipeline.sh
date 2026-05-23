#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$PROJECT_DIR"

ENV_NAME="${CONDA_ENV_NAME:-rnaseq_enrich}"
LOG_DIR="logs"
CONFIG_FILE="config/config.yaml"

mkdir -p "$LOG_DIR"

timestamp() {
  date '+%Y-%m-%d %H:%M:%S'
}

log_info() {
  echo "[$(timestamp)] [INFO] $*"
}

log_warn() {
  echo "[$(timestamp)] [WARN] $*"
}

log_error() {
  echo "[$(timestamp)] [ERROR] $*" >&2
}

activate_conda() {
  if command -v conda >/dev/null 2>&1; then
    eval "$(conda shell.bash hook)"
  elif [[ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]]; then
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
  elif [[ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]]; then
    source "$HOME/anaconda3/etc/profile.d/conda.sh"
  else
    log_error "未找到 conda，请先确认 conda 已安装"
    exit 1
  fi

  conda activate "$ENV_NAME"
  log_info "已激活 conda 环境: $ENV_NAME"
}

check_file() {
  local f="$1"
  if [[ ! -f "$f" ]]; then
    log_error "缺少文件: $f"
    exit 1
  fi
}

run_stage() {
  local stage_name="$1"
  local script_path="$2"
  local log_file="$3"
  local critical="${4:-yes}"

  log_info "开始执行: $stage_name"
  log_info "脚本: $script_path"

  if [[ ! -f "$script_path" ]]; then
    if [[ "$critical" == "yes" ]]; then
      log_error "脚本不存在: $script_path"
      exit 1
    else
      log_warn "脚本不存在，跳过: $script_path"
      return 0
    fi
  fi

  set +e
  Rscript "$script_path" 2>&1 | tee "$log_file"
  local exit_code=${PIPESTATUS[0]}
  set -e

  if [[ $exit_code -ne 0 ]]; then
    if [[ "$critical" == "yes" ]]; then
      log_error "$stage_name 执行失败，退出码: $exit_code"
      exit $exit_code
    else
      log_warn "$stage_name 执行失败，但设为非关键步骤，继续后续流程"
    fi
  else
    log_info "$stage_name 执行完成"
  fi
}

verify_output() {
  local output_file="$1"
  local label="$2"
  if [[ -f "$output_file" ]]; then
    log_info "$label 输出检查通过: $output_file"
  else
    log_error "$label 输出缺失: $output_file"
    exit 1
  fi
}

main() {
  log_info "🚀 启动总控流程"

  check_file "$CONFIG_FILE"
  activate_conda

  # Stage 1: DESeq2
  run_stage \
    "Stage 1 - DESeq2差异分析" \
    "scripts/01_deseq2_analysis.R" \
    "$LOG_DIR/01_deseq2_analysis.log" \
    "yes"

  verify_output \
    "results/01_DESeq2/DESeq2_Annotated_Results.csv" \
    "DESeq2"

  # Stage 2: Enrichment
  # 如果你希望富集失败不影响后面的机器学习，可把 yes 改成 no
  run_stage \
    "Stage 2 - 富集分析" \
    "scripts/02_enrichment.R" \
    "$LOG_DIR/02_enrichment.log" \
    "yes"

  # Stage 3: Mouse LASSO
  run_stage \
    "Stage 3 - 小鼠特征筛选" \
    "scripts/03_lasso_mouse.R" \
    "$LOG_DIR/03_lasso_mouse.log" \
    "yes"

  MOUSE_FS_DIR="results/03_MachineLearning/mouse_feature_selection"
  CORE_1SE="$MOUSE_FS_DIR/Core_Biomarkers_1se.csv"
  CORE_MIN="$MOUSE_FS_DIR/Core_Biomarkers_min.csv"
  CAND_FILE="$MOUSE_FS_DIR/Candidate_Genes_for_ML.csv"
  SKIP_FILE="$MOUSE_FS_DIR/ML_SKIP_REASON.txt"

  if [[ -f "$CORE_1SE" ]]; then
    log_info "小鼠特征筛选输出检查通过: $CORE_1SE"
    export MOUSE_GENE_MODE="core"
    export MOUSE_GENE_FILE="$CORE_1SE"
  elif [[ -f "$CORE_MIN" ]]; then
      log_info "小鼠特征筛选输出检查通过: $CORE_MIN"
      export MOUSE_GENE_MODE="core"
      export MOUSE_GENE_FILE="$CORE_MIN"
  elif [[ -f "$CAND_FILE" && -f "$SKIP_FILE" ]]; then
      log_warn "小鼠样本量不足，Stage 3 进入降级模式；将使用候选基因池继续执行人类验证"
      log_info "候选池文件: $CAND_FILE"
      export MOUSE_GENE_MODE="candidate"
      export MOUSE_GENE_FILE="$CAND_FILE"
  else
    log_error "小鼠特征筛选输出缺失：既没有 Core_Biomarkers，也没有 Candidate_Genes_for_ML.csv"
    exit 1 
  fi

  # Stage 4: Human validation
  run_stage \
    "Stage 4 - 人类验证" \
    "scripts/04_human_validation.R" \
    "$LOG_DIR/04_human_validation.log" \
    "yes"

  log_info "✅ 全部流程执行完成"
  log_info "结果目录:"
  log_info "  - results/01_DESeq2"
  log_info "  - results/02_Enrichment"
  log_info "  - results/03_MachineLearning"
}

main "$@"
