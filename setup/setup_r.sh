#!/usr/bin/env bash
set -euo pipefail

echo "======================================"
echo "🚀 创建 RNAseq R 环境"
echo "======================================"

ENV_NAME="rnaseq_r"

# 初始化 conda
source "$(conda info --base)/etc/profile.d/conda.sh"

# 删除旧环境
conda env remove -n ${ENV_NAME} -y || true

# 强制严格 channel 优先级（关键）
conda config --set channel_priority strict

echo "📦 创建基础环境..."

# 直接用 conda-forge 创建
conda create -n ${ENV_NAME} \
    -c conda-forge \
    -c bioconda \
    python=3.11 \
    r-base=4.3 \
    -y

conda activate ${ENV_NAME}

echo "📦 安装 RNAseq / Bioconductor 包..."

conda install -y \
    -c conda-forge \
    -c bioconda \
    bioconductor-deseq2 \
    bioconductor-clusterprofiler \
    bioconductor-org.mm.eg.db \
    bioconductor-enrichplot \
    r-tidyverse \
    r-readr \
    r-data.table \
    r-pheatmap

echo "🧪 测试 R 包..."

Rscript -e "
library(DESeq2)
library(clusterProfiler)
library(org.Mm.eg.db)

cat('✅ DESeq2 OK\n')
cat('✅ clusterProfiler OK\n')
cat('✅ org.Mm.eg.db OK\n')
"

echo "======================================"
echo "✅ RNAseq R 环境安装完成"
echo "======================================"