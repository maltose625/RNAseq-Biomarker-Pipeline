#!/usr/bin/env bash

mkdir -p envs

echo "📦 导出 rnaseq_py ..."
conda env export -n rnaseq_py > envs/rnaseq_py.yml

echo "📦 导出 rnaseq_r ..."
conda env export -n rnaseq_r > envs/rnaseq_r.yml

echo "✅ 所有环境已导出"


conda install -c conda-forge \
  tidyverse \
  data.table \
  ggplot2 \
  r-reshape2 \
  r-viridis \
  -y