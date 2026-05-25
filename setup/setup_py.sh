#!/usr/bin/env bash

echo "======================================"
echo "🚀 创建 Python RNAseq 环境"
echo "======================================"

conda create -n rnaseq_py python=3.11 -y

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate rnaseq_py

echo "📦 安装 Python 包..."

conda install -c conda-forge -y 
pandas 
numpy 
scipy 
matplotlib 
seaborn 
scikit-learn 
openpyxl 
pyyaml 
python-dateutil

echo "🧪 测试 pandas..."

python -c "
import pandas
import numpy
import scipy
print('✅ Python 环境正常')
print('pandas version:', pandas.**version**)
"

echo "======================================"
echo "✅ rnaseq_py 环境安装完成"
echo "======================================"
