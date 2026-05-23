#!/usr/bin/env bash
# =================================================================
# 脚本名称: 00_upstream_align_quant.sh
# 核心功能: HISAT2 流式比对 + featureCounts 全局定量 (极限 I/O 优化版)
# =================================================================

# 开启严格报错模式
set -euo pipefail
shopt -s nullglob

echo "================================================="
echo "🚀 启动上游自动化流水线 (比对 & 定量) ..."
echo "================================================="

# 1. 定义动态相对路径 (基于我们新建的 Linux 目录结构)
# 获取当前脚本所在目录的上一级，即项目根目录 ~/RNAseq_Pipeline
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# 所有的输入输出全部指向 Linux 原生系统 (软链接也在这里)
QC_DIR="${PROJECT_DIR}/data/clean_counts"  # 假设质控后的 .gz 文件软链接在此
ALIGN_DIR="${PROJECT_DIR}/results/00_Alignment"
QUANT_DIR="${PROJECT_DIR}/results/00_Quantification"

mkdir -p "${ALIGN_DIR}" "${QUANT_DIR}"

# 2. 外部数据库路径 (这个可以指向 D 盘的软链接，或者 Linux 里的路径)
INDEX="$HOME/Databases/GRCm39/GRCm39_hisat2"
GTF="/mnt/d/A/WSL_Microbiome_Project/Databases/Mus_musculus/GRCm39/Mus_musculus.GRCm39.111.gtf" # 建议后续也改成软链接

# 3. 动态获取样本列表 (告别手动写 SRR 名字，自动识别文件夹里的配对文件)
# 这行代码会自动寻找 QC_DIR 下所有 _clean_1.fq.gz 文件，并提取出样本名
SAMPLES=($(ls "${QC_DIR}"/*_clean_1.fq.gz | awk -F"/" '{print $NF}' | sed 's/_clean_1.fq.gz//'))

if [ ${#SAMPLES[@]} -eq 0 ]; then
    echo "❌ 错误: 在 ${QC_DIR} 未找到任何 clean_1.fq.gz 文件！"
    exit 1
fi

echo "📦 检测到待处理样本: ${SAMPLES[*]}"

# =================================================================
# 模块一：批量流式比对 (全程不落盘临时文件)
# =================================================================
for SAMPLE in "${SAMPLES[@]}"; do
    echo "-------------------------------------------------"
    echo "🧬 正在处理样本: ${SAMPLE} ..."
    
    # 智能防呆机制
    if [ -f "${ALIGN_DIR}/${SAMPLE}_sorted.bam" ]; then
        echo "⏭️ [跳过] 发现 ${SAMPLE}_sorted.bam 已存在。"
        continue
    fi

    # 极限流式比对核心代码：直接读 .gz -> 内存比对 -> 内存转 BAM -> 内存排序 -> 仅输出最终 BAM
    echo "⚙️ 启动 HISAT2 流式比对 (直接读取 .gz，杜绝磁盘暴涨)..."
    hisat2 -p 4 \
        -x "${INDEX}" \
        -1 "${QC_DIR}/${SAMPLE}_clean_1.fq.gz" \
        -2 "${QC_DIR}/${SAMPLE}_clean_2.fq.gz" \
        --summary-file "${ALIGN_DIR}/${SAMPLE}_summary.txt" \
        | samtools view -bS - \
        | samtools sort -@ 4 -m 1G \
        # 新增：建立 BAM 索引
    echo "🔍 正在为 BAM 文件建立索引..."
    samtools index "${ALIGN_DIR}/${SAMPLE}_sorted.bam"
          
    echo "✅ 样本 ${SAMPLE} BAM 构建与索引完成！"
          -o "${ALIGN_DIR}/${SAMPLE}_sorted.bam"
          
    echo "✅ 样本 ${SAMPLE} BAM 构建完成！"
done

# =================================================================
# 模块二：全局定量
# =================================================================
echo "================================================="
echo "📊 开始生成全局表达矩阵 (featureCounts) ..."
echo "================================================="

# 动态获取所有刚刚生成的 BAM 文件
BAM_FILES=("${ALIGN_DIR}"/*_sorted.bam)

# 喂给 featureCounts 一次性统计
featureCounts -T 4 -p \
    -t exon -g gene_id \
    -a "${GTF}" \
    -o "${QUANT_DIR}/All_Samples_counts.txt" \
    "${BAM_FILES[@]}"

# 提取核心数据列 (去除特征信息，只保留 GeneID 和 Counts) 方便后续 DESeq2 直接使用
cut -f1,7- "${QUANT_DIR}/All_Samples_counts.txt" | grep -v "^#" > "${PROJECT_DIR}/data/clean_counts/Final_Counts_Matrix.txt"

echo "🎉 [SUCCESS] 上游分析全部竣工！最终表达矩阵已输出至: data/clean_counts/Final_Counts_Matrix.txt"