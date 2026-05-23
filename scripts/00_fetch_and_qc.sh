#!/usr/bin/env bash
# =============================================================================
# 脚本名称：00_fetch_and_qc.sh
# 功能描述：统一测序数据处理流水线 (数据下载 -> fastp 质控，极限 I/O 优化版)
# =============================================================================

# 开启严格错误检测模式 (工业级标配)
set -euo pipefail
shopt -s nullglob

echo "================================================="
echo "🚀 启动自动化数据获取与质控流水线..."
echo "================================================="

# --- 1. 动态相对路径配置 ---
# 自动定位项目根目录，彻底告别 D 盘硬编码
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# 定义 Linux 原生数据存放路径
RAW_DIR="${PROJECT_DIR}/data/raw_fastq"
QC_DIR="${PROJECT_DIR}/data/clean_counts"
LIST_FILE="${PROJECT_DIR}/config/srr_list.txt" # 建议把样本列表放在 config 目录下

# 确保工作目录存在
mkdir -p "${RAW_DIR}" "${QC_DIR}"

# --- 2. 命令行参数解析 ---
STEP="all" # 默认执行全部
METHOD="direct" # 默认下载方式

usage() {
    echo "用法: bash $0 [-s <all|fetch|qc>] [-m <direct|kingfisher>]"
    echo "  -s : 执行步骤 [all|fetch|qc] (默认: all)"
    echo "  -m : 下载模式 [direct|kingfisher] (默认: direct)"
    exit 1
}

while getopts "s:m:h" opt; do
    case $opt in
        s) STEP="$OPTARG" ;;
        m) METHOD="$OPTARG" ;;
        h) usage ;;
        ?) usage ;;
    esac
done

echo "⚙️ 当前执行阶段: [ $STEP ] | 下载模式: [ $METHOD ]"

# =============================================================================
# 模块一：数据获取 (Fetch)
# =============================================================================
if [[ "$STEP" == "all" || "$STEP" == "fetch" ]]; then
    echo "-------------------------------------------------"
    echo "🌐 [模块一: 数据获取] 启动..."
    
    if [[ ! -f "$LIST_FILE" ]]; then
        echo "❌ 错误：未找到样本列表文件: $LIST_FILE"
        echo "请在 config 目录下创建 srr_list.txt 并填入样本 SRR 号。"
        exit 1
    fi

    cd "${RAW_DIR}"

    if [[ "$METHOD" == "direct" ]]; then
        echo "⬇️ [wget 直下模式]..."
        # 假设 srr_list.txt 里面是直接的下载链接
        wget -i "$LIST_FILE" -c -q --show-progress
    elif [[ "$METHOD" == "kingfisher" ]]; then
        echo "🦅 [Kingfisher 级联爬取模式]..."
        while IFS= read -r id || [[ -n "$id" ]]; do
            [[ -z "$id" || "$id" =~ ^# ]] && continue
            
            # 智能防呆：检查是否已存在压缩包
            if [[ -f "${id}_1.fastq.gz" || -f "${id}.fastq.gz" ]]; then
                 echo "⏭️ [跳过] $id 原始数据已存在。"
                 continue
            fi
            
            echo "⬇️ 正在获取: $id ..."
            kingfisher get -r "$id" -m ena-ftp aws-http prefetch -f fastq.gz --download-threads 8
        done < "$LIST_FILE"
    else
        echo "❌ 未知下载模式：$METHOD"
        exit 1
    fi
    echo "✅ 数据获取阶段完成！"
fi

# =============================================================================
# 模块二：数据质控 (QC) - 极限 I/O 优化
# =============================================================================
if [[ "$STEP" == "all" || "$STEP" == "qc" ]]; then
    echo "-------------------------------------------------"
    echo "🧪 [模块二: 自动化 QC 流水线] 启动..."
    
    cd "${RAW_DIR}"
    
    # 自动识别所有的 R1 压缩包
    files=(*_1.fastq.gz)
    
    if [[ ${#files[@]} -eq 0 ]]; then
        echo "⚠️ 警告：${RAW_DIR} 中没有找到任何 _1.fastq.gz 结尾的文件，跳过质控。"
    else
        for r1 in "${files[@]}"; do
            sample_name="${r1%_1.fastq.gz}"
            r2="${sample_name}_2.fastq.gz"
            
            # 质控输出依然保持 .gz 压缩状态！绝不产生庞大的未压缩 .fq
            clean_r1="${QC_DIR}/${sample_name}_clean_1.fq.gz"
            clean_r2="${QC_DIR}/${sample_name}_clean_2.fq.gz"

            if [[ -f "${r2}" ]]; then
                if [[ -f "${clean_r1}" && -f "${clean_r2}" ]]; then
                    echo "⏭️ [跳过] ${sample_name} 干净数据已存在。"
                else
                    echo "📊 正在质控: ${sample_name} ..."
                    # fastp 原生支持直接读取和输出 .gz 文件
                    fastp -i "${r1}" -I "${r2}" \
                          -o "${clean_r1}" -O "${clean_r2}" \
                          -j "${QC_DIR}/${sample_name}_fastp.json" \
                          -h "${QC_DIR}/${sample_name}_fastp.html" \
                          --thread 4 --detect_adapter_for_pe --cut_right \
                          --length_required 50 --qualified_quality_phred 15
                fi
            else
                echo "⚠️ 警告: 找不到 ${sample_name} 的 R2 配对文件！"
            fi
        done

        # 汇总质控报告
        if command -v multiqc &> /dev/null; then
            echo "📈 生成 MultiQC 全局质控报告..."
            multiqc "${QC_DIR}" -o "${QC_DIR}/MultiQC_Report"
        else
            echo "⚠️ 警告: 未安装 MultiQC，跳过全局报告生成。建议在 conda 环境中安装。"
        fi
    fi
    echo "✅ 质控阶段完成！"
fi

echo "================================================="
echo "🎉 数据获取与质控流水线圆满结束！"
echo "================================================="