#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def parse_arguments():
    """使用 argparse 处理命令行参数 (技能点 1: 命令行编程)"""
    parser = argparse.ArgumentParser(description="RNA-seq Pipeline 结果自动化汇总与出图工具")
    parser.add_argument("-i", "--input_dir", type=str, default="results", 
                        help="分析结果根目录的路径 (默认: results)")
    parser.add_argument("-o", "--output_dir", type=str, default="results/00_Summary", 
                        help="汇总报告输出路径 (默认: results/00_Summary)")
    return parser.parse_args()

def extract_deg_stats(results_dir):
    """使用 pandas 处理差异分析表格 (技能点 2&3: os 路径管理 & pandas)"""
    deseq_file = os.path.join(results_dir, "01_DESeq2", "DESeq2_Annotated_Results.csv")
    
    if not os.path.exists(deseq_file):
        print(f"[WARN] 未找到差异分析文件: {deseq_file}")
        return None
    
    # 用 pandas 读取 csv
    df = pd.read_csv(deseq_file)
    
    # 过滤出显著上调和下调的基因数
    if "Significance" in df.columns:
        up_count = len(df[df["Significance"] == "Up-regulated"])
        down_count = len(df[df["Significance"] == "Down-regulated"])
        return {"Up-regulated": up_count, "Down-regulated": down_count}
    return None

def extract_ml_metrics(results_dir):
    """提取机器学习模块的 AUC 结果"""
    metrics_file = os.path.join(results_dir, "03_MachineLearning", "human_validation", "model_metrics.csv")
    
    if not os.path.exists(metrics_file):
        print(f"[WARN] 未找到机器学习指标文件: {metrics_file}")
        return None
    
    # 用 pandas 读取并转换为字典
    df = pd.read_csv(metrics_file)
    # 将 Metric 列作为键，Value 列作为值
    metrics_dict = dict(zip(df["Metric"], df["Value"]))
    
    # 安全地提取我们要的 AUC 分数
    try:
        dev_auc = float(metrics_dict.get("DevelopmentApparentAUC", 0))
        ext_auc = float(metrics_dict.get("ExternalAUC", 0))
        return {"Development AUC": dev_auc, "External AUC": ext_auc}
    except ValueError:
        return None

def plot_summary_dashboard(deg_stats, auc_stats, output_dir):
    """使用 seaborn 和 matplotlib 生成数据面板图 (技能点 4: 数据可视化)"""
    os.makedirs(output_dir, exist_ok=True)
    
    # 设置图形样式与大小
    sns.set_theme(style="whitegrid")
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # 图 1: 差异基因条形图
    if deg_stats:
        deg_labels = list(deg_stats.keys())
        deg_values = list(deg_stats.values())
        sns.barplot(x=deg_labels, y=deg_values, ax=axes[0], palette=["#d62728", "#1f77b4"])
        axes[0].set_title("Mouse Differentially Expressed Genes (DEGs)", fontsize=14, pad=10)
        axes[0].set_ylabel("Number of Genes")
        # 在柱子上打上数字标签
        for i, v in enumerate(deg_values):
            axes[0].text(i, v + (max(deg_values)*0.02), str(v), ha='center', fontweight='bold')
    
    # 图 2: 模型 AUC 条形图
    if auc_stats:
        auc_labels = list(auc_stats.keys())
        auc_values = list(auc_stats.values())
        sns.barplot(x=auc_labels, y=auc_values, ax=axes[1], palette=["#2ca02c", "#ff7f0e"])
        axes[1].set_title("Human Cohort Validation Performance (AUC)", fontsize=14, pad=10)
        axes[1].set_ylabel("AUC Score")
        axes[1].set_ylim(0, 1.1)
        # 在柱子上打上数字标签
        for i, v in enumerate(auc_values):
            axes[1].text(i, v + 0.02, f"{v:.3f}", ha='center', fontweight='bold')

    plt.tight_layout()
    
    # 保存出图
    plot_path = os.path.join(output_dir, "Project_Final_Dashboard.pdf")
    plt.savefig(plot_path, format="pdf", dpi=300)
    print(f"[INFO] 📊 项目终极汇总图表已生成: {plot_path}")

def main():
    args = parse_arguments()
    print(f"--- 启动 Python 自动化汇总工具 ---")
    print(f"[INFO] 扫描结果目录: {args.input_dir}")
    
    deg_stats = extract_deg_stats(args.input_dir)
    auc_stats = extract_ml_metrics(args.input_dir)
    
    if deg_stats or auc_stats:
        plot_summary_dashboard(deg_stats, auc_stats, args.output_dir)
    else:
        print("[ERROR] 未能提取到足够的数据用于出图，请检查 results 目录的完整性。")

if __name__ == "__main__":
    main()