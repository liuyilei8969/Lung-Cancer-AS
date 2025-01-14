import scipy.stats as stats
import math
import pandas as pd

def read_data(fg_file, bg_file):
    """
    读取前景和背景文件，返回数据框。
    """
    fg_data = pd.read_csv(fg_file, sep="\t")
    bg_data = pd.read_csv(bg_file, sep="\t")
    return fg_data, bg_data

def calculate_enrichment(fg_data, bg_data):
    """
    计算富集分析结果。
    参数：
        fg_data: 前景数据，包含Motif和Count列。
        bg_data: 背景数据，包含Motif和Count列。
        test_method: 使用的统计方法（binom或hyperg）。
    返回：
        结果数据框，包含Motif、fg_count、bg_count、log2FC、pval、zscore列。
    """
    # 合并fg和bg数据
    merged_data = pd.merge(fg_data, bg_data, on="RBP", how="outer", suffixes=("_fg", "_bg")).fillna(0)
    fg_sum = merged_data["Count_fg"].sum()
    bg_sum = merged_data["Count_bg"].sum()

    results = []

    for _, row in merged_data.iterrows():
        RBP = row["RBP"]
        fg_count = row["Count_fg"]
        bg_count = row["Count_bg"]

        # 计算log2FC
        log2FC = (fg_sum == 0 or bg_count - fg_count == 0 or fg_count == 0 or bg_sum - fg_sum == 0)
        log2FC = math.log2(fg_count * (bg_sum - fg_sum) / (fg_sum * (bg_count - fg_count))) if not log2FC else 'NA'

        # 计算zscore
        r = fg_sum / bg_sum
        zscore = (fg_count - bg_count * r) / math.sqrt(bg_count * r * (1 - r))

        # 计算p值

        pval_hyperg = stats.hypergeom.sf(fg_count - 1, bg_sum, bg_count, fg_sum)


        results.append({
            "RBP": RBP,
            "fg_count": fg_count,
            "bg_count": bg_count,
            "log2FC": log2FC,
            "zscore": zscore,
            "hyperg": pval_hyperg
        })

    return pd.DataFrame(results)

def save_results(results, output_file):
    """
    保存结果到文件。
    """
    results.to_csv(output_file, sep="\t", index=False)

if __name__ == "__main__":
    fg_file = "DS_RBP_exon_results.txt"  # 前景文件路径，或为DS_RBP_downstream_results.txt，DS_RBP_upstream_results.txt
    bg_file = "All_RBP_exon_results.txt"  # 背景文件路径，或为All_RBP_downstream_results.txt，All_RBP_upstream_results.txt
    output_file = "DS_exon_enrich_RBP.txt"  #DS_downstream_enrich_RBP.txt，DS_upstream_enrich_RBP.txt

    # 读取数据
    fg_data, bg_data = read_data(fg_file, bg_file)

    # 计算富集分析结果
    results = calculate_enrichment(fg_data, bg_data)

    # 保存结果
    save_results(results, output_file)

    print(f"Enrichment analysis results saved to {output_file}")
