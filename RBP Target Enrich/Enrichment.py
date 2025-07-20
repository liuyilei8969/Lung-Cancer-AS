import pandas as pd
from scipy.stats import hypergeom
import numpy as np
from statsmodels.stats.multitest import multipletests

# 读取数据
foreground = pd.read_csv('DE_DS_count.txt', sep='\t', index_col=0)
background = pd.read_csv('E_S_count.txt', sep='\t', index_col=0)
rbps = foreground.index.intersection(background.index)

# 设置总数
total_foreground = 127  # 对于Subtype-Specific，替换为178
total_background = 10931  # ExS matrix中的总行数

# 计算超几何分布检验
def hypergeometric_test(fg_count, bg_count, total_fg, total_bg):
    M = total_fg + total_bg  # 总体大小
    n = fg_count + bg_count  # 成功数
    N = total_fg  # 目标群体大小
    x = fg_count  # 成功事件数
    p_value = hypergeom.sf(x - 1, M, N, n)  # 计算上尾概率
    return p_value

# 计算fold change分析
def fold_change(fg_count, bg_count):
    return fg_count * total_background / ((bg_count if bg_count > 0 else np.nan) * total_foreground)  # 避免除以零

# 计算每个RBP的分析结果
results = []

for rbp in rbps:
    fg_count = foreground.loc[rbp, 'total']
    bg_count = background.loc[rbp, 'total']

    # 超几何检验
    p_value_hypergeometric = hypergeometric_test(fg_count, bg_count, total_foreground, total_background)
    
    # Fold change分析
    fc = fold_change(fg_count, bg_count)

    # 存储结果
    results.append({
        'RBP': rbp,
        'DScount': fg_count,
        'count': bg_count,
        'hypergeo_p_value': p_value_hypergeometric,
        'fold_change': fc
    })

# 将结果保存为DataFrame
results_df = pd.DataFrame(results)

# 多重检验校正：使用Benjamini-Hochberg方法
pvals = results_df['hypergeo_p_value']
_, corrected_pvals, _, _ = multipletests(pvals, method='fdr_bh')

# 将校正后的p值加入结果
results_df['corrected_p_value'] = corrected_pvals

# 可选择将结果保存到文件
results_df.to_csv('RBPactivity_corrected.txt', index=False, sep='\t')

print("Analysis complete, results saved to 'RBPactivity_corrected.txt'.")
