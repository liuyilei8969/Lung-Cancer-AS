import pandas as pd
from scipy.stats import hypergeom, binom
import numpy as np
from scipy.stats import chi2_contingency
from statsmodels.stats.multitest import multipletests

# 读取数据
foreground = pd.read_csv('DE_DS_count.txt', sep='\t', index_col=0)
background = pd.read_csv('../result/count.txt', sep='\t', index_col=0)
rbps = foreground.index.intersection(background.index)
# 设置总数
total_foreground = 178
total_background = 10931

# 计算超几何分布检验
def hypergeometric_test(fg_count, bg_count, total_fg, total_bg):
    M = total_fg + total_bg  # 总体大小
    n = fg_count + bg_count  # 成功数
    N = total_fg  # 目标群体大小
    x = fg_count  # 成功事件数
    p_value = hypergeom.sf(x - 1, M, N, n)  # 计算上尾概率
    return p_value

# 计算二项分布检验
def binomial_test(fg_count, total_fg, prob):
    p_value = binom.sf(fg_count - 1, total_fg, prob)  # 上尾概率
    return p_value

# 计算fold change分析
def fold_change(fg_count, bg_count):
    return fg_count * total_background / ((bg_count if bg_count > 0 else np.nan) * total_foreground ) # 避免除以零

# 计算对数似然比检验 (Log-Likelihood Ratio Test)
def log_likelihood_ratio(fg_count, bg_count, total_fg, total_bg):
    # 如果fg_count或bg_count为零，避免计算
    if fg_count == 0 or bg_count == 0:
        return np.nan  # 返回NaN，表示无法计算
    fg_repression = total_fg - fg_count
    bg_repression = total_bg - bg_count
    observed = np.array([[fg_count, fg_repression], [bg_count, bg_repression]])
    chi2, p_val, _, _ = chi2_contingency(observed)
    return p_val


# 计算每个RBP的分析结果
results = []

for rbp in rbps:
    fg_count = foreground.loc[rbp, 'total']
    bg_count = background.loc[rbp, 'total']

    # 超几何检验
    p_value_hypergeometric = hypergeometric_test(fg_count, bg_count, total_foreground, total_background)

    # 二项分布检验
    binom_p_value = binomial_test(fg_count, total_foreground, fg_count / total_foreground)

    # Fold change分析
    fc = fold_change(fg_count, bg_count)

    # 对数似然比检验
    log_lr_p_value = log_likelihood_ratio(fg_count, bg_count, total_foreground, total_background)

    # 存储结果
    results.append({
        'RBP': rbp,
        'DScount': fg_count,
        'count': bg_count,
        'hypergeo_p_value': p_value_hypergeometric,
        'binom_p_value': binom_p_value,
        'fold_change': fc,
        'log_lr_p_value': log_lr_p_value,
    })

# 将结果保存为DataFrame
results_df = pd.DataFrame(results)



# 可选择将结果保存到文件
results_df.to_csv('../result/subtype/RBPactivity.txt', index=False, sep = '\t')
