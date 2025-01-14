import pandas as pd

# 读取第一个表 (Motif - RBP)，根据mCross生成
motif_rbp = pd.read_csv('motif-RBP.txt', sep='\t', header=None, names=['Motif', 'RBP'])

# 读取第二个表 (Motif - Count)
motif_count = pd.read_csv('DS_exon_motif_results.txt', sep='\t', header=None, names=['Motif', 'Count']) #或为DS_upstream_motif_results.txt，DS_downstream_motif_results.txt
motif_count['Count'] = pd.to_numeric(motif_count['Count'], errors='coerce')
# 合并两个表格，找到Motif对应的Count
merged = pd.merge(motif_rbp, motif_count, on='Motif', how='inner')

# 计算每个RBP对应的Count总和
rbp_count = merged.groupby('RBP')['Count'].sum().reset_index()

# 输出结果
rbp_count.to_csv('DS_RBP_exon_results.txt', sep='\t')   #DS_RBP_upstream_results.txt，DS_RBP_downstream_results.txt
