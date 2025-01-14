import os
import pandas as pd
from collections import defaultdict

# 定义文件夹路径
folder_path = "/path/to/DS/result"

# 初始化存储 RBP 和 Target 的数据结构
rbp_target_map = defaultdict(lambda: defaultdict(list))
rbp_cellline_count = defaultdict(set)  # 记录每个 RBP 出现的细胞系

# 遍历文件夹，筛选目标文件
for file_name in os.listdir(folder_path):
    if file_name.endswith("-01.dn.txt") or file_name.endswith("-01.up.txt"):
        # 提取细胞系名称、RBP名称和调控方向
        parts = file_name.split(".")
        cell_line = parts[0]
        rbp_name = parts[1]
        direction = "dn" if "dn.txt" in file_name else "up"

        # 记录 RBP 所在的细胞系
        rbp_cellline_count[rbp_name].add(cell_line)

        # 读取文件并获取第五列（调控的 target）
        file_path = os.path.join(folder_path, file_name)
        df = pd.read_csv(file_path, sep="\t", header=None)
        targets = df.iloc[:, 4].tolist()

        # 将 target 加入数据结构
        for target in targets:
            rbp_target_map[target][rbp_name].append(direction)

# 初始化结果矩阵
all_targets = sorted(rbp_target_map.keys())
all_rbps = sorted({rbp for target in rbp_target_map for rbp in rbp_target_map[target]})

# 填充矩阵
result_matrix = pd.DataFrame(0, index=all_targets, columns=all_rbps)

for target, rbps in rbp_target_map.items():
    for rbp, directions in rbps.items():
        if len(rbp_cellline_count[rbp]) == 1:  # 如果 RBP 只出现在一个细胞系中

            if str(directions[0]) == 'up':
                result_matrix.loc[target, rbp] = 1
            else:
                result_matrix.loc[target, rbp] = -1
        elif len(set(directions)) == 1:  # 如果 RBP 出现在两个细胞系中，需要调控方向一致
            if str(directions[0]) == 'up':
                result_matrix.loc[target, rbp] = 1
            else:
                result_matrix.loc[target, rbp] = -1

bed_file = pd.read_csv('../data/Hs.seq.all.cass.chrom.can.exon.bed', sep='\t', header=None)  # 读取BED文件，假设以tab分隔

# 提取BED文件中的剪接事件名称
bed_events = bed_file[3].unique()  # 提取（剪接事件名称）

# 从RBP行为矩阵中筛选出共有的剪接事件行
filtered_matrix = result_matrix[result_matrix.index.isin(bed_events)]

# 保存结果为 tab 分隔的文件
output_path = "ExS.txt"
filtered_matrix.to_csv(output_path, sep="\t")
