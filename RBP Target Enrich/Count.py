import pandas as pd

# 读取ExS矩阵
rbp_matrix = pd.read_csv('DExDS.txt', sep = '\t') # 或者其他的ExS矩阵
# 创建一个空的字典来存储每个RBP的activation和repression计数
activation_repression_count = {}

# 遍历计数
for rbp in rbp_matrix.columns[1:]:  # 跳过第一列（剪接事件名称）
    # 统计每个RBP列的 repression 和 activation 的数量
    repression_count = (rbp_matrix[rbp] == -1).sum()
    activation_count = (rbp_matrix[rbp] == 1).sum()

    # 将统计结果存储到字典中
    activation_repression_count[rbp] = {'activation': activation_count, 'repression': repression_count}
    #activation_repression_count[rbp] = {'total': activation_count + repression_count}

# 将字典转换为DataFrame
activation_repression_df = pd.DataFrame.from_dict(activation_repression_count, orient='index')

# 输出结果
activation_repression_df.to_csv('DE_DS_count.txt', sep = '\t')
