import pandas as pd

# 读取相关文件
ExS = pd.read_csv('ExS.txt', sep='\t') 

# Correctly reading down-regulated files for events and RBPs
up_events = set(pd.read_csv('up_events.txt', header=None)[0])  
dn_events = set(pd.read_csv('dn_events.txt', header=None)[0])  # Fixed to read 'dn_events.txt'
up_rbps = set(pd.read_csv('up_RBP.txt', header=None)[0]) 
dn_rbps = set(pd.read_csv('dn_RBP.txt', header=None)[0])  # Fixed to read 'dn_RBP.txt'

# 筛选出矩阵中关注的剪接事件和RBP
filtered_events = up_events | dn_events
filtered_rbps = up_rbps | dn_rbps

# 保留矩阵中出现在filtered_events中的行，filtered_rbps中的列
filtered_matrix = ExS[ExS['name'].isin(filtered_events)]
filtered_matrix = filtered_matrix[['name'] + [rbp for rbp in filtered_matrix.columns[1:] if rbp in filtered_rbps]]

# 修改矩阵值
corrected_matrix = filtered_matrix.copy()

# Vectorized approach to modify values in the matrix based on conditions
for event in filtered_matrix['name']:
    # Mask for up-regulated events and RBPs
    up_mask = filtered_matrix['name'] == event
    for rbp in filtered_matrix.columns[1:]:
        if event in up_events and rbp in up_rbps:
            corrected_matrix.loc[up_mask, rbp] = filtered_matrix.loc[up_mask, rbp].where(filtered_matrix[rbp] == 1, 0)
        # Mask for down-regulated events and RBPs
        elif event in dn_events and rbp in dn_rbps:
            corrected_matrix.loc[up_mask, rbp] = filtered_matrix.loc[up_mask, rbp].where(filtered_matrix[rbp] == 1, 0)
        elif event in up_events and rbp in dn_rbps:
            corrected_matrix.loc[up_mask, rbp] = filtered_matrix.loc[up_mask, rbp].where(filtered_matrix[rbp] == -1, 0)
        elif event in dn_events and rbp in up_rbps:
            corrected_matrix.loc[up_mask, rbp] = filtered_matrix.loc[up_mask, rbp].where(filtered_matrix[rbp] == -1, 0)

# 保存修改后的矩阵
filtered_matrix.to_csv('DExDS.txt', index=False, sep='\t')
corrected_matrix.to_csv('correctDExDS.txt', index=False, sep='\t')

print("Matrices have been successfully saved.")
