import pandas as pd

# 读取相关文件
ExS = pd.read_csv('ExS.txt', sep = '\t') 
up_events = set(pd.read_csv('up_events.txt', header=None)[0])  
dn_events = set(pd.read_csv('up_events.txt', header=None)[0])  
up_rbps = set(pd.read_csv('up_RBP.txt', header=None)[0]) 
dn_rbps = set(pd.read_csv('up_RBP.txt', header=None)[0]) 

# 筛选出矩阵中关注的剪接事件和RBP
filtered_events = up_events | dn_events
filtered_rbps = up_rbps | dn_rbps

# 保留矩阵中出现在filtered_events中的行，filtered_rbps中的列
filtered_matrix = rbp_matrix[rbp_matrix['name'].isin(filtered_events)]
filtered_matrix = filtered_matrix[['name'] + [rbp for rbp in filtered_matrix.columns[1:] if rbp in filtered_rbps]]

# 修改矩阵值
corrected_matrix = filtered_matrix.copy()
for event in filtered_matrix['name']:
    for rbp in filtered_matrix.columns[1:]:
        value = filtered_matrix.loc[filtered_matrix['name'] == event, rbp].values[0]
        if event in up_events and rbp in up_rbps:
            corrected_matrix.loc[filtered_matrix['name'] == event, rbp] = value if value == 1 else 0
        elif event in dn_events and rbp in dn_rbps:
            corrected_matrix.loc[filtered_matrix['name'] == event, rbp] = value if value == 1 else 0
        elif event in up_events and rbp in dn_rbps:
            corrected_matrix.loc[filtered_matrix['name'] == event, rbp] = value if value == -1 else 0
        elif event in dn_events and rbp in up_rbps:
            corrected_matrix.loc[filtered_matrix['name'] == event, rbp] = value if value == -1 else 0

# 保存修改后的矩阵
filtered_matrix.to_csv('DExDS.txt', index=False, sep = '\t')
corrected_matrix.to_csv('correctDExDS.txt', index=False, sep = '\t')
