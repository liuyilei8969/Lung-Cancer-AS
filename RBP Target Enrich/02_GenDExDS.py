import pandas as pd

# Read relevant files
ExS = pd.read_csv('ExS.txt', sep='\t') 

# Correctly reading down-regulated files for events and RBPs
up_events = set(pd.read_csv('up_events.txt', header=None)[0])  
dn_events = set(pd.read_csv('dn_events.txt', header=None)[0])  # Fixed to read 'dn_events.txt'
up_rbps = set(pd.read_csv('up_RBP.txt', header=None)[0]) 
dn_rbps = set(pd.read_csv('dn_RBP.txt', header=None)[0])  # Fixed to read 'dn_RBP.txt'

# Filter the matrix to include only the events and RBPs of interest
filtered_events = up_events | dn_events
filtered_rbps = up_rbps | dn_rbps

# Retain rows from the matrix that are in filtered_events and columns that are in filtered_rbps
filtered_matrix = ExS[ExS['name'].isin(filtered_events)]
filtered_matrix = filtered_matrix[['name'] + [rbp for rbp in filtered_matrix.columns[1:] if rbp in filtered_rbps]]

# Modify matrix values
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

# Save the modified matrices
filtered_matrix.to_csv('DExDS.txt', index=False, sep='\t')
corrected_matrix.to_csv('correctDExDS.txt', index=False, sep='\t')

print("Matrices have been successfully saved.")
