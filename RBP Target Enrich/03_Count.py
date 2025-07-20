import pandas as pd

# Read the ExS matrix
rbp_matrix = pd.read_csv('DExDS.txt', sep='\t')  # Or use other ExS matrices

# Initialize an empty dictionary to store activation and repression counts for each RBP
activation_repression_count = {}

# Loop through each RBP column (skip the first column, which contains splicing event names)
for rbp in rbp_matrix.columns[1:]:  # Skipping the first column ('Splicing Event')
    # Count the number of activations (1) and repressions (-1) for each RBP
    repression_count = (rbp_matrix[rbp] == -1).sum()
    activation_count = (rbp_matrix[rbp] == 1).sum()

    # Store the counts in the dictionary
    activation_repression_count[rbp] = {'activation': activation_count, 'repression': repression_count}

# Convert the dictionary to a DataFrame
activation_repression_df = pd.DataFrame.from_dict(activation_repression_count, orient='index')

# Save the result to a tab-separated file
activation_repression_df.to_csv('DE_DS_count.txt', sep='\t')

print("Activation and repression counts have been saved to 'DE_DS_count.txt'")
