import os
import pandas as pd
from collections import defaultdict

# Define folder path
folder_path = "/path/to/ENCODE_DS/result"

# Initialize data structures to store RBP and target data
rbp_target_map = defaultdict(lambda: defaultdict(list))
rbp_cellline_count = defaultdict(set)  # Track cell lines each RBP appears in

# Loop through files in the folder
for file_name in os.listdir(folder_path):
    if file_name.endswith("-01.dn.txt") or file_name.endswith("-01.up.txt"):
        # Extract cell line, RBP name, and direction
        parts = file_name.split(".")
        cell_line = parts[0]
        rbp_name = parts[1]
        direction = "dn" if "dn.txt" in file_name else "up"

        # Track the cell lines for each RBP
        rbp_cellline_count[rbp_name].add(cell_line)

        # Read the file and extract the target column (5th column, zero-indexed)
        file_path = os.path.join(folder_path, file_name)
        df = pd.read_csv(file_path, sep="\t", header=None)
        targets = df.iloc[:, 4].tolist()

        # Add targets to the map for the RBP
        for target in targets:
            rbp_target_map[target][rbp_name].append(direction)

# Initialize the result matrix for all targets and RBPs
all_targets = sorted(rbp_target_map.keys())
all_rbps = sorted({rbp for target in rbp_target_map for rbp in rbp_target_map[target]})

# Create the result matrix with default value 0
result_matrix = pd.DataFrame(0, index=all_targets, columns=all_rbps)

# Fill the result matrix based on the RBP-target mapping
for target, rbps in rbp_target_map.items():
    for rbp, directions in rbps.items():
        if len(rbp_cellline_count[rbp]) == 1:  # If the RBP appears in only one cell line
            result_matrix.loc[target, rbp] = 1 if directions[0] == 'up' else -1
        elif len(set(directions)) == 1:  # If the RBP appears in multiple cell lines but with consistent regulation
            result_matrix.loc[target, rbp] = 1 if directions[0] == 'up' else -1

# Load the BED file and extract the splicing event names
bed_file = pd.read_csv('../data/Hs.seq.all.cass.chrom.can.exon.bed', sep='\t', header=None)

# Extract unique splicing event names (column index 3 for BED format)
bed_events = bed_file[3].unique()

# Filter the result matrix to keep only the events present in the BED file
filtered_matrix = result_matrix[result_matrix.index.isin(bed_events)]

# Save the result matrix to a tab-separated file
output_path = "ExS.txt"
filtered_matrix.to_csv(output_path, sep="\t")

print(f"Filtered result matrix saved to {output_path}")
