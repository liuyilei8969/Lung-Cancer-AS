import pandas as pd
from scipy.stats import hypergeom
import numpy as np
from statsmodels.stats.multitest import multipletests

# Read data
foreground = pd.read_csv('DE_DS_count.txt', sep='\t', index_col=0)
background = pd.read_csv('E_S_count.txt', sep='\t', index_col=0)
rbps = foreground.index.intersection(background.index)

# Set total counts
total_foreground = 127  # For Subtype-Specific, replace with 178
total_background = 10931  # Total number of rows in ExS matrix

# Perform hypergeometric test
def hypergeometric_test(fg_count, bg_count, total_fg, total_bg):
    M = total_fg + total_bg  # Total population size
    n = fg_count + bg_count  # Number of successes
    N = total_fg  # Size of the target population
    x = fg_count  # Number of successful events
    p_value = hypergeom.sf(x - 1, M, N, n)  # Calculate upper tail probability
    return p_value

# Perform fold change analysis
def fold_change(fg_count, bg_count):
    return fg_count * total_background / ((bg_count if bg_count > 0 else np.nan) * total_foreground)  # Avoid division by zero

# Perform analysis for each RBP
results = []

for rbp in rbps:
    fg_count = foreground.loc[rbp, 'total']
    bg_count = background.loc[rbp, 'total']

    # Perform hypergeometric test
    p_value_hypergeometric = hypergeometric_test(fg_count, bg_count, total_foreground, total_background)
    
    # Perform fold change analysis
    fc = fold_change(fg_count, bg_count)

    # Store results
    results.append({
        'RBP': rbp,
        'DScount': fg_count,
        'count': bg_count,
        'hypergeo_p_value': p_value_hypergeometric,
        'fold_change': fc
    })

# Convert results to DataFrame
results_df = pd.DataFrame(results)

# Multiple testing correction: Use Benjamini-Hochberg method
pvals = results_df['hypergeo_p_value']
_, corrected_pvals, _, _ = multipletests(pvals, method='fdr_bh')

# Add corrected p-values to results
results_df['corrected_p_value'] = corrected_pvals

# Optionally, save results to a file
results_df.to_csv('RBPactivity_corrected.txt', index=False, sep='\t')

print("Analysis complete, results saved to 'RBPactivity_corrected.txt'.")
