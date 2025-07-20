import scipy.stats as stats
import math
import pandas as pd

def read_data(fg_file, bg_file):
    """
    Reads the foreground and background files and returns dataframes.
    """
    fg_data = pd.read_csv(fg_file, sep="\t")
    bg_data = pd.read_csv(bg_file, sep="\t")
    return fg_data, bg_data

def calculate_enrichment(fg_data, bg_data):
    """
    Performs enrichment analysis.
    Parameters:
        fg_data: Foreground data containing Motif and Count columns.
        bg_data: Background data containing Motif and Count columns.
        test_method: Statistical method used (binom or hyperg).
    Returns:
        A dataframe containing Motif, fg_count, bg_count, log2FC, pval, zscore columns.
    """
    # Merge foreground and background data
    merged_data = pd.merge(fg_data, bg_data, on="RBP", how="outer", suffixes=("_fg", "_bg")).fillna(0)
    fg_sum = merged_data["Count_fg"].sum()
    bg_sum = merged_data["Count_bg"].sum()

    results = []

    for _, row in merged_data.iterrows():
        RBP = row["RBP"]
        fg_count = row["Count_fg"]
        bg_count = row["Count_bg"]

        # Calculate log2FC (log fold change)
        log2FC = (fg_sum == 0 or bg_count - fg_count == 0 or fg_count == 0 or bg_sum - fg_sum == 0)
        log2FC = math.log2(fg_count * (bg_sum - fg_sum) / (fg_sum * (bg_count - fg_count))) if not log2FC else 'NA'

        # Calculate z-score
        r = fg_sum / bg_sum
        zscore = (fg_count - bg_count * r) / math.sqrt(bg_count * r * (1 - r))

        # Calculate p-value using hypergeometric distribution
        pval_hyperg = stats.hypergeom.sf(fg_count - 1, bg_sum, bg_count, fg_sum)

        results.append({
            "RBP": RBP,
            "fg_count": fg_count,
            "bg_count": bg_count,
            "log2FC": log2FC,
            "zscore": zscore,
            "hyperg": pval_hyperg
        })

    return pd.DataFrame(results)

def save_results(results, output_file):
    """
    Saves the results to a file.
    """
    results.to_csv(output_file, sep="\t", index=False)

if __name__ == "__main__":
    fg_file = "DS_RBP_exon_results.txt"  # Foreground file path (e.g., DS_RBP_downstream_results.txt, DS_RBP_upstream_results.txt)
    bg_file = "All_RBP_exon_results.txt"  # Background file path (e.g., All_RBP_downstream_results.txt, All_RBP_upstream_results.txt)
    output_file = "DS_exon_enrich_RBP.txt"  # Output file (e.g., DS_downstream_enrich_RBP.txt, DS_upstream_enrich_RBP.txt)

    # Read data
    fg_data, bg_data = read_data(fg_file, bg_file)

    # Perform enrichment analysis
    results = calculate_enrichment(fg_data, bg_data)

    # Save the results
    save_results(results, output_file)

    print(f"Enrichment analysis results saved to {output_file}")
