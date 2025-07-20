from Bio import SeqIO
from collections import defaultdict

# Set the range for motif length
min_len, max_len = 7, 7

# Initialize a dictionary to count motif frequencies
motif_counts = defaultdict(int)

# Read the fasta file and convert sequences to uppercase
fasta_file = "upstream.fasta"  # Or exon.fasta, downstream.fasta
for record in SeqIO.parse(fasta_file, "fasta"):
    sequence = str(record.seq).upper()  # Convert to uppercase
    # Iterate through all substrings of each length in the sequence
    for length in range(min_len, max_len + 1):
        for i in range(len(sequence) - length + 1):
            motif = sequence[i:i+length]
            motif_counts[motif] += 1

# Sort the motifs by frequency
sorted_motifs = sorted(motif_counts.items(), key=lambda x: x[1], reverse=True)

# Save the results to a txt file
output_file = "motif_upstream_results.txt"  # motif_exon_results.txt, or motif_downstream_results.txt
with open(output_file, "w") as f:
    f.write("Motif\tCount\n")
    for motif, count in sorted_motifs:
        f.write(f"{motif}\t{count}\n")

print(f"Motif counts saved to {output_file}")
