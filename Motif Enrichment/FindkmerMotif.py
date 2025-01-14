from Bio import SeqIO
from collections import defaultdict

# 设置motif长度范围
min_len, max_len = 7, 7

# 初始化字典来统计子序列频率
motif_counts = defaultdict(int)

# 读取fasta文件并将序列转换为大写
fasta_file = "upstream.fasta"  # 或为exon.fasta，downstream.fasta
for record in SeqIO.parse(fasta_file, "fasta"):
    sequence = str(record.seq).upper()  # 转换为大写
    # 遍历序列中每种长度的所有子序列
    for length in range(min_len, max_len + 1):
        for i in range(len(sequence) - length + 1):
            motif = sequence[i:i+length]
            motif_counts[motif] += 1

# 将motif按频率排序
sorted_motifs = sorted(motif_counts.items(), key=lambda x: x[1], reverse=True)

# 将结果保存到txt文件
output_file = "motif_upstream_results.txt" # motif_exon_results.txt，或为motif_downstream_results.txt
with open(output_file, "w") as f:
    f.write("Motif\tCount\n")
    for motif, count in sorted_motifs:
        f.write(f"{motif}\t{count}\n")

print(f"Motif counts saved to {output_file}")
