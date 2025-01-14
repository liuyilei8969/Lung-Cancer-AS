from Bio import SeqIO

def read_fasta(fasta_file):
    """
    读取FASTA文件，并以字典形式返回染色体序列。
    """
    fasta_dict = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        fasta_dict[record.id] = record.seq
    return fasta_dict

def read_bed(bed_file):
    """
    读取BED文件，并返回包含id和坐标信息的列表。
    """
    bed_entries = []
    with open(bed_file, "r") as file:
        for line in file:
            fields = line.strip().split()
            chr_name, start, end, seq_id, _, strand = fields[:6]
            start = int(start)
            end = int(end)
            bed_entries.append((chr_name, start, end, seq_id, strand))
    return bed_entries

def extract_sequences(fasta_dict, bed_entries, exon_file, upstream_file, downstream_file):
    """
    根据BED文件中的坐标，从FASTA中提取对应的序列并输出到三个新的FASTA文件。
    """
    with open(exon_file, "w") as exon_out, open(upstream_file, "w") as upstream_out, open(downstream_file, "w") as downstream_out:
        for chr_name, start, end, seq_id, strand in bed_entries:
            if chr_name in fasta_dict:
                exon_seq = fasta_dict[chr_name][start:end]

                if strand == "+":
                    # 正链上下游
                    upstream_seq = fasta_dict[chr_name][start-200:start]
                    downstream_seq = fasta_dict[chr_name][end:end+200]
                elif strand == "-":
                    # 负链上下游
                    upstream_seq = fasta_dict[chr_name][end:end+200]      
                    downstream_seq = fasta_dict[chr_name][start-200:start]  

                    # 反向互补处理
                    exon_seq = exon_seq.reverse_complement()
                    upstream_seq = upstream_seq.reverse_complement()
                    downstream_seq = downstream_seq.reverse_complement()

                # 格式化输出FASTA
                exon_out.write(f">{seq_id}_exon\n")
                for i in range(0, len(exon_seq), 80):
                    exon_out.write(f"{exon_seq[i:i+80]}\n")

                upstream_out.write(f">{seq_id}_upstream\n")
                for i in range(0, len(upstream_seq), 80):
                    upstream_out.write(f"{upstream_seq[i:i+80]}\n")

                downstream_out.write(f">{seq_id}_downstream\n")
                for i in range(0, len(downstream_seq), 80):
                    downstream_out.write(f"{downstream_seq[i:i+80]}\n")
            else:
                print(f"Warning: {chr_name} not found in the FASTA file.")

if __name__ == "__main__":
    fasta_file = "hg19.fa"  # 输入的FASTA文件
    bed_file = "DSposition.txt"    # 输入的BED文件
    exon_file = "exon.fasta"     # 输出的外显子FASTA文件
    upstream_file = "upstream.fasta"  # 输出的上游序列FASTA文件
    downstream_file = "downstream.fasta"  # 输出的下游序列FASTA文件

    fasta_dict = read_fasta(fasta_file)      # 读取FASTA文件
    bed_entries = read_bed(bed_file)         # 读取BED文件
    extract_sequences(fasta_dict, bed_entries, exon_file, upstream_file, downstream_file)  # 提取并输出序列
