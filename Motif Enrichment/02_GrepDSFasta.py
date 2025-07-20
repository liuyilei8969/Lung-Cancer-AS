from Bio import SeqIO

def read_fasta(fasta_file):
    """
    Reads a FASTA file and returns a dictionary of chromosome sequences.
    """
    fasta_dict = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        fasta_dict[record.id] = record.seq
    return fasta_dict

def read_bed(bed_file):
    """
    Reads a BED file and returns a list containing id and coordinate information.
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
    Extracts sequences based on coordinates in the BED file from the FASTA file 
    and outputs them to three new FASTA files.
    """
    with open(exon_file, "w") as exon_out, open(upstream_file, "w") as upstream_out, open(downstream_file, "w") as downstream_out:
        for chr_name, start, end, seq_id, strand in bed_entries:
            if chr_name in fasta_dict:
                exon_seq = fasta_dict[chr_name][start:end]

                if strand == "+":
                    # For the positive strand, extract upstream and downstream sequences
                    upstream_seq = fasta_dict[chr_name][start-200:start]
                    downstream_seq = fasta_dict[chr_name][end:end+200]
                elif strand == "-":
                    # For the negative strand, extract upstream and downstream sequences
                    upstream_seq = fasta_dict[chr_name][end:end+200]      
                    downstream_seq = fasta_dict[chr_name][start-200:start]  

                    # Reverse complement the sequences
                    exon_seq = exon_seq.reverse_complement()
                    upstream_seq = upstream_seq.reverse_complement()
                    downstream_seq = downstream_seq.reverse_complement()

                # Format and write sequences to FASTA files
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
    fasta_file = "hg19.fa"  # Input FASTA file
    bed_file = "DSposition.txt"    # Input BED file
    exon_file = "exon.fasta"     # Output exon FASTA file
    upstream_file = "upstream.fasta"  # Output upstream sequence FASTA file
    downstream_file = "downstream.fasta"  # Output downstream sequence FASTA file

    fasta_dict = read_fasta(fasta_file)      # Read the FASTA file
    bed_entries = read_bed(bed_file)         # Read the BED file
    extract_sequences(fasta_dict, bed_entries, exon_file, upstream_file, downstream_file)  # Extract and output sequences
