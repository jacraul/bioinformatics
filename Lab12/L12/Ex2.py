import math
import matplotlib.pyplot as plt
import os


motifs = [
    "GAGGTAAAC",
    "TCCGTAAGT",
    "CAGGTTGGA",
    "ACAGTCAGT",
    "TAGGTCATT",
    "TAGGTACTG",
    "ATGGTAACT",
    "CAGGTATAC",
    "TGTGTGAGT",
    "AAGGTAAGT"
]

file_names = [
    "H5N1_Goose_Guangdong_1996.fasta",
    "H7N9_Anhui_2013.fasta",
    "B_Brisbane_2008.fasta",
    "H1N1_PR8_1934.fasta",
    "H1N1_SpanishFlu_1918.fasta",
    "B_Yamagata_1988.fasta",
    "B_Lee_1940.fasta",
    "B_Victoria_1987.fasta",
    "H1N1_California_2009.fasta",
    "H3N2_NewYork_2004.fasta"
]

motif_length = 9
num_sequences = len(motifs)
bases = ['A', 'C', 'G', 'T']
background_freq = 0.25


count_matrix = {b: [0] * motif_length for b in bases}
prob_matrix = {b: [0.0] * motif_length for b in bases}
log_likelihood_matrix = {b: [0.0] * motif_length for b in bases}

pseudocount = 1
for b in bases:
    for i in range(motif_length):
        count_matrix[b][i] = pseudocount

for seq in motifs:
    for i, char in enumerate(seq):
        if char in count_matrix:
            count_matrix[char][i] += 1


denom = num_sequences + (4 * pseudocount)

for b in bases:
    for i in range(motif_length):
        prob_matrix[b][i] = count_matrix[b][i] / denom

for b in bases:
    for i in range(motif_length):
        p_obs = prob_matrix[b][i]
        score = math.log(p_obs / background_freq)
        log_likelihood_matrix[b][i] = score


def parse_fasta(filename):
    """
    Reads a FASTA file and concatenates all sequences into one string.
    Returns the concatenated sequence and a short ID.
    """
    sequence_data = ""
    try:
        with open(filename, 'r') as f:
            lines = f.readlines()
            
        for line in lines:
            line = line.strip()
            if not line: continue
            if line.startswith(">"):
                continue # Skip headers, we just want the raw sequence stream
            else:
                # Filter to ensure only valid bases ACGT are processed
                clean_line = "".join([c for c in line.upper() if c in bases])
                sequence_data += clean_line
                
        return sequence_data
    except FileNotFoundError:
        print(f"Warning: File {filename} not found.")
        return None

def scan_genome(sequence, matrix):
    """
    Scans the sequence with the PWM matrix.
    Returns a list of scores corresponding to each position.
    """
    scores = []
    # Sliding window
    for i in range(len(sequence) - motif_length + 1):
        window = sequence[i : i + motif_length]
        current_score = 0
        
        for pos, char in enumerate(window):
            if char in matrix:
                current_score += matrix[char][pos]
            else:
                current_score += -99.0 # Heavy penalty for non-base characters
        
        scores.append(current_score)
    return scores

def plot_genome_chart(filename, scores):
    """
    Creates a chart for the genome showing signal strength.
    """
    # Create figure
    plt.figure(figsize=(12, 5))
    
    # Plot data
    plt.plot(scores, color='#1f77b4', linewidth=0.6, label='Motif Score')
    
    # Add threshold line (Score > 0 implies likely match)
    plt.axhline(y=0, color='red', linestyle='--', linewidth=1, label='Signal Threshold (0)')
    
    # Styling
    plt.title(f"Genome Scan: {filename}", fontsize=14)
    plt.xlabel("Genomic Position (bp)", fontsize=12)
    plt.ylabel("Log-Likelihood Score", fontsize=12)
    plt.legend(loc="upper right")
    plt.grid(True, linestyle=':', alpha=0.6)
    
    # Show plot
    plt.tight_layout()
    plt.show()

# --- 4. Main Execution Loop ---

print(f"{'File Name':<35} | {'Genome Length':<15} | {'Max Score':<10}")
print("-" * 70)

for f_name in file_names:
    genome_seq = parse_fasta(f_name)
    
    if genome_seq:
        # Scan the genome
        scores = scan_genome(genome_seq, log_likelihood_matrix)
        
        # Print summary statistics
        max_score = max(scores) if scores else 0
        print(f"{f_name:<35} | {len(genome_seq):<15} | {max_score:.2f}")
        
        # Generate the chart
        plot_genome_chart(f_name, scores)

print("\nProcessing Complete.")