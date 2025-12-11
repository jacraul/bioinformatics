import sys

# ==========================================
# 1. READ DATA
# ==========================================
def read_fasta_robust(filename):
    sequence = []
    try:
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith(">") or line.startswith("["):
                    continue
                sequence.append(line.upper())
    except FileNotFoundError:
        print(f"Error: File {filename} is missing.")
        return ""
    return "".join(sequence)

# ==========================================
# 2. THE 3 SCORING EQUATIONS
# ==========================================

def score_1_identity(seq1, seq2):
    """
    EQUATION 1: Simple Identity
    Calculates the percentage of identical characters at the same position.
    """
    length = min(len(seq1), len(seq2))
    if length == 0: return 0.0, 0
    
    matches = 0
    for i in range(length):
        if seq1[i] == seq2[i]:
            matches += 1
            
    percentage = (matches / length) * 100
    return percentage, matches

def score_2_weighted(seq1, seq2, match_reward=1, mismatch_penalty=1):
    """
    EQUATION 2: Weighted Score
    Awards points for matches and DEDUCTS points for mismatches.
    The result is then normalized to a 0-100% scale.
    """
    length = min(len(seq1), len(seq2))
    if length == 0: return 0.0, 0
    
    raw_score = 0
    # Calculate theoretical max and min scores for normalization
    max_possible_score = length * match_reward
    min_possible_score = length * -mismatch_penalty
    
    for i in range(length):
        if seq1[i] == seq2[i]:
            raw_score += match_reward
        else:
            raw_score -= mismatch_penalty
            
    # Normalization (Min-Max Scaling) to bring score between 0 and 100
    # Formula: (val - min) / (max - min) * 100
    normalized_score = (raw_score - min_possible_score) / (max_possible_score - min_possible_score) * 100
    
    return normalized_score, raw_score

def score_3_jaccard_kmers(seq1, seq2, k=3):
    """
    EQUATION 3: Jaccard Similarity (on K-mers)
    Breaks sequences into 'words' of k letters (trigrams) and checks overlap.
    Useful for detecting structural similarity even if there are small shifts.
    """
    if len(seq1) < k or len(seq2) < k: return 0.0, 0, 0

    # Generate sets of k-mers (e.g., "ATCG" -> {"ATC", "TCG"})
    set1 = set(seq1[i : i+k] for i in range(len(seq1) - k + 1))
    set2 = set(seq2[i : i+k] for i in range(len(seq2) - k + 1))
    
    intersection = len(set1.intersection(set2))
    union = len(set1.union(set2))
    
    if union == 0: return 0.0, 0, 0
    
    jaccard_index = (intersection / union) * 100
    return jaccard_index, intersection, union

# ==========================================
# 3. FIND BEST WINDOW (Intermediate Layer)
# ==========================================
def find_best_alignment_window(s_flu, s_cov, window_size=50):
    """
    Scans the COVID genome with a chunk of Influenza to find 
    a region with relevant similarity to apply the calculations on.
    """
    print(f"Searching for best local alignment (Window: {window_size}bp)...")
    
    best_score = -1
    best_flu_seq = ""
    best_cov_seq = ""
    best_positions = (0, 0) # (flu_start, cov_start)

    # Optimization: scan Influenza every 100bp
    for i in range(0, len(s_flu) - window_size, 100):
        chunk_flu = s_flu[i : i+window_size]
        
        # Scan Covid every 500bp
        for j in range(0, len(s_cov) - window_size, 500):
            chunk_cov = s_cov[j : j+window_size]
            
            # Use simple score for scanning
            perc, _ = score_1_identity(chunk_flu, chunk_cov)
            
            if perc > best_score:
                best_score = perc
                best_flu_seq = chunk_flu
                best_cov_seq = chunk_cov
                best_positions = (i, j)
                
    return best_flu_seq, best_cov_seq, best_positions


if __name__ == "__main__":
    f_flu = 'Influenza.fasta'
    f_cov = 'Covid19.fasta'

    seq_flu = read_fasta_robust(f_flu)
    seq_cov = read_fasta_robust(f_cov)

    if seq_flu and seq_cov:
        w_flu, w_cov, pos = find_best_alignment_window(seq_flu, seq_cov, window_size=60)
        
        print("\n" + "="*60)
        print(f"CALCULATION RESULTS FOR BEST LOCAL ALIGNMENT")
        print("="*60)
        print(f"Influenza Position: {pos[0]} | COVID Position: {pos[1]}")
        print("-" * 60)
        print(f"Sequence 1 (Flu):   {w_flu}")
        print(f"Sequence 2 (Covid): {w_cov}")
        print("-" * 60)
        
        res1, matches = score_1_identity(w_flu, w_cov)
        print(f"\n1. EQUATION: IDENTITY SCORE")
        print(f"   Formula: (Matches / Total_Len) * 100")
        print(f"   Calc:    ({matches} / {len(w_flu)}) * 100")
        print(f"   RESULT:  {res1:.2f}%")
        
        res2, raw = score_2_weighted(w_flu, w_cov, match_reward=1, mismatch_penalty=1)
        print(f"\n2. EQUATION: WEIGHTED SCORE")
        print(f"   Rule:    Match = +1 | Mismatch = -1")
        print(f"   Raw Score: {raw} (out of max possible {len(w_flu)})")
        print(f"   Normalization: (Raw - Min) / (Max - Min) * 100")
        print(f"   RESULT:  {res2:.2f}%")
        
        res3, inter, union = score_3_jaccard_kmers(w_flu, w_cov, k=3)
        print(f"\n3. EQUATION: JACCARD SIMILARITY (K-mers, k=3)")
        print(f"   Formula: |A intersect B| / |A union B| * 100")
        print(f"   Common Trigrams: {inter}")
        print(f"   Total Unique Trigrams: {union}")
        print(f"   RESULT:  {res3:.2f}%")
        print("="*60)
        
    else:
        print("Could not read input files.")