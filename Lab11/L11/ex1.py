import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches


def needleman_wunsch(seq1, seq2, match=1, mismatch=-1, gap=0):
    n = len(seq1) 
    m = len(seq2) 
    
    score_matrix = np.zeros((m + 1, n + 1))
    
    for i in range(m + 1): score_matrix[i][0] = i * gap
    for j in range(n + 1): score_matrix[0][j] = j * gap
    
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if seq2[i-1] == seq1[j-1]:
                diagonal = score_matrix[i-1][j-1] + match
            else:
                diagonal = score_matrix[i-1][j-1] + mismatch
            
            up = score_matrix[i-1][j] + gap
            left = score_matrix[i][j-1] + gap
            
            score_matrix[i][j] = max(diagonal, up, left)
            
    align1 = ""
    align2 = ""
    matches_count = 0
    path_coords = [] 
    
    i, j = m, n
    path_coords.append((i, j))
    
    while i > 0 or j > 0:
        score = score_matrix[i][j]
        
        diag_val = score_matrix[i-1][j-1] if (i>0 and j>0) else -999
        up_val   = score_matrix[i-1][j]   if i>0 else -999
        left_val = score_matrix[i][j-1]   if j>0 else -999
        
        is_match = (seq2[i-1] == seq1[j-1]) if (i>0 and j>0) else False
        target_diag = diag_val + (match if is_match else mismatch)
        
        if i > 0 and j > 0 and score == target_diag:
            align1 = seq1[j-1] + align1
            align2 = seq2[i-1] + align2
            if is_match: matches_count += 1
            i -= 1; j -= 1
        elif i > 0 and score == up_val + gap:
            align1 = "-" + align1
            align2 = seq2[i-1] + align2
            i -= 1
        else:
            align1 = seq1[j-1] + align1
            align2 = "-" + align2
            j -= 1
            
        path_coords.append((i, j))

    return align1, align2, matches_count, score_matrix, path_coords


def print_console_results(s1, s2, as1, as2, matches):
    match_line = ""
    for c1, c2 in zip(as1, as2):
        match_line += "|" if (c1 == c2 and c1 != '-') else " "

    total_len = len(as1)
    similarity = int((matches / total_len) * 100)

    print("\nShow Alignment:")
    print("-" * 25)
    print(as1)
    print(match_line)
    print(as2)
    print("")
    print(f"Maches = {matches}")
    print(f"Lenght = {total_len}")
    print(f"Similarity = {similarity} %")
    print(f"Tracing back: M[{len(s1)},{len(s2)}]")
    print("-" * 25)


def show_figures(matrix, path, s1, s2):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle(f"Needleman-Wunsch Alignment (Gap=0, Match=1, Mismatch=-1)", fontsize=14)

    cax = ax1.imshow(matrix, cmap='magma', aspect='auto')
    ax1.set_title("Graphic representation of the alignment matrix")
    ax1.set_xlabel("Sequence 1 (Cols)")
    ax1.set_ylabel("Sequence 2 (Rows)")
    fig.colorbar(cax, ax=ax1, fraction=0.046, pad=0.04)

    rows, cols = matrix.shape
    ax2.set_title("Traceback path deviation")
    ax2.set_xlim(0, cols)
    ax2.set_ylim(0, rows)
    
    ax2.invert_yaxis()
    
    ax2.set_xticks(np.arange(cols) + 0.5)
    ax2.set_yticks(np.arange(rows) + 0.5)
    ax2.set_xticklabels(['-'] + list(s1))
    ax2.set_yticklabels(['-'] + list(s2))
    ax2.xaxis.tick_top() 
    ax2.set_facecolor('#FFFFE0')

    for x in range(cols + 1):
        ax2.axvline(x, color='black', linewidth=1)
    for y in range(rows + 1):
        ax2.axhline(y, color='black', linewidth=1)

    for r, c in path:
        rect = patches.Rectangle((c, r), 1, 1, linewidth=1, edgecolor='black', facecolor='#D32F2F')
        ax2.add_patch(rect)

    plt.tight_layout()
    plt.show()


S1 = "ACCGTGAAGCCAATAC"
S2 = "AGCGTGCAGCCAATAC"

aligned_s1, aligned_s2, match_count, mat, path = needleman_wunsch(S1, S2)

print_console_results(S1, S2, aligned_s1, aligned_s2, match_count)

print("Opening graphic visualization...")
show_figures(mat, path, S1, S2)