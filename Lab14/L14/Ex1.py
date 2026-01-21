import math

# --- 1. Datele de intrare ---
s1_island = "ATCGATTCGATATCATACACGTAT"
s2_non_island = "CTCGACTAGTATGAAGTCCACGCTTG"
s_test = "CAGGTTGGAAACGTAA"
bases = ['A', 'C', 'G', 'T']

# --- Funcție: Calculare Probabilități de Tranziție ---
def get_transition_probabilities(sequence, bases):
    # Inițializare cu 1 (Laplace Smoothing)
    counts = {b1: {b2: 1.0 for b2 in bases} for b1 in bases}
    
    # Numărare tranziții
    for i in range(len(sequence) - 1):
        current_base = sequence[i]
        next_base = sequence[i+1]
        counts[current_base][next_base] += 1
        
    # Calcul probabilități
    probs = {b1: {b2: 0.0 for b2 in bases} for b1 in bases}
    for b1 in bases:
        total_transitions = sum(counts[b1].values())
        for b2 in bases:
            probs[b1][b2] = counts[b1][b2] / total_transitions
            
    return probs

# --- Antrenare Modele ---
prob_plus = get_transition_probabilities(s1_island, bases)
prob_minus = get_transition_probabilities(s2_non_island, bases)

# --- AFISARE TRANSITION MATRICES (NOU) ---

# Funcție helper pentru printare tabele
def print_matrix(matrix, title):
    print(f"\n--- {title} ---")
    print(f"{'':<5} {'A':<8} {'C':<8} {'G':<8} {'T':<8}")
    for b1 in bases:
        row_str = f"{b1:<5}"
        for b2 in bases:
            row_str += f"{matrix[b1][b2]:<8.3f} "
        print(row_str)

# Printam cele doua matrici de probabilitati
print_matrix(prob_plus, "Transition Probabilities: CpG Island (+)")
print_matrix(prob_minus, "Transition Probabilities: Non-Island (-)")

# --- Calculare Matrice Log-Likelihood ---
log_likelihood_matrix = {b1: {b2: 0.0 for b2 in bases} for b1 in bases}

print("\n--- Log-Likelihood Matrix (Beta) ---")
print(f"{'':<5} {'A':<8} {'C':<8} {'G':<8} {'T':<8}")

for b1 in bases:
    row_str = f"{b1:<5}"
    for b2 in bases:
        p_plus = prob_plus[b1][b2]
        p_minus = prob_minus[b1][b2]
        
        # Formula: log2( P+(t|s) / P-(t|s) )
        score = math.log2(p_plus / p_minus)
        log_likelihood_matrix[b1][b2] = score
        
        row_str += f"{score:<8.3f} "
    print(row_str)

# --- Scorare Secvență S ---
total_score = 0
print("\n--- Scoring Sequence S ---")
print(f"Sequence: {s_test}")

for i in range(len(s_test) - 1):
    curr = s_test[i]
    nxt = s_test[i+1]
    step_score = log_likelihood_matrix[curr][nxt]
    total_score += step_score

print(f"Total Log-Likelihood Score: {total_score:.4f}")

if total_score > 0:
    print("Result: The sequence is likely a CpG ISLAND (+).")
else:
    print("Result: The sequence is likely NON-ISLAND (-).")