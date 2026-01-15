import json
import numpy as np
import random

def load_dna_matrix(filename="dna_matrix.json"):
    with open(filename, 'r') as f:
        matrix = json.load(f)
    return np.array(matrix)

def synthesize_dna(matrix, length=50, start_state_index=None):
    states = ['A', 'C', 'G', 'T']
    n_states = len(states)
    
    # If no start state is given, pick a random one
    if start_state_index is None:
        current_idx = random.randint(0, n_states - 1)
    else:
        current_idx = start_state_index
        
    sequence = [states[current_idx]]
    
    print(f"Start State: {states[current_idx]}")

    for _ in range(length - 1):
        # Extract the probability column for the current state
        # (Column 'current_idx' contains probabilities of going TO other states)
        probs = matrix[:, current_idx]
        
        # Check if the column sums to ~1 (valid stochastic matrix)
        if np.sum(probs) == 0:
            print("Warning: Absorbing state reached (no outgoing transitions).")
            break
            
        # Normalize to ensure sum is exactly 1.0 for np.random.choice
        probs = probs / np.sum(probs)
        
        # Sample the next state based on these probabilities
        next_idx = np.random.choice(range(n_states), p=probs)
        
        sequence.append(states[next_idx])
        current_idx = next_idx
        
    return "".join(sequence)

if __name__ == "__main__":
    try:
        matrix = load_dna_matrix()
        
        print("--- Synthesizing New DNA Sequence ---")
        new_dna = synthesize_dna(matrix, length=50)
        print("\nGenerated Sequence:")
        print(new_dna)
        
    except FileNotFoundError:
        print("Error: 'dna_matrix.json' not found. Please run the Exercise 1 GUI first.")