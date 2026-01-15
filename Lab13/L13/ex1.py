import tkinter as tk
from tkinter import messagebox
import numpy as np
import json
import random


class StatePredictor:
    def __init__(self, transition_matrix, initial_vector):
        self.matrix = np.array(transition_matrix)
        self.current_state = np.array(initial_vector)

    def predict(self, steps=5):
        predictions = []
        for _ in range(steps):
            next_state = np.dot(self.matrix, self.current_state)
            predictions.append(next_state)
            self.current_state = next_state
        return predictions

class SequenceSynthesizer:
    def __init__(self, transition_matrix, states):
        self.matrix = np.array(transition_matrix)
        self.states = states
        self.n_states = len(states)

    def generate(self, length=20):
        current_idx = random.randint(0, self.n_states - 1)
        sequence = [self.states[current_idx]]
        
        for _ in range(length - 1):
            probs = self.matrix[:, current_idx]
            
            if np.sum(probs) == 0:
                break
                
            probs = probs / np.sum(probs)
            
            next_idx = np.random.choice(range(self.n_states), p=probs)
            
            sequence.append(self.states[next_idx])
            current_idx = next_idx
            
        return "".join(sequence)


class DnaTool:
    def __init__(self, root):
        self.root = root
        self.root.title("Ex1: DNA Analyzer & Synthesizer")
        self.root.geometry("650x700")
        self.matrix = None
        self.states = ['A', 'C', 'G', 'T']

        lbl_frame = tk.LabelFrame(root, text="1. Input Data", padx=10, pady=10)
        lbl_frame.pack(fill="x", padx=10, pady=5)
        
        tk.Button(lbl_frame, text="Generate Random DNA (50 chars)", command=self.generate_random_dna).pack(anchor="w")
        self.txt_seq = tk.Text(lbl_frame, height=2, width=70)
        self.txt_seq.pack(pady=5)

        mat_frame = tk.LabelFrame(root, text="2. Transition Matrix (JSON)", padx=10, pady=10)
        mat_frame.pack(fill="x", padx=10, pady=5)
        
        tk.Button(mat_frame, text="Compute & Save Matrix", command=self.compute_matrix).pack(anchor="w")
        self.txt_matrix = tk.Text(mat_frame, height=5, width=70)
        self.txt_matrix.pack(pady=5)

        pred_frame = tk.LabelFrame(root, text="3. State Prediction (Distribution)", padx=10, pady=10)
        pred_frame.pack(fill="x", padx=10, pady=5)
        
        tk.Button(pred_frame, text="Predict Next 5 Steps (Start: 100% A)", command=self.run_prediction).pack(anchor="w")
        self.txt_pred = tk.Text(pred_frame, height=5, width=70)
        self.txt_pred.pack(pady=5)

        syn_frame = tk.LabelFrame(root, text="4. Synthesize New Sequence", padx=10, pady=10)
        syn_frame.pack(fill="both", expand=True, padx=10, pady=5)
        
        tk.Button(syn_frame, text="Generate New DNA from Matrix", command=self.run_synthesis).pack(anchor="w")
        self.txt_syn = tk.Text(syn_frame, height=3, width=70, fg="blue")
        self.txt_syn.pack(pady=5)

    def generate_random_dna(self):
        seq = "".join(random.choices("ACGT", k=50))
        self.txt_seq.delete(1.0, tk.END)
        self.txt_seq.insert(tk.END, seq)

    def compute_matrix(self):
        seq = self.txt_seq.get(1.0, tk.END).strip()
        if not seq: return
        
        state_map = {char: i for i, char in enumerate(self.states)}
        counts = np.zeros((4, 4))
        
        for i in range(len(seq) - 1):
            curr, next_ = seq[i], seq[i+1]
            c, r = state_map[curr], state_map[next_]
            counts[r][c] += 1
            
        self.matrix = np.zeros((4, 4))
        for col in range(4):
            total = np.sum(counts[:, col])
            if total > 0:
                self.matrix[:, col] = counts[:, col] / total
            else:
                self.matrix[col, col] = 1.0 

        with open("dna_matrix.json", "w") as f:
            json.dump(self.matrix.tolist(), f)
        
        self.txt_matrix.delete(1.0, tk.END)
        for row in self.matrix:
            self.txt_matrix.insert(tk.END, str([round(x, 2) for x in row]) + "\n")
        messagebox.showinfo("Success", "Matrix computed and saved to dna_matrix.json")

    def run_prediction(self):
        if self.matrix is None: return
        initial = [1.0, 0.0, 0.0, 0.0] 
        predictor = StatePredictor(self.matrix, initial)
        results = predictor.predict(5)
        
        self.txt_pred.delete(1.0, tk.END)
        for i, res in enumerate(results, 1):
            self.txt_pred.insert(tk.END, f"Step {i}: {[round(x, 3) for x in res]}\n")

    def run_synthesis(self):
        if self.matrix is None: return
        synth = SequenceSynthesizer(self.matrix, self.states)
        new_seq = synth.generate(length=50)
        self.txt_syn.delete(1.0, tk.END)
        self.txt_syn.insert(tk.END, new_seq)

if __name__ == "__main__":
    root = tk.Tk()
    DnaTool(root)
    root.mainloop()