import tkinter as tk
from tkinter import messagebox
import numpy as np
import json
import random
import re

# --- Core Logic (Reused) ---

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

class TextSynthesizer:
    def __init__(self, transition_matrix, word_list):
        self.matrix = np.array(transition_matrix)
        self.words = word_list # List of words indexed 0..N
        self.n = len(word_list)

    def generate(self, length=15):
        curr_idx = random.randint(0, self.n - 1)
        output = [self.words[curr_idx]]
        
        for _ in range(length - 1):
            probs = self.matrix[:, curr_idx]
            if np.sum(probs) == 0: break
            probs = probs / np.sum(probs)
            
            curr_idx = np.random.choice(range(self.n), p=probs)
            output.append(self.words[curr_idx])
            
        return " ".join(output)


class TextTool:
    def __init__(self, root):
        self.root = root
        self.root.title("Ex2: English Text Markov Chain")
        self.root.geometry("700x750")
        self.matrix = None
        self.unique_words = []
        self.mapping = {} 

        lbl_frame = tk.LabelFrame(root, text="1. English Text Input (~300 chars)", padx=10, pady=10)
        lbl_frame.pack(fill="x", padx=10, pady=5)
        
        tk.Button(lbl_frame, text="Load Default Text", command=self.load_text).pack(anchor="w")
        self.txt_input = tk.Text(lbl_frame, height=4, width=80)
        self.txt_input.pack(pady=5)

        mat_frame = tk.LabelFrame(root, text="2. Transition Matrix & Symbols", padx=10, pady=10)
        mat_frame.pack(fill="x", padx=10, pady=5)
        
        tk.Button(mat_frame, text="Compute & Save JSON", command=self.compute).pack(anchor="w")
        self.txt_legend = tk.Text(mat_frame, height=3, width=80, bg="#f0f0f0")
        self.txt_legend.pack(pady=2)
        self.txt_matrix = tk.Text(mat_frame, height=5, width=80)
        self.txt_matrix.pack(pady=2)

        pred_frame = tk.LabelFrame(root, text="3. State Prediction (Probabilities)", padx=10, pady=10)
        pred_frame.pack(fill="x", padx=10, pady=5)
        
        tk.Button(pred_frame, text="Predict Distribution (5 Steps)", command=self.run_prediction).pack(anchor="w")
        self.txt_pred = tk.Text(pred_frame, height=5, width=80)
        self.txt_pred.pack(pady=5)

        syn_frame = tk.LabelFrame(root, text="4. Synthesize New Text", padx=10, pady=10)
        syn_frame.pack(fill="both", expand=True, padx=10, pady=5)
        
        tk.Button(syn_frame, text="Generate New Sentence", command=self.run_synthesis).pack(anchor="w")
        self.txt_syn = tk.Text(syn_frame, height=4, width=80, fg="blue", font=("Arial", 10, "bold"))
        self.txt_syn.pack(pady=5)

    def load_text(self):
        text = ("The cat sat on the mat. The dog sat on the log. "
                "The cat saw the dog. The dog saw the cat. "
                "They played on the mat. The mat was blue. "
                "The log was brown. The cat ran to the log. "
                "The dog ran to the mat. It was a good day for the cat and the dog.")
        self.txt_input.delete(1.0, tk.END)
        self.txt_input.insert(tk.END, text)

    def compute(self):
        raw = self.txt_input.get(1.0, tk.END).strip()
        if not raw: return
        
        clean = re.sub(r'[^\w\s]', '', raw).lower()
        words = clean.split()
        self.unique_words = sorted(list(set(words)))
        n = len(self.unique_words)
        
        self.mapping = {w: chr(65+i) for i, w in enumerate(self.unique_words)}
        
        counts = np.zeros((n, n))
        for i in range(len(words)-1):
            curr_i = self.unique_words.index(words[i])
            next_i = self.unique_words.index(words[i+1])
            counts[next_i][curr_i] += 1 
            
        self.matrix = np.zeros((n, n))
        for c in range(n):
            s = np.sum(counts[:, c])
            if s > 0: self.matrix[:, c] = counts[:, c] / s
            else: self.matrix[c, c] = 1.0

        data = {"mapping": self.mapping, "matrix": self.matrix.tolist()}
        with open("text_matrix.json", "w") as f:
            json.dump(data, f)

        self.txt_legend.delete(1.0, tk.END)
        self.txt_legend.insert(tk.END, f"Mapping: {self.mapping}")
        
        self.txt_matrix.delete(1.0, tk.END)
        for row in self.matrix:
            self.txt_matrix.insert(tk.END, " ".join([f"{x:.1f}" for x in row]) + "\n")

    def run_prediction(self):
        if self.matrix is None: return
        n = len(self.unique_words)
        v0 = np.zeros(n); v0[0] = 1.0 
        
        pred = StatePredictor(self.matrix, v0)
        results = pred.predict(5)
        
        self.txt_pred.delete(1.0, tk.END)
        self.txt_pred.insert(tk.END, f"Start: '{self.unique_words[0]}'\n")
        for i, vec in enumerate(results, 1):
            top_idx = np.argmax(vec)
            self.txt_pred.insert(tk.END, f"Step {i}: Most likely -> '{self.unique_words[top_idx]}' ({vec[top_idx]:.2f})\n")

    def run_synthesis(self):
        if self.matrix is None: return
        synth = TextSynthesizer(self.matrix, self.unique_words)
        text = synth.generate(length=15)
        self.txt_syn.delete(1.0, tk.END)
        self.txt_syn.insert(tk.END, text + "...")

if __name__ == "__main__":
    root = tk.Tk()
    TextTool(root)
    root.mainloop()