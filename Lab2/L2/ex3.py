# Design an application by using the AI, which contains a GUI that allows the user to select a FASTA file. 
# The content of the FASTA file should be analysed by using a sliding window of 30  positions. 
# the content of each sliding window should be used in order to compute the relative frequencies
# of the nucleotides found in the alphabet of the sequence. 
# the output of the application should be a chart containing 4 signals: 
# one signal for each of the symbols of the alphabet of the sequence. 
# In order to plot this chart, first, the vectors that contain the frequency of 
# each symbol must be computed and only then plotted.

import tkinter as tk
from tkinter import filedialog, messagebox
import matplotlib.pyplot as plt


def load_fasta(filepath):
    with open(filepath, 'r') as f:
        lines = f.readlines()
    seq = ''.join(line.strip() for line in lines if not line.startswith('>')).upper()
    return seq

def compute_frequencies(seq, window_size=30):
    alphabet = sorted(set(seq))  # e.g., ['A', 'C', 'G', 'T']
    freq_vectors = {base: [] for base in alphabet}

    for i in range(len(seq) - window_size + 1):
        window = seq[i:i+window_size]
        for base in alphabet:
            freq = window.count(base) / window_size
            freq_vectors[base].append(freq)
    return freq_vectors


def open_and_analyze():
    filepath = filedialog.askopenfilename(
        title="Select FASTA File",
        filetypes=(("FASTA files", "*.fasta *.fa *.txt"), ("All files", "*.*"))
    )
    if not filepath:
        return
    
    try:
        seq = load_fasta(filepath)
        if len(seq) < 30:
            messagebox.showwarning("Warning", "Sequence too short for window size 30.")
            return
        
        freqs = compute_frequencies(seq, 30)
        
        plt.figure(figsize=(10, 6))
        for base, values in freqs.items():
            plt.plot(values, label=f"{base} frequency")
        
        plt.title("Relative Nucleotide Frequencies (Sliding Window = 30)")
        plt.xlabel("Window Position")
        plt.ylabel("Relative Frequency")
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.show()

    except Exception as e:
        messagebox.showerror("Error", str(e))



#gui
root = tk.Tk()
root.title("FASTA Nucleotide Frequency Analyzer")
root.geometry("400x150")

label = tk.Label(root, text="Select a FASTA file to analyze nucleotide frequencies:")
label.pack(pady=10)

button = tk.Button(root, text="Open FASTA File", command=open_and_analyze)
button.pack(pady=10)

root.mainloop()
