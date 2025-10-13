import tkinter as tk
from itertools import product
from collections import Counter

bases = "ATGC"

def generate_combinations(length):
    return [''.join(p) for p in product(bases, repeat=length)]

def calculate_percentages(seq, length):
    combos = generate_combinations(length)
    counts = Counter(seq[i:i+length] for i in range(len(seq)-length+1))
    total = sum(counts.values())
    percentages = {k: (counts.get(k, 0)/total)*100 for k in combos}
    return percentages

def format_grid(percent_dict, cols=5):
    items = list(percent_dict.items())
    lines = []
    for i in range(0, len(items), cols):
        row = items[i:i+cols]
        line = " | ".join(f"{k}: {v:.2f}%" for k, v in row)
        lines.append(line)
    return "\n".join(lines)

def calculate():
    seq = entry_seq.get().upper().replace(" ", "")
    length = int(entry_length.get())
    if length < 1:
        output_text.set("Length must be >= 1")
        return
    percentages = calculate_percentages(seq, length)
    output_text.set(format_grid(percentages))

# GUI 
root = tk.Tk()
root.title("Nucleotide Percentage Calculator")

tk.Label(root, text="Enter DNA Sequence:").grid(row=0, column=0, sticky="w")
entry_seq = tk.Entry(root, width=50)
entry_seq.grid(row=0, column=1, columnspan=4)

tk.Label(root, text="Combination Length (2 or 3):").grid(row=1, column=0, sticky="w")
entry_length = tk.Entry(root, width=5)
entry_length.grid(row=1, column=1, sticky="w")
entry_length.insert(0, "2")

tk.Button(root, text="Calculate Percentages", command=calculate).grid(row=1, column=2, columnspan=2)

output_text = tk.StringVar()
output_label = tk.Label(root, textvariable=output_text, justify="left", font=("Courier", 10))
output_label.grid(row=2, column=0, columnspan=5, sticky="w")

entry_seq.insert(0, "ATTGTCCCAATCTGTTG")

root.mainloop()