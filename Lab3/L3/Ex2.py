import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import re
import math
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

def calculate_tms(subsequence, na_concentration):
    a_count = subsequence.count('A')
    t_count = subsequence.count('T')
    g_count = subsequence.count('G')
    c_count = subsequence.count('C')
    
    total_bases = len(subsequence)
    if total_bases == 0 or total_bases != (a_count + t_count + g_count + c_count):
        return -1, -1

    # Formula 1: Wallace Rule (Basic)
    tm_wallace = (2 * (a_count + t_count)) + (4 * (g_count + c_count))
    
    # Formula 2: Marmur & Doty Formula (Advanced)
    gc_content = (g_count + c_count) / total_bases * 100
    # Tm = 81.5 + 16.6 * log10([Na+]) + 0.41 * (%GC) - 600 / length
    tm_advanced = 81.5 + (16.6 * math.log10(na_concentration)) + (0.41 * gc_content) - (600 / total_bases)
    
    return tm_wallace, round(tm_advanced, 2)

class TmAnalyzerApp:
    def __init__(self, root):
        self.root = root
        self.root.title("DNA Melting Temperature Analyzer")
        self.root.geometry("800x700")
        self.root.minsize(600, 600)

        style = ttk.Style(self.root)
        style.theme_use("clam")
        style.configure("TButton", padding=6, relief="flat", background="#0078D7", foreground="white")
        style.map("TButton", background=[('active', '#005a9e')])
        style.configure("TLabel", padding=5, font=('Helvetica', 10))
        style.configure("Treeview.Heading", font=('Helvetica', 10, 'bold'))
        style.configure("TEntry", padding=5)

        top_frame = ttk.Frame(self.root, padding="10")
        top_frame.pack(fill=tk.X, padx=10, pady=5)

        self.filepath_label = ttk.Label(top_frame, text="No file selected.", wraplength=400)
        self.filepath_label.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(0, 10))
        
        self.browse_button = ttk.Button(top_frame, text="Browse...", command=self.browse_file)
        self.browse_button.pack(side=tk.LEFT)

        control_frame = ttk.Frame(self.root, padding="10")
        control_frame.pack(fill=tk.X, padx=10, pady=5)
        
        ttk.Label(control_frame, text="Window Size:").pack(side=tk.LEFT, padx=(0, 5))
        self.window_size_var = tk.StringVar(value="8")
        self.window_size_entry = ttk.Entry(control_frame, textvariable=self.window_size_var, width=5)
        self.window_size_entry.pack(side=tk.LEFT, padx=(0, 20))

        ttk.Label(control_frame, text="Na+ Conc. (M):").pack(side=tk.LEFT, padx=(0, 5))
        self.na_conc_var = tk.StringVar(value="0.1")
        self.na_conc_entry = ttk.Entry(control_frame, textvariable=self.na_conc_var, width=7)
        self.na_conc_entry.pack(side=tk.LEFT, padx=(0, 20))

        self.analyze_button = ttk.Button(control_frame, text="Analyze Sequence", command=self.run_analysis)
        self.analyze_button.pack(side=tk.LEFT)

        results_frame = ttk.Frame(self.root, padding="10")
        results_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)

        self.tree = ttk.Treeview(results_frame, columns=("Position", "Sequence", "Tm_Basic", "Tm_Advanced"), show="headings")
        self.tree.heading("Position", text="Position")
        self.tree.heading("Sequence", text="Sequence")
        self.tree.heading("Tm_Basic", text="Tₘ Basic (°C)")
        self.tree.heading("Tm_Advanced", text="Tₘ Advanced (°C)")
        
        self.tree.column("Position", width=80, anchor=tk.CENTER)
        self.tree.column("Sequence", width=120, anchor=tk.CENTER)
        self.tree.column("Tm_Basic", width=120, anchor=tk.CENTER)
        self.tree.column("Tm_Advanced", width=150, anchor=tk.CENTER)

        scrollbar = ttk.Scrollbar(results_frame, orient=tk.VERTICAL, command=self.tree.yview)
        self.tree.configure(yscroll=scrollbar.set)
        
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        self.tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        
        plot_frame = ttk.Frame(self.root, padding="10")
        plot_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        self.fig = Figure(figsize=(6, 4), dpi=100)
        self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.fig, master=plot_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        self.status_var = tk.StringVar(value="Ready.")
        self.status_bar = ttk.Label(self.root, textvariable=self.status_var, relief=tk.SUNKEN, anchor=tk.W, padding=5)
        self.status_bar.pack(side=tk.BOTTOM, fill=tk.X)
        
        self.selected_filepath = ""
        self.analysis_results = []

    def browse_file(self):
        filepath = filedialog.askopenfilename(
            title="Select a FASTA file",
            filetypes=(("FASTA files", "*.fasta *.fa *.fna"), ("All files", "*.*"))
        )
        if filepath:
            self.selected_filepath = filepath
            display_path = filepath if len(filepath) < 60 else "..." + filepath[-57:]
            self.filepath_label.config(text=display_path)
            self.status_var.set(f"File selected: {display_path}")

    def plot_results(self):
        self.ax.clear()
        if not self.analysis_results:
            self.canvas.draw()
            return
            
        positions = [res[0] for res in self.analysis_results]
        tm_basic_values = [res[2] for res in self.analysis_results]
        tm_advanced_values = [res[3] for res in self.analysis_results]

        self.ax.plot(positions, tm_basic_values, label='Tₘ Basic (Wallace Rule)', marker='o', linestyle='-', color='blue')
        self.ax.plot(positions, tm_advanced_values, label='Tₘ Advanced (Marmur & Doty)', marker='x', linestyle='--', color='red')
        
        self.ax.set_title("Melting Temperature Analysis")
        self.ax.set_xlabel("Window Start Position")
        self.ax.set_ylabel("Temperature (°C)")
        self.ax.legend()
        self.ax.grid(True)
        self.fig.tight_layout()
        self.canvas.draw()

    def run_analysis(self):
        for item in self.tree.get_children():
            self.tree.delete(item)
        self.analysis_results = []

        if not self.selected_filepath:
            messagebox.showwarning("Input Missing", "Please select a FASTA file first.")
            return
            
        try:
            window_size = int(self.window_size_var.get())
            if window_size <= 0:
                raise ValueError
        except ValueError:
            messagebox.showerror("Invalid Input", "Window size must be a positive integer.")
            return

        try:
            na_concentration = float(self.na_conc_var.get())
            if na_concentration <= 0:
                raise ValueError
        except ValueError:
            messagebox.showerror("Invalid Input", "Na+ concentration must be a positive number.")
            return

        self.status_var.set("Analyzing... Please wait.")
        self.root.update_idletasks()

        sequence = self.parse_fasta_internal()
        if not sequence:
            self.status_var.set("Analysis failed. Please check the file.")
            self.plot_results()
            return

        if len(sequence) < window_size:
            messagebox.showerror("Error", "The sequence is shorter than the specified window size.")
            self.status_var.set("Analysis failed.")
            self.plot_results()
            return
        
        for i in range(len(sequence) - window_size + 1):
            subsequence = sequence[i : i + window_size]
            tm_wallace, tm_advanced = calculate_tms(subsequence, na_concentration)
            
            if tm_wallace != -1:
                self.analysis_results.append((i + 1, subsequence, tm_wallace, tm_advanced))
                self.tree.insert("", tk.END, values=(i + 1, subsequence, f"{tm_wallace}", f"{tm_advanced}"))
            else:
                self.tree.insert("", tk.END, values=(i + 1, subsequence, "Invalid", "Invalid"))
        
        self.status_var.set(f"Analysis complete. Found {len(self.analysis_results)} valid windows.")
        self.plot_results()

    def parse_fasta_internal(self):
        if not self.selected_filepath:
            return None
        try:
            with open(self.selected_filepath, 'r') as file:
                sequence_lines = []
                for line in file:
                    if not line.startswith('>'):
                        sequence_lines.append(re.sub(r'[^a-zA-Z]', '', line))
                return "".join(sequence_lines).upper()
        except Exception:
            return None

if __name__ == "__main__":
    root = tk.Tk()
    app = TmAnalyzerApp(root)
    root.mainloop()

