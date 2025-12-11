import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import sys


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
        print(f"Eroare: Fisierul {filename} lipseste.")
        return ""
    return "".join(sequence)


def get_relative_alignment_data(seq_flu, seq_cov, samples=600):

    data = []
    len_f = len(seq_flu)
    len_c = len(seq_cov)
    
    print(f"Generare hartă relativă ({samples} puncte)...")
    
    for p in range(samples):
        idx_f = int((p / samples) * len_f)
        idx_c = int((p / samples) * len_c)
        
        # Limite
        idx_f = min(idx_f, len_f - 1)
        idx_c = min(idx_c, len_c - 1)
        
        char_f = seq_flu[idx_f]
        char_c = seq_cov[idx_c]
        is_match = (char_f == char_c)
        
        data.append({
            'x': p,
            'flu_char': char_f,
            'cov_char': char_c,
            'is_match': is_match,
            'flu_idx': idx_f,
            'cov_idx': idx_c
        })
        
    return data


class ZoomableGenomeMap:
    def __init__(self, data, width):
        self.data = data
        self.width = width
        
        self.fig, self.ax = plt.subplots(figsize=(15, 8))
        self.text_artists_flu = []
        self.text_artists_cov = []
        
        self.zoom_threshold = 120 
        self.text_visible = False

    def draw(self):
        self.ax.set_xlim(-10, self.width + 10)
        self.ax.set_ylim(-0.5, 1.5)
        self.ax.axis('off')

        # --- LINII DE FUNDAL ---
        self.ax.hlines(1, 0, self.width, colors='blue', linewidth=4, alpha=0.2)
        self.ax.text(-20, 1, "INFLUENZA", va='center', ha='right', fontsize=12, fontweight='bold', color='blue')
        
        self.ax.hlines(0, 0, self.width, colors='red', linewidth=4, alpha=0.2)
        self.ax.text(-20, 0, "COVID-19", va='center', ha='right', fontsize=12, fontweight='bold', color='red')

        segments = []
        for d in self.data:
            if d['is_match']:
                segments.append([(d['x'], 0), (d['x'], 1)])
        
        lc = LineCollection(segments, colors='green', alpha=0.5, linewidths=2)
        self.ax.add_collection(lc)

        print("Se pregătesc literele (ascunse)...")
        for d in self.data:
            color = 'green' if d['is_match'] else 'darkred'
            weight = 'bold' if d['is_match'] else 'normal'
            
            t1 = self.ax.text(d['x'], 1, d['flu_char'], ha='center', va='center', 
                              fontsize=10, color=color, fontweight=weight, visible=False, clip_on=True)
            self.text_artists_flu.append(t1)
            
            t2 = self.ax.text(d['x'], 0, d['cov_char'], ha='center', va='center', 
                              fontsize=10, color=color, fontweight=weight, visible=False, clip_on=True)
            self.text_artists_cov.append(t2)

        self.ax.set_title("Harta Aliniere Genomică\nFolosește Lupa (Zoom Box) pe o zonă mică pentru a vedea LITERELE", fontsize=14)
        
        self.ax.callbacks.connect('xlim_changed', self.on_xlims_change)
        
        plt.tight_layout()
        plt.show()

    def on_xlims_change(self, event_ax):
        xlim = event_ax.get_xlim()
        visible_range = xlim[1] - xlim[0]
        
        should_be_visible = visible_range < self.zoom_threshold

        if should_be_visible != self.text_visible:
            self.text_visible = should_be_visible
            print(f"Zoom detectat: {int(visible_range)} puncte. Afisare text: {self.text_visible}")
            
            for t in self.text_artists_flu:
                t.set_visible(self.text_visible)
            for t in self.text_artists_cov:
                t.set_visible(self.text_visible)
            
            self.fig.canvas.draw_idle()


if __name__ == "__main__":
    f_flu = 'Influenza.fasta'
    f_cov = 'Covid19.fasta'

    s_flu = read_fasta_robust(f_flu)
    s_cov = read_fasta_robust(f_cov)

    if s_flu and s_cov:
        data_points = get_relative_alignment_data(s_flu, s_cov, samples=600)
        
        if data_points:
            viz = ZoomableGenomeMap(data_points, width=600)
            viz.draw()
        else:
            print("Eroare: Nu sunt date.")
    else:
        print("Lipsesc fisierele FASTA.")