import matplotlib.pyplot as plt
import collections
import os
import re

def calculate_cg_content(sequence):

    sequence = sequence.upper()
    length = len(sequence)
    if length == 0:
        return 0.0
    
    c_count = sequence.count('C')
    g_count = sequence.count('G')
    
    cg_percent = ((c_count + g_count) / length) * 100
    return cg_percent

def calculate_kappa_ic(sequence):
    sequence = sequence.upper()
    N = len(sequence)
    if N < 2:
        return 0.0
    
    total_percentage = 0.0
    
    for shift in range(1, N):
        matches = 0
        length_of_overlap = N - shift
        
        for i in range(length_of_overlap):
            if sequence[i] == sequence[i + shift]:
                matches += 1
        
        if length_of_overlap > 0:
            pct = (matches / length_of_overlap) * 100.0
            total_percentage += pct
            
    ic_value = total_percentage / (N - 1)
    return ic_value

def analyze_sequence(sequence, window_size=30):

    cg_values = []
    ic_values = []
    
    for i in range(len(sequence) - window_size + 1):
        window = sequence[i : i + window_size]
        
        cg = calculate_cg_content(window)
        ic = calculate_kappa_ic(window)
        
        cg_values.append(cg)
        ic_values.append(ic)
        
    return cg_values, ic_values

def read_fasta(file_path):

    header = None
    sequence = []
    
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            if line.startswith(">"):
                if header:
                    yield header, "".join(sequence)
                header = line[1:] # Remove >
                sequence = []
            else:
                sequence.append(line)
        
        if header and sequence:
            yield header, "".join(sequence)

def process_promoters_combined(fasta_file, output_folder="ODS_Patterns"):
 
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
        print(f"Created output folder: {output_folder}")
    
    print(f"Processing {fasta_file}...")
    
    plt.figure(figsize=(12, 10))
    
    count = 0
    
    for header, sequence in read_fasta(fasta_file):
        count += 1
        
        cg_points, ic_points = analyze_sequence(sequence, window_size=30)
        
        if not cg_points:
            print(f"[{count}] Warning: Sequence too short. Skipped.")
            continue

        plt.scatter(cg_points, ic_points, alpha=0.4, s=15, label=f"Seq {count}" if count <= 5 else "") 

        if count % 10 == 0:
            print(f"Processed {count} sequences...")

    print(f"Finished processing {count} sequences.")

    plt.title(f'Combined ODS Patterns (N={count} Sequences)')
    plt.xlabel('C+G %')
    plt.ylabel('Kappa IC')
    plt.grid(True, linestyle='--', alpha=0.3)
    
    output_filename = f"{output_folder}/Combined_ODS_Pattern.png"
    plt.savefig(output_filename)
    plt.close()
    
    print(f"Merged ODS pattern saved to: {output_filename}")

def main():
    FASTA_FILE = "promoters.fasta"
    
    if not os.path.exists(FASTA_FILE):
        print(f"Error: {FASTA_FILE} not found. Please upload the file.")
    else:
        process_promoters_combined(FASTA_FILE)

if __name__ == "__main__":
    main()