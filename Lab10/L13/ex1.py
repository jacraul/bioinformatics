import matplotlib.pyplot as plt
import collections

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

def main():
    # 1. Use the sequence S
    S = "CGGACTGATCTATCTAAAAAAAAAAAAAAAAAAAAAAAAAAACGTAGCATCTATCGATCTATCTAGCGATCTATCTACTACG"
    WINDOW_SIZE = 30
    
    print(f"Analyzing Sequence (Length: {len(S)})")
    print("-" * 40)

    
    total_cg = calculate_cg_content(S)
    total_ic = calculate_kappa_ic(S)
    
    print(f"Verification on Total Sequence:")
    print(f"Calculated CG%: {total_cg:.2f} (Target: 29.27)")
    print(f"Calculated IC:  {total_ic:.2f} (Target: 27.53)")
    
    if abs(total_ic - 27.53) < 2.0:
        print(">> VERIFICATION SUCCESSFUL: Values are within acceptable range (PromKappa Algorithm).")
    else:
        print(">> VERIFICATION WARNING: Values deviate.")
    
    print("-" * 40)

    print(f"Generating Pattern with Sliding Window (Size: {WINDOW_SIZE})...")
    cg_points, ic_points = analyze_sequence(S, WINDOW_SIZE)
    
    avg_cg = sum(cg_points) / len(cg_points)
    avg_ic = sum(ic_points) / len(ic_points)
    
    print(f"Center of Weight (Centroid):")
    print(f"CG: {avg_cg:.2f}")
    print(f"IC: {avg_ic:.2f}")
    
    plt.figure(figsize=(12, 6))
    
    plt.subplot(1, 2, 1)
    plt.scatter(cg_points, ic_points, alpha=0.6, c='blue', label='Window (30b)')
    plt.plot(cg_points, ic_points, alpha=0.2, c='blue') 
    plt.title('Promoter Pattern (Sliding Window)')
    plt.xlabel('C+G %')
    plt.ylabel('Kappa IC (PromKappa)')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.legend()

    plt.subplot(1, 2, 2)
    plt.scatter(avg_cg, avg_ic, s=200, c='red', marker='X', label='Center of Weight')
    plt.scatter(cg_points, ic_points, alpha=0.1, c='gray')
    
    plt.title('Center of Weight (Centroid)')
    plt.xlabel('Average C+G %')
    plt.ylabel('Average Kappa IC')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.legend()

    plt.tight_layout()
    
    plt.savefig('promoter_pattern_analysis.png')
    print(">> Charts generated and saved to 'promoter_pattern_analysis.png'")
    plt.show()

if __name__ == "__main__":
    main()