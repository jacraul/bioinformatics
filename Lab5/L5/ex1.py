import random
import time

# --- Step 1: Take the DNA sequence from the provided file ---

file_name = 'C:\\Users\\RAUL\\OneDrive\\Desktop\\Politehnica\\Anul4\\Bioinformatics\\Project_L5\\L5\\seq.fasta'
sequence_lines = []
original_sequence = ""

# Read the FASTA file
try:
    with open(file_name, 'r') as f:
        for line in f:
            if not line.startswith('>'):
                sequence_lines.append(line.strip())

    original_sequence = "".join(sequence_lines)
    seq_len = len(original_sequence)

    print(f"--- Step 1: Loaded Sequence ---")
    print(f"Successfully loaded {seq_len} bp from '{file_name}'.")
    print("-" * 30)

    # --- Step 2 & 3: Take 2000 samples and store in a list ---

    print("--- Step 2 & 3: Sampling and Storing ---")
    samples = []
    num_samples = 2000
    min_len = 100
    max_len = 150

    for i in range(num_samples):
        sample_length = random.randint(min_len, max_len)
        max_start = seq_len - sample_length
        start_index = random.randint(0, max_start)
        end_index = start_index + sample_length
        sample = original_sequence[start_index:end_index]
        samples.append(sample)

    print(f"Generated and stored {len(samples)} samples in the 'samples' list.")
    print("-" * 30)
    
    # --- Step 4: Rebuild using a Greedy Overlap Algorithm ---
    
    print("--- Step 4: Attempting to Rebuild with Overlap Logic ---")
    
    # WARNING: Running this on 2000 samples is computationally massive (O(n^2)).
    # We will use a small subset of 100 samples to demonstrate the logic.
    if not samples:
        print("No samples to rebuild.")
        exit()

    rebuild_list = samples[:100]
    print(f"Using a small subset of {len(rebuild_list)} samples to rebuild.")

    # Set our parameters
    min_overlap = 70  # The smallest overlap we will trust
    max_sample_len = 150 # The max possible length of a sample

    # Start with a "seed" fragment from our list
    reconstructed_sequence = rebuild_list.pop(0)
    start_time = time.time()

    while len(rebuild_list) > 0:
        found_match = False
        best_match_sample = None
        best_match_index = -1
        best_match_k = -1 # This will store the length of the overlap

        # Find the *longest possible* overlap
        # We search from the max length down to our minimum
        for k in range(max_sample_len - 1, min_overlap - 1, -1):
            if k >= len(reconstructed_sequence):
                continue # Skip if our sequence is shorter than the overlap
            
            # Get the end of our growing sequence
            end_of_sequence = reconstructed_sequence[-k:]
            
            # Search the *entire* remaining list for a match
            for i, sample in enumerate(rebuild_list):
                if sample.startswith(end_of_sequence):
                    # Found a match! This is our new best.
                    best_match_sample = sample
                    best_match_index = i
                    best_match_k = k
                    found_match = True
                    break  # Stop searching the list
            
            if found_match:
                break # Stop searching for shorter overlaps

        # --- Now, stitch the best match we found ---
        if found_match:
            # Add the new, non-overlapping part
            reconstructed_sequence += best_match_sample[best_match_k:]
            
            # Remove the sample we just used
            rebuild_list.pop(best_match_index)
            
            print(f"  ...Stitched! New length: {len(reconstructed_sequence)} bp. ({len(rebuild_list)} samples left)")
        else:
            # If we went through all fragments and all 'k' lengths and
            # found no match, our assembly is stuck.
            print("  ...Assembly STUCK. No overlap found for this fragment.")
            break # Exit the 'while' loop

    # --- Final Report ---
    end_time = time.time()
    print("-" * 30)
    print("Assembly Finished.")
    print(f"Time taken: {end_time - start_time:.2f} seconds.")
    print(f"Final reconstructed length: {len(reconstructed_sequence)} bp (from 100 fragments)")
    print(f"Original sequence length: {seq_len} bp")
    
    print("\nNote: This is a 'contig' (a continuous fragment).")
    print("It's incomplete because we only used 100 samples and got stuck.")
    print("This demonstrates the *logic* required for Step 4.")


except FileNotFoundError:
    print(f"Error: The file '{file_name}' was not found.")
except Exception as e:
    print(f"An error occurred: {e}")