import matplotlib.pyplot as plt
from collections import Counter

with open('Project_L7\\L7\\dna.fasta', 'r') as f:
    sequence = ''.join(line.strip() for line in f if not line.startswith('>')).upper()

repetitions = Counter()
for length in range(6, 11):
    for i in range(len(sequence) - length + 1):
        pattern = sequence[i:i + length]
        repetitions[pattern] += 1

repetitions = {p: c for p, c in repetitions.items() if c > 1}

top_20 = sorted(repetitions.items(), key=lambda x: x[1], reverse=True)[:20]
patterns, frequencies = zip(*top_20) if top_20 else ([], [])

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 10))

ax1.bar(range(len(patterns)), frequencies, color='steelblue')
ax1.set_xlabel('DNA Pattern')
ax1.set_ylabel('Frequency')
ax1.set_title('Top 20 DNA Repetitions (6-10bp)')
ax1.set_xticks(range(len(patterns)))
ax1.set_xticklabels(patterns, rotation=45, ha='right')
for i, v in enumerate(frequencies):
    ax1.text(i, v, str(v), ha='center', va='bottom')

length_freq = Counter()
for p, c in repetitions.items():
    length_freq[len(p)] += c
lengths = sorted(length_freq.keys())
counts = [length_freq[l] for l in lengths]

ax2.bar(lengths, counts, color='coral')
ax2.set_xlabel('Pattern Length (bp)')
ax2.set_ylabel('Total Frequency')
ax2.set_title('Repetition Frequency by Length')
ax2.set_xticks(lengths)
for l, c in zip(lengths, counts):
    ax2.text(l, c, str(c), ha='center', va='bottom')

plt.tight_layout()
plt.savefig('screenshot1.png', dpi=300)
print(f"Found {len(repetitions)} repetitive patterns")
print(f"Top pattern: {top_20[0][0]} appears {top_20[0][1]} times")
plt.show()