S = "ATTGTCCCAATCTGTTG"

dinucleotides = []
trinucleotides = []

for i in range(len(S)):
    if i + 1 < len(S):
        di = S[i:i+2]
        if di not in dinucleotides:
            dinucleotides.append(di)
    if i + 2 < len(S):
        tri = S[i:i+3]
        if tri not in trinucleotides:
            trinucleotides.append(tri)

print("Existing Dinucleotides:")
print(", ".join(dinucleotides))

print("\nExisting Trinucleotides:")
print(", ".join(trinucleotides))
