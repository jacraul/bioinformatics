def read_fasta(filename):
    seq = ""
    with open(filename, "r") as f:
        for line in f:
            if line.startswith(">"): 
                continue
            seq += line.strip()       
    return seq


s = read_fasta("input.fasta")

aparitii = {}
sumall = 0

for ch in s:
    sumall += 1
    if ch in aparitii:
        aparitii[ch] += 1
    else:
        aparitii[ch] = 1

print("Simbols and percentage:\n")
for lit, cnt in aparitii.items():
    p = cnt / sumall * 100
    print(f"{lit} : {p:.2f}%")
