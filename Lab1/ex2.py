import numpy as np
s = "ACGGGGCATATGCGC"
aparitii = {}
sumall = 0

for ch in s:
    sumall += 1
    if ch != " ": 
        if ch in aparitii:
            aparitii[ch] += 1
        else:
            aparitii[ch] = 1

for lit, cnt in aparitii.items():
    p=cnt/sumall*100
    print(lit, ":", p, "%")
