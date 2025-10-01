s = "ana are mere"
aparitii = {}

for ch in s:
    if ch != " ": 
        if ch in aparitii:
            aparitii[ch] += 1
        else:
            aparitii[ch] = 1

for lit, cnt in aparitii.items():
    print(lit, "apare de", cnt, "ori")
