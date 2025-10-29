def translate_cds(sequence): 
    codon_table = {
        'UUU': 'Phe', 'UUC': 'Phe', 'UUA': 'Leu', 'UUG': 'Leu',
        'UCU': 'Ser', 'UCC': 'Ser', 'UCA': 'Ser', 'UCG': 'Ser',
        'UAU': 'Tyr', 'UAC': 'Tyr', 'UAA': 'Stop', 'UAG': 'Stop',
        'UGU': 'Cys', 'UGC': 'Cys', 'UGA': 'Stop', 'UGG': 'Trp',
        'CUU': 'Leu', 'CUC': 'Leu', 'CUA': 'Leu', 'CUG': 'Leu',
        'CCU': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
        'CAU': 'His', 'CAC': 'His', 'CAA': 'Gln', 'CAG': 'Gln',
        'CGU': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg',
        'AUU': 'Ile', 'AUC': 'Ile', 'AUA': 'Ile', 'AUG': 'Met',
        'ACU': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
        'AAU': 'Asn', 'AAC': 'Asn', 'AAA': 'Lys', 'AAG': 'Lys',
        'AGU': 'Ser', 'AGC': 'Ser', 'AGA': 'Arg', 'AGG': 'Arg',
        'GUU': 'Val', 'GUC': 'Val', 'GUA': 'Val', 'GUG': 'Val',
        'GCU': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
        'GAU': 'Asp', 'GAC': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu',
        'GGU': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly'
    }

    protein = []
    
    mrna = sequence.upper().replace('T', 'U')

    for i in range(0, len(mrna), 3):
        codon = mrna[i:i+3]

        if len(codon) < 3:
            break
        amino_acid = codon_table.get(codon, '?')
        if amino_acid == 'Stop':
            break
        protein.append(amino_acid)

    return "".join(protein)

dna_sequence = "ATGGCCGCTTTCTAA"
print(f"DNA Sequence: {dna_sequence}")
print(f"Protein:      {translate_cds(dna_sequence)}")
print("-" * 20)

rna_sequence = "AUGGCCGCUUUCUAA"
print(f"RNA Sequence: {rna_sequence}")
print(f"Protein:      {translate_cds(rna_sequence)}")
print("-" * 20)

unknown_sequence = "AUGNNNCAACUGUAG"
print(f"Unknown Seq:  {unknown_sequence}")
print(f"Protein:      {translate_cds(unknown_sequence)}")
print("-" * 20)

short_sequence = "UUCAUGGGU"
print(f"Short Seq:    {short_sequence}")
print(f"Protein:      {translate_cds(short_sequence)}") 