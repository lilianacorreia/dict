# dict

#1. Write a Python script to transcribe this sequence to mRNA

sequence = ('ATGTCTAAAGGTGAAGAATTATTCACTGGTGTTGTCCCAATTTTGGTTGAATTAGATGGTGATGTTAATGGT'
            'CACAAATTTTCTGTCTCCGGTGAAGGTGAAGGTGATGCTACTTACGGTAAATTGACCTTAAAATTTATTTGT'
            'ACTACTGGTAAATTGCCAGTTCCATGGCCAACCTTAGTCACTACTTTCGGTTATGGTGTTCAATGTTTTGCT'
'AGATACCCAGATCATATGAAACAACATGACTTTTTCAAGTCTGCCATGCCAGAAGGTTATGTTCAAGAAAGA'
'ACTATTTTTTTCAAAGATGACGGTAACTACAAGACCAGAGCTGAAGTCAAGTTTGAAGGTGATACCTTAGTT'
'AATAGAATCGAATTAAAAGGTATTGATTTTAAAGAAGATGGTAACATTTTAGGTCACAAATTGGAATACAAC'
'TATAACTCTCACAATGTTTACATCATGGCTGACAAACAAAAGAATGGTATCAAAGTTAACTTCAAAATTAGA'
'CACAACATTGAAGATGGTTCTGTTCAATTAGCTGACCATTATCAACAAAATACTCCAATTGGTGATGGTCCA'
'GTCTTGTTACCAGACAACCATTACTTATCCACTCAATCTGCCTTATCCAAAGATCCAAACGAAAAGAGAGAC'
'CACATGGTCTTGTTAGAATTTGTTACTGCTGCTGGTATTACCCATGGTATGGATGAATTGTACAAATAA')
dna2cdna = {"A": "T", "T": "A", "C": "G", "G": "C"}
cDNA_sequence = str() # Creating an empty string.

for nucleotide in sequence:
    cDNA_sequence += dna2cdna[nucleotide]
print('DNA:', sequence, '\n\ncDNA:', cDNA_sequence )

dna2rna = {"A": "U", "T": "A", "C": "G", "G": "C"}
rna_sequence = str() # Creating an empty string.

for nucleotide in cDNA_sequence:
    rna_sequence += dna2rna[nucleotide]
print('\nRNA:', rna_sequence)

#2. obtain the translation (protein sequence) of the Green Fluorescent Protein 

codons = []
for i in range(0, len(rna_sequence), 3):
    codons.append(rna_sequence[i:i+3])

# translation
protein = []
for i in codons:
    protein.append(codon2aa[i])

protein = ''.join(protein)
protein


#2. obtain the translation (protein sequence) of the Green Fluorescent Protein 

codon2aa = {
"UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L", "CUU": "L",
"CUC": "L", "CUA": "L", "CUG": "L", "AUU": "I", "AUC": "I",
"AUA": "I", "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
"UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S", "AGU": "S",
"AGC": "S", "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
"ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCU": "A",
"GCC": "A", "GCA": "A", "GCG": "A", "UAU": "Y", "UAC": "Y",
"CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "AAU": "N",
"AAC": "N", "AAA": "K", "AAG": "K", "GAU": "D", "GAC": "D",
"GAA": "E", "GAG": "E", "UGU": "C", "UGC": "C", "UGG": "W",
"CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R",
"AGG": "R", "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G",
"AUG": "<Met>", "UAA": "<STOP>", "UAG": "<STOP>", "UGA": "<STOP>"
}

codons = []
for i in range(0, len(rna_sequence), 3):
    codons.append(rna_sequence[i:i+3])

# translation
protein = []
for i in codons:
    protein.append(codon2aa[i])

protein = ''.join(protein)
print(protein)
#Analyse the proteins involved in the canonical Wnt signaling pathway

canonical = ['Wnt', 'LPR5/6', 'Frizzled', 'Dsh', 'APC', 'Axin', 'GSK3', 'CK1', 'Beta-Cathenin']
noncanonical = ['Wnt', 'Frizzled', 'Dsh','DAAM1', 'Profilin', 'Rho','ROCK']

# 1. What protein(s), if any, do the two pathways have in common?
common = []
for i in canonical:
    for j in noncanonical:
        if i == j:
            common.append(i)
print(common)

#2. Which proteins are unique to only one of the pathways (it is not important which pathway they belong to)?
unique = []
for i in (canonical + noncanonical):
    if i not in common:
        unique.append(i)
        
print(unique)

#3. Is there a protein called "APC" in either collection (it is not important which one)?
answer_3 = "APC" in (canonical + noncanonical)

print(answer_3)

# 4. Provide a non-repetitive collection of the proteins involved in both pathways; and ensure that 
#    the collection cannot be accidentally modified later

non_repetitive = canonical.copy()
for i in noncanonical:
    if i not in canonical:
        non_repetitive.append(i)
    
non_repetitive = tuple(non_repetitive)
print(non_repetitive)


##########################################################################################

#answer 1.

numbers = [7, 16, 0.3, 0, 15, -4, 5, 3, 15]

minimum = numbers[0]

for value in numbers:
    if value < minimum:
        minimum = value


print('The minimum value is:', minimum, 'The minimum was set to the first number of the array because any given defined number could not be included in the numbers list.',
'for example minimum = -20. The printed result would be -20 even though -20 is not in the numbers list.',
'By the way the for loop is defined, we could had chosen any',
'given number of the list that the printed result would always be -4.')

#As example for the question 1 if minimum would be minimum = -20
numbers = [7, 16, 0.3, 0, 15, -4, 5, 3, 15]

minimum = -20

for value in numbers:
    if value < minimum:
        minimum = value

print('if minimum=', minimum, 'it would be wrong')

#answer 2
numbers = [0, -2.1, 1.5, 3]

for item in numbers:
    new_value= sum(numbers)
print('The sum of numbers in the array is:',new_value)

#answer 3
numbers= [2, 1, 3]

for value in numbers:
    print(str(value)*value)

#answer 4
#write a script that using at most two for loops, finds the variance of the numbers, and display the mean, and the variance.
#Note that you will need to calculate the mean as a part of your calculations to find the variance.

numbers = [7, 16, 0.3, 0, 15, -4, 5, 3, 15]


mean_numbers = sum(numbers)/ len(numbers)

x = []
for i in numbers:
    x.append((i - mean_numbers)**2)

variance = sum(x)/len(numbers)

print("Mean: ", round(mean_numbers, 3), "\nVariance: ", round(variance, 3))


