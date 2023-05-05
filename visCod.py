import re
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter

# Codon table for translation
rna2aa = {"UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
          "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
          "UAU": "Y", "UAC": "Y", "UAA": "STOP", "UAG": "STOP",
          "UGU": "C", "UGC": "C", "UGA": "STOP", "UGG": "W",
          "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
          "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
          "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
          "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
          "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
          "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
          "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
          "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
          "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
          "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
          "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
          "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G", }

# Open codon table file
with open("eColiCodon.txt", "r") as codonTable:
    # Parse frequencies from table
    codTab = {}
    for line in codonTable:
        tempLine = line.split(">")
        for i in range(0, len(tempLine), 3):
            codTab[tempLine[i]] = float(tempLine[i + 1].replace("(", ""))

# Open DNA sequence file
with open("Ecoli.txt", "r") as seqFile:
    geneSeq = seqFile.read()

# Translate DNA sequence to RNA and verify sequence is valid
rnaGene = geneSeq.replace("T", "U")
rnaGeneSeq = rnaGene.replace("\n", "")

remove_lower = lambda text: re.sub('[a-z]', '', text)
rnaGene = remove_lower(rnaGene)

# Get frequencies of each codon for later plotting
aaSeq = ""
locCodFreq = []
counterList = []
for i in range(0, len(geneSeq), 3):
    a = rnaGeneSeq[i:i + 3]
    counterList.append(a)
    locCodFreq.append(codTab.get(a))
    aaSeq += str(rna2aa.get(a))

aaSeq = remove_lower(aaSeq)

# Remove none type objects caused by STOP in codon table
aaSeq = aaSeq.replace("None", "")  # Remove "None" strings
aaSeq = aaSeq.replace("STOP", "")  # Remove "STOP" strings
aaSeq = aaSeq.upper()  # Convert all characters to uppercase
lenAASeq = len(aaSeq)  # Get the length of the amino acid sequence
print(aaSeq)

# Count the frequency of codons and amino acids
print(Counter(counterList))  # Counter of codon frequency
print(Counter(aaSeq))  # Counter of amino acid frequency
print(len(Counter(aaSeq)))  # Number of unique amino acids
print(len(aaSeq))  # Total number of amino acids
print("# of old codons used:{}\n# of new codons used:{}".format(len(codTab), len(Counter(counterList))))

# Create lists to store amino acid sequence and codon frequency
aaSeqList = []
locFreq = []

# Get the amino acid sequence and corresponding codon frequencies
lenCodFreq = len(locCodFreq)
const = lenCodFreq - lenAASeq
for i in range(0, len(aaSeq)):
    aaSeqList.append(aaSeq[i])
for i in range(0, len(locCodFreq) - const):
    locFreq.append(locCodFreq[i])

print(len(locCodFreq))

# Plot the histogram of codon frequency vs. amino acid sequence
plt.xlabel('AA')
plt.ylabel('Freq')
plt.title("Histogram of Freq")

# Plot the histogram using matplotlib bar()
indexes = np.arange(len(aaSeq))
plt.bar(indexes, locFreq)
plt.xticks(indexes, aaSeqList)
plt.show()
