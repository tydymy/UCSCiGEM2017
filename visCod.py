import re
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from collections import Counter
'''
Program for plotting realtive RNA codon frequency vs amino acid. Gene optimization. 
'''

#Codon table for translation
rna2aa = {"UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
    "UCU":"S", "UCC":"s", "UCA":"S", "UCG":"S",
    "UAU":"Y", "UAC":"Y", "UAA":"STOP", "UAG":"STOP",
    "UGU":"C", "UGC":"C", "UGA":"STOP", "UGG":"W",
    "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
    "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
    "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
    "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G",}

#import codon table
codonTable = open("eColiCodon.txt","r")

#Parse frequncies from table
codTab = {}
for lines in codonTable:
    tempLine = lines.split(">")
    for i in range(0,len(tempLine),3):
        codTab[tempLine[i]] = float(tempLine[i+1].replace("(",""))
        
#import DNA sequence
with open("Ecoli.txt", "r") as seqFile:
    geneSeq = seqFile.read()
    
#translate to rna and verify sequence is valid
rnaGene = geneSeq.replace("T", "U")
rnaGeneSeq = rnaGene.replace("\n", "")

remove_lower = lambda text: re.sub('[a-z]', '', text)

rnaGene = remove_lower(rnaGene)

#get frequencies of each codon for later plotting
aaSeq = ""
locCodFreq = []
counterList = []
for i in range(0, len(geneSeq),3):
    a = rnaGeneSeq[i:i+3]
    counterList.append(a)
    locCodFreq.append(codTab.get(a))
    aaSeq += str(rna2aa.get(a))

aaSeq = remove_lower(aaSeq)

#remove none type objects caused by STOP in codon table
aaSeq = aaSeq.replace("None","")
aaSeq = aaSeq.replace("STOP", "")
aaSeq = aaSeq.upper()
lenAASeq = len(aaSeq)
print(aaSeq)



print(Counter(counterList))
print(Counter(aaSeq))
print(len(Counter(aaSeq)))
print(len(aaSeq))
print("# of old codons used:{}\n# of new codons used:{}".format(len(codTab), len(Counter(counterList))))
aaSeqList = []
locFreq = []

#get aa sequence, plot vs codon frequency list
lenCodFreq = len(locCodFreq)
const = lenCodFreq - lenAASeq
for i in range(0, len(aaSeq)):
    aaSeqList.append(aaSeq[i])
for i in range(0, len(locCodFreq)-const):#-8):
    locFreq.append(locCodFreq[i])

print(len(locCodFreq))
#print(locCodFreq)
#syntax for plotting
# the histogram of the data
plt.xlabel('AA')
plt.ylabel('Freq')
plt.title("Histogram of Freq")

# Plot histogram using matplotlib bar().
indexes = np.arange(len(aaSeq))
plt.bar(indexes, locFreq)
plt.xticks(indexes, aaSeqList)
plt.show()
    





# UCSCiGEM2017
