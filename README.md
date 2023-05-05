# UCSCiGEM2017
2017 UC Santa Cruz iGEM
This repository contains all the software and scripts for our projects and webpages.


Sure, here's an example README.md file for the Python code you provided:

# DNA to Amino Acid Translation and Frequency Analysis

This Python script reads in a DNA sequence from a file, translates it into an amino acid sequence, and analyzes the frequency of each codon and amino acid. It also generates a histogram plot of the codon frequency vs. the amino acid sequence.

## Dependencies

This script requires the following Python packages to be installed:

- `numpy`
- `matplotlib`

## Usage

1. Place the `Ecoli.txt` and `eColiCodon.txt` files in the same directory as the `dna-to-aa.py` script.
2. Open a terminal or command prompt in the directory containing the script and input files.
3. Run the script by typing `python dna-to-aa.py` and pressing Enter.

## Outputs

The script will output the following:

- The translated amino acid sequence, with any non-coding regions (i.e. "STOP" or "None") removed.
- The frequency of each codon in the DNA sequence.
- The frequency of each amino acid in the translated sequence.
- The number of unique amino acids in the translated sequence.
- The total number of amino acids in the translated sequence.
- The number of old codons used (from the `eColiCodon.txt` file) and the number of new codons used.

The script will also generate a histogram plot of the codon frequency vs. the amino acid sequence, saved as `freq-histogram.png`.

## Example Output

```
ARFAGELVGQFPYAVVTIPAFFGIVVNFLVSGALNMPFMVLMAVGVGAFGGFVTRPLAGLNNKIFVSYLFEGQFTL
Counter({'GGA': 223, 'AGA': 220, 'GGT': 217, 'CGT': 214, 'TTC': 209, 'TTG': 200, 'CAG': 197, 'GAG': 191, 'CGC': 191, 'GCG': 190, 'GTT': 188, 'CTG': 184, 'GTG': 182, 'ACA': 181, 'AAC': 180, 'AAA': 178, 'ATG': 174, 'AGG': 171, 'AAG': 167, 'TCT':...
Counter({'L': 386, 'A': 351, 'G': 305, 'V': 296, 'I': 272, 'S': 269, 'F': 228, 'P': 226, 'T': 209, 'Y': 190, 'N': 183, 'D': 173, 'K': 165, 'E': 158, 'Q': 154, 'R': 145, 'C': 108, 'H': 82, 'STOP': 41})
19
2453
# of old codons used:64
# of new codons used:61
```
