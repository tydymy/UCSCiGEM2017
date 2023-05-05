#include <iostream>
#include <fstream>
#include <unordered_map>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;

unordered_map<string, float> codTab;
unordered_map<string, string> rna2aa = {{"UUU", "F"}, {"UUC", "F"}, {"UUA", "L"}, {"UUG", "L"}, {"UCU", "S"}, {"UCC", "S"}, {"UCA", "S"}, {"UCG", "S"}, {"UAU", "Y"}, {"UAC", "Y"}, {"UAA", "STOP"}, {"UAG", "STOP"}, {"UGU", "C"}, {"UGC", "C"}, {"UGA", "STOP"}, {"UGG", "W"}, {"CUU", "L"}, {"CUC", "L"}, {"CUA", "L"}, {"CUG", "L"}, {"CCU", "P"}, {"CCC", "P"}, {"CCA", "P"}, {"CCG", "P"}, {"CAU", "H"}, {"CAC", "H"}, {"CAA", "Q"}, {"CAG", "Q"}, {"CGU", "R"}, {"CGC", "R"}, {"CGA", "R"}, {"CGG", "R"}, {"AUU", "I"}, {"AUC", "I"}, {"AUA", "I"}, {"AUG", "M"}, {"ACU", "T"}, {"ACC", "T"}, {"ACA", "T"}, {"ACG", "T"}, {"AAU", "N"}, {"AAC", "N"}, {"AAA", "K"}, {"AAG", "K"}, {"AGU", "S"}, {"AGC", "S"}, {"AGA", "R"}, {"AGG", "R"}, {"GUU", "V"}, {"GUC", "V"}, {"GUA", "V"}, {"GUG", "V"}, {"GCU", "A"}, {"GCC", "A"}, {"GCA", "A"}, {"GCG", "A"}, {"GAU", "D"}, {"GAC", "D"}, {"GAA", "E"}, {"GAG", "E"}, {"GGU", "G"}, {"GGC", "G"}, {"GGA", "G"}, {"GGG", "G"}};

void parseCodonTable(string fileName) {
    ifstream file(fileName);
    if (!file) {
        cout << "Error opening codon table file!" << endl;
        exit(1);
    }
    string line;
    while (getline(file, line)) {
        size_t start = 0;
        while (start < line.length()) {
            size_t end = line.find(">", start);
            string codon = line.substr(start, end - start);
            start = end + 1;
            end = line.find("(", start);
            string freqStr = line.substr(start, end - start - 1);
            float freq = stof(freqStr);
            codTab[codon] = freq;
            start = end + 1;
        }
    }
}

string translateToRNA(string dna) {
    string rna = "";
    for (char c : dna) {
        if (c == 'T') {
            rna += 'U';
        }
        else {
            rna += c;
        }
    }
    return rna;
}

string removeLowerCase(string str) {
    string result = "";
    for (char c : str) {
        if (isupper(c)) {
            result += c;
        }
    }
    return result;
}

void
