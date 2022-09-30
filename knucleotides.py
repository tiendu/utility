import argparse
import copy
import re
import random
from itertools import groupby

parser = argparse.ArgumentParser(description='k-nucleotide clustering.')
parser.add_argument('-i', type = str, required = True)
parser.add_argument('-k', type = int, required = True)
parser.add_argument('-p', type = int, required = True)
args = parser.parse_args()

def fasta_iter(filename):
    fin = open(filename, 'rb')
    faiter = (x[1] for x in groupby(fin, lambda line: str(line, 'utf-8')[0] == ">"))
    for header in faiter:
        headerstr = str(header.__next__(), 'utf-8')
        long_name = headerstr.strip().replace('>', '')
        name = long_name.split()[0]
        seq = "".join(str(s, 'utf-8').strip() for s in faiter.__next__())
        yield name, seq
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
class IdSequence:
    def __init__(self, id, sequence):
        self.id = id
        self.sequence = sequence
    def base_freq(self):
        bases = ["A", "T", "G", "C"]
        bases = {i:0 for i in bases}
        for i in range(0, len(self.sequence), 1):
            bases[self.sequence[i:i + 1:1]] += 1
        for i in bases.keys():
            bases[i] = bases[i] / len(self.sequence)
        return bases
    def kmer_freq(self, k):
        kmers = generate_kmer(k)
        kmers = {i:0 for i in kmers}
        for i in range(0, len(self.sequence) - 3, 1):
            kmers[self.sequence[i:i + 4:1]] += 1
        total = 0
        for i in kmers.keys():
            total += kmers[i]
        for i in kmers.keys():
            kmers[i] = kmers[i] / total 
        return kmers
    def norm_kmer_freq(self, k):
        base_freq = self.base_freq()
        kmer_freq = self.kmer_freq(k)
        norm_freq = {}
        for i in kmer_freq.keys():
            norm_freq[i] = kmer_freq[i]
            nu = [*i]
        for j in nu:
            if base_freq[j] != 0:
                norm_freq[i] = norm_freq[i] / base_freq[j]
        return norm_freq
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
def generate_kmer(k):
    bases_1 = ["A", "T", "G", "C"]
    bases_2 = bases_1.copy()
    i = 0
    while i < k - 1:
        i += 1
        temp = []
        for m in bases_1:
            for n in bases_2:
                temp.append(m + n)
        bases_2 = None
        bases_2 = temp
    return bases_2

def fix_sequence(sequence):
    temp = ""
    for i in range(0, len(sequence)):
        nu = [*sequence]
        if nu[i] == "R":
            nu[i] = random.choice(["A", "G"])
        elif nu[i] == "Y":
            nu[i] = random.choice(["C", "T"])
        elif nu[i] == "S":
            nu[i] = random.choice(["G", "C"])
        elif nu[i] == "W":
            nu[i] = random.choice(["A", "T"])
        elif nu[i] == "K":
            nu[i] = random.choice(["G", "T"])
        elif nu[i] == "M":
            nu[i] = random.choice(["A", "C"])
        elif nu[i] == "B":
            nu[i] = random.choice(["C", "G", "T"])
        elif nu[i] == "D":
            nu[i] = random.choice(["A", "G", "T"])
        elif nu[i] == "H":
            nu[i] = random.choice(["A", "C", "T"])
        elif nu[i] == "V":
            nu[i] = random.choice(["A", "C", "G"])
        elif nu[i] == "N":
            nu[i] = random.choice(["A", "T", "G", "C"])
        temp += nu[i]
    return temp
# A.................Adenine
# C.................Cytosine
# G.................Guanine
# T (or U)..........Thymine (or Uracil)
# R.................A or G
# Y.................C or T
# S.................G or C
# W.................A or T
# K.................G or T
# M.................A or C
# B.................C or G or T
# D.................A or G or T
# H.................A or C or T
# V.................A or C or G
# N.................any base
# . or -............gap

def remove_tandem_repeat(sequence, k, c):
    for i in generate_kmer(k):
        sequence = re.sub("(" + i + ")" + "{" + str(c) + ",}", "", sequence)
        return sequence
def create_reverse_complement(sequence):
    rev = ""
    for i in range(len(sequence) - 1, -1, -1):
        nu = [*sequence]
        if nu[i] == "A":
            rev += "T"
        elif nu[i] == "T":
            rev += "A"
        elif nu[i] == "G":
            rev += "C"
        elif nu[i] == "C":
            rev += "G"
    return rev

fasta = fasta_iter(args.i)
labels = []
features = []
for i in fasta:
    n = IdSequence(i[0], remove_tandem_repeat(fix_sequence(i[1]), args.k, 3))
    labels.append(n.id)
    if args.p == 1:
        palindromes = []
        for j in n.norm_kmer_freq(args.k).keys():
            if j == create_reverse_complement(j):
                palindromes.append(n.norm_kmer_freq(args.k)[j])
        features.append(palindromes)
    else:
        features.append(list(n.norm_kmer_freq(args.k).values()))
features = np.array(features)

