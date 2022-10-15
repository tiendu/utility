import argparse
from itertools import product, groupby

class IdSequence:
    def __init__(self, id, sequence):
        self.id = id
        self.sequence = sequence

def delineate(string):
    conversion = {
        "A": ["A"],
        "C": ["C"],
        "G": ["G"],
        "T": ["T"],
        "R": ["A","G"],
        "Y": ["C","T"],
        "S": ["G","C"],
        "W": ["A","T"],
        "K": ["G","T"],
        "M": ["A","C"],
        "B": ["C","G","T"],
        "D": ["A","G","T"],
        "H": ["A","C","T"],
        "V": ["A","C","G"],
        "N": ["A","T","G","C"],
    }
    res = []
    i = 0
    for sub in [zip(conversion.keys(), char) for char in product(*conversion.values())]:
        tmp = string
        for repls in sub:
            tmp = tmp.replace(*repls)
        res.append(tmp)
        i += 1
    return list(set(res))
    

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


def fasta_iter(filename):
    fin = open(filename, 'rb')
    faiter = (x[1] for x in groupby(fin, lambda line: str(line, 'utf-8')[0] == ">"))
    for header in faiter:
        headerstr = str(header.__next__(), 'utf-8')
        long_name = headerstr.strip().replace('>', '')
        name = long_name.split()[0]
        seq = "".join(str(s, 'utf-8').strip() for s in faiter.__next__())
        yield name, seq

parser = argparse.ArgumentParser(description='k-nucleotide clustering.')
parser.add_argument('-i', type = str, required = True)
args = parser.parse_args()

fasta = fasta_iter(args.i)
for i in fasta:
    res = delineate(i[1])
    for j in range(len(res)):
        print(f'>{i[0]}_{j}')
        print(f'{res[j]}')


