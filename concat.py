import argparse
from itertools import groupby

class IdSequence:
    def __init__(self, id, sequence):
        self.id = id
        self.sequence = sequence

def fasta_iter(filename):
    fin = open(filename, 'rb')
    faiter = (x[1] for x in groupby(fin, lambda line: str(line, 'utf-8')[0] == ">"))
    for header in faiter:
        headerstr = str(header.__next__(), 'utf-8')
        long_name = headerstr.strip().replace('>', '')
        name = long_name.split()[0]
        seq = "".join(str(s, 'utf-8').strip() for s in faiter.__next__())
        yield name, seq

def overlap_concat(s1, s2):
    l = min(len(s1), len(s2))
    for i in range(l, 0, -1):
        if s1.endswith(s2[:i]):
            return s1 + s2[i:]
    return ""
        
parser = argparse.ArgumentParser()
parser.add_argument('-i', type = str, required = True)
args = parser.parse_args()

sequences = []
for i in fasta_iter(args.i):
    n = IdSequence(i[0], i[1])
    sequences.append(n)

sequences = sorted(sequences, key = lambda i: len(i.sequence), reverse = True)

for i in range(len(sequences)):
    for j in range(i + 1, len(sequences)):
        concat = overlap_concat(sequences[i].sequence, sequences[j].sequence)
        if (concat != "" and len(concat) < (len(sequences[i].sequence) + len(sequences[j].sequence))):
            print(f'>{sequences[i].id}_{sequences[j].id}')
            print(f'{concat}')
