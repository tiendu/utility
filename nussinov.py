from dataclasses import dataclass
from itertools import groupby
import sys

@dataclass(frozen=True)
class Seq:
    id: str
    sequence: str
    def transcribe(self):
        '''Convert DNA sequence to RNA sequence.'''
        return Seq(self.id, self.sequence.replace('T', 'U'))
    def reverse_complement(self):
        '''Calculate the reverse complement of the DNA sequence.'''
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        rev_comp_seq = ''.join(complement.get(base, base) for base in reversed(self.sequence))
        return Seq(self.id, rev_comp_seq)

def can_pair(base1, base2):
    '''Check if two RNA bases can form a pair.'''
    pairs = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G'}
    return pairs.get(base1) == base2

def nussinov_half_matrix(seq):
    '''Predict RNA secondary structure using the Nussinov algorithm with half-matrix representation.'''
    n = len(seq)
    M = [[0] * n for _ in range(n)]

    # Fill the DP table
    for k in range(1, n):
        for i in range(n - k):
            j = i + k
            if can_pair(seq[i], seq[j]):
                M[i][j] = M[i+1][j-1] + 1
            for l in range(i, j):
                M[i][j] = max(M[i][j], M[i][l] + M[l+1][j])
            M[i][j] = max(M[i][j], M[i+1][j], M[i][j-1])

    def traceback(i, j, structure):
        '''Enhanced traceback to find optimal RNA secondary structure.'''
        tb_stack = [(i, j)]
        last_direction = None

        while tb_stack:
            i, j = tb_stack.pop()
            if i >= j:
                continue
            candidate_stack = []
            if M[i][j] == M[i+1][j-1] + (1 if can_pair(seq[i], seq[j]) else 0):
                candidate_stack.append(('DIAGONAL', i+1, j-1))
            if M[i][j] == M[i][j-1]:
                candidate_stack.append(('LEFT', i, j-1))
            if M[i][j] == M[i+1][j]:
                candidate_stack.append(('DOWN', i+1, j))

            if not candidate_stack:
                last_direction = None
                for k in range(i+1, j):
                    if M[i][j] == M[i][k] + M[k+1][j]:
                        tb_stack.append((k+1, j))
                        tb_stack.append((i, k))
                        break
            else:
                pushed = False
                while candidate_stack:
                    d, m, n = candidate_stack.pop()
                    if last_direction is None or d == last_direction:
                        tb_stack.append((m, n))
                        if can_pair(seq[i], seq[j]):
                            structure[i] = '('
                            structure[j] = ')'
                        pushed = True
                        break
                if not pushed:
                    tb_stack.append((m, n))
                last_direction = d

    structure = ['.' for _ in range(n)]
    traceback(0, n-1, structure)
    return ''.join(structure)

def read_fasta(filename):
    '''Read sequences from a FASTA file and return a list of Seq dataclass instances.'''
    sequences = []
    try:
        with open(filename, 'rt') as fin:
            faiter = (x[1] for x in groupby(fin, lambda line: line[0] == '>'))
            for header in faiter:
                headerstr = next(header).strip()
                name = headerstr[1:]  # Removing ">" character
                seq = ''.join(s.strip().upper() for s in next(faiter))
                sequences.append(Seq(name, seq))
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
        sys.exit(1)
    return sequences

def main():
    '''Main function to process FASTA file and predict RNA structure.'''
    if len(sys.argv) != 2:
        print(f'Usage: python {sys.argv[0]} <fasta_file>')
        sys.exit(1)

    fasta_file = sys.argv[1]
    sequences = read_fasta(fasta_file)

    for sequence in sequences:
        rna_sequence = sequence.transcribe()
        structure = nussinov_half_matrix(rna_sequence.sequence)
        print(f'>{sequence.id}')
        print(rna_sequence)
        print(structure)

if __name__ == '__main__':
    main()
