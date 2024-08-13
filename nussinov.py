from dataclasses import dataclass
from itertools import groupby
import sys
import re
import concurrent.futures
from typing import List, Tuple

'''
Nussinov algorithm is a fundamental tool in computational biology used
to predict RNA secondary structures by identifying base-pairing interactions
within RNA sequences.
The algorithm detects regions of high structural stability or folding propensity
in RNA molecules.
This script implements the Nussinov algorithm to predict RNA secondary structures
from input RNA sequences, providing insights into potential functional regions
within the RNA construct.
'''

@dataclass(frozen=True)
class Seq:
    id: str
    sequence: str

    def transcribe(self) -> 'Seq':
        '''Convert DNA sequence to RNA sequence.'''
        return Seq(self.id, self.sequence.replace('T', 'U'))

    def reverse_complement(self) -> 'Seq':
        '''Calculate the reverse complement of the DNA sequence.'''
        complement = str.maketrans('ATGC', 'TACG')
        return Seq(self.id, self.sequence[::-1].translate(complement))

    def extract_subsequence(self, start: int, end: int) -> 'Seq':
        '''Extract a subsequence given a start and end position.'''
        if end > len(self.sequence):
            return Seq(self.id + f' {start+1}..{len(self.sequence)}', self.sequence[start:end])
        return Seq(self.id + f' {start+1}..{end+1}', self.sequence[start:end])

    def locate_subsequence(self, subsequence: str) -> Tuple[int, int]:
        '''Locate the start and end positions of a subsequence within the original sequence.'''
        match = re.search(subsequence, self.sequence)
        if match:
            return match.start(), match.end()
        else:
            return -1, -1

    def length(self) -> int:
        return len(self.sequence)

def can_pair(base1: str, base2: str) -> bool:
    '''Check if two bases can form a pair.'''
    pairs = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G'}
    return pairs.get(base1) == base2

def nussinov_half_matrix(transcript: str) -> str:
    '''Predict RNA secondary structure using the Nussinov algorithm with half-matrix representation.'''
    n = len(transcript)
    M = [[0] * n for _ in range(n)]

    # Fill the DP table
    for k in range(1, n):
        for i in range(n - k):
            j = i + k
            if can_pair(transcript[i], transcript[j]):
                M[i][j] = M[i+1][j-1] + 1
            for l in range(i, j):
                M[i][j] = max(M[i][j], M[i][l] + M[l+1][j])
            M[i][j] = max(M[i][j], M[i+1][j], M[i][j-1])

    def traceback(i: int, j: int, structure: List[str]):
        '''Enhanced traceback to find optimal RNA secondary structure.'''
        tb_stack = [(i, j)]
        last_direction = None

        while tb_stack:
            i, j = tb_stack.pop()
            if i >= j:
                continue
            candidate_stack = []
            if M[i][j] == M[i+1][j-1] + (1 if can_pair(transcript[i], transcript[j]) else 0):
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
                        if can_pair(transcript[i], transcript[j]):
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

def read_fasta(filename: str) -> List[Seq]:
    '''Read sequences from a FASTA file.'''
    sequences = []
    with open(filename, 'rt') as fin:
        faiter = (x[1] for x in groupby(fin, lambda line: line[0] == '>'))
        for header in faiter:
            headerstr = next(header).strip()
            name = headerstr[1:]  # Removing ">" character
            seq = ''.join(s.strip().upper() for s in next(faiter))
            sequences.append(Seq(name, seq))
    return sequences

def nussinov(sequence: Seq) -> List[Tuple[str, str, str, str]]:
    forward_transcript = sequence.transcribe()
    forward_structure = nussinov_half_matrix(forward_transcript.sequence)

    reverse_transcript = sequence.reverse_complement().transcribe()
    reverse_structure = nussinov_half_matrix(reverse_transcript.sequence)

    results = []

    results.append((sequence.id, forward_transcript.sequence, 'FW', forward_structure))
    results.append((sequence.id, reverse_transcript.sequence, 'RV', reverse_structure))

    return results

def main():
    '''Main function to process FASTA file and predict RNA structure.'''
    if len(sys.argv) != 3:
        print(f'Usage: python {sys.argv[0]} <input_file> <output_file>')
        return

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    sequences = read_fasta(input_file)

    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures = [executor.submit(nussinov, sequence) for sequence in sequences]

        with open(output_file, 'w') as fout:
            for future in concurrent.futures.as_completed(futures):
                for result in future.result():
                    seqid, seq, direction, structure = result
                    fout.write(f'>{seqid}|{direction}\n')
                    fout.write(f'{seq}\n')
                    fout.write(f'{structure}\n')

if __name__ == '__main__':
    main()
