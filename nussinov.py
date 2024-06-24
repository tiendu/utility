from dataclasses import dataclass
from itertools import groupby
import sys
import re
import numpy as np

'''
Nussinov algorithm is a fundamental tool in computational biology used
 to predict RNA secondary structures by identifying base-pairing inter-
actions within RNA sequences. 
The algorithm detects regions of high structural stability or folding 
propensity in RNA molecules. 
This script implements the Nussinov algorithm to predict RNA secondary
 structures from input RNA sequences, providing insights into potentia-
l functional regions within the RNA construct.
'''

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
    def extract_subsequence(self, start, end):
        '''Extract a subsequence given a start and end position.'''
        if end > len(self.sequence):
            return Seq(self.id + f' {start+1}..{len(self.sequence)}', self.sequence[start:end])
        return Seq(self.id + f' {start+1}..{end+1}', self.sequence[start:end])
    def locate_subsequence(self, subsequence):
        '''Locate the start and end positions of a subsequence within the original sequence.'''
        match = re.search(subsequence, self.sequence)
        if match:
            return match.start(), match.end()
        else:
            return -1, -1

def can_pair(base1, base2):
    '''Check if two bases can form a pair.'''
    pairs = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G'}
    return pairs.get(base1) == base2

def nussinov_half_matrix(seq):
    '''Predict RNA secondary structure using the Nussinov algorithm with half-matrix representation.'''
    n = len(seq)
    M = np.zeros((n, n), dtype=int)

    for k in range(1, n):
        for i in range(n - k):
            j = i + k
            M[i, j] = M[i+1, j]
            for l in range(i, j):
                M[i, j] = max(M[i, j], M[i, l] + M[l+1, j])
            if can_pair(seq[i], seq[j]):
                M[i, j] = max(M[i, j], M[i+1, j-1] + 1)

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

def sum_consecutive_pairs(structure):
    '''Sum the lengths of all consecutive pairs longer than 6 in the structure string.'''
    total_sum = 0
    current_length = 0
    stack = []
    threshold = 6  # Minimum length threshold for consecutive pairs

    for char in structure:
        if char == '(':
            stack.append(char)
        elif char == ')' and stack:
            stack.pop()
            current_length += 1
        else:
            if current_length > threshold:
                total_sum += current_length
            current_length = 0

    # Check if last consecutive pairs were longer than threshold
    if current_length > threshold:
        total_sum += current_length

    return total_sum

def read_fasta(filename):
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

def is_stem_loop(structure):
    '''Check if a structure contains a stem-loop pattern.'''
    in_stem = False
    stem_length = 0
    loop_length = 0
    min_stem_length = 12  # Minimum length for a stem
    min_loop_length = 0  # Minimum length for a loop

    for char in structure:
        if char == '(':
            if not in_stem:
                in_stem = True
                stem_length = 1
            else:
                stem_length += 1
        elif char == '.':
            if in_stem:
                if stem_length >= min_stem_length and loop_length >= min_loop_length:
                    return True
                in_stem = False
                stem_length = 0
                loop_length = 0
            else:
                loop_length += 1
        else:  # Reset on non-stem, non-loop characters
            in_stem = False
            stem_length = 0
            loop_length = 0
    
    # Check at the end of the structure
    if in_stem and stem_length >= min_stem_length and loop_length >= min_loop_length:
        return True

    return False

def main():
    '''Main function to process FASTA file and predict RNA structure.'''
    if len(sys.argv) != 2:
        print(f'Usage: python {sys.argv[0]} <fasta_file>')
        return

    fasta_file = sys.argv[1]
    sequences = read_fasta(fasta_file)
    chunk_size = 200

    for sequence in sequences:
        for i in range(0, len(sequence.sequence), chunk_size):
            subsequence = sequence.extract_subsequence(i, chunk_size+i)
            
            forward_transcript = subsequence.transcribe().sequence
            forward_structure = nussinov_half_matrix(forward_transcript)
            forward_max_pairs = sum_consecutive_pairs(forward_structure)

            reverse_transcript = subsequence.reverse_complement().transcribe().sequence
            reverse_structure = nussinov_half_matrix(reverse_transcript)
            reverse_max_pairs = sum_consecutive_pairs(reverse_structure)

            if is_stem_loop(forward_structure):
                print(f'{sequence.id}\tFW\t{i+1}..{chunk_size+i+1}')
                print(f'Forward structure: {forward_structure}')
                print(f'Forward max pairs: {forward_max_pairs}')
                print(f'{"#" * 20}')
            elif is_stem_loop(reverse_structure):
                print(f'{sequence.id}\tRV\t{i+1}..{chunk_size+i+1}')
                print(f'Reverse structure: {reverse_structure}')
                print(f'Reverse max pairs: {reverse_max_pairs}')
                print(f'{"#" * 20}')

if __name__ == '__main__':
    main()
