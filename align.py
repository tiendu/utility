import re
import sys
from dataclasses import dataclass
from typing import List, Tuple

@dataclass
class SequenceAlignment:
    sequence1: str  # First input sequence
    sequence2: str  # Second input sequence
    match: int = 2  # Score for a match
    mismatch: int = -1  # Penalty for a mismatch
    gap_opening: int = -2  # Penalty for opening a gap
    gap_extension: int = -1  # Penalty for extending a gap
    aligned_sequence1: str = ''  # Aligned sequence 1
    aligned_sequence2: str = ''  # Aligned sequence 2
    score: int = 0  # Alignment score

    def set_params(self, 
                   match: int = None, 
                   mismatch: int = None, 
                   gap_opening: int = None, 
                   gap_extension: int = None
                   ) -> None:
        """Set custom parameters for the alignment."""
        if match is not None:
            self.match = match
        if mismatch is not None:
            self.mismatch = mismatch
        if gap_opening is not None:
            self.gap_opening = gap_opening
        if gap_extension is not None:
            self.gap_extension = gap_extension

@dataclass
class GlobalAlignment(SequenceAlignment):
    """Perform global sequence alignment."""
    def align(self) -> Tuple[str, str, str, int]:
        # Initialize the matrix
        rows = len(self.sequence1) + 1
        cols = len(self.sequence2) + 1
        matrix = [[0] * cols for _ in range(rows)]

        # Fill the first row and column with gap penalties
        for i in range(rows):
            matrix[i][0] = i * self.gap_opening
        for j in range(cols):
            matrix[0][j] = j * self.gap_opening

        # Fill the rest of the matrix
        for i in range(1, rows):
            for j in range(1, cols):
                match_score = matrix[i - 1][j - 1] + (self.match if self.sequence1[i - 1] == self.sequence2[j - 1] else self.mismatch)
                gap_open_score = max(matrix[i - 1][j] + self.gap_opening, matrix[i][j - 1] + self.gap_opening)
                gap_extend_score = matrix[i - 1][j - 1] + self.gap_extension
                matrix[i][j] = max(match_score, gap_open_score, gap_extend_score)

        # Traceback to find the alignment
        alignment_symbols = ''
        i, j = rows - 1, cols - 1
        while i > 0 and j > 0:
            if matrix[i][j] == matrix[i - 1][j - 1] + (self.match if self.sequence1[i - 1] == self.sequence2[j - 1] else self.mismatch):
                self.aligned_sequence1 = self.sequence1[i - 1] + self.aligned_sequence1
                self.aligned_sequence2 = self.sequence2[j - 1] + self.aligned_sequence2
                alignment_symbols = '|' + alignment_symbols if self.sequence1[i - 1] == self.sequence2[j - 1] else ' ' + alignment_symbols
                i -= 1
                j -= 1
            elif matrix[i][j] == matrix[i - 1][j] + self.gap_opening:
                self.aligned_sequence1 = self.sequence1[i - 1] + self.aligned_sequence1
                self.aligned_sequence2 = '-' + self.aligned_sequence2
                alignment_symbols = ' ' + alignment_symbols
                i -= 1
            elif matrix[i][j] == matrix[i][j - 1] + self.gap_opening:
                self.aligned_sequence1 = '-' + self.aligned_sequence1
                self.aligned_sequence2 = self.sequence2[j - 1] + self.aligned_sequence2
                alignment_symbols = ' ' + alignment_symbols
                j -= 1
            elif matrix[i][j] == matrix[i - 1][j - 1] + self.gap_extension:
                self.aligned_sequence1 = self.sequence1[i - 1] + self.aligned_sequence1
                self.aligned_sequence2 = self.sequence2[j - 1] + self.aligned_sequence2
                alignment_symbols = '|' + alignment_symbols if self.sequence1[i - 1] == self.sequence2[j - 1] else ' ' + alignment_symbols
                i -= 1
                j -= 1

        # Fill in the rest of the alignment if one sequence is longer than the other
        while i > 0:
            self.aligned_sequence1 = self.sequence1[i - 1] + self.aligned_sequence1
            self.aligned_sequence2 = '-' + self.aligned_sequence2
            alignment_symbols = ' ' + alignment_symbols
            i -= 1
        while j > 0:
            self.aligned_sequence1 = '-' + self.aligned_sequence1
            self.aligned_sequence2 = self.sequence2[j - 1] + self.aligned_sequence2
            alignment_symbols = ' ' + alignment_symbols
            j -= 1

        self.score = matrix[rows - 1][cols - 1]

        return self.aligned_sequence1, alignment_symbols, self.aligned_sequence2, self.score
    
@dataclass
class LocalAlignment(SequenceAlignment):
    """Perform local sequence alignment."""
    def align(self) -> Tuple[str, str, str, int]:
        rows = len(self.sequence1) + 1
        cols = len(self.sequence2) + 1
        matrix = [[0] * cols for _ in range(rows)]
        
        max_i = 0
        max_j = 0

        # Fill the matrix and track maximum score position
        for i in range(1, rows):
            for j in range(1, cols):
                match_score = matrix[i - 1][j - 1] + (self.match if self.sequence1[i - 1] == self.sequence2[j - 1] else self.mismatch)
                gap_open_score = max(matrix[i - 1][j] + self.gap_opening, matrix[i][j - 1] + self.gap_opening, 0)
                gap_extend_score = matrix[i - 1][j - 1] + self.gap_extension
                matrix[i][j] = max(match_score, gap_open_score, gap_extend_score, 0)

                if matrix[i][j] > self.score:
                    self.score = matrix[i][j]
                    max_i = i
                    max_j = j
        
        alignment_symbols = ''
        i, j = max_i, max_j

        # Traceback to find the local alignment
        while i > 0 and j > 0 and matrix[i][j] != 0:
            if matrix[i][j] == matrix[i - 1][j - 1] + (self.match if self.sequence1[i - 1] == self.sequence2[j - 1] else self.mismatch):
                self.aligned_sequence1 = self.sequence1[i - 1] + self.aligned_sequence1
                self.aligned_sequence2 = self.sequence2[j - 1] + self.aligned_sequence2
                alignment_symbols = '|' + alignment_symbols if self.sequence1[i - 1] == self.sequence2[j - 1] else ' ' + alignment_symbols
                i -= 1
                j -= 1
            elif matrix[i][j] == matrix[i - 1][j] + self.gap_opening:
                self.aligned_sequence1 = self.sequence1[i - 1] + self.aligned_sequence1
                self.aligned_sequence2 = '-' + self.aligned_sequence2
                alignment_symbols = ' ' + alignment_symbols
                i -= 1
            elif matrix[i][j] == matrix[i][j - 1] + self.gap_opening:
                self.aligned_sequence1 = '-' + self.aligned_sequence1
                self.aligned_sequence2 = self.sequence2[j - 1] + self.aligned_sequence2
                alignment_symbols = ' ' + alignment_symbols
                j -= 1
            elif matrix[i][j] == matrix[i - 1][j - 1] + self.gap_extension:
                self.aligned_sequence1 = self.sequence1[i - 1] + self.aligned_sequence1
                self.aligned_sequence2 = self.sequence2[j - 1] + self.aligned_sequence2
                alignment_symbols = '|' + alignment_symbols if self.sequence1[i - 1] == self.sequence2[j - 1] else ' ' + alignment_symbols
                i -= 1
                j -= 1

        return self.aligned_sequence1, alignment_symbols, self.aligned_sequence2, self.score

    def find_local_alignment_location(self) -> List[int]:
        """Find the positions of local alignment in the original sequences."""
        # Remove gap characters from the aligned sequences
        aligned_sequence1_without_gaps = self.aligned_sequence1.replace('-', '')
        aligned_sequence2_without_gaps = self.aligned_sequence2.replace('-', '')

        # Find the starting position of the aligned sequence in the original sequence
        match1 = re.search(aligned_sequence1_without_gaps, self.sequence1)
        match2 = re.search(aligned_sequence2_without_gaps, self.sequence2)
        
        if match1 and match2:
            start1 = match1.start()
            start2 = match2.start()
            return [start1 + 1, start2 + 1]

def print_alignment(seq1: str, symbols: str, seq2: str, width: int=60) -> None:
    '''Print sequences and alignment symbols in chunks with annotations.'''
    max_len = max(len(seq1), len(seq2))

    for i in range(0, max_len, width):
        before_chunk = seq1[i:i+width] if i < len(seq1) else ''
        after_chunk = seq2[i:i+width] if i < len(seq2) else ''
        annotation_chunk = symbols[i:i+width] if i < len(symbols) else ''

        print(f'{i+1:4d} {before_chunk}')
        print(f'     {annotation_chunk}')
        print(f'{i+1:4d} {after_chunk}\n')

def main(sequence1: str, sequence2: str) -> None:
    # Perform global alignment
    global_alignment = GlobalAlignment(sequence1, sequence2)
    global_alignment.set_params(match=3, mismatch=-2, gap_opening=-3, gap_extension=-2)
    global_aligned_sequence1, global_alignment_symbols, global_aligned_sequence2, global_score = global_alignment.align()

    # Perform local alignment
    local_alignment = LocalAlignment(sequence1, sequence2)
    local_alignment.set_params(match=3, mismatch=-2, gap_opening=-3, gap_extension=-2)
    local_aligned_sequence1, local_alignment_symbols, local_aligned_sequence2, local_score = local_alignment.align()
    local_alignment_position1, local_alignment_position2 = local_alignment.find_local_alignment_location()

    # Print global alignment results
    print('Global Alignment\n')
    print(global_aligned_sequence1)
    print(global_alignment_symbols)
    print(global_aligned_sequence2)
    print('Alignment Score:', global_score)

    print('=' * 40)  # Separator

    # Print local alignment results
    print('Local Alignment\n')
    print(local_aligned_sequence1)
    print(local_alignment_symbols)
    print(local_aligned_sequence2)
    print('Alignment Score:', local_score)
    print('Aligned Sequence 1 Position:', local_alignment_position1)
    print('Aligned Sequence 2 Position:', local_alignment_position2)

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print(f'Usage: python {sys.argv[0]} sequence1 sequence2')
    else:
        main(sys.argv[1], sys.argv[2])
