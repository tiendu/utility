import sys
from dataclasses import dataclass
from typing import Tuple

@dataclass
class SequenceAlignment:
    sequence1: str
    sequence2: str
    match: int = 2
    mismatch: int = -1
    gap_opening: int = -2
    gap_extension: int = -1

    def set_params(self, 
                   match: int = None, 
                   mismatch: int = None, 
                   gap_opening: int = None, 
                   gap_extension: int = None
                   ) -> None:
        if match is not None:
            self.match = match
        if mismatch is not None:
            self.mismatch = mismatch
        if gap_opening is not None:
            self.gap_opening = gap_opening
        if gap_extension is not None:
            self.gap_extension = gap_extension

    def align(self) -> Tuple[str, str, str, int]:
        # Initialize the matrix
        rows = len(self.sequence1) + 1
        cols = len(self.sequence2) + 1
        matrix = [[0] * cols for _ in range(rows)]

        for i in range(rows):
            matrix[i][0] = i * self.gap_opening
        for j in range(cols):
            matrix[0][j] = j * self.gap_opening

        # Fill the matrix
        for i in range(1, rows):
            for j in range(1, cols):
                match_score = matrix[i - 1][j - 1] + (self.match if self.sequence1[i - 1] == self.sequence2[j - 1] else self.mismatch)
                gap_open_score = max(matrix[i - 1][j] + self.gap_opening, matrix[i][j - 1] + self.gap_opening)
                gap_extend_score = matrix[i - 1][j - 1] + self.gap_extension
                matrix[i][j] = max(match_score, gap_open_score, gap_extend_score)

        # Traceback to find the alignment
        alignment1 = ''
        alignment2 = ''
        alignment_symbols = ''
        i, j = rows - 1, cols - 1
        while i > 0 and j > 0:
            if matrix[i][j] == matrix[i - 1][j - 1] + (self.match if self.sequence1[i - 1] == self.sequence2[j - 1] else self.mismatch):
                alignment1 = self.sequence1[i - 1] + alignment1
                alignment2 = self.sequence2[j - 1] + alignment2
                alignment_symbols = '|' + alignment_symbols if self.sequence1[i - 1] == self.sequence2[j - 1] else ' ' + alignment_symbols
                i -= 1
                j -= 1
            elif matrix[i][j] == matrix[i - 1][j] + self.gap_opening:
                alignment1 = self.sequence1[i - 1] + alignment1
                alignment2 = '-' + alignment2
                alignment_symbols = ' ' + alignment_symbols
                i -= 1
            elif matrix[i][j] == matrix[i][j - 1] + self.gap_opening:
                alignment1 = '-' + alignment1
                alignment2 = self.sequence2[j - 1] + alignment2
                alignment_symbols = ' ' + alignment_symbols
                j -= 1
            elif matrix[i][j] == matrix[i - 1][j - 1] + self.gap_extension:
                alignment1 = self.sequence1[i - 1] + alignment1
                alignment2 = self.sequence2[j - 1] + alignment2
                alignment_symbols = '|' + alignment_symbols if self.sequence1[i - 1] == self.sequence2[j - 1] else ' ' + alignment_symbols
                i -= 1
                j -= 1

        # Fill in the rest of the alignment if one sequence is longer than the other
        while i > 0:
            alignment1 = self.sequence1[i - 1] + alignment1
            alignment2 = '-' + alignment2
            alignment_symbols = ' ' + alignment_symbols
            i -= 1
        while j > 0:
            alignment1 = '-' + alignment1
            alignment2 = self.sequence2[j - 1] + alignment2
            alignment_symbols = ' ' + alignment_symbols
            j -= 1

        score = matrix[rows - 1][cols - 1]

        return alignment1, alignment_symbols, alignment2, score
    
@dataclass
class LocalAlignment(SequenceAlignment):
    def align(self) -> Tuple[str, str, str, int]:
        rows = len(self.sequence1) + 1
        cols = len(self.sequence2) + 1
        matrix = [[0] * cols for _ in range(rows)]
        
        max_score = 0
        max_i = 0
        max_j = 0

        for i in range(1, rows):
            for j in range(1, cols):
                match_score = matrix[i - 1][j - 1] + (self.match if self.sequence1[i - 1] == self.sequence2[j - 1] else self.mismatch)
                gap_open_score = max(matrix[i - 1][j] + self.gap_opening, matrix[i][j - 1] + self.gap_opening, 0)
                gap_extend_score = matrix[i - 1][j - 1] + self.gap_extension
                matrix[i][j] = max(match_score, gap_open_score, gap_extend_score, 0)

                if matrix[i][j] > max_score:
                    max_score = matrix[i][j]
                    max_i = i
                    max_j = j
        
        alignment1 = ''
        alignment2 = ''
        alignment_symbols = ''
        i, j = max_i, max_j

        while i > 0 and j > 0 and matrix[i][j] != 0:
            if matrix[i][j] == matrix[i - 1][j - 1] + (self.match if self.sequence1[i - 1] == self.sequence2[j - 1] else self.mismatch):
                alignment1 = self.sequence1[i - 1] + alignment1
                alignment2 = self.sequence2[j - 1] + alignment2
                alignment_symbols = '|' + alignment_symbols if self.sequence1[i - 1] == self.sequence2[j - 1] else ' ' + alignment_symbols
                i -= 1
                j -= 1
            elif matrix[i][j] == matrix[i - 1][j] + self.gap_opening:
                alignment1 = self.sequence1[i - 1] + alignment1
                alignment2 = '-' + alignment2
                alignment_symbols = ' ' + alignment_symbols
                i -= 1
            elif matrix[i][j] == matrix[i][j - 1] + self.gap_opening:
                alignment1 = '-' + alignment1
                alignment2 = self.sequence2[j - 1] + alignment2
                alignment_symbols = ' ' + alignment_symbols
                j -= 1
            elif matrix[i][j] == matrix[i - 1][j - 1] + self.gap_extension:
                alignment1 = self.sequence1[i - 1] + alignment1
                alignment2 = self.sequence2[j - 1] + alignment2
                alignment_symbols = '|' + alignment_symbols if self.sequence1[i - 1] == self.sequence2[j - 1] else ' ' + alignment_symbols
                i -= 1
                j -= 1

        return alignment1, alignment_symbols, alignment2, max_score

def main(sequence1: str, sequence2: str) -> None:
    global_alignment = SequenceAlignment(sequence1, sequence2)
    global_alignment.set_params(match=3, mismatch=-2, gap_opening=-3, gap_extension=-2)
    aligned_sequence1, alignment_symbols, aligned_sequence2, global_score = global_alignment.align()

    local_alignment = LocalAlignment(sequence1, sequence2)
    local_alignment.set_params(match=3, mismatch=-2, gap_opening=-3, gap_extension=-2)
    local_aligned_sequence1, local_alignment_symbols, local_aligned_sequence2, local_score = local_alignment.align()

    print('Global Alignment\n')
    print(aligned_sequence1)
    print(alignment_symbols)
    print(aligned_sequence2)
    print('Alignment Score:', global_score)
    print('=' * 40)  # Separator
    print('Local Alignment\n')
    print(local_aligned_sequence1)
    print(local_alignment_symbols)
    print(local_aligned_sequence2)
    print('Alignment Score:', local_score)

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print(f'Usage: python {sys.argv[0]} sequence1 sequence2')
    else:
        main(sys.argv[1], sys.argv[2])
