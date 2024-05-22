import sys
from dataclasses import dataclass
from typing import List, Tuple
from concurrent.futures import ProcessPoolExecutor, as_completed


@dataclass
class SequenceAlignment:
    sequence1: str
    sequence2: str
    match: int = 2
    mismatch: int = -1
    gap_opening: int = -2
    gap_extension: int = -1
    aligned_sequence1: str = ''
    aligned_sequence2: str = ''
    score: int = 0

    def set_params(self, match: int = None, mismatch: int = None, gap_opening: int = None, gap_extension: int = None) -> None:
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
    def align(self) -> Tuple[str, str, str, int]:
        rows = len(self.sequence1) + 1
        cols = len(self.sequence2) + 1
        matrix = [[0] * cols for _ in range(rows)]

        for i in range(rows):
            matrix[i][0] = i * self.gap_opening
        for j in range(cols):
            matrix[0][j] = j * self.gap_opening

        for i in range(1, rows):
            for j in range(1, cols):
                match_score = matrix[i - 1][j - 1] + (self.match if self.sequence1[i - 1] == self.sequence2[j - 1] else self.mismatch)
                gap_open_score1 = matrix[i - 1][j] + self.gap_opening
                gap_open_score2 = matrix[i][j - 1] + self.gap_opening
                matrix[i][j] = max(match_score, gap_open_score1, gap_open_score2)

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

def compute_pairwise_distance(pair: Tuple[str, str]) -> Tuple[int, int, int]:
    seq1, seq2, i, j = pair
    alignment = GlobalAlignment(seq1, seq2)
    alignment.set_params(match=3, mismatch=-2, gap_opening=-3, gap_extension=-2)
    _, _, _, score = alignment.align()
    return i, j, -score

def compute_pairwise_distances(sequences: List[str]) -> List[List[int]]:
    num_sequences = len(sequences)
    distances = [[0] * num_sequences for _ in range(num_sequences)]
    pairs = [(sequences[i], sequences[j], i, j) for i in range(num_sequences) for j in range(i + 1, num_sequences)]

    with ProcessPoolExecutor() as executor:
        futures = [executor.submit(compute_pairwise_distance, pair) for pair in pairs]
        for future in as_completed(futures):
            i, j, score = future.result()
            distances[i][j] = distances[j][i] = score

    return distances

def align_two(seq1: str, seq2: str) -> Tuple[str, str]:
    alignment = GlobalAlignment(seq1, seq2)
    alignment.set_params(match=3, mismatch=-2, gap_opening=-3, gap_extension=-2)
    aligned_seq1, _, aligned_seq2, _ = alignment.align()
    return aligned_seq1, aligned_seq2

def align_progressively(sequences: List[str]) -> List[str]:
    num_sequences = len(sequences)
    aligned_sequences = sequences.copy()
    distance_matrix = compute_pairwise_distances(sequences)
    unaligned = set(range(num_sequences))

    while len(unaligned) > 1:
        min_dist = float('inf')
        closest_pair = (0, 1)
        for i in unaligned:
            for j in unaligned:
                if i != j and distance_matrix[i][j] < min_dist:
                    min_dist = distance_matrix[i][j]
                    closest_pair = (i, j)

        i, j = closest_pair
        aligned_seq1, aligned_seq2 = align_two(aligned_sequences[i], aligned_sequences[j])
        aligned_sequences[i] = aligned_seq1
        aligned_sequences[j] = aligned_seq2
        for k in unaligned:
            if k != i and k != j:
                seq1 = aligned_sequences[i]
                seq2 = aligned_sequences[k]
                _, _, _, score = GlobalAlignment(seq1, seq2).align()
                distance_matrix[i][k] = distance_matrix[k][i] = -score

                seq1 = aligned_sequences[j]
                seq2 = aligned_sequences[k]
                _, _, _, score = GlobalAlignment(seq1, seq2).align()
                distance_matrix[j][k] = distance_matrix[k][j] = -score

        unaligned.remove(j)
    max_len = max(len(seq) for seq in aligned_sequences)
    for i in range(num_sequences):
        aligned_sequences[i] = aligned_sequences[i].ljust(max_len, '-')

    return aligned_sequences

def main(sequences: List[str]) -> None:
    aligned_sequences = align_progressively(sequences)

    for seq in aligned_sequences:
        print(seq)

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print(f'Usage: python {sys.argv[0]} sequence1 sequence2 [sequence3 ...]')
    else:
        main(sys.argv[1:])
