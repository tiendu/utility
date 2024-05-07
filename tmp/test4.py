import sys
from dataclasses import dataclass
from typing import List
from itertools import groupby


@dataclass(frozen=True)
class Seq:
    id: str
    sequence: str
    def __hash__(self):
        return hash((self.id, self.sequence))

@dataclass
class SequenceAlignment:
    sequence1: str  # First input sequence
    sequence2: str  # Second input sequence
    match: int = 2  # Score for a match
    mismatch: int = -1  # Penalty for a mismatch
    gap_opening: int = -2  # Penalty for opening a gap
    gap_extension: int = -1  # Penalty for extending a gap
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
    def align(self) -> int:
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

        self.score = matrix[rows - 1][cols - 1]
        self.alignment_length = rows + cols - 2
        return self.score


def pairwise_comparison(sequences: List[str]) -> List[List[int]]:
    similarity_scores = []
    for i in range(len(sequences)):
        similarity_scores.append([0] * len(sequences))

    for i in range(len(sequences)):
        for j in range(i + 1, len(sequences)):
            global_alignment = GlobalAlignment(sequences[i].sequence, sequences[j].sequence)
            global_alignment.set_params(match=3, mismatch=-2, gap_opening=-3, gap_extension=-2)
            score = global_alignment.align()
            similarity_scores[i][j] = similarity_scores[j][i] = score / global_alignment.alignment_length

    return similarity_scores

def read_sequences_from_fasta(file_path: str) -> List[Seq]:
    sequences = []
    with open(file_path, 'rt') as fin:
        faiter = (x[1] for x in groupby(fin, lambda line: line[0] == ">"))
        for header in faiter:
            headerstr = next(header).strip()
            name = headerstr[1:]  # Removing ">" character
            seq = "".join(s.strip().upper() for s in next(faiter))
            sequences.append(Seq(name, seq))
    return sequences

def main(file_path: str) -> None:
    sequences = read_sequences_from_fasta(file_path)
    similarity_scores = pairwise_comparison(sequences)

    print("Pairwise Similarity Scores:")
    for i in range(len(sequences)):
        for j in range(i + 1, len(sequences)):
            similarity_score = similarity_scores[i][j]
            print(f"{sequences[i].id}\t{sequences[j].id}\t{similarity_score:.2f}")

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print(f'Usage: python {sys.argv[0]} file.fasta')
    else:
        main(sys.argv[1])
