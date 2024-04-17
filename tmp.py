import sys
import re
from dataclasses import dataclass, field
from typing import List


@dataclass
class ScoringMatrix:
    match: int
    mismatch: int
    residue1: str = ''
    residue2: str = ''
    
    def get_score(self, residue1, residue2):
        if residue1 == residue2:
            return self.match
        else:
            return self.mismatch

# @dataclass
# class Identity(ScoringMatrix):
#     pass

@dataclass
class ScoreSet:
    scoring_matrix: ScoringMatrix
    gap: int
    begin_gap: int 
    end_gap: int
    use_begin_gap_top: bool = True
    use_begin_gap_left: bool = True
    use_end_gap_bottom: bool = True
    use_end_gap_right: bool = True
    
    def set_param(self, scoring_matrix, gap, begin_gap, end_gap):
        self.scoring_matrix = scoring_matrix
        self.gap = gap
        self.begin_gap = begin_gap
        self.end_gap = end_gap
        
    def get_score(self, residue1, residue2):
        return self.scoring_matrix.get_score(residue1, residue2)

@dataclass
class Node:
    value: int = 0
    traceback_i: int = 0
    traceback_j: int = 0

@dataclass
class AlignPairQuad:
    sequence1: str
    sequence2: str
    score_set: ScoreSet
    nodes: List[Node] = field(default_factory=list)
    score: int = 0
    aligned_sequence1: str = ''
    aligned_sequence2: str = ''

    def initialize_matrix(self):
        self.nodes = [[Node() for _ in range(len(self.sequence2) + 1)] for _ in range(len(self.sequence1) + 1)]
        self.nodes[0][0].value = 0

        # Initialize first column
        for i in range(1, len(self.nodes)):
            if self.score_set.use_begin_gap_left:
                self.nodes[i][0].value = self.nodes[i - 1][0].value - self.score_set.begin_gap
            else:
                self.nodes[i][0].value = self.nodes[i - 1][0].value - self.score_set.gap
            self.nodes[i][0].traceback_i = i - 1
            self.nodes[i][0].traceback_j = 0

        # Initialize first row
        for j in range(1, len(self.nodes[0])):
            if self.score_set.use_begin_gap_top:
                self.nodes[0][j].value = self.nodes[0][j - 1].value - self.score_set.begin_gap
            else:
                self.nodes[0][j].value = self.nodes[0][j - 1].value - self.score_set.gap
            self.nodes[0][j].traceback_i = 0
            self.nodes[0][j].traceback_j = j - 1

    def fill_matrix(self):
        for i in range(1, len(self.nodes)):
            for j in range(1, len(self.nodes[0])):
                a, b, c = 0, 0, 0
                
                # Handle end gaps
                if i == len(self.nodes) - 1 and j == len(self.nodes[0]) - 1:
                    if self.score_set.use_end_gap_right:
                        a = self.nodes[i - 1][j].value - self.score_set.end_gap
                    else:
                        a = self.nodes[i - 1][j].value - self.score_set.gap

                    if self.score_set.use_end_gap_bottom:
                        b = self.nodes[i][j - 1].value - self.score_set.end_gap
                    else:
                        b = self.nodes[i][j - 1].value - self.score_set.gap
                elif i == len(self.nodes) - 1:
                    a = self.nodes[i - 1][j].value - self.score_set.gap
                    if self.score_set.use_end_gap_bottom:
                        b = self.nodes[i][j - 1].value - self.score_set.end_gap
                    else:
                        b = self.nodes[i][j - 1].value - self.score_set.gap
                elif j == len(self.nodes[0]) - 1:
                    if self.score_set.use_end_gap_right:
                        a = self.nodes[i - 1][j].value - self.score_set.end_gap
                    else:
                        a = self.nodes[i - 1][j].value - self.score_set.gap
                    b = self.nodes[i][j - 1].value - self.score_set.gap
                else:
                    a = self.nodes[i - 1][j].value - self.score_set.gap
                    b = self.nodes[i][j - 1].value - self.score_set.gap

                c = self.nodes[i - 1][j - 1].value + self.score_set.get_score(self.sequence1[i - 1], self.sequence2[j - 1])

                if a >= b and a >= c:
                    self.nodes[i][j].value = a
                    self.nodes[i][j].traceback_i = i - 1
                    self.nodes[i][j].traceback_j = j
                elif b >= c and b >= a:
                    self.nodes[i][j].value = b
                    self.nodes[i][j].traceback_i = i
                    self.nodes[i][j].traceback_j = j - 1
                else:
                    self.nodes[i][j].value = c
                    self.nodes[i][j].traceback_i = i - 1
                    self.nodes[i][j].traceback_j = j - 1
        
        self.score = self.nodes[-1][-1].value

    def align(self):
        current_i = len(self.nodes) - 1
        current_j = len(self.nodes[0]) - 1
        current_node = self.nodes[len(self.nodes) - 1][len(self.nodes[0]) - 1]
        
        while current_node.traceback_i != 0 and current_node.traceback_j != 0:
            if current_node.traceback_i == current_i - 1 and current_node.traceback_j == current_j - 1:
                self.aligned_sequence1 += self.sequence1[current_i - 1]
                self.aligned_sequence2 += self.sequence2[current_j - 1]
            elif current_node.traceback_j == current_j - 1:
                self.aligned_sequence1 += '-'
                self.aligned_sequence2 += self.sequence2[current_j - 1]
            else:
                self.aligned_sequence1 += self.sequence1[current_i - 1]
                self.aligned_sequence2 += '-'
            
            current_i = current_node.traceback_i
            current_j = current_node.traceback_j

            current_node = self.nodes[current_node.traceback_i][current_node.traceback_j]
        
        self.aligned_sequence1 = self.aligned_sequence1[::-1] if self.aligned_sequence1 != '' else self.aligned_sequence1
        self.aligned_sequence2 = self.aligned_sequence2[::-1] if self.aligned_sequence2 != '' else self.aligned_sequence2

@dataclass
class AlignPairLinear:
    sequence1: str
    sequence2: str
    score_set: ScoreSet
    score: int = 0
    aligned_sequence1: str = ''
    aligned_sequence2: str = ''
    Sn: List[int] = field(default_factory=list)
    Sp: List[int] = field(default_factory=list)
    
    def align(self):
        if len(self.sequence1) == 0:
            for i in range(1, len(self.sequence2) + 1):
                self.aligned_sequence1 += '-'
                self.aligned_sequence2 += self.sequence2[i - 1]
                self.score += self.score_set.gap
        elif len(self.sequence2) == 0:
            for i in range(1, len(self.sequence1) + 1):
                self.aligned_sequence1 += self.sequence1[i - 1]
                self.aligned_sequence2 += '-'
                self.score += self.score_set.gap
        else:
            self.Sn = [0] * (len(self.sequence2) + 1)
            self.Sp = [0] * (len(self.sequence2) + 1)
            self.path(0, 0, len(self.sequence1), len(self.sequence2))
            
    def path(self, i1, j1, i2, j2):
        if i1 + 1 == i2 or j1 == j2:
            sub_sequence1 = self.sequence1[i1:i2]
            sub_sequence2 = self.sequence2[j1:j2]
            sub_score_set = ScoreSet(
                scoring_matrix=self.score_set.scoring_matrix,
                gap=self.score_set.gap,
                begin_gap=self.score_set.begin_gap,
                end_gap=self.score_set.end_gap
                )
            
            if j1 == j2:
                if j1 == 0:
                    sub_score_set.set_param(self.score_set.scoring_matrix,
                                            self.score_set.begin_gap,
                                            self.score_set.begin_gap,
                                            self.score_set.begin_gap
                                            )
                elif j1 == len(self.sequence2):
                    sub_score_set.set_param(self.score_set.scoring_matrix,
                                            self.score_set.end_gap,
                                            self.score_set.end_gap,
                                            self.score_set.end_gap
                                            )
                else:
                    sub_score_set.set_param(self.score_set.scoring_matrix,
                                            self.score_set.gap,
                                            self.score_set.gap,
                                            self.score_set.gap
                                            )
            else:
                sub_score_set.set_param(self.score_set.scoring_matrix,
                                        self.score_set.gap,
                                        self.score_set.begin_gap,
                                        self.score_set.end_gap
                                        )
                sub_score_set.use_begin_gap_top = False
                sub_score_set.use_begin_gap_left = False
                sub_score_set.use_end_gap_bottom = False
                sub_score_set.use_end_gap_right = False
                if i1 == 0:
                    sub_score_set.use_begin_gap_top = True
                if j1 == 0:
                    sub_score_set.use_begin_gap_left = True
                if j2 == len(self.sequence2):
                    sub_score_set.use_end_gap_right = True
                if i2 == len(self.sequence1):
                    sub_score_set.use_end_gap_bottom = True
            alignment = AlignPairQuad(
                sub_sequence1, 
                sub_sequence2, 
                sub_score_set
                )
            alignment.initialize_matrix()
            alignment.fill_matrix()
            alignment.align()
            self.aligned_sequence1 += alignment.aligned_sequence1
            self.aligned_sequence2 += alignment.aligned_sequence2
            self.score += alignment.score
        else:
            middle = (i1 + i2) // 2
            self.Sn = [0] * (len(self.sequence2) + 1)
            self.Sp = [0] * (len(self.sequence2) + 1)
            
            # Linear-space computation of alignment score to middle row (forward pass)
            if i1 == 0:
                self.Sn = [-self.score_set.begin_gap * j for j in range(len(self.Sn))]
            else:
                self.Sn = [-self.score_set.gap * j for j in range(len(self.Sn))]

            diag = 0
            for i in range(i1 + 1, middle + 1):
                left = self.Sn[j1]
                if j1 == 0:
                    left -= self.score_set.begin_gap
                else:
                    left -= self.score_set.gap

                self.Sn[j1] = left

                for j in range(j1 + 1, j2 + 1):
                    if j == len(self.sequence2) and i == len(self.sequence1):
                        left = max(
                            self.Sn[j] - self.score_set.end_gap,
                            max(left - self.score_set.end_gap, diag + self.score_set.get_score(self.sequence1[i - 1], self.sequence2[j - 1]))
                            )
                    elif i == len(self.sequence1):
                        left = max(
                            self.Sn[j] - self.score_set.gap,
                            max(left - self.score_set.end_gap, diag + self.score_set.get_score(self.sequence1[i - 1], self.sequence2[j - 1]))
                            )
                    elif j == len(self.sequence2):
                        left = max(
                            self.Sn[j] - self.score_set.end_gap,
                            max(left - self.score_set.gap, diag + self.score_set.get_score(self.sequence1[i - 1], self.sequence2[j - 1]))
                            )
                    else:
                        left = max(
                            self.Sn[j] - self.score_set.gap,
                            max(left - self.score_set.gap, diag + self.score_set.get_score(self.sequence1[i - 1], self.sequence2[j - 1]))
                            )
                    diag = self.Sn[j]
                    self.Sn[j] = left

            # Linear-space computation of alignment score to middle row (reverse pass)
            if i2 == len(self.sequence1):
                self.Sp = [-self.score_set.end_gap * j for j in range(len(self.Sp))]
            else:
                self.Sp = [-self.score_set.gap * j for j in range(len(self.Sp))]

            diag = 0
            for i in range(i2 - 1, middle - 1, -1):
                right = self.Sp[j2]
                if j2 == len(self.sequence2):
                    right -= self.score_set.end_gap
                else:
                    right -= self.score_set.gap

                self.Sp[j2] = right

                for j in range(j2 - 1, j1 - 1, -1):
                    if j == 0 and i == 0:
                        right = max(
                            self.Sp[j] - self.score_set.begin_gap,
                            max(right - self.score_set.begin_gap, diag + self.score_set.get_score(self.sequence1[i], self.sequence2[j]))
                        )
                    elif j == 0:
                        right = max(
                            self.Sp[j] - self.score_set.begin_gap,
                            max(right - self.score_set.gap, diag + self.score_set.get_score(self.sequence1[i], self.sequence2[j]))
                        )
                    elif i == 0:
                        right = max(
                            self.Sp[j] - self.score_set.gap,
                            max(right - self.score_set.begin_gap, diag + self.score_set.get_score(self.sequence1[i], self.sequence2[j]))
                        )
                    else:
                        right = max(
                            self.Sp[j] - self.score_set.gap,
                            max(right - self.score_set.gap, diag + self.score_set.get_score(self.sequence1[i], self.sequence2[j]))
                        )
                    diag = self.Sp[j]
                    self.Sp[j] = right

            # Find the value of j that maximizes self.Sn[j] + self.Sp[j]
            max_value = self.Sn[j1] + self.Sp[j1]
            max_j = j1

            for j in range(j1 + 1, j2 + 1):
                if self.Sn[j] + self.Sp[j] >= max_value:
                    max_value = self.Sn[j] + self.Sp[j]
                    max_j = j

            self.path(i1, j1, middle, max_j)
            self.path(middle, max_j, i2, j2)
    
    def get_aligned_sequence1(self):
        return self.aligned_sequence1

    def get_aligned_sequence2(self):
        return self.aligned_sequence2

# Example usage
sequence1 = "ATCGATCG"
sequence2 = "ATCATCGATCG"
scoring_matrix = ScoringMatrix(match=2, mismatch=-1)
score_set = ScoreSet(scoring_matrix=scoring_matrix, gap=2, begin_gap=0, end_gap=0)

# Create an instance of AlignPairQuad
aligner = AlignPairLinear(sequence1, sequence2, score_set, [], 0)

# Perform alignment
aligner.align()
print("Linear Alignment:")
print("Aligned Sequence 1:", aligner.get_aligned_sequence1())
print("Aligned Sequence 2:", aligner.get_aligned_sequence2())
print("Alignment Score:", aligner.score)
