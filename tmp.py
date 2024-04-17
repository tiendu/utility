import sys
import re
from dataclasses import dataclass
from typing import List


@dataclass
class ScoringMatrix:
    match: int
    mismatch: int

@dataclass
class Identity(ScoringMatrix):
    pass

@dataclass
class ScoreSet:
    scoring_matrix: ScoringMatrix
    gap: int
    begin_gap: int 
    end_gap: int
    use_begin_gap_top = True
    use_begin_gap_left = True
    use_end_gap_bottom = True
    use_end_gap_right = True
    
    def set_param(self, scoring_matrix, gap, begin_gap, end_gap):
        self.scoring_matrix = scoring_matrix
        self.gap = gap
        self.begin_gap = begin_gap
        self.end_gap = end_gap

@dataclass
class Node:
    value: int
    traceback_i: int
    trackback_j: int

@dataclass
class AlignPairQuad:
    sequence1: str
    sequence2: str
    score_set: ScoreSet
    nodes: List[Node]
    aligned_sequence1: str
    aligned_sequence2: str
    score: int
    def initialize_matrix(self):
        self.score_set = score_set
        
        

@dataclass
class AlignPairLinear:
    sequence1: str
    sequence2: str
    score_set: ScoreSet
    aligned_sequence1: str
    aligned_sequence2: str
    score: int
    
    def align(self):
        self.aligned_sequence1 = ''
        self.aligned_sequence2 = ''
        score = 0
        if len(self.sequence1) == 0:
            for i in range(1, len(self.sequence2) + 1):
                aligned_sequence1 += '-'
                aligned_sequence2 += self.sequence2[i - 1]
                score += self.score_set.gap
        elif len(self.sequence2) == 0:
            for i in range(1, len(self.sequence1) + 1):
                aligned_sequence1 += self.sequence1[i - 1]
                aligned_sequence2 += '-'
                score += self.score_set.gap
        else:
            self.path(0, 0, len(self.sequence1), len(self.sequence2))
            
    def path(self, i1, j1, i2, j2):
        if i1 + 1 == i2 or j1 == j2:
            sub_sequence1 = ''
            sub_sequence2 = ''
            for i in range(i1 + 1, i2):
                sub_sequence1 += self.sequence1[i - 1]
            for j in range(j1 + 1, j2):
                sub_sequence2 += self.sequence2[j - 1]
            sub_score_set = ScoreSet()
            if j1 == j2:
                if j1 == 0:
                    sub_score_set.set_param(self.score_set.scoring_matrix,
                                            self.score_set.begin_gap,
                                            self.score_set.begin_gap,
                                            self.score_set.begin_gap)
                elif j1 == len(self.sequence2):
                    sub_score_set.set_param(self.score_set.scoring_matrix,
                                            self.score_set.end_gap,
                                            self.score_set.end_gap,
                                            self.score_set.end_gap)
                else:
                    sub_score_set.set_param(self.score_set.scoring_matrix,
                                            self.score_set.gap,
                                            self.score_set.gap,
                                            self.score_set.gap)
            else:
                sub_score_set.set_param(self.score_set.scoring_matrix,
                                        self.score_set.gap,
                                        self.score_set.begin_gap,
                                        self.score_set.end_gap)
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

# Example usage
scoring_matrix = ScoringMatrix(match=2, mismatch=-1)
score_set = ScoreSet(scoring_matrix=scoring_matrix, gap=2, begin_gap=0, end_gap=0)
aligner = AlignPairLinear(sequence1="ATCGATCGATCG", sequence2="ATCATCGATCG", score_set=score_set)
aligner.align()

print(aligner)
