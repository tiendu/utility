from typing import List, Tuple
from dataclasses import dataclass
import concurrent.futures
import os
from functools import partial
from itertools import product, groupby
import sys
import csv
import gzip
from pathlib import Path
import argparse

# Constants
FASTQ_EXTENSIONS = ['.fastq', '.fq']
FASTA_EXTENSIONS = ['.fasta', '.fa', '.fna', 'faa']

@dataclass(frozen=True)
class Seq:
    id: str
    sequence: str
    quality: str = ''
    
    def __post_init__(self):
        object.__setattr__(self, 'sequence', self.sequence.upper())

def read_sequences(file_path: Path) -> List[Seq]:
    sequences = []
    file_type = ''

    if any(file_path.suffix == ext for ext in FASTQ_EXTENSIONS):
        file_type = 'FASTQ'
    elif any(file_path.suffix == ext for ext in FASTA_EXTENSIONS):
        file_type = 'FASTA'
    else:
        raise ValueError(f'Unrecognized file extension for {file_path}. Expected FASTA {FASTA_EXTENSIONS} or FASTQ {FASTQ_EXTENSIONS}')

    opener = gzip.open if file_path.suffix == '.gz' else open

    if file_type == 'FASTQ':
        with opener(file_path, 'rt') as fin:
            groups = groupby(enumerate(fin), key=lambda x: x[0] // 4)
            for _, group in groups:
                header_line, sequence_line, _, quality_line = [line.strip() for _, line in group]
                name = header_line[1:]
                seq = sequence_line.upper()
                qual = quality_line.upper()
                sequences.append(Seq(name, seq, qual))
    elif file_type == 'FASTA':
        with opener(file_path, 'rt') as fin:
            faiter = (x[1] for x in groupby(fin, lambda line: line[0] == '>'))
            for header in faiter:
                headerstr = next(header).strip()
                name = headerstr[1:]
                seq = ''.join(s.strip().upper() for s in next(faiter))
                sequences.append(Seq(name, seq))

    return sequences

def reverse_translate(sequence: str) -> str:
    aa_to_nu = {
        'W': 'TGG', 'Y': 'TAY', 'C': 'TGY', 'E': 'GAR',
        'K': 'AAR', 'Q': 'CAR', 'S': 'WSN', 'L': 'YTN',
        'R': 'MGN', 'G': 'GGN', 'F': 'TTY', 'D': 'GAY',
        'H': 'CAY', 'N': 'AAY', 'M': 'ATG', 'A': 'GCN',
        'P': 'CCN', 'T': 'ACN', 'V': 'GTN', 'I': 'ATH',
        '*': 'TRR', 'X': 'NNN'
    }
    trans_table = str.maketrans(aa_to_nu)
    try:
        rev_trans = sequence.translate(trans_table)
    except KeyError as e:
        raise ValueError(f"Invalid amino acid '{e.args[0]}' in sequence. Valid amino acids: {list(aa_to_nu.keys())}") from e
    return rev_trans

def get_minimizer(sequence: str, k: int, w: int) -> List[Tuple[str, int]]:
    minimizers = []
    n = len(sequence)
    
    for i in range(n - w + 1):
        window = sequence[i:i + w]
        k_mers = [(window[j:j + k], i + j) for j in range(w - k + 1)]
        min_kmer = min(k_mers)
        minimizers.append(min_kmer)

    # Deduplicate consecutive identical minimizers
    unique_minimizers = []
    for m in minimizers:
        if not unique_minimizers or m[0] != unique_minimizers[-1][0]:
            unique_minimizers.append(m)
    
    return unique_minimizers

def match_minimizers(minimizers1: List[Tuple[str, int]], minimizers2: List[Tuple[str, int]], similarity_threshold: float) -> float:
    set1 = set(minimizer[0] for minimizer in minimizers1)
    set2 = set(minimizer[0] for minimizer in minimizers2)
    
    intersection = set1.intersection(set2)
    union = set1.union(set2)
    
    if not union:
        return 0.0
    
    similarity = len(intersection) / len(union)
    
    if similarity > similarity_threshold:
        return similarity

def match_pair(sequence1: Seq, sequence2: Seq, sequence1_type: str, sequence2_type: str, k: int, w: int, similarity_threshold: float) -> List:
    matches = []

    if sequence1_type == 'aa':
        sequence1.sequence = reverse_translate(sequence1.sequence)
    elif sequence1_type != 'nu':
        raise ValueError('Incorrect input 1 type!')

    if sequence2_type == 'aa':
        sequence2.sequence = reverse_translate(sequence2.sequence)
    elif sequence2_type != 'nu':
        raise ValueError('Incorrect input 2 type!')

    minimizers1 = get_minimizer(sequence1.sequence, k, w)
    minimizers2 = get_minimizer(sequence2.sequence, k, w)

    similarity = match_minimizers(minimizers1, minimizers2, similarity_threshold)
    if similarity:
        matches.append((sequence1.id, sequence2.id, f'{similarity * 100:.1f}%'))

    return matches

def calculate_dynamic_params(sequence: str, default_k: int, default_w: int) -> Tuple[int, int]:
    length = len(sequence)
    
    if length < default_w:
        w = max(3, (length // 2) | 1)  # Ensure w is odd and at least 3
        k = max(5, (w // 2) | 1)  # Ensure k is odd and at least 5
    else:
        w = default_w
        k = default_k
    
    return k, w

def process_concurrently(sequences1: List[Seq], sequences2: List[Seq], sequences1_type: str, sequences2_type: str, k: int, w: int, similarity_threshold: float, output_file: str) -> None:
    results = []

    max_workers = os.cpu_count() or 1
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        func = partial(match_pair, sequence1_type=sequences1_type, sequence2_type=sequences2_type, k=k, w=w, similarity_threshold=similarity_threshold)
        futures = [executor.submit(func, sequence1, sequence2) for sequence1, sequence2 in product(sequences1, sequences2)]
        for future in concurrent.futures.as_completed(futures):
            result = future.result()
            if result:
                results.extend(result)

    if results:
        headers = ['Query_ID', 'Reference_ID', 'Details']
        with open(output_file, 'w', newline='') as csvfile:
            csv_writer = csv.writer(csvfile)
            csv_writer.writerow(headers)
            csv_writer.writerows(results)

def main():
    parser = argparse.ArgumentParser(description="Sequence comparison using minimizers.")
    parser.add_argument("input1", help="Path to the first input file")
    parser.add_argument("input2", help="Path to the second input file")
    parser.add_argument("input1_type", choices=["nu", "aa"], help="Type of the first input file (nu for nucleotide, aa for amino acid)")
    parser.add_argument("input2_type", choices=["nu", "aa"], help="Type of the second input file (nu for nucleotide, aa for amino acid)")
    parser.add_argument("-k", type=int, default=31, help="Length of k-mers for minimizer calculation (default: 31)")
    parser.add_argument("-w", type=int, default=51, help="Window size for minimizer calculation (default: 51)")
    parser.add_argument("-t", "--threshold", type=float, default=0.8, help="Similarity threshold for sequence matching (default: 0.8)")
    parser.add_argument("-o", "--output", help="Path to the output CSV file")
    
    args = parser.parse_args()

    sequences1 = read_sequences(args.input1)
    sequences2 = read_sequences(args.input2)

    if not sequences1 or not sequences2:
        print("Error: One or both input files are empty.")
        sys.exit(1)
    
    # Sort sequences by length in descending order
    sequences1.sort(key=lambda x: len(x.sequence), reverse=True)
    sequences2.sort(key=lambda x: len(x.sequence), reverse=True)

    # Determine dynamic parameters for both sequences
    k1, w1 = calculate_dynamic_params(sequences1[0].sequence, args.k, args.w)
    k2, w2 = calculate_dynamic_params(sequences2[0].sequence, args.k, args.w)
    
    # Use the smallest k and w values from both sequences
    k = min(k1, k2)
    w = min(w1, w2)
    
    process_concurrently(sequences1, sequences2, args.input1_type, args.input2_type, k, w, args.threshold, args.output)

if __name__ == "__main__":
    main()
