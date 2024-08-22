from typing import List, Tuple, Dict, Set
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
from collections import defaultdict

# Constants
FASTQ_EXTENSIONS = ['.fastq', '.fq']
FASTA_EXTENSIONS = ['.fasta', '.fa', '.fna', '.faa']
K = 11  # Fixed length for k-mers
W = 3 * K  # Fixed window size

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

    if any(file_path.suffix in ext for ext in FASTQ_EXTENSIONS):
        file_type = 'FASTQ'
    elif any(file_path.suffix in ext for ext in FASTA_EXTENSIONS):
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
    rev_trans = []
    for aa in sequence:
        if aa in aa_to_nu:
            rev_trans.append(aa_to_nu[aa])
        else:
            raise ValueError(f"Unexpected amino acid '{aa}' encountered in sequence.")
    return ''.join(rev_trans)

def delineate(dna: str) -> List[str]:
    conversion = {
        'A': ['A'],
        'C': ['C'],
        'G': ['G'],
        'T': ['T'],
        'R': ['A', 'G'],
        'Y': ['C', 'T'],
        'S': ['G', 'C'],
        'W': ['A', 'T'],
        'K': ['G', 'T'],
        'M': ['A', 'C'],
        'B': ['C', 'G', 'T'],
        'D': ['A', 'G', 'T'],
        'H': ['A', 'C', 'T'],
        'V': ['A', 'C', 'G'],
        'N': ['A', 'T', 'G', 'C'],
    }

    all_variants = ['']
    for nucleotide in dna:
        all_variants = [prefix + base for prefix in all_variants for base in conversion.get(nucleotide, [nucleotide])]
    
    return all_variants

def get_minimizers(sequence: str, k: int, w: int) -> Dict[str, Set[int]]:
    minimizers = defaultdict(set)
    n = len(sequence)
    
    for i in range(n - w + 1):
        window = sequence[i:i + w]
        k_mers = [(window[j:j + k], i + j) for j in range(w - k + 1)]
        min_kmer, min_pos = min(k_mers)
        minimizers[min_kmer].add(min_pos)
    
    return minimizers

def match_minimizers(query_minimizers: Dict[str, Set[int]], reference_minimizers: Dict[str, Set[int]], threshold: float) -> float:
    query_keys = set(query_minimizers.keys())
    reference_keys = set(reference_minimizers.keys())
    
    intersection = query_keys & reference_keys
    union = query_keys | reference_keys
    
    if not union:
        return 0.0
    
    similarity = len(intersection) / len(union)
    return similarity

def map_query_to_reference(query: Seq, reference: Seq, k: int, w: int, threshold: float) -> List[Tuple[str, str, int, float]]:
    matches = []
    
    # Reverse translate if the query or reference is an amino acid sequence
    query_seq = reverse_translate(query.sequence) if query.id.endswith("aa") else query.sequence
    reference_seq = reverse_translate(reference.sequence) if reference.id.endswith("aa") else reference.sequence
    
    # Delineate sequences
    query_variants = delineate(query_seq)
    reference_variants = delineate(reference_seq)
    
    for query_variant in query_variants:
        for reference_variant in reference_variants:
            reference_minimizers = get_minimizers(reference_variant, k, w)
            query_minimizers = get_minimizers(query_variant, k, w)
            
            similarity = match_minimizers(query_minimizers, reference_minimizers, threshold)
            
            if similarity >= threshold:
                for kmer in query_minimizers:
                    if kmer in reference_minimizers:
                        for ref_pos in reference_minimizers[kmer]:
                            matches.append((query.id, reference.id, ref_pos, similarity))
    
    return matches

def process_concurrently(query_sequences: List[Seq], reference_sequences: List[Seq], query_type: str, reference_type: str, k: int, w: int, similarity_threshold: float, output_file: str) -> None:
    results = []

    max_workers = os.cpu_count() or 1
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        func = partial(map_query_to_reference, k=k, w=w, threshold=similarity_threshold)
        futures = [executor.submit(func, query, reference) for query, reference in product(query_sequences, reference_sequences)]
        for future in concurrent.futures.as_completed(futures):
            result = future.result()
            if result:
                results.extend(result)

    if results:
        headers = ['Query_ID', 'Reference_ID', 'Reference_Pos', 'Similarity']
        with open(output_file, 'w', newline='') as csvfile:
            csv_writer = csv.writer(csvfile)
            csv_writer.writerow(headers)
            csv_writer.writerows(results)

def main():
    parser = argparse.ArgumentParser(description="Sequence comparison using minimizers.")
    parser.add_argument("--query", required=True, help="Path to the query input file")
    parser.add_argument("--reference", required=True, help="Path to the reference input file")
    parser.add_argument("--query_type", choices=["nu", "aa"], required=True, help="Type of the query input file (nu for nucleotide, aa for amino acid)")
    parser.add_argument("--reference_type", choices=["nu", "aa"], required=True, help="Type of the reference input file (nu for nucleotide, aa for amino acid)")
    parser.add_argument("-t", "--threshold", type=float, default=0.8, help="Similarity threshold for sequence matching (default: 0.8)")
    parser.add_argument("-o", "--output", required=True, help="Path to the output CSV file")
    
    args = parser.parse_args()

    query_sequences = read_sequences(Path(args.query))
    reference_sequences = read_sequences(Path(args.reference))

    if not query_sequences or not reference_sequences:
        print("Error: One or both input files are empty.")
        sys.exit(1)
    
    # Sort sequences by length in descending order
    query_sequences.sort(key=lambda x: len(x.sequence), reverse=True)
    reference_sequences.sort(key=lambda x: len(x.sequence), reverse=True)
    
    process_concurrently(query_sequences, reference_sequences, args.query_type, args.reference_type, K, W, args.threshold, args.output)

if __name__ == "__main__":
    main()
