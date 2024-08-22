
import hashlib
from typing import List, Tuple, Dict, Generator, Set
from dataclasses import dataclass
import concurrent.futures
import os
from functools import partial
from itertools import product, groupby
import sys
import csv
import gzip
from pathlib import Path
from collections import Counter
from math import sqrt
import argparse

# Constants
FASTQ_EXTENSIONS = ['.fastq', '.fq']
FASTA_EXTENSIONS = ['.fasta', '.fa', '.fna', '.faa']
K = 11  # Fixed length for k-mers

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

def hash_string(string: str, hash_function=hashlib.sha3_256) -> str:
    return hash_function(string.encode()).hexdigest()

def generate_hashed_kmers(string: str, k: int) -> Generator[str, None, None]:
    for i in range(len(string) - k + 1):
        kmer = string[i:i + k]
        yield hash_string(kmer)

def calculate_normalized_similarity(query_kmers: Set[str], reference_kmers: Set[str], query_len: int, reference_len: int) -> float:
    """Calculate normalized similarity based on k-mer sets."""
    intersection = query_kmers & reference_kmers
    union = query_kmers | reference_kmers
    if not union:
        return 0.0
    jaccard_index = len(intersection) / len(union)
    
    # Normalize based on the length of the sequences
    query_normalization = len(query_kmers) / query_len
    reference_normalization = len(reference_kmers) / reference_len
    
    # Return the average of Jaccard index and the normalized values
    return (jaccard_index + query_normalization + reference_normalization) / 3

def cosine_similarity(query_kmers: Counter, reference_kmers: Counter) -> float:
    """Calculate cosine similarity between k-mer frequency vectors."""
    intersection = set(query_kmers) & set(reference_kmers)
    dot_product = sum(query_kmers[kmer] * reference_kmers[kmer] for kmer in intersection)
    norm_query = sqrt(sum(count**2 for count in query_kmers.values()))
    norm_reference = sqrt(sum(count**2 for count in reference_kmers.values()))
    
    if norm_query == 0 or norm_reference == 0:
        return 0.0
    
    return dot_product / (norm_query * norm_reference)

def map_query_to_reference(query: Seq, reference: Seq, k: int, threshold: float) -> List[Tuple[str, str, float]]:
    query_kmers = Counter(generate_hashed_kmers(query.sequence, k))
    reference_kmers = Counter(generate_hashed_kmers(reference.sequence, k))
    
    # similarity = cosine_similarity(query_kmers, reference_kmers)
    similarity = calculate_normalized_similarity(query_kmers, reference_kmers, len(query.sequence), len(reference.sequence))
    
    if similarity >= threshold:
        return [(query.id, reference.id, similarity)]
    
    return []

def process_concurrently(query_sequences: List[Seq], 
                         reference_sequences: List[Seq], 
                         query_type: str, reference_type: str, 
                         k: int, similarity_threshold: float, 
                         output_file: str) -> None:
    results = []

    max_workers = os.cpu_count() or 1
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        func = partial(map_query_to_reference, k=k, threshold=similarity_threshold)
        futures = [executor.submit(func, query, reference) for query, reference in product(query_sequences, reference_sequences)]
        for future in concurrent.futures.as_completed(futures):
            result = future.result()
            if result:
                results.extend(result)

    if results:
        headers = ['Query_ID', 'Reference_ID', 'Similarity']
        with open(output_file, 'w', newline='') as csvfile:
            csv_writer = csv.writer(csvfile)
            csv_writer.writerow(headers)
            csv_writer.writerows(results)

def main():
    parser = argparse.ArgumentParser(description="Sequence comparison using k-mers.")
    parser.add_argument("--query", required=True, help="Path to the query input file")
    parser.add_argument("--reference", required=True, help="Path to the reference input file")
    parser.add_argument("--query_type", choices=["nu", "aa"], required=True, help="Type of the query input file (nu for nucleotide, aa for amino acid)")
    parser.add_argument("--reference_type", choices=["nu", "aa"], required=True, help="Type of the reference input file (nu for nucleotide, aa for amino acid)")
    parser.add_argument("-t", "--threshold", type=float, default=0.01, help="Similarity threshold for sequence matching (default: 0.8)")
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
    
    process_concurrently(query_sequences, reference_sequences, args.query_type, args.reference_type, K, args.threshold, args.output)

if __name__ == "__main__":
    main()
