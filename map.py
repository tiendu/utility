import hashlib
from typing import List, Tuple, Generator
from collections import Counter
from dataclasses import dataclass
import concurrent.futures
from functools import partial
from itertools import product
from math import sqrt
import os
import csv
import gzip
from pathlib import Path
import argparse
import sys

# Constants
FASTQ_EXTENSIONS = ['.fastq', '.fq']
FASTA_EXTENSIONS = ['.fasta', '.fa', '.fna', '.faa']

@dataclass(frozen=False)
class Seq:
    id: str
    sequence: str
    quality: str = ''
    
    def __post_init__(self):
        object.__setattr__(self, 'sequence', self.sequence.upper())

def read_sequences(file_path: Path) -> List[Seq]:
    sequences = []
    file_type = None

    if any(file_path.suffix in ext for ext in FASTQ_EXTENSIONS):
        file_type = 'FASTQ'
    elif any(file_path.suffix in ext for ext in FASTA_EXTENSIONS):
        file_type = 'FASTA'
    else:
        raise ValueError(f'Unrecognized file extension for {file_path}. Expected FASTA {FASTA_EXTENSIONS} or FASTQ {FASTQ_EXTENSIONS}')

    opener = gzip.open if file_path.suffix == '.gz' else open

    if file_type == 'FASTQ':
        with opener(file_path, 'rt') as fin:
            while True:
                header_line = fin.readline().strip()
                if not header_line:
                    break
                sequence_line = fin.readline().strip()
                _ = fin.readline().strip()  # Plus line
                quality_line = fin.readline().strip()
                sequences.append(Seq(header_line[1:], sequence_line, quality_line))
    elif file_type == 'FASTA':
        with opener(file_path, 'rt') as fin:
            seq_id, seq = None, []
            for line in fin:
                line = line.strip()
                if line.startswith('>'):
                    if seq_id:
                        sequences.append(Seq(seq_id, ''.join(seq)))
                    seq_id, seq = line[1:], []
                else:
                    seq.append(line.upper())
            if seq_id:
                sequences.append(Seq(seq_id, ''.join(seq)))

    return sequences

def hash_string(string: str, hash_function=hashlib.sha3_256) -> str:
    return hash_function(string.encode()).hexdigest()

def generate_hashed_kmers(string: str, k: int) -> Generator[str, None, None]:
    for i in range(len(string) - k + 1):
        yield hash_string(string[i:i + k])

def cosine_similarity(query_kmers: list[str], reference_kmers: list[str]) -> float:
    query_kmers = Counter(query_kmers)
    reference_kmers = Counter(reference_kmers)
    intersection = set(query_kmers) & set(reference_kmers)
    dot_product = sum(query_kmers[kmer] * reference_kmers[kmer] for kmer in intersection)
    norm_query = sqrt(sum(count**2 for count in query_kmers.values()))
    norm_reference = sqrt(sum(count**2 for count in reference_kmers.values()))
    
    if norm_query == 0 or norm_reference == 0:
        return 0.0

    return dot_product / (norm_query * norm_reference)

def map_query_to_reference(query: Seq, 
                           reference: Seq, 
                           similarity_threshold: float, 
                           coverage_threshold: float) -> List[Tuple[str, str, str, float, float]]:
    k = max(len(query.sequence) // 5 | 1, 3)
    query_kmers = list(generate_hashed_kmers(query.sequence, k))
    results = []

    for i in range(len(reference.sequence) - len(query.sequence) + 1):
        subref = reference.sequence[i:i + len(query.sequence)]
        subref_kmers = list(generate_hashed_kmers(subref, k))
        similarity = cosine_similarity(query_kmers, subref_kmers)
        coverage = len(subref) / len(reference.sequence)

        if similarity >= similarity_threshold and coverage >= coverage_threshold:
            results.append((query.id, reference.id, f'{i}..{i+len(query.sequence)}', round(similarity, 3), round(coverage, 3)))

    return results

def process_concurrently(query_sequences: List[Seq], 
                         reference_sequences: List[Seq], 
                         similarity_threshold: float, 
                         coverage_threshold: float, 
                         output_file: str,
                         threads: int) -> None:
    results = []

    max_workers = threads or 1
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        func = partial(map_query_to_reference, similarity_threshold=similarity_threshold, coverage_threshold=coverage_threshold)
        futures = [executor.submit(func, query, reference) for query, reference in product(query_sequences, reference_sequences)]
        for future in concurrent.futures.as_completed(futures):
            result = future.result()
            if result:
                results.extend(result)

    if results:
        headers = ['Query_ID', 'Reference_ID', 'Region', 'Similarity', 'Coverage']
        with open(output_file, 'w', newline='') as csvfile:
            csv_writer = csv.writer(csvfile)
            csv_writer.writerow(headers)
            csv_writer.writerows(results)

def main():
    parser = argparse.ArgumentParser(description="Sequence comparison using k-mers.")
    parser.add_argument("--query", required=True, help="Path to the query input file")
    parser.add_argument("--reference", required=True, help="Path to the reference input file")
    parser.add_argument("-s", "--similarity", type=float, default=0.8, help="Similarity threshold for sequence matching (default: 0.5)")
    parser.add_argument("-c", "--coverage", type=float, default=0.0, help="Coverage threshold for sequence matching (default: 0.0)")
    parser.add_argument("-t", "--threads", type=int, help="Number of threads")
    parser.add_argument("-o", "--output", required=True, help="Path to the output CSV file")
    
    args = parser.parse_args()

    query_sequences = read_sequences(Path(args.query))
    reference_sequences = read_sequences(Path(args.reference))

    if not query_sequences or not reference_sequences:
        print("Error: One or both input files are empty.")
        sys.exit(1)

    process_concurrently(query_sequences, reference_sequences, args.similarity, args.coverage, args.output, args.threads)

if __name__ == "__main__":
    main()
