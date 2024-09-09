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

def reverse_complement(dna: str) -> str:
    complement = str.maketrans('ATGCRYSWKMBDHVN', 'TACGYRSWMKVHDBN')
    return dna[::-1].translate(complement)

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

def euclidean_similarity(query_kmers: List[str], reference_kmers: List[str]) -> float:
    query_kmers = Counter(query_kmers)
    reference_kmers = Counter(reference_kmers)
    all_kmers = set(query_kmers.keys()).union(set(reference_kmers.keys()))    
    sum_of_squares = sum((query_kmers[kmer] - reference_kmers[kmer]) ** 2 for kmer in all_kmers)
    euclidean_dist = sqrt(sum_of_squares)
    max_possible_dist = sqrt(len(all_kmers))
    if max_possible_dist == 0:
        return 0.0
    norm_dist = euclidean_dist / max_possible_dist  # Normalize the distance to be in the range [0, 1]
    norm_similarity = 1 - norm_dist  # Convert to a similarity measure: similarity = 1 - normalized distance
    return norm_similarity

def truncate_string(string: str, width: int) -> str:
    return string if len(string) <= width else string[:width - 3] + '...'

def print_stdout_table(headers: List[str], rows: List[List[str]], col_widths: dict):
    def format_cell(value: str, width: int) -> str:
        truncated = value if len(value) <= width else value[:width - 3] + '...'
        return truncated.ljust(width)  # Left-justify to ensure even spacing

    header_row = " | ".join([format_cell(header, col_widths[header]) for header in headers])
    print(f"| {header_row} |")
    divider_row = "-" * (sum(col_widths.values()) + 3 * (len(headers) - 1) + 2)
    print(divider_row)
    for row in rows:
        formatted_row = " | ".join([format_cell(str(value), col_widths[header]) for value, header in zip(row, headers)])
        print(f"| {formatted_row} |")

def map_query_to_reference(query: Seq, 
                           reference: Seq, 
                           similarity_threshold: float, 
                           coverage_threshold: float,
                           is_nucleotide: bool,
                           similarity_func=euclidean_similarity) -> List[Tuple[str, str, str, float, float, str]]:
    if len(query.sequence) > len(reference.sequence):
        query, reference = reference, query
    k = max(len(query.sequence) // 7 | 1, 3)
    query_kmers = list(generate_hashed_kmers(query.sequence, k))
    results = []
    for i in range(len(reference.sequence) - len(query.sequence) + 1):
        subref = reference.sequence[i:i + len(query.sequence)]
        subref_kmers = list(generate_hashed_kmers(subref, k))
        if is_nucleotide:
            rev_complement_seq = reverse_complement(subref)
            rev_complement_kmers = list(generate_hashed_kmers(rev_complement_seq, k))
            fwd_similarity = similarity_func(query_kmers, subref_kmers)
            rev_similarity = similarity_func(query_kmers, rev_complement_kmers)
            if fwd_similarity > rev_similarity:
                best_similarity = fwd_similarity
                strand = '+'
            else:
                best_similarity = rev_similarity
                strand = '-'
        else:  # Amino acid mode: only forward strand
            best_similarity = similarity_func(query_kmers, subref_kmers)
            strand = '+'  # No reverse complement for amino acids
        coverage = len(subref) / len(reference.sequence)
        if best_similarity >= similarity_threshold and coverage >= coverage_threshold:
            results.append((query.id, reference.id, f'{i}..{i+len(query.sequence)}', round(best_similarity, 3), round(coverage, 3), strand))
    return results

def process_concurrently(query_sequences: List[Seq], 
                         reference_sequences: List[Seq], 
                         similarity_threshold: float, 
                         coverage_threshold: float, 
                         is_nucleotide: bool,
                         output_file: str,
                         threads: int,
                         method: str) -> None:
    results = []
    if method == 'euclid':
        similarity_func = euclidean_similarity
    elif method == 'cosine':
        similarity_func = cosine_similarity
    max_workers = threads or 1
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        func = partial(map_query_to_reference, 
                       similarity_threshold=similarity_threshold, 
                       coverage_threshold=coverage_threshold, 
                       is_nucleotide=is_nucleotide,
                       similarity_func=similarity_func)
        futures = [executor.submit(func, query, reference) for query, reference in product(query_sequences, reference_sequences)]
        for future in concurrent.futures.as_completed(futures):
            result = future.result()
            if result:
                results.extend(result)
    headers = ['Query_ID', 'Reference_ID', 'Region', 'Similarity', 'Coverage', 'Strand']
    col_widths = {
        'Query_ID': 20,
        'Reference_ID': 20,
        'Region': 15,
        'Similarity': 10,
        'Coverage': 10,
        'Strand': 5
    }
    if results:
        print_stdout_table(headers, results, col_widths)
        with open(output_file, 'w', newline='') as csvfile:
            csv_writer = csv.writer(csvfile)
            csv_writer.writerow(headers)
            csv_writer.writerows(results)

def main():
    parser = argparse.ArgumentParser(description='Sequence comparison using k-mers.')
    parser.add_argument('--query', required=True, help='Path to the query input file')
    parser.add_argument('--reference', required=True, help='Path to the reference input file')
    parser.add_argument('-s', '--similarity', type=float, default=0.8, help='Similarity threshold for sequence matching (default: 0.8)')
    parser.add_argument('-c', '--coverage', type=float, default=0.0, help='Coverage threshold for sequence matching (default: 0.0)')
    parser.add_argument('-t', '--threads', type=int, default=4, help='Number of threads')
    parser.add_argument('-o', '--output', required=True, help='Path to the output CSV file')
    parser.add_argument('--mode', choices=['nu', 'aa'], default='nu', 
                        help="Comparison mode: 'nu' for DNA/RNA, 'aa' for proteins (default: 'nu')")
    parser.add_argument('--method', choices=['euclid', 'cosine'], default='cosine', help='Method for calculating similarity')
    args = parser.parse_args()
    query_sequences = read_sequences(Path(args.query))
    reference_sequences = read_sequences(Path(args.reference))
    if not query_sequences or not reference_sequences:
        print('Error: One or both input files are empty.')
        sys.exit(1)
    is_nucleotide = args.mode == 'nu'
    process_concurrently(query_sequences, reference_sequences, args.similarity, args.coverage, is_nucleotide, args.output, args.threads, args.method)

if __name__ == '__main__':
    main()
