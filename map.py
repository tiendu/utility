import hashlib
from typing import List, Tuple, Generator
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

def reverse_translate(sequence: str) -> str:
    aa_to_nu = {
        'W': 'TGG', 'Y': 'TAY', 'C': 'TGY', 'E': 'GAR',
        'K': 'AAR', 'Q': 'CAR', 'S': 'WSN', 'L': 'YTN',
        'R': 'MGN', 'G': 'GGN', 'F': 'TTY', 'D': 'GAY',
        'H': 'CAY', 'N': 'AAY', 'M': 'ATG', 'A': 'GCN',
        'P': 'CCN', 'T': 'ACN', 'V': 'GTN', 'I': 'ATH',
        '*': 'TRR', 'X': 'NNN'
    }

    return ''.join(aa_to_nu.get(aa, 'NNN') for aa in sequence)

def delineate(dna: str) -> List[str]:
    conversion = {
        'A': ['A'], 'C': ['C'], 'G': ['G'], 'T': ['T'],
        'R': ['A', 'G'], 'Y': ['C', 'T'], 'S': ['G', 'C'],
        'W': ['A', 'T'], 'K': ['G', 'T'], 'M': ['A', 'C'],
        'B': ['C', 'G', 'T'], 'D': ['A', 'G', 'T'],
        'H': ['A', 'C', 'T'], 'V': ['A', 'C', 'G'],
        'N': ['A', 'T', 'G', 'C'],
    }

    variants = ['']
    for nucleotide in dna:
        variants = [prefix + base for prefix in variants for base in conversion.get(nucleotide, [nucleotide])]
    
    return variants

def expand_sequences(sequences: List[Seq]) -> List[Seq]:
    expanded_sequences = []

    for seq in sequences:
        delineated_seqs = delineate(seq.sequence)
        expanded_sequences.extend([Seq(seq.id, variant) for variant in delineated_seqs])

    return expanded_sequences

def hash_string(string: str, hash_function=hashlib.sha3_256) -> str:
    return hash_function(string.encode()).hexdigest()

def generate_hashed_kmers(string: str, k: int) -> Generator[str, None, None]:
    for i in range(len(string) - k + 1):
        yield hash_string(string[i:i + k])

def cosine_similarity_ordered(query_kmers: List[str], reference_kmers: List[str]) -> float:
    dot_product = sum(1 for q_kmer, r_kmer in zip(query_kmers, reference_kmers) if q_kmer == r_kmer)
    norm_query = sqrt(len(query_kmers))
    norm_reference = sqrt(len(reference_kmers))
    
    if norm_query == 0 or norm_reference == 0:
        return 0.0
    
    return dot_product / (norm_query * norm_reference)

def map_query_to_reference(query: Seq, reference: Seq, similarity_threshold: float, coverage_threshold: float) -> List[Tuple[str, str, str, float, float]]:
    k = max(len(query.sequence) // 5 | 1, 3)
    query_kmers = list(generate_hashed_kmers(query.sequence, k))
    results = []

    for i in range(len(reference.sequence) - len(query.sequence) + 1):
        subref = reference.sequence[i:i + len(query.sequence)]
        subref_kmers = list(generate_hashed_kmers(subref, k))
        similarity = cosine_similarity_ordered(query_kmers, subref_kmers)
        coverage = len(subref) / len(reference.sequence)

        if similarity >= similarity_threshold and coverage >= coverage_threshold:
            results.append((query.id, reference.id, f'{i}..{i+len(query.sequence)}', similarity, coverage))

    return results 

def process_concurrently(query_sequences: List[Seq], 
                         reference_sequences: List[Seq], 
                         similarity_threshold: float, 
                         coverage_threshold: float, 
                         output_file: str) -> None:
    results = []

    max_workers = os.cpu_count() - 1 or 1
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
    parser.add_argument("--query_type", choices=["nu", "aa"], required=True, help="Type of the query input file (nu for nucleotide, aa for amino acid)")
    parser.add_argument("--reference_type", choices=["nu", "aa"], required=True, help="Type of the reference input file (nu for nucleotide, aa for amino acid)")
    parser.add_argument("-t", "--threshold", type=float, default=0.5, help="Similarity threshold for sequence matching (default: 0.5)")
    parser.add_argument("-c", "--coverage", type=float, default=0.0, help="Coverage threshold for sequence matching (default: 0.0)")
    parser.add_argument("-o", "--output", required=True, help="Path to the output CSV file")
    
    args = parser.parse_args()

    query_sequences = read_sequences(Path(args.query))
    reference_sequences = read_sequences(Path(args.reference))

    if not query_sequences or not reference_sequences:
        print("Error: One or both input files are empty.")
        sys.exit(1)

    if args.query_type == 'aa':
        query_sequences = [Seq(seq.id, reverse_translate(seq.sequence)) for seq in query_sequences]
    if args.reference_type == 'aa':
        reference_sequences = [Seq(seq.id, reverse_translate(seq.sequence)) for seq in reference_sequences]

    query_sequences = expand_sequences(query_sequences)
    reference_sequences = expand_sequences(reference_sequences)
    process_concurrently(query_sequences, reference_sequences, args.threshold, args.coverage, args.output)

if __name__ == "__main__":
    main()
