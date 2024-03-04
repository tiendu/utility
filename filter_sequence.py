import argparse
import os
import gzip
import logging
from collections import Counter
from itertools import groupby
from dataclasses import dataclass
from typing import List, Dict
import concurrent.futures

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Constants
FASTQ_EXTENSIONS = ['.fastq', '.fq']
FASTA_EXTENSIONS = ['.fasta', '.fa', '.fna']
LINE_LENGTH = 80

# Data class to represent a sequence.
@dataclass(frozen=True)
class Sequence:
    id: str
    sequence: str
    quality: str = ''

def read_sequences_from_file(file_path: str, file_type: str) -> List[Sequence]:
    sequences = []

    # Determine the appropriate file opening mode based on the file extension.
    opener = gzip.open if file_path.endswith('.gz') else open
    
    if file_type == "FASTQ":
        with opener(file_path, 'rt') as f:
            fqiter = groupby(enumerate(f), key=lambda x: x[0] // 4)
            for _, fq in fqiter:
                header, sequence, _, quality = [line.strip() for _, line in fq]
                sequence_id = header[1:]  # remove "@" character
                sequence = sequence.upper()
                quality = quality
                sequences.append(Sequence(sequence_id, sequence, quality))
    elif file_type == "FASTA":
        with opener(file_path, 'rt') as f:
            faiter = (x[1] for x in groupby(f, lambda line: line[0] == ">"))
            for header in faiter:
                header = next(header).strip()
                sequence_id = header[1:]  # remove ">" character
                sequence = "".join(string.strip().upper() for string in next(faiter))
                sequences.append(Sequence(sequence_id, sequence))
    return sequences

def write_sequences_to_file(sequences: List[Sequence], file_path: str) -> None:
    with open(file_path, 'w') as f:
        for sequence in sequences:
            if sequence.quality == '':
                f.write(f">{sequence.id}\n")
                # Write sequence with a defined line length
                for i in range(0, len(sequence.sequence), LINE_LENGTH):
                    f.write(sequence.sequence[i:i+80] + "\n")
            else:
                f.write(f"@{sequence.id}\n{sequence.sequence}\n+\n{sequence.quality}\n")

def generate_kmers(sequence: str, k: int) -> List[str]:
    return [sequence[i:i+k] for i in range(len(sequence) - k + 1)]

def generate_kmer_counts(sequence: str, k: int) -> Dict[str, int]:
    kmers = generate_kmers(sequence, k)
    return Counter(kmers)

def generate_kmer_matrix_for_chunk(chunk: List[Sequence], k: int) -> Dict[str, Dict[str, int]]:
    kmer_matrix: Dict[str, Dict[str, int]] = {}
    for sequence in chunk:
        kmer_counts = generate_kmer_counts(sequence.sequence, k)
        kmer_matrix[sequence.id] = kmer_counts
    return kmer_matrix

def filter_sequences(sequences: List[Sequence], k: int, num_threads: int) -> List[Sequence]:
    total_size = len(sequences)
    chunk_size = max(1, total_size // num_threads)
    chunks = [sequences[i:i + chunk_size] for i in range(0, total_size, chunk_size)]

    results: List[Dict[str, Dict[str, int]]] = []

    with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as executor:
        futures = [executor.submit(generate_kmer_matrix_for_chunk, chunk, k) for chunk in chunks]
        concurrent.futures.wait(futures)
        for future in futures:
            results.append(future.result())

    # Combine k-mer counts from all chunks
    kmer_counts: Dict[str, Dict[str, int]] = {}
    for result in results:
        for sequence_id, counts in result.items():
            for kmer, count in counts.items():
                kmer_counts.setdefault(kmer, {}).setdefault(sequence_id, 0)
                kmer_counts[kmer][sequence_id] += count

    # Calculate total counts for each k-mer
    total_kmer_counts = {kmer: sum(counts.values()) for kmer, counts in kmer_counts.items()}

    # Filter out k-mers with total count <= 1
    kmer_matrix = {kmer: counts for kmer, counts in kmer_counts.items() if total_kmer_counts[kmer] > 1 and len(counts) > 1}

    # Extract filtered sequences based on k-mer matrix
    filtered_sequence_ids = set.union(*[set(counts.keys()) for counts in kmer_matrix.values()])
    filtered_sequences = [sequence for sequence in sequences if sequence.id in filtered_sequence_ids]

    return filtered_sequences

def main():
    parser = argparse.ArgumentParser(description='Generate k-mer matrix from FASTA file.')
    parser.add_argument('-i', '--input_file', required=True, help='Path to the input FASTA file.')
    parser.add_argument('-o', '--output_file', required=True, help='Path to the output file for the k-mer matrix.')
    parser.add_argument('-k', '--k_mer', type=int, default=11, help='Length of k-mer.')
    parser.add_argument('-t', '--num_threads', type=int, default=4, help='Number of threads to use.')
    args = parser.parse_args()

    if not os.path.exists(args.input_file):
        logging.error(f"Input file '{args.input_file}' not found.")
        return
    
    # Check file type.
    file_type = ''
    if any(ext in args.input_file for ext in FASTQ_EXTENSIONS):
        file_type = "FASTQ"
    elif any(ext in args.input_file for ext in FASTA_EXTENSIONS):
        file_type = "FASTA"
    else:
        raise ValueError(f"Unrecognized file extension for {args.input_file}. Expected FASTA (.fasta, .fa, .fna) or FASTQ (.fastq, .fq).")

    # Validate CPUs.
    max_threads = os.cpu_count()
    if args.num_threads < 1 or args.num_threads > max_threads:
        logging.warning(f"Invalid number of threads. Adjusting thread count to be between 1 and {max_threads}")
        args.num_threads = min(max(args.num_threads, 1), max_threads)
    if args.num_threads >= len(read_sequences_from_file(args.input_file, file_type)) * 0.1:
        logging.warning(f"Number of sequences too low, thread count adjusted to 1")
        args.num_threads = 1

    sequences = read_sequences_from_file(args.input_file, file_type)
    filtered_sequences = filter_sequences(sequences, args.k_mer, args.num_threads)
    write_sequences_to_file(filtered_sequences, args.output_file)
    
if __name__ == '__main__':
    main()
