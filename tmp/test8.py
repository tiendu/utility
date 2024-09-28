
import argparse
import os
import gzip
import logging
import hashlib
import random
from functools import partial
from itertools import groupby
from dataclasses import dataclass
from typing import List, Generator, Dict, Set, Callable, Any
from concurrent.futures import ProcessPoolExecutor, wait, as_completed

# Constants
FASTQ_EXTENSIONS = ['.fastq', '.fq']
FASTA_EXTENSIONS = ['.fasta', '.fa', '.fna', '.faa']
LINE_LENGTH = 80
CHUNK_SIZE = 100_000

# Setup logging to display messages with INFO level and above.
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

@dataclass(frozen=True)
class Seq:
    id: str
    sequence: str
    quality: str = ''
    def __hash__(self) -> int:
        return hash((self.id, self.sequence, self.quality))

    def length(self) -> int:
        return len(self.sequence)

def read_sequences_from_file(file_path: str, file_type: str) -> List[Seq]:
    sequences = []
    opener = gzip.open if file_path.endswith('.gz') else open
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

def write_sequences_to_file(sequences: List[Seq], file_path: str) -> None:
    opener = gzip.open if file_path.endswith('.gz') else open
    with opener(file_path, 'wt') as f:
        for seq in sequences:
            if seq.quality == '':
                f.write(f'>{seq.id}\n')
                for i in range(0, len(seq.sequence), LINE_LENGTH):
                    f.write(seq.sequence[i:i+LINE_LENGTH] + '\n')
            else:
                f.write(f'@{seq.id}\n{seq.sequence}\n+\n{seq.quality}\n')

def kmerize(string: str, k: int, modifier=None) -> list[str]:
    if len(string) < k:
        return []
    if modifier:
        return [modifier(string[i:i+k]) for i in range(len(string) - k + 1)]
    return [string[i:i+k] for i in range(len(string) - k + 1)]

def hash_string(string: str, hash_function=hashlib.sha3_256) -> str:
    return hash_function(string.encode()).hexdigest()

def sequence_to_bits(dna: str) -> list:
    nucl_to_bits = {
        'A': 0b0001,  # A = 1
        'T': 0b0010,  # T = 2
        'G': 0b0100,  # G = 4
        'C': 0b1000,  # C = 8
        'W': 0b0011,  # W = A | T = 3
        'R': 0b0101,  # R = A | G = 5
        'Y': 0b1010,  # Y = C | T = 10
        'S': 0b1100,  # S = G | C = 12
        'K': 0b0110,  # K = G | T = 6
        'M': 0b1001,  # M = A | C = 9
        'B': 0b1110,  # B = C | G | T = 14
        'D': 0b0111,  # D = A | G | T = 7
        'H': 0b1011,  # H = A | C | T = 11
        'V': 0b1101,  # V = A | C | G = 13
        'N': 0b1111   # N = A | C | G | T = 15 (any nucleotide)
    }
    return [nucl_to_bits[nu.upper()] for nu in dna]

def compare_two_strings(str1: str, str2: str, k: int=3) -> bool:
    if not k or k > min(len(str1), len(str2)):
        k = min(len(str1), len(str2))
    bits1 = kmerize(str1, k, sequence_to_bits)
    bits2 = kmerize(str2, k, sequence_to_bits)
    def compare_bit_sets(set1: list, set2: list) -> bool:
        return all(bit1 & bit2 for bit1, bit2 in zip(set1, set2))

    if not bits1 or not bits2:
        return False  # Return False for consistency in boolean values
    return any(compare_bit_sets(bit1, bit2)
               for bit1 in bits1 for bit2 in bits2)

def round_robin_divide(items: List[Any], 
                       chunk_size: int, 
                       num_threads: int, 
                       key: Callable[[Any], Any], 
                       is_descending: bool) -> List[List[Any]]:
    def is_sorted(lst: List[Any], comparison_func: Callable[[Any, Any], bool]) -> bool:
        return all(comparison_func(lst[i], lst[i + 1]) for i in range(len(lst) - 1))
    
    if is_descending:
        order = lambda x, y: x >= y
    else:
        order = lambda x, y: x <= y
    if not is_sorted([key(item) for item in items], order):
        items = sorted(items, key=key, reverse=is_descending)
    chunks = [[] for _ in range(num_threads)]
    for i, item in enumerate(items):
        chunks[i % num_threads].append(item)
    divided_chunks = []
    for chunk in chunks:
        for i in range(0, len(chunk), chunk_size):
            divided_chunks.append(chunk[i:i + chunk_size])
    return divided_chunks

def deduplicate_chunk(sequences: List[str], mode: str) -> List[str]:
    sequences = sorted(sequences, key=lambda seq: seq.length(), reverse=True)
    k = sequences[-1].length()  # Shortest sequence length
    uniq_seqs = {}
    uniq_kmers = set()
    for seq in sequences:
        kmer_hashes = set(hashed_kmer for hashed_kmer in kmerize(seq.sequence, k, hash_string))
        if all(kmer_hash in uniq_kmers for kmer_hash in kmer_hashes):
            continue
        if mode == 'nu':
             if not any(compare_two_strings(seq.sequence, unique_seq.sequence, k)
                        for unique_seq in uniq_seqs.values()):
                uniq_seqs[seq.sequence] = seq
                uniq_kmers.update(kmer_hashes)
        elif mode == 'aa':
            uniq_seqs[seq.sequence] = seq
            uniq_kmers.update(kmer_hashes)
    return uniq_seqs

def deduplicate_concurrently(sequences: List[Seq], num_threads: int, mode: str) -> List[Seq]:
    deduped_seqs = dict()
#    while sequences:
#        logging.info(f'Current number of sequences: {len(sequences)}')
#        size = max(1, len(sequences) // num_threads)  # Calculate chunk size
#        chunks = [sequences[i:i + size] for i in range(0, len(sequences), size)]  # Proper chunking
#        with ProcessPoolExecutor(max_workers=num_threads) as executor:
#            func = partial(deduplicate_chunk,
#                           mode=mode)
#            futures = [executor.submit(func, chunk) for chunk in chunks]
#            wait(futures)
#            for future in as_completed(futures):
#                uniq_seqs = future.result()
#                for seq, obj in uniq_seqs.items():
#                    if not any(seq in deduped_seq or deduped_seq in seq for deduped_seq in deduped_seqs):
#                        deduped_seqs[seq] = obj
#        if len(deduped_seqs) == len(sequences):
#            break
#        sequences = list(deduped_seqs.values())
#        random.shuffle(sequences)
#    logging.info(f'Final number of sequences: {len(sequences)}')
    sequences = list(deduplicate_chunk(sequences, mode).values())
    return sequences

def main() -> None:
    parser = argparse.ArgumentParser(description='Deduplicate FASTA/FASTQ sequences.')
    parser.add_argument('-i', '--input_file', required=True, help='Path to the input file.')
    parser.add_argument('-o', '--output_file', required=True, help='Path to the output file.')
    parser.add_argument('-t', '--num_threads', type=int, default=4, help='Number of threads to use. Default is 4.')
    parser.add_argument('-m', '--mode', choices=['nu', 'aa'], default='nu',
                        help="Comparison mode: 'nu' for DNA/RNA, 'aa' for proteins (default: 'nu')")
    args = parser.parse_args()
    if not os.path.exists(args.input_file):
        raise ValueError(f'The input file "{args.input_file}" does not exist.')
    try:
        with open(args.output_file, 'w'):
            pass
    except Exception as e:
        raise ValueError(f'Error creating output file "{args.output_file}": {e}')
    file_type = ''
    if any(ext in args.input_file for ext in FASTQ_EXTENSIONS):
        file_type = 'FASTQ'
    elif any(ext in args.input_file for ext in FASTA_EXTENSIONS):
        file_type = 'FASTA'
    else:
        raise ValueError(f'Unrecognized file extension for {args.input_file}. Expected FASTA (.fasta, .fa, .fna) or FASTQ (.fastq, .fq).')
    max_threads = os.cpu_count()
    if args.num_threads < 1 or args.num_threads > max_threads:
        logging.warning(f'Invalid number of threads. Adjusting thread count to be between 1 and {max_threads}')
        args.num_threads = min(max(args.num_threads, 1), max_threads)
    sequences = read_sequences_from_file(args.input_file, file_type)
    if not sequences:
        raise ValueError('No sequences detected!')
    if args.num_threads >= len(sequences) * 0.1:
        logging.warning(f'Number of sequences too low, thread count adjusted to 1')
        args.num_threads = 1
    deduplicated_sequences = deduplicate_concurrently(sequences, args.num_threads, args.mode)
    write_sequences_to_file(deduplicated_sequences, args.output_file)

if __name__ == '__main__':
    main()
